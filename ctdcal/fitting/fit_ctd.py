"""
A module for fitting CTD discrete/bottle and continuous data.
"""
import logging
from pathlib import Path

import gsw
import numpy as np
import pandas as pd

from ctdcal.common import validate_dir
from ctdcal.fitting.fit_common import get_node, NodeNotFoundError, multivariate_fit, apply_polyfit
from ctdcal import get_ctdcal_config
from ctdcal import process_ctd as process_ctd
from ctdcal.processors.functions_oxy import calculate_dV_dt
from ctdcal.processors.functions_salt import CR_to_cond
from ctdcal.plotting.plot_fit import _intermediate_residual_plot
# from ctdcal.fitting.fit_legacy import load_fit_yaml
from ctdcal.flagging.flag_common import _flag_btl_data, outliers, nan_values, by_residual

cfg = get_ctdcal_config()
log = logging.getLogger(__name__)


def _prepare_fit_data(df, param, ref_param, zRange=None):
    """Remove non-finite data, trim to desired zRange, and remove extreme outliers"""

    good_data = df[np.isfinite(df[ref_param]) & np.isfinite(df[param])].copy()

    if zRange is not None:
        zMin, zMax = zRange.split(":")
        good_data = good_data[
            (good_data["CTDPRS"] > int(zMin)) & (good_data["CTDPRS"] < int(zMax))
        ]

    good_data["Diff"] = good_data[ref_param] - good_data[param]
    good_data["Flag"] = outliers(good_data["Diff"], n_sigma2=4)

    df_good = good_data[good_data["Flag"] == 2].copy()
    df_bad = good_data[good_data["Flag"] == 4].copy()

    return df_good, df_bad


def calibrate_ct(sensor, fit_groups, btl):
    # unimplemented... coming back to this
    #
    T_flag, T_fit_coefs = pd.DataFrame(), pd.DataFrame()
    for group in fit_groups:
        btl_rows = btl["cast_id"].isin(group).values
        good_rows = btl_rows & (btl_df["REFTMP_FLAG_W"] == 2)
        time_rows = time_df["SSSCC"].isin(ssscc_sublist).values


def calibrate_temp(btl_df, time_df, fit_groups, fit_coeffs, outdir, report_dir, cast_id='SSSCC'):
    """
    Least-squares fit CTD temperature data against reference data.

    Parameters
    -----------
    btl_df : DataFrame
        CTD data at bottle stops
    time_df : DataFrame
        Continuous CTD data

    Returns
    --------

    """
    outdir = validate_dir(outdir, create=True)
    sensor_cols = ['CTDTMP1', 'CTDTMP2']
    for sensor_idx, sensor in sorted(enumerate(fit_coeffs)):  # ['primary', 'secondary'] because they are alphabetically in the order we want
        t_flag, t_fit_coefs = pd.DataFrame(), pd.DataFrame()
        for coeffs_idx, coeffs in enumerate(fit_coeffs[sensor]):
            casts = fit_groups[coeffs_idx]
            btl_rows = btl_df[cast_id].isin(casts).values
            good_rows = btl_rows & (btl_df['REFTMP_FLAG_W'] == 2)
            time_rows = time_df[cast_id].isin(casts).values

            # 1) plot pre-fit residual
            _intermediate_residual_plot(
                btl_df.loc[btl_rows, 'REFTMP']
                - btl_df.loc[btl_rows, sensor_cols[sensor_idx]],
                btl_df.loc[btl_rows, 'CTDPRS'],
                btl_df.loc[btl_rows, cast_id],
                xlabel=f'{sensor_cols[sensor_idx].upper()} Residual (T90 C)',
                f_out=Path(outdir, 'temp_residual_%s_%s_prefit.pdf' % (sensor, coeffs_idx)),
            )

            # 2) prepare data for fitting
            # NOTE: df_bad will be overwritten during post-fit data flagging but is
            # left here for future debugging (if necessary)
            df_good, df_bad = _prepare_fit_data(
                btl_df[good_rows],
                sensor_cols[sensor_idx],
                'REFTMP',
                zRange=coeffs['%s_%s' % (sensor, coeffs_idx)].zRange
            )
            _intermediate_residual_plot(
                df_good['Diff'],
                df_good['CTDPRS'],
                df_good[cast_id],
                xlabel=f"{sensor_cols[sensor_idx].upper()} Residual (T90 C)",
                f_out=Path(outdir, 'temp_residual_%s_%s_fit_data.pdf' % (sensor, coeffs_idx)),
            )

            # 3) calculate fit coefs
            P_order = coeffs['%s_%s' % (sensor, coeffs_idx)].P_order
            T_order = coeffs['%s_%s' % (sensor, coeffs_idx)].T_order
            coef_dict = multivariate_fit(
                df_good['Diff'],
                (df_good['CTDPRS'], P_order),
                (df_good[sensor_cols[sensor_idx]], T_order),
                coef_names=["cp", "ct"],
            )

            # 4) apply fit
            P_coefs = tuple(coef_dict[f"cp{n}"] for n in np.arange(1, P_order + 1))
            T_coefs = tuple(coef_dict[f"ct{n}"] for n in np.arange(1, T_order + 1))
            btl_df.loc[btl_rows, sensor_cols[sensor_idx]] = apply_polyfit(
                btl_df.loc[btl_rows, sensor_cols[sensor_idx]],
                (coef_dict["c0"],) + T_coefs,
                (btl_df.loc[btl_rows, 'CTDPRS'], P_coefs),
            )
            time_df.loc[time_rows, sensor_cols[sensor_idx]] = apply_polyfit(
                time_df.loc[time_rows, sensor_cols[sensor_idx]],
                (coef_dict["c0"],) + T_coefs,
                (time_df.loc[time_rows, 'CTDPRS'], P_coefs),
            )

            # 4.5) flag CTDTMP and make residual plots
            df_ques, df_bad = _flag_btl_data(
                btl_df[btl_rows],
                param=sensor_cols[sensor_idx],
                ref='REFTMP',
                cast_id=cast_id,
                f_out=Path(outdir, 'temp_residual_%s_%s.pdf' % (sensor, coeffs_idx)),
            )

            # 5) handle quality flags
            t_flag = pd.concat([t_flag, df_bad, df_ques])

            # 6) handle fit params
            coef_df = pd.DataFrame()
            coef_df[cast_id] = casts
            coef_names = ["cp2", "cp1", "ct2", "ct1", "c0"]
            coef_df[coef_names] = 0.0
            for k, v in coef_dict.items():
                coef_df[k] = v

            t_fit_coefs = pd.concat([t_fit_coefs, coef_df])

        # one more fig with all cuts
        _intermediate_residual_plot(
            btl_df['REFTMP'] - btl_df[sensor_cols[sensor_idx]],
            btl_df['CTDPRS'],
            btl_df[cast_id],
            xlabel=f"{sensor_cols[sensor_idx]} Residual (T90 C)",
            show_thresh=True,
            f_out=Path(outdir, 'temp_residual_%s_all_postfit.pdf' % sensor)
        )

        # export temp quality flags
        report_file = Path(report_dir, 'qual_flag_%s.csv' % sensor_cols[sensor_idx])
        t_flag.sort_index().to_csv(report_file, index=False)

        # export temp fit params (formated to 5 sig figs, scientific notation)
        t_fit_coefs[coef_names] = t_fit_coefs[coef_names].applymap(
            lambda x: np.format_float_scientific(x, precision=4, exp_digits=1)
        )
        report_file = Path(report_dir, 'fit_coef_%s.csv' % sensor_cols[sensor_idx])
        t_fit_coefs.to_csv(report_file, index=False)

    # flag temperature data
    time_df["CTDTMP_FLAG_W"] = 2

    return True

def calibrate_cond(btl_df, time_df, fit_groups, fit_coeffs, outdir, report_dir, flag_file, cast_id='SSSCC', ref_node='salt'):
    """
    Least-squares fit CTD conductivity data against bottle salts.

    Parameters
    -----------
    btl_df : DataFrame
        CTD data at bottle stops
    time_df : DataFrame
        Continuous CTD data
    user_cfg : Munch object
        Dictionary of user configuration parameters
    ref_node : str
        Name of reference parameter

    Returns
    --------

    """
    log.info("Calibrating conductivity")
    # calculate BTLCOND values from autosal data
    btl_df['BTLCOND'] = CR_to_cond(
        btl_df["CRavg"],
        btl_df["BathTEMP"],
        btl_df['CTDTMP1'],
        btl_df['CTDPRS'],
    )

    # merge in handcoded salt flags
    flag_file = Path(flag_file)
    salt_flags_manual = None
    if flag_file.exists():
        try:
            salt_flags_manual = get_node(flag_file, ref_node)
        except NodeNotFoundError:
            log.info("No previously flagged values for %s found in flag file." % ref_node)
    else:
        log.info("No pre-existing flag file found.")

    if salt_flags_manual is not None:
        log.info("Merging previously flagged values for %s." % ref_node)
        salt_flags_manual_df = pd.DataFrame.from_dict(salt_flags_manual)
        salt_flags_manual_df = salt_flags_manual_df.rename(
                columns={"cast_id": cast_id, "bottle_num": "btl_fire_num", "value": "SALNTY_FLAG_W"},
        ).drop(columns=["notes"])
        btl_df = btl_df.merge(salt_flags_manual_df, on=[cast_id, "btl_fire_num"], how="left")
        btl_df["SALNTY_FLAG_W"] = nan_values(
            btl_df["SALNTY"], old_flags=btl_df["SALNTY_FLAG_W"]
        )
    else:
        btl_df["SALNTY_FLAG_W"] = nan_values(btl_df["SALNTY"])

    # here lies the calibrating
    outdir = validate_dir(outdir, create=True)
    t_cols = ['CTDTMP1', 'CTDTMP2']
    c_cols = ['CTDCOND1', 'CTDCOND2']

    for sensor_idx, sensor in sorted(enumerate(fit_coeffs)):  # ['primary', 'secondary'] because they are alphabetically in the order we want
        c_flag, c_fit_coefs = pd.DataFrame(), pd.DataFrame()
        for coeffs_idx, coeffs in enumerate(fit_coeffs[sensor]):
            casts = fit_groups[coeffs_idx]
            btl_rows = btl_df[cast_id].isin(casts).values
            good_rows = btl_rows & (btl_df['SALNTY_FLAG_W'] == 2)
            time_rows = time_df[cast_id].isin(casts).values

            # 1) plot pre-fit residual
            _intermediate_residual_plot(
                btl_df.loc[btl_rows, 'BTLCOND']
                - btl_df.loc[btl_rows, c_cols[sensor_idx]],
                btl_df.loc[btl_rows, 'CTDPRS'],
                btl_df.loc[btl_rows, cast_id],
                xlabel=f'{c_cols[sensor_idx].upper()} Residual (mS/cm)',
                f_out=Path(outdir, 'cond_residual_%s_%s_prefit.pdf' % (sensor, coeffs_idx)),
            )

            # 2) prepare data for fitting
            # NOTE: df_bad will be overwritten during post-fit data flagging
            # but is left here for future debugging (if necessary)
            df_good, df_bad = _prepare_fit_data(
                    btl_df[good_rows],
                    c_cols[sensor_idx],
                    'BTLCOND',
                    zRange=coeffs['%s_%s' % (sensor, coeffs_idx)].zRange
            )
            _intermediate_residual_plot(
                    df_good['Diff'],
                    df_good['CTDPRS'],
                    df_good[cast_id],
                    xlabel=f"{c_cols[sensor_idx].upper()} Residual (mS/cm)",
                    f_out=Path(outdir, 'cond_residual_%s_%s_fit_data.pdf' % (sensor, coeffs_idx)),
            )
            # 3) calculate fit coefs
            P_order = coeffs['%s_%s' % (sensor, coeffs_idx)].P_order
            T_order = coeffs['%s_%s' % (sensor, coeffs_idx)].T_order
            C_order = coeffs['%s_%s' % (sensor, coeffs_idx)].C_order
            coef_dict = multivariate_fit(
                    df_good['Diff'],
                    (df_good['CTDPRS'], P_order),
                    (df_good[t_cols[sensor_idx]], T_order),
                    (df_good[c_cols[sensor_idx]], C_order),
                    coef_names=['cp', 'ct', 'cc'],
            )

            # 4) apply fit
            P_coefs = tuple(coef_dict[f"cp{n}"] for n in np.arange(1, P_order + 1))
            T_coefs = tuple(coef_dict[f"ct{n}"] for n in np.arange(1, T_order + 1))
            C_coefs = tuple(coef_dict[f"cc{n}"] for n in np.arange(1, C_order + 1))
            btl_df.loc[btl_rows, c_cols[sensor_idx]] = apply_polyfit(
                btl_df.loc[btl_rows, c_cols[sensor_idx]],
                (coef_dict["c0"],) + C_coefs,
                (btl_df.loc[btl_rows, 'CTDPRS'], P_coefs),
                (btl_df.loc[btl_rows, t_cols[sensor_idx]], T_coefs),
            )
            time_df.loc[time_rows, c_cols[sensor_idx]] = apply_polyfit(
                time_df.loc[time_rows, c_cols[sensor_idx]],
                (coef_dict["c0"],) + C_coefs,
                (time_df.loc[time_rows, 'CTDPRS'], P_coefs),
                (time_df.loc[time_rows, t_cols[sensor_idx]], T_coefs),
            )

            # 4.5) flag CTDCOND and make residual plots
            df_ques, df_bad = _flag_btl_data(
                btl_df[btl_rows],
                param=c_cols[sensor_idx],
                ref='BTLCOND',
                cast_id=cast_id,
                f_out=Path(outdir, 'cond_residual_%s_%s.pdf' % (sensor, coeffs_idx)),
            )

            # 5) handle quality flags
            c_flag = pd.concat([c_flag, df_bad, df_ques])

            # 6) handle fit params
            coef_df = pd.DataFrame()
            coef_df[cast_id] = casts
            coef_names = ['cp2', 'cp1', 'ct2', 'ct1', 'cc2', 'cc1', 'c0']
            coef_df[coef_names] = 0.0
            for k, v in coef_dict.items():
                coef_df[k] = v

            c_fit_coefs = pd.concat([c_fit_coefs, coef_df])

        # one more fig with all cuts
        _intermediate_residual_plot(
            btl_df['BTLCOND'] - btl_df[c_cols[sensor_idx]],
            btl_df['CTDPRS'],
            btl_df[cast_id],
            xlabel=f"{c_cols[sensor_idx]} Residual (mS/cm)",
            show_thresh=True,
            f_out=Path(outdir, 'cond_residual_%s_all_postfit.pdf' % sensor)
        )

        # export cond quality flags
        report_file = Path(report_dir, 'qual_flag_%s.csv' % c_cols[sensor_idx])
        c_flag.sort_index().to_csv(report_file, index=False)

        # export cond fit params
        c_fit_coefs[coef_names] = c_fit_coefs[coef_names].applymap(
            lambda x: np.format_float_scientific(x, precision=4, exp_digits=1)
        )
        report_file = Path(report_dir, 'fit_coef_%s.csv' % c_cols[sensor_idx])
        c_fit_coefs.to_csv(report_file, index=False)

    # recalculate salinity with calibrated C/T
    time_df['CTDSAL'] = gsw.SP_from_C(
        time_df['CTDCOND1'],
        time_df['CTDTMP1'],
        time_df['CTDPRS'],
    )
    btl_df['CTDSAL'] = gsw.SP_from_C(
        btl_df['CTDCOND1'],
        btl_df['CTDTMP1'],
        btl_df['CTDPRS'],
    )

    # flag salinity data
    time_df['CTDSAL_FLAG_W'] = 2
    btl_df['CTDSAL_FLAG_W'] = by_residual(
        btl_df['CTDSAL'],
        btl_df["SALNTY"],
        btl_df['CTDPRS'],
    )
    bad_rows = btl_df["SALNTY_FLAG_W"].isin([3, 4])
    btl_df.loc[bad_rows, 'CTDSAL_FLAG_W'] = 2  # bad salts not used for QC
    btl_df['CTDSAL_FLAG_W'] = nan_values(
        btl_df['CTDSAL'], old_flags=btl_df['CTDSAL_FLAG_W']
    )

    return btl_df, time_df


def _get_pressure_offset(start_vals, end_vals):
    """
    Finds unique values and calculate mean for pressure offset

    Parameters
    ----------
    start_vals : array_like
        Array of initial ondeck pressure values

    end_vals : array_like
        Array of ending ondeck pressure values
    Returns
    -------
    p_off : float
         Average pressure offset

    """
    log.warning("Use of _get_pressure_offset() is deprecated. Use calibrate_pressure() instead.")
    p_start = pd.Series(np.unique(start_vals))
    p_end = pd.Series(np.unique(end_vals))
    p_start = p_start[p_start.notnull()]
    p_end = p_end[p_end.notnull()]
    p_off = p_start.mean() - p_end.mean()

    # JACKSON THINKS THIS METHOD SHOULD BE USED TO KEEP START END PAIRS
    #    p_df = pd.DataFrame()
    #    p_df['p_start'] = p_start
    #    p_df['p_end'] = p_end
    #    p_df = p_df[p_df['p_end'].notnull()]
    #    p_df = p_df[p_df['p_start'].notnull()]
    #    p_off = p_df['p_start'].mean() - p_df['p_end'].mean()
    ##########################################################

    p_off = np.around(p_off, decimals=4)

    return p_off


def calibrate_pressure(btl_data, time_data, fit_groups, report_dir):
    offsets_file = Path(report_dir, 'ondeck_pressure.csv')
    p = pd.read_csv(offsets_file, dtype={"cast_id": str}, na_values="Started in Water")

    for group in fit_groups:
        # get mean offset by fit group
        offsets = p[['pressure_start', 'pressure_end']][p['cast_id'].isin(group)].mean().round(4)
        # TODO: save these to a pressure coeffs file?
        print('Pressure offset for casts %s-%s: start %s, end %s'
              % (group[0], group[-1], offsets['pressure_start'], offsets['pressure_end'])
              )

        # this is how we have historically computed the final offset
        # TODO: this doesn't seem scientifically sound to me. is it??
        offset = offsets['pressure_start'] - offsets['pressure_end']

        # should we do something like this instead?
        # offset = np.mean([offsets['pressure_start'], offsets['pressure_end']])

        # apply offsets to btl and time
        # btl_data['CTDPRS'] += offset
        # time_data['CTDPRS'] += offset

        # i think we should do this instead
        btl_data['CTDPRS'] -= offsets['pressure_start']
        time_data['CTDPRS'] -= offsets['pressure_start']


        # add flag cols
        btl_data['CTDPRS_FLAG_W'] = 2
        time_data['CTDPRS_FLAG_W'] = 2



def apply_pressure_offset(df, p_col="CTDPRS"):
    """
    Calculate pressure offset using deck pressure log and apply it to the data.
    Pressure flag column is added with value 2, indicating the data are calibrated.

    Parameters
    ----------
    df : DataFrame
        DataFrame containing column with pressure values
    p_col : str, optional
        Pressure column name in DataFrame (defaults to CTDPRS)

    Returns
    -------
    df : DataFrame
        DataFrame containing updated pressure values and a new flag column

    """
    log.warning("Use of apply_pressure_offset() is deprecated. Use calibrate_pressure() instead.")
    p_log = pd.read_csv(
        cfg.dirs["logs"] + "ondeck_pressure.csv",
        dtype={"cast_id": str},
        na_values="Started in Water",
    )
    p_offset = _get_pressure_offset(p_log['pressure_start'], p_log['pressure_end'])
    df[p_col] += p_offset
    df[p_col + "_FLAG_W"] = 2

    return df


def load_time_all(casts, time_dir):
    time_dir = Path(time_dir)
    pkl_all = sorted(tuple(time_dir.glob('*_time.pkl')))
    for i, pkl in enumerate(pkl_all):
        if str(pkl.stem)[:-5] not in casts:
            pkl_all.pop(i)
    # for cast in casts:
    #     log.info("Loading TIME data for station: " + cast + "...")
    #     time_file = Path(time_dir, "%s_time.pkl" % cast)
    #     time_data = pd.read_pickle(time_file)
    #     df_list.append(time_data)
    df_data_all = pd.concat(
            [pd.read_pickle(pkl) for pkl in pkl_all],
            axis=0,
            sort=False,
    )

    df_data_all['cast_id'] = df_data_all['cast_id'].astype('category')
    df_data_all.reset_index(inplace=True)

    return df_data_all

def load_all_ctd_files(ssscc_list):
    """
    Load CTD files for station/cast list and merge into a dataframe.

    Parameters
    ----------
    ssscc_list : list of str
        List of stations to load

    Returns
    -------
    df_data_all : DataFrame
        Merged dataframe containing all loaded data

    """
    log.warning("Use of load_all_ctd_files() is deprecated. Use load_time_all() instead.")
    df_list = []
    for ssscc in ssscc_list:
        log.info("Loading TIME data for station: " + ssscc + "...")
        time_file = cfg.dirs["time"] + ssscc + "_time.pkl"
        time_data = pd.read_pickle(time_file)
        time_data["SSSCC"] = str(ssscc)
        time_data["dv_dt"] = calculate_dV_dt(
            time_data["CTDOXYVOLTS"], time_data["scan_datetime"]
        )
        df_list.append(time_data)
        # print("** Finished TIME data station: " + ssscc + " **")
    df_data_all = pd.concat(df_list, axis=0, sort=False)

    df_data_all["master_index"] = range(len(df_data_all))

    return df_data_all


