import sys
from pathlib import Path

import numpy as np
import pandas as pd
from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import (
    BoxSelectTool,
    Button,
    ColumnDataSource,
    DataTable,
    Div,
    MultiSelect,
    Range1d,
    Select,
    StringFormatter,
    TableColumn,
    TextInput,
)
from bokeh.plotting import figure

from ctdcal import get_ctdcal_config, io
from ctdcal.common import load_user_config, validate_file
from ctdcal.fitting.common import df_node_to_BottleFlags, get_node, save_node

# cfg = get_ctdcal_config()
USERCONFIG = '/Users/als026/data/ices/ices.yaml'
user_cfg = load_user_config(validate_file(USERCONFIG))
INST = 'ctd'
FLAGFILE = Path(user_cfg.datadir, "fit", INST, user_cfg.bottleflags_man)

# TODO: abstract parts of this to a separate file
# TODO: following above, make parts reusable?

# load continuous CTD data and make into a dict (only ~20MB)
file_list = sorted(Path(user_cfg.datadir, 'export', INST).glob("*ct1.csv"))
ssscc_list = [ssscc.stem[:-4] for ssscc in file_list]
ctd_data = []
for f in file_list:
    print(f"Loading {f}")
    header, df = io.load_exchange_ctd(f)
    df["SSSCC"] = header["CASTNO"]
    # df["SSSCC"] = header["STNNBR"].zfill(3) + header["CASTNO"].zfill(2)
    ctd_data.append(df)
ctd_data = pd.concat(ctd_data, axis=0, sort=False)

# load bottle file
fname = list(Path(user_cfg.datadir, 'export', INST).glob("*hy1.csv"))[0]
btl_data = io.load_exchange_btl(fname).replace(-999, np.nan)
btl_data["SSSCC"] = btl_data["CASTNO"]

#   Create salinity residuals
btl_data["Residual"] = btl_data["SALNTY"] - btl_data["CTDSAL"]
btl_data[["CTDPRS", "Residual"]] = btl_data[["CTDPRS", "Residual"]].round(4)
btl_data["Comments"] = ""
btl_data["New Flag"] = btl_data["SALNTY_FLAG_W"].copy()

#   Temperature
btl_data["t_res"] = (btl_data["REFTMP"] - btl_data["CTDTMP"]).round(4)
btl_data["New T Flag"] = btl_data["REFTMP_FLAG_W"].copy()

#   Aaaand oxygen
if "OXYGEN" not in btl_data.columns:
    btl_data["OXYGEN"] = np.nan
    btl_data["OXYGEN_FLAG_W"] = 9
btl_data["o_res"] = (btl_data["OXYGEN"] - btl_data["CTDOXY"]).round(4)
btl_data["New O Flag"] = btl_data["OXYGEN_FLAG_W"].copy()

deltas_d = {"CTDSAL": "Residual", "CTDTMP": "t_res", "CTDOXY": "o_res"}

# update with old handcoded flags if file exists
if FLAGFILE.exists():
    salt_flags_manual = get_node(FLAGFILE, "salt")
    salt_flags_manual_df = pd.DataFrame.from_dict(salt_flags_manual)
    salt_flags_manual_df = salt_flags_manual_df.rename(
        columns={
            "value": "New Flag",
            "cast_id": "SSSCC",
            "bottle_num": "SAMPNO",
            "notes": "Comments",
        }
    )

    # there's gotta be a better way... but this is good enough for now
    btl_data = btl_data.merge(salt_flags_manual_df, on=["SSSCC", "SAMPNO"], how="left")
    merge_rows = ~btl_data["New Flag_y"].isnull() | ~btl_data["Comments_y"].isnull()
    btl_data.loc[merge_rows, "New Flag_x"] = btl_data.loc[merge_rows, "New Flag_y"]
    btl_data.loc[merge_rows, "Comments_x"] = btl_data.loc[merge_rows, "Comments_y"]
    btl_data = btl_data.rename(
        columns={"New Flag_x": "New Flag", "Comments_x": "Comments"}
    ).drop(columns=["New Flag_y", "Comments_y"])

# make downcast data point by interpolating bottle points on CTD data
downcast_data = []
for ssscc in ssscc_list:
    btl_rows = btl_data["SSSCC"] == ssscc
    ctd_rows = ctd_data["SSSCC"] == ssscc
    df = pd.DataFrame()
    for v in ["CTDTMP", "CTDSAL", "CTDOXY"]:
        df["CTDPRS"] = btl_data.loc[btl_rows, "CTDPRS"]
        df[v] = np.interp(
            btl_data.loc[btl_rows, "CTDPRS"],
            ctd_data.loc[ctd_rows, "CTDPRS"],
            ctd_data.loc[ctd_rows, v],
        )
    downcast_data.append(df)

downcast_data = pd.concat(downcast_data)

# intialize widgets
parameter = Select(
    title="Parameter", options=["CTDSAL", "CTDTMP", "CTDOXY"], value="CTDSAL"
)
ref_dict = {"CTDSAL": "SALNTY", "CTDTMP": "REFTMP", "CTDOXY": "OXYGEN"}
ref_param = Select(
    title="Reference",
    options=["SALNTY", "REFTMP", "OXYGEN"],
    value=ref_dict[parameter.value],
    disabled=True,
)

station = Select(title="Station", options=ssscc_list, value=ssscc_list[0])

# explanation of flags:
# https://cchdo.github.io/hdo-assets/documentation/manuals/pdf/90_1/chap4.pdf
flag_list = MultiSelect(
    title="Plot data flagged as:",
    value=["1", "2", "3", "9"],
    options=[
        ("1", "1 [Uncalibrated]"),
        ("2", "2 [Acceptable]"),
        ("3", "3 [Questionable]"),
        ("4", "4 [Bad]"),
        ("9", "9 [No Data]"),
    ],
    min_height=98,
)
# returns list of select options, e.g., ['2'] or ['1','2']
flag_input = Select(
    title="Flag:",
    options=[
        ("1", "1 [Uncalibrated]"),
        ("2", "2 [Acceptable]"),
        ("3", "3 [Questionable]"),
        ("4", "4 [Bad]"),
    ],
    value="3",
)
comment_box = TextInput(value="", title="Comment:")

# button_type: default, primary, success, warning or danger
flag_button = Button(label="Apply to selected", button_type="primary")
comment_button = Button(label="Apply to selected", button_type="warning")
save_button = Button(label="Save flagged data", button_type="success")
exit_button = Button(label="Exit flagging tool", button_type="danger")

vspace = Div(text=""" """, width=200, height=50)
residual_guide_text = Div(
    text="""<br><br>
    <u>Questionable Limits</u><br>
    <b>0.020</b> [0-500 m]<br>
    <b>0.010</b> [500-1000 m]<br>
    <b>0.005</b> [1000-2000 m]<br>
    <b>0.002</b> [2000-6000 m]""",
    width=150,
    height=135,
)
bulk_flag_text = Div(
    text="""<br><br>
    <b>Bulk Bottle Flagging:</b><br>
    Select multiple bottles using the table
    (with shift/control) or the 'Box Select' tool on the plot.""",
    width=150,
    height=135,
)

# set up datasources
src_table = ColumnDataSource(data=dict())
src_table_changes = ColumnDataSource(data=dict())
src_plot_trace = ColumnDataSource(data=dict(x=[], y=[]))
src_plot_downcast = ColumnDataSource(data=dict(x=[], y=[]))
src_plot_upcast = ColumnDataSource(data=dict(x=[], y=[]))
src_plot_btl = ColumnDataSource(data=dict(x=[], y=[]))

# set up plots
fig = figure(
    height=900,
    width=400,
    title="{} vs CTDPRS [Station {}]".format(parameter.value, station.value),
    tools="pan,box_zoom,wheel_zoom,box_select,reset",
    y_axis_label="Pressure (dbar)",
)
fig.line(
    "x",
    "y",
    line_color="#000000",
    line_width=2,
    source=src_plot_trace,
    legend_label="CTD Trace",
)
btl_sal = fig.asterisk(
    "x",
    "y",
    size=12,
    line_width=1.5,
    color="#0033CC",
    source=src_plot_btl,
    legend_label="Bottle sample",
)
ctd_sal = fig.circle(
    "x",
    "y",
    size=7,
    color="#BB0000",
    source=src_plot_downcast,
    legend_label="Downcast CTD sample",
)
upcast_sal = fig.triangle(
    "x",
    "y",
    size=7,
    color="#00BB00",
    source=src_plot_upcast,
    legend_label="Upcast CTD sample",
)
fig.select(BoxSelectTool).continuous = False
# fig.y_range = Range1d(max_y_value, 0)
fig.y_range.flipped = True  # invert y-axis
fig.legend.location = "bottom_right"
fig.legend.border_line_width = 3
fig.legend.border_line_alpha = 1
btl_sal.nonselection_glyph.line_alpha = 0.2
ctd_sal.nonselection_glyph.fill_alpha = 1  # makes CTDSAL *not* change on select
upcast_sal.nonselection_glyph.fill_alpha = 1  # makes CTDSAL *not* change on select

threshes = {
    "CTDSAL": np.array([0.002, 0.005, 0.010, 0.020]),
    "CTDTMP": np.array([0.002, 0.005, 0.010, 0.020]),
    "CTDOXY": np.array([0.625, 1.250, 2.500, 5.000]),   #   Not advised to follow this - biology happens
}

#   Residuals plot
src_plot_btl_del = ColumnDataSource(data=dict(x=[], y=[]))
fig2 = figure(
    height=900,
    width=400,
    title="{} residual vs CTDPRS [Station {}]".format(parameter.value, station.value),
    tools="pan,box_zoom,wheel_zoom,box_select,reset",
    # y_axis_label="Pressure (dbar)",
    y_range=fig.y_range,
)
# thresh = np.array([0.002, 0.005, 0.010, 0.020])
thresh = threshes[parameter.value]
p_range = np.array([6000, 2000, 1000, 500])
thresh = np.append(thresh, thresh[-1])
p_range = np.append(p_range, 0)
btl_sal2 = fig2.asterisk(
    "x",
    "y",
    size=12,
    line_width=1.5,
    color="#0033CC",
    source=src_plot_btl_del,
)
fig2.step(thresh, p_range)
fig2.step(-thresh, p_range)
fig2.select(BoxSelectTool).continuous = False
fig2.y_range = fig.y_range
# fig2.y_range.flipped = True  # invert y-axis

# define callback functions


def update_selectors():

    print("exec update_selectors()")
    ctd_rows = ctd_data["SSSCC"] == station.value
    table_rows = btl_data["SSSCC"] == station.value
    btl_rows = (btl_data["New Flag"].isin([int(v) for v in flag_list.value])) & (
        btl_data["SSSCC"] == station.value
    )

    # update table data
    current_table = btl_data[table_rows].reset_index()

    # print("Parameter: ", parameter.value, current_table[parameter.value][0])
    ref_param.value = ref_dict[parameter.value]
    # print("Reference parameter: ", ref_param.value, current_table[ref_param.value][0])

    src_table.data = {  # this causes edit_flag() to execute
        "SSSCC": current_table["SSSCC"],
        "SAMPNO": current_table["SAMPNO"],
        "CTDPRS": current_table["CTDPRS"],
        "t_res": current_table["t_res"],
        "o_res": current_table["o_res"],
        "CTD Param": current_table[parameter.value].round(4),
        "Reference": current_table[ref_param.value].round(4),
        "diff": current_table["Residual"],
        "flag": current_table["New Flag"],
        "Comments": current_table["Comments"],
    }

    # update plot data
    src_plot_trace.data = {
        "x": ctd_data.loc[ctd_rows, parameter.value],
        "y": ctd_data.loc[ctd_rows, "CTDPRS"],
    }
    src_plot_upcast.data = {
        "x": btl_data.loc[table_rows, parameter.value],
        "y": btl_data.loc[table_rows, "CTDPRS"],
    }
    src_plot_downcast.data = {
        "x": downcast_data.loc[table_rows, parameter.value],
        "y": downcast_data.loc[table_rows, "CTDPRS"],
    }
    src_plot_btl.data = {
        "x": btl_data.loc[btl_rows, ref_param.value],
        "y": btl_data.loc[btl_rows, "CTDPRS"],
    }
    src_plot_btl_del.data = {
        "x": btl_data.loc[btl_rows, deltas_d[parameter.value]],
        "y": btl_data.loc[btl_rows, "CTDPRS"],
    }

    # update plot labels/axlims
    fig.title.text = "{} vs CTDPRS [Station {}]".format(parameter.value, station.value)
    fig.xaxis.axis_label = parameter.value

    fig2.title.text = "{} Residual".format(parameter.value)
    fig2.xaxis.axis_label = parameter.value

    # Set the y-range from 0 to the maximum value
    max_y = ctd_data.loc[ctd_rows, "CTDPRS"].max()
    fig.y_range = Range1d(max_y + 0.05 * max_y, 0)
    fig2.y_range = fig.y_range

    # Reset the x-axes when plots refresh
    max_x_fig = ctd_data.loc[ctd_rows, parameter.value].max()
    min_x_fig = ctd_data.loc[ctd_rows, parameter.value].min()
    fig.x_range = Range1d(min_x_fig, max_x_fig)

    max_x_fig2 = max(btl_data.loc[btl_rows, deltas_d[parameter.value]].max(), 0.025)
    min_x_fig2 = min(btl_data.loc[btl_rows, deltas_d[parameter.value]].min(), -0.025)
    fig2.x_range = Range1d(min_x_fig2, max_x_fig2)

    # deselect all datapoints
    btl_sal.data_source.selected.indices = []
    src_table.selected.indices = []


def edit_flag():

    print("exec edit_flag()")

    btl_data.loc[
        btl_data["SSSCC"] == src_table.data["SSSCC"].values[0],
        "New Flag",
    ] = src_table.data["flag"].values
    btl_data.loc[
        btl_data["SSSCC"] == src_table.data["SSSCC"].values[0],
        "Comments",
    ] = src_table.data["Comments"].values

    edited_rows = (btl_data["New Flag"].isin([3, 4])) | (btl_data["Comments"] != "")

    src_table_changes.data = {
        "SSSCC": btl_data.loc[edited_rows, "SSSCC"],
        "SAMPNO": btl_data.loc[edited_rows, "SAMPNO"],
        "diff": btl_data.loc[edited_rows, "Residual"],
        "flag_old": btl_data.loc[edited_rows, "SALNTY_FLAG_W"],
        "flag_new": btl_data.loc[edited_rows, "New Flag"],
        "Comments": btl_data.loc[edited_rows, "Comments"],
    }


def apply_flag():

    print("Applying flags")

    table_rows = btl_data["SSSCC"] == station.value
    selected_rows = src_table.selected.indices

    # update table data
    current_table = btl_data[table_rows].reset_index()
    current_table.loc[selected_rows, "New Flag"] = int(flag_input.value)
    src_table.data = {  # this causes edit_flag() to execute
        "SSSCC": current_table["SSSCC"],
        "SAMPNO": current_table["SAMPNO"],
        "CTDPRS": current_table["CTDPRS"],
        "CTDSAL": current_table["CTDSAL"],
        "SALNTY": current_table["SALNTY"],
        "diff": current_table["Residual"],
        "flag": current_table["New Flag"],
        "Comments": current_table["Comments"],
        "CTD Param": current_table[parameter.value].round(4),
        "Reference": current_table[ref_param.value].round(4),
        "o_res": current_table["o_res"],
        "t_res": current_table["t_res"],
    }


def apply_comment():

    print("Applying Comments")

    table_rows = btl_data["SSSCC"] == station.value
    selected_rows = src_table.selected.indices

    # update table data
    current_table = btl_data[table_rows].reset_index()
    current_table.loc[selected_rows, "Comments"] = comment_box.value
    src_table.data = {  # this causes edit_flag() to execute
        "SSSCC": current_table["SSSCC"],
        "SAMPNO": current_table["SAMPNO"],
        "CTDPRS": current_table["CTDPRS"],
        "CTDSAL": current_table["CTDSAL"],
        "SALNTY": current_table["SALNTY"],
        "diff": current_table["Residual"],
        "flag": current_table["New Flag"],
        "Comments": current_table["Comments"],
        "CTD Param": current_table[parameter.value].round(4),
        "Reference": current_table[ref_param.value].round(4),
        "o_res": current_table["o_res"],
        "t_res": current_table["t_res"],
    }


def save_data():

    print("Saving flagged data...")

    # get data from table
    df_out = pd.DataFrame.from_dict(src_table_changes.data)

    # minor changes to columns/names/etc.
    df_out = df_out.rename(
        columns={
            "flag_new": "value",
            "SSSCC": "cast_id",
            "SAMPNO": "bottle_num",
            "Comments": "notes",
        }
    ).drop(columns=["flag_old", "diff"])

    # save it
    salt = df_node_to_BottleFlags(df_out)
    flagfile = validate_file(FLAGFILE, create=True)
    save_node(flagfile, salt, "salt", create_new=True)


def exit_bokeh():
    print("Stopping flagging software")
    sys.exit()


def selected_from_plot(attr, old, new):

    # update using bottle number, not index
    # currently there is a bug if not all data from a cast are plotted
    src_table.selected.indices = new


def selected_from_table(attr, old, new):

    # update using bottle number, not index
    # currently there is a bug if not all data from a cast are plotted
    btl_sal.data_source.selected.indices = new


# set up change callbacks
parameter.on_change("value", lambda attr, old, new: update_selectors())
station.on_change("value", lambda attr, old, new: update_selectors())
flag_list.on_change("value", lambda attr, old, new: update_selectors())
flag_button.on_click(apply_flag)
comment_button.on_click(apply_comment)
save_button.on_click(save_data)
exit_button.on_click(exit_bokeh)
src_table.on_change("data", lambda attr, old, new: edit_flag())
src_table.selected.on_change("indices", selected_from_table)
btl_sal.data_source.selected.on_change("indices", selected_from_plot)


# build data tables
columns = []
fields = [
    "SSSCC",
    "SAMPNO",
    "CTDPRS",
    "CTD Param",
    "Reference",
    "t_res",
    "diff",
    "o_res",
    "flag",
    "Comments",
]
ref_dict
titles = [
    "SSSCC",
    "Bottle",
    "CTDPRS",
    "CTD Param",
    "Reference",
    "t_res",
    "s_res",
    "o_res",
    "Flag",
    "Comments",
]
widths = [50, 40, 65, 65, 65, 65, 65, 65, 15, 200]
for field, title, width in zip(fields, titles, widths):
    if field == "flag":
        strfmt_in = {"text_align": "center", "font_style": "bold"}
    elif field == "Comments":
        strfmt_in = {}
    else:
        strfmt_in = {"text_align": "right"}
    columns.append(
        TableColumn(
            field=field,
            title=title,
            width=width,
            formatter=StringFormatter(**strfmt_in),
        )
    )

columns_changed = []
fields = ["SSSCC", "SAMPNO", "diff", "flag_old", "flag_new", "Comments"]
titles = ["SSSCC", "Bottle", "Residual", "Old", "New", "Comments"]
widths = [50, 40, 65, 15, 15, 375]
for field, title, width in zip(fields, titles, widths):
    if field == "flag_old":
        strfmt_in = {"text_align": "center", "font_style": "bold"}
    elif field == "flag_new":
        strfmt_in = {"text_align": "center", "font_style": "bold", "text_color": "red"}
    elif field == "Comments":
        strfmt_in = {}
    else:
        strfmt_in = {"text_align": "right"}
    columns_changed.append(
        TableColumn(
            field=field,
            title=title,
            width=width,
            formatter=StringFormatter(**strfmt_in),
        )
    )

data_table = DataTable(
    source=src_table,
    columns=columns,
    index_width=20,
    #   width is originally 565 + 50 + 20 = 635
    width=565 + 150 + 20,  # sum of col widths + fudge factor + idx width
    height=600,
    editable=True,
    fit_columns=False,
    sortable=False,
)
data_table_changed = DataTable(
    source=src_table_changes,
    columns=columns_changed,
    index_width=20,
    width=565 + 50 + 20,  # sum of col widths + fudge factor + idx width
    height=200,
    editable=False,
    fit_columns=False,
    sortable=False,
)
data_table_title = Div(text="""<b>All Station Data:</b>""", width=200, height=15)
data_table_changed_title = Div(
    text="""<b>Flagged Data:</b>""", width=200, height=15
)

controls = column(
    parameter,
    ref_param,
    station,
    flag_list,
    residual_guide_text,
    bulk_flag_text,
    flag_input,
    flag_button,
    comment_box,
    comment_button,
    vspace,
    save_button,
    exit_button,
    width=170,
)
tables = column(
    data_table_title, data_table, data_table_changed_title, data_table_changed
)

curdoc().add_root(row(controls, tables, fig, fig2))
curdoc().title = "CTDO Data Flagging Tool"

update_selectors()
