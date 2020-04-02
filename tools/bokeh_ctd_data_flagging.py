import pandas as pd
import numpy as np
import glob

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.plotting import figure
from bokeh.models import (
    Button,
    ColumnDataSource,
    # CustomJS,
    DataTable,
    TableColumn,
    Select,
    MultiSelect,
    StringFormatter,  # https://docs.bokeh.org/en/latest/_modules/bokeh/models/widgets/tables.html#DataTable
    Div,
    TextInput,
    BoxSelectTool,
)

# TODO: load in existing handcoded salts if it exists

# load continuous CTD data and make into a dict (only ~20MB)
file_list = sorted(glob.glob("../data/pressure/*.csv"))
ssscc_list = [ssscc.strip("../data/pressure/")[:5] for ssscc in file_list]
ctd_data = []
for f in file_list:
    df = pd.read_csv(f, header=12, skiprows=[13], skipfooter=1, engine="python")
    df["SSSCC"] = f.strip("../data/pressure/")[:5]
    ctd_data.append(df)
ctd_data = pd.concat(ctd_data, axis=0, sort=False)


# adapted from compare_salinities.ipynb
# load salt file
file_list = sorted(glob.glob("../data/salt/*.csv"))
ssscc_list = [ssscc.strip("../data/salt/")[:5] for ssscc in file_list]
salt_data = []
for f in file_list:
    df = pd.read_csv(f, usecols=["STNNBR", "CASTNO", "SAMPNO", "SALNTY"])
    df["SSSCC"] = f.strip("../data/salt/")[:5]
    salt_data.append(df)
salt_data = pd.concat(salt_data, axis=0, sort=False)
salt_data["SALNTY"] = salt_data["SALNTY"].round(4)
if "SALNTY_FLAG_W" not in salt_data.columns:
    salt_data["SALNTY_FLAG_W"] = 2

# load ctd btl data
df_ctd_btl = pd.read_csv(
    "../data/scratch_folder/ctd_to_bottle.csv",
    skiprows=[1],
    skipfooter=1,
    engine="python",
)
df_btl_all = pd.merge(df_ctd_btl, salt_data, on=["STNNBR", "CASTNO", "SAMPNO"])
btl_data = df_btl_all.loc[
    :,
    [
        "SSSCC",
        "SAMPNO",
        "CTDPRS",
        "CTDTMP",
        "REFTMP",
        "CTDSAL",
        "SALNTY",
        "SALNTY_FLAG_W",
    ],
]
btl_data["Residual"] = btl_data["SALNTY"] - btl_data["CTDSAL"]
btl_data[["CTDPRS", "Residual"]] = btl_data[["CTDPRS", "Residual"]].round(4)
btl_data["Comments"] = ""
btl_data["New Flag"] = btl_data["SALNTY_FLAG_W"].copy()


# intialize widgets
save_button = Button(label="Save flagged data", button_type="success")
parameter = Select(title="Parameter", options=["CTDSAL", "CTDTMP"], value="CTDSAL")
ref_param = Select(title="Reference", options=["SALNTY"], value="SALNTY")
# ref_param.options = ["foo","bar"]  # can dynamically change dropdowns
station = Select(title="Station", options=ssscc_list, value="00101")
# explanation of flags:
# https://cchdo.github.io/hdo-assets/documentation/manuals/pdf/90_1/chap4.pdf
flag_list = MultiSelect(
    title="Plot data flagged as:",
    value=["1", "2", "3"],
    options=[
        ("1", "1 [Uncalibrated]"),
        ("2", "2 [Acceptable]"),
        ("3", "3 [Questionable]"),
        ("4", "4 [Bad]"),
    ],
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

vspace = Div(text=""" """, width=200, height=95)

# set up datasources
src_table = ColumnDataSource(data=dict())
src_table_changes = ColumnDataSource(data=dict())
src_plot_trace = ColumnDataSource(data=dict(x=[], y=[]))
src_plot_ctd = ColumnDataSource(data=dict(x=[], y=[]))
src_plot_btl = ColumnDataSource(data=dict(x=[], y=[]))

# set up plots
fig = figure(
    plot_height=800,
    plot_width=400,
    title="{} vs CTDPRS [Station {}]".format(parameter.value, station.value),
    tools="pan,box_zoom,wheel_zoom,box_select,reset",
    y_axis_label="Pressure (dbar)",
)
fig.select(BoxSelectTool).select_every_mousemove = False
fig.y_range.flipped = True  # invert y-axis
fig.line(
    "x", "y", line_color="#000000", line_width=2, source=src_plot_trace,
)
btl_sal = fig.asterisk(
    "x", "y", size=10, line_width=1.5, color="#0033CC", source=src_plot_btl
)
ctd_sal = fig.circle("x", "y", size=7, color="#BB0000", source=src_plot_ctd)
btl_sal.nonselection_glyph.line_alpha = 0.2
ctd_sal.nonselection_glyph.fill_alpha = 1  # makes CTDSAL change on select

parameter.on_change("value", lambda attr, old, new: update_selectors())
station.on_change("value", lambda attr, old, new: update_selectors())
flag_list.on_change("value", lambda attr, old, new: update_selectors())


def update_selectors():

    print("exec update_selectors()")
    ctd_rows = ctd_data["SSSCC"] == station.value
    table_rows = btl_data["SSSCC"] == station.value
    btl_rows = (btl_data["New Flag"].isin(flag_list.value)) & (
        btl_data["SSSCC"] == station.value
    )

    # update table data
    current_table = btl_data[table_rows].reset_index()
    src_table.data = {  # this causes update_flag() to execute
        "SSSCC": current_table["SSSCC"],
        "SAMPNO": current_table["SAMPNO"],
        "CTDPRS": current_table["CTDPRS"],
        "CTDSAL": current_table["CTDSAL"],
        "SALNTY": current_table["SALNTY"],
        "diff": current_table["Residual"],
        "flag": current_table["New Flag"],
        "comments": current_table["Comments"],
    }

    # update plot data
    src_plot_trace.data = {
        "x": ctd_data.loc[ctd_rows, parameter.value],
        "y": ctd_data.loc[ctd_rows, "CTDPRS"],
    }
    src_plot_ctd.data = {
        "x": btl_data.loc[table_rows, parameter.value],
        "y": btl_data.loc[table_rows, "CTDPRS"],
    }
    src_plot_btl.data = {
        "x": btl_data.loc[btl_rows, "SALNTY"],
        "y": btl_data.loc[btl_rows, "CTDPRS"],
    }

    # update plot labels/axlims
    fig.title.text = "{} vs CTDPRS [Station {}]".format(parameter.value, station.value)
    fig.xaxis.axis_label = parameter.value

    # deselect all datapoints
    # breakpoint()
    btl_sal.data_source.selected.indices = []
    src_table.selected.indices = []


def update_flag():

    print("exec update_flag()")

    # breakpoint()

    btl_data.loc[
        btl_data["SSSCC"] == src_table.data["SSSCC"].values[0], "New Flag",
    ] = src_table.data["flag"].values
    btl_data.loc[
        btl_data["SSSCC"] == src_table.data["SSSCC"].values[0], "Comments",
    ] = src_table.data["comments"].values

    edited_rows = (btl_data["New Flag"].isin([3, 4])) | (btl_data["Comments"] != "")

    src_table_changes.data = {
        "SSSCC": btl_data.loc[edited_rows, "SSSCC"],
        "SAMPNO": btl_data.loc[edited_rows, "SAMPNO"],
        "diff": btl_data.loc[edited_rows, "Residual"],
        "flag_old": btl_data.loc[edited_rows, "SALNTY_FLAG_W"],
        "flag_new": btl_data.loc[edited_rows, "New Flag"],
        "comments": btl_data.loc[edited_rows, "Comments"],
    }


def apply_flag():

    print("Applying flags")

    table_rows = btl_data["SSSCC"] == station.value
    selected_rows = src_table.selected.indices

    # update table data
    current_table = btl_data[table_rows].reset_index()
    current_table.loc[selected_rows, "New Flag"] = int(flag_input.value)
    src_table.data = {  # this causes update_flag() to execute
        "SSSCC": current_table["SSSCC"],
        "SAMPNO": current_table["SAMPNO"],
        "CTDPRS": current_table["CTDPRS"],
        "CTDSAL": current_table["CTDSAL"],
        "SALNTY": current_table["SALNTY"],
        "diff": current_table["Residual"],
        "flag": current_table["New Flag"],
        "comments": current_table["Comments"],
    }


flag_button.on_click(apply_flag)


def apply_comment():

    print("Applying comments")

    table_rows = btl_data["SSSCC"] == station.value
    selected_rows = src_table.selected.indices

    # update table data
    current_table = btl_data[table_rows].reset_index()
    current_table.loc[selected_rows, "Comments"] = comment_box.value
    src_table.data = {  # this causes update_flag() to execute
        "SSSCC": current_table["SSSCC"],
        "SAMPNO": current_table["SAMPNO"],
        "CTDPRS": current_table["CTDPRS"],
        "CTDSAL": current_table["CTDSAL"],
        "SALNTY": current_table["SALNTY"],
        "diff": current_table["Residual"],
        "flag": current_table["New Flag"],
        "comments": current_table["Comments"],
    }


comment_button.on_click(apply_comment)


def save_data():

    print("Saving flagged data...")

    # get data from table
    df_out = pd.DataFrame.from_dict(src_table_changes.data)

    # minor changes to columns/names/etc.
    df_out = df_out.rename(columns={"flag_new": "salinity_flag"}).drop(
        columns="flag_old"
    )

    # save it
    df_out.to_csv("salt_flags_handcoded.csv", index=None)


save_button.on_click(save_data)


columns = [
    TableColumn(
        field="SSSCC",
        title="SSSCC",
        width=40,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="SAMPNO",
        title="Bottle",
        width=20,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="CTDPRS",
        title="CTDPRS",
        width=75,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="CTDSAL",
        title="CTDSAL",
        width=65,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="SALNTY",
        title="SALNTY",
        width=65,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="diff",
        title="Residual",
        width=65,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="flag",
        title="Flag",
        width=15,
        formatter=StringFormatter(text_align="center", font_style="bold"),
    ),
    TableColumn(
        field="comments",
        title="Comments",
        width=135,
        formatter=StringFormatter(text_align="left"),
    ),
]
columns_changed = [
    TableColumn(
        field="SSSCC",
        title="SSSCC",
        width=40,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="SAMPNO",
        title="Bottle",
        width=20,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="diff",
        title="Residual",
        width=40,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="flag_old",
        title="Old",
        width=20,
        formatter=StringFormatter(text_align="center", font_style="bold"),
    ),
    TableColumn(
        field="flag_new",
        title="New",
        width=20,
        formatter=StringFormatter(
            text_align="center", font_style="bold", text_color="red"
        ),
    ),
    TableColumn(
        field="comments",
        title="Comments",
        width=200,
        formatter=StringFormatter(text_align="left"),
    ),
]
data_table = DataTable(
    source=src_table,
    columns=columns,
    index_width=20,
    width=480 + 20,  # sum of col widths + idx width
    height=600,
    editable=True,
    fit_columns=True,
    sortable=False,
)
data_table_changed = DataTable(
    source=src_table_changes,
    columns=columns_changed,
    index_width=20,
    width=480 + 20,  # sum of col widths + idx width
    height=200,
    editable=False,
    fit_columns=True,
    sortable=False,
)
data_table_title = Div(text="""<b>All Station Data:</b>""", width=200, height=15)
data_table_changed_title = Div(text="""<b>Edited Data:</b>""", width=200, height=15)

# run update() when user selects new column (may indicate new flag value)
src_table.on_change("data", lambda attr, old, new: update_flag())


def selected_from_plot(attr, old, new):

    src_table.selected.indices = new


def selected_from_table(attr, old, new):

    btl_sal.data_source.selected.indices = new


src_table.selected.on_change("indices", selected_from_table)
btl_sal.data_source.selected.on_change("indices", selected_from_plot)

controls = column(
    parameter,
    ref_param,
    station,
    flag_list,
    vspace,
    flag_input,
    flag_button,
    comment_box,
    comment_button,
    vspace,
    save_button,
    width=170,
)
tables = column(
    data_table_title, data_table, data_table_changed_title, data_table_changed
)

curdoc().add_root(row(controls, tables, fig))
curdoc().title = "CTDO Data Flagging Tool"

update_selectors()
