import pandas as pd
import glob
import pickle
import gsw

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.plotting import figure
from bokeh.models import (
    Button,
    ColumnDataSource,
    DataTable,
    TableColumn,
    Select,
    MultiSelect,
    StringFormatter,  # https://docs.bokeh.org/en/latest/_modules/bokeh/models/widgets/tables.html#DataTable
    Div,
    TextInput,
    BoxSelectTool,
)

# TODO: abstract parts of this to a separate file
# TODO: following above, make parts reusable?

# load continuous CTD data and make into a dict (only ~20MB)
file_list = sorted(glob.glob("../../data/pressure/*ct1.csv"))
ssscc_list = [ssscc.strip("../../data/pressure/")[:5] for ssscc in file_list]
ctd_data = []
for f in file_list:
    print(f"Loading {f}")
    df = pd.read_csv(f, header=12, skiprows=[13], skipfooter=1, engine="python")
    df["SSSCC"] = f.strip("../../data/pressure/")[:5]
    ctd_data.append(df)
ctd_data = pd.concat(ctd_data, axis=0, sort=False)

# load bottle trip file
file_list = sorted(glob.glob("../../data/bottle/*.pkl"))
ssscc_list = [ssscc.strip("../../data/bottle/")[:5] for ssscc in file_list]
upcast_data = []
for f in file_list:
    with open(f, "rb") as x:
        df = pickle.load(x)
        df["SSSCC"] = f.strip("../../data/bottle/")[:5]
        # change to secondary if that is what's used
        upcast_data.append(df[["SSSCC", "CTDCOND1", "CTDTMP1", "CTDPRS"]])
upcast_data = pd.concat(upcast_data, axis=0, sort=False)
upcast_data["CTDSAL"] = gsw.SP_from_C(
    upcast_data["CTDCOND1"], upcast_data["CTDTMP1"], upcast_data["CTDPRS"]
)

# load salt file (adapted from compare_salinities.ipynb)
file_list = sorted(glob.glob("../../data/salt/*.csv"))
ssscc_list = [ssscc.strip("../../data/salt/")[:5] for ssscc in file_list]
salt_data = []
for f in file_list:
    df = pd.read_csv(f, usecols=["STNNBR", "CASTNO", "SAMPNO", "SALNTY"])
    df["SSSCC"] = f.strip("../../data/salt/")[:5]
    salt_data.append(df)
salt_data = pd.concat(salt_data, axis=0, sort=False)
salt_data["SALNTY"] = salt_data["SALNTY"].round(4)
if "SALNTY_FLAG_W" not in salt_data.columns:
    salt_data["SALNTY_FLAG_W"] = 2

# load ctd btl data
df_ctd_btl = pd.read_csv(
    "../../data/scratch_folder/ctd_to_bottle.csv",
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

# update with old handcoded flags if file exists
handcoded_file = "salt_flags_handcoded.csv"
if glob.glob(handcoded_file):
    handcodes = pd.read_csv(handcoded_file, dtype={"SSSCC": str}, keep_default_na=False)
    handcodes = handcodes.rename(columns={"salinity_flag": "New Flag"}).drop(
        columns="diff"
    )
    # there's gotta be a better way... but this is good enough for now
    btl_data = btl_data.merge(handcodes, on=["SSSCC", "SAMPNO"], how="left")
    merge_rows = ~btl_data["New Flag_y"].isnull() | ~btl_data["Comments_y"].isnull()
    btl_data.loc[merge_rows, "New Flag_x"] = btl_data.loc[merge_rows, "New Flag_y"]
    btl_data.loc[merge_rows, "Comments_x"] = btl_data.loc[merge_rows, "Comments_y"]
    btl_data = btl_data.rename(
        columns={"New Flag_x": "New Flag", "Comments_x": "Comments"}
    ).drop(columns=["New Flag_y", "Comments_y"])

# intialize widgets
save_button = Button(label="Save flagged data", button_type="success")
parameter = Select(title="Parameter", options=["CTDSAL", "CTDTMP"], value="CTDSAL")
ref_param = Select(title="Reference", options=["SALNTY"], value="SALNTY")
# ref_param.options = ["foo","bar"]  # can dynamically change dropdowns
station = Select(title="Station", options=ssscc_list, value=ssscc_list[0])
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

vspace = Div(text=""" """, width=200, height=65)
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
src_plot_ctd = ColumnDataSource(data=dict(x=[], y=[]))
src_plot_upcast = ColumnDataSource(data=dict(x=[], y=[]))
src_plot_btl = ColumnDataSource(data=dict(x=[], y=[]))

# set up plots
fig = figure(
    plot_height=800,
    plot_width=400,
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
    source=src_plot_ctd,
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
fig.select(BoxSelectTool).select_every_mousemove = False
fig.y_range.flipped = True  # invert y-axis
fig.legend.location = "bottom_right"
fig.legend.border_line_width = 3
fig.legend.border_line_alpha = 1
btl_sal.nonselection_glyph.line_alpha = 0.2
ctd_sal.nonselection_glyph.fill_alpha = 1  # makes CTDSAL *not* change on select
upcast_sal.nonselection_glyph.fill_alpha = 1  # makes CTDSAL *not* change on select

# define callback functions


def update_selectors():

    print("exec update_selectors()")
    ctd_rows = ctd_data["SSSCC"] == station.value
    table_rows = btl_data["SSSCC"] == station.value
    btl_rows = (btl_data["New Flag"].isin([int(v) for v in flag_list.value])) & (
        btl_data["SSSCC"] == station.value
    )
    # breakpoint()

    # update table data
    current_table = btl_data[table_rows].reset_index()
    src_table.data = {  # this causes edit_flag() to execute
        "SSSCC": current_table["SSSCC"],
        "SAMPNO": current_table["SAMPNO"],
        "CTDPRS": current_table["CTDPRS"],
        "CTDSAL": current_table["CTDSAL"],
        "SALNTY": current_table["SALNTY"],
        "diff": current_table["Residual"],
        "flag": current_table["New Flag"],
        "Comments": current_table["Comments"],
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
    src_plot_upcast.data = {
        "x": upcast_data.loc[upcast_data["SSSCC"] == station.value, "CTDSAL"],
        "y": upcast_data.loc[upcast_data["SSSCC"] == station.value, "CTDPRS"],
    }
    src_plot_btl.data = {
        "x": btl_data.loc[btl_rows, "SALNTY"],
        "y": btl_data.loc[btl_rows, "CTDPRS"],
    }

    # update plot labels/axlims
    fig.title.text = "{} vs CTDPRS [Station {}]".format(parameter.value, station.value)
    fig.xaxis.axis_label = parameter.value

    # deselect all datapoints
    btl_sal.data_source.selected.indices = []
    src_table.selected.indices = []


def edit_flag():

    print("exec edit_flag()")

    btl_data.loc[
        btl_data["SSSCC"] == src_table.data["SSSCC"].values[0], "New Flag",
    ] = src_table.data["flag"].values
    btl_data.loc[
        btl_data["SSSCC"] == src_table.data["SSSCC"].values[0], "Comments",
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
    }


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


def selected_from_plot(attr, old, new):

    src_table.selected.indices = new


def selected_from_table(attr, old, new):

    btl_sal.data_source.selected.indices = new


# set up change callbacks
parameter.on_change("value", lambda attr, old, new: update_selectors())
station.on_change("value", lambda attr, old, new: update_selectors())
flag_list.on_change("value", lambda attr, old, new: update_selectors())
flag_button.on_click(apply_flag)
comment_button.on_click(apply_comment)
save_button.on_click(save_data)
src_table.on_change("data", lambda attr, old, new: edit_flag())
src_table.selected.on_change("indices", selected_from_table)
btl_sal.data_source.selected.on_change("indices", selected_from_plot)


# build data tables
columns = []
fields = ["SSSCC", "SAMPNO", "CTDPRS", "CTDSAL", "SALNTY", "diff", "flag", "Comments"]
titles = ["SSSCC", "Bottle", "CTDPRS", "CTDSAL", "SALNTY", "Residual", "Flag", "Comments"]
widths = [40, 20, 75, 65, 65, 65, 15, 135]
for (field, title, width) in zip(fields, titles, widths):
    if field == "flag":
        strfmt_in = {"text_align": "center", "font_style": "bold"}
    elif field == "Comments":
        strfmt_in = {}
    else:
        strfmt_in = {"text_align": "right"}
    columns.append(TableColumn(
        field=field,
        title=title,
        width=width,
        formatter=StringFormatter(**strfmt_in)
    ))

columns_changed = []
fields = ["SSSCC", "SAMPNO", "diff", "flag_old", "flag_new", "Comments"]
titles = ["SSSCC", "Bottle", "Residual", "Old", "New", "Comments"]
widths = [40, 20, 40, 20, 20, 200]
for (field, title, width) in zip(fields, titles, widths):
    if field == "flag_old":
        strfmt_in = {"text_align": "center", "font_style": "bold"}
    elif field == "flag_new":
        strfmt_in = {"text_align": "center", "font_style": "bold", "text_color": "red"}
    elif field == "Comments":
        strfmt_in = {}
    else:
        strfmt_in = {"text_align": "right"}
    columns_changed.append(TableColumn(
        field=field,
        title=title,
        width=width,
        formatter=StringFormatter(**strfmt_in)
    ))

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
data_table_changed_title = Div(text="""<b>Flagged Data:</b>""", width=200, height=15)

controls = column(
    parameter,
    ref_param,
    station,
    flag_list,
    bulk_flag_text,
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
