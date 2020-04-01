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
)

# load continuous CTD data and make into a dict (only ~20MB)
file_list = sorted(glob.glob("../data/pressure/*.csv"))
ssscc_list = [ssscc.strip("../data/pressure/")[:5] for ssscc in file_list]
ctd_data = []
for f in file_list:
    df = pd.read_csv(f, header=12, skiprows=[13], skipfooter=1, engine="python")
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
if "SALNTY_FLAG_W" not in salt_data.columns:
    salt_data["SALNTY_FLAG_W"] = 2

breakpoint()

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
        "STNNBR",
        "CASTNO",
        "SAMPNO",
        "CTDPRS",
        "CTDTMP",
        "REFTMP",
        "CTDSAL",
        "SALNTY",
        "SALNTY_FLAG_W",
    ],
]

### TODO: pick up from here


# intialize widgets
button = Button(label="Save flagged data", button_type="success")
parameter = Select(
    title="Parameter",
    options=["CTDTMP", "CTDSAL", "CTDOXY", "CTDRINKO"],
    value="CTDSAL",
)
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

source_table = ColumnDataSource(data=dict())
source_table_changes = ColumnDataSource(data=dict())
source_plot_ssscc = ColumnDataSource(data=dict(x=[], y=[]))
source_plot_all = ColumnDataSource(data=dict(x=[], y=[]))

### set up plots
plot_ssscc = figure(
    plot_height=800,
    plot_width=400,
    title="{} vs CTDPRS [Station {}]".format(parameter.value, station.value),
    tools="crosshair,pan,reset,wheel_zoom",
    y_axis_label="Pressure (dbar)",
)
plot_all = figure(
    plot_height=800,
    plot_width=400,
    title="{} vs CTDPRS [Station {}]".format(parameter.value, station.value),
    tools="crosshair,reset,wheel_zoom",
)
plot_ssscc.y_range.flipped = True  # invert y-axis
plot_all.y_range.flipped = True  # invert y-axis
plot_ssscc.line(
    "x", "y", line_color="#000000", line_width=2, source=source_plot_ssscc,
)
plot_all.line(
    "x", "y", line_color="#000000", line_width=2, source=source_plot_all,
)

parameter.on_change("value", lambda attr, old, new: update_selectors())
station.on_change("value", lambda attr, old, new: update_selectors())
flag_list.on_change("value", lambda attr, old, new: update_selectors())


def update_selectors():

    print("exec update_selectors()")

    df_edited = ctd_dict[station.value]

    # breakpoint()
    flag_col = parameter.value + "_FLAG_W"
    rows = df_edited[flag_col].isin(flag_list.value)
    source_plot_ssscc.data = {
        "x": df_edited.loc[rows, parameter.value],
        "y": df_edited.loc[rows, "CTDPRS"],
    }
    source_plot_all.data = {
        "x": df_edited.loc[rows, parameter.value],
        "y": df_edited.loc[rows, "CTDPRS"],
    }
    plot_ssscc.title.text = "{} vs CTDPRS [Station {}]".format(
        parameter.value, station.value
    )
    plot_all.title.text = "{} vs CTDPRS [Station {}]".format(
        parameter.value, station.value
    )
    plot_ssscc.x_range = plot_all.x_range
    plot_ssscc.y_range = plot_all.y_range
    plot_ssscc.xaxis.axis_label = parameter.value
    plot_all.xaxis.axis_label = parameter.value


def update_flag():

    print("exec update_flag()")

    # breakpoint()

    df_edited.loc[
        df_edited["SSSCC_rinko"] == source_table.data["SSSCC"].values[0], "Flag New",
    ] = source_table.data["flag"].values
    df_edited.loc[
        df_edited["SSSCC_rinko"] == source_table.data["SSSCC"].values[0], "Comments",
    ] = source_table.data["comments"].values

    edited_rows = (df_edited["Flag"] - df_edited["Flag New"]) != 0

    source_table_changes.data = {
        "SSSCC": df_edited.loc[edited_rows, "SSSCC_rinko"],
        "CTDPRS": df_edited.loc[edited_rows, "CTDPRS_rinko_ctd"],
        "diff": df_edited.loc[edited_rows, "Residual"],
        "flag_old": df_edited.loc[edited_rows, "Flag"],
        "flag_new": df_edited.loc[edited_rows, "Flag New"],
        "comments": df_edited.loc[edited_rows, "Comments"],
    }


def save_data():

    print("Save button clicked...")

    # df_out =


button.on_click(save_data)


columns = [
    TableColumn(
        field="SSSCC",
        title="SSSCC",
        width=40,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="CTDPRS",
        title="CTDPRS",
        width=80,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="REFOXY",
        title="REFOXY",
        width=80,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="CTDRINKO",
        title="CTDRINKO",
        width=80,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="diff",
        title="Residual",
        width=80,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="flag",
        title="Flag",
        width=20,
        formatter=StringFormatter(text_align="center", font_style="bold"),
    ),
    TableColumn(
        field="comments",
        title="Comments",
        width=100,
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
        field="CTDPRS",
        title="CTDPRS",
        width=80,
        formatter=StringFormatter(text_align="right"),
    ),
    TableColumn(
        field="diff",
        title="Residual",
        width=80,
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
        width=100,
        formatter=StringFormatter(text_align="left"),
    ),
]
data_table = DataTable(
    source=source_table,
    columns=columns,
    index_width=20,
    width=480 + 20,  # sum of col widths + idx width
    height=600,
    editable=True,
    fit_columns=True,
    sortable=False,
)
data_table_changed = DataTable(
    source=source_table_changes,
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
# source_table.selected.on_change("indices", lambda attr, old, new: update_flag())
# source_table.on_change("data", lambda attr, old, new: update_flag())

# TODO: highlight scatterpoint based on selected row
# source_table.selected.on_change("indices", lambda attr, old, new: update_highlight())
# source_table.selected.indices -> could likely be used to highlight point

controls = column(parameter, station, flag_list, button, width=170)
tables = column(
    data_table_title, data_table, data_table_changed_title, data_table_changed
)

curdoc().add_root(row(controls, tables, plot_ssscc, plot_all))
curdoc().title = "CTDO Data Flagging Tool"

update_selectors()
