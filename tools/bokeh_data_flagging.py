import pandas as pd
import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.plotting import figure, show
from bokeh.models import (
    Button,
    ColumnDataSource,
    # CustomJS,
    DataTable,
    Slider,
    TableColumn,
    Select,
    MultiSelect,
)


### intialize data
df = pd.read_csv("../all_rinko_merged.csv")
df["Residual"] = df["REFOXY_rinko"] - df["CTDRINKO"]
df["Flag"] = 1
ssscc_list = df["SSSCC_rinko"].unique()

### intialize widgets
button = Button(label="Save flagged data", button_type="success")
# button.js_on_click(CustomJS(args=dict(source=source),
#                             code=open(join(dirname(__file__), "download.js")).read()))
parameter = Select(
    title="Parameter", options=["CTDTMP", "CTDCOND", "SBE43", "RINKO"], value="RINKO"
)
station = Select(title="Station", options=[str(x) for x in ssscc_list], value="00101")
station.on_change("value", lambda attr, old, new: update())
# explanation of flags:
# https://cchdo.github.io/hdo-assets/documentation/manuals/pdf/90_1/chap4.pdf
flag_list = MultiSelect(
    title="Display data flagged as:",
    value=["1", "2", "3"],
    options=[
        ("1", "1 [Uncalibrated]"),
        ("2", "2 [Acceptable]"),
        ("3", "3 [Questionable]"),
        ("4", "4 [Bad]"),
    ],
)

source_table = ColumnDataSource(data=dict())
source_plot = ColumnDataSource(data=dict(x=[], y=[]))

### set up plots
x_var = "RINKO"
ssscc = "00101"
plot = figure(
    plot_height=800,
    plot_width=500,
    title="{} vs CTDPRS [Station {}]".format(x_var, ssscc),
    tools="crosshair,pan,reset,save,wheel_zoom",
    y_axis_label="Pressure (dbar)",
)
plot.y_range.flipped = True  # invert y-axis
plot.scatter(
    "x",
    "y",
    fill_color="#999999",
    line_color="#000000",
    size=10,
    line_width=2,
    source=source_plot,
)
plot.scatter(df["REFOXY_rinko"], df["CTDPRS_rinko_ctd"])


def update():
    current = df[df["SSSCC_rinko"] == int(station.value)]
    source_table.data = {
        "CTDPRS": current["CTDPRS_rinko_ctd"],
        "REFOXY": current["REFOXY_rinko"],
        "CTDRINKO": current["CTDRINKO"],
        "diff": current["Residual"],
        "flag": current["Flag"],
    }
    source_plot.data = {
        "x": current["REFOXY_rinko"],
        "y": current["CTDPRS_rinko_ctd"],
    }


columns = [
    TableColumn(field="CTDPRS", title="CTDPRS"),
    TableColumn(field="REFOXY", title="REFOXY"),
    TableColumn(field="CTDRINKO", title="CTDRINKO"),
    TableColumn(field="diff", title="Residual"),
    TableColumn(field="flag", title="Flag"),
]

data_table = DataTable(
    source=source_table, columns=columns, width=500, height=800, editable=True,
)

controls = column(parameter, station, flag_list, button, width=170)

curdoc().add_root(row(controls, data_table, plot))
curdoc().title = "CTDO Data Flagging Tool"

update()
