import pandas as pd
import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import (
    Button,
    ColumnDataSource,
    # CustomJS,
    DataTable,
    Slider,
    TableColumn,
    Select,
)

df = pd.read_csv("../all_rinko_merged.csv")
df["Residual"] = df["REFOXY_rinko"] - df["CTDRINKO"]
df["Flag"] = 1
ssscc_list = df["SSSCC_rinko"].unique()

parameter = Select(
    title="Parameter", options=["CTDTMP", "CTDCOND", "SBE43", "RINKO"], value="RINKO"
)
station = Select(title="Station", options=[str(x) for x in ssscc_list], value="00101")
source = ColumnDataSource(data=dict())


def update():
    current = df[df["SSSCC_rinko"] == ssscc_list[slider.value]]
    source.data = {
        "CTDPRS": current["CTDPRS_rinko_ctd"],
        "REFOXY": current["REFOXY_rinko"],
        "CTDRINKO": current["CTDRINKO"],
        "diff": current["Residual"],
        "flag": current["Flag"],
    }


slider = Slider(
    title="Station", start=1, end=len(ssscc_list), value=1, step=1, format="0"
)
slider.on_change("value", lambda attr, old, new: update())

button = Button(label="Save flagged data", button_type="success")
# button.js_on_click(CustomJS(args=dict(source=source),
#                             code=open(join(dirname(__file__), "download.js")).read()))

columns = [
    TableColumn(field="CTDPRS", title="CTDPRS"),
    TableColumn(field="REFOXY", title="REFOXY"),
    TableColumn(field="CTDRINKO", title="CTDRINKO"),
    TableColumn(field="diff", title="Residual"),
    TableColumn(field="flag", title="Flag"),
]

data_table = DataTable(
    source=source, columns=columns, width=500, height=800, editable=True
)

controls = column(parameter, station, slider, button)

curdoc().add_root(row(controls, data_table))
curdoc().title = "CTDO Data Flagging Tool"

update()
