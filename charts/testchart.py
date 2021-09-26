#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge

fig = plt.figure(dpi=300)
ax = fig.gca()

radius_mults = (
    1.115486238,
    1.307838899,
    1.781859099,
    1.101604669,
    1.0,
    1.582363521,
    1.340079679,
    1.033863188,
    1.140432036,
)
color_labels = (
    "XRN2",
    "SND1/TSN",
    "XRN2",
    "XRN2",
    "",
    "XRN2",
    "XRN2",
    "SND1/TSN",
    "XRN2",
)
wedge_labels = (
    "let-7a",
    "let-7b",
    "let-7c",
    "let-7d",
    "let-7e",
    "let-7f",
    "let-7g",
    "let-7i",
    "miR-98",
)

items = min(map(len, (radius_mults, color_labels, wedge_labels)))
angle_inc = 360 / items
start_angle = 90
base_radius = 0.23

cmap = {
    "XRN2"    : "tab:blue"  ,
    "SND1/TSN": "lightcoral",
    "ZSWIM8"  : "lightgreen",
    ""        : "gainsboro" ,
}

ax.set_aspect("equal")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

label_r_mult = 0.7
legend_r_mult, legend_r_offset = 1.07, 0.05

for i in range(items):
    w = Wedge(
        (0.5, 0.5), base_radius * radius_mults[i], i * angle_inc, (i + 1) * angle_inc
    )
    ang = (w.theta2 - w.theta1) / 2 + w.theta1
    cx, cy = w.center
    raw_x, raw_y = np.cos(np.deg2rad(ang)), np.sin(np.deg2rad(ang))
    w_x, w_y = w.r * raw_x, w.r * raw_y
    lab_y = cy + label_r_mult * w_y
    lab_x = cx + label_r_mult * w_x
    x_sgn, y_sgn = np.sign(raw_x), np.sign(raw_y)
    leg_y = cy + legend_r_mult * w_y
    leg_x = cx + legend_r_mult * w_x
    label_text = f"{wedge_labels[i]}\n{radius_mults[i]-1:.1%}"
    legend_text = color_labels[i]
    print(label_text, legend_text, ang)
    ax.add_patch(w)
    ax.annotate(
        label_text, (lab_x, lab_y),
        ha="center",
        va="center",
        fontsize="xx-small"
    )
    ax.annotate(
        legend_text,
        (leg_x, leg_y),
        ha="right"
        if (x_sgn < 0) and (ang >= 110 and ang <= 250)
        else ("center"
        if (x_sgn < 0) and (ang <= 110 or ang >= 250)
        else "left"),
        va="center",
        fontsize="xx-small",
    )
    w.set_color(cmap[color_labels[i]])
    w.set_edgecolor("black")
    w.set_linewidth(0.2)

ax.set_title("HPCC")
fig.savefig("Figures/Raw Images/hpcc_pie.svg", bbox_inches="tight")
# %%
