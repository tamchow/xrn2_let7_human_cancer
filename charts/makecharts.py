#%%
# Hack to get autocompletion working in VSCode
try:
    from .charts import *
except ImportError:
    from charts import *

raw_data = pd.read_excel("Raw data excel.xlsx", engine="openpyxl", sheet_name=None)

a549_pre_to_pri, _ = make_bar_figure(
    df=raw_data["A549_pre2pri"],
    x="miRNA",
    xlabel=False,
    y="Ratio of pre-miRNA to pri-miRNA",
    hue="Group",
    figsize=(2, 3),
    filename="A549_pre_to_pri",
    title="A549",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=(
        lambda fig, ax: (
            fig,
            rescale_yaxis(move_legend(ax, loc="bottom_out", offset=0.2), ymax=1.5),
        )
    ),
)
huh7_pre_to_pri, _ = make_bar_figure(
    df=raw_data["HUH7_pre2pri"],
    x="miRNA",
    xlabel=False,
    y="Ratio of pre-miRNA to pri-miRNA",
    hue="Group",
    figsize=(2, 3),
    filename="HUH7_pre_to_pri",
    title="HUH7",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=(
        lambda fig, ax: (
            fig,
            rescale_yaxis(move_legend(ax, loc="bottom_out", offset=0.2), ymax=1.6),
        )
    ),
)
u87_pre_to_pri, _ = make_bar_figure(
    df=raw_data["U87_pre2pri"],
    x="miRNA",
    xlabel=False,
    y="Ratio of pre-miRNA to pri-miRNA",
    hue="Group",
    figsize=(2, 3),
    filename="U87_pre_to_pri",
    title="U87",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=(
        lambda fig, ax: (
            fig,
            rescale_yaxis(move_legend(ax, loc="bottom_out", offset=0.2), ymax=1.6),
        )
    ),
)
a549_targets_1, _ = make_bar_figure(
    df=raw_data["A549_targets_1"],
    x="Target",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Group",
    figsize=(5, 2),
    stat_high=False,
    title="A549",
    filename="A549_targets_1",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (fig, move_legend(ax, loc="bottom_out", ncol=2)),
)
a549_targets_2, _ = make_bar_figure(
    df=raw_data["A549_targets_2"],
    x="Target",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Group",
    figsize=(3, 3),
    stat_high=False,
    title="A549",
    filename="A549_targets_2",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (fig, move_legend(ax, loc="bottom_out", ncol=2)),
)
huh7_targets_1, _ = make_bar_figure(
    df=raw_data["HUH7_targets_1"],
    x="Target",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Group",
    figsize=(7, 2.5),
    stat_high=False,
    title="HUH7",
    filename="HUH7_targets_1",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (
        fig,
        rescale_yaxis(
            move_legend(ax, loc="bottom_out", ncol=2), ymax=1.5, tick_gap=0.3
        ),
    ),
)
huh7_targets_2, _ = make_bar_figure(
    df=raw_data["HUH7_targets_2"],
    x="Target",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Group",
    figsize=(4, 2.5),
    stat_high=False,
    title="HUH7",
    filename="HUH7_targets_2",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (
        fig,
        rescale_yaxis(
            move_legend(ax, loc="bottom_out", ncol=2), ymax=1.5, tick_gap=0.3
        ),
    ),
)
u87_targets_1, _ = make_bar_figure(
    df=raw_data["U87_targets_1"],
    x="Target",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Group",
    figsize=(5, 2),
    stat_high=False,
    title="U87",
    filename="U87_targets_1",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (fig, move_legend(ax, loc="bottom_out", ncol=2)),
)
u87_targets_2, _ = make_bar_figure(
    df=raw_data["U87_targets_2"],
    x="Target",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Group",
    figsize=(3, 3),
    stat_high=False,
    title="U87",
    filename="U87_targets_2",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (fig, move_legend(ax, loc="bottom_out", ncol=2)),
)
A549_colonies, _ = make_box_figure(
    df=raw_data["A549_colonies"],
    x="Treatment",
    xlabel=False,
    y="Colonies formed/1000 cells",
    figsize=(3, 3),
    filename="A549_colonies",
    palette=("tab:orange", "tab:green", "tab:blue"),
)
HUH7_colonies, _ = make_box_figure(
    df=raw_data["HUH7_colonies"],
    x="Treatment",
    xlabel=False,
    y="Colonies formed/1000 cells",
    figsize=(3, 3),
    filename="HUH7_colonies",
    palette=("tab:orange", "tab:green", "tab:blue"),
)
a549_pri_miRNAs, _ = make_bar_figure(
    df=raw_data["A549_pri"],
    x="pri-miRNA",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Group",
    figsize=(2, 3),
    title="A549",
    filename="A549_pri_miRNAs",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (
        fig,
        rescale_yaxis(move_legend(ax, loc="bottom_out"), ymax=1.5),
    ),
)
huh7_pri_miRNAs, _ = make_bar_figure(
    df=raw_data["HUH7_pri"],
    x="pri-miRNA",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Group",
    figsize=(2, 3),
    title="HUH7",
    filename="HUH7_pri_miRNAs",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (
        fig,
        rescale_yaxis(move_legend(ax, loc="bottom_out"), ymax=1.8, tick_gap=0.2),
    ),
)
u87_pri_miRNAs, _ = make_bar_figure(
    df=raw_data["U87_pri"],
    x="pri-miRNA",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Group",
    figsize=(2, 3),
    title="U87",
    filename="U87_pri_miRNAs",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (
        fig,
        rescale_yaxis(move_legend(ax, loc="bottom_out"), ymax=1.8, tick_gap=0.2),
    ),
)
miRNAs_a549, _ = make_bar_figure(
    df=raw_data["A549_miRNA"],
    x="miRNA",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Group",
    figsize=(8, 3),
    filename="A549_miRNAs",
    title="A549",
    palette=("tab:orange", "tab:green", "tab:blue"),
    yaxis_break=((0, 1.6), (1.6, 4.5)),
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="bottom_out", ncol=2),
    ),
)
miRNAs_huh7, _ = make_bar_figure(
    df=raw_data["HUH7_miRNA"],
    x="miRNA",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Group",
    figsize=(8, 3),
    filename="HUH7_miRNAs",
    title="HUH7",
    palette=("tab:orange", "tab:green", "tab:blue"),
    yaxis_break=((0, 1.6), (1.6, 4.5)),
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="bottom_out", ncol=2),
    ),
)
miRNAs_u87, _ = make_bar_figure(
    df=raw_data["U87_miRNA"],
    x="miRNA",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Group",
    figsize=(8, 3),
    filename="U87_miRNAs",
    title="U87",
    yaxis_break=((0, 1.6), (1.6, 4.5)),
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="bottom_out", ncol=2),
    ),
)
stability, _ = make_line_figure(
    df=raw_data["stability"],
    x="Time (hours)",
    y="Relative Expression",
    hue="miRNA",
    figsize=(4, 4),
    filename="stability",
    title="",
)
miRNAs_ago2ip, _ = make_bar_figure(
    df=raw_data["ago2ip"],
    x="miRNA",
    xlabel=False,
    y="Relative miRNA level\n(Normalized)",
    hue="Group",
    figsize=(7, 3),
    filename="AGO2IP_miRNAs",
    title="AGO2 IP from A549 cells",
    significance_groups=("Control", "XRN2 shRNA"),
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="bottom_out", ncol=3, offset=0.075),
    ),
    palette=("tab:orange", "tab:blue", "tab:green"),
)
a549_scratch_assay, _ = make_line_figure(
    df=raw_data["A549_wound"],
    x="Time (hours)",
    y="Gap Size (%)",
    hue="Treatment",
    figsize=(4, 4),
    filename="A549_scratch_assay",
    transform=lambda y: y * 100,
    title="",
    palette=("tab:orange", "tab:green"),
)
HUH7_scratch_assay, _ = make_line_figure(
    df=raw_data["HUH7_wound"],
    x="Time (hours)",
    y="Gap Size (%)",
    hue="Treatment",
    figsize=(4, 4),
    filename="HUH7_scratch_assay",
    transform=lambda y: y * 100,
    title="",
    palette=("tab:orange", "tab:green"),
)
a549_addtl_miRNAs, _ = make_bar_figure(
    df=raw_data["A549_addtl_miRNAs"],
    x="miRNA",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Treatment",
    figsize=(2, 3),
    filename="A549_addtl_miRNAs",
    title="A549",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (
        fig,
        rescale_yaxis(
            move_legend(ax, loc="bottom_out", ncol=1, offset=0.3), ymax=5.5, tick_gap=1
        ),
    ),
)
huh7_addtl_miRNAs, _ = make_bar_figure(
    df=raw_data["HUH7_addtl_miRNAs"],
    x="miRNA",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Treatment",
    figsize=(2, 3),
    filename="HUH7_addtl_miRNAs",
    title="HUH7",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (
        fig,
        rescale_yaxis(
            move_legend(ax, loc="bottom_out", ncol=1, offset=0.3), ymax=5.5, tick_gap=1
        ),
    ),
)
u87_addtl_miRNAs, _ = make_bar_figure(
    df=raw_data["U87_addtl_miRNAs"],
    x="miRNA",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Treatment",
    figsize=(2, 3),
    filename="U87_addtl_miRNAs",
    title="U87",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (
        fig,
        rescale_yaxis(
            move_legend(ax, loc="bottom_out", ncol=1, offset=0.3), ymax=5.5, tick_gap=1
        ),
    ),
)
a549_emt, _ = make_bar_figure(
    df=raw_data["A549_emt"],
    x="mRNA",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Treatment",
    stat_high=None,
    figsize=(4, 2),
    filename="A549_emt",
    title="A549",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="bottom_out", ncol=2, offset=0.2),
    ),
)
a549_mmp, _ = make_bar_figure(
    df=raw_data["A549_mmp"],
    x="mRNA",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Treatment",
    stat_high=None,
    figsize=(3, 2),
    filename="A549_mmp",
    title="A549",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (
        fig,
        rescale_yaxis(
            move_legend(ax, loc="bottom_out", ncol=2, offset=0.2),
            ymax=1.2,
            tick_gap=0.3,
        ),
    ),
)
huh7_emt, _ = make_bar_figure(
    df=raw_data["HUH7_emt"],
    x="mRNA",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Treatment",
    stat_high=None,
    figsize=(4, 2),
    filename="HUH7_emt",
    title="HUH7",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="bottom_out", ncol=2, offset=0.2),
    ),
)
huh7_mmp, _ = make_bar_figure(
    df=raw_data["HUH7_mmp"],
    x="mRNA",
    xlabel=False,
    y="Normalized Fold Change",
    hue="Treatment",
    stat_high=None,
    figsize=(3, 2),
    filename="HUH7_mmp",
    title="HUH7",
    palette=("tab:orange", "tab:green", "tab:blue"),
    plot_postprocess=lambda fig, ax: (
        fig,
        rescale_yaxis(
            move_legend(ax, loc="bottom_out", ncol=2, offset=0.2),
            ymax=1.2,
            tick_gap=0.3,
        ),
    ),
)
release_rna_1, _ = make_bar_figure(
    df=raw_data["rna_release_1"],
    x="miRNA",
    # xlabel='Target RNA to let-7a',
    xlabel=False,
    y="Fraction of AGO-bound miRNA\n(Normalized)",
    hue="Group",
    skip_statannot=True,
    use_orig_order=True,
    figsize=(10, 4),
    filename="release_rna_1",
    palette=consistent_palette,
    title="",
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="bottom_out", ncol=3, expand=True, offset=0.25),
    ),
)
release_rna_2, _ = make_bar_figure(
    df=raw_data["rna_release_2"],
    x="miRNA",
    # xlabel='Target RNA to let-7a',
    xlabel=False,
    y="Fraction of AGO-bound miRNA\n(Normalized)",
    hue="Group",
    skip_statannot=True,
    use_orig_order=True,
    figsize=(7, 4),
    filename="release_rna_2",
    palette=consistent_palette,
    title="",
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="bottom_out", ncol=3, expand=True, offset=0.25),
    ),
)
release_rna_3, _ = make_bar_figure(
    df=raw_data["rna_release_3"],
    x="miRNA",
    # xlabel='Target RNA to let-7a',
    xlabel=False,
    y="Fraction of AGO-bound miRNA\n(Normalized)",
    hue="Group",
    skip_statannot=True,
    use_orig_order=True,
    figsize=(4, 4),
    filename="release_rna_3",
    palette=consistent_palette,
    title="",
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="right_out", offset=0.5),
    ),
)
release_protein_1, _ = make_bar_figure(
    df=raw_data["protein_release_2"],
    x="Treatment",
    xlabel=False,
    y="Fraction of AGO-bound miRNA\n(Normalized)",
    hue="Group",
    skip_statannot=True,
    use_orig_order=True,
    figsize=(10, 4),
    filename="release_protein_1",
    palette=consistent_palette,
    title="",
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="bottom_out", ncol=4, offset=0.075),
    ),
)
release_protein_2, _ = make_bar_figure(
    df=raw_data["protein_release_3"],
    x="Treatment",
    xlabel=False,
    y="Fraction of AGO-bound miRNA\n(Normalized)",
    hue="Group",
    skip_statannot=True,
    use_orig_order=True,
    figsize=(6, 3),
    filename="release_protein_2",
    palette=consistent_palette,
    title="",
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="bottom_out", ncol=2, offset=0.075),
    ),
)
release_protein_3, _ = make_bar_figure(
    df=raw_data["protein_release_4"],
    x="Treatment",
    xlabel=False,
    y="Fraction of AGO-bound miRNA\n(Normalized)",
    hue="Group",
    skip_statannot=True,
    use_orig_order=True,
    figsize=(6, 3),
    filename="release_protein_3",
    palette=consistent_palette,
    title="",
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="bottom_out", ncol=2, offset=0.075),
    ),
)
release_protein_4, _ = make_bar_figure(
    df=raw_data["protein_release_5"],
    x="Treatment",
    xlabel=False,
    y="Fraction of AGO-bound miRNA\n(Normalized)",
    hue="Group",
    skip_statannot=True,
    use_orig_order=True,
    figsize=(4, 4),
    filename="release_protein_4",
    palette=consistent_palette,
    title="",
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="bottom_out", ncol=2, expand=True, offset=0.075),
    ),
)
release_protein_5, _ = make_bar_figure(
    df=raw_data["protein_release_6"],
    x="Treatment",
    xlabel=False,
    y="Fraction of AGO-bound miRNA\n(Normalized)",
    hue="Group",
    skip_statannot=True,
    use_orig_order=True,
    figsize=(3, 4),
    filename="release_protein_5",
    palette=consistent_palette,
    title="",
    plot_postprocess=lambda fig, ax: (
        fig,
        move_legend(ax, loc="bottom_out", ncol=2, expand=True, offset=0.075),
    ),
)
# release_comp, _ = make_bar_figure(
#     df=raw_data["release_comp (2)"],
#     x="Treatment",
#     xlabel=False,
#     y="Fraction of AGO-bound miRNA\n(Normalized)",
#     hue="Group",
#     skip_statannot=True,
#     figsize=(8, 4),
#     filename="release_comparison",
#     title="",
#     palette=("tab:blue",),
#     plot_postprocess=lambda fig, ax: (
#         fig,
#         move_legend(ax, loc="bottom_out", ncol=1),
#     ),
# )
# %%
