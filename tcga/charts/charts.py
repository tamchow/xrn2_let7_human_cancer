import builtins
import itertools as it
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import scipy.stats as stats
import statannot as sa

debug = True


def print(*args, **kwargs):
    if debug:
        builtins.print(args, kwargs)


def dedup_col_name(name):
    if "." in name:
        return name[: name.index(".")]
    else:
        return name


def safe_dropna(df):
    df.dropna(how="all", axis="columns", inplace=True)
    df.dropna(how="all", axis="index", inplace=True)
    return df


def clean(df, id_name, variable_name, value_name):
    id_col = df.columns[0]
    df = safe_dropna(df)
    orig_x = df[id_col]

    if id_name:
        df.rename(columns={id_col: id_name}, inplace=True)
    else:
        id_name = id_col
    new_df = df.melt(id_vars=id_name, var_name=variable_name, value_name=value_name)
    new_df[variable_name] = new_df[variable_name].map(dedup_col_name)
    new_df.sort_values(by=variable_name, inplace=True)
    if id_name:
        new_df.sort_values(by=id_name, inplace=True)
    orig_cols = list(
        dict.fromkeys(dedup_col_name(col_name) for col_name in df.columns[1:])
    )
    return new_df, orig_cols, orig_x


def chunk(iterable, size):
    iterable = iter(iterable)
    return iter(lambda: tuple(it.islice(iterable, size)), ())


def get_box_pairs(var_col, group_col):
    if isinstance(var_col, pd.Series):
        var_col = var_col.unique()
    if isinstance(group_col, pd.Series):
        group_col = group_col.unique()

    return list(chunk(it.product(var_col, group_col), size=len(group_col)))


default_palette = list(mcolors.TABLEAU_COLORS.keys())

consistent_palette = ["tab:orange", "tab:green", "tab:blue"]
consistent_palette = (
    consistent_palette[:-1]
    + [color for color in default_palette if color not in consistent_palette]
    + [consistent_palette[-1]]
)


def set_labels(fig, axes, x, xlabel, y, ylabel, title, yaxis_break=None):
    if ylabel is None and ylabel is not False:
        ylabel = y
    elif ylabel is False:
        ylabel = None
    if xlabel is None and xlabel is not False:
        xlabel = x
    elif xlabel is False:
        xlabel = None
    if yaxis_break:
        fig.supxlabel(xlabel, y=0)
        fig.supylabel(ylabel, x=0.05)
        fig.suptitle(title, y=1.02)
    else:
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)
    return fig, axes


def savefig(fig, filename):
    if filename:
        fig.savefig(f"./Figures/Raw Images/{filename}.svg", bbox_inches="tight")
    return fig


def make_bar_figure(
    df,
    x,
    y,
    hue,
    figsize=None,
    filename=None,
    xlabel=None,
    ylabel=None,
    title=None,
    significance_groups=None,
    stat_high=True,
    hline=False,
    use_stripplot=False,
    skip_statannot=False,
    order=None,
    use_orig_order=False,
    palette=default_palette,
    plot_postprocess=None,
    yaxis_break=None,
):
    df, orig_order, orig_x = clean(df, x, hue, y)
    order = order or np.sort(df[hue].unique()).tolist()
    if use_orig_order:
        order = orig_order
    fig = plt.figure(dpi=300, figsize=figsize)
    ax1 = None
    if use_stripplot:
        ax1 = sns.stripplot(
            ax=fig.gca(),
            data=df,
            x=x,
            y=y,
            hue=hue,
            order=orig_x,
            hue_order=order,
            edgecolor="black",
            linewidth=1,
            size=10,
            palette=palette,
            jitter=0.25,
        )
    else:
        if yaxis_break:
            height_ratios = [ylims[1] - ylims[0] for ylims in yaxis_break]
            gridspec_kws = {
                "ncols": 1,
                "nrows": 2,
                "hspace": 0.05,
                "height_ratios": height_ratios,
            }
            gs = fig.add_gridspec(**gridspec_kws)
            ax_top = fig.add_subplot(gs[0])
            ax_bottom = fig.add_subplot(gs[1], sharex=ax_top)
            ax_top.set_ylim(yaxis_break[1])
            ax_top = sns.barplot(
                ax=ax_top,
                data=df,
                x=x,
                y=y,
                hue=hue,
                order=orig_x,
                hue_order=order,
                errwidth=1,
                capsize=0.05,
                ci="sd",
                edgecolor="white",
                linewidth=2,
                palette=palette,
            )
            ax_bottom.set_ylim(yaxis_break[0])
            ax_bottom = sns.barplot(
                ax=ax_bottom,
                data=df,
                x=x,
                y=y,
                hue=hue,
                order=orig_x,
                hue_order=order,
                errwidth=1,
                capsize=0.05,
                ci="sd",
                edgecolor="white",
                linewidth=2,
                palette=palette,
            )
            ax_bottom.set_ylabel(None)
            ax_bottom.set_xlabel(None)
            ax_top.set_ylabel(None)
            ax_top.set_xlabel(None)
            ax_top.get_xaxis().set_visible(False)
            ax_bottom.get_legend().remove()
            d = 0.015  # how big to make the diagonal lines in axes coordinates
            # arguments to pass to plot, just so we don't keep repeating them
            kwargs = dict(
                transform=ax_top.transAxes, color="k", linewidth=0.75, clip_on=False
            )
            ax_top.plot((-d, +d), (0, 0), **kwargs)  # top-left break line
            kwargs.update(transform=ax_bottom.transAxes)
            ax_bottom.plot((-d, +d), (1, 1), **kwargs)  # bottom-left break line
            ax1 = ax_top
        else:
            ax1 = sns.barplot(
                ax=fig.gca(),
                data=df,
                x=x,
                y=y,
                hue=hue,
                order=orig_x,
                hue_order=order,
                errwidth=1,
                capsize=0.05,
                ci="sd",
                edgecolor="white",
                linewidth=2,
                palette=palette,
            )

    if not skip_statannot:
        box_pairs = get_box_pairs(
            df[x], significance_groups if significance_groups is not None else df[hue]
        )

        def single_sided_ttest():
            pvalues = []
            for pair in box_pairs:
                data1 = df[df[x] == pair[0][0]][df[hue] == pair[0][1]][y]
                data2 = df[df[x] == pair[1][0]][df[hue] == pair[1][1]][y]
                stat, p = stats.ttest_ind(data2, data1, nan_policy="omit")
                if stat_high is not None:
                    if stat_high:
                        p = p if stat > 0 else 1
                    else:
                        p = p if stat < 0 else 1
                # print("T-Test(independent) stat={:.2e} p-value={:.2e}".format(stat, p))
                pvalues.append(p)
            return pvalues

        if yaxis_break:
            ax1, test_results = sa.add_stat_annotation(
                ax1,
                data=df,
                x=x,
                y=y,
                hue=hue,
                hue_order=order,
                box_pairs=box_pairs,
                perform_stat_test=False,
                pvalues=single_sided_ttest(),
                linewidth=1,
                loc="outside",
            )
        else:
            ax1, test_results = sa.add_stat_annotation(
                ax1,
                data=df,
                x=x,
                y=y,
                hue=hue,
                hue_order=order,
                box_pairs=box_pairs,
                perform_stat_test=False,
                pvalues=single_sided_ttest(),
                linewidth=1,
            )

    ax1_legend = ax1.get_legend()
    ax1_legend.set_title(None)
    ax1.set_xticklabels(
        [label.get_text().replace("\\n", "\n") for label in ax1.get_xticklabels()]
    )
    fig, ax1 = set_labels(fig, ax1, x, xlabel, y, ylabel, title, yaxis_break)
    if isinstance(hline, (int, float)):
        ax1.axhline(hline, color="black", linewidth="0.5")
    if plot_postprocess:
        fig, ax1 = plot_postprocess(fig, ax1)
    sns.despine(fig=fig)
    if yaxis_break:
        sns.despine(ax=ax1, bottom=True)
    fig = savefig(fig, filename)
    return df, ax1


def make_box_figure(
    df,
    x,
    y,
    figsize=None,
    filename=None,
    xlabel=None,
    ylabel=None,
    title=None,
    palette=default_palette,
):
    df, _, _ = clean(df, None, x, y)
    fig = plt.figure(dpi=300, figsize=figsize)
    ax1 = sns.boxplot(ax=fig.gca(), data=df, x=x, y=y, palette=palette)
    box_pairs = list(chunk(df[x].unique(), size=2))
    ax1, test_results = sa.add_stat_annotation(
        ax1,
        data=df,
        x=x,
        y=y,
        box_pairs=box_pairs,
        test="t-test_ind",
    )
    fig, ax1 = set_labels(fig, ax1, x, xlabel, y, ylabel, title)
    sns.despine(fig=fig)
    fig = savefig(fig, filename)
    return df, ax1


def make_line_figure(
    df,
    x=None,
    y=None,
    hue=None,
    figsize=None,
    filename=None,
    xlabel=None,
    ylabel=None,
    transform=None,
    title=None,
    palette=None,
):
    df, markers, _ = clean(df, x, hue, y)
    markers = list(dict.fromkeys(markers))
    if transform:
        df[y] = df[y].apply(transform)
    fig = plt.figure(dpi=300, figsize=figsize)
    ax1 = fig.gca()
    from matplotlib.lines import Line2D

    ax1 = sns.lineplot(
        data=df,
        x=x,
        y=y,
        hue=hue,
        hue_order=markers,
        style=hue,
        style_order=markers,
        palette=palette,
        err_style="bars",
        err_kws={"capsize": 4.0, "capthick": 0.5, "elinewidth": 0.5},
        dashes=False,
        markers=dict(zip(markers, Line2D.filled_markers)),
    )
    fig, ax1 = set_labels(fig, ax1, x, xlabel, y, ylabel, title)
    sns.despine(fig=fig)
    fig = savefig(fig, filename)
    return df, ax1


def volcano_plot(
    df,
    pval_col="pval",
    quant_col="measure",
    label_col="label",
    figsize=None,
    filename=None,
    log2_quant=False,
    xlabel=None,
    ylabel=None,
    colors=("red", "gray", "blue"),
    pval_cutoff=(0.01, 0.01),
    quant_cutoff=(-0.05, 0.05),
    clip_limits=(1e-10, 10),
    lines=False,
    labels=None,
    label_neg_outliers=0,
    label_pos_outliers=0,
    plot_postprocess=None,
    title_item=''
):
    def mark_change(q, p):
        return (
            "-ve"
            if q <= quant_cutoff[0] and p <= pval_cutoff[0]
            else "+ve"
            if q >= quant_cutoff[1] and p <= pval_cutoff[0]
            else "/"
        )

    df["change"] = df.apply(lambda x: mark_change(x[quant_col], x[pval_col]), axis=1)

    df["mlog10pval"] = -np.log10(np.clip(df[pval_col], clip_limits[0], clip_limits[1]))

    if log2_quant:
        df["log2quant"] = np.log2(df[quant_col])
        quant_col = "log2quant"

    dist_scale_factor = np.mean(df["mlog10pval"]) / np.mean(abs(df[quant_col]))
    df["distance"] = df.apply(
        lambda x: np.hypot(dist_scale_factor * abs(x[quant_col]), x["mlog10pval"]),
        axis=1,
    )

    fig = plt.figure(dpi=300, figsize=figsize)
    ax = fig.gca()
    ax = sns.scatterplot(
        x=quant_col,
        y="mlog10pval",
        hue="change",
        data=df,
        palette={"-ve": colors[0], "/": colors[1], "+ve": colors[2]},
        linewidth=0.2,
        edgecolor="black",
        s=24,
    )

    grouped_by_dir = df.groupby(by="change")
    neg = (
        grouped_by_dir.get_group("-ve").sort_values(by="distance", ascending=False)
        if "-ve" in grouped_by_dir.groups
        else None
    )
    pos = (
        grouped_by_dir.get_group("+ve").sort_values(by="distance", ascending=False)
        if "+ve" in grouped_by_dir.groups
        else None
    )

    current_handles, current_labels = ax.get_legend_handles_labels()
    new_labels = {current_label: "" for current_label in current_labels}

    def get_count(df, col, val=None):
        if val is None:
            return df[col].count()
        return df.loc[df[col] == val, col].count()

    total_cnt = get_count(df, "change")

    if "-ve" in new_labels:
        new_labels["-ve"] = f"Negative ({get_count(df, 'change', '-ve')})"
    if "/" in new_labels:
        new_labels["/"] = f"Unaltered ({get_count(df, 'change', '/')})"
    if "+ve" in new_labels:
        new_labels["+ve"] = f"Positive ({get_count(df, 'change', '+ve')})"

    new_labels_list = list(new_labels.values())
    ax.legend(
        current_handles,
        new_labels_list,
        title=f"Of {total_cnt}{title_item}:",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.5),
        ncol=len(new_labels_list),
        handletextpad=0.1,
        columnspacing=0.3
    )

    if lines:
        for xc in quant_cutoff:
            ax.axvline(
                x=np.log2(max(0, xc)) if log2_quant else xc,
                linestyle="--",
                color="black",
                linewidth=0.5,
            )
        for yc in pval_cutoff:
            ax.axhline(y=-np.log10(yc), linestyle="--", color="black", linewidth=0.5)

    if labels == None:
        labels = []

    if label_neg_outliers > 0 and (neg is not None):
        print(neg.head(label_neg_outliers))
        labels += neg.head(label_neg_outliers)[label_col].tolist()
    if label_pos_outliers > 0 and (pos is not None):
        print(neg.head(label_pos_outliers))
        labels += pos.head(label_pos_outliers)[label_col].tolist()

    all_labels = set(df[label_col].to_list())
    if labels:
        for label in labels:
            if label not in all_labels:
                continue
            text_x = df.loc[
                    df[label_col] == label, "log2quant" if log2_quant else quant_col
                ].values[0]
            text_y = df.loc[df[label_col] == label, "mlog10pval"].values[0]
            ax.text(
                x=text_x,
                y=text_y,
                s=label,
                fontsize="xx-small",
                color="black",
            )

    import matplotlib.ticker as mticker

    ax.grid(False)
    ax.yaxis.set_major_locator(mticker.FixedLocator(ax.get_yticks().tolist()))
    ax.set_yticklabels([f"{ts:.1f}" for ts in ax.get_yticks()])

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    if plot_postprocess:
        fig, ax = plot_postprocess(fig, ax)

    fig = savefig(fig, filename)
    return df, ax


def grouped_boxplots(
    df, id_col, target_col, value_col, group_col, figsize=None, filename=None
):
    fig = plt.figure(dpi=300, figsize=figsize)
    sns.set(font_scale=1.25, style="ticks")
    ax = fig.gca()
    order = df[target_col].unique().tolist()
    groups = list(reversed(np.sort(df[group_col].unique())))
    palette = ("tab:orange", "tab:green", "tab:blue")
    ax = sns.boxplot(
        data=df,
        x=target_col,
        y=value_col,
        hue=group_col,
        order=order,
        ax=ax,
        dodge=True,
        width=0.8,
        fliersize=0,
        zorder=5,
        palette=palette,
    )
    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, 0.8))
    ax = sns.stripplot(
        data=df,
        x=target_col,
        y=value_col,
        hue=group_col,
        order=order,
        ax=ax,
        dodge=True,
        linewidth=0.5,
        edgecolor="gray",
        zorder=2,
        marker=".",
        palette=palette,
    )
    ax.set_xlabel(None)
    handles, labels = ax.get_legend_handles_labels()
    labels = [
        f"{group} ({len(df[df[group_col] == group][id_col].unique())}/{len(df[id_col].unique())})"
        for group in labels
    ]
    plt.legend(
        handles[0 : len(groups)],
        labels,
        loc="center",
        bbox_to_anchor=(0.5, -0.2),
        ncol=2,
    )

    box_pairs = list(chunk(it.product(order, groups), len(groups)))

    test_results = sa.add_stat_annotation(
        ax,
        data=df,
        x=target_col,
        y=value_col,
        hue=group_col,
        order=order,
        box_pairs=box_pairs,
        test="t-test_ind",
        text_format="star",
        loc="inside",
        verbose=2,
    )
    sns.despine(ax=ax)
    fig, savefig(fig, filename)
    return ax


def survival_analysis_plot(
    ax, survival_data, pval, labels, palette=("tab:green", "tab:orange")
):
    for survival_datum, label, color in zip(survival_data, labels, palette):
        series_x, series_y = zip(*((x[2] / 30, x[1]) for x in survival_datum))
        ax = sns.lineplot(
            x=series_x,
            y=series_y,
            label=label,
            drawstyle="steps-post",
            ax=ax,
            color=color,
        )

    ax.legend(loc="upper right")
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Proportion Surviving")
    ax.set_title(f"P-Value = {round(pval, 3)}")
    return ax


def make_violinplots(
    df,
    mirna_cols,
    mrna_cols,
    figsize=(10, 5),
    filename=None,
    miRNA_label="miRNA (Reads Per Million)",
    mRNA_label="mRNA (FPKM)",
):
    fig = plt.figure(dpi=300, figsize=figsize)
    mirna_ax = fig.add_axes((0, 0, 0.46, 1))
    mrna_ax = fig.add_axes((0.54, 0, 0.46, 1))
    mirna_ax = sns.violinplot(data=df[mirna_cols], ax=mirna_ax, scale="width")
    mirna_ax.yaxis.grid(True)
    mirna_ax.set_ylabel(miRNA_label)
    mrna_ax = sns.violinplot(data=df[mrna_cols], ax=mrna_ax, scale="width")
    mrna_ax.yaxis.grid(True)
    mrna_ax.set_ylabel(mRNA_label)
    fig = savefig(fig, filename)
    return (mirna_ax, mrna_ax)


def make_heatmap(df, basis_col=None, columns=None, cutoff=None, clip_outliers=(0, 0)):
    def argmedian(x):
        return np.argpartition(x, len(x) // 2)[len(x) // 2]

    df = df.copy()
    df = df.sort_values(by=basis_col)
    low_group = df.iloc[: argmedian(df[basis_col])]
    high_group = df.iloc[argmedian(df[basis_col]) :]
    sns.set(font_scale=0.4)
    fig = plt.figure(dpi=600, figsize=(3.5, 2))
    hmap = fig.gca()

    if basis_col is None or not any(columns):
        raise ValueError("No columns provided for heatmap generation")
    else:
        if clip_outliers[0] > 0 or clip_outliers[1] > 0:
            for column in columns:
                low_group[column] = (
                    low_group[column].sort_values().iloc[clip_outliers[0] :]
                )
                high_group[column] = (
                    high_group[column].sort_values().iloc[clip_outliers[1] :]
                )
            df = pd.concat((low_group, high_group)).sort_values(by=basis_col)
            df.interpolate(axis=0, inplace=True)
        if cutoff and 0 < cutoff < len(columns):
            print(f"Dropping {len(columns) - cutoff} cols to achieve {cutoff}")
            columns = columns[:cutoff]
        hmap = sns.heatmap(
            pd.concat((df[basis_col], df[columns]), axis=1),
            ax=hmap,
            robust=True,
            cmap="viridis",
            cbar=True,
        )
        hmap.axvline(x=1, color="white", linewidth=2)
        hmap.axhline(y=argmedian(df[basis_col]), color="white", linewidth=1)
    hmap.set(yticks=[])
    hmap.set(ylabel="")
    return hmap


def gse_plot(gse_file, size_func, size_col_name=None):
    gse_data = pd.read_csv(gse_file, index_col="term_id")
    gse_data["percent_mapped"] = 100 * (
        gse_data["intersection_size"] / gse_data["query_size"]
    )
    gse_data["point_size"] = [
        size_func(intersection.split(",")) for intersection in gse_data["intersections"]
    ]
    new_term_col = "Enriched GO Terms"
    new_pt_size_col = "point_size" if size_col_name is None else size_col_name
    gse_data.rename(
        columns={"term_name": new_term_col, "point_size": new_pt_size_col}, inplace=True
    )
    sns.set_theme(style="ticks", font_scale=2.5)
    fig = plt.figure(dpi=600, figsize=(10, 10))
    default_markersize = plt.rcParams["lines.markersize"]
    plt.rcParams["lines.markersize"] = default_markersize * 3
    ax = fig.gca()
    ax = sns.scatterplot(
        data=gse_data,
        x="percent_mapped",
        y="negative_log10_of_adjusted_p_value",
        hue=new_term_col,
        size=new_pt_size_col,
        ax=ax,
        palette=sns.color_palette("husl", len(gse_data)),
        linewidth=1,
        edgecolor="black",
    )
    plt.rcParams["lines.markersize"] = default_markersize
    handles, labels = ax.get_legend_handles_labels()

    last_term_legend_idx = len(gse_data[new_term_col]) + 1  # +1 for title
    title_1_mm = labels[0].replace(" ", "\\ ")
    labels[0] = f"$\\hspace{{8.0}}\\bf{{\\mathrm{{{title_1_mm}}}}}$"
    title_2_mm = labels[last_term_legend_idx].replace(" ", "\\ ")
    labels[last_term_legend_idx] = f"$\\hspace{{8.0}}\\bf{{\\mathrm{{{title_2_mm}}}}}$"
    ax.legend(
        handles, labels, bbox_to_anchor=(1.025, 1), loc="upper left", labelspacing=0.5
    )

    import matplotlib.ticker as mticker

    ax.xaxis.set_major_locator(mticker.FixedLocator(ax.get_xticks().tolist()))
    ax.set_xticklabels([f"{xt}%" for xt in ax.get_xticks()])

    ax.set_xlabel("Target genes mapped to GO term")
    ax.set_ylabel("$-\\log_{10}\\left(p_{adj}\\right)$")
    sns.set_theme(style="ticks", font_scale=1)
    return ax


def rescale_yaxis(ax, ymin=None, ymax=None, tick_gap=None):
    axymin, axymax = ax.get_ylim()
    if ymin is None:
        ymin = axymin
    if ymax is None:
        ymax = axymax
    axyticks = ax.get_yticks()
    if tick_gap is None:
        tick_gap = abs(axyticks[1] - axyticks[0])
    ax.set_ylim((ymin, ymax))
    ax.set_yticks(np.arange(ymin, ymax + tick_gap, tick_gap))
    return ax


def flip(items, ncol):
    return it.chain(*[items[i::ncol] for i in range(ncol)])


def move_legend(
    ax, bbox=None, loc=None, expand=False, ncol=1, offset=None, figure=False
):
    old_ax = ax
    data = ax.get_legend_handles_labels()

    if expand:
        handles, labels = data
        data = flip(handles, ncol), flip(labels, ncol)

    if figure and loc is None:
        loc = "lower center"
    elif loc is None:
        loc = "center right"
    if loc == "bottom_out":
        if offset is None:
            offset = 0.15
        bbox = (0.5, 0 - offset)
        loc = "lower center"
        figure = True
    elif loc == "right_out":
        if offset is None:
            offset = 0.1
        bbox = (1 + offset, 0.5)
        loc = "center right"
        figure = True

    if figure and not loc == "best":
        fig = ax.figure
        ax.get_legend().remove()
        ax = fig

    if bbox:
        ax.legend(*data, bbox_to_anchor=bbox, loc=loc, ncol=ncol)
    elif loc:
        ax.legend(*data, loc=loc, ncol=ncol)
    return old_ax
