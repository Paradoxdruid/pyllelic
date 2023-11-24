#!/usr/bin/env python3
"""Utilities to visualize data for use in pyllelic."""

from typing import Dict, List, Optional, Union

import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.subplots as sp
import seaborn as sns
from matplotlib.figure import Figure


def _create_histogram(
    data: pd.DataFrame, cell_line: str, position: str, backend: str
) -> Union[go.Figure, Figure]:
    """Generate a graph figure showing fractional methylation in
    a given cell line at a given site.

    Args:
        data (pd.DataFrame): dataframe of individual data
        cell_line (str): name of cell line
        position (str): genomic position
        backend (str): which plotting backend to use

    Returns:
        Union[go.Figure, Figure]: plotly or matplotlib figure object

    Raises:
        ValueError: invalid plotting backend provided
    """
    if backend == "plotly":
        fig: go.Figure = go.Figure()
        fig.add_trace(
            go.Histogram(
                x=data.loc[cell_line, position],
                xbins={
                    "start": -0.1,
                    "end": 1.1,
                    "size": 0.2,
                },  # offset bins to center displayed bars
            )
        )
        fig.update_layout(
            title_text=f"Methylation at pos: {position} for cell line: {cell_line}",
            xaxis_title_text="Fraction of Sites Methylated in Read",
            yaxis_title_text="Read Count",
            bargap=0.2,
            template="seaborn",
        )

        fig.update_xaxes(range=[-0.1, 1.1])

        return fig

    if backend == "matplotlib":
        sns.set_theme(style="white", rc={"figure.figsize": (6, 4)})
        ax = sns.histplot(
            data.loc[cell_line, position], color="black", discrete=True, shrink=0.2
        )
        ax.set(xlabel="Methylation", ylabel="Count")

        return ax

    raise ValueError("Invalid plotting backend")


def _create_heatmap(
    df: pd.DataFrame,
    min_values: int,
    width: int,
    height: int,
    title_type: str,
    backend: str,
) -> Union[go.Figure, Figure]:
    """Generate a graph figure showing heatmap of mean methylation across
    cell lines.

    Args:
        df (pd.DataFrame): dataframe of mean methylation
        min_values (int): minimum number of points data must exist at a position
        width (int): figure width
        height (int): figure height
        title_type (str): type of figure being plotted
        backend (str): which plotting backend to use

    Returns:
        Union[go.Figure, Figure]: plotly or matplotlib figure object

    Raises:
        ValueError: invalid plotting backend
    """

    values_needed = len(df.index) - min_values
    df = df.loc[:, (df.isnull().sum(axis=0) <= values_needed)]

    if backend == "plotly":
        fig: go.Figure = go.Figure()
        fig.add_trace(
            go.Heatmap(
                z=df,
                x=df.columns,
                y=df.index,
                hoverongaps=False,
            ),
        )
        fig.update_traces(
            {"showscale": False, "coloraxis": None, "colorscale": "RdBu_r"},
            selector={"type": "heatmap"},
        )

        fig.update_layout(
            title_text=f"{title_type} Methylation Heatmap",
            xaxis_title_text="Position",
            yaxis_title_text="Cell Line",
            # aspect="equal",
            template="seaborn",
            autosize=False,
            width=width,
            height=height,
        )

        return fig

    if backend == "matplotlib":
        mplwidth: float = width / 66.667
        mplheight: float = height / 66.667
        df = df.apply(pd.to_numeric)
        sns.set_theme(style="white", rc={"figure.figsize": (mplwidth, mplheight)})
        ax = sns.heatmap(df, cmap="vlag", cbar=False)
        ax.set(xlabel="Position", ylabel="Cell Line")

        return ax

    raise ValueError("Invalid plotting backend")


def _create_methylation_diffs_bar_graph(
    df: pd.DataFrame, backend: str
) -> Union[go.Figure, Figure]:
    """Generate a graph figure showing bar graph of significant methylation across
    cell lines.

    Args:
        df (pd.DataFrame): dataframe of significant methylation positions
        backend (str): which plotting backend to use

    Returns:
        Union[go.Figure, Figure]: plotly or matplotlib figure object

    Raises:
        ValueError: invalid plotting backend
    """
    data = df.pivot_table(index="position", columns="cellLine", values="ad_stat")
    data = data.dropna(axis=1, how="all").count(axis=1).to_frame()

    if backend == "plotly":
        fig = go.Figure()
        fig.add_trace(go.Bar(x=data.index, y=data.values.flatten()))
        fig.update_layout(
            xaxis_type="linear",
            showlegend=False,
            title="Significant Methylation Differences",
            template="seaborn",
        )
        fig.update_xaxes(
            tickformat="r",
            tickangle=45,
            nticks=40,
            title="Position",
            range=[int(data.index.min()), int(data.index.max())],
        )
        fig.update_yaxes(title="# of significant differences")
        fig.update_traces(width=10)

        return fig
    if backend == "matplotlib":
        data = data.set_index(pd.to_numeric(data.index, errors="coerce"))
        ax = data.plot(
            kind="bar",
            legend=False,
            figsize=(12, 8),
        )
        ax.set(xlabel="Position", ylabel="# of significant differences")
        return ax

    raise ValueError("Invalid plotting backend")


def _make_stacked_fig(df: pd.DataFrame, backend: str) -> Union[go.Figure, Figure]:
    """Generate a graph figure showing methylated and unmethylated reads across
    cell lines.

    Args:
        df (pd.DataFrame): dataframe of individual read data
        backend (str): plotting backend to use

    Returns:
        Union[go.Figure, Figure]: plotly or matplotlib figure

    Raises:
        ValueError: invalid plotting backend
    """

    if backend == "plotly":
        return _make_stacked_plotly_fig(df)

    if backend == "matplotlib":
        return _make_stacked_mpl_fig(df)

    raise ValueError("Invalid plotting backend")


def _make_binary(data: Optional[List[int]]) -> List[int]:
    if isinstance(data, list):
        new = [data.count(0), data.count(1)]
    else:
        new = [0, 0]
    return new


def _make_methyl_df(df: pd.DataFrame, row: str) -> pd.DataFrame:
    new_df = pd.DataFrame(
        {each[0]: each[1].values[0] for each in df.loc[row].to_frame().iterrows()}
    )
    new_df = new_df.rename(index={0: "unmethylated", 1: "methylated"})
    return new_df


def _make_stacked_plotly_fig(df: pd.DataFrame) -> go.Figure:
    """Generate a graph figure showing methylated and unmethylated reads across
    cell lines.

    Args:
        df (pd.DataFrame): dataframe of individual read data

    Returns:
        go.Figure: plotly figure
    """

    def make_it_all(df: pd.DataFrame) -> Dict[str, go.Figure]:
        def make_methyl_bar(df: pd.DataFrame) -> go.Figure:
            fig = px.bar(
                df.T,
                height=20,
                width=800,
                color_discrete_sequence=["white", "black"],
                labels={"value": "reads", "variable": "status"},
            )
            fig.update_layout(showlegend=False, margin={"l": 0, "r": 0, "t": 0, "b": 0})
            fig.update_xaxes(visible=False)
            fig.update_yaxes(visible=False)
            return fig

        df2 = df.applymap(_make_binary)
        fig_dict = {}
        for each in df2.index:
            methyl_df = _make_methyl_df(df2, each)
            new_fig = make_methyl_bar(methyl_df)
            fig_dict[each] = new_fig
        return fig_dict

    fig_dict = make_it_all(df)

    this_figure = sp.make_subplots(
        rows=len(fig_dict),
        cols=2,
        column_widths=[0.1, 0.9],
        horizontal_spacing=0,
        vertical_spacing=0.005,
        shared_xaxes=True,
    )

    for i, (key, each) in enumerate(fig_dict.items()):
        this_figure.append_trace(each["data"][0], row=i + 1, col=2)
        this_figure.append_trace(each["data"][1], row=i + 1, col=2)
        this_figure.add_annotation(
            text=key,
            xref="x",
            yref="y",
            x=1,
            y=1,
            showarrow=False,
            font={"size": 12},
            row=i + 1,
            col=1,
        )

    this_figure.update_layout(
        showlegend=False,
        margin={"l": 0, "r": 0, "t": 40, "b": 0},
        barmode="stack",
        width=800,
        height=(len(fig_dict) * 40) + 40,
        template="ggplot2",
        paper_bgcolor="rgba(132,132,132,1)",
        title_text="Methylated and Unmethylated Reads by Cell Line",
        title={"font": {"color": "white"}},
    )
    this_figure.update_xaxes(
        showgrid=False,
        nticks=11,
        tickangle=45,
        color="white",
    )
    this_figure.update_xaxes(showticklabels=False, visible=False)
    this_figure.update_xaxes(
        showticklabels=True, tickcolor="white", visible=True, row=len(fig_dict), col=2
    )
    this_figure.update_yaxes(visible=False)

    return this_figure


def _make_stacked_mpl_fig(df: pd.DataFrame) -> Figure:
    """Generate a graph figure showing methylated and unmethylated reads across
    cell lines.

    Args:
        df (pd.DataFrame): dataframe of individual read data

    Returns:
        Figure: matplotlib figure
    """
    df2 = df.applymap(_make_binary)
    df_dict = {}
    for each in df2.index:
        methyl_df = _make_methyl_df(df2, each)
        df_dict[each] = methyl_df
    new_df = pd.DataFrame()
    for k, v in df_dict.items():
        int_df = v.T.reset_index()
        int_df["cell"] = k
        new_df = pd.concat([new_df, int_df])

    sns.set_theme(style="dark")
    fig, ax = plt.subplots(
        len(new_df["cell"].unique()),
        1,
        sharex="all",
        sharey="all",
        squeeze=False,
    )
    for i, each in enumerate(new_df["cell"].unique()):
        tmp_df = new_df[new_df["cell"] == each].loc[
            :, ["index", "unmethylated", "methylated"]
        ]
        ax[i] = tmp_df.plot(
            kind="bar",
            legend=False,
            width=1,
            stacked=True,
            color=["white", "black"],
            xticks=[],
            title="",
            yticks=[],
            xlabel="",
            ylabel="",
            figsize=(12, 0.5),
            ax=ax[i],
        )
        ax[i].set_ylabel(each, rotation=0)
        ax[i].yaxis.set_label_coords(-0.05, 0.3)
    fig.set_size_inches(12, 0.6 * len(new_df["cell"].unique()))
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.suptitle("Methylation of reads per cell line", y=1.02)
    plt.tight_layout(pad=0)

    return fig
