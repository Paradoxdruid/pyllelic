#!/usr/bin/env python3
"""Utilities to visualize data for use in pyllelic."""

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import plotly.subplots as sp


def _create_histogram(data: pd.DataFrame, cell_line: str, position: str) -> go.Figure:
    """Generate a graph figure showing fractional methylation in
    a given cell line at a given site.

    Args:
        data (pd.DataFrame): dataframe of individual data
        cell_line (str): name of cell line
        position (str): genomic position

    Returns:
        go.Figure: plotly figure object
    """
    fig: go.Figure = go.Figure()
    fig.add_trace(
        go.Histogram(
            x=data.loc[cell_line, position],
            xbins=dict(
                start=-0.1,
                end=1.1,
                size=0.2,
            ),  # offset bins to center displayed bars
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


def _create_heatmap(
    df: pd.DataFrame, min_values: int, width: int, height: int, title_type: str
) -> go.Figure:
    """Generate a graph figure showing heatmap of mean methylation across
    cell lines.

    Args:
        df (pd.DataFrame): dataframe of mean methylation
        min_values (int): minimum number of points data must exist at a position
        width (int): figure width
        height (int): figure height
        title_type (str): type of figure being plotted

    Returns:
        go.Figure: plotly figure object
    """

    values_needed = len(df.index) - min_values
    df = df.loc[:, (df.isnull().sum(axis=0) <= values_needed)]

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
        dict(showscale=False, coloraxis=None, colorscale="RdBu_r"),
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


def _create_methylation_diffs_bar_graph(df: pd.DataFrame) -> go.Figure:
    """Generate a graph figure showing bar graph of significant methylation across
    cell lines.

    Args:
        df (pd.DataFrame): dataframe of significant methylation positions

    Returns:
        go.Figure: plotly figure object
    """
    data = df.pivot(index="position", columns="cellLine", values="ad_stat")
    data = data.dropna(axis=1, how="all").count(axis=1).to_frame()

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
    fig.update_traces(width=50)

    return fig


def _make_stacked_fig(df: pd.DataFrame) -> go.Figure:
    """Generate a graph figure showing methylated and unmethylated reads across
    cell lines.

    Args:
        df (pd.DataFrame): dataframe of individual read data

    Returns:
        go.Figure: plotly figure option
    """

    def make_it_all(df):
        def make_binary(data):
            if isinstance(data, list):
                new = [data.count(0), data.count(1)]
            else:
                new = [0, 0]
            return new

        def make_methyl_df(df, row):
            new_df = pd.DataFrame(
                {
                    each[0]: each[1].values[0]
                    for each in df.loc[row].to_frame().iterrows()
                }
            )
            new_df = new_df.rename(index={0: "unmethylated", 1: "methylated"})
            return new_df

        def make_methyl_bar(df):
            fig = px.bar(
                df.T,
                height=20,
                width=800,
                color_discrete_sequence=["white", "black"],
                labels={"value": "reads", "variable": "status"},
            )
            fig.update_layout(showlegend=False, margin=dict(l=0, r=0, t=0, b=0))
            fig.update_xaxes(visible=False)
            fig.update_yaxes(visible=False)
            return fig

        df2 = df.applymap(make_binary)
        fig_dict = {}
        for each in df2.index:
            methyl_df = make_methyl_df(df2, each)
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
            font=dict(
                size=12,
            ),
            row=i + 1,
            col=1,
        )

    this_figure.update_layout(
        showlegend=False,
        margin=dict(l=0, r=0, t=40, b=0),
        barmode="stack",
        width=800,
        height=len(fig_dict) * 40,
        template="ggplot2",
        paper_bgcolor="rgba(132,132,132,1)",
        title_text="Methylated and Unmethylated Reads by Cell Line",
        title={"font": {"color": "white"}},
    )
    this_figure.update_xaxes(visible=False)
    this_figure.update_yaxes(visible=False)

    return this_figure
