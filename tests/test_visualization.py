#!/usr/bin/env python3
"""pytest unit tests for visualization."""

# Testing
import pytest  # noqa  # pylint: disable=unused-import

# Required libraries for test data
import pandas as pd
from .inputs import EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA

# Module to test
import pyllelic.visualization as viz  # noqa  # pylint: disable=unused-import


TEST_INDIVIDUAL_DATA = pd.DataFrame.from_dict(
    {
        "cellLine": ["1", "2"],
        "position": ["1", "2"],
        "ad_stat": [0, 1],
        "p_crit": [1, 2],
        "diff": [1, 2],
        "raw": [1, 2],
    }
).astype("object")


def test__create_heatmap(mocker):
    mocked_go = mocker.patch("pyllelic.visualization.go")

    _ = viz._create_heatmap(
        TEST_INDIVIDUAL_DATA,
        min_values=1,
        height=600,
        width=600,
        title_type="means",
        backend="plotly",
    )

    mocked_go.Figure.assert_called_once()
    mocked_go.Heatmap.assert_called_once()


def test__create_heatmap_mpl(mocker):
    mocked_mpl = mocker.patch("pyllelic.visualization.sns")

    _ = viz._create_heatmap(
        TEST_INDIVIDUAL_DATA,
        min_values=1,
        height=600,
        width=600,
        title_type="means",
        backend="matplotlib",
    )

    mocked_mpl.heatmap.assert_called_once()


def test__create_heatmap_backend_error():
    with pytest.raises(ValueError):
        _ = viz._create_heatmap(
            TEST_INDIVIDUAL_DATA,
            min_values=1,
            height=600,
            width=600,
            title_type="means",
            backend="other",
        )


def test__create_histogram(mocker):
    mocked_go = mocker.patch("pyllelic.visualization.go")
    TEST_CELL_LINE = "TEST1"
    TEST_POSITION = "1"
    intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
    TEST_DATA = intermediate.astype("object")

    _ = viz._create_histogram(
        TEST_DATA, TEST_CELL_LINE, TEST_POSITION, backend="plotly"
    )

    mocked_go.Figure.assert_called_once()
    mocked_go.Histogram.assert_called_once_with(
        x=TEST_DATA.loc[TEST_CELL_LINE, TEST_POSITION],
        xbins=dict(
            start=-0.1,
            end=1.1,
            size=0.2,
        ),
    )


def test__create_histogram_mpl(mocker):
    mocked_mpl = mocker.patch("pyllelic.visualization.sns")
    TEST_CELL_LINE = "TEST1"
    TEST_POSITION = "1"
    intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
    TEST_DATA = intermediate.astype("object")

    _ = viz._create_histogram(
        TEST_DATA, TEST_CELL_LINE, TEST_POSITION, backend="matplotlib"
    )

    mocked_mpl.histplot.assert_called_once()


def test__create_histogram_backend_error():
    TEST_CELL_LINE = "TEST1"
    TEST_POSITION = "1"
    intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
    TEST_DATA = intermediate.astype("object")

    with pytest.raises(ValueError):
        _ = viz._create_histogram(
            TEST_DATA, TEST_CELL_LINE, TEST_POSITION, backend="other"
        )


def test__create_methylation_diffs_bar_graph(mocker):
    mocked_go = mocker.patch("pyllelic.visualization.go")

    _ = viz._create_methylation_diffs_bar_graph(TEST_INDIVIDUAL_DATA, backend="plotly")

    mocked_go.Figure.assert_called_once()
    mocked_go.Bar.assert_called_once()


def test__create_methylation_diffs_bar_graph_mpl(mocker):
    mocked_mpl = mocker.patch("pyllelic.visualization.pd.DataFrame.plot")

    _ = viz._create_methylation_diffs_bar_graph(
        TEST_INDIVIDUAL_DATA, backend="matplotlib"
    )

    mocked_mpl.assert_called_once()


def test__create_methylation_diffs_bar_graph_invalid():
    with pytest.raises(ValueError):
        _ = viz._create_methylation_diffs_bar_graph(
            TEST_INDIVIDUAL_DATA, backend="FAKE"
        )


def test__make_stacked_fig(mocker):
    mocked_px = mocker.patch("pyllelic.visualization.px")
    mocked_sp = mocker.patch("pyllelic.visualization.sp")

    intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
    TEST_DATA = intermediate.astype("object")

    _ = viz._make_stacked_fig(TEST_DATA, backend="plotly")

    mocked_px.bar.assert_called()
    mocked_sp.make_subplots.assert_called_once()


def test__make_stacked_fig_mpl(mocker):
    mocked_df_plot = mocker.patch("pyllelic.visualization.pd.DataFrame.plot")
    mocked_mpl = mocker.patch("pyllelic.visualization.plt")
    mocked_mpl.subplots.return_value = (mocker.MagicMock(), mocker.MagicMock())

    intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
    TEST_DATA = intermediate.astype("object")

    _ = viz._make_stacked_fig(TEST_DATA, backend="matplotlib")

    mocked_df_plot.assert_called()
    mocked_mpl.subplots.assert_called()


def test__make_stacked_fig_invalid():
    intermediate = EXPECTED_INTERMEDIATE_INDIVIDUAL_DATA
    TEST_DATA = intermediate.astype("object")
    with pytest.raises(ValueError):
        _ = viz._make_stacked_fig(TEST_DATA, backend="FAKE")
