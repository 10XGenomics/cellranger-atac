"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

Layout handling tools for plotly plot data
"""

from __future__ import division, print_function
import numpy as np
import copy

# Subplot layout definitions
MAXIMUM_SUBPLOT_COLUMNS = 4
SUBPLOT_GRID_HEIGHT = 350
SUBPLOT_GRID_WIDTH = 380
SUBPLOT_GRID_HEIGHT_PADDING = 40
SUBPLOT_GRID_WIDTH_PADDING = 40

# Default config kwargs for plotly plots
PLOT_CONFIG_KWARGS = {
    "staticPlot": False,
    "displayModeBar": True,
    "modeBarButtons": [["toImage"]],
    "showAxisDragHandles": True,
    "scrollZoom": False,
}


def layout_samples_into_grid(num_samples, sample_index=None):
    """If sample_index is provided, give the (zero-indexed) row and column indices of the sample.  Otherwise,
    give the total number of rows and columns in the needed grid.
    """
    total_cols = min(MAXIMUM_SUBPLOT_COLUMNS, num_samples)
    total_rows = num_samples // total_cols + (1 if num_samples % total_cols else 0)
    if sample_index is None:
        return total_rows, total_cols
    if sample_index >= num_samples:
        raise ValueError('Sample index {} is larger than the total number of samples {}'.format(sample_index, num_samples))
    # Invert row indices so that subplots fill in from the top instead of the bottom
    row = sample_index // total_cols
    row = total_rows - row - 1
    col = sample_index % total_cols
    return row, col


def initialize_empty_subplot(num_samples, title=""):
    rows, cols = layout_samples_into_grid(num_samples)
    return {
        "layout": {
            "showlegend": False,
            "height": SUBPLOT_GRID_HEIGHT * rows + SUBPLOT_GRID_HEIGHT_PADDING * (rows - 1),
            "width": SUBPLOT_GRID_WIDTH * cols + SUBPLOT_GRID_WIDTH_PADDING * (cols - 1),
            "annotations": [],
            "title": title,
        },
        "data": [],
    }


def set_common_axis_range(plot_data, log_zero_sub_vals=None):
    """ Sets a common x and y range for a series of plotly charts contained in plot_data.

        log_zero_sub_vals is a dictionary {"x": xsub, "y": ysub}
        and specifies how we handle log(0) if encountered in a log axis.
        For example, if the x value of any data point is 0, the min value on
        a log scale is technically -infinity, but instead we set it to 
        log(xsub).
    """
    bounds = {
        "x": {"min": np.inf, "max": -np.inf},
        "y": {"min": np.inf, "max": -np.inf}
    }
    ops = {
        "max": max,
        "min": min
    }
    if log_zero_sub_vals is None:
        log_zero_sub_vals = {
            "x": 1,
            "y": 1,
        }

    for dataset in plot_data["data"]:
        for axis in ["x", "y"]:
            axis_type = plot_data["layout"]["{}axis".format(axis)].get("type", "-")
            for op_type, op in ops.items():
                local_bound = op(dataset[axis])
                if axis_type == "log":
                    local_bound = np.log10(local_bound)
                    if np.isinf(local_bound):
                        local_bound = np.log10(log_zero_sub_vals[axis])
                bounds[axis][op_type] = op(local_bound, bounds[axis][op_type])

    for i in xrange(len(plot_data["data"])):
        axis_id = "" if i == 0 else str(i + 1)
        for axis in ["x", "y"]:
            key = "{}axis{}".format(axis, axis_id)
            if key not in plot_data["layout"]:
                continue
            plot_data["layout"][key]["range"] = [bounds[axis]["min"], bounds[axis]["max"]]


def place_chart_in_grid(plot_data, subplot_data, num_subplots, subplot_index, subplot_title,
                        xlabel=None, ylabel=None, xscale=None, yscale=None):
    """ Given a plotly chart ("subplot_data") defined on a single axis, place it
        correctly on a grid of subplots ("plot_data"). The other parameters are
        xlabel = x axis label
        ylabel = y axis label
        subplot_title = subplot title for the chart in the grid
        xscale = any plotly axis type (e.g. log)
        yscale = any plotly axis type (e.g. log)
    """
    total_rows, total_cols = layout_samples_into_grid(num_subplots)
    row, col = layout_samples_into_grid(num_subplots, subplot_index)

    x_padding = SUBPLOT_GRID_WIDTH_PADDING / (SUBPLOT_GRID_WIDTH * total_cols + SUBPLOT_GRID_WIDTH_PADDING * (total_cols - 1)) if total_cols > 1 else 0.
    y_padding = SUBPLOT_GRID_HEIGHT_PADDING / (SUBPLOT_GRID_HEIGHT * total_rows + SUBPLOT_GRID_HEIGHT_PADDING * (total_rows - 1)) if total_rows > 1 else 0.

    axis_id = "" if subplot_index == 0 else str(subplot_index + 1)

    extents = {"y": np.linspace(0, 1, total_rows + 1),
               "x": np.linspace(0, 1, total_cols + 1)}

    layout = plot_data["layout"]

    if any(("type" in item and item["type"] == "sankey") for item in subplot_data["data"]):
        # Sankey diagrams have to be handled specially to get layout and display working OK
        for item in subplot_data["data"]:
            item["domain"] = {
                'x': [extents["x"][col] + x_padding, extents["x"][col + 1] - x_padding],
                'y': [extents["y"][row] + y_padding, extents["y"][row + 1] - y_padding],
            }
            item["node"]["pad"] = 10
            item["node"]["thickness"] = 20
            item["node"]["label"] = [l.replace(" Fragments", "").replace(" Barcodes", "")
                                     for l in item["node"]["label"]]
    else:
        # Other plots need separate axes defined
        layout["xaxis{}".format(axis_id)] = copy.deepcopy(subplot_data["layout"]["xaxis"])
        layout["xaxis{}".format(axis_id)].update({
            "domain": [extents["x"][col] + x_padding, extents["x"][col + 1] - x_padding],
            "anchor": "y{}".format(axis_id),
            "showline": True,
            "zeroline": False,
        })
        xlabel = xlabel if xlabel is not None else subplot_data["layout"]["xaxis"].get("title")
        xscale = xscale if xscale is not None else subplot_data["layout"]["xaxis"].get("type")
        if xlabel is not None:
            layout["xaxis{}".format(axis_id)]["title"] = xlabel
        if xscale is not None:
            layout["xaxis{}".format(axis_id)]["type"] = xscale

        layout["yaxis{}".format(axis_id)] = copy.deepcopy(subplot_data["layout"]["yaxis"])
        layout["yaxis{}".format(axis_id)].update({
            "domain": [extents["y"][row] + y_padding, extents["y"][row + 1] - y_padding],
            "anchor": "x{}".format(axis_id),
            "showline": True,
            "zeroline": False,
        })
        ylabel = ylabel if ylabel is not None else subplot_data["layout"]["yaxis"].get("title")
        yscale = yscale if yscale is not None else subplot_data["layout"]["yaxis"].get("type")
        if ylabel is not None:
            layout["yaxis{}".format(axis_id)]["title"] = ylabel
        if yscale is not None:
            layout["yaxis{}".format(axis_id)]["type"] = yscale

        # Attach the subplot's datasets to the correct axes
        for axis in ["x", "y"]:
            for d in subplot_data["data"]:
                d["{}axis".format(axis)] = "{}{}".format(axis, axis_id)

    plot_data["data"].extend(subplot_data["data"])

    layout['annotations'].append({
        "text": subplot_title,
        "x": (col + 0.5) / total_cols,
        "y": (1 + row) / total_rows - y_padding / 2,
        "xref": "paper",
        "yref": "paper",
        "xanchor": "center",
        "yanchor": "top",
        "font": {"size": 14},
        "visible": True,
        "showarrow": False,
    })

