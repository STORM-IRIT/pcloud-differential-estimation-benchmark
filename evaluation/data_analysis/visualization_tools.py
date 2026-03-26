import plotly.graph_objects as go
import plotly.offline as pyo
import plotly.subplots as sp
import numpy as np 
import random
import math

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import methods as m
import sys
import os
import re

random_color = ['rgba(243,101,81,1)', 'rgba(131,22,23,1)', 'rgba(96,128,64,1)', 'rgba(171,157,84,1)', 'rgba(171,117,84,1)', 'rgba(81,100,145,1)', 'rgba(81,129,242,1)', 'rgba(157,45,243,1)', 'rgba(107,0,243,1)', 'rgba(18,100,235,1)','rgba(18,230,7,1)']

properties = {
    "Mean curvature (mean)": m.mean_curvature_estimation,
    "Gaussian curvature (mean)": m.principal_curvature_estimation,
    "K1 mean": m.principal_curvature_estimation,
    "K2 mean": m.principal_curvature_estimation,
    "Timings (mean)": m.all_methods,
    "D1 mean": m.curvature_direction,
    "D2 mean": m.curvature_direction,
    "normal mean": m.normal_estimation,
    "iShape mean": m.principal_curvature_estimation,
    "pos mean": m.all_methods,
    "Nb neighbors (mean)" : m.all_methods,
    "non_stable_ratio" : m.all_methods,
}

def generate_random_color(alpha='1'):
    '''Generate a random color in RGBA format.'''
    r, g, b = random.sample(range(0, 256), 3)
    return f'rgba({r},{g},{b},{alpha})'


dash_patterns = [
    'solid',
    'dot',
    'dash',
    'dashdot',
    'longdashdot',
    '5px,10px',  # 5px dash, 10px space
    '10px,5px,2px,5px',  # More complex pattern
    'longdash',
]

def transpose_subplots(subplots_dict):
    subplots = subplots_dict['data']

    for plot in subplots :
        plot['row'], plot['col'] = plot['col'], plot['row']
    
    subplots_dict = {
        'data': subplots,
        'transpose': True,
    }
    return subplots_dict

def create_subplots(dataset, rows, cols, to_legend, colors, line_width=1.5, split_by=None):
    """
    Rows are configured as : 
    [
        ( "constraint", [the constrained row's param], "x", "x_label")
    ]
    Cols are configured as : 
    [
        ( "constraint", [the constrained col's param], "y", "y_label", "log")
    ]
    """

    subplot_data = []

    displayed_legend = set()

    for i, row in enumerate(rows) :
        data_row = dataset[dataset[row[0]].isin(row[1])]
        for j, col in enumerate(cols) :
            data_col = data_row[data_row[col[0]].isin(col[1])]
            traces = []
            for legend in data_col[to_legend].unique() :
                data_legend = data_col[data_col[to_legend] == legend]
                
                showlegend = legend not in displayed_legend
                if showlegend :
                    displayed_legend.add(legend)
                
                x_datas = []
                y_datas = []
                lines = []
                line_names = []
                
                if split_by is not None :
                    
                    is_unique_split = len(data_legend[split_by].unique()) == 1

                    for s_idx, split in enumerate( data_legend[split_by].unique() ) :
                        data_split = data_legend[data_legend[split_by] == split]
                        x_data = sorted(data_split[row[2]].unique())
                        y_data = data_split[col[2]].groupby(data_split[row[2]]).mean()
                        x_datas.append(x_data)
                        y_datas.append(y_data)
                        lines.append(dict(width=line_width*1.4, color=colors[legend], dash=dash_patterns[s_idx % len(dash_patterns)]) if (s_idx > 0) else dict(color=colors[legend]))
                        if is_unique_split :
                            line_names.append(f"{legend}")
                        else :
                            line_names.append(f"{legend}_{split}")
                else :
                    x_datas = [sorted(data_legend[row[2]].unique())]
                    y_datas = [data_legend[col[2]].groupby(data_legend[row[2]]).mean()]
                    lines = [dict(color=colors[legend], width=line_width)]
                    line_names = [legend]

                for x_data, y_data, line, name in zip(x_datas, y_datas, lines, line_names) :
                    traces.append(
                        dict(
                            x=x_data, 
                            y=y_data, 
                            mode='lines', 
                            name=name, 
                            legendgroup=name,
                            showlegend=showlegend,
                            line=line,
                            marker=dict(color=colors[legend])
                        )
                    )

            subplot_data.append(
            {
                'title': "",
                'traces': traces,
                'x_title': row[3],
                'y_title': col[3],
                'x_basename' : row[2],
                'y_basename' : col[2],
                'row': i+1,
                'col': j+1, 
                'y_axis_type': "linear" if len(col) < 5 else col[4],
            })
    subplot_data = {
        'data': subplot_data,
        'transpose': False,
    }
    return subplot_data


def is_shared_axis(subplot_data, max_rows, max_cols):
    x_axis = [
        "" for _ in range(max_rows)
    ]
    y_axis = [
        "" for _ in range(max_cols)
    ]
    for data in subplot_data :
        x_axis[data['row']-1] = data['x_basename']
        y_axis[data['col']-1] = data['y_basename']
    shared_x = True
    shared_y = True
    for i in range(max_rows-1) :
        if x_axis[i] != x_axis[i+1] :
            shared_x = False
            break
    for i in range(max_cols-1) :
        if y_axis[i] != y_axis[i+1] :
            shared_y = False
            break
    return shared_x, shared_y

base_fonts = {
    'title_font': 18,
    'axis_title_font': 15,
    'axis_tick_font': 10,
    'annotation_font': 15,
    'legend_font' : 12,
    'font_factor': 1.35,
    'font_family': 'Times New Roman',
}

def get_legend_height (list_legend, width, fonts, font_heights, show_titles):
    char_width = fonts['legend_font'] * 0.6
    char_height = font_heights['legend_height']
    bar_width = 30
    spacing = 30 if show_titles else 0
    total_length = 0
    for legend in list_legend :
        total_length += ( len(legend) * char_width ) + spacing + bar_width
    return ( ( math.ceil(total_length / width) ) * char_height ) + bar_width * 2 

def get_font_height (fonts, showtitles) :
    font_heights = {
        'title_height': fonts['title_font'] * fonts['font_factor'] if showtitles else fonts['annotation_font'] * fonts['font_factor'],
        'axis_title_height': fonts['axis_title_font'] * fonts['font_factor'] if showtitles else 0,
        'axis_tick_height': fonts['axis_tick_font'] * fonts['font_factor'],
        'annotation_height': fonts['annotation_font'] * fonts['font_factor'],
        'legend_height': fonts['legend_font'] * fonts['font_factor'] if showtitles else 0,
    }
    return font_heights

def find_titles (subplot_data, max_rows, max_cols, show_titles) :
    x_names = [ ["" for _ in range(max_cols)] for _ in range(max_rows) ]
    y_names = [ ["" for _ in range(max_cols)] for _ in range(max_rows) ]

    titles = []
    for data in subplot_data:
        row, col = data.get('row', 1), data.get('col', 1)
        x_names[row-1][col-1] = data['x_title'] if show_titles else ""
        y_names[row-1][col-1] = data['y_title'] if show_titles else ""
        
        if show_titles:
            titles.append(data['title'])
        else : 
            if data.get('row', 1) == 1 :
                titles.append(data['y_title'])
            else :
                titles.append("")

    return x_names, y_names, titles

def find_titles_transposed (subplot_data, max_rows, max_cols, show_titles) :
    x_names = [ ["" for _ in range(max_cols)] for _ in range(max_rows) ]
    y_names = [ ["" for _ in range(max_cols)] for _ in range(max_rows) ]

    titles = []
    for data in subplot_data:
        row, col = data.get('row', 1), data.get('col', 1)
        x_names[row-1][col-1] = data['x_title'] if show_titles else ""
        y_names[row-1][col-1] = data['y_title'] if show_titles else ""
        
        if show_titles:
            titles.append(data['title'])
        else : 
            if data.get('col', 1) == 1 :
                y_names[row-1][col-1] = data['y_title']
            else :
                titles.append("")

    return x_names, y_names, titles

def simplify_axis_name (fig, x_names, y_names, max_rows, max_cols, collapse_x, collapse_y) :
    for col in range(max_cols) :
        same = False
        for row in range(max_rows) :
            if x_names[row][col] == x_names[0][col] :
                same = True
            else :
                same = False
                break
        if same and collapse_x :
            for row in range(max_rows-1) :
                fig.update_xaxes(title_text="", row=row+1, col=col+1)

    for row in range(max_rows) :
        same = False
        for col in range(max_cols) :
            if y_names[row][col] == y_names[row][0] :
                same = True
            else :
                same = False
                break
        if same and collapse_y :
            for col in range(1, max_cols) :
                fig.update_yaxes(title_text="", row=row+1, col=col+1)
    return fig

def fix_log_ticks(fig, y_data, row=1, col=1, show_minor=True, max_ticks=6, min_ticks=2):
    """
    Correct the log ticks to ensure major ticks are always visible, 
    by adjusting the range to include the upper and lower ticks.
    """
    y_valid = np.array([v for v in y_data if v is not None and np.isfinite(v) and v > 0])
    if len(y_valid) == 0:
        return fig

    y_min, y_max = y_valid.min(), y_valid.max()

    pmin = int(np.floor(np.log10(y_min)))
    pmax = int(np.ceil(np.log10(y_max)))
    if pmax == pmin:
        pmin -= 1

    p_range = pmax - pmin
    step = max(1, int(np.ceil(p_range / max_ticks)))
    major_exponents = list(range(pmin, pmax + 1, step))
    if len(major_exponents) < min_ticks:
        major_exponents = [pmin, pmax]

    major_ticks = [10 ** i for i in major_exponents]
    major_labels = [f"10<sup>{i}</sup>" for i in major_exponents]

    # Ticks mineurs
    minor_ticks = []
    if show_minor and p_range <= 3:
        for i in range(pmin, pmax):
            minor_ticks.extend([m * 10 ** i for m in [2, 5]])

    axis_min = y_min
    axis_max = y_max

    # if len(major_ticks) > 1 and p_range <= 3:

    #     # Distances logarithmiques
    #     log_y_min, log_y_max = np.log10(y_min), np.log10(y_max)
    #     log_major_ticks = np.log10(major_ticks)

    #     tick_below_min = max([t for t in major_ticks if t <= y_min], default=y_min)
    #     tick_above_max = min([t for t in major_ticks if t >= y_max], default=y_max)

    #     dist_min_below = log_y_min - np.log10(tick_below_min)
    #     dist_max_above = np.log10(tick_above_max) - log_y_max

    #     if dist_min_below < dist_max_above:
    #         axis_min = tick_below_min  # étendre vers le tick inférieur
    #         axis_max = y_max
    #     else:
    #         axis_min = y_min
    #         axis_max = tick_above_max  # étendre vers le tick supérieur


    # axis_max = max(y_max, major_ticks[-1]) if p_range <= 3 else y_max

    fig.update_yaxes(
        row=row, col=col,
        type="log",
        tickvals=major_ticks + minor_ticks,
        ticktext=major_labels + [""] * len(minor_ticks),
        
        showline=True,
        mirror=False,
        linewidth=1.75,
        linecolor='rgb(115,115,115)',
        range=[np.log10(axis_min), np.log10(axis_max)]
    )

    return fig

def create_modular_subplots(
    subplots, 
    title=None, 
    page_width=718.75, 
    graph_ratio=0.8, 
    fonts=base_fonts, 
    margin_percent=4, 
    show_titles=True, 
    collapse_x=False, 
    collapse_y=False, 
    show_legend=True,
    row_ylabels=None
):
    isTransposed = subplots["transpose"]
    subplot_data = subplots["data"]
    font_heights = get_font_height (fonts, show_titles)
    max_rows = max(data.get('row', 1) for data in subplot_data)
    max_cols = max(data.get('col', 1) for data in subplot_data)

    margin_px = page_width * (margin_percent / 100)
    available_width = page_width - 2 * margin_px
    column_widths = [1/max_cols] * max_cols
    subgraph_height = (available_width / max_cols) * graph_ratio
    total_graph_height = subgraph_height * max_rows

    subgraph_height = (available_width / max_cols) * graph_ratio
    subplot_px_height = total_graph_height / max_rows

    max_ticks_allowed = max(3, int(subplot_px_height / 20))

    total_height = (
        2 * margin_px +
        font_heights['title_height'] +
        total_graph_height +
        ( ( font_heights['axis_title_height'] + font_heights['axis_tick_height'] ) * ( max_rows + 1) )
    )
    row_heights = [ ( subgraph_height ) / total_height] * max_rows
    
    x_names, y_names, titles = find_titles(subplot_data, max_rows, max_cols, show_titles)
    if isTransposed :
        x_names, y_names, titles = find_titles_transposed(subplot_data, max_rows, max_cols, show_titles)
    
    is_shared_x, is_shared_y = is_shared_axis(subplot_data, max_rows, max_cols)

    fig = sp.make_subplots( rows=max_rows, cols=max_cols,
        column_widths=column_widths, row_heights=row_heights,
        shared_xaxes=is_shared_x & collapse_x,
        shared_yaxes=is_shared_y & collapse_y,
        horizontal_spacing=0.08 if show_titles else 0.05, 
        vertical_spacing=0.05 if show_titles else 0.05,
        subplot_titles=titles,
    )
    legend_names = []
    
    for data in subplot_data:
        row, col = data.get('row', 1), data.get('col', 1)
        y_values = []
            
        for trace in data['traces']:
            fig.add_trace(go.Scatter(**trace), row=row, col=col)
            
            if trace['name'] not in legend_names:
                legend_names.append(trace['name'])

            if 'y' in trace:
                y_values.extend(trace['y'])


        fig.update_xaxes(title_text=x_names[row-1][col-1], row=row, col=col,
                        title_font=dict(size=fonts['axis_title_font']),
                        gridcolor="lightgray",
                        gridwidth=1.25,
                        showgrid=True,
                        
                        showline=True,
                        linewidth=1.75,
                        linecolor='rgb(115,115,115)',
                        mirror=False,
                        
                        tickfont=dict(size=fonts['axis_tick_font']))
        
        fig.update_yaxes(title_text=y_names[row-1][col-1], row=row, col=col,
            type=data.get('y_axis_type', 'linear'),
            title_font=dict(size=fonts['axis_title_font']),
            gridcolor="lightgray",
            gridwidth=1.25,
            showgrid=True,
            
            showline=True,
            linewidth=1.75,
            linecolor='rgb(115,115,115)',
            mirror=False,

            tickfont=dict(size=fonts['axis_tick_font']),
        )
                
        for i in fig['layout']['annotations']:
            i['font'] = dict(size=fonts['annotation_font'], family=fonts['font_family'])
        
        if data.get('y_axis_type','linear') == 'log' and len(y_values) > 0:
                fig = fix_log_ticks(
                    fig, y_values, row=row, col=col, 
                    show_minor=True, 
                    max_ticks=max_ticks_allowed
                )

        if row_ylabels is not None and row_ylabels[row-1]:
                if col == 1:
                    fig.update_yaxes(title_text=row_ylabels[row-1], row=row, col=col,
                                        title_font=dict(size=fonts['axis_title_font']))
                else:
                    fig.update_yaxes(title_text="", row=row, col=col)

    fig = simplify_axis_name (fig, x_names, y_names, max_rows, max_cols, collapse_x, collapse_y)

    legend_height = get_legend_height (legend_names, available_width, fonts, font_heights, show_titles)
    total_height += legend_height

    fig.update_layout(
        font=dict(family=fonts['font_family'], color='black'),
        plot_bgcolor='white',
        
        title = dict( text = title if show_titles else "", font = dict(size=fonts['title_font'], family=fonts['font_family']) ),
        width=page_width,
        height=total_height,
        showlegend=show_legend,
        legend=dict(
            orientation='h',
            itemsizing='constant',
            y= - ( font_heights['axis_title_height'] + font_heights['axis_tick_height'] ) / ( total_graph_height ),
            x=0,
            font=dict(size=fonts['legend_font'], family=fonts['font_family']),
        ),
        margin=dict(l=margin_px, r=margin_px, t= ( show_titles * margin_px ) + font_heights['title_height'], b=margin_px + legend_height),
    )

    return fig