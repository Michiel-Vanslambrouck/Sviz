# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 17:24:26 2019
Python 3.6
@author: mcmic
Requirements:   this file must be run in a folder that also contains all output folders from SNPeffects
                AND the program MutPosFilter.py must have ran to generate the correct mutation positions
                any external reports are expected in a subfolder called External_data, under which all dataset
                names have a folder in which there is a pph2-short.txt and dataset_provean.tsv
"""

import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from pathlib import Path
import os
import numpy as np
import dash_table
from dash.dependencies import Input, Output
import warnings

print("////////////// Sviz \\\\\\\\\\\\\\\\\\\\\\\\\\\\")


# ----------------------------------- Data preprocessing -----------------------------------------------

warnings.simplefilter("ignore")


def find_files():  # goes through all folders in the same dir and returns a list of finalreports
    report = {}
    print("Auto detecting files")
    dir_name = os.path.dirname(os.path.realpath(__file__))
    path_list = Path(dir_name).glob('**/finalreport2.tab')
    for path in path_list:
        folder_name = os.path.split(os.path.split(str(path))[0])[1]
        report[folder_name] = pd.read_csv(str(path), sep='\t')
    return report, dir_name


def preprocess(df_in, folder):  # reformat columns, returns datasets after duplicate pdb_file+mutationstring removal
    df_in.columns = \
        df_in.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
    cols = df_in.columns.tolist()
    cols = cols[:4] + cols[-3:] + cols[4:-3]
    df_in = df_in[cols]
    df_out = df_in.sort_values('alignment_score', ascending=False).drop_duplicates(['pdb_file', 'mutationstring'],
                                                                                   keep='first').sort_index()
    print("Report found in <" + folder + "> preprocessing... Lines in: " + str(df_in.variant.count()) + ". Lines out: "
          + str(df_out.variant.count()) + ".")
    df_out = df_out.round(5)
    return df_out


def fixtotalenergy(df_in):  # return true total energy after dividing by the amount of sites in the mutationstring
    # THIS STEP IS PRESENT TO FIX AN ERROR IN THE PIPELINE AND WILL LIKELY NEED TO BE REMOVED
    amount_mutations = [(str(df_in['mutationstring'][i])).count(",") + 1 for i in df_in.index.values]
    df_in.total_energy /= amount_mutations


def addexternaldata(df_in, name):  # to add sift and polyfen data, will check if present
    df_out = df_in
    if os.path.exists(headfolder + "/External_data/" + name + "/" + "pph2-short.txt"):
        polyphen_df = pd.read_csv(headfolder + "/External_data/" + name + "/" + "pph2-short.txt",
                                  sep='\t',  skipinitialspace=True)
        polyphen_df.columns = [str.strip(i) for i in polyphen_df.columns]
        polyphen_min = polyphen_df[['#o_acc', 'pos', 'pph2_prob', 'prediction']]
        polyphen_min.columns = ['uniprot_id', 'abs_mutpos', 'pph2_probability', 'pph2_prediction']
        for i in ['uniprot_id', 'pph2_prediction']:
            polyphen_min[i] = polyphen_min[i].str.strip()
        df_out = pd.merge(df_out, polyphen_min, how='left',
                          left_on=['uniprot_id', 'abs_mutpos'], right_on=['uniprot_id', 'abs_mutpos'])
        df_out = df_out.round(5)
        print('   Polyphen integrated, new size ' + str(df_out.variant.count()))
    if os.path.exists(headfolder + "/External_data/" + name + "/" + name + "_provean.tsv"):
        sift_df = pd.read_csv(headfolder + "/External_data/" + name + "/" + name + "_provean.tsv", sep='\t')
        sift_df.columns = \
            sift_df.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
        sift_min = sift_df[
            ['protein_id', 'position', 'score', 'prediction_cutoff=-2.5', 'score.1', 'prediction_cutoff=0.05']]
        sift_min.columns = \
            ['uniprot_id', 'abs_mutpos', 'provean_score', 'provean_prediction', 'sift_score', 'sift_prediction']
        df_out = pd.merge(df_out, sift_min, how='left',
                          left_on=['uniprot_id', 'abs_mutpos'], right_on=['uniprot_id', 'abs_mutpos'])
        df_out = df_out.round(5)
        print('   Sift integrated, new size: ' + str(df_out.variant.count()))
    return df_out


# actually performing the functions for preprocessing
reports, headfolder = find_files()
for x in reports:
    reports[x] = preprocess(reports[x], x)
    fixtotalenergy(reports[x])
    reports[x] = addexternaldata(reports[x], x)
    reports[x].to_csv('FinalFinalReport/' + x + '.csv')

# ----------------------------------- DASH plots processing -----------------------------------------------

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']  # from dash
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

L = []  # this is a list for adding massplots, to later sort into an html container
for r in reports:
    rep = reports[r]
    L.append(
        html.Div([
            html.H4(children=r, style={'textAlign': 'center'}),
            dcc.Graph(
                id='MASS plot' + r,
                figure={
                    'data': [
                        go.Scatter(
                            ids=rep['variant'],
                            x=rep['delta_tango_score'],
                            y=rep['total_energy'],
                            text=[rep['variant'].tolist()[i] + "<br>" + rep['gene_name'].tolist()[i]
                                  + "<br>" + rep['pdb_file'].tolist()[i] + "<br>" + rep['mutationstring'].tolist()[i]
                                  for i in range(len(rep['variant'].tolist()))],
                            mode='markers',
                            opacity=0.9,
                            marker=dict(
                                size=10,
                                color=rep['delta_waltz_score'],  # set color equal to a variable
                                showscale=True,
                                colorbar=dict(title="dWALTZ", tickvals=[-200, -100, 0, 100, 200], ticks='outside'),
                                symbol=list(map(lambda x: 3 if x == '+' else 0, rep['mutation_in_tm'])),
                                colorscale=[[0, 'rgba(140,40,220, 0.7)'],
                                            [0.25, 'rgba(30,144,255, 0.7)'],
                                            [0.5, 'rgba(50,205,50, 0.7)'],
                                            [0.65, 'rgba(250,220,1, 0.7)'],
                                            [0.85, 'rgba(255,140,0, 0.7)'],
                                            [1, 'rgba(178,34,34, 0.7)']],
                                cmin=-200, cmax=200
                            ),
                        )
                    ],
                    'layout': go.Layout(
                        xaxis={'title': 'dTANGO'},
                        yaxis={'title': 'ddG FoldX'},
                        margin={'l': 100, 'b': 100, 't': 20, 'r': 100},
                        paper_bgcolor='rgba(0,0,0,0)',
                        plot_bgcolor='#eeeefa',
                        showlegend=False,
                        hovermode='closest',
                    )
                }
            )
        ], style={'backgroundColor': '#fafafe'})
    )
L_half = int(np.round(len(L) / 2))  # want to know half length of list to sort into two columns


def differential(df1, df2, colx):  # return the intersection of two datasets
    col1 = reports[df1][['gene_name', 'mutationstring', colx, 'variant']]
    col2 = reports[df2][['gene_name', 'mutationstring', colx, 'variant']]
    # an inner join ensures that we get the genes present in both sets and each ones values
    merge_df = pd.merge(col1, col2, how='inner', left_on=['gene_name'], right_on=['gene_name'])
    out_df = merge_df[['gene_name', 'mutationstring_x', 'mutationstring_y',
                       colx + '_x', colx + '_y', 'variant_x', 'variant_y']]
    # we are not interested in nonexitent differences between identical mutations
    # nonsense_rows = out_df[out_df.mutationstring_x == out_df.mutationstring_y]
    # out_df = out_df.drop(nonsense_rows.index)
    # out_df.to_csv('testo.csv')
    return out_df


# ----------------------------------- DASH html structure -----------------------------------------------

app.layout = \
    html.Div(children=[  # light gray page border

        # Page head
        html.Div([  # title bar
            html.H1('Sviz'),
            html.H4('visualization for SNPeffect pipeline'),
        ], style={'textAlign': 'center', 'color': '#404040', 'backgroundColor': '#eeeefa', 'height': 120}),

        # MASS plots
        html.Div(html.H5(children="Individual MASS plots", style={'textAlign': 'center', 'marginTop': 10})),

        html.Div(children=L[0:L_half],  # left column of mass plots
                 style={'textAlign': 'center',
                        'width': '49%', 'display': 'inline-block', 'vertical-align': 'top', 'align': 'left', }),

        html.Div(children=L[L_half:len(L)],  # right column of mass plots
                 style={'textAlign': 'center', 'width': '49%', 'display': 'inline-block', 'align': 'right', }),

        # Dataset in detail
        html.Div(html.H5(children="Choose any dataset view", style={'textAlign': 'center', 'marginTop': 10})),

        html.Div(  # buttons for dataset view
            id='detail container',
            children=[
                html.Div(
                    children=[
                        dcc.Dropdown(
                            options=[{'label': name, 'value': name} for name in reports],
                            placeholder="Select dataset",
                            value=list(reports.keys())[0],
                            id="dropdownA",
                            style={'marginLeft': 10, 'marginRight': 10}
                        )], style={'width': '20%', 'display': 'inline-block', 'align': 'center', 'textAlign': 'center'}
                ),
                html.Div(
                    children=[
                        dcc.Dropdown(
                            options=[{'label': name, 'value': name} for name in list(reports.values())[0].columns
                                     if np.issubdtype(list(reports.values())[0][name].dtype, np.number)],
                            placeholder="Select variable for x",
                            value='abs_mutpos',
                            id="dropdownB",
                            style={'marginLeft': 10, 'marginRight': 10}
                        )], style={'width': '20%', 'display': 'inline-block', 'align': 'center', 'textAlign': 'center'}
                ),
                html.Div(
                    children=[
                        dcc.Dropdown(
                            options=[{'label': name, 'value': name} for name in list(reports.values())[0].columns
                                     if np.issubdtype(list(reports.values())[0][name].dtype, np.number)],
                            placeholder="Select variable for y",
                            value='position_mutation',
                            id="dropdownC",
                            style={'marginLeft': 10, 'marginRight': 10}
                        )], style={'width': '20%', 'display': 'inline-block', 'align': 'center', 'textAlign': 'center'}
                ),
                html.Div(
                    children=[
                        dcc.Dropdown(
                            options=[{'label': name, 'value': name} for name in list(reports.values())[0].columns
                                     if np.issubdtype(list(reports.values())[0][name].dtype, np.number)],
                            placeholder="Select var for color",
                            value='alignment_score',
                            id="dropdownD",
                            style={'marginLeft': 10, 'marginRight': 10}
                        )], style={'width': '20%', 'display': 'inline-block', 'align': 'center', 'textAlign': 'center'}
                ),
                html.Div(
                    children=[
                        dcc.Dropdown(
                            options=[{'label': 'linear', 'value': 'linear'}, {'label': 'log', 'value': 'log'}],
                            placeholder="Select axis type",
                            value='linear',
                            id="dropdownE",
                            style={'marginLeft': 10, 'marginRight': 10}
                        )], style={'width': '20%', 'display': 'inline-block', 'align': 'center', 'textAlign': 'center'}
                ),
            ], style={'marginLeft': 200, 'marginRight': 200, 'marginBottom': 10, 'backgroundColor': '#fafafe'}),

        # plot for dataset detail
        html.Div(
            dcc.Graph(id="det",
                      figure={'data': [go.Scatter(x=[0, 0], y=[0, 0])],
                              'layout': go.Layout(plot_bgcolor='#eeeefa', paper_bgcolor='rgba(0,0,0,0)',
                                                  margin=dict(t=20, r=0, b=20, l=20))}, style={'height': 600}),
            style={'marginLeft': 200, 'marginRight': 200, 'backgroundColor': '#fafafe'}),

        # Differential analysis block
        html.Div(html.H5(children="Choose any differential analysis", style={'textAlign': 'center', 'marginTop': 60})),

        html.Div(  # buttons for differential score
            id='differential container',
            children=[
                html.Div(
                    children=[
                        dcc.Dropdown(
                            options=[{'label': name, 'value': name} for name in reports],
                            placeholder="Select dataset",
                            value=list(reports.keys())[0],
                            id="dropdown1",
                            style={'marginLeft': 20, 'marginRight': 20}
                        )], style={'width': '33%', 'display': 'inline-block', 'align': 'center', 'textAlign': 'center'}
                ),
                html.Div(
                    children=[
                        dcc.Dropdown(
                            options=[{'label': name, 'value': name} for name in list(reports.values())[0].columns
                                     if np.issubdtype(list(reports.values())[0][name].dtype, np.number)],
                            placeholder="Select variable for x",
                            value='total_energy',
                            id="dropdown3a",
                            style={'marginLeft': 20, 'marginRight': 20}
                        )], style={'width': '33%', 'display': 'inline-block', 'align': 'center', 'textAlign': 'center'}
                ),
                html.Div(
                    children=[
                        dcc.Dropdown(
                            options=[{'label': name, 'value': name} for name in reports],
                            placeholder="Select dataset",
                            value=list(reports.keys())[1],
                            id="dropdown2",
                            style={'marginLeft': 20, 'marginRight': 20}
                        )], style={'width': '33%', 'display': 'inline-block', 'align': 'center', 'textAlign': 'center'}
                ),
            ], style={'marginLeft': 300, 'marginRight': 300, 'marginBottom': 10, 'backgroundColor': '#fafafe'}),

        # plot for differential score
        html.Div(
            dcc.Graph(id="difff",
                      figure={'data': [go.Scatter(x=[0, 0], y=[0, 0])],
                              'layout': go.Layout(plot_bgcolor='#eeeefa', paper_bgcolor='rgba(0,0,0,0)',
                                                  margin=dict(t=20, r=0, b=20, l=20))}, style={'height': 500}),
            style={'marginLeft': 300, 'marginRight': 300, 'backgroundColor': '#fafafe'}),

        # whole dataset table view
        html.H5(children="Full table",
                style={'textAlign': 'center', 'marginTop': 80, 'marginBottom': 15}),

        html.Div(
            children=[
                dcc.Dropdown(
                    options=[{'label': name, 'value': name} for name in reports],
                    placeholder="Select dataset",
                    id="dropdown4",
                    style={'align': 'center'}
                )], style={'width': '33%', 'display': 'inline-block', 'align': 'center', 'textAlign': 'center'}
        ),

        dash_table.DataTable(
            id='table',
            columns=[{"name": i, "id": i} for i in sorted(list(reports.values())[0].columns)],
            sort_action="custom",
            sort_mode="multi",
            sort_by=[],
            filter_action="native",
            style_table={'overflowX': 'scroll',
                         'maxHeight': '800px',
                         'overflowY': 'scroll',
                         'border': 'thin lightgrey solid',
                         },
            style_header={'backgroundColor': '#eeeefa'},
            style_cell={
                # all three widths are needed
                'minWidth': '60px', 'maxWidth': '300px',
                'overflow': 'hidden',
                'textOverflow': 'clip',
                'textAlign': 'left',
                'backgroundColor': '#fafafe',
            },
            style_as_list_view=False,
        ),
        html.Div(id="table-filter-container")
    ], style={'marginLeft': 150, 'marginRight': 150, 'marginTop': 10, 'marginBottom': 10,
              'backgroundColor': '#fafafe'})


@app.callback( #interactive segment for the detail plot
    Output("det", "figure"),
    [Input("dropdownA", "value"),
     Input("dropdownB", "value"),
     Input("dropdownC", "value"),
     Input("dropdownD", "value"),
     Input("dropdownE", "value")],
)
def update_graph(dropdownA, dropdownB, dropdownC, dropdownD, dropdownE):
    if (dropdownA is None) | (dropdownB is None) | (dropdownC is None) | (dropdownD is None):
        figure = {'data': [go.Scatter(x=[0, 0], y=[0, 0])],
                  'layout': go.Layout(plot_bgcolor='#eeeefa', paper_bgcolor='rgba(0,0,0,0)',
                                      margin=dict(t=20, r=0, b=40, l=40))}
        return figure

    repo = reports[dropdownA]
    figure = dict(
        data=[go.Scatter(
            x=repo[str(dropdownB)],
            y=repo[str(dropdownC)],
            text=[dropdownA + " variant: " + str(repo['variant'].tolist()[i])
                  + "<br>" + "Gene name: " + str(repo['gene_name'].tolist()[i])
                  + "<br>" + dropdownB + ": " + str(repo[dropdownB].tolist()[i])
                  + "<br>" + dropdownC + ": " + str(repo[dropdownC].tolist()[i])
                  + "<br>" + dropdownD + ": " + str(repo[dropdownD].tolist()[i])
                  for i in range(len(repo['gene_name'].tolist()))],
            mode='markers',
            opacity=0.9,
            marker={'size': 10,
                    'color': repo[str(dropdownD)],
                    'showscale': True,
                    'colorscale': [[0, 'rgba(140,40,220, 0.7)'],
                                   [0.04, 'rgba(30,144,255, 0.7)'],
                                   [0.1, 'rgba(50,205,50, 0.7)'],
                                   [0.3, 'rgba(250,220,1, 0.7)'],
                                   [0.7, 'rgba(255,140,0, 0.7)'],
                                   [1, 'rgba(178,34,34, 0.7)']],
                    'colorbar': dict(title=str(dropdownD), ticks='outside')})],
        layout=go.Layout(
            xaxis={'type': dropdownE, 'title': dropdownB},
            yaxis={'type': dropdownE, 'title': dropdownC},
            margin=dict(t=20, r=0, b=40, l=40),
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='#eeeefa',
            hovermode='closest',
        )
    )
    return figure


@app.callback(  # interactive segment for theh differential plot
    Output("difff", "figure"),
    [Input("dropdown1", "value"),
     Input("dropdown2", "value"),
     Input("dropdown3a", "value"), ],
)
def update_graph(dropdown1, dropdown2, dropdown3a):
    if (dropdown1 is None) | (dropdown2 is None) | (dropdown3a is None):
        figure = {'data': [go.Scatter(x=[0, 0], y=[0, 0])],
                  'layout': go.Layout(plot_bgcolor='#eeeefa', paper_bgcolor='rgba(0,0,0,0)',
                                      margin=dict(t=20, r=0, b=40, l=40))}
        return figure

    diff = differential(dropdown1, dropdown2, dropdown3a)
    figure = dict(
        data=[go.Scatter(
            x=diff[str(dropdown3a + '_x')],
            y=diff[str(dropdown3a + '_y')],
            text=["Gene name: " + str(diff.iloc[i, 0])
                  + "<br>------------------------------"
                  + "<br>" + dropdown1 + " variant: " + str(diff.iloc[i, 5])
                  + "<br> Mut " + str(diff.iloc[i, 1])
                  + "<br>" + str(dropdown3a) + ' : ' + str(diff.iloc[i, 3])
                  + "<br>------------------------------"
                  + "<br>" + dropdown2 + " variant: " + str(diff.iloc[i, 6])
                  + "<br> Mut " + str(diff.iloc[i, 2])
                  + "<br>" + str(dropdown3a) + ' : ' + str(diff.iloc[i, 4])
                  for i in range(len(diff['gene_name'].tolist()))],
            mode='markers',
            opacity=0.9,
            marker={'size': [9 if diff.iloc[i, 1] == diff.iloc[i, 2] else 12 for i in
                             range(len(diff['gene_name'].tolist()))],
                    'color': (diff[str(dropdown3a + '_x')] - diff[str(dropdown3a + '_y')]) /
                             (np.abs(diff[str(dropdown3a + '_x')]) + np.abs(diff[str(dropdown3a + '_y')])),
                    'showscale': True,
                    'symbol': [100 if diff.iloc[i, 1] == diff.iloc[i, 2] else 0 for i in
                               range(len(diff['gene_name'].tolist()))],
                    'colorscale': [[0, 'rgba(140,40,220, 0.7)'],
                                   [0.25, 'rgba(30,144,255, 0.7)'],
                                   [0.5, 'rgba(50,205,50, 0.7)'],
                                   [0.65, 'rgba(250,220,1, 0.7)'],
                                   [0.85, 'rgba(255,140,0, 0.7)'],
                                   [1, 'rgba(178,34,34, 0.7)']],
                    'colorbar': dict(title="ratio"
                                     # , tickvals=[0.1, 1, 10], ticks='outside'), 'cmin': 0.01, 'cmax': 10
                                     )})],
        layout=go.Layout(
            xaxis={'type': 'linear', 'title': dropdown1 + ' ' + dropdown3a},
            yaxis={'type': 'linear', 'title': dropdown2 + ' ' + dropdown3a},
            margin=dict(t=20, r=0, b=40, l=40),
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='#eeeefa',
            hovermode='closest',
        )
    )
    return figure


@app.callback( # interactive segment for the datatable
    [Output("table", "columns"), Output("table", "data")],
    [Input("dropdown4", "value"), Input("table", "sort_by")],
)
def update_table(dropdown4, sort_by):
    if dropdown4 is None:
        columns = [{"name": i, "id": i} for i in list(reports.values())[0].columns]
        data = list(reports.values())[0].to_dict("rows")
        return columns, [data[0]]
    else:
        columns = [{"name": i, "id": i} for i in reports[dropdown4].columns]
        data = reports[dropdown4].to_dict("rows")
        if sort_by is not None and len(sort_by):
            dfObj = pd.DataFrame(data)
            for col in sort_by:
                dfObj.sort_values(
                    by=col["column_id"],
                    ascending=(col["direction"] == "asc"),
                    inplace=True,
                )
            return columns, dfObj.to_dict("rows")
        return columns, data


@app.callback(
    Output('table-filter-container', "children"),
    [Input('table', "data")])
def update_graph(rows):
    if rows is None:
        dff = reports.values()
    else:
        dff = pd.DataFrame(rows)
    return html.Div()


if __name__ == '__main__':
    app.run_server(debug=False)
