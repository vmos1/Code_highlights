import dash
from dash import Dash, dcc, html, State, dash_table
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd 
from dash import dash_table 


from modules import *


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = html.Div([

    html.Div([
#        dcc.Markdown('Calculate '),
        dbc.Row(children=[
            dbc.Col(children=[
                html.Div([
                    html.H4("Simple Interest"),
                    html.Label("Principal Amount:"),
                    dcc.Input(id="principal", type="number", value=1000),
                    html.Br(), html.Br(),
                    html.Label("Rate of Interest (% per term):"),
                    dcc.Input(id="rate", type="number", value=5),
                    
                    html.Br(), html.Br(),
                    html.Label("Time Period (terms):"),
                    # dcc.Slider(id='time',min=1,max=10,step=1,value=3),
                    dcc.Dropdown(id='time',options = np.arange(100),value=3),

                    html.Br(), html.Br(),
                    html.H6("Final amount"),
                    html.Div(id="output-interest"),
                    html.Br(), html.Br(),
                ],),

                html.Details([
                    html.Summary("Show table "),
                    html.Div([
                        dash_table.DataTable(
                            id='table1',
                            data=[],  # Empty data at first, will be filled by callback

                            style_table={'overflowX': 'auto'},  # Make table scrollable horizontally if needed
                            style_cell={'textAlign': 'left'},  # Align text to the left
                        )
                    ]),
                ],open=False),
            ], style={"padding":"20px"}, md = 6),
            dbc.Col(children=[
               html.Div([
                    html.H4("Compound Interest"),
                    html.Label("Principal Amount:"),
                    dcc.Input(id="principal2", type="number", value=1000),
                    html.Br(), html.Br(),
                    html.Label("Rate of Interest (% per term):"),
                    dcc.Input(id="rate2", type="number", value=5),
                    
                    html.Br(), html.Br(),
                    html.Label("Time Period (terms):"),
                    dcc.Dropdown(id='time2',options = np.arange(100),value=3),
                    html.Br(), html.Br(),

                    html.H6("Final amount"),
                    html.Div(id="output-interest2"),
                    html.Br(), html.Br(),
                ],),
                html.Details([
                    html.Summary("Show table "),
                    html.Div([
                        dash_table.DataTable(
                            id='table2',
                            data=[],  # Empty data at first, will be filled by callback
                            style_table={'overflowX': 'auto'},  # Make table scrollable horizontally if needed
                            style_cell={'textAlign': 'left'},  # Align text to the left
                        )
                    ]),
                ]),
            ], style={"padding":"20px"}, md=6),
        ]),
        html.Hr(),
        dbc.Row(children=[
            dbc.Col(children=[
                html.Div([                    
                    html.Br(), html.Br(),
                    
                    html.H2("SIP amount"),

                    html.Label("Periodic contribution"),
                    dcc.Input(id='contrib',type='number',value=1000),

                    html.Br(), html.Br(),
                    
                    html.Label("Rate of Interest (% per term):"),
                    dcc.Input(id="rate3", type="number", value=5),
                    
                    html.Br(), html.Br(),
                    
                    html.Label("Time Period (terms):"),
                    dcc.Dropdown(id='n3',options = np.arange(100),value=10),

                    html.H6("Final amount"),
                    html.Div(id="sip_output"),
                ]),
                html.Details([
                    html.Summary("Show table "),
                    html.Div([
                        dash_table.DataTable(
                            id='table3',
                            data=[],  # Empty data at first, will be filled by callback
                            style_table={'overflowX': 'auto'},  # Make table scrollable horizontally if needed
                            style_cell={'textAlign': 'left'},  # Align text to the left
                        )
                    ]),
                ],open=False),
            ], style={"padding":"20px"}, md = 4),
            dbc.Col(children=[
                html.Div([                    
                    html.Br(), html.Br(),
                    html.H3("EMI payment amount"),

                    html.Label("Total loan"),
                    dcc.Input(id='loan',type='number',value=1000),

                    html.Br(), html.Br(),
                    
                    html.Label("Rate of Interest (% per term):"),
                    dcc.Input(id="rate4", type="number", value=5),
                    
                    html.Br(), html.Br(),
                    
                    html.Label("Time Period (terms):"),
                    dcc.Dropdown(id='n4',options = np.arange(1,100),value=10),

                    html.H6("Term (monthly) payment"),
                    html.Div(id="emi_term"),
                ]),
                html.Details([
                    html.Summary("Show table "),
                    html.Div([
                        dash_table.DataTable(
                            id='table4',
                            data=[],  # Empty data at first, will be filled by callback
                            style_table={'overflowX': 'auto'},  # Make table scrollable horizontally if needed
                            style_cell={'textAlign': 'left'},  # Align text to the left
                        )
                    ]),
                ],open=False),
            ], style={"padding":"20px"}, md=4),
            dbc.Col(children=[
                html.Div([
                ]),
                html.Div([                    
                    # html.Br(), html.Br(),
                    html.H3("EMI period"),

                    html.Label("Total loan"),
                    dcc.Input(id='loan2',type='number',value=1000),

                    html.Br(), html.Br(),
                    
                    html.Label("Rate of Interest (% per term):"),
                    dcc.Input(id="rate5", type="number", value=5),
                    
                    html.Br(), html.Br(),
                    
                    html.Label("EMI:"),
                    dcc.Input(id='emi2', type="number",value=500),

                    html.H6("Number of Terms "),
                    html.Div(id="terms"),
                ]),
                html.Details([
                    html.Summary("Show table "),
                    html.Div([
                        dash_table.DataTable(
                            id='table5',
                            data=[],  # Empty data at first, will be filled by callback
                            style_table={'overflowX': 'auto'},  # Make table scrollable horizontally if needed
                            style_cell={'textAlign': 'left'},  # Align text to the left
                        )
                    ]),
                ],open=False),
            ], style={"padding":"20px"}, md=4),
        ])

    ]),
])

@app.callback(
    Output("output-interest", "children"),
    Output("table1", "data"),
    Input("principal", "value"),
    Input("rate", "value"),
    Input("time", "value")
)
def calculate_interest(principal, rate, time):
    amount,df = f_interest(principal, rate/100, time, 'simple')
    return amount,df.to_dict('records')

@app.callback(
    Output("output-interest2", "children"),
    Output("table2", "data"),
    Input("principal2", "value"),
    Input("rate2", "value"),
    Input("time2", "value")
)
def calculate_interest(principal, rate, time):
    amount,df = f_interest(principal, rate/100, time, 'compound')
    return amount ,df.to_dict('records')

@app.callback(
    Output('sip_output',"children"),
    Output("table3", "data"),
    Input("contrib","value"),
    Input("rate3","value"),
    Input("n3","value"),
)
def calc_sip(contrib,rate3,n2):
    amount,df = f_SIP(contrib,rate3/100,n2)
    return amount,df.to_dict('records')


@app.callback(
    Output('emi_term','children'),
    Output("table4", "data"),
    Input("loan","value"),
    Input("rate4","value"),
    Input("n4","value"),
)
def calc_emi_amount(loan,rate,terms):
    amount,df= f_EMI(loan,rate/100,terms,mode='amount')
    return amount,df.to_dict('records')

@app.callback(
    Output('terms','children'),
    Output("table5", "data"),
    Input("loan2","value"),
    Input("rate5","value"),
    Input("emi2","value"),
)
def calc_emi_terms(loan, rate, emi):

    terms,df = f_EMI(loan, rate/100, emi, mode='period')
    return terms, df.to_dict('records')

if __name__ == '__main__':
    app.run_server(debug=True)
