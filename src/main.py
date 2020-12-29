"""
    Genome Analysis Project
"""
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd
import dash_bootstrap_components as dbc
import re

from dash.dependencies import Input, Output
from string import digits
from Bio.Seq import Seq


""" 
    During transcription, the RNA polymerase read the template DNA strand in the 3′→5′ direction, but the mRNA is formed in the 5′ to 3′ direction.
    The codons of the mRNA reading frame are translated in the 5′→3′ direction into amino acids by a ribosome to produce a polypeptide chain.  

    for graphs use plotly or seaborn

    Show content of each nucleotide for coding, template and mRNA strand, GC content
    Percent expected mutations
    Expected amino acide chain -> create 3d structure?

    ~10−3 per generation for satellite DNA expansion/contraction
    human genome mutation rate is similarly estimated to be ~1.1×10−8 per site per generation

    http://prody.csb.pitt.edu/
"""

external_stylesheets = [dbc.themes.CYBORG]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

colors = {
    'background': '#040404',
    'text': '#7FDBFF'
}

app.layout = html.Div(children=[
    
    dbc.Row(dbc.Col(dcc.Input(id="text", type="text", placeholder="enter coding strand", value=''), width={'size':6, "offset":0, 'order':1}),
    ),


    dbc.Row([dbc.Col(dcc.Graph(id="bar-chart-coding"), width={'size':3, "offset":0, 'order':0}),
    dbc.Col(dcc.Graph(id="bar-chart-template"), width={'size':3, "offset":0, 'order':1}),
    dbc.Col(dcc.Graph(id="bar-chart-mRNA"), width={'size':3, "offset":0, 'order':2}),]
    ),
    
    html.H4("Amino Acid Sequence: "),
    html.Div(id='output'),
])

@app.callback(
    Output('bar-chart-coding', 'figure'),
    Output('bar-chart-template', 'figure'),
    Output('bar-chart-mRNA', 'figure'),
    Output('output', 'children'),
    Input('text', 'value')
)
def update_graph (text_input):
    text_input = text_input.upper()
    if len(re.findall("[^ATCG]+", text_input)) > 0:
        text_input = ""

    coding_strand = Seq(text_input)
    template_strand = coding_strand.complement()
    mRNA = coding_strand.transcribe()
    amino_acid_sequence = mRNA.translate()

    coding_strand_counts = [coding_strand.count('A'), coding_strand.count('T'), coding_strand.count('C'), coding_strand.count('G')]
    template_strand_counts = [template_strand.count('A'), template_strand.count('T'), template_strand.count('C'), template_strand.count('G')]
    mRNA_strand_counts = [mRNA.count('A'), mRNA.count('U'), mRNA.count('C'), mRNA.count('G')]

    coding_strand_data = {'Base': ['A','T','C','G'],
                        'Count': coding_strand_counts
                        }
    template_strand_data = {'Base': ['A','T','C','G'],
                        'Count': template_strand_counts
                        }
    mRNA_strand_data = {'Base': ['A','U','C','G'],
                        'Count': mRNA_strand_counts
                        }

    coding_strand_df = pd.DataFrame(coding_strand_data, columns = ['Base','Count'])
    template_strand_df = pd.DataFrame(template_strand_data, columns = ['Base','Count'])
    mRNA_strand_df = pd.DataFrame(mRNA_strand_data, columns = ['Base','Count'])

    coding_strand_fig = px.bar(coding_strand_df, x='Base', y='Count', color='Base', color_discrete_map={'A':'lightcyan',
                'T':'cyan',
                'C':'royalblue',
                'G':'darkblue'}, title='Coding Strand')
    template_strand_fig = px.bar(template_strand_df, x='Base', y='Count', color='Base', color_discrete_map={'A':'lightcyan',
                'T':'cyan',
                'C':'royalblue',
                'G':'darkblue'}, title='Template Strand')
    mRNA_strand_fig = px.bar(mRNA_strand_df, x='Base', y='Count', color='Base', color_discrete_map={'A':'lightcyan',
                'U':'cyan',
                'C':'royalblue',
                'G':'darkblue'}, title='mRNA Strand')

    coding_strand_fig.update_layout(
    plot_bgcolor=colors['background'],
    paper_bgcolor=colors['background'],
    font_color=colors['text']
        )
    template_strand_fig.update_layout(
        plot_bgcolor=colors['background'],
        paper_bgcolor=colors['background'],
        font_color=colors['text']
        )
    mRNA_strand_fig.update_layout(
        plot_bgcolor=colors['background'],
        paper_bgcolor=colors['background'],
        font_color=colors['text']
        )

    return coding_strand_fig, template_strand_fig, mRNA_strand_fig, f"{amino_acid_sequence}"

if __name__ == '__main__':
    app.run_server(debug=True)