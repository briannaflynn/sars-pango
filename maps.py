#!/usr/bin/python

import plotly.graph_objects as go
# version 5.1.0
import pandas as pd
# version 1.1.5
from PyPDF2 import PdfFileWriter, PdfFileReader
# version version 1.26.0
import sys
import toytree
import json

id = sys.argv[1]

def get_clade_dict(json_path, tree):

  def clade_getter(json_path):
    file = open(json_path)
    data = json.load(file)
    file.close()
    return data

  clade_dict = clade_getter(json_path)
  nodes_dict = clade_dict['nodes']

  tip_keys = tree.get_tip_labels()
  
  clade_dict = {k: nodes_dict[k] for k in tip_keys if k != "SAMPLE"}

  clade_membs = list(clade_dict.values())
  clades = [i['clade_membership'] for i in clade_membs]
  nodes = list(clade_dict.keys())

  clade_dict = dict(zip(nodes, clades))
  clade_dict['SAMPLE'] = 'Sample'
  
  return clade_dict
  
tree = toytree.tree(f"{id}-results/tree/tree_raw.nwk", tree_format=1)
clade_dict = get_clade_dict(f"{id}-results/tree/clades.json", tree)

names = list(clade_dict.keys())

def get_countries(names:list):
  countries = []
  for i in range(len(names)):
    n = names[i]
    n = n.split('/')
    n = n[0]
    countries.append(n)

  return countries

countries = get_countries(names)

df = pd.DataFrame()
df['Virus name'] = names
df['clade'] = list(clade_dict.values())
df['countries'] = countries

clade_df = df

clade_df.replace({'countries':'USA'}, 'United States', inplace = True)
clade_df.replace({'countries':'Wuhan'}, 'China', inplace = True)
clade_df.replace({'countries':'mink'}, 'Netherlands', inplace = True)
clade_df.replace({'countries':'England'}, 'United Kingdom', inplace = True)
clade_df.replace({'countries':'Scotland'}, 'United Kingdom', inplace = True)

sampleClade = pd.read_csv(f'{id}-results/clade_assignment.tsv', delimiter = "\t")
clade = dict(sampleClade['clade'])
clade = clade[0]

df = clade_df
df = df.loc[df['clade'] == clade]
df = pd.value_counts(df.countries).to_frame().reset_index()
df.columns = ['COUNTRY','count']

wrld_df = pd.read_csv('world_codes.csv')

counts = wrld_df.merge(df, on = "COUNTRY")

def perc(dataframe):
  counts = list(dataframe['count'])
  s = sum(counts)
  percs = [round(c/s * 100) for c in counts]
  percs_df = pd.DataFrame({'count':percs})
  dataframe.update(percs_df)
  return dataframe

df = perc(counts)

title = 'Global Prevalence of SARS-CoV-2 Strain ' + clade
fig = go.Figure(data=go.Choropleth(
    
    locations = df['CODE'],
    z = df['count'],
    text = df['COUNTRY'],
    colorscale = 'YlOrRd',
    autocolorscale=False,
    reversescale=False,
    marker_line_color='darkgray',
    marker_line_width=0.5,
    colorbar_ticksuffix = '%',
    colorbar_title = 'Percentage of total<br>' + clade + ' samples',
))

fig.update_layout(
    title_text=title,
    geo=dict(
        showframe=False,
        showcoastlines=False,
        projection_type='equirectangular'
    ),
    annotations = [dict(
        x=0.55,
        y=0.1,
        xref='paper',
        yref='paper',
        text="",
        showarrow = False
    )]
)
fig.update_geos(fitbounds="locations")

#uncomment to populate an interactive figure in web browser
#fig.show()

fig.write_image(f'{id}_map.pdf')

output = PdfFileWriter() 
input = PdfFileReader(open(f'{id}_map.pdf', 'rb')) 

n = input.getNumPages()

for i in range(n):
  page = input.getPage(i)
  page.cropBox.upperLeft = (0,60)
  page.cropBox.lowerRight = (800,1200)
  output.addPage(page) 
  
outputStream = open(f'{id}_result.pdf','wb') 
output.write(outputStream) 
outputStream.close()
