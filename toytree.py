#!/usr/bin/python

import numpy as np
import toytree
import toyplot
import json
import toyplot.pdf
import sys

sars = sys.argv[1]
jso = sys.argv[2]

def get_clade_dict(json_path, tree):

  def clade_getter(json_path):
    file = open(json_path)
    data = json.load(file)
    file.close()
    return data

  clade_dict = clade_getter(json_path)
  nodes_dict = clade_dict['nodes']

  tip_keys = tree.get_tip_labels()
  
  clade_dict = {k: nodes_dict[k] for k in tip_keys}

  clade_membs = list(clade_dict.values())
  clades = [i['clade_membership'] for i in clade_membs]
  nodes = list(clade_dict.keys())

  clade_dict = dict(zip(nodes, clades))
  
  return clade_dict

def get_emap(dict, clades, color_list=['#b4d4f2', '#008080', '#0ff1ce', '#ff99ff', '#f0b456', '#5956f0', '#ff0019']):
  
  def get_keys_from_values(d, val):
    return [k for k, v in d.items() if v == val]

  keylist = [get_keys_from_value(dict, i) for i in clades]
  
  color_list = color_list[:len(keylist)]

  emap = {}
  for i in range(len(keylist)):
    j = tuple(keylist[i])
    v = color_list[i]
    emap[j] = v
    
  return emap

def get_figure(tree, emap, pth = "./", name = "tree-plot", type = "pdf", height = 1400, width = 800, tip_labels = False, tip_labels_align = False, ew = 5):
  
  ecolors = tree.get_edge_values_mapped(emap)
  elabels =  tree.get_edge_values('idx')
  ewidths = [ew for i in elabels]
  
  canvas, _, _ = tree.draw(height = height, width = width, tip_labels = tip_labels, tip_labels_align = tip_labels_align, edge_colors = ecolors, edge_widths = ewidths)
  
  fname = pth + name + "." + type
  if type == "pdf":
    toyplot.pdf.render(canvas, fname)
  elif type == "svg":
    toyplot.svg.render(canvas, fname)
  
  return tree.draw(height = height, width = width, tip_labels = tip_labels, tip_labels_align = tip_labels_align, edge_colors = ecolors, edge_widths = ewidths)

tre1 = toytree.tree(sars, tree_format=1)

clade_dict = get_clade_dict(jso, tre1)

c_uni = list(np.unique(np.array(list(clade_dict.values()))))

emap = get_emap(clade_dict, c_uni)

j = get_figure(tre1, emap)
