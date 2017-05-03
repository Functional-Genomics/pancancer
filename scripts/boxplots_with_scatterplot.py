#!/usr/bin/env python
import scipy as sp
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import sys
import seaborn as sns
sns.set_style("whitegrid")
#plt.style.use('ggplot')


fig = plt.figure(figsize=(15,15))

font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 30,
        }

if len(sys.argv[1:])<8:
	sys.stderr.write('Error: missing argument\nUsage:\nscatterplot.py <file.tsv> <log> <x> <y> <hue> <file_colors> <title> <plot.png> [label_order_file.txt]\nlog = [ t | f ]\nlabel_order_file is optional\n')
	sys.exit(1)

if len(sys.argv[1:])==8:
	table,log,xval,yval,hueval,file_colors,title,outfile=sys.argv[1:]
else:
	table,log,xval,yval,hueval,file_colors,title,outfile,gene_order=sys.argv[1:]
	gene_order = pd.read_csv(gene_order,sep='\t')
	order = gene_order[gene_order.columns[0]].values.tolist()
	
table = pd.read_csv(table,sep='\t',index_col=[0], na_filter=True, na_values=0)
table = table.fillna(0)
log = str(log)
xval = str(xval)
yval = str(yval)
hueval = str(hueval)
file_colors = pd.read_csv(file_colors,sep='\t',header=None)
title = str(title)
#gene_order = pd.read_csv(sys.argv[8],sep='\t')
#order = gene_order[gene_order.columns[0]].values.tolist()
colors_dict = file_colors.set_index(0)[1].to_dict()

if log == 't':
	table[table.columns[1]]=sp.log10(table[table.columns[1]])
#ax = sns.swarmplot(x=table.columns[0], y=table.columns[1], hue=table.columns[2],data=table,palette=colors_dict,size=10)
#bp=sns.boxplot(x=xval, y=yval,data=table,fliersize=0)
if len(sys.argv[1:])==8:
	ax=sns.stripplot(x=xval, y=yval, hue=hueval,data=table,palette=colors_dict,size=10,jitter=0.1,edgecolor="black",alpha=0.4)
	bp=sns.boxplot(x=xval, y=yval,data=table,fliersize=0,palette=colors_dict)
else:
	ax=sns.stripplot(x=xval, y=yval, hue=hueval,data=table,palette=colors_dict,size=10,jitter=0.1,edgecolor="black",alpha=0.4,order=order)
	bp=sns.boxplot(x=xval, y=yval,data=table,order=order,fliersize=0,palette=colors_dict)
plt.setp(bp.artists, alpha=0.5)
bp.artists[0].set_edgecolor('black')
bp.artists[0].set_linewidth(2)
ax.set_ylabel(yval,fontdict=font)
ax.set_xlabel('',fontdict=font)
bp.set_title(title,fontdict=font)
sns.set_style("whitegrid")
plt.xticks(rotation=90,fontsize=14)
plt.yticks(rotation=0,fontsize=24)
plt.legend(bbox_to_anchor=(1.07,-0.2), loc=1,ncol=5, borderaxespad=0.,fontsize=21)
plt.subplots_adjust(bottom=0.3)
plt.savefig(outfile,format='pdf', dpi=400)

