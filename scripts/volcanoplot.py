#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import sys
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()

#if len(sys.argv[1:]) < 4:
#	sys.stderr.write('Error: missing argument\nUsage:\nvolcanoplot.py file.tsv log color_key plot.png\n\nOptional: labels_file with coordinates and text that will be overlayed to the scatterplot\nlog [t | f | negative; default is f]\n')
	#sys.exit(1)
if len(sys.argv[1:]) == 4:
	file1,log_scale,color_key,outfile = sys.argv[1:]
elif len(sys.argv[1:]) == 5:
	file1,log_scale,color_key,outfile,labels_file =  sys.argv[1:]
else:
	sys.stderr.write('Error: missing argument\nUsage:\nvolcanoplot.py <file.tsv> <log> <color_key> <plot.png> <labels_file.tsv>\n\nlog:  [t | f | negative; default is f]\ncolor_key: column header to select colors for the plot\nOptional: labels_file: file with coordinates and text that will be overlayed to the scatterplot\n')
	sys.exit(1)
file1 = pd.read_csv(file1,sep='\t',index_col=[0], na_filter=True, na_values=0)
file1 = file1.fillna(0)
log_scale = str(log_scale)
x=file1[file1.columns[0]].astype(float).values
y=file1[file1.columns[1]].astype(float).values
color = file1[color_key].values.tolist()
if log_scale == 'True' or log_scale == 'T' or log_scale == 't' or log_scale == 'TRUE':
	x = (np.log10(x))
	y = (np.log10(y))
if log_scale == 'negative':
	y = (-np.log10(y))

ax1.scatter(x,y,color=color,s=80, edgecolor='none',alpha=0.7)

ax1.set_xlabel(file1.columns[0],fontsize=12)
ax1.set_ylabel(file1.columns[1],fontsize=12)
xmax = float(np.abs(x).max())
xmargin = float(xmax)/50
xlim = xmax+xmargin
ax1.set_xlim(-xlim,xlim)

ax2.tick_params(right='off', labelright='off')
ax1.tick_params(direction='out')
ax1.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_position('zero')
ax2.spines['left'].set_linestyle('--')
ax2.spines['left'].set_color('gray')

ax1.spines['top'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')

#create legend
set_labels = set(file1['histotype'].values.tolist())
handles = []
for l in set_labels:
	color_x = pd.unique(file1[file1['histotype']==l][color_key].values)[0]
	patch = mpatches.Patch(color=color_x, label=l)
	handles.append(patch)

plt.legend(handles=handles,markerscale=3,borderpad=0.5,labelspacing=0.5,loc='upper center',fontsize=10,ncol=3,bbox_to_anchor=(0.5,1.1),fancybox=True,shadow=False,frameon=False)


axes = plt.gca()
if len(sys.argv[1:]) == 5:
	print 'adding text...'
	labels_file = pd.read_csv(labels_file,sep='\t')
	if log_scale=='negative':
		labels_file[labels_file.columns[2]]=-np.log10(labels_file[labels_file.columns[2]].astype(float).values)
		labels_file[labels_file.columns[4]]=-np.log10(labels_file[labels_file.columns[4]].astype(float).values)
	ax1.scatter(labels_file[labels_file.columns[1]],labels_file[labels_file.columns[2]],color='orange',s=80, edgecolor='none',alpha=0.7)
	for i in xrange(labels_file.shape[0]):
		ax1.annotate(labels_file[labels_file.columns[0]][i], xy=(labels_file[labels_file.columns[1]][i],labels_file[labels_file.columns[2]][i]),xytext=(labels_file[labels_file.columns[3]][i],labels_file[labels_file.columns[4]][i]), arrowprops=dict(facecolor='black',arrowstyle='->'), fontsize=12)
	print 'done.'
#plt.plot([llim,ulim],[llim,ulim],'--',color='black')
plt.savefig(outfile,format='png', dpi=300)
plt.close()
