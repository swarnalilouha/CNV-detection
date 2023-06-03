##Scrip for plotting read depths across specific genomic regions
#!/usr/bin/env python

import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
%matplotlib inline 
import glob
import itertools
from matplotlib.backends.backend_pdf import PdfPages
from functools import reduce

##Takes in folder and filenames ending with ".txt"
filenames= glob.glob("shasta_1800bp/*.txt")
df = []
for name in filenames:
	if name.endswith(".txt"):
		file_name=name.partition(".")
		df.append(pd.read_table(name, delim_whitespace=True, names=('chromosome', 'start', 'end', (file_name[0]))))
		file_name = ()

df_merge = reduce(lambda x, y: pd.merge(x, y, on = ['chromosome', 'start', 'end']), df)
df_merge.shape

def select_scaffold(scaffold, start, end):
	is_chromosome = df_merge['chromosome'].str.contains(scaffold)
	df_chromosome = df_merge[is_chromosome]
	is_range = df_chromosome['start'].between(start, end)
	df_range = df_chromosome[is_range]
	return df_range

## Enter scaffold name, start pos, and end pos of genomic region to be plotted
list = select_scaffold("NC_036633",0,692922)
list.head()

## Plots genomic regions with matplotlib
fig = plt.figure(figsize=(20,20))
fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True) ##Depending on the number of sample, make changes in nrows and ncols
fig.subplots_adjust(hspace=1,wspace=0.4)

fig.text(0.5, 0.001, 'window_size (bp)', ha='center')
fig.text(0.001, 0.5, 'Read depth', va='center', rotation='vertical')
pdf_filename = 'Scaffold_22_amp.pdf'
pdf = PdfPages(pdf_filename)

## The user needs to input the axis for every sample plotted. The user also needs to specify the names of the files, and set the title for every sample on the plot
axs[0,0].plot(list['start'], list['shasta_1800bp/_1_output_win_7000358660_S4_R2_001'])
axs[0,0].set_title('7000358660')

axs[0,1].plot(list['start'], list['shasta_1800bp/_1_output_win_7000358662_S5_R2_001'])
axs[0,1].set_title('7000358662')

axs[1,0].plot(list['start'], list['shasta_1800bp/_2_output_win__1_9008537374_S2_R1_001'])
axs[1,0].set_title('9008537374')

## Code used to export the plot out in a pdf file
pdf.savefig(fig)
plt.close(fig)
pdf.close()
           
