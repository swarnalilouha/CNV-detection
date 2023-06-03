## Script for the detection of CNV's in the genome by simulation of read depthS across the genome using a gaussian distribution.
#!/usr/bin/env python

import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
%matplotlib inline 
import sys

## Reading in the file with read depth across many scaffolds
data = pd.read_table('_3_output_win__3_1_DS73500_TAAGGCGATATCCTCT_L001_R1_001_AHA4WPADXX.filt.fastq_intersect.txt', delim_whitespace=True, names=('chromosome', 'start', 'end', 'coverage'))
data.head()

# Calculates the window size from the start and end positions
data['window_size'] = data['end'] - data['start'].astype(int)
data_df = data[["chromosome", "start", "end", "window_size", "coverage"]].copy()

# Selecting the window size from the first row of the file and getting rid of rows with different window sizes at the end of each scaffold
is_window_size = data_df['window_size']==data_df.iloc[0,2]
data_df_window_size = data_df[is_window_size]

## Function for calculating statistics of distribution of read depths
def calculate_statistics(dataframe):
	stats = dataframe["coverage"].agg(['count', 'mean', 'min', 'max', np.std])
	stats = stats.round(2)
	stats_df = pd.DataFrame({'Mean': [stats["mean"]], 'Std dev':[stats["std"]], 'Count':[stats["count"]], 'Min':[stats["min"]], 'Max':[stats["max"]]}, index=[0])
	return stats_df

# Calculates the statistics of the raw distribution of read depths
stats1 = calculate_statistics(data_df_window_size)
stats1

# Function for plotting the distribution and its statistics
def make_plots(df, stats):
	fig = plt.figure(figsize=(10,10))
	fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
	fig.subplots_adjust(hspace=1)
	sns.distplot(df["coverage"], ax=ax1)

	ax2.axis('off')
	table = ax2.table(cellText=stats.values, colLabels=stats.columns, bbox=[0,0,1,1], loc='right')
	table.auto_set_font_size(False)
	table.set_fontsize(12)

	ax1.set_title("Distribution")
	ax2.set_title("Statistics")
	fig.savefig("Distribution_read_depth.pdf")
	plt.close(fig)

make_plots(data_df_window_size, stats1)

# Function for filtering out high read depths arising from high GC content/repeat regions
def remove_outliers(data, s):
	return data[abs( data["coverage"] - np.mean(data["coverage"]) < s * np.std(data["coverage"]))]

## Removes genomic regions which are 6 std dev away from the mean of the distribution
remove_repeats = remove_outliers(data_df_window_size, s=6)

## Calculates statistics of distribution after removal of genomic regions possibly containing repeats
stats2 = calculate_statistics(remove_repeats)
stats2

make_plots(remove_repeats, stats2)


##Code for CWL2
## Function for detecting amplifications 
def CNV(df1, m):
	return df1[abs(df1["coverage"] - np.mean(df1["coverage"]) > m * np.std(df1["coverage"]))]

# Detecting amplifications in genomic regions 2 std dev away from mean on the right of the distribution after removing repeats
copy_number_variations = CNV(remove_repeats,2)
copy_number_variations.head()


amplification = copy_number_variations[["chromosome", "start", "end", "coverage"]].copy()
#amplification['Color'] = 'green'
amplification.head(10)
#file1 = amplification.sort_values(by=['chromosome','start'])
#file1.to_csv("amplification.txt", sep="\t", index=False, header=False)
#amplification.index[1]

# Function for detecting deletions
def deletion(df2,n):
	return df2[abs(df2["coverage"] < np.mean(df2["coverage"]) - (n * np.std(df2["coverage"])))]

## Detecting amplifications in genomic regions 2 std dev away from mean on the left of the distribution after removing repeats
Deletions = deletion(remove_repeats,2)

drops = Deletions[["chromosome", "start", "end", "coverage"]].copy()
#drops['Color'] = 'red'
drops.head()

#def nothing(df3,m,n):
   # a =  df3[abs(df3["coverage"] > np.mean(df3["coverage"]) - m * np.std(df3["coverage"]))]
   # b =  df3[abs(df3["coverage"] > np.mean(df3["coverage"]) + m * np.std(df3["coverage"]))]
   # frame1 = [a,b]
   # frame1_concat = pd.concat(frame1)
   # frame1_remove_duplicates = frame1_concat.drop_duplicates(keep=False)
    #return frame1_remove_duplicates
    
   # c =  df3[abs(df3["coverage"] < np.mean(df3["coverage"]) + n * np.std(df3["coverage"]))]
   # d =  df3[abs(df3["coverage"] < np.mean(df3["coverage"]) - n * np.std(df3["coverage"]))]
   # frame2 = [c,d]
   # frame2_concat = pd.concat(frame2)
   # frame2_remove_duplicates = frame2_concat.drop_duplicates(keep=False)
   #return frame2_remove_duplicates
    
    #frame3 = [frame1_remove_duplicates,frame2_remove_duplicates]
    #result = pd.concat(frame3)
   # result = result.drop_duplicates(keep='first')
   # return result 

#Nothing = nothing(remove_repeats,2,2)
#non_significant = Nothing[["chromosome", "start", "end", "coverage"]].copy()
#non_significant['Color'] = 'blue'
#non_significant.head()
#non_significant.shape

# Final dataframe containing amplifications and deletions
all_frames = [amplification,drops]
plot_file = pd.concat(all_frames)
plot_file.shape

file1 = plot_file.sort_values(by=['chromosome','start'])
file1.head()

file1.to_csv("final.txt", sep="\t", index=False, header=False)
