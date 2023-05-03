# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-04-27 23:12:24
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-04-27 23:40:56

import matplotlib.pyplot as plt
import statistics
import numpy as np
import re
import scipy
import sys
import pandas as pd
import seaborn as sns
import matplotlib as mpl

def create_figure_df(coverage_file):
    coverages = {}
    xaxis = []
    yaxis = []
    
    with open(coverage_file, "r") as cf:
        for line in cf:
            cov_list = line.split(", ")
            for cov in cov_list:
                nums = re.search('([0-9]+)\: ([0-9]+)', cov)
                if nums is None:
                    break
                xaxis.append(int(nums.group(1)) + 1)
                yaxis.append(int(nums.group(2)))
                coverages.update({nums.group(1): nums.group(2)})
        avg = round(statistics.mean(yaxis), 2)
    x_y_spline = scipy.interpolate.make_interp_spline(xaxis, yaxis)
    x_ = np.linspace(min(xaxis), max(xaxis), 500)
    y_ = x_y_spline(x_)
    df_fig = pd.DataFrame()
    df_fig['nt'] = xaxis
    df_fig['count'] = yaxis
    return df_fig

### code_Coverage Histograms ###

def list_element_count(LIST_IN=[]):
    '''
    ### Count the number of every element in a list
    '''
    OUTPUT = {}
    for element in LIST_IN:
        if element in OUTPUT:
            OUTPUT[element] += 1
        else:
            OUTPUT[element]  = 1
    return OUTPUT

def draw(virus_name, raw_coverage_data, correct_coverage_data):
    df_fig_cor= create_figure_df(correct_coverage_data)
    df_fig_raw= create_figure_df(raw_coverage_data)

    xaxis1 = df_fig_cor['nt'].tolist()
    yaxis1 = df_fig_cor['count'].tolist()
    xaxis2 = df_fig_raw['nt'].tolist()
    yaxis2 = df_fig_raw['count'].tolist()

    print (f"Raw_mean: {np.mean(yaxis2)}, Cor_mean: {np.mean(yaxis1)}")

    # set default line width and font size
    mpl.rcParams['axes.linewidth'] = 2  # set the default line width of the axis
    mpl.rcParams['axes.labelsize'] = 24  # set the default font size of the axis labels

    fig, ax = plt.subplots(figsize=(24, 10))
    # customize the plot margins
    plt.subplots_adjust(top=0.95, left=0.1, right=0.95)

    fig = sns.lineplot(data=df_fig_cor, x='nt', y='count', color='#e41a1c', lw=3, label="Corrected")
    fig = sns.lineplot(data=df_fig_raw, x='nt', y='count', color='#000000', lw=3, label="Original") #style=True, dashes=[(2,2)], 
    # axes = plt.gca()
    fig.legend(fontsize=24, frameon=False) #labels=['Original', 'Corrected']

    # customize the axis lines
    # axes = plt.gca()  # get the current axes
    ax.spines['left'].set_position(('outward', 8))  # set the position of the left axis
    ax.spines['bottom'].set_position(('outward', 8))  # set the position of the bottom axis
    ax.spines['right'].set_visible(False)  # hide the right axis
    ax.spines['top'].set_visible(False)  # hide the top axis
    plt.xlim(left=0)  # set the lower limit of the x axis to 0
    plt.ylim(bottom=0)  # set the lower limit of the y axis to 0
    # set the x-axis range to [0, 200000]
    plt.xlim([0, 200000])
    plt.ylim([0, 6500])
    plt.yticks([0, 1000, 2000, 3000, 4000, 5000, 6000, 6500]) 
    # ax.spines['left'].set_linewidth(10)

    # customize the tick label font size
    plt.xticks(fontsize=24)  # set the font size of the x axis tick labels to 14
    plt.yticks(fontsize=24)  # set the font size of the y axis tick labels to 14
    plt.xlabel('Location (nt)', fontsize=32, fontname="serif")  # set the x axis label to "Time (s)" with font size 16
    plt.ylabel('Coverage', fontsize=32, fontname="serif")  # set the y axis label to "Voltage (V)" with font size 16

    plt.savefig(virus_name + "_coverage.png")
    plt.close()

    ##########################################################################################################################
    ll = []
    for i in range(len(yaxis1)):
        difference = yaxis1[i] - yaxis2[i]
        ll.append(difference)
        if difference < 0:
            print(i)
            
    df_fig = pd.DataFrame()
    df_fig['nt'] = xaxis1
    df_fig['count'] = ll

    fig, ax = plt.subplots(figsize=(24, 12))
    # customize the plot margins
    plt.subplots_adjust(bottom=0.1, top=0.3)

    fig = sns.lineplot(data=df_fig, x='nt', y='count', color='#e41a1c', lw=2, zorder=1)
    plt.fill_between(df_fig['nt'], df_fig['count'], color='#e41a1c')
    # fig.legend(fontsize=14, frameon=False) #labels=['Original', 'Corrected']

    ax.spines['left'].set_position(('outward', 8))  # set the position of the left axis
    ax.spines['bottom'].set_position(('outward', 8))  # set the position of the bottom axis
    ax.spines['right'].set_visible(False)  # hide the right axis
    ax.spines['top'].set_visible(False)  # hide the top axis
    plt.xlim(left=0)  # set the lower limit of the x axis to 0
    plt.ylim(bottom=0)  # set the lower limit of the y axis to 0
    # set the x-axis range to [0, 200000]
    plt.xlim([0, 200000])

    plt.ylim([0, 2500])
    plt.yticks([0, 1000, 2000, 2500])  # set the y-axis tick values to 0, 4, and 8

    # customize the tick label font size
    plt.xticks(fontsize=24)  # set the font size of the x axis tick labels to 14
    plt.yticks(fontsize=24)  # set the font size of the y axis tick labels to 14
    plt.xlabel('Location (nt)', fontsize=32, fontname="serif")  # set the x axis label to "Time (s)" with font size 16
    plt.ylabel('Coverage\nDifference', fontsize=32, fontname="serif")  # set the y axis label to "Voltage (V)" with font size 16

    plt.savefig(virus_name + ".difference.png")
    plt.close()

    dd = list_element_count(df_fig_cor['count'].tolist())

    dfdf = pd.DataFrame(dd.items(), columns=['keys', 'values'])
    dfdf = dfdf.sort_values('keys', ascending=False, ignore_index=True)

    # fig, ax = ss.new_fig(10, 10)
    fig, ax = plt.subplots(figsize=(14, 10))

    fig = sns.histplot(data=df_fig_raw, x='count', color='#000000', alpha=0.3, 
                    kde=True, line_kws={'linewidth':3}, label="Original")
    fig = sns.histplot(data=df_fig_cor, x='count', color='#ff0000', alpha=0.3, 
                    kde=True, line_kws={'linewidth':3}, label="Corrected")

    fig.legend(fontsize=24, frameon=False) #labels=['Original', 'Corrected']

    ax.spines['left'].set_position(('outward', 8))  # set the position of the left axis
    ax.spines['bottom'].set_position(('outward', 8))  # set the position of the bottom axis
    ax.spines['right'].set_visible(False)  # hide the right axis
    ax.spines['top'].set_visible(False)  # hide the top axis
    plt.xlim(left=0)  # set the lower limit of the x axis to 0
    plt.ylim(bottom=0)  # set the lower limit of the y axis to 0
    # set the x-axis range to [0, 200000]
    plt.xlim([0, 6500])
    plt.ylim([0, 4500])
    plt.yticks([0, 1000, 2000, 3000, 4000, 4500]) 
    # ax.spines['left'].set_linewidth(10)

    # customize the tick label font size
    plt.xticks(fontsize=24)  # set the font size of the x axis tick labels to 14
    plt.yticks(fontsize=24)  # set the font size of the y axis tick labels to 14
    plt.xlabel('Coverage', fontsize=32, fontname="serif")  # set the x axis label to "Time (s)" with font size 16
    plt.ylabel('Nucleotide Count', fontsize=32, fontname="serif")  # set the y axis label to "Voltage (V)" with font size 16
    plt.savefig(virus_name + ".histogram.png")
    plt.close()

if __name__ == '__main__':        
    virus_name = sys.argv[1]
    raw_coverage_data = sys.argv[2]
    correct_coverage_data = sys.argv[3]
    draw(virus_name, raw_coverage_data, correct_coverage_data)