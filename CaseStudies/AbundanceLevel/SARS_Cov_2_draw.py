# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-05-03 16:46:29
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-05-03 16:51:13

import matplotlib.pyplot as plt
import statistics
import numpy as np
import re
import scipy
import os
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import sys

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

    # customize the tick label font size
    plt.xticks(fontsize=24)  # set the font size of the x axis tick labels to 14
    plt.yticks(fontsize=24)  # set the font size of the y axis tick labels to 14
    plt.xlabel('Location (nt)', fontsize=32, fontname="serif")  # set the x axis label to "Time (s)" with font size 16
    plt.ylabel('Coverage', fontsize=32, fontname="serif")  # set the y axis label to "Voltage (V)" with font size 16

    # ylimlim = ax.get_ylim()

    # ss.save_file(fig, f"{virus_name}_total", 'fig')
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

    # customize the tick label font size
    plt.xticks(fontsize=24)  # set the font size of the x axis tick labels to 14
    plt.yticks(fontsize=24)  # set the font size of the y axis tick labels to 14
    plt.xlabel('Location (nt)', fontsize=32, fontname="serif")  # set the x axis label to "Time (s)" with font size 16
    plt.ylabel('Coverage\nDifference', fontsize=32, fontname="serif")  # set the y axis label to "Voltage (V)" with font size 16

    plt.savefig(virus_name + ".difference.png")
    plt.close()

    ########################################################################################################

    df_fig_neg = df_fig.iloc[18845:18895]
    df_fig_pos_color = df_fig.iloc[18889:18895]
    df_fig_neg_color = df_fig.iloc[18855:18888]

    fig, ax = plt.subplots(figsize=(10, 8))
    # customize the plot margins
    plt.subplots_adjust(top=0.5)

    fig = sns.lineplot(data=df_fig_neg, x='nt', y='count', color='#ffffff', lw=2, zorder=1)
    plt.fill_between(df_fig_neg['nt'], df_fig_neg['count'], color='#ffffff')
    fig = sns.lineplot(data=df_fig_neg_color, x='nt', y='count', color='#377eb8', lw=2, zorder=1)
    plt.fill_between(df_fig_neg_color['nt'], df_fig_neg_color['count'], color='#377eb8')
    fig = sns.lineplot(data=df_fig_pos_color, x='nt', y='count', color='#e41a1c', lw=2, zorder=1)
    plt.fill_between(df_fig_pos_color['nt'], df_fig_pos_color['count'], color='#e41a1c')

    ax.spines['left'].set_position(('outward', 8))  # set the position of the left axis
    ax.spines['bottom'].set_position(('outward', 8))  # set the position of the bottom axis
    ax.spines['right'].set_visible(False)  # hide the right axis
    ax.spines['top'].set_visible(False)  # hide the top axis

    plt.xticks(list(range(18845,18896,5)))
    plt.ylim([-2, 5])
    # ss.fig_axis_line(y=0, zorder=5, ls='--', lw=3)
    plt.axhline(y=0, zorder=5, ls='--', lw=3)

    # customize the tick label font size
    plt.xticks(fontsize=14)  # set the font size of the x axis tick labels to 14
    plt.yticks(fontsize=14)  # set the font size of the y axis tick labels to 14
    plt.xlabel('Location (nt)', fontsize=24, fontname="serif")  # set the x axis label to "Time (s)" with font size 16
    plt.ylabel('Coverage\nDifference', fontsize=24, fontname="serif")  # set the y axis label to "Voltage (V)" with font size 16

    plt.savefig('highlight_negative_SARS-Cov-2.png')

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