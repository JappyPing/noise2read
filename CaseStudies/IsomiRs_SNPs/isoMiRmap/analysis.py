# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-05-03 15:56:56
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-05-03 16:19:33

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from matplotlib.lines import Line2D
import sys

def add_mature_miRNA(df_data):
    data = df_data.copy()
    mirna_list = []
    if 'Hairpin locations (comma delimited)' in data.columns.tolist():
        target_col = 'Hairpin locations (comma delimited)'
    else:
        target_col = 'Hairpin(s)'
    for i in data[target_col].tolist():
        try:
            mirna = i.split('|')[1].split('&')[0]
        except:
            mirna = 'None'
        mirna_list.append(mirna)
    data.insert(loc = 0,
                column = 'Mature miRNA',
                value = mirna_list)
    return data

def df_col_log10(Df_In, Col_List=[]):
    '''
    ### Convert dataframe columns with log10 transform
    '''
    df_data = Df_In.copy()
    if isinstance(Col_List, list) is True:
        col_list = Col_List[:]
    else:
        col_list = [Col_List]
    for col in col_list:
        df_data[col] = list(np.log10(df_data[col].tolist()))
    return df_data

# raw_dir = "../isomimap_result/raw_result/"
# correct_dir = "../isomimap_result/correction/"

def analysis(raw_dir, correct_dir, out_dir):

    file_names = os.listdir(raw_dir)
    fn_lst =  [name for name in file_names if name.endswith('.txt') and "countsmeta.txt" not in name and "IsoMiRmap_v5" in name]

    for fn in fn_lst:
        raw_fn = os.path.join(raw_dir, fn)
        base_name = fn.split('-IsoMiRmap_v5')
        cor_fn = base_name[0].split('.')[0] + "_corrected-IsoMiRmap_v5" + base_name[1]
        correct_fn = os.path.join(correct_dir, cor_fn)
        
        print(correct_fn, raw_fn)

        data_cor = pd.read_csv(correct_fn, sep='\t', skiprows=6)
        data_raw = pd.read_csv(raw_fn, sep='\t', skiprows=6)
        
        data_cor = add_mature_miRNA(data_cor)
        data_raw = add_mature_miRNA(data_raw)

        mirna_list = list(set(data_raw['Mature miRNA']))
        
        res = {}

        for loc, mirna in enumerate(mirna_list):

            df_t_cor = data_cor[data_cor['Mature miRNA'] == mirna]
            df_t_cor = df_t_cor.reset_index(drop=True)
            read_count_l_cor = df_t_cor['Unnormalized read counts'].tolist()

            df_t_raw = data_raw[data_raw['Mature miRNA'] == mirna]
            df_t_raw = df_t_raw.reset_index(drop=True)
            read_count_l_raw = df_t_raw['Unnormalized read counts'].tolist()

            res[loc] = {'Mature miRNA'          :mirna,
                        'n_isoform_cor'         :df_t_cor.shape[0],
                        'read_count_cor'        :sum(read_count_l_cor),
                        'read_count_cor_log10'  :sum(read_count_l_cor),
                        'n_isoform_raw'         :df_t_raw.shape[0],
                        'read_count_raw'        :sum(read_count_l_raw),
                        'read_count_raw_log10'  :sum(read_count_l_raw),}
        
        df_res = pd.DataFrame.from_dict(res, 'index')
        df_res = df_col_log10(df_res, ['read_count_raw_log10', 'read_count_cor_log10'])
        df_res['n_isoform_difference'] = list(np.array(df_res['n_isoform_raw']) - np.array(df_res['n_isoform_cor']))
        df_res['read_count_difference'] = list(np.array(df_res['read_count_raw']) - np.array(df_res['read_count_cor']))
        df_res['n_isoform_difference_percent'] = list(np.array((df_res['n_isoform_raw']) - np.array(df_res['n_isoform_cor'])) / np.array(df_res['n_isoform_raw']))
        df_res['read_count_difference_percent'] = list((np.array(df_res['read_count_raw']) - np.array(df_res['read_count_cor'])) / np.array(df_res['read_count_raw']))

        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(1, 1, 1)
        
        fig = sns.kdeplot(data=df_res, x='n_isoform_raw', y='read_count_raw_log10', color='#000000', 
                            alpha=0.4, fill=True)
        fig = sns.kdeplot(data=df_res, x='n_isoform_cor', y='read_count_cor_log10', color='#e41a1c', 
                            alpha=0.4, fill=True)
                            
        fig = sns.scatterplot(data=df_res, x='n_isoform_raw', y='read_count_raw_log10', 
                                s=100, color='#000000')#, hue="Identical miRNA",  palette=palette, ,legend='auto'
        fig = sns.scatterplot(data=df_res, x='n_isoform_cor', y='read_count_cor_log10', 
                                s=100, color='#e41a1c')#legend='auto'
        
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)

        plt.ylabel('log\u2081\u2080 (IsomiRs count)', fontname="serif", fontsize=28)
        plt.xlabel("Identical isomiRs' Number", fontname="serif", fontsize=28)

        
        ax.spines['left'].set_position(('outward', 8))  # set the position of the left axis
        ax.spines['bottom'].set_position(('outward', 8))  # set the position of the bottom axis
        ax.spines['right'].set_visible(False)  # hide the right axis
        ax.spines['top'].set_visible(False)  # hide the top axis

        custom = [Line2D([], [], marker='o', color='#000000', linestyle='None'),
                Line2D([], [], marker='o', color='#e41a1c', linestyle='None')]

        lgnd = plt.legend(custom, ['Original', 'Corrected'], loc='lower right', fontsize=24, frameon=False)#prop={'size': 22}

        plt.savefig(os.path.join(out_dir, fn+'.png'), dpi=150, format='png')

if __name__ == '__main__':
    raw_dir = sys.argv[1]
    correct_dir = sys.argv[2] 
    out_dir = sys.argv[3]   
    analysis(raw_dir, correct_dir, out_dir)     
