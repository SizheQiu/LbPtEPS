from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

def volcano( table, lfc_col, pv_col, lfc_cutoff, pv_cutoff, ax):
    lfc_list = list(table[lfc_col])
    mlg10pv = [ -np.log10(x) for x in list(table[pv_col]) ]
    x_cut1 = lfc_cutoff; x_cut2 = - lfc_cutoff;
    y_cut = -np.log10( pv_cutoff )
    color_list = []
    for i in range(len(lfc_list) ):
        if mlg10pv[i] > y_cut:
            if lfc_list[i] > x_cut1:
                color_list.append('tab:red')
            elif lfc_list[i] < x_cut2:
                color_list.append('tab:blue')
            else:
                color_list.append('grey')
        else:
            color_list.append('grey')
    ax.scatter(lfc_list, mlg10pv, c=color_list, marker='o',linewidth=0.5,edgecolor='black',alpha=0.5)
    ax.axhline(y=y_cut, color='grey', linestyle='--')
    ax.axvline(x=x_cut1, color='grey', linestyle='--')
    ax.axvline(x=x_cut2, color='grey', linestyle='--')
    return ax

def annot_volcano( table, name_dict, lfc_col, pv_col, lfc_cutoff, pv_cutoff, fontsize, ax ):
    x_cut1 = lfc_cutoff; x_cut2 = - lfc_cutoff;
    y_cut = -np.log10( pv_cutoff )
    
    sig_table = table[ ( (table[lfc_col]>x_cut1) | (table[lfc_col]<x_cut2) ) & (table[pv_col]<pv_cutoff)  ]
    sig_table = sig_table.reset_index().drop(['index'],axis=1)
    x_list = list(sig_table[lfc_col])
    y_list = [ -np.log10(x) for x in list(sig_table[pv_col]) ]
    names = [ name_dict[gid] for gid in list(sig_table['ID']) ]
    
    for i in range(len(x_list)):
        ax.text(1.01*x_list[i], y_list[i], s=names[i], fontsize=fontsize)
    return ax





