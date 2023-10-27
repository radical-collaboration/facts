import glob
import os
import shutil
import re
import fnmatch
#
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
#
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.basemap import Basemap
import cartopy
#
PD=os.getcwd(); 

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create a list of file paths.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def filePATH(ssp,path): 
    
    # files=[path+f'wf1e/{ssp}/coupling.{ssp}.emuAIS.emulandice.AIS_globalsl_quantiles.nc',
    #        path+f'wf2e/{ssp}/coupling.{ssp}.larmip.larmip.AIS_globalsl_quantiles.nc',
    #        path+f'wf3e/{ssp}/coupling.{ssp}.deconto21.deconto21.AIS_AIS_globalsl_quantiles.nc',
    #        path+f'wf4/{ssp}/coupling.{ssp}.bamber19.bamber19.icesheets_AIS_globalsl_quantiles.nc']
    # return files

    files=[path+f'wf1e/{ssp}/coupling.{ssp}.total.workflow.wf1e.global_quantiles.nc',
           path+f'wf2e/{ssp}/coupling.{ssp}.total.workflow.wf2e.global_quantiles.nc',
           path+f'wf3e/{ssp}/coupling.{ssp}.total.workflow.wf3e.global_quantiles.nc',
           path+f'wf4/{ssp}/coupling.{ssp}.total.workflow.wf4.global_quantiles.nc']
    return files

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def nc2var(files,loc):
    sls     = np.empty((0, 0))
    Qvals   = np.empty((0, 0))
    label   = []
    #
    for fi0,fi1 in enumerate(files):
        #
        lab=fi1.split('/')[-1].split('.')[4]
        # lab=fi1.split('/')[-1].split('.')[2]
        # ................................................................
        # Exract Data.
        dataset = xr.open_dataset(fi1)
        time    = dataset['years'].values
        subt    = np.where(time == 2100)[0][0]
        quant   = dataset['quantiles'].values
        qlevs = np.arange(0.01, 1, 0.01)
        subq = np.where((np.around(qlevs,decimals=2) == 0.17) | (np.around(qlevs,decimals=2) == 0.83))[0]
        #
        sl      = dataset['sea_level_change'].values[:,subt,loc].T
        qvals   = np.percentile(sl, qlevs * 100)
        if sls.size == 0: 
            sls = sl
            label=lab
            Qvals=qvals
        else: 
            sls=np.vstack([sls,sl])
            label=np.vstack([label,lab])
            Qvals=np.vstack([Qvals,qvals])
    
    return sls, Qvals, label, quant, qlevs, subq 



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
colors = {
    'ssp119' : np.array([0, 173, 207]) / 255,
    'ssp126' : np.array([23, 60, 102]) / 255,
    'ssp245' : np.array([247, 148, 32]) / 255,
    'ssp370' : np.array([231, 29, 37]) / 255,
    'ssp585' : np.array([149, 27, 30]) / 255
}
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def ax_properties(ax,axes,ss1,fnt,xlab,ylab, x_min, x_max, y_min=None, y_max=None):
    ax.set_xlabel(xlab, fontsize=fnt+2)    
    ax.set_ylabel(ylab, fontsize=fnt+2)    
    #
    ax.set_xlim(x_min, x_max)  
    x_ticks= np.around(np.arange(x_min, x_max+0.1, .5), decimals=1)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticks,fontsize=fnt, rotation=0)
    #
    if y_min is not None and y_max is not None:
        ax.set_ylim(y_min, y_max)
        if ax==axes[0,0] or ax==axes[0,1]:
            y_ticks = np.around(np.arange(y_min, y_max, 1), decimals=1)
        else: 
            y_ticks = np.around(np.arange(y_min, y_max, 0.2), decimals=1)
        #
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_ticks, fontsize=fnt, rotation=0) 
    #
    ax.tick_params(direction='in', length=3.5, width=1, axis='both')
    #
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1)
    #
    # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # Title and legend.
    if ax==axes[0,0] or ax==axes[0,1]:
        title=ss1[:4].upper()+'-'+ss1[4]+'.'+ss1[5];  ax.set_title(title,fontsize=fnt+10)
    if ax==axes[0,0]: ax.legend(fontsize=fnt+1) 
    #
    # 0 Horizontal line. 
    if ax==axes[1,0] or ax==axes[1,1]:
        ax.axhline(0, color='black', linestyle='-', linewidth=1)
        # ax.text(0.65, 0.15, '17th-83rd percentile ranges', transform=ax.transAxes, verticalalignment='top', fontsize='5', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        ax.text(0.5, 0.15, '17th-83rd percentile ranges', transform=ax.transAxes, verticalalignment='top', fontsize=7)
    #    
    # 0 Horizontal line. 
    if ax==axes[2,0] or ax==axes[2,1]:
        ax.plot([0, 4], [0, 0], color='black', linestyle='-')
        # if ax==axes[2,0]: ax.legend(fontsize=fnt+1) 
    return(None)
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def plot(loc,ssp,path):
    #
    # fig, axes = plt.subplots(3, 2, figsize=(10*2, 3*3)); 
    # fig, axes = plt.subplots(3, 2, figsize=(12, 8))
    fig, axes = plt.subplots(3, 2, figsize=(10, 8))
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
    #
    for ss0,ss1 in enumerate(ssp):
        #
        files=filePATH(ss1,path)   
        sls, Qvals, label, quant, qlevs, subq = nc2var(files,loc)
        #
        fnt=8
        lws = [2, 2, 1, 1] * 2
        color_list = list(colors.values())[:5]; 
        hp = []  
        binwidth = 50
        #
        for sss in range(4):
            # 
            current_color = color_list[sss]
            if sss < 2: linewidth = 3.5  
            else: linewidth = 1.5
            # ==================================================================================
            # PLOT  Pnl 1
            ax= axes[0,ss0]
            # =================================================================================
            # n, edges = np.histogram(sls[sss, :], np.arange(0, 4001, binwidth), density=True)
            bin_edges = np.arange(0, 4000 + binwidth, binwidth)
            n, edges = np.histogram(sls[sss, :], bins=bin_edges)
            #
            hist_sum = np.sum(n)
            scaling_factor = binwidth/1000
            n = n / hist_sum / scaling_factor
            #
            xx=.5*(edges[1:]+edges[:-1])/1000
            ax.plot(xx, n, color=current_color,label=label[sss, 0])
            # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
            # Axis Properties
            xlab = 'GMSL rise in 2100 (m)'
            ylab = 'Probability density' 
            x_min = 0;  x_max = 3.5
            y_min = 0;  y_max = 5
            #
            ax_properties(ax,axes,ss1,fnt,xlab,ylab, x_min, x_max,y_min,y_max)
            #
            #
            # ====================================================================================
            # PLOT  Pnl 2
            ax= axes[1,ss0]
            # =====================================================================================
            #
            xx=sls[sss,:]/1000
            xx1=sls[sss,subq]/1000
            yy=quant
            yy1 = [-0.05 * sss-0.05] * 2
            ax.plot(xx, yy, color=current_color, linewidth=linewidth,label=label[sss, 0])
            ax.plot(xx1,yy1,color=current_color, linewidth=linewidth)
            # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
            # Axis Properties
            xlab = 'GMSL rise in 2100 (m)'
            ylab = 'CDF quantile' 
            x_min = -0.0;  x_max = 3.5
            y_min = -0.25;  y_max = 1
            #
            ax_properties(ax,axes,ss1,fnt,xlab,ylab, x_min, x_max,y_min,y_max)
            #            
            #
            # ====================================================================================================================================
            # PLOT  Pnl 3
            ax= axes[2,ss0]
            # ====================================================================================================================================
            #
            # P-box 
            pbox1l = np.min(Qvals[:2, :], axis=0)
            pbox1h = np.max(Qvals[:2, :], axis=0)
            pbox2l = np.min(Qvals, axis=0)
            pbox2h = np.max(Qvals, axis=0)
            #
            # ax.fill_betweenx(np.concatenate([qlevs,qlevs[::-1]]), np.concatenate([pbox2l, pbox2h[::-1]]) / 1000,  color=[0.8, 0.8, 1],label='Med. conf.')
            # ax.fill_betweenx(np.concatenate([qlevs,qlevs[::-1]]), np.concatenate([pbox1l, pbox1h[::-1]]) / 1000,  color=[0.4, 0.4, 1],label='Low. conf.')
            ax.fill_betweenx(np.concatenate([qlevs,qlevs[::-1]]), np.concatenate([pbox2l, pbox2h[::-1]]) / 1000,  color=[0.8, 0.8, 1])
            ax.fill_betweenx(np.concatenate([qlevs,qlevs[::-1]]), np.concatenate([pbox1l, pbox1h[::-1]]) / 1000,  color=[0.4, 0.4, 1])
            #
            ax.plot([pbox2l[subq[0]] / 1000, pbox2h[subq[1]] / 1000], [-0.1, -0.1], color=[0.8, 0.8, 1], linewidth=6)
            ax.plot([pbox1l[subq[0]] / 1000, pbox1h[subq[1]] / 1000], [-0.1, -0.1], color=[0.4, 0.4, 1], linewidth=6)
            # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
            # Axis Properties
            xlab = 'GMSL rise in 2100 (m)'
            ylab = 'p-box quantile' 
            x_min = -0.0;  x_max = 3.5
            y_min = -0.25;  y_max = 1
            #
            ax_properties(ax,axes,ss1,fnt,xlab,ylab, x_min, x_max,y_min,y_max)
    plt.show()