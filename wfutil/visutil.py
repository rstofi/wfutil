"""Utility functions for visualization.

This library contains wrapper functions to wisualize `WFdata` objects.
Mostly waterfall plots
"""

__all__ = ['quick_and_dirty_plot']

import os
import copy
import logging
import numpy as np
import numpy.ma as ma

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.patheffects import withStroke

import wfutil as wf

#=== Globals ===
#RCparams for plotting
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

matplotlib.rcParams['xtick.major.size'] = 9
matplotlib.rcParams['ytick.major.size'] = 9

matplotlib.rcParams['xtick.major.width'] = 3
matplotlib.rcParams['ytick.major.width'] = 3

matplotlib.rcParams['xtick.minor.size'] = 6
matplotlib.rcParams['ytick.minor.size'] = 6

matplotlib.rcParams['xtick.minor.width'] = 2
matplotlib.rcParams['ytick.minor.width'] = 2

matplotlib.rcParams['axes.linewidth'] = 2

plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16


#4 sampled colors from viridis
c0 = '#440154';#Purple
c1 = '#30678D';#Blue
c2 = '#35B778';#Greenish
c3 = '#FDE724';#Yellow


#*******************************************************************************
#=== Functions I ===
def quick_and_dirty_plot(WFD, plot_phase=False):
    """Quick&Dirty plot function for a `WFdata` object. It dumps the time-chan
    waterfall plot for the first polarisation

    For complex numbers it plots the amplitude by default but also the phase can
    be plotted... results in a unizorm zeros plot if the data is not complex

    NOTE: the created waterfall plots are forces to be rectangular shape

    Parameters:
    ===========
    WFD: <WFdata>
        The iunput `WFdata` object to plot
        
    plot_phase: bool, opt
        If True the code try to plot the phase of the data if possible

    Returns:
    ========
    Create a plot to show
    """

    fig, ax = plt.subplots(1, 1, figsize=(10,9.))
    
    if plot_phase:
        img = ax.matshow(np.angle(WFD.data_array[0,...])) #Phase

    else:
        img = ax.matshow(np.absolute(WFD.data_array[0,...])) #Amp

    cb = plt.colorbar(img, aspect=30, fraction=0.04975, pad=0)

    cb.ax.yaxis.get_offset_text().set_fontsize(18)
    cb.ax.tick_params(labelsize=18)
    cb.ax.tick_params(direction='in', length=6, width=2)

    cb.ax.set_ylabel(r'Data', fontsize = 18)

    ax.set_aspect('auto') #This forces the rectangular shape

    plt.xlabel(r'Channel', fontsize = 18)
    plt.ylabel(r'Time', fontsize = 18)

    plt.show()
    #plt.savefig(fname, bbox_inches='tight')

    plt.close()


#*******************************************************************************
#=== MAIN ===
if __name__ == "__main__":
    pass

    #exit()
