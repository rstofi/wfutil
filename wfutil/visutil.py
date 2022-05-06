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

from matplotlib.backends.backend_pdf import PdfPages

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
#=== Functions ===
def quick_and_dirty_plot(WFD,
                        plot_phase=False,
                        savefig=False,
                        fname='./test_flag_fraction_matrix.png',
                        ptitle=None,
                        ctitle=None):
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
    
    if ptitle != None:
        fig.suptitle('{0:s}'.format(ptitle), fontsize = 18)

    if plot_phase:
        data = np.angle(WFD.data_array[0,...])
        #img = ax.matshow(np.angle(WFD.data_array[0,...])) #Phase

    else:
        data = np.absolute(WFD.data_array[0,...])
        #img = ax.matshow(np.absolute(WFD.data_array[0,...])) #Amp

    #Plot the data

    #TO DO add an option to use the masked data as the min max color not only the visible values

    img = ax.matshow(data, cmap='viridis', vmin=np.min(data), vmax=np.max(data))

    cb = plt.colorbar(img, aspect=30, fraction=0.04975, pad=0)

    cb.ax.yaxis.get_offset_text().set_fontsize(18)
    cb.ax.tick_params(labelsize=18)
    cb.ax.tick_params(direction='in', length=6, width=2)

    if ctitle != None:
        cb.ax.set_ylabel('{0:s}'.format(ctitle), fontsize = 18)
    else:
        cb.ax.set_ylabel(r'Data', fontsize = 18)

    ax.set_aspect('auto') #This forces the rectangular shape

    plt.xlabel(r'Channel', fontsize = 18)
    plt.ylabel(r'Time', fontsize = 18)
    
    if savefig:
        plt.savefig(fname, bbox_inches='tight')
    else:
        plt.show()

    plt.close()

def plot_flag_fraction_for_sant(WFD, 
                                save_fig=False,
                                fname='./test_flag_fraction_matrix.png',
                                ptitle='Baseline flag fraction plot for example antenna',
                                recalculate_mask=False,
                                flag_treshold=0.66):
    """Create diagnostics plots for the crosscol flag fractions or the derived flags
    """
    #Get polarisation axis
    N_pol = np.size(WFD.pol_array)

    #If new mask needs to be calculated
    if recalculate_mask:
        mask = np.zeros(np.shape(WFD.data_array), dtype=bool)
        mask[WFD.data_array > flag_treshold] = True

        #Apply new mask
        WFD.apply_mask(mask)

    #Set up figure based on the number of polarisations
    fig, ax = plt.subplots(N_pol,1,figsize=(15,1+(N_pol * 5)))
    
    #Add Suptitle
    fig.suptitle('{0:s}'.format(ptitle), fontsize = 18)
    

    for i in range(0,N_pol):
        ax[i].set_title('{}'.format(WFD.pol_array[i]), fontsize = 16)
        pdata = ax[i].imshow(WFD.data_array[i,...], aspect='auto',
                            vmin=0.0, vmax=1.0)
        ax[i].set_ylabel('Time', fontsize = 16)

    
        if i == (N_pol - 1) :
            ax[i].set_xlabel('Channel',fontsize = 16)
    
    #Add colorbar
    #if not only_flags:
    #cax = fig.colorbar(pdata, ax=ax, fraction=0.046, pad=0.04)
    cax = fig.colorbar(pdata, ax=ax, aspect=40)
    cax.set_label('Flag fraction', fontsize=18)
    
    if save_fig:
        plt.savefig(fname, bbox_inches="tight")
        plt.close()
    else:
        #plt.tight_layout()
        plt.show()

#*******************************************************************************
#=== MAIN ===
if __name__ == "__main__":
    pass

    #exit()
