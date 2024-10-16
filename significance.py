import ROOT as r
import os
import numpy as np
import math
from Optimizing import OptimizeCuts
import mplhep as hep
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
import numpy as np


base_path='/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/templates_v5/templates_'

def Plot_Sig_VAE(vaes,significances,opt_vae_cut_val,sig_optimal,process):
    plt.style.use(hep.style.CMS)
    plt.figure(figsize=(12,9))

    plt.plot(vaes,significances,'o',color='k')
    plt.plot(np.array(opt_vae_cut_val),np.array(sig_optimal),'*',color='r',label='Optimal Cut',markersize=12)

    
    yTitle='Signficance'
    xTitle='VAE Loss'

    plt.xlabel(xTitle, horizontalalignment='right', x=1.0)
    plt.ylabel(yTitle,horizontalalignment='right', y=1.0)
    plt.legend(loc='upper right',ncol=2)

    plt.title(process+" Significance vs. VAE Loss Run 2")

    plt.savefig(process+'_Sig_VAE.png')
    plt.savefig(process+'_Sig_VAE.pdf')


    return

if __name__=='__main__':
    processes=['MX2200_MY400']
    years=['2016','2016APV','2017','2018']

    lowX,highX,lowY,highY=16,27,5,16 #2200, 90
    #lowX,highX,lowY,highY=16,27,16,27 #2200, 190
    #lowX,highX,lowY,highY=16,27,22,33 #2200, 250
    lowX,highX,lowY,highY=16,27,37,48 #2200, 400

    vaes=[]
    significances=[]

    for process in processes:
        optimal_vae_cut=-1
        sig_optimal=-1
        for vae_cut in range(100): # 100 vae bins
            P_Sig=OptimizeCuts(process,vae_cut,lowX,highX,lowY,highY,years)
            print('vae bin index: ',vae_cut,', significance: ',P_Sig)
            significances.append(P_Sig)

            filepath=base_path+process+'_2016.root'
            file=r.TFile.Open(filepath)
            keyname='mjj_my_vaeloss_'+process+'_2016_IR_Pass_nom'
            h=file.Get(keyname)
            cut_val=h.GetZaxis().GetBinLowEdge(vae_cut)
            vaes.append(cut_val)

            if P_Sig>sig_optimal:
                optimal_vae_cut=vae_cut
                sig_optimal=P_Sig

        opt_vae_cut_val=vaes[optimal_vae_cut]

        Plot_Sig_VAE(vaes,significances,opt_vae_cut_val,sig_optimal,process)
