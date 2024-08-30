import ROOT as r
import os
import numpy as np
import math
import mplhep as hep
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
import numpy as np
import matplotlib.cm as cm

base_path='/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/CMSSW_14_1_0_pre4/src/AnomalousSearchFits/SR_run2_new_signals/'
MXs=[1400,1600,1800,2200,2600,3000]
MYs=[90,125,190,250,300,400]

fb_conversion=5

def get_limit(filepath):

    file=r.TFile.Open(filepath)
    keyname='limit'

    tree=file.Get(keyname)

    tree.GetEntry(2)

    limit=tree.limit*fb_conversion

    return limit

def Plot(Limits):

    plt.style.use(hep.style.CMS)

    x_edges=MXs+[3400]
    y_edges=MYs+[425]

    plt.pcolormesh(x_edges, y_edges, Limits, shading='auto', cmap='viridis')

    plt.colorbar() #label='Efficiency')

    plt.ylabel(r'$M_{Y}$', fontsize=24)
    plt.xlabel(r'$M_{X}$', fontsize=24)

    plt.title("Mass vs Expected Limits")
    plt.savefig('Limit_vs_Mass.png')

def MakePlots():
    
    limits=[]
    for MY in MYs:
        limits.append([])    
        for MX in MXs:    
            filepath=base_path+'MX'+str(MX)+'_MY'+str(MY)+'-0_area/higgsCombineTest.AsymptoticLimits.mH120.root'
            if os.path.exists(filepath):
                print('Getting limit with MX= '+str(MX)+' MY= '+str(MY))
                limit=get_limit(filepath)
                limits[-1].append(limit)
            else:
                print("No file found")
    print(limits)
    Plot(limits)
    return
    

if __name__=='__main__':
    
    MakePlots()

