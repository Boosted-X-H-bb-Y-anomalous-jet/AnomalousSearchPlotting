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

base_path='/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/templates_v4/templates_'
MXs=[1400,1600,1800,2200,2600,3000]#[1200,1400,1600,1800,2000,2200,2400,2500,2600,2800,3000,3500,4000]
MYs=[90,125,190,250,300,400]
years=['2016','2016APV','2017','2018']

generated_events=5*138

def get_efficiency(MX,MY):
    SR_Pass_events=0
    for year in years:
        filepath=base_path+'MX'+MX+'_MY'+MY+'_'+year+'.root'

        file=r.TFile.Open(filepath)
        keyname='mjj_my_MX'+MX+'_MY'+MY+'_'+year+'_SR_Pass_nom'

        h=file.Get(keyname)

        SR_Pass_events+=h.Integral()


    efficiency=SR_Pass_events/generated_events

    return efficiency

def Plot(Efficiencies):

    plt.style.use(hep.style.CMS)

    x_edges=MXs+[3400]
    y_edges=MYs+[425]

    plt.pcolormesh(x_edges, y_edges, Efficiencies, shading='auto', cmap='viridis')
    #plt.grid(True)

    #fig, ax = plt.subplots()

    #MX,MY=np.meshgrid(MXs, MYs)

    #im = ax.imshow(Efficiencies, interpolation='bilinear', origin='lower', cmap=cm.viridis, extent=(MXs[0],MXs[-1],MYs[0],MYs[-1]))
    #ConPlot = ax.contour(MX,MY,Efficiencies,[-5],colors='k', extent=(MXs[0],MXs[-1],MYs[0],MYs[-1]))
    
    #CB = fig.colorbar(im, shrink=0.8)
    #CB.ax.set_title(r'Efficiency',fontsize=14)

    plt.colorbar(label='Efficiency')

    plt.ylabel(r'$M_{Y}$', fontsize=24)
    plt.xlabel(r'$M_{X}$', fontsize=24)

    #plt.aspect((MXs[-1]-MXs[0])/(MYs[-1]-MYs[0]))

    plt.title("Selection Efficiency vs Mass")
    plt.savefig('Selection_Efficiency_vs_Mass.png')

    """plt.style.use(hep.style.CMS)
    plt.figure(figsize=(12,9))

    plt.plot(Masses,Efficiencies,'o',color='r')
    
    yTitle='Efficiency'
    xTitle=MXorMY

    plt.xlabel(xTitle, horizontalalignment='right', x=1.0)
    plt.ylabel(yTitle,horizontalalignment='right', y=1.0)
    plt.legend(loc='upper right',ncol=2)\
    
    plt.title("MX vs. Efficiency with MY= "+ConstantVal)
    plt.savefig('MX_vs_Eff_MY'+ConstantVal+'.png')
    print('Saving plot as MX_vs_Eff_MY'+ConstantVal+'.png')"""

def MakePlots():
    
    efficiencies=[]
    for MY in MYs:
        efficiencies.append([])    
        for MX in MXs:    
            filepathcheck=base_path+'MX'+str(MX)+'_MY'+str(MY)+'_2016.root'
            if os.path.exists(filepathcheck):
                print('Calculating Efficiency with MX= '+str(MX)+' MY= '+str(MY))
                efficiency=get_efficiency(str(MX),str(MY))
                efficiencies[-1].append(efficiency)
            else:
                print("No file found")
    
    Plot(efficiencies)
    return
    

if __name__=='__main__':
    
    MakePlots()

