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

base_path='/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/templates_v8/templates_'
MXs=[1400,1600,1800,2000,2200,2600,3000]#[1200,1400,1600,1800,2000,2200,2400,2500,2600,2800,3000,3500,4000]
MYs=[90,125,190,250,300,400]

generated_events=5*138

def get_efficiency(MX,MY):
    SR_Pass_events=0
    filepath=base_path+'MX'+MX+'_MY'+MY+'_run2.root'

    file=r.TFile.Open(filepath)
    keyname='mjj_my_MX'+MX+'_MY'+MY+'_SR_Pass_nom'

    h=file.Get(keyname)

    SR_Pass_events+=h.Integral()


    efficiency=SR_Pass_events/generated_events

    return efficiency

def Plot(Efficiencies, title="Selection Efficiency vs Mass", outputFile="plots/Selection_Efficiency_vs_Mass"):

    plt.style.use(hep.style.CMS)

    x_edges = MXs + [3400]
    y_edges = MYs + [425]
    z_min = 0.
    z_max = 1.0

    Efficiencies = np.array(Efficiencies)

    mesh = plt.pcolormesh(x_edges, y_edges, Efficiencies, shading='auto', cmap='viridis', vmin=z_min, vmax=z_max)
    plt.colorbar(mesh, label='Efficiency')

    for i in range(Efficiencies.shape[0]):
        for j in range(Efficiencies.shape[1]):
            plt.text(x_edges[j] + 0.5 * (x_edges[j+1] - x_edges[j]), 
                     y_edges[i] + 0.5 * (y_edges[i+1] - y_edges[i]), 
                     f'{Efficiencies[i, j]:.2f}', 
                     ha='center', va='center', color='white', fontsize=18, fontweight='bold')

    plt.ylabel(r'$M_{Y}$', fontsize=24)
    plt.xlabel(r'$M_{X}$', fontsize=24)
    #plt.title(title)

    lumi = 138
    lumiText = str(lumi)+ " $fb^{-1} (13 TeV)$"    
    hep.cms.lumitext(lumiText)
    hep.cms.text("Preliminary",loc=0)

    plt.savefig(f'{outputFile}.png')
    plt.savefig(f'{outputFile}.pdf')
    
    plt.cla()
    plt.clf()


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
    
def get_efficiencies_from_pkl(process):
    import pickle
    lumis = {"2016":16.8,"2016APV":19.5,"2017":41.5,"2018":59.8}
    with open("efficiencies.pkl", "rb") as f:
        gen_reco_efficiencies = pickle.load(f)
    reco = sum(gen_reco_efficiencies["reco"][process][year] * lumis[year] for year in lumis)/sum(lumis[year] for year in lumis)
    gen = sum(gen_reco_efficiencies["gen"][process][year] * lumis[year] for year in lumis) / sum(lumis[year] for year in lumis)
    return reco, gen


def MakeGenRecoPlots():
    efficiencies_reco=[]
    efficiencies_gen=[]
    for MY in MYs:
        efficiencies_reco.append([])    
        efficiencies_gen.append([])    
        for MX in MXs:    
            reco, gen =get_efficiencies_from_pkl(f"MX{MX}_MY{MY}")
            efficiencies_reco[-1].append(reco)
            efficiencies_gen[-1].append(gen)

    
    Plot(efficiencies_gen,title="Generator-level efficiency",outputFile="plots/gen_efficiency")
    Plot(efficiencies_reco,title="Reconstruction-level efficiency",outputFile="plots/reco_efficiency")
    return

def MakeDelphesPlot():
    efficiencies = []
    import csv

    with open("delphes_anomaly_tagging_efficiencies.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        data = {row[0]: float(row[2]) for row in reader}

    for MY in MYs:
        efficiencies.append([])
        for MX in MXs:
            key = f"MX{MX}_MY{MY}"
            if key in data:
                efficiencies[-1].append(data[key])
            else:
                efficiencies[-1].append(0.0)

    Plot(efficiencies, title="Delphes-based reduced selection efficiency", outputFile="plots/delphes_efficiency")
    return

if __name__=='__main__':
    
    #MakePlots()
    #MakeGenRecoPlots()
    MakeDelphesPlot()
