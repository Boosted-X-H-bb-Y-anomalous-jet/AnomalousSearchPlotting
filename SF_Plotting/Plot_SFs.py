import mplhep as hep
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
import numpy as np

def plot(years,Mxs,SFs,Uncs_up,Uncs_down,colors):

    for i in range(len(years)):
        plt.style.use(hep.style.CMS)
        plt.figure(figsize=(12,9))

        plt.errorbar(Mxs,SFs[i],yerr=(Uncs_down[i],Uncs_up[i]),fmt='o',color=colors[i],label=years[i])

    
        yTitle='Scale Factor'
        xTitle='MX'
        plt.ylim(ymin=0.5, ymax=1.5)

        plt.xlabel(xTitle, horizontalalignment='right', x=1.0)
        plt.ylabel(yTitle,horizontalalignment='right', y=1.0)
        plt.legend(loc='upper right',ncol=2)

        plt.savefig('vae_SFs'+years[i]+'.pdf')

    return

if __name__ == '__main__':
    with open("SFs.txt", 'r') as file:
        lines=file.readlines()
    
    years=['2016','2016APV','2017','2018']
    colors=['k','g','b','r']

    SFs2016=[]
    SFs2016APV=[]
    SFs2017=[]
    SFs2018=[]
    MXs=[]

    Uncs_up=[[],[],[],[]]
    Uncs_down=[[],[],[],[]]
    for line in lines:
        content = line.split(",")
        if '2016APV' in content[0]:
            SFs2016APV.append(float(content[2]))
            MXs.append(float(content[1][2:6]))
            unc_up=float(content[3])
            unc_down=float(content[4])
            Uncs_up[1].append(unc_up)
            Uncs_down[1].append(unc_down)
        elif '2016' in content[0]:
            SFs2016.append(float(content[2]))
            unc_up=float(content[3])
            unc_down=float(content[4])
            Uncs_up[0].append(unc_up)
            Uncs_down[0].append(unc_down)
        elif '2017' in content[0]:
            SFs2017.append(float(content[2]))
            unc_up=float(content[3])
            unc_down=float(content[4])
            Uncs_up[2].append(unc_up)
            Uncs_down[2].append(unc_down)
        elif '2018' in content[0]:
            SFs2018.append(float(content[2]))
            unc_up=float(content[3])
            unc_down=float(content[4])
            Uncs_up[3].append(unc_up)
            Uncs_down[3].append(unc_down)
    
    SFs=[SFs2016,SFs2016APV,SFs2017,SFs2018]
    plot(years,MXs,SFs,Uncs_up,Uncs_down,colors)
            



