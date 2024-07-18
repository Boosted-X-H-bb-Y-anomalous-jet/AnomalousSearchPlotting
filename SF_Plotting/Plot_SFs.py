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

        plt.xlabel(xTitle, horizontalalignment='right', x=1.0)
        plt.ylabel(yTitle,horizontalalignment='right', y=1.0)
        plt.legend(loc='upper right',ncol=2)

        plt.savefig('SFs'+years[i]+'.png')

    return

if __name__ == '__main__':
    with open("SFs.txt", 'r') as file:
        lines=file.readlines()

    for i in range(len(lines)):
        lines[i]=lines[i].split(', ')
    
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
        if '2016APV' in line[0]:
            SFs2016APV.append(float(line[3][0:4]))
            MXs.append(float(line[0][10:14]))
            unc_up=float(line[1][17:21])/float(line[0][-4:])
            unc_down=float(line[1][23:])/float(line[0][-4:])
            Uncs_up[1].append(unc_up)
            Uncs_down[1].append(unc_down)
        elif '2016' in line[0]:
            SFs2016.append(float(line[3][0:4]))
            unc_up=float(line[1][17:21])/float(line[0][-4:])
            unc_down=float(line[1][23:])/float(line[0][-4:])
            Uncs_up[0].append(unc_up)
            Uncs_down[0].append(unc_down)
        elif '2017' in line[0]:
            SFs2017.append(float(line[3][0:4]))
            unc_up=float(line[1][17:21])/float(line[0][-4:])
            unc_down=float(line[1][23:])/float(line[0][-4:])
            Uncs_up[2].append(unc_up)
            Uncs_down[2].append(unc_down)
        elif '2018' in line[0]:
            SFs2018.append(float(line[3][0:4]))
            unc_up=float(line[1][17:21])/float(line[0][-4:])
            unc_down=float(line[1][23:])/float(line[0][-4:])
            Uncs_up[3].append(unc_up)
            Uncs_down[3].append(unc_down)
    
    SFs=[SFs2016,SFs2016APV,SFs2017,SFs2018]
    plot(years,MXs,SFs,Uncs_up,Uncs_down,colors)
            



