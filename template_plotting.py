import matplotlib
matplotlib.use('Agg')

import ROOT as r
from optparse import OptionParser
from time import sleep
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mplhep as hep
from root_numpy import hist2array
import ctypes
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import os 
#import misc.plotRPFWithUnc as plotRPFWithUnc

matplotlib.use('Agg')
r.gROOT.SetBatch(True)
r.gStyle.SetOptFit(111)

base_path='/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/templates_v2/'

log=True

def plot(histosData,edgesData,colorsData,labelsData,histosSig,edgesSig,colorsSig,labelsSig,Pass):
    QCD=[]
    #QCDedges=[]
    nonQCDhistos=[]
    #nonQCDedges=[]
    nonQCDcolors=[]
    nonQCDlabels=[]
    outlinecolors=[]

    
    for i in range(len(histosSig)):
        if 'QCD' in labelsSig[i]:
            QCD.append(histosSig[i])
            #QCDedges.append(edgesSig[i])
        else:
            nonQCDhistos.append(histosSig[i])
            #nonQCDedges.append(edgesSig[i])
            nonQCDcolors.append(colorsSig[i])
            nonQCDlabels.append(labelsSig[i])
            outlinecolors.append('dimgray')
            print(labelsSig[i])
    


    QCDhistos=[]
    for j in range(len(QCD[0])):
        val=0
        for histo in QCD:
            val+=histo[j]
        QCDhistos.append(val)

    if histosData!=None:
        errors=[]
        errorsCenters=[]
        for i in range(len(histosData[0])):
            errors.append(np.sqrt(histosData[0][i]))
            errorsCenters.append((edgesData[0][i]+edgesData[0][i+1])/2)
        
        Total_QCD_Sim=np.sum(QCDhistos)
        Total_Data=np.sum(histosData)
        Total_QCD_Data=Total_Data
        for data_set in nonQCDhistos:
            Total_QCD_Data=Total_QCD_Data-np.sum(data_set)
        ratio=Total_QCD_Data/Total_QCD_Sim
    
        for i in range(len(QCDhistos)):
            QCDhistos[i]=QCDhistos[i]*ratio

    
    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()

    hep.histplot(QCDhistos,edgesSig[0],stack=False,ax=ax,label=['QCD'],linewidth=3,histtype="fill",color=['burlywood'])
    hep.histplot(QCDhistos,edgesSig[0],stack=False,ax=ax,linewidth=3,histtype="step",color=['dimgray'])
    hep.histplot(nonQCDhistos,edgesSig[0],stack=False,ax=ax,label=nonQCDlabels,linewidth=3,histtype="fill",color=nonQCDcolors)
    hep.histplot(nonQCDhistos,edgesSig[0],stack=False,ax=ax,linewidth=3,histtype="step",color=outlinecolors)
    if histosData!=None:
        ax.errorbar(errorsCenters,histosData[0],yerr=errors,fmt='o', color='black')
        hep.histplot(histosData,edgesData[0],stack=False,ax=ax,label=labelsData,linewidth=3,histtype="errorbar",color=colorsData)

    if(log):
        ax.set_yscale("log")

    yTitle='Events'
    xTitle='$M_{jj}$ [GeV]'

    ax.set_ylabel(yTitle)
    ax.set_xlabel(xTitle)
    plt.xlabel(xTitle, horizontalalignment='right', x=1.0)
    plt.ylabel(yTitle,horizontalalignment='right', y=1.0)
    if histosData!=None:
        if Pass:
            regionText = "CR Pass"
            plt.ylim(1., 10**3)
            plt.xlim(1000, 3000)
        else:
            regionText = "CR Fail"
            plt.ylim(10., 10**6)
            plt.xlim(1000, 3000)
    else:
        if Pass:
            regionText = "SR Pass"
            plt.ylim(1., 10**3)
            plt.xlim(1000, 3000)
        else:
            regionText = "SR Fail"
            plt.ylim(10., 10**6)
            plt.xlim(1000, 3000)

    plt.text(0.8, 0.8, regionText,horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
    hep.cms.text("Simulation WiP",loc=0)
    #lumiText =  $XYZ fb^{-1} (13 TeV)$"    
    lumiText = f"{year} (13 TeV)"    
    hep.cms.lumitext(lumiText)
    plt.legend(loc='upper right',ncol=2)#loc="best",ncol=2)#loc = 'best'    
    plt.tight_layout()
    if histosData!=None:
        if Pass:
            outFile='CRPass_'+year+'_mjj'
        else:
            outFile='CRFail_'+year+'_mjj'
    else:
        if Pass:
            outFile='SRPass_'+year+'_mjj'
        else:
            outFile='SRFail_'+year+'_mjj'
    print("Saving {0}".format(outFile))

    plt.savefig(outFile)
    plt.savefig(outFile.replace(".png",".pdf"))

    plt.clf()

    return

def make_histograms_CR(mcpaths,datapath,year):
    
    histosDataPass = []
    edgesDataPass  = []
    colorsDataPass = []
    labelsDataPass = []
    histosSigPass  = []
    edgesSigPass   = []
    labelsSigPass  = []
    colorsSigPass  = []

    histosDataFail = []
    edgesDataFail  = []
    colorsDataFail = []
    labelsDataFail = []
    histosSigFail  = []
    edgesSigFail   = []
    labelsSigFail  = []
    colorsSigFail  = []

    print(year)
    for item in mcpaths:

        #print(mcpaths)
        mcpath=item[0]
        process=item[1]
        mcsig=r.TFile.Open(mcpath)
        mckeys=mcsig.GetListOfKeys()
        for key in mckeys:
            keyname=key.GetName()
            if '_nom' not in keyname:
                continue
            if 'SR' in keyname:
                continue
            if 'Pass' in keyname:
                pf_flag=True
                #index=1
            elif 'Fail' in keyname:
                pf_flag=False
                #index=0
            else:
                print('pf_flag error')

            if process=='TTToHadronic':
                color='cornflowerblue'
            elif process=='TTToSemiLeptonic':
                color='darkblue'
            elif 'QCD' in process:
                color='burlywood'
                

            print('Creating histogram for mc ',process,keyname)
            h=mcsig.Get(keyname)
            projection = h.ProjectionX("proj_name")
            hist, edges = hist2array(projection,return_edges=True)
            if pf_flag:
                histosSigPass.append(hist)
                labelsSigPass.append(process)
                colorsSigPass.append(color)
                edgesSigPass.append(edges[0])
            else:
                histosSigFail.append(hist)
                labelsSigFail.append(process)
                colorsSigFail.append(color)
                edgesSigFail.append(edges[0])


    #print(datapath)
    data=r.TFile.Open(datapath[0][0])
    process=datapath[0][1]
    datakeys=data.GetListOfKeys()

    for key in datakeys:

        color='k'

        keyname=key.GetName()
        if 'SR' in keyname:
            print('No histogram created')
            continue
        if 'Pass' in keyname:
            pf_flag=True
        elif 'Fail' in keyname:
            pf_flag=False
        else:
            print('pf_flag error')
        print('Creating histogram for data ',keyname)
        h=data.Get(keyname)
        projection = h.ProjectionX("proj_name")
        hist, edges = hist2array(projection,return_edges=True)
        if pf_flag:
            histosDataPass.append(hist)
            labelsDataPass.append(process)
            colorsDataPass.append(color)
            edgesDataPass.append(edges[0])
        else:
            histosDataFail.append(hist)
            labelsDataFail.append(process)
            colorsDataFail.append(color)
            edgesDataFail.append(edges[0])

    plot(histosDataPass,edgesDataPass,colorsDataPass,labelsDataPass,histosSigPass,edgesSigPass,colorsSigPass,labelsSigPass,True)
    plot(histosDataFail,edgesDataFail,colorsDataFail,labelsDataFail,histosSigFail,edgesSigFail,colorsSigFail,labelsSigFail,False)

    return

def make_histograms_SR(mcpaths,datapath,year):
    
    histosSigPass  = []
    edgesSigPass   = []
    labelsSigPass  = []
    colorsSigPass  = []

    histosSigFail  = []
    edgesSigFail   = []
    labelsSigFail  = []
    colorsSigFail  = []

    print(year)
    for item in mcpaths:

        #print(mcpaths)
        mcpath=item[0]
        process=item[1]
        mcsig=r.TFile.Open(mcpath)
        mckeys=mcsig.GetListOfKeys()
        for key in mckeys:
            keyname=key.GetName()
            if '_nom' not in keyname:
                continue
            if 'CR' in keyname:
                continue
            if 'Pass' in keyname:
                pf_flag=True
                #index=1
            elif 'Fail' in keyname:
                pf_flag=False
                #index=0
            else:
                print('pf_flag error')

            if process=='TTToHadronic':
                color='cornflowerblue'
            elif process=='TTToSemiLeptonic':
                color='darkblue'
            elif 'QCD' in process:
                color='burlywood'
                

            print('Creating histogram for mc ',process,keyname)
            h=mcsig.Get(keyname)
            projection = h.ProjectionX("proj_name")
            hist, edges = hist2array(projection,return_edges=True)
            if pf_flag:
                histosSigPass.append(hist)
                labelsSigPass.append(process)
                colorsSigPass.append(color)
                edgesSigPass.append(edges[0])
            else:
                histosSigFail.append(hist)
                labelsSigFail.append(process)
                colorsSigFail.append(color)
                edgesSigFail.append(edges[0])


    

    plot(None,None,None,None,histosSigPass,edgesSigPass,colorsSigPass,labelsSigPass,True)
    plot(None,None,None,None,histosSigFail,edgesSigFail,colorsSigFail,labelsSigFail,False)

    return


def plotrun2(path):
    data=r.TFile.Open(path)
    process='run2'
    datakeys=data.GetListOfKeys()

    histosData = []
    edgesData  = []
    colorsData = []
    labelsData = []
    errors=[]
    errorsCenters=[]

    for key in datakeys:

        keyname=key.GetName()
        if ('CR'in keyname) and ('Pass' in keyname):
            label='CRPass'
            color='k'
        elif ('CR'in keyname) and ('Fail' in keyname):
            label='CRFail'
            color='r'
        elif ('SR'in keyname) and ('Pass' in keyname):
            label='SRPass'
            color='g'
        elif ('SR'in keyname) and ('Fail' in keyname):
            label='SRFail'
            color='b'
        else:
            print('pf_flag error')
        print('Creating histogram for run 2 ',label)
        h=data.Get(keyname)
        projection = h.ProjectionX("proj_name")
        hist, edges = hist2array(projection,return_edges=True)
        histosData.append(hist)
        labelsData.append(label)
        colorsData.append(color)
        edgesData.append(edges[0])
        newerrors=[]
        for i in range(len(hist)):
            newerrors.append(np.sqrt(hist[i]))
        errors.append(newerrors)
        newCenters=[]
        for i in range(len(histosData[0])):
            newCenters.append((edgesData[0][i]+edgesData[0][i+1])/2)
        errorsCenters.append(newCenters)
    
    
    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()
    
    for i in range(len(histosData)):
        ax.errorbar(errorsCenters[i],histosData[i],yerr=errors[i],fmt='o', color=colorsData[i])
    hep.histplot(histosData,edgesData[0],stack=False,ax=ax,label=labelsData,linewidth=3,histtype="errorbar",color=colorsData)

    if(log):
        ax.set_yscale("log")

    yTitle='Events'
    xTitle='$M_{jj}$ [GeV]'

    ax.set_ylabel(yTitle)
    ax.set_xlabel(xTitle)
    plt.xlabel(xTitle, horizontalalignment='right', x=1.0)
    plt.ylabel(yTitle,horizontalalignment='right', y=1.0)
    
    regionText = "Run 2 Data"
    plt.ylim(10., 10**6)
    plt.xlim(1000, 3000)

    plt.text(0.8, 0.8, regionText,horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
    hep.cms.text("Simulation WiP",loc=0)
    #lumiText =  $XYZ fb^{-1} (13 TeV)$"    
    lumiText = f"Run2 (13 TeV)"    
    hep.cms.lumitext(lumiText)
    plt.legend(loc='upper right',ncol=2)#loc="best",ncol=2)#loc = 'best'    
    plt.tight_layout()
    
    outFile='run2_mjj'
    
    print("Saving {0}".format(outFile))

    plt.savefig(outFile)
    plt.savefig(outFile.replace(".png",".pdf"))

    plt.clf()

    return


def getfilepaths(process,year,data):
    if data:
        filepaths=[]
        newfile=base_path+'templates_data_obs_'+year+'.root'
        if not os.path.exists(newfile):
            print('File does not exist: ',newfile)
        else:
            filepaths.append([newfile,year+' data'])
        return filepaths
    else:
        if process=='MXMY':
            filepaths=[]
            MXs=['1200','1400','1600','1800','2000','2200','2400','2500','2600','2800','3000','3500','4000']
            for MX in MXs:
                newfile=base_path+'templates_MX'+MX+'_MY90_'+year+'.root'
                if not os.path.exists(newfile):
                    print('File does not exist: ',newfile)
                else:
                    filepaths.append([newfile,'MX '+MX+' MY 90'])
            return filepaths
        elif process=='QCD':
            filepaths=[]
            QCD_ranges=['700to1000','1000to1500','1500to2000','2000toInf']
            for range in QCD_ranges:
                newfile=base_path+'templates_QCD_HT'+range+'_'+year+'.root'    
                if not os.path.exists(newfile):
                    print('File does not exist: ',newfile)
                else:
                    filepaths.append([newfile,'QCD '+range])
            return filepaths
        elif process=='TTToHadronic':
            filepaths=[]
            newfile=base_path+'templates_TTToHadronic_'+year+'.root'
            if not os.path.exists(newfile):
                print('File does not exist: ',newfile)
            else:
                filepaths.append([newfile,process])
            return filepaths
        elif process=='TTToSemiLeptonic':
            filepaths=[]
            newfile=base_path+'templates_TTToSemiLeptonic_'+year+'.root'
            if not os.path.exists(newfile):
                print('File does not exist: ',newfile)
            else:
                filepaths.append([newfile,process])
            return filepaths
        elif process=='run2':
            return base_path+'templates_data_obs_run2.root'
        else:
            print("Invalid process name")
            return
    


if __name__ == '__main__':
    wp = "tight_medium"

    years=['2016','2016APV','2017','2018']
    processes=['QCD','TTToHadronic','TTToSemiLeptonic']#['MXMY','QCD','TTToHadronic','TTToSemiLeptonic']
    #years=['2016']

    for year in years:
        filepaths=[]
        for process in processes:
            mcpaths=getfilepaths(process,year,False)
            for item in mcpaths:
                filepaths.append(item)
        data=getfilepaths('x',year,True)
        make_histograms_CR(filepaths,data,year)
        make_histograms_SR(filepaths,data,year)

    newprocess='run2'
    path=getfilepaths(newprocess,None,False)
    plotrun2(path)