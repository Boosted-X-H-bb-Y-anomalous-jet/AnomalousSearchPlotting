import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import ROOT as r
plt.style.use(hep.style.CMS)
from PyHist import PyHist
import copy
import matplotlib

def get_binning_x(hLow,hSig,hHigh):
    bins = []
    for i in range(1,hLow.GetNbinsX()+1):
        bins.append(hLow.GetXaxis().GetBinLowEdge(i))
    for i in range(1,hSig.GetNbinsX()+1):
        bins.append(hSig.GetXaxis().GetBinLowEdge(i))
    for i in range(1,hHigh.GetNbinsX()+2):#low edge of overflow is high edge of last bin
        bins.append(hHigh.GetXaxis().GetBinLowEdge(i))
    bins = np.array(bins,dtype='float64')
    return bins

def get_binning_y(hLow,hSig,hHigh):
    #histos should have same binning in Y
    bins = []
    for i in range(1,hLow.GetNbinsY()+2):
        bins.append(hLow.GetYaxis().GetBinLowEdge(i))
    bins = np.array(bins,dtype='float64')
    return bins

def get2DPostfitPlot(file, process, region, prefit=False):
    if not r.TFile.Open(file):
        raise FileNotFoundError(f"The file '{file}' does not exist or cannot be opened.")
    
    f = r.TFile.Open(file)
    fitStatus = "prefit" if prefit else "postfit"
    
    hLow = f.Get(f"{region}_LOW_{fitStatus}/{process}")
    if not hLow:
        raise ValueError(f"Histogram '{region}_LOW_{fitStatus}/{process}' does not exist in the file '{file}'.")
    
    hSig = f.Get(f"{region}_SIG_{fitStatus}/{process}")
    if not hSig:
        raise ValueError(f"Histogram '{region}_SIG_{fitStatus}/{process}' does not exist in the file '{file}'.")
    
    hHigh = f.Get(f"{region}_HIGH_{fitStatus}/{process}")
    if not hHigh:
        raise ValueError(f"Histogram '{region}_HIGH_{fitStatus}/{process}' does not exist in the file '{file}'.")
    
    h2 = merge_low_sig_high(hLow, hSig, hHigh, hName=f"h2_{process}_{region}")
    h2.SetDirectory(0)
    return h2

def merge_low_sig_high(hLow,hSig,hHigh,hName="temp"):
    n_x_low     = hLow.GetNbinsX()
    n_x_sig     = hSig.GetNbinsX()
    n_x_high    = hHigh.GetNbinsX()
    n_x         = n_x_low + n_x_sig + n_x_high
    n_y         = hLow.GetNbinsY()#assumes Y bins are the same
    bins_x      = get_binning_x(hLow,hSig,hHigh)
    bins_y      = get_binning_y(hLow,hSig,hHigh)
    h_res       = r.TH2F(hName,"",n_x,bins_x,n_y,bins_y)
    for i in range(1,n_x_low+1):
        for j in range(1,n_y+1):
            h_res.SetBinContent(i+0,j,hLow.GetBinContent(i,j))
            h_res.SetBinError(i+0,j,hLow.GetBinError(i,j))

    for i in range(1,n_x_sig+1):
        for j in range(1,n_y+1):
            h_res.SetBinContent(i+n_x_low,j,hSig.GetBinContent(i,j))
            h_res.SetBinError(i+n_x_low,j,hSig.GetBinError(i,j))

    for i in range(1,n_x_high+1):
        for j in range(1,n_y+1):
            h_res.SetBinContent(i+n_x_sig+n_x_low,j,hHigh.GetBinContent(i,j))
            h_res.SetBinError(i+n_x_sig+n_x_low,j,hHigh.GetBinError(i,j))
    return h_res

def get_hists(input_file,region,processes):
    histos_region_dict = {}
    for process in processes:
        h2 = get2DPostfitPlot(input_file,process,region)
        histos_region_dict[process]=h2
    return histos_region_dict

def plotShapesWithRatioAndBand(hData,hMC,hTotalBkg,labelsMC,colorsMC,xlabel,outputFile,xRange=[],yRange=[],projectionText="",yRangeLog=[]):
    #Adapting this to new class
    edges = hData.bin_edges  
    centresData = hData.get_bin_centers()

    plt.style.use([hep.style.CMS])
    matplotlib.rcParams.update({'font.size': 35})
    f, axs = plt.subplots(2,1, sharex=True, sharey=False,gridspec_kw={'height_ratios': [4, 1],'hspace': 0.00}, figsize=(10, 11),constrained_layout=True)
    axs = axs.flatten()
    plt.sca(axs[0])

    mc_yields = [histo.bin_values for histo in hMC]
    hep.histplot(mc_yields,edges,stack=True,label = labelsMC, histtype="fill",facecolor=colorsMC)

    yerr = hData.get_error_pairs()
    xerrorsData = [width / 2 for width in hData.bin_widths]
    plt.errorbar(centresData,hData.bin_values,xerr=xerrorsData,yerr=yerr, fmt='o',color="k",label = "Observed",zorder=11,markersize=10)

    def calcRatio(hData,hMC,dataErrs):
        ratioVals=[]
        ratioErrs = copy.deepcopy(dataErrs)
        for i in range(len(hData)):
            print(hData.get_error_pairs()[0][0])
            data      = hData[i]
            mc        = hMC[i]
            ratioVal     = data/(mc+0.000001)#Protect division by zero
            ratioErrs[0][i] = ratioErrs[0][i]/(mc+0.000001)
            ratioErrs[1][i] = ratioErrs[1][i]/(mc+0.000001)
            ratioVals.append(ratioVal)

        return ratioVals, ratioErrs

    def calcSystBand(hMC,uncBand):
        systBand = [[],[]]
        for i in range(len(hMC)):
            systBandUp = uncBand[0][i]/hMC[i]
            systBandDn = uncBand[1][i]/hMC[i]
            systBand[0].append(1+systBandUp)
            systBand[1].append(1-systBandDn)
        return systBand

    def calcPulls(hData,hMC,dataErrs,uncBand):
        pulls = []
        for i in range(len(hMC)):
            num = hData[i]-hMC[i]
            if(num>0):
                err_data = dataErrs[0][i]#error low
                err_mc = uncBand[1][i]#error high
            else:
                err_data = dataErrs[1][i]
                err_mc = uncBand[0][i]
            den = np.sqrt(err_data**2+err_mc**2)
            if(den==0):
                den = 0.0000000001#Avoid division by zero
            pull = num/den
            pulls.append(pull)
        return pulls

    axs[0].legend()
    plt.ylabel("Events / GeV",horizontalalignment='center', y=0.5)
    if(xRange):
        axs[0].set_xlim(xRange)
    if(yRange):
        axs[0].set_ylim(yRange)
    else:
        yMaximum = max(hData.bin_values)*1.4
        axs[0].set_ylim([0.,yMaximum])
    
    
    lumi = 138
    lumiText = str(lumi)+ " $fb^{-1} (13 TeV)$"    
    hep.cms.lumitext(lumiText)
    hep.cms.text("WiP",loc=0)
    plt.legend(loc="best",ncol=1)

    if(projectionText):
        plt.text(0.80, 0.60, projectionText, horizontalalignment='center',verticalalignment='center',transform=axs[0].transAxes)
    
    plt.sca(axs[1])#switch to lower pad
    #If plotting ratios
    if False:
        axs[1].set_ylabel("Data/bkg.",horizontalalignment='center', y=0.5)
        ratioVals, ratioErrs = calcRatio(hData.bin_values,hTotalBkg.bin_values,hData.get_error_pairs())
        systBand = calcSystBand(hTotalBkg.bin_values,hTotalBkg.get_error_pairs())
        axs[1].axhline(y=1.0, xmin=0, xmax=1, color="grey",linestyle="--",alpha=0.5)
        axs[1].set_ylim([0.0,2.3])
        plt.xlabel(xlabel,horizontalalignment='center', x=0.5)
        plt.errorbar(centresData,ratioVals,yerr=ratioErrs,xerr=xerrorsData, fmt='o',color="k",markersize=10) 

        axs[1].tick_params(axis='x', pad=10)  #Fix overlap between numbers on x and y axis of ratio plot
        systBand_lower = list(systBand[0]) + [systBand[0][-1]]
        systBand_upper = list(systBand[1]) + [systBand[1][-1]]
        plt.fill_between(edges, systBand_lower, systBand_upper, color='lightgrey', alpha=0.8, label="Bkg. uncert.",step='pre')
        plt.legend(bbox_to_anchor=(-0.01, 1.1),loc='upper left')
    else:#Plotting pulls  
        pulls = calcPulls(hData.bin_values,hTotalBkg.bin_values,hData.get_error_pairs(),hTotalBkg.get_error_pairs())  
        axs[1].set_ylabel("Pulls",horizontalalignment='center', y=0.5)
        axs[1].axhline(y=0.0, xmin=0, xmax=1, color="grey",linestyle="--",alpha=0.5)
        axs[1].set_ylim([-2.3,2.3])
        plt.xlabel(xlabel,horizontalalignment='center', x=0.5)
        hep.histplot(pulls,edges,linewidth=1,histtype="fill",facecolor="grey",edgecolor='black')

    print("Saving ", outputFile)
    #plt.tight_layout()
    plt.savefig(outputFile,bbox_inches="tight")
    plt.savefig(outputFile.replace("png","pdf"),bbox_inches="tight")

    axs[0].set_yscale("log")
    if yRangeLog:
        axs[0].set_ylim(yRangeLog)
    
    log_outputFile = outputFile.replace(".png", "_log.png")
    log_outputFile_pdf = log_outputFile.replace(".png", ".pdf")
    print(f"Saving {log_outputFile}")
    plt.savefig(log_outputFile,bbox_inches="tight")
    plt.savefig(log_outputFile_pdf,bbox_inches="tight")

    plt.clf()
    plt.cla()

def plot_projection(histos_dict, region, processes, labels_dict, colors_dict, axis="X",yRange=[],yRangeLog=[1,10**5]):
    if axis not in ["X", "Y"]:
        raise ValueError("Invalid axis. Choose 'X' or 'Y'.")

    axis_label = "$M_{JJ} [GeV]$" if axis == "X" else "$M_{J}^{Y} [GeV]$"
    file_suffix = "mjj" if axis == "X" else "mjy"
    x_range = [1300, 3500] if axis == "X" else [40,500]
    projection_method = f"Projection{axis}"

    # Extract histograms
    histos = histos_dict[region]

    h_data = getattr(histos["data_obs"], projection_method)(f"data_{axis}")
    h_data.SetBinErrorOption(1)
    h_data = PyHist(h_data)
    h_data.divide_by_bin_width()

    h_mc = []
    labels_mc = []
    colors_mc = []
    for process in processes:
        h_temp = getattr(histos[process], projection_method)(f"{process}_{axis}")
        h_temp = PyHist(h_temp)
        h_temp.divide_by_bin_width()
        h_mc.append(h_temp)
        labels_mc.append(labels_dict[process])
        colors_mc.append(colors_dict[process])

    h_bkg = getattr(histos["TotalBkg"], projection_method)(f"bkg_{axis}")
    h_bkg = PyHist(h_bkg)
    h_bkg.divide_by_bin_width()

    plotShapesWithRatioAndBand(h_data, h_mc, h_bkg,labels_mc, colors_mc,axis_label,f"{region}_{file_suffix}.png",xRange=x_range,yRange=yRange,projectionText=region.replace("_", " "),yRangeLog=yRangeLog)

if __name__ == "__main__":
    if False:
        input_file = "~/nobackup/el8_anomalous/el9_fitting/CMSSW_14_1_0_pre4/src/AnomalousSearchFits/CR_run2/MX1400_MY90-0_area/postfitshapes_b.root"
        processes_CR_Pass =["WJets800","ZJets800","TTToHadronic","TTToSemiLeptonic","Background_0","TotalBkg","data_obs"]
        processes_CR_Fail =["WJets800","ZJets800","TTToHadronic","TTToSemiLeptonic","Background","TotalBkg","data_obs"]
        histos_dict={}
        histos_dict["CR_Pass"] = get_hists(input_file,"CR_Pass",processes_CR_Pass)
        histos_dict["CR_Fail"] = get_hists(input_file,"CR_Fail",processes_CR_Fail)
        labels_dict={"WJets800":"__nolabel__","ZJets800":"V+jets","TTToHadronic":r"$t\bar t$","TTToSemiLeptonic":"__nolabel__","Background_0":"Multijet","Background":"Multijet"}
        colors_dict={"WJets800":"#1f77b4","ZJets800":"#1f77b4","TTToHadronic":"#d42e12","TTToSemiLeptonic":"#d42e12","Background_0":"#f39c12","Background":"#f39c12"}
        
        plot_projection(histos_dict,"CR_Pass",["Background_0","TTToHadronic","TTToSemiLeptonic","WJets800","ZJets800"],labels_dict,colors_dict,axis="X",yRange=[0, 8],yRangeLog=[0.01,10**2])    
        plot_projection(histos_dict,"CR_Pass",["Background_0","TTToHadronic","TTToSemiLeptonic","WJets800","ZJets800"],labels_dict,colors_dict,axis="Y",yRange=[0, 20],yRangeLog=[0.01,10**2])
        plot_projection(histos_dict,"CR_Fail",["Background","TTToHadronic","TTToSemiLeptonic","WJets800","ZJets800"],labels_dict,colors_dict,axis="X",yRange=[],yRangeLog=[1,10**5])
        plot_projection(histos_dict,"CR_Fail",["Background","TTToHadronic","TTToSemiLeptonic","WJets800","ZJets800"],labels_dict,colors_dict,axis="Y",yRange=[],yRangeLog=[1,10**5])
    elif True:
        input_file = "~/nobackup/el8_anomalous/el9_fitting/CMSSW_14_1_0_pre4/src/AnomalousSearchFits/SR_run2_one_signal/MX2200_MY250-2_area/postfitshapes_b.root"
        processes_SR_Pass =["TTToHadronic","TTToSemiLeptonic","Background_2","TotalBkg","data_obs"]
        processes_SR_Fail =["TTToHadronic","TTToSemiLeptonic","Background","TotalBkg","data_obs"]
        histos_dict={}
        histos_dict["SR_Pass"] = get_hists(input_file,"SR_Pass",processes_SR_Pass)
        histos_dict["SR_Fail"] = get_hists(input_file,"SR_Fail",processes_SR_Fail)
        labels_dict={"TTToHadronic":r"$t\bar t$","TTToSemiLeptonic":"__nolabel__","Background_2":"Multijet","Background":"Multijet"}
        colors_dict={"TTToHadronic":"#d42e12","TTToSemiLeptonic":"#d42e12","Background_2":"#f39c12","Background":"#f39c12"}

        plot_projection(histos_dict,"SR_Pass",["Background_2","TTToHadronic","TTToSemiLeptonic"],labels_dict,colors_dict,axis="X",yRange=[0, 13],yRangeLog=[0.01,10**2])    
        plot_projection(histos_dict,"SR_Pass",["Background_2","TTToHadronic","TTToSemiLeptonic"],labels_dict,colors_dict,axis="Y",yRange=[0, 30],yRangeLog=[0.01,10**2])
        plot_projection(histos_dict,"SR_Fail",["Background","TTToHadronic","TTToSemiLeptonic"],labels_dict,colors_dict,axis="X",yRange=[],yRangeLog=[1,10**5])
        plot_projection(histos_dict,"SR_Fail",["Background","TTToHadronic","TTToSemiLeptonic"],labels_dict,colors_dict,axis="Y",yRange=[],yRangeLog=[1,10**5])

    else:
        exit()