import matplotlib
matplotlib.use('Agg')

import ROOT as r
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mplhep as hep
from root_numpy import hist2array, root2array
import ctypes
from pathlib import Path
r.gErrorIgnoreLevel = r.kFatal

matplotlib.use('Agg')
r.gROOT.SetBatch(True)
r.gStyle.SetOptFit(111)

import csv

def write_limits_to_csv(filename, MX_values, limits):
    fieldnames = ["MX", "expected", "minus2sigma", "minus1sigma", "plus1sigma", "plus2sigma"]
    if "observed" in limits and len(limits["observed"]) > 0:
        fieldnames.append("observed")

    with open(filename, mode='w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, MX in enumerate(MX_values):
            row = {
                "MX": MX,
                "expected": f"{limits['expected'][i]:.3f}",
                "minus2sigma": f"{limits['minus2sigma'][i]:.3f}",
                "minus1sigma": f"{limits['minus1sigma'][i]:.3f}",
                "plus1sigma": f"{limits['plus1sigma'][i]:.3f}",
                "plus2sigma": f"{limits['plus2sigma'][i]:.3f}",
            }
            if "observed" in limits and len(limits["observed"]) > 0:
                row["observed"] = f"{limits['observed'][i]:.3f}"
            writer.writerow(row)

def plot_1d_limits_WW(path_template,MX_values,my,output_name,prefit_xsec_in_fb=1.,observed=False):
    limits = {
        "expected": [],
        "observed": [],
        "minus2sigma": [],
        "minus1sigma": [],
        "plus1sigma": [],
        "plus2sigma": []
    }

    for MX in MX_values:
        path = path_template.format(MX,my)
        data = root2array(path, treename='limit')

        limits["expected"].append(data['limit'][2]*prefit_xsec_in_fb)
        limits["minus2sigma"].append(data['limit'][0]*prefit_xsec_in_fb)
        limits["minus1sigma"].append(data['limit'][1]*prefit_xsec_in_fb)
        limits["plus1sigma"].append(data['limit'][3]*prefit_xsec_in_fb)
        limits["plus2sigma"].append(data['limit'][4]*prefit_xsec_in_fb)
        if observed:
            limits["observed"].append(data['limit'][5]*prefit_xsec_in_fb)

    for key in limits:
        limits[key] = np.array(limits[key])

    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()

    hep.cms.text("WiP",loc=0)
    lumiText = "138 $fb^{-1}$ (13 TeV)"    
    hep.cms.lumitext(lumiText)
    ax.set_xlim(MX_values[0],MX_values[-1])

    plt.fill_between(MX_values, limits["minus2sigma"], limits["plus2sigma"], color='darkorange', label='Expected ±2$\sigma$')
    plt.fill_between(MX_values, limits["minus1sigma"], limits["plus1sigma"], color='forestgreen', label='Expected ±1$\sigma$')
    plt.plot(MX_values, limits["expected"], color='black', linestyle='--', label='Expected')
    if observed:
        plt.plot(MX_values, limits["observed"], color='black', marker='o', linestyle='-', label='Observed')

    plt.xlabel("$M_{X} [GeV]$")
    plt.ylabel(r'$\sigma(pp \rightarrow X \rightarrow H(bb)Y(WW(4q))})\,[fb]$',horizontalalignment='right', y=1.0)
    plt.xlabel("$M_{X} [GeV]$", horizontalalignment='right', x=1.0)
    plt.text(1700, 40, '$M_{Y}$ = '+f'{my} GeV')
    plt.yscale('log')
    plt.ylim(0.1, 10**2)
    plt.legend(loc=(0.10,0.60)
      , title='95% CL upper limits'
      ,ncol=2
      ,title_fontsize=17
      ,fontsize=17
      )

    plt.legend()

    print(f"Saving plots/{output_name}.png")
    plt.savefig(f"plots/{output_name}.pdf")
    plt.savefig(f"plots/{output_name}.png")
def plot_1d_limits_TPrime(path_template, MX_values, output_name, prefit_xsec_in_fb=1., observed=False):
    limits = {
        "expected": [],
        "observed": [],
        "minus2sigma": [],
        "minus1sigma": [],
        "plus1sigma": [],
        "plus2sigma": []
    }

    for MX in MX_values:
        path = path_template.format(MX)
        data = root2array(path, treename='limit')

        limits["expected"].append(data['limit'][2] * prefit_xsec_in_fb)
        limits["minus2sigma"].append(data['limit'][0] * prefit_xsec_in_fb)
        limits["minus1sigma"].append(data['limit'][1] * prefit_xsec_in_fb)
        limits["plus1sigma"].append(data['limit'][3] * prefit_xsec_in_fb)
        limits["plus2sigma"].append(data['limit'][4] * prefit_xsec_in_fb)
        if observed:
            limits["observed"].append(data['limit'][5] * prefit_xsec_in_fb)

    for key in limits:
        limits[key] = np.array(limits[key])

    # ---- Plotting the expected line with projections ----
    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()

    hep.cms.text("", loc=1)
    lumiText = "138 $fb^{-1}$ (13 TeV)"
    hep.cms.lumitext(lumiText)
    ax.set_xlim(MX_values[0], MX_values[-1])

    plt.plot(MX_values, limits["expected"], color='black', linestyle='--', label='Expected', linewidth=2)
    plt.xlabel("$M_{X} [GeV]$")
    plt.ylabel(r"$\sigma(pp \rightarrow \mathrm{X} \rightarrow \mathrm{H(bb)\,Y(bqq)})})\,[fb]$", horizontalalignment='right', y=1.0)
    plt.yscale('log')
    plt.ylim(0.2, 100)

    # # Load predicted limits from pickle file
    # import pickle
    # with open("delphes_pred_limits.pkl", "rb") as f:
    #     pred_limits = pickle.load(f)
    #     print(pred_limits)
    # delphes_predicted = []
    # for MX in MX_values:
    #     key = f"TPrime_MX-{MX}_MY-125"
    #     if key in pred_limits:
    #         delphes_predicted.append(pred_limits[key]["delphesPred"])

    # plt.plot(MX_values, delphes_predicted, color='blue', marker='s', linestyle='-', linewidth=2, label='Delphes-predicted')

    plt.legend(loc=(0.01, 0.65), title='95% CL upper limits', ncol=2, title_fontsize=26, fontsize=26)

    print(f"Saving plots/{output_name}.png")
    plt.savefig(f"plots/{output_name}.pdf")
    plt.savefig(f"plots/{output_name}.png")
    
    # ---- Create the Brazilian flag plot (shaded regions) ----
    f_brazilian, ax_brazilian = plt.subplots()

    hep.cms.text("", loc=1)
    hep.cms.lumitext(lumiText)
    ax_brazilian.set_xlim(MX_values[0], MX_values[-1]+6)
    ax_brazilian.set_xticks(np.arange(1400, 3200, 200))
    #ax_brazilian.set_xticks(np.arange(1400, 3000, 200), minor=True)
    #ax_brazilian.grid(True, which='both', axis='x', linestyle='--', linewidth=0.5)

    plt.fill_between(MX_values, limits["minus2sigma"], limits["plus2sigma"], color='#F5BB54', label=r'95% CL exp.')
    plt.fill_between(MX_values, limits["minus1sigma"], limits["plus1sigma"], color="#607641", label=r'68% CL exp.')
    if observed:
        plt.plot(MX_values, limits["observed"], color='black', marker='o', linestyle='-', label='Observed')

    # Plot the expected line
    plt.plot(MX_values, limits["expected"], color='black', linestyle='--', label='Expected', linewidth=2)

    plt.xlabel("$M_{X} [GeV]$")
    plt.ylabel(r"$\sigma(pp \rightarrow \mathrm{X} \rightarrow \mathrm{H(bb)\,Y(bqq)})})$ [fb]", horizontalalignment='right', y=1.0)
    plt.yscale('log')
    plt.ylim(0.2, 100)

    # Add the legend
    handles, labels = plt.gca().get_legend_handles_labels()

    # Swap the legend order for 68% and 95% expected
    handles[2], handles[3] = handles[3], handles[2]
    labels[2], labels[3] = labels[3], labels[2]
    plt.legend(handles=handles, labels=labels,loc="upper right", title='95% CL upper limits', ncol=1, title_fontsize=26, fontsize=26)

    # Save the Brazilian flag plot separately
    print(f"Saving plots/brazilian_{output_name}_plot.pdf")
    plt.savefig(f"plots/brazilian_flag_{output_name}.pdf")
    plt.savefig(f"plots/brazilian_flag_{output_name}.png")

def plot_1d_limits_generic(path_template, MX_values, MY_value, output_name, key_prefix, ylabel, prefit_xsec_in_fb=1., observed=False,label=""):
    limits = {
        "expected": [],
        "observed": [],
        "minus2sigma": [],
        "minus1sigma": [],
        "plus1sigma": [],
        "plus2sigma": []
    }

    for MX in MX_values:
        path = path_template.format(MX, MY_value)
        data = root2array(path, treename='limit')

        limits["expected"].append(data['limit'][2] * prefit_xsec_in_fb)
        limits["minus2sigma"].append(data['limit'][0] * prefit_xsec_in_fb)
        limits["minus1sigma"].append(data['limit'][1] * prefit_xsec_in_fb)
        limits["plus1sigma"].append(data['limit'][3] * prefit_xsec_in_fb)
        limits["plus2sigma"].append(data['limit'][4] * prefit_xsec_in_fb)
        if observed:
            limits["observed"].append(data['limit'][5] * prefit_xsec_in_fb)

    for k in limits:
        limits[k] = np.array(limits[k])

    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()
    hep.cms.text("", loc=1)
    lumiText = "138 $fb^{-1}$ (13 TeV)"
    hep.cms.lumitext(lumiText)
    ax.set_xlim(MX_values[0], MX_values[-1])
    plt.fill_between(MX_values, limits["minus2sigma"], limits["plus2sigma"], color='#F5BB54', label=r'95% CL exp.')
    plt.fill_between(MX_values, limits["minus1sigma"], limits["plus1sigma"], color='#607641', label=r'68% CL exp.')
    plt.plot(MX_values, limits["expected"], color='black', linestyle='--', label='Median expected', linewidth=2)
    if observed:
        plt.plot(MX_values, limits["observed"], color='black', marker='o', linestyle='-', label='Observed')

    plt.xlabel("$M_{X} [GeV]$")
    plt.ylabel(ylabel, horizontalalignment='right', y=1.0)
    plt.yscale('log')
    plt.ylim(0.1, 100)

    # import pickle
    # with open("delphes_pred_limits.pkl", "rb") as f:
    #     pred_limits = pickle.load(f)
    # delphes_predicted = []
    # for MX in MX_values:
    #     dict_key = f"{key_prefix}_MX-{MX}_MY-{MY_value}"
    #     if dict_key in pred_limits:
    #         delphes_predicted.append(pred_limits[dict_key]["delphesPred"])

    # plt.plot(MX_values, delphes_predicted, color='blue', linestyle='-.', linewidth=2, label='Delphes-predicted median expected')

    # # Swap the legend order for 68% and Delphes-predicted expected
    # handles, labels = plt.gca().get_legend_handles_labels()
    # handles[2], handles[4] = handles[4], handles[2]
    # labels[2], labels[4] = labels[4], labels[2]
   
    #plt.legend(handles=handles, labels=labels,loc='upper right', bbox_to_anchor=(1.0, 1.0), title='95% CL upper limits', ncol=1, title_fontsize=26, fontsize=26)
    plt.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0), title='95% CL upper limits', ncol=1, title_fontsize=26, fontsize=26)
    
    plt.text(0.70, 0.55, label, horizontalalignment='center',verticalalignment='center',transform=ax.transAxes)

    print(f"Saving plots/{output_name}.png")
    plt.savefig(f"plots/{output_name}.pdf")
    plt.savefig(f"plots/{output_name}.png")

    csv_filename = f"plots/{output_name}.csv"
    write_limits_to_csv(csv_filename, MX_values, limits)
    print(f"Saved CSV: {csv_filename}")



def MakePlots():
    MX = [1400, 1600, 1800, 2000, 2200, 2600, 3000]
    MY = [90, 125, 190, 250, 300, 400]
    xsec = 5.0#fb
    bbww_path_template = '/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/CMSSW_14_1_0_pre4/src/AnomalousSearchFits/SR_run2/MX{}_MY{}-2_area/higgsCombineTest.AsymptoticLimits.mH120.root'
    for my in MY:
        continue
        plot_1d_limits_WW(bbww_path_template,MX,my,f"bbWW_limits_MY{my}",observed=True,prefit_xsec_in_fb=xsec)
    tprime_path_template = '/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/CMSSW_14_1_0_pre4/src/AnomalousSearchFits/SR_run2/TPrime_MX-{}_MY-125-2_area/higgsCombineTest.AsymptoticLimits.mH120.root'
    yuu_path_template = '/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/CMSSW_14_1_0_pre4/src/AnomalousSearchFits/SR_run2/XToYH_HTo2BYTo2Up_MX-{}_MY-200-2_area/higgsCombineTest.AsymptoticLimits.mH120.root'
    ytt_path_template = '/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/CMSSW_14_1_0_pre4/src/AnomalousSearchFits/SR_run2/XToYH_HTo2BYTo2T_Hadronic_MX-{}_MY-400-2_area/higgsCombineTest.AsymptoticLimits.mH120.root'
    yww_path_template = '/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/CMSSW_14_1_0_pre4/src/AnomalousSearchFits/SR_run2/MX{}_MY{}-2_area/higgsCombineTest.AsymptoticLimits.mH120.root'
    #plot_1d_limits_TPrime(tprime_path_template,[1400,1600,1800,2000,2400,3000],f"TPrime_limits",prefit_xsec_in_fb=xsec,observed=True)

    plot_1d_limits_generic(
            yww_path_template, MX, 250, 
            "Yww_limits", "XToYH_HTo2BYTo4q",r"$\sigma(pp \rightarrow \mathrm{X} \rightarrow \mathrm{H} \mathrm{Y}) \mathcal{B}(\mathrm{H}\rightarrow\mathrm{bb})\mathcal{B}(\mathrm{Y}\rightarrow\mathrm{WW}\rightarrow\mathrm{4q})$ [fb]",
            prefit_xsec_in_fb=xsec, observed=True, label = r"$\mathrm{Y}\rightarrow\mathrm{WW}\rightarrow\mathrm{4q}$"+"\n"+r"$M_Y = 250\ \mathrm{GeV}$"
        )
    
    plot_1d_limits_generic(
            yuu_path_template, MX, 200, 
            "Yuu_limits", "XToYH_HTo2BYTo2Up",r"$\sigma(pp \rightarrow \mathrm{X} \rightarrow \mathrm{H} \mathrm{Y}) \mathcal{B}(\mathrm{H}\rightarrow\mathrm{bb})\mathcal{B}(\mathrm{Y}\rightarrow\mathrm{uu})$ [fb]",
            prefit_xsec_in_fb=xsec, observed=True, label=r"$\mathrm{Y}\rightarrow\mathrm{uu}$" + "\n" + r"$M_Y = 200\ \mathrm{GeV}$"
        )

    plot_1d_limits_generic(
            ytt_path_template, MX, 400, 
            "Ytt_limits", "XToYH_HTo2BYTo2T_Hadronic", r"$\sigma(pp \rightarrow \mathrm{X} \rightarrow \mathrm{H} \mathrm{Y}) \mathcal{B}(\mathrm{H}\rightarrow\mathrm{bb})\mathcal{B}(\mathrm{Y}\rightarrow tt\rightarrow bqq bqq)$ [fb]",
            prefit_xsec_in_fb=xsec, observed=True, label=r"$\mathrm{Y}\rightarrow tt\rightarrow bqq\ bqq$" + "\n" + r"$M_Y = 400\ \mathrm{GeV}$"
        )

    plot_1d_limits_generic(
        tprime_path_template, [1400,1600,1800,2000,2400,3000], 125,
        "TPrime_limits", "TPrime",r"$\sigma(pp \rightarrow \mathrm{X} \rightarrow \mathrm{H} \mathrm{Y}) \mathcal{B}(\mathrm{H}\rightarrow\mathrm{bb})\mathcal{B}(\mathrm{Y}\rightarrow\mathrm{bqq})$ [fb]",
        prefit_xsec_in_fb=xsec, observed=True, label=r"$\mathrm{Y}\rightarrow\mathrm{bqq}$" + "\n" + r"$M_Y = 170\ \mathrm{GeV}$"
        )

if __name__=='__main__':
    MakePlots()