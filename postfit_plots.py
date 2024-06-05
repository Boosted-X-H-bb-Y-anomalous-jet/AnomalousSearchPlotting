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

matplotlib.use('Agg')
r.gROOT.SetBatch(True)
r.gStyle.SetOptFit(111)

def plot_limits(path_template,MX_values,output_name,prefit_xsec_in_fb=1.,observed=False):
    limits = {
        "expected": [],
        "observed": [],
        "minus_2sigma": [],
        "minus_1sigma": [],
        "plus_1sigma": [],
        "plus_2sigma": []
    }

    for MX in MX_values:
        path = path_template.format(MX)
        data = root2array(path, treename='limit')

        limits["expected"].append(data['limit'][2]*prefit_xsec_in_fb)
        limits["minus_2sigma"].append(data['limit'][0]*prefit_xsec_in_fb)
        limits["minus_1sigma"].append(data['limit'][1]*prefit_xsec_in_fb)
        limits["plus_1sigma"].append(data['limit'][3]*prefit_xsec_in_fb)
        limits["plus_2sigma"].append(data['limit'][4]*prefit_xsec_in_fb)
        if observed:
            limits["observed"].append(data['limit'][5]*prefit_xsec_in_fb) #Don't have observed yet

    for key in limits:
        limits[key] = np.array(limits[key])

    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()

    hep.cms.text("WiP",loc=0)
    lumiText = "138 $fb^{-1} (13 TeV)$"    
    hep.cms.lumitext(lumiText)
    ax.set_xlim(MX_values[0],MX_values[-1])

    plt.fill_between(MX_values, limits["minus_2sigma"], limits["plus_2sigma"], color='darkorange', label='Expected ±2$\sigma$')
    plt.fill_between(MX_values, limits["minus_1sigma"], limits["plus_1sigma"], color='forestgreen', label='Expected ±1$\sigma$')
    plt.plot(MX_values, limits["expected"], color='black', linestyle='--', label='Expected')
    if observed:
        plt.plot(MX_values, limits["observed"], color='red', marker='o', linestyle='-', label='Observed')

    plt.xlabel("$M_{X} [GeV]$")
    plt.ylabel(r'$\sigma(pp \rightarrow X \rightarrow H(bb)Y(WW(4q))})\,[fb]$',horizontalalignment='right', y=1.0)
    plt.xlabel("$M_{X} [GeV]$", horizontalalignment='right', x=1.0)
    plt.text(1400, 40, r'$M_{Y} = 90\,GeV$')
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

MX_values = [1200, 1400, 1600, 2000, 2500, 3000]
path_template = '/uscms_data/d3/roguljic/el8_anomalous/fitting/AnomalousSearchFits/SR_run2/MX{}_MY90-1_area/higgsCombineTest.AsymptoticLimits.mH120.root'
xsec = 5.0#fb
plot_limits(path_template,MX_values,"bbWW_limits_MY90",observed=False,prefit_xsec_in_fb=xsec)