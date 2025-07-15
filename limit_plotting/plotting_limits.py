import ROOT as r
import os
import numpy as np
import math
import mplhep as hep
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import csv

def write_limits_csv(filename, MXs, MYs, base_path, fb_conversion):
    rows = []

    for MY in MYs:
        for MX in MXs:
            filepath = os.path.join(
                base_path,
                f"MX{MX}_MY{MY}-2_area/higgsCombineTest.AsymptoticLimits.mH120.root"
            )
            if not os.path.exists(filepath):
                print(f"Missing file: MX={MX}, MY={MY}")
                continue

            limits = get_limit_all(filepath, fb_conversion)
            row = {
                "MX": MX,
                "MY": MY,
            }
            # Format floats to .3f
            for key, value in limits.items():
                row[key] = f"{value:.3f}" if value is not None else ""
            rows.append(row)

    with open(filename, mode='w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["MX", "MY", "expected", "minus2sigma", "minus1sigma", "plus1sigma", "plus2sigma", "observed"])
        writer.writeheader()
        writer.writerows(rows)

# Suppress ROOT warnings
r.gROOT.SetBatch(True)
r.gErrorIgnoreLevel = r.kBreak

base_path = '/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/CMSSW_14_1_0_pre4/src/AnomalousSearchFits/SR_run2/'
second_base_path = '/eos/uscms/store/user/rkansal/bbVV/cards/25Feb19ResUnblinded/NMSSM_XToYHTo2W2BTo4Q2B_MX-{MX}_MY-{MY}/'
MXs = [1400, 1600, 1800, 2000, 2200, 2600, 3000]
MYs = [90, 125, 190, 250, 300, 400]

fb_conversion_first = 5
fb_conversion_second = 1

def get_limit_all(filepath, conversion_factor):
    file = r.TFile.Open(filepath)
    tree = file.Get("limit")

    results = {
        "minus2sigma": None,
        "minus1sigma": None,
        "expected": None,
        "plus1sigma": None,
        "plus2sigma": None,
        "observed": None,
    }

    quantile_indices = {
        0: "minus2sigma",
        1: "minus1sigma",
        2: "expected",
        3: "plus1sigma",
        4: "plus2sigma",
        5: "observed",
    }

    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        key = quantile_indices.get(i)
        if key is not None:
            results[key] = tree.limit * conversion_factor

    return results


def get_limit(filepath, conversion_factor, expected=True):
    file = r.TFile.Open(filepath)
    keyname = 'limit'

    tree = file.Get(keyname)
    if expected:
        exp_idx = 2  #Expected: 2, Observed: 5
    else:
        exp_idx = 5
    tree.GetEntry(exp_idx)
    limit = tree.limit * conversion_factor

    return limit

def Plot_v2(Limits, Ratio, expected=True):
    plt.style.use(hep.style.CMS)

    # Plot for the expected limits (first plot)
    plt.figure(figsize=(15, 10))
    
    # Adjust tick locations to multiples of 200
    plt.xticks(range(0, 3500, 200))
    plt.yticks(range(0, 500, 100))
    
    # Set axis limits with some padding
    plt.xlim(1300, 3100)
    plt.ylim(50, 450)

    # Flatten Limits for color normalization
    flat_limits = np.array([val for row in Limits for val in row if not np.isnan(val)])

    # Create color normalization
    norm = LogNorm(vmin=0.1, vmax=20.)
    cmap = plt.cm.viridis

    # Plot colored boxes for Limits
    for i in range(len(MYs)):
        for j in range(len(MXs)):
            limit_value = Limits[i][j]
            
            if not np.isnan(limit_value):
                # Map limit value to color
                color = cmap(norm(limit_value))
                
                # Create rectangle
                rect = plt.Rectangle((MXs[j]-50, MYs[i]-15), 100, 30, 
                                     facecolor=color, 
                                     edgecolor='black', 
                                     linewidth=1)
                plt.gca().add_patch(rect)
                
                # Add text annotation
                plt.text(MXs[j], MYs[i],
                         f'{limit_value:.1f}',
                         color='white',
                         ha='center',
                         va='center',
                         fontsize=17,
                         fontweight="bold")

    # Manually create ScalarMappable
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_clim(vmin=0.1, vmax=20.)
    
    # Add colorbar
    cbar = plt.colorbar(sm, pad=0.03)
    if expected:
        cbar.set_label('Expected 95% CL upper limits [fb]')
    else:
        cbar.set_label('Observed 95% CL upper limits [fb]')
        
    plt.ylabel(r'$M_{Y}$ [GeV]', fontsize=24)
    plt.xlabel(r'$M_{X}$ [GeV]', fontsize=24)

    hep.cms.text("", loc=0)
    lumi = 138
    lumiText = str(lumi)+ " $fb^{-1}$ (13 TeV)"    
    hep.cms.lumitext(lumiText)

    plt.tight_layout()
    if expected:
        plt.savefig('plots/Limit_logscale.png', bbox_inches='tight')
        plt.savefig('plots/Limit_logscale.pdf', bbox_inches='tight')
    else:
        plt.savefig('plots/Limit_logscale_obs.png', bbox_inches='tight')
        plt.savefig('plots/Limit_logscale_obs.pdf', bbox_inches='tight')

    # Plot for the ratio (second plot)
    plt.figure(figsize=(15, 10))
    
    # Adjust tick locations to multiples of 200
    plt.xticks(range(0, 3500, 200))
    plt.yticks(range(0, 500, 100))
    
    # Set axis limits with some padding
    plt.xlim(1200, 3300)
    plt.ylim(50, 450)

    # Flatten Ratio for color normalization
    flat_ratio = np.array([val for row in Ratio for val in row if not np.isnan(val)])

    # Create color normalization for Ratio
    norm_ratio = LogNorm(vmin=flat_ratio.min(), vmax=flat_ratio.max())
    cmap_ratio = plt.cm.coolwarm

    # Plot colored boxes for Ratio
    for i in range(len(MYs)):
        for j in range(len(MXs)):
            value = Ratio[i][j]
            
            if not np.isnan(value):
                # Map ratio value to color
                color = cmap_ratio(norm_ratio(value))
                
                # Create rectangle
                rect = plt.Rectangle((MXs[j]-50, MYs[i]-25), 100, 50, 
                                     facecolor=color, 
                                     edgecolor='black', 
                                     linewidth=1)
                plt.gca().add_patch(rect)
                
                # Add text annotation
                plt.text(MXs[j], MYs[i],
                         f'{value:.1f}',
                         color='black',
                         ha='center',
                         va='center',
                         fontsize=8,
                         fontweight="bold")

    # Manually create ScalarMappable for Ratio
    sm_ratio = plt.cm.ScalarMappable(cmap=cmap_ratio, norm=norm_ratio)
    sm_ratio.set_clim(flat_ratio.min(), flat_ratio.max())
    
    # Add colorbar
    cbar = plt.colorbar(sm_ratio, pad=0.03)
    cbar.set_label('Ratio of anomalous versus targeted search limits')
    
    plt.ylabel(r'$M_{Y}$ [GeV]', fontsize=24)
    plt.xlabel(r'$M_{X}$ [GeV]', fontsize=24)

    plt.tight_layout()
    if expected:
        plt.savefig('Limit_Ratio_logscale.png', bbox_inches='tight')
        plt.savefig('Limit_Ratio_logscale.pdf', bbox_inches='tight')
    else:
        plt.savefig('Limit_Ratio_logscale_obs.png', bbox_inches='tight')
        plt.savefig('Limit_Ratio_logscale_obs.pdf', bbox_inches='tight')

def Plot(Limits, Ratio, expected=True):
    plt.style.use(hep.style.CMS)

    x_edges = MXs + [3400]
    y_edges = MYs + [450]  # Changed from 425 to 450 for more space

    # Plot for the expected limits (first plot)
    plt.figure(figsize=(12, 8))
    mesh = plt.pcolormesh(x_edges, y_edges, Limits, shading='auto', cmap='viridis', 
                          norm=LogNorm(vmin=0.1, vmax=20))
    cbar = plt.colorbar(mesh, pad=0.03)
    if expected:
        cbar.set_label('Expected 95% CL upper limits [fb]')
    else:
        cbar.set_label('Observed 95% CL upper limits [fb]')
        
    plt.ylabel(r'$M_{Y}$ [GeV]', fontsize=24)
    plt.xlabel(r'$M_{X}$ [GeV]', fontsize=24)

    # Constant offset for text positioning
    #offset = 14  # Pixels in figure coordinates
    offset = 0  # Pixels in figure coordinates

    # Add text annotations for Limits with Ratio in parentheses below
    for i in range(len(MYs)):
        for j in range(len(MXs)):
            limit_value = Limits[i][j]
            ratio_value = Ratio[i][j]
            
            if not np.isnan(limit_value):
                # Calculate center of cell
                x_center = x_edges[j] + (x_edges[j+1] - x_edges[j]) / 2
                y_center = y_edges[i] + (y_edges[i+1] - y_edges[i]) / 2
                
                # Display limit value above center
                plt.text(x_center, y_center + offset * 0.5,
                         f'{limit_value:.1f}',
                         color='white',
                         ha='center',
                         va='center',
                         fontsize=17,
                         fontweight="bold",
                         bbox=dict(facecolor='none', alpha=0.6, edgecolor='none'))
                
                # # Display ratio value in parentheses below center if available
                # if not np.isnan(ratio_value):
                #     plt.text(x_center, y_center - offset * 0.5,
                #              f'({ratio_value:.1f})',
                #              color='white',
                #              ha='center',
                #              va='center',
                #              fontsize=14,
                #              bbox=dict(facecolor='none', alpha=0.6, edgecolor='none'))
    
    hep.cms.text("", loc=0)

    plt.tight_layout()
    if expected:
        plt.savefig('plots/Limit_logscale.png', bbox_inches='tight')
        plt.savefig('plots/Limit_logscale.pdf', bbox_inches='tight')
    else:
        plt.savefig('plots/Limit_logscale_obs.png', bbox_inches='tight')
        plt.savefig('plots/Limit_logscale_obs.pdf', bbox_inches='tight')

    # Plot for the ratio (second plot)
    plt.figure(figsize=(12, 8))
    mesh = plt.pcolormesh(x_edges, y_edges, Ratio, shading='auto', cmap='coolwarm', 
                          norm=LogNorm(vmin=0.1, vmax=10))
    cbar = plt.colorbar(mesh, pad=0.03)
    cbar.set_label('Ratio of anomalous versus targeted search limits')
    plt.ylabel(r'$M_{Y}$ [GeV]', fontsize=24)
    plt.xlabel(r'$M_{X}$ [GeV]', fontsize=24)

    # Add text annotations for Ratio
    for i in range(len(MYs)):
        for j in range(len(MXs)):
            value = Ratio[i][j]
            if not np.isnan(value):
                plt.text(x_edges[j] + (x_edges[j+1] - x_edges[j]) / 2,
                         y_edges[i] + (y_edges[i+1] - y_edges[i]) / 2,
                         f'{value:.1f}',
                         color='black',
                         ha='center',
                         va='center',
                         fontsize=8,
                         bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    plt.tight_layout()
    if expected:
        plt.savefig('Limit_Ratio_logscale.png', bbox_inches='tight')
        plt.savefig('Limit_Ratio_logscale.pdf', bbox_inches='tight')
    else:
        plt.savefig('Limit_Ratio_logscale_obs.png', bbox_inches='tight')
        plt.savefig('Limit_Ratio_logscale_obs.pdf', bbox_inches='tight')

def MakePlots(expected=True):
    limits = []
    ratios = []

    for MY in MYs:
        limits.append([])
        ratios.append([])    
        for MX in MXs:    
            # First set of files
            filepath_first = base_path + f'MX{MX}_MY{MY}-2_area/higgsCombineTest.AsymptoticLimits.mH120.root'
            
            if os.path.exists(filepath_first):
                limit_first = get_limit(filepath_first, fb_conversion_first, expected=expected)
                limits[-1].append(limit_first)
            else:
                print(f"No file found for MX= {MX}, MY= {MY} in first set")
                limits[-1].append(np.nan)

            # Second set of files
            filepath_second = second_base_path.format(MX=MX, MY=MY) + 'higgsCombine.AsymptoticLimits.mH125.42.root'
            
            if os.path.exists(filepath_second):
                limit_second = get_limit(filepath_second, fb_conversion_second, expected=expected)
            else:
                print(f"No file found for MX= {MX}, MY= {MY} in second set")
                limit_second = np.nan

            # Calculate the ratio
            if not math.isnan(limit_first) and not math.isnan(limit_second) and limit_second != 0:
                ratio = limit_first / limit_second
            else:
                ratio = np.nan
            ratios[-1].append(ratio)

    # Plot both limits and the ratio
    Plot_v2(limits, ratios, expected=expected)
    return
    

if __name__=='__main__':
    #MakePlots()
    #MakePlots(expected=False)
    write_limits_csv(filename="Limits_2Dgrid.csv",MXs=MXs,MYs=MYs,base_path=base_path,fb_conversion=5.)