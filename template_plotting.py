import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

base_path = '/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/templates_v5/'

def get_histogram(process, year, region, variable):
    def add_histograms(counts_list, edges_list):
        total_counts = np.sum(counts_list, axis=0).tolist()
        return total_counts, edges_list[0]  # Assuming all histograms have the same binning

    years = ["2016APV", "2016", "2017", "2018"] if year == "run2" else [year]
    processes = ["QCD_HT700to1000", "QCD_HT1000to1500", "QCD_HT1500to2000", "QCD_HT2000toInf"] if process == "QCD" else [process]

    all_counts = []
    all_edges = None
    
    for yr in years:
        for proc in processes:
            file_name = f"{base_path}/templates_{proc}_{yr}.root"
            hist_name = f"mjj_my_{proc}_{yr}_{region}_nom"
            
            with uproot.open(file_name) as file:
                hist = file[hist_name]
                if variable == "mjj":
                    projection = hist.to_hist().project(0)
                elif variable == "mjy":
                    projection = hist.to_hist().project(1)
                else:
                    raise ValueError("Invalid variable: choose 'mjj' or 'mjy'.")

                counts = projection.values().tolist()
                edges = projection.axes.edges.tolist()

                all_counts.append(counts)
                if all_edges is None:
                    all_edges = edges

    if len(all_counts) > 1:
        return add_histograms(all_counts, all_edges)
    else:
        return all_counts[0], all_edges

def plot_histograms(processes, year, region, variable, colors=None):
    fig, (ax1, ax2) = plt.subplots(2, 1,sharex=True,sharey=False, figsize=(10, 8), gridspec_kw={'height_ratios': [4, 1], 'hspace': 0.1})
    hist_data = []
    bin_edges = None
    stacked_counts = None
    data_counts = None

    for process in processes:
        counts, edges = get_histogram(process, year, region, variable)
        
        if bin_edges is None:
            bin_edges = edges
        
        if process == "data_obs":
            data_counts = counts
        else:
            if stacked_counts is None:
                stacked_counts = np.array(counts)
            else:
                stacked_counts += np.array(counts)
            color = colors.get(process, None) if colors else None
            hist_data.append((np.array(edges), np.array(counts), color, process))

    if data_counts is not None:
        bin_centers = (np.array(bin_edges[:-1]) + np.array(bin_edges[1:])) / 2
        yerr = np.sqrt(np.array(data_counts))

    if hist_data:
        edges_list = [data[0] for data in hist_data]
        counts_list = [data[1] for data in hist_data]
        colors_list = [data[2] for data in hist_data]
        labels_list = [data[3] for data in hist_data]

        hep.histplot(counts_list,edges_list[0],color=colors_list,label=labels_list,histtype='fill',stack=True,ax=ax1)

    if data_counts is not None and region!="SR_Pass":
        ax1.errorbar(bin_centers, data_counts, yerr=yerr, fmt='o', color='black', label="data_obs")

    if variable == "mjj":
        ax1.set_xlim(1200, 3500)
        ax1.set_xlabel(r'$M_{jj}$ [GeV]')
    elif variable == "mjy":
        ax1.set_xlabel(r'$M_{j}^{Y}$ [GeV]')
        ax1.set_xlim(0, 500)

    ax1.set_yscale('log')
    ax1.set_ylabel("Events")
    ax1.set_xlabel("")
    ax1.legend(loc="upper right")
    ax1.text(0.15, 0.95, region.replace("_", " "), transform=ax1.transAxes, fontsize="large", verticalalignment='top')

    if "Pass" in region:
        ax1.set_ylim(1e-1, 1e4)
    else:
        ax1.set_ylim(1e-1, 1e6)

    if data_counts is not None and stacked_counts is not None:
        ratio = np.array(data_counts) / np.array(stacked_counts)
        ratio_error = np.sqrt(np.array(data_counts)) / np.array(stacked_counts)
        if region=="SR_Pass":
            ratio = np.zeros(len(ratio))
            ratio_error = np.zeros(len(ratio))
        ax2.errorbar(bin_centers, ratio, yerr=ratio_error, fmt='o', color='black')
        ax2.set_ylabel('Data / MC')
        ax2.set_xlabel(r'$M_{jj}$ [GeV]' if variable == "mjj" else r'$M_{j}^{Y}$ [GeV]')
        ax2.set_ylim(0, 2)
        ax2.axhline(y=1, color='grey', linestyle='--')

    ax1.tick_params(labelbottom=False)
    hep.cms.label("WiP", loc=0, ax=ax1)
    lumiText = f"{year if year != 'run2' else 'Run 2'} (13 TeV)"
    hep.cms.lumitext(lumiText, ax=ax1)

    plt.tight_layout()
    print(f"Saving plots/{region}_{year}_{variable}.png")
    plt.savefig(f"plots/{region}_{year}_{variable}.png")
    plt.savefig(f"plots/{region}_{year}_{variable}.pdf")


processes = ["TTToSemiLeptonic","TTToHadronic","QCD","data_obs"]
colors = {"QCD":"burlywood","TTToHadronic":"cornflowerblue","TTToSemiLeptonic":"darkblue","data_obs":"black"}
for region in ["CR_Pass","CR_Fail","SR_Pass","SR_Fail"]:
    plot_histograms(processes, "run2",region, "mjj", colors=colors)
    plot_histograms(processes, "run2",region, "mjy", colors=colors)