import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

base_path = '/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/templates_v6/'

def add_histograms(counts_list, edges_list):
    total_counts = np.sum(counts_list, axis=0).tolist()
    return total_counts, edges_list#[0]

def get_histogram(process, year, region, variable):
    years = ["2016APV", "2016", "2017", "2018"] if year == "run2" else [year]
    processes = ["QCD_HT700to1000", "QCD_HT1000to1500", "QCD_HT1500to2000", "QCD_HT2000toInf"] if process == "QCD" else [process]

    all_counts = []
    all_edges = None
    
    for yr in years:
        for proc in processes:
            file_name = f"{base_path}/templates_{proc}_{yr}.root"
            if variable in ["mjj", "mjy"]:
                hist_name = f"mjj_my_{proc}_{region}_nom"
            elif variable == "phi":
                hist_name = f"phi_{proc}_{region}_nom"
            elif variable == "eta":
                hist_name = f"eta_{proc}_{region}_nom"
            elif variable == "eta_phi":
                hist_name = f"eta_phi_{proc}_{region}_nom"
            else:
                raise ValueError("Invalid variable: choose 'mjj', 'mjy', 'phi', 'eta' or 'eta_phi'.")

            with uproot.open(file_name) as file:
                hist = file[hist_name]
                
                #For mjj/mjy we need to project from 2D histogram
                if variable in ["mjj", "mjy"]:
                    projection = hist.to_hist().project(0) if variable == "mjj" else hist.to_hist().project(1)
                    counts = projection.values().tolist()
                    edges = projection.axes.edges.tolist()
                else:
                    counts = hist.to_hist().values().tolist()
                    edges = hist.to_hist().axes.edges.tolist()
                all_counts.append(counts)
                if all_edges is None:
                    all_edges = edges

    if len(all_counts) > 1:
        return add_histograms(all_counts, all_edges)
    else:
        return all_counts[0], all_edges

def plot_histograms(processes, year, region, variable, colors=None):
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=False, figsize=(10, 8), gridspec_kw={'height_ratios': [4, 1], 'hspace': 0.1})
    hist_data = []
    bin_edges = None
    stacked_counts = None
    data_counts = None

    for process in processes:
        counts, edges = get_histogram(process, year, region, variable)
        edges = edges[0]
        
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

        hep.histplot(counts_list, edges_list[0], color=colors_list, label=labels_list, histtype='fill', stack=True, ax=ax1)

    if data_counts is not None and region != "SR_Pass":
        ax1.errorbar(bin_centers, data_counts, yerr=yerr, fmt='o', color='black', label="data_obs")

    # Set labels and limits based on the variable being plotted
    if variable == "mjj":
        ax1.set_xlim(1300, 3500)
        ax1.set_xlabel(r'$M_{jj}$ [GeV]')
    elif variable == "mjy":
        ax1.set_xlim(0, 500)
        ax1.set_xlabel(r'$M_{j}^{Y}$ [GeV]')
    elif variable == "phi":
        ax1.set_xlim(-3.2, 3.2)
        ax1.set_xlabel(r'$\phi$')
    elif variable == "eta":
        ax1.set_xlim(-2.5, 2.5)
        ax1.set_xlabel(r'$\eta$')

    ax1.set_yscale('log')
    ax1.set_ylabel("Events")
    ax1.set_xlabel("")
    ax1.legend(loc="upper right")
    ax1.text(0.15, 0.95, region.replace("_", " "), transform=ax1.transAxes, fontsize="large", verticalalignment='top')

    if region=="IR_Pass":
        ax1.set_ylim(1, 1e5)
    elif region=="IR_Fail":
        ax1.set_ylim(10, 1e8)
    elif "Pass" in region:
        ax1.set_ylim(1e-1, 1e4)
    else:
        ax1.set_ylim(1e-1, 1e7)

    if data_counts is not None and stacked_counts is not None:
        ratio = np.array(data_counts) / np.array(stacked_counts)
        ratio_error = np.sqrt(np.array(data_counts)) / np.array(stacked_counts)
        if region == "SR_Pass":
            ratio = np.zeros(len(ratio))
            ratio_error = np.zeros(len(ratio))
        ax2.errorbar(bin_centers, ratio, yerr=ratio_error, fmt='o', color='black')
        ax2.set_ylabel('Data / MC')
        ax2.set_xlabel(r'$M_{jj}$ [GeV]' if variable == "mjj" else r'$M_{j}^{Y}$ [GeV]' if variable == "mjy" else r'$\phi$' if variable == "phi" else r'$\eta$')
        ax2.set_ylim(0, 2)
        ax2.axhline(y=1, color='grey', linestyle='--')

    ax1.tick_params(labelbottom=False)
    hep.cms.label("WiP", loc=0, ax=ax1,data=True)
    #lumiText = f"{year if year != 'run2' else 'Run 2'} (13 TeV)"
    lumiText = f"{year if year != 'run2' else 'Run 2'} (13 TeV)"
    hep.cms.lumitext(lumiText, ax=ax1)

    plt.tight_layout()
    print(f"Saving plots/{region}_{year}_{variable}.png")
    plt.savefig(f"plots/{region}_{year}_{variable}.png")
    plt.savefig(f"plots/{region}_{year}_{variable}.pdf")


def plot_2d_histogram(processes, year, region):
    for process in processes:
        fig, ax = plt.subplots(figsize=(10, 8))
        #print(process)
        counts, edges = get_histogram(process, year, region, variable='eta_phi')
        counts_array = np.array(counts)
        #print(np.shape(counts_array))
        #print(np.shape(edges))
        #This will make bins with zero entries appear as white
        counts_array_masked = np.ma.masked_where(counts_array == 0, counts_array)
        
        edges_eta = [edge for sublist in edges[0] for edge in sublist]
        edges_phi = [edge for sublist in edges[1] for edge in sublist]
            
        c = ax.pcolormesh(edges_eta, edges_phi,counts_array_masked.T, shading='auto', cmap='viridis')
        fig.colorbar(c, ax=ax, label='Counts')

        ax.set_ylabel(r'$\phi$')
        ax.set_xlabel(r'$\eta$')
        data_flag = "False"
        if("data" in process):
            data_flag = True
        hep.cms.label("WiP", loc=0, ax=ax,data=data_flag)
        lumiText = f"{year if year != 'run2' else 'Run 2'} (13 TeV)"
        hep.cms.lumitext(lumiText, ax=ax)

        plt.tight_layout()
        print(f"Saving plots/{process}_{region}_{year}_eta_phi.png")
        plt.savefig(f"plots/{process}_{region}_{year}_eta_phi.png")
        plt.savefig(f"plots/{process}_{region}_{year}_eta_phi.pdf")

        ax.cla()
        plt.clf() 


def get_histogram_ratio(process, year, region1, region2, variable):
    import ROOT
    years = ["2016APV", "2016", "2017", "2018"] if year == "run2" else [year]
    processes = ["QCD_HT700to1000", "QCD_HT1000to1500", "QCD_HT1500to2000", "QCD_HT2000toInf"] if process == "QCD" else [process]

    total_pass = None
    total_fail = None
    
    for yr in years:
        for proc in processes:
            file_name = f"{base_path}/templates_{proc}_{yr}.root"

            if variable in ["mjj", "mjy"]:
                hist_name1 = f"mjj_my_{proc}_{region1}_nom"
                hist_name2 = f"mjj_my_{proc}_{region2}_nom"
            elif variable == "phi":
                hist_name1 = f"phi_{proc}_{region1}_nom"
                hist_name2 = f"phi_{proc}_{region2}_nom"
            elif variable == "eta":
                hist_name1 = f"eta_{proc}_{region1}_nom"
                hist_name2 = f"eta_{proc}_{region2}_nom"
            else:
                raise ValueError("Invalid variable: choose 'mjj', 'mjy', 'phi', 'eta'.")

            file = ROOT.TFile.Open(file_name)

            hist_pass = file.Get(hist_name1)
            hist_fail = file.Get(hist_name2)

            if variable == "mjj":
                hist_pass = hist_pass.ProjectionX()
                hist_fail = hist_fail.ProjectionX()
            elif variable == "mjy":
                hist_pass = hist_pass.ProjectionY()
                hist_fail = hist_fail.ProjectionY()

            if total_pass is None:
                total_pass = hist_pass.Clone(f"total_{region1}")
                total_fail = hist_fail.Clone(f"total_{region2}")
                total_pass.SetDirectory(0)
                total_fail.SetDirectory(0)
            else:
                total_pass.Add(hist_pass)
                total_fail.Add(hist_fail)

            file.Close()

    hist_ratio = total_pass.Clone(f"ratio_{region1}_{region2}")
    hist_ratio.Divide(total_fail)
    
    return hist_ratio


def plot_qcd_ratio(year, variable, region1="CR_Pass", region2="CR_Fail"):
    hist_ratio = get_histogram_ratio("QCD", year, region1, region2, variable)
    
    bin_centers = np.array([hist_ratio.GetBinCenter(i) for i in range(1, hist_ratio.GetNbinsX() + 1)])
    ratio_values = np.array([hist_ratio.GetBinContent(i) for i in range(1, hist_ratio.GetNbinsX() + 1)])
    ratio_errors = np.array([hist_ratio.GetBinError(i) for i in range(1, hist_ratio.GetNbinsX() + 1)])
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.errorbar(bin_centers, ratio_values, yerr=ratio_errors, fmt='o', label=f'QCD {region1}/{region2}')
    
    if variable == "mjj":
        ax.set_xlim(1300, 3500)
        ax.set_xlabel(r'$M_{jj}$ [GeV]')
    elif variable == "mjy":
        ax.set_xlim(0, 500)
        ax.set_xlabel(r'$M_{j}^{Y}$ [GeV]')
        
    ax.set_ylabel(r'$R_{P/F}$')
    ax.set_ylim(0, 0.02)
    
    hep.cms.label("WiP", loc=0, ax=ax, data=False)
    lumiText = f"{year if year != 'run2' else 'Run 2'} (13 TeV)"
    hep.cms.lumitext(lumiText, ax=ax)
    
    plt.tight_layout()
    plt.savefig(f"plots/QCD_ratio_{region1}_{region2}_{year}_{variable}.png")
    plt.savefig(f"plots/QCD_ratio_{region1}_{region2}_{year}_{variable}.pdf")
    plt.clf()


# Example calls
plot_qcd_ratio("2018", "mjj")
plot_qcd_ratio("2018", "mjy")

exit()





processes = ["TTToSemiLeptonic", "TTToHadronic", "QCD", "data_obs"]
colors = {"QCD": "burlywood", "TTToHadronic": "cornflowerblue", "TTToSemiLeptonic": "darkblue", "data_obs": "black"}


for region in ["CR_Pass", "CR_Fail", "SR_Pass", "SR_Fail"]:
    # Plot mjj and mjy histograms
    plot_histograms(processes, "run2", region, "mjj", colors=colors)
    plot_histograms(processes, "run2", region, "mjy", colors=colors)

for region in ["IR_Pass", "IR_Fail"]:
    # Plot phi and eta histograms
    plot_histograms(processes, "run2", region, "phi", colors=colors)
    plot_histograms(processes, "run2", region, "eta", colors=colors)
    plot_2d_histogram(processes, "2018", "IR_Fail")
    plot_2d_histogram(["MX2600_MY250"], "2018", region)
