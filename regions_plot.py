import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('agg')  # or 'pdf'
import matplotlib.pyplot as plt
import mplhep as hep
import ROOT as r
r.gROOT.SetBatch(True) # Prevents ROOT from trying to open a window

def regions_plot():
    font_size = 20
    plt.rcParams.update({'font.size': font_size, 'axes.titlesize': font_size, 'axes.labelsize': font_size, 'xtick.labelsize': font_size, 'ytick.labelsize': font_size, 'legend.fontsize': font_size})
    vae_loss_range = np.linspace(0, 1, 6)  # from 0 to 1 scaled by 10^-4

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.fill_betweenx([0, 0.5], 0.05, 0.2, color='skyblue', alpha=0.5)
    ax.fill_betweenx([0.5, 1], 0.05, 0.2, color='dodgerblue', alpha=0.5)
    ax.fill_betweenx([0, 0.5], 0.2, 1, color='lightcoral', alpha=0.5)
    ax.fill_betweenx([0.5, 1], 0.2, 1, color='red', alpha=0.5)

    ax.text(0.12, 0.25, 'CR Fail', horizontalalignment='center', verticalalignment='center', fontsize=font_size)
    ax.text(0.12, 0.75, 'CR Pass', horizontalalignment='center', verticalalignment='center', fontsize=font_size)
    ax.text(0.45, 0.25, 'SR Fail', horizontalalignment='center', verticalalignment='center', fontsize=font_size)
    ax.text(0.45, 0.75, 'SR Pass', horizontalalignment='center', verticalalignment='center', fontsize=font_size)

    ax.set_xlabel('VAE Loss ($\\times 10^{-4}$)')
    ax.set_ylabel('ParticleNet Score')

    ax.set_xticks(vae_loss_range)
    ax.set_xticklabels([f'{x:.1f}' for x in vae_loss_range])

    ax.set_yticks([0.25, 0.75])
    ax.set_yticklabels(['Fail', 'Pass'])

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig("plots/regions.pdf")

def vae_plot():

    def vae_hist_from_file(file_name,hist_name,rebin=1,norm=False):
        f = r.TFile.Open(file_name)
        h = f.Get(hist_name)
        h = h.ProjectionZ(hist_name+"_projz")
        h.RebinX(rebin)
        if norm:
            h.Scale(1./h.Integral())
        values,errors,edges = hist_to_numpy(h)
        f.Close()
        return values,errors,edges

    def hist_to_numpy(h):
        values = []
        errors = []
        edges  = []
        n_bins = h.GetNbinsX()
        for i in range(n_bins):
            values.append(h.GetBinContent(i+1))        
            errors.append(h.GetBinError(i+1))
            edges.append(h.GetBinLowEdge(i+1))
        
        edges.append(h.GetBinLowEdge(n_bins)+h.GetBinWidth(n_bins))
        return values, errors, edges

    def non_stack_plot(histos,edges,labels,out_file,colors=[],x_title="X",y_title="Y", alpha=0.7,x_lim=[None,None],y_lim=[None,None], year='', category=''):
        plt.style.use([hep.style.CMS])
        hep.histplot(histos,edges[0],stack=False,label=labels,histtype="step",color=colors,alpha=alpha,linewidth=3)
        plt.xlabel(x_title, horizontalalignment='center', x=0.5)
        plt.ylabel(y_title, horizontalalignment='center', y=0.5)
        plt.xlim(x_lim)
        plt.ylim(y_lim)
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.legend()
        plt.text(0.55, 0.95, f"{year} PNetXbb {category}", transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='right')
        print(f"Saving {out_file}")
        plt.savefig(out_file)
        plt.cla()
        plt.clf()

    base_path='/uscms_data/d3/roguljic/el8_anomalous/el9_fitting/templates_v3/templates_'
    for year in ["2016APV","2016","2017","2018"]:
        for category in ["Pass","Fail"]:
            labels = [r"t$\bar{t}$","Data","Signal"]
            histo_names = []
            histos = []
            edges_all = []
            rebin = 1
            for process in ["TTToHadronic","data_obs","MX1600_MY90"]:
                filepath=base_path+process+'_'+year+'.root'
                keyname=f'mjj_my_vaeloss_{process}_{year}_IR_{category}_nom'
                values,errors,edges = vae_hist_from_file(filepath,keyname,rebin=rebin,norm=True)
                histos.append(values)
                edges_all.append(edges)

            vae_plot_name = f"plots/VAE_loss_IR_{category}_{year}.png"
            non_stack_plot(histos,edges_all,labels,vae_plot_name,x_title="VAE loss",y_title="Normalized yield", alpha=0.7,x_lim=[0.,0.0001], year=year, category=category)

vae_plot()
regions_plot()