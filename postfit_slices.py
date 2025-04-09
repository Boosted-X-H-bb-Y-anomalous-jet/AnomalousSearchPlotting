from plots_for_paper import *

def plot_slices(histos_dict, region, processes, labels_dict, colors_dict, axis="X", yRangeLog=[1, 10**5]):
    if axis not in ["X", "Y"]:
        raise ValueError("Invalid axis. Choose 'X' or 'Y'.")

    axis_label = "$M_{JJ} [GeV]$" if axis == "X" else "$M_{J}^{Y} [GeV]$"
    file_suffix = "mjj" if axis == "X" else "mjy"
    slice_var = "$M_{J}^{Y}$" if axis == "X" else "$M_{JJ}$"
    x_range = [1300, 3500] if axis == "X" else [40, 500]

    histos = histos_dict[region]

    slicing_axis = histos["data_obs"].GetYaxis() if axis == "X" else histos["data_obs"].GetXaxis()
    projection_method = f"Projection{axis}"
    n_bins = slicing_axis.GetNbins()

    h_bkg = histos["TotalBkg"]

    # Slice and plot
    for bin_idx in range(1, n_bins + 1):
        h_data_slice = getattr(histos["data_obs"], projection_method)(f"data_{axis}_slice_{bin_idx}", bin_idx, bin_idx)
        for i in range(1,h_data_slice.GetNbinsX()+1):
            #Sometimes float precision messes up the counts of data events so we force them to be int
            h_data_slice.SetBinContent(i,round(h_data_slice.GetBinContent(i)))
        h_data_slice.SetBinErrorOption(1)
        h_data_slice = PyHist(h_data_slice)
        h_data_slice.divide_by_bin_width()

        yRange=[0,1.5*max(h_data_slice.bin_values)]

        h_mc_slices = []
        labels_mc = []
        colors_mc = []
        for process in processes:
            h_temp = getattr(histos[process], projection_method)(f"{process}_{axis}_slice_{bin_idx}", bin_idx, bin_idx)
            h_temp = PyHist(h_temp)
            h_temp.divide_by_bin_width()
            h_mc_slices.append(h_temp)
            labels_mc.append(labels_dict[process])
            colors_mc.append(colors_dict[process])

        h_bkg_slice = getattr(h_bkg, projection_method)(f"bkg_{axis}_slice_{bin_idx}", bin_idx, bin_idx)
        h_bkg_slice = PyHist(h_bkg_slice)
        h_bkg_slice.divide_by_bin_width()

        bin_low_edge = slicing_axis.GetBinLowEdge(bin_idx)
        bin_up_edge = slicing_axis.GetBinUpEdge(bin_idx)

        # Format the bin range as text
        projectionText = f"{region}\n{bin_low_edge:.0f}<{slice_var}<{bin_up_edge:.0f} GeV"

        plotShapesWithRatioAndBand(h_data_slice, h_mc_slices, h_bkg_slice, labels_mc, colors_mc,axis_label, f"slices/{region}_{file_suffix}_slice_{bin_idx}.png",xRange=x_range, yRange=yRange, projectionText=projectionText, yRangeLog=yRangeLog)


if __name__ == '__main__':
        input_file = "~/nobackup/el8_anomalous/el9_fitting/CMSSW_14_1_0_pre4/src/AnomalousSearchFits/SR_run2_one_signal/MX2200_MY250-2_area/postfitshapes_b.root"
        processes_SR_Pass =["GluGluHToBB","VBFHToBB","WplusH_HToBB_WToQQ","WminusH_HToBB_WToQQ","ZH_HToBB_ZToQQ","ggZH_HToBB_ZToQQ","TTToHadronic","TTToSemiLeptonic","Background_2","TotalBkg","data_obs"]
        processes_SR_Fail =["GluGluHToBB","VBFHToBB","WplusH_HToBB_WToQQ","WminusH_HToBB_WToQQ","ZH_HToBB_ZToQQ","ggZH_HToBB_ZToQQ","TTToHadronic","TTToSemiLeptonic","Background","TotalBkg","data_obs"]
        histos_dict={}
        histos_dict["SR_Pass"] = get_hists(input_file,"SR_Pass",processes_SR_Pass)
        #histos_dict["SR_Fail"] = get_hists(input_file,"SR_Fail",processes_SR_Fail)
        labels_dict={"GluGluHToBB":"SM Higgs","VBFHToBB":"__nolabel__","WplusH_HToBB_WToQQ":"__nolabel__","WminusH_HToBB_WToQQ":"__nolabel__","ZH_HToBB_ZToQQ":"__nolabel__","ggZH_HToBB_ZToQQ":"__nolabel__","TTToHadronic":r"$t\bar t$","TTToSemiLeptonic":"__nolabel__","Background_2":"Multijet","Background":"Multijet"}
        colors_dict={"GluGluHToBB":"forestgreen","VBFHToBB":"forestgreen","WplusH_HToBB_WToQQ":"forestgreen","WminusH_HToBB_WToQQ":"forestgreen","ZH_HToBB_ZToQQ":"forestgreen","ggZH_HToBB_ZToQQ":"forestgreen","TTToHadronic":"#d42e12","TTToSemiLeptonic":"#d42e12","Background_2":"#f39c12","Background":"#f39c12"}

        plot_slices(histos_dict,"SR_Pass",["GluGluHToBB","VBFHToBB","WplusH_HToBB_WToQQ","WminusH_HToBB_WToQQ","ZH_HToBB_ZToQQ","ggZH_HToBB_ZToQQ","Background_2","TTToHadronic","TTToSemiLeptonic"],labels_dict,colors_dict,axis="X",yRangeLog=[0.0001,10**1])    
        #plot_slices(histos_dict,"SR_Pass",["GluGluHToBB","VBFHToBB","WplusH_HToBB_WToQQ","WminusH_HToBB_WToQQ","ZH_HToBB_ZToQQ","ggZH_HToBB_ZToQQ","Background_2","TTToHadronic","TTToSemiLeptonic"],labels_dict,colors_dict,axis="Y",yRangeLog=[0.0001,10**1])
        #plot_slices(histos_dict,"SR_Fail",["GluGluHToBB","VBFHToBB","WplusH_HToBB_WToQQ","WminusH_HToBB_WToQQ","ZH_HToBB_ZToQQ","ggZH_HToBB_ZToQQ","Background","TTToHadronic","TTToSemiLeptonic"],labels_dict,colors_dict,axis="X",yRangeLog=[0.0001,10**1])
        #plot_slices(histos_dict,"SR_Fail",["GluGluHToBB","VBFHToBB","WplusH_HToBB_WToQQ","WminusH_HToBB_WToQQ","ZH_HToBB_ZToQQ","ggZH_HToBB_ZToQQ","Background","TTToHadronic","TTToSemiLeptonic"],labels_dict,colors_dict,axis="Y",yRangeLog=[0.0001,10**1])
