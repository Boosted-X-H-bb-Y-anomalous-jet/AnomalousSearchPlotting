import pickle
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kFatal

def get_efficiencies_from_pkl(process):
    lumis = {"2016":16.8,"2016APV":19.5,"2017":41.5,"2018":59.8}
    with open("/uscms/home/roguljic/nobackup/el8_anomalous/systematics/gen_efficiencies/efficiencies.pkl", "rb") as f:
        gen_reco_efficiencies = pickle.load(f)
    reco = sum(gen_reco_efficiencies["reco"][process][year] * lumis[year] for year in lumis)/sum(lumis[year] for year in lumis)
    gen = sum(gen_reco_efficiencies["gen"][process][year] * lumis[year] for year in lumis) / sum(lumis[year] for year in lumis)
    return reco, gen

def eff_ratios(processA,processB):
    #Returns reco- and gen-efficiency A/B ratios
    recoA, genA = get_efficiencies_from_pkl(processA)
    recoB, genB = get_efficiencies_from_pkl(processB)
    recoRatio = recoA/recoB
    genRatio = genA/genB
    return recoRatio,genRatio


def eff_ratios_from_tpl(processA,processB):
    #Returns reco- and gen-efficiency A/B ratios
    effA = get_full_reco_efficiency(processA)
    effB = get_full_reco_efficiency(processB)
    ratio = effA/effB
    return ratio

def get_limit(process):
    base_dir = "/uscms/home/roguljic/nobackup/el8_anomalous/el9_fitting/CMSSW_14_1_0_pre4/src/AnomalousSearchFits/SR_run2"
    file_path = f"{base_dir}/{process}-2_area/higgsCombineTest.AsymptoticLimits.mH120.root"
    f = ROOT.TFile.Open(file_path)
    limit_tree = f.Get("limit")
    limit_tree.GetEntry(2)
    median_expected_limit = limit_tree.limit
    return median_expected_limit

def predict_limit(limitA,recoRatio,genRatio,fullRecoRatio):
    #Predict limitB from reco and genRatios (A/B) and limitA
    recoPred = limitA*recoRatio
    genPred = limitA*genRatio
    fullRecoPred = limitA*fullRecoRatio
    return recoPred,genPred,fullRecoPred


def get_full_reco_efficiency(process):
    file_path = f"/uscms/home/roguljic/nobackup/el8_anomalous/el9_fitting/templates_v8/templates_{process}_run2.root"
    hist_name = f"mjj_my_{process}_SR_Pass_nom"
    int_lumi = 138.
    xsec = 5.
    
    root_file = ROOT.TFile.Open(file_path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"File not found or corrupted: {file_path}")
        return None
    
    hist = root_file.Get(hist_name)
    if not hist:
        print(f"Histogram not found: {hist_name}")
        return None
    
    total_yield = hist.Integral(7, 50, 5, 50)
    efficiency = total_yield / (int_lumi * xsec)
    
    root_file.Close()
    return efficiency

process_pairs = {"XToYH_HTo2BYTo2Up_MX-2000_MY-200":"MX2000_MY190","XToYH_HTo2BYTo2T_MX-2000_MY-400":"MX2000_MY400","TPrime_MX-1400_MY-125":"MX1400_MY190","TPrime_MX-1600_MY-125":"MX1600_MY190","TPrime_MX-1800_MY-125":"MX1800_MY190","TPrime_MX-2000_MY-125":"MX2000_MY190","TPrime_MX-2400_MY-125":"MX2600_MY190","TPrime_MX-3000_MY-125":"MX3000_MY190"}
import pickle
results = {}
xsec = 5.#fb
print("Process | Y->WW limit | Predicted from Gen eff | Predicted from Reco eff | Predicted from full Reco eff | True limit|")
for procB, procA in process_pairs.items():
    limitA = xsec*get_limit(procA)
    limitB = xsec*get_limit(procB)
    recoRatio,genRatio = eff_ratios(procA,procB)
    fullRecoRatio = eff_ratios_from_tpl(procA,procB)
    recoPred,genPred,fullRecoPred=predict_limit(limitA,recoRatio,genRatio,fullRecoRatio)
    print(f"{procB} | {limitA:.3f} | {recoPred:.3f} | {genPred:.3f} | {fullRecoPred:.3f} | {limitB:.3f}")
    results[procB] = {"recoPred": recoPred, "genPred": genPred, "limit": limitB}

with open("limit_plotting/pred_limits.pkl", "wb") as f:
    pickle.dump(results, f)
