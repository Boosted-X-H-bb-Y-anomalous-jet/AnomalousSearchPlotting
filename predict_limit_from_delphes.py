import pickle
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kFatal
import csv

def get_limit(process):
    base_dir = "/uscms/home/roguljic/nobackup/el8_anomalous/el9_fitting/CMSSW_14_1_0_pre4/src/AnomalousSearchFits/SR_run2"
    file_path = f"{base_dir}/{process}-2_area/higgsCombineTest.AsymptoticLimits.mH120.root"
    f = ROOT.TFile.Open(file_path)
    limit_tree = f.Get("limit")
    limit_tree.GetEntry(2)
    median_expected_limit = limit_tree.limit
    return median_expected_limit

def eff_ratio(procA, procB, filename="delphes_anomaly_tagging_efficiencies.csv"):
    #Returns efficiency of procB/ efficiency of procA
    efficiencies = {}
    with open(filename, 'r', newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        
        for row in reader:
            process_name = row[0].strip()
            efficiency = float(row[2])
            efficiencies[process_name] = efficiency
    
    return efficiencies[procB] / efficiencies[procA]

process_pairs = {}

for mx in ["1400","1600","1800","2000","2200","2600","3000"]:
    process_Yuu = f"XToYH_HTo2BYTo2Up_MX-{mx}_MY-200"
    process_Yuu_pair = f"MX{mx}_MY190"
    process_pairs[process_Yuu]=process_Yuu_pair

    process_Ytt = f"XToYH_HTo2BYTo2T_Hadronic_MX-{mx}_MY-400"
    process_Ytt_pair = f"MX{mx}_MY400"
    process_pairs[process_Ytt]=process_Ytt_pair 


for mx in ["1400","1600","1800","2000","2400","3000"]:
    process_Ybqq = f"TPrime_MX-{mx}_MY-125"
    process_Ybqq_pair = f"MX{mx}_MY190"
    if mx=="2400":
        process_Ybqq_pair = f"MX2200_MY190"
    process_pairs[process_Ybqq]=process_Ybqq_pair

import pickle
results = {}
xsec = 5.#fb
print("Process | Y->WW limit | Predicted from Delphes eff | True limit|")
for procB, procA in process_pairs.items():
    limitA = xsec*get_limit(procA)
    limitB = xsec*get_limit(procB)
    delphesRatio = eff_ratio(procA,procB)
    delphesPred = limitA*delphesRatio
    results[procB] = {"delphesPred": delphesPred, "limit": limitB}
    print(f"{procB} | {limitA:.3f} | {delphesPred:.3f} | {limitB:.3f}")

with open("limit_plotting/delphes_pred_limits.pkl", "wb") as f:
    pickle.dump(results, f)
