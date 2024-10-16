import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use(hep.style.CMS)

def parse_data(filename):
    data = {}
    with open(filename) as f:
        for line in f:
            parts = line.split(',')
            year, process = parts[0], parts[1]
            nominal, up, down = map(float, parts[2:5])
            if year not in data:
                data[year] = []
            data[year].append((process, nominal, up, down))
    return data

def plot_MX_MY90(data, year):
    MX_vals, noms, ups, downs = [], [], [], []
    for process, nominal, up, down in data[year]:
        if 'MY90' in process:
            MX = int(process.split('MX')[-1].split('_')[0])
            MX_vals.append(MX)
            noms.append(nominal)
            ups.append(up)
            downs.append(down)
    MX_vals, noms, ups, downs = map(np.array, (MX_vals, noms, ups, downs))
    plt.errorbar(MX_vals, noms, yerr=[downs, ups], fmt='o')
    plt.xlabel('MX')
    plt.ylabel('Scale Factor')
    hep.cms.text("WiP",loc=0)
    lumiText = f"{year} (13 TeV)"    
    hep.cms.lumitext(lumiText)
    plt.ylim(0.5, 1.5)

    plt.savefig(f'plot_{year}_MX_MY90.png')
    plt.savefig(f'plot_{year}_MX_MY90.pdf')
    plt.clf()

def plot_MY_MX2200(data, year):
    MY_vals, noms, ups, downs = [], [], [], []
    for process, nominal, up, down in data[year]:
        if 'MX2200' in process:
            MY = int(process.split('MY')[-1])
            MY_vals.append(MY)
            noms.append(nominal)
            ups.append(up)
            downs.append(down)
    MY_vals, noms, ups, downs = map(np.array, (MY_vals, noms, ups, downs))
    plt.errorbar(MY_vals, noms, yerr=[downs, ups], fmt='o')
    plt.xlabel('MY')
    plt.ylabel('Scale Factor')
    hep.cms.text("WiP",loc=0)
    lumiText = f"{year} (13 TeV)"    
    hep.cms.lumitext(lumiText)
    plt.ylim(0.5, 1.5)
    plt.savefig(f'plot_{year}_MY_MX2200.png')
    plt.savefig(f'plot_{year}_MY_MX2200.pdf')
    plt.clf()

def plot_TT(process_name, data):
    years = ['2016APV', '2016', '2017', '2018']
    years_ordered = np.arange(len(years))
    noms, ups, downs = [], [], []
    for year in years:
        for process, nominal, up, down in data[year]:
            if process == process_name:
                noms.append(nominal)
                ups.append(up)
                downs.append(down)
    noms, ups, downs = map(np.array, (noms, ups, downs))
    plt.errorbar(years_ordered, noms, yerr=[downs, ups], fmt='o')
    plt.xlabel('Year')
    plt.ylabel('Scale Factor')
    hep.cms.text("WiP",loc=0)
    lumiText = f"{year} (13 TeV)"    
    hep.cms.lumitext(lumiText)
    plt.xticks(ticks=years_ordered, labels=years)
    plt.ylim(0.5, 1.5)
    plt.savefig(f'plot_{process_name}.png')
    plt.savefig(f'plot_{process_name}.pdf')
    plt.clf()

data = parse_data('SFs.txt')

for year in ['2016APV', '2016', '2017', '2018']:
    plot_MX_MY90(data, year)
    plot_MY_MX2200(data, year)

plot_TT('TTToHadronic', data)
plot_TT('TTToSemiLeptonic', data)
