import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

mpl.rc('text', usetex=True)
mpl.rcParams.update({'font.size': 16})

main_path = "parallel_runs/bifurcation_file.txt"
df = pd.read_csv(main_path, sep='\t', names=['Gamma', 'Ra_c'])

def jeffrey(gamma_values):
    Ra_values = []
    for gamma in gamma_values:
        Ra_values.append(1708)
    return Ra_values

def jeffrey_corrected(gamma_values):
    Ra_values = []
    for gamma in gamma_values:
        Ra_values.append(1708*(1+0.52*gamma**(-2))**2)
    return Ra_values

def olga_approach(gamma_values):
    Ra_values = []
    for gamma in gamma_values:
        Ra_values.append((2*math.pi)**4*(1+gamma**(-2))*(1+gamma**(-2)/4))
    return Ra_values

fig, ax = plt.subplots()
ax.plot(df['Gamma'].values, df["Ra_c"].values, label="Present sims.", color='blue')
ax.plot(df['Gamma'].values, jeffrey_corrected(df['Gamma'].values), label="1708(1+0.52$\Gamma^{-2}$)$^{-2}$", color='black', linestyle='--')
ax.plot(df['Gamma'].values, olga_approach(df['Gamma'].values), label="(2$\pi$)$^4$(1+$\Gamma^{-2}$)(1+$\Gamma^{-2}$/4))", color='red')
ax.plot(df['Gamma'].values, jeffrey(df['Gamma'].values), label="1708", color='lightblue', linestyle='--')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('$\Gamma$')
ax.set_ylabel('Ra$_{c,\Gamma}$')
plt.legend()
plt.show()