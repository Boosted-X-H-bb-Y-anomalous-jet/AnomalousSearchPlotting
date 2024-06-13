import numpy as np
import matplotlib.pyplot as plt
font_size = 20
plt.rcParams.update({'font.size': font_size, 'axes.titlesize': font_size, 'axes.labelsize': font_size, 'xtick.labelsize': font_size, 'ytick.labelsize': font_size, 'legend.fontsize': font_size})
vae_loss_range = np.linspace(0, 1, 6)  # from 0 to 1 scaled by 10^-4

fig, ax = plt.subplots(figsize=(10, 6))

ax.fill_betweenx([0, 0.5], 0.25, 0.4, color='skyblue', alpha=0.5)
ax.fill_betweenx([0.5, 1], 0.25, 0.4, color='dodgerblue', alpha=0.5)
ax.fill_betweenx([0, 0.5], 0.5, 1, color='lightcoral', alpha=0.5)
ax.fill_betweenx([0.5, 1], 0.5, 1, color='red', alpha=0.5)

ax.text(0.32, 0.25, 'CR Fail', horizontalalignment='center', verticalalignment='center', fontsize=font_size)
ax.text(0.32, 0.75, 'CR Pass', horizontalalignment='center', verticalalignment='center', fontsize=font_size)
ax.text(0.75, 0.25, 'SR Fail', horizontalalignment='center', verticalalignment='center', fontsize=font_size)
ax.text(0.75, 0.75, 'SR Pass', horizontalalignment='center', verticalalignment='center', fontsize=font_size)

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
