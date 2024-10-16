import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Arc
import numpy as np

# Create a new figure
fig, ax = plt.subplots()

# Set limits
ax.set_xlim(-1, 3)
ax.set_ylim(-2, 2)

# Remove axes
ax.axis('off')

# Scalar X decays at origin (0, 0)
ax.text(0, 0, 'X', fontsize=15, ha='center', va='center')

# Scalar Y, H outgoing
ax.text(2, 1, 'Y', fontsize=15, ha='center', va='center')
ax.text(2, -1, 'H', fontsize=15, ha='center', va='center')

# Draw lines for particles
ax.add_patch(FancyArrow(0, 0, 2, 1, width=0.05, head_width=0.1, head_length=0.1, color='k'))
ax.add_patch(FancyArrow(0, 0, 2, -1, width=0.05, head_width=0.1, head_length=0.1, color='k'))

# Save and show the diagram
plt.savefig("scalar_decay.png")
