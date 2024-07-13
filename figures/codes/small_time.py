import matplotlib.pyplot as plt
import numpy as np

species = ("W26L31", "W32L45", "W14L31", "W18L31", "W22L31", "W10L15")

penguin_means = {
    'ISSH': (169.13, 912.43, 170.77, 473.27, 284.4, 25.72),
    'Set Covering (CPLEX)': (96.67, 95.03, 96.6, 96.76, 96.57, 98.17),

}

x = np.arange(len(species))  # the label locations
width = 0.3
padding = 0.1  # padding between groups of bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in penguin_means.items():
    offset = width * multiplier + padding * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=3)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Execution Time (milliseconds)')
ax.set_xticks(x + width + padding / 2, species)
ax.legend(loc='upper right', ncols=2)
ax.set_ylim(0, 1000)

plt.show()
