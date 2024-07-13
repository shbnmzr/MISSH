import matplotlib.pyplot as plt
import numpy as np

species = ("W26L31", "W32L45", "W14L31", "W18L31", "W22L31", "W10L15")

penguin_means = {
    'ISSH': (2308.27, 3754.86, 2514.23, 3110.68, 2808.64, 1693.48),
    'Set Covering (CPLEX)': (125.82, 130.98, 122.88, 123.16, 124.01, 116.56),
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
ax.set_ylim(0, 4000)

plt.show()
