import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# +30% on top of previous 1.30 => 1.69 total
sns.set_context("paper", font_scale=1.69)

df = pd.read_csv('differentially_present_binary.csv', index_col=0)

# 0 -> lavender, 1 -> sea-foam green
cmap = sns.color_palette(["#C5B4E3", "#9FD8C0"])

# Increase figure WIDTH by 20% to widen cells (16 -> 19.2)
plt.figure(figsize=(19.2, len(df) * 2.4), dpi=600)

ax = sns.heatmap(
    df.T,
    cmap=cmap,
    vmin=0, vmax=1,
    cbar=False,
    linewidths=0.4
)

plt.xlabel('')
plt.ylabel('')
ax.xaxis.tick_top()
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

plt.subplots_adjust(left=0.2, right=0.85, top=0.95, bottom=0.12)
plt.tight_layout()
plt.savefig('differentially_present_binary.jpeg', format='jpeg', dpi=600, bbox_inches='tight')
plt.show()
