"""Shared figure style: L-frame axes, outward ticks, no grid, Okabe-Ito colours."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 8, 'axes.labelsize': 9, 'legend.fontsize': 7, 'axes.titlesize': 8,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.8,
                     'xtick.direction': 'out', 'ytick.direction': 'out', 'axes.grid': False,
                     'savefig.dpi': 300, 'savefig.bbox': 'tight'})
C = {'black': '#000000', 'orange': '#E69F00', 'blue': '#0072B2', 'green': '#009E73',
     'vermillion': '#D55E00', 'sky': '#56B4E9', 'purple': '#CC79A7'}


def subplots(*a, **k):
    return plt.subplots(*a, **k)


def panel(ax, letter, x=-0.14, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=11, fontweight='bold', va='bottom')


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
