import sys
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

from plot_utils import load_matrix

def plot_matrix_heatmaps(paths, figsize_per_plot: tuple = (4, 4)) -> plt.Figure:
    if isinstance(paths, str):
        paths = [paths]

    n_plots = len(paths)
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(figsize_per_plot[0] * n_cols + 1, figsize_per_plot[1] * n_rows),
        squeeze=False
    )

    # Load all matrices first so we can compute the shared color range
    data = np.abs(np.array([load_matrix(p) for p in paths]))

    temp = np.ma.masked_less_equal(data, 0.0)
    temp = np.ma.masked_where(~np.isfinite(temp), temp)
    temp = np.log10(temp)
    min_val = np.nanmin(temp)
    max_val = np.nanmax(temp) 

    cmap = mpl.cm.get_cmap('cividis')
    normalizer = mpl.colors.Normalize(min_val, max_val)

    for i, path in enumerate(paths):
        ax = axes[i // n_cols][i % n_cols]
        plot_data = np.log10(np.ma.masked_less_equal(data[i], 0.0))
        im = ax.imshow(plot_data, cmap=cmap, norm=normalizer)
        
        parts = str(Path(path).parent).split('/')[2:]
        parts = [p if len(p) <= 35 else p[:18] + ' [...] ' + p[-17:]for p in parts]
        title = '\n'.join(parts)
        ax.set_title(title, fontsize=7, loc='left')
        
        overlay = np.zeros((*data[i].shape, 4))  # fully transparent by default
        overlay[np.isnan(data[i])] = [1, 0, 0, 1]  # red for nan
        overlay[np.isinf(data[i])] = [0, 0, 0, 1]  # black for +/- infinity
        overlay[data[i] == 0]          = [1, 1, 1, 1]  # white for zero
        ax.imshow(overlay)

        ax.axes.xaxis.set_ticks([])
        ax.axes.yaxis.set_ticks([])

    # Hide unused subplots
    for idx in range(n_plots, n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].set_visible(False)

    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6, format=mpl.ticker.StrMethodFormatter('$10^{{{x:.0f}}}$'))

    return fig

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <path1> [path2 ...]")
        sys.exit(1)

    paths = sys.argv[1:]

    # Paths whose directory is one of these should come first, in this order
    priority = ["A", "B", "C_ref", "C_ref_abs"]

    def sort_key(path):
        parts = Path(path).parts
        for i, name in enumerate(priority):
            if name in parts:
                return (0, i, path)
        return (1, 0, path)

    paths = sorted(paths, key=sort_key)

    # Within the "remaining" group, sort_key already orders alphabetically
    # by path because of the trailing `path` element in the tuple.

    fig = plot_matrix_heatmaps(paths)
    output_path = "results/figure.png"
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

if __name__ == "__main__":
    main()