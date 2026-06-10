import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import defaultdict
from pathlib import Path

from plot_utils import load_matrix, parse_path_metadata

def collect_max_vals(paths):
    """Group max values by impl and scale. Returns (grouped, value_descriptor) where
    grouped is dict[impl][s] -> list of (d, max)."""
    grouped = defaultdict(lambda: defaultdict(list))
    descriptors = set()
    for p in paths:
        meta = parse_path_metadata(p)
        matrix = load_matrix(p)
        max_val = np.nanmax(matrix)
        grouped[meta["impl"]][meta["s"]].append((meta["d"], max_val))
        descriptors.add(meta["value_descriptor"])
    if len(descriptors) > 1:
        raise ValueError(f"Expected a single value descriptor across all paths, got: {sorted(descriptors)}")
    value_descriptor = descriptors.pop() if descriptors else "value"
    return grouped, value_descriptor

def plot_lines(grouped, axes_row, scale_style, value_descriptor):
    """Draw one line plot per impl into the given row of axes."""
    impls = sorted(grouped.keys())
    for i, impl in enumerate(impls):
        ax = axes_row[i]
        for s in sorted(grouped[impl]):
            points = sorted(grouped[impl][s])
            ds = [pt[0] for pt in points]
            ys = [pt[1] for pt in points]
            style = scale_style[s]
            ax.plot(
                ds, ys,
                label=rf"$\pm 2^{{{s}}}$",
                marker=style["marker"],
                markersize=16,
                linewidth=1.5,
            )
        ax.tick_params(labelsize=FONTSIZE)
        ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True))
        ax.yaxis.get_offset_text().set_fontsize(FONTSIZE)
        ax.set_xlabel("dimension (m = n = k)", fontsize=FONTSIZE)
        ax.set_ylabel(r"$f_A(n)$", fontsize=FONTSIZE)
        ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.6)
        ax.legend(fontsize=FONTSIZE, loc='upper left', bbox_to_anchor=(1.01, 1),
                  borderaxespad=0., title="scale factor", title_fontsize=FONTSIZE)

def build_heatmap_grid(grouped):
    """
    Reshape the grouped data into a dict[impl] -> (ds, scales, Z) where
    Z[i, j] = max value at scale scales[i], dimension ds[j].
    Missing combinations become NaN.
    """
    all_ds = sorted({d for impl_data in grouped.values()
                       for points in impl_data.values()
                       for d, _ in points})
    all_scales = sorted({s for impl_data in grouped.values() for s in impl_data})

    grids = {}
    for impl, impl_data in grouped.items():
        Z = np.full((len(all_scales), len(all_ds)), np.nan)
        for s, points in impl_data.items():
            i = all_scales.index(s)
            for d, val in points:
                j = all_ds.index(d)
                Z[i, j] = val
        grids[impl] = Z
    return all_ds, all_scales, grids


def plot_heatmaps(grids, all_ds, all_scales, axes_row, cmap, value_descriptor):
    """Draw one heatmap per impl into the given row of axes."""
    impls = sorted(grids.keys())
    for i, impl in enumerate(impls):
        ax = axes_row[i]
        Z = grids[impl]
        finite = Z[np.isfinite(Z)]
        vmin, vmax = (finite.min(), finite.max()) if finite.size > 0 else (0.0, 1.0)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        im = ax.imshow(
            Z,
            cmap=cmap,
            norm=norm,
            aspect='auto',
            origin='lower',
            interpolation='nearest',
        )
        labels = {(r, c): f"{Z[r, c]:.2g}"
                  for r in range(Z.shape[0]) for c in range(Z.shape[1])
                  if np.isfinite(Z[r, c])}
        max_len = max((len(s) for s in labels.values()), default=5)
        cell_fontsize = FONTSIZE * min(1.0, 5 / max_len)

        for (row, col), label in labels.items():
            rgba = cmap(norm(Z[row, col]))
            luminance = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
            text_color = "white" if luminance < 0.5 else "black"
            ax.text(col, row, label, ha="center", va="center",
                    fontsize=cell_fontsize, color=text_color)
        ax.set_xlabel("dimension (m = n = k)", fontsize=FONTSIZE)
        ax.set_ylabel("scale factor", fontsize=FONTSIZE)
        ax.set_xticks(range(len(all_ds)))
        ax.set_xticklabels(all_ds, fontsize=FONTSIZE, rotation=45, ha='right')
        ax.set_yticks(range(len(all_scales)))
        ax.set_yticklabels([rf"$\pm 2^{{{s}}}$" for s in all_scales], fontsize=FONTSIZE)
        cb = ax.figure.colorbar(im, ax=ax, shrink=0.75)
        cb.set_label(r"$f_A(n)$", fontsize=FONTSIZE)
        cb.ax.tick_params(labelsize=FONTSIZE)
        cb.ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True))
        cb_ot = cb.ax.yaxis.get_offset_text()
        cb_ot.set_fontsize(FONTSIZE)
        cb_ot.set_ha('left')
        cb_ot.set_x(0)

def plot_max_vals(paths, figsize_per_plot: tuple = (9, 6)) -> plt.Figure:
    if isinstance(paths, str):
        paths = [paths]

    grouped, value_descriptor = collect_max_vals(paths)
    impls = sorted(grouped.keys())
    n_impls = len(impls)

    mpl.rcParams['font.family'] = 'serif'
    global FONTSIZE
    FONTSIZE = 18
    fig, axes = plt.subplots(
        n_impls, 2,
        figsize=(figsize_per_plot[0] * 2, figsize_per_plot[1] * n_impls),
        squeeze=False,
    )

    # Shared scale-factor styling for the line plots
    markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', '<', '>', 'h', 'p']
    all_scales_for_lines = sorted({s for impl_data in grouped.values() for s in impl_data})
    scale_style = {
        s: {"marker": markers[i % len(markers)]}
        for i, s in enumerate(all_scales_for_lines)
    }

    all_ds, all_scales, grids = build_heatmap_grid(grouped)
    heatmap_cmap = mpl.cm.get_cmap('cividis')

    line_axes = [axes[i][0] for i in range(n_impls)]
    heat_axes = [axes[i][1] for i in range(n_impls)]

    plot_lines(grouped, line_axes, scale_style, value_descriptor)
    plot_heatmaps(grids, all_ds, all_scales, heat_axes,
                  cmap=heatmap_cmap, value_descriptor=value_descriptor)

    fig.tight_layout(w_pad=6, h_pad=4)
    return fig


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <unique_dir>")
        sys.exit(1)

    unique_dir = Path(sys.argv[1])
    paths = sorted((unique_dir / "raw_data").glob("**/f_of_n.fp64"))
    if not paths:
        print(f"No f_of_n.fp64 files found under {unique_dir / 'raw_data'}")
        sys.exit(1)

    fig = plot_max_vals(paths)
    output_path = unique_dir / "figure.png"
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()