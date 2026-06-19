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

def plot_lines(grouped, ax, scale_style, y_label):
    """Draw one line per scale into ax."""
    impl = next(iter(grouped))
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
    all_ds_flat = sorted({pt[0] for points in grouped[impl].values() for pt in points})
    if all_ds_flat:
        ax.plot(all_ds_flat, all_ds_flat, color='black', linestyle=':', linewidth=1.5)
        ax.text(all_ds_flat[0], all_ds_flat[0], r"  $y=x$",
                ha='left', va='top', fontsize=FONTSIZE * 0.8, color='black')
    ax.set_xscale('log', base=2)
    ax.xaxis.set_major_formatter(mpl.ticker.LogFormatterMathtext(2.0))
    ax.set_yscale('log', base=10)
    ax.yaxis.set_major_formatter(mpl.ticker.LogFormatterMathtext(10.0))
    ax.tick_params(labelsize=FONTSIZE)
    ax.set_xlabel("dimension (m = n = k)", fontsize=FONTSIZE)
    ax.set_ylabel(y_label, fontsize=FONTSIZE)
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


def plot_heatmap(grids, all_ds, all_scales, ax, cmap, y_label):
    """Draw the heatmap for the single impl into ax."""
    impl = next(iter(grids))
    Z = grids[impl]
    finite = Z[np.isfinite(Z)]
    vmin, vmax = (finite.min(), finite.max()) if finite.size > 0 else (0.0, 1.0)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    im = ax.imshow(
        np.where(np.isfinite(Z), Z, np.nan),
        cmap=cmap,
        norm=norm,
        aspect='auto',
        origin='lower',
        interpolation='nearest',
    )
    labels = {}
    for r in range(Z.shape[0]):
        for c in range(Z.shape[1]):
            val = Z[r, c]
            if np.isposinf(val):
                labels[(r, c)] = "+INF"
            elif np.isneginf(val):
                labels[(r, c)] = "-INF"
            elif np.isfinite(val):
                labels[(r, c)] = f"{val:.2g}"
    max_len = max((len(s) for s in labels.values()), default=5)
    cell_fontsize = FONTSIZE * min(1.0, 5 / max_len)

    for (row, col), label in labels.items():
        val = Z[row, col]
        if not np.isfinite(val):
            ax.add_patch(mpl.patches.Rectangle(
                (col - 0.5, row - 0.5), 1, 1, color='red', zorder=2))
            text_color = "white"
        else:
            rgba = cmap(norm(val))
            luminance = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
            text_color = "white" if luminance < 0.5 else "black"
        ax.text(col, row, label, ha="center", va="center",
                fontsize=cell_fontsize, color=text_color, zorder=3)
    ax.set_xlabel("dimension (m = n = k)", fontsize=FONTSIZE)
    ax.set_ylabel("scale factor", fontsize=FONTSIZE)
    ax.set_xticks(range(len(all_ds)))
    ax.set_xticklabels(all_ds, fontsize=FONTSIZE, rotation=45, ha='right')
    ax.set_yticks(range(len(all_scales)))
    ax.set_yticklabels([rf"$\pm 2^{{{s}}}$" for s in all_scales], fontsize=FONTSIZE)
    cb = ax.figure.colorbar(im, ax=ax, shrink=0.75)
    cb.set_label(y_label, fontsize=FONTSIZE)
    cb.ax.tick_params(labelsize=FONTSIZE)
    cb.ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True))
    cb_ot = cb.ax.yaxis.get_offset_text()
    cb_ot.set_fontsize(FONTSIZE)
    cb_ot.set_ha('left')
    cb_ot.set_x(0)

def plot_max_vals(paths, y_label="f(n)", figsize_per_plot: tuple = (9, 6)) -> plt.Figure:
    if isinstance(paths, str):
        paths = [paths]

    grouped, value_descriptor = collect_max_vals(paths)
    impls = sorted(grouped.keys())
    if len(impls) != 1:
        raise ValueError(f"Expected exactly one impl in paths, got: {impls}")

    all_ds, all_scales, grids = build_heatmap_grid(grouped)
    heatmap_cmap = mpl.cm.get_cmap('cividis')

    mpl.rcParams['font.family'] = 'serif'
    global FONTSIZE
    FONTSIZE = 18
    line_w = figsize_per_plot[0]
    heat_w = max(6.0, len(all_ds) * 1.5)
    plot_h = max(6.0, len(all_scales) * 1.0)
    fig, axes = plt.subplots(
        1, 2,
        figsize=(line_w + heat_w, plot_h),
        gridspec_kw={'width_ratios': [line_w, heat_w]},
    )

    markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', '<', '>', 'h', 'p']
    all_scales_for_lines = sorted({s for impl_data in grouped.values() for s in impl_data})
    scale_style = {
        s: {"marker": markers[i % len(markers)]}
        for i, s in enumerate(all_scales_for_lines)
    }

    plot_lines(grouped, axes[0], scale_style, y_label)
    plot_heatmap(grids, all_ds, all_scales, axes[1], cmap=heatmap_cmap, y_label=y_label)

    fig.tight_layout(w_pad=6)
    return fig


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_type_dir> <impl>")
        print(f"  e.g. {sys.argv[0]} results/raw_data/test2c openblas")
        sys.exit(1)

    input_type_dir = Path(sys.argv[1])
    impl = sys.argv[2]
    input_type = input_type_dir.name

    y_label_map = {
        "test2c":                r"$f_A(n)$",
        "test3":                 r"$f_B(n)$",
        "test2c_w_inner_permute": "max relative error",
        "test2c_w_outer_permute": "max relative error",
    }
    y_label = y_label_map.get(input_type, r"$f(n)$")

    filename = "rel_errors.fp64" if "permute" in input_type else "f_of_n.fp64"

    WIDTH = 60
    prefix = f"━━━ Generating figure for {input_type} / {impl} "
    print("\n" + prefix + "━" * max(0, WIDTH - len(prefix)))

    paths = sorted(input_type_dir.glob(f"**/{impl}/{filename}"))
    if not paths:
        print(f"No f_of_n.fp64 files found for impl '{impl}' under {input_type_dir}")
        sys.exit(1)

    fig = plot_max_vals(paths, y_label=y_label)

    output_dir = Path("results/figures")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{input_type}_{impl}.png"
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

    sep = "━" * 60
    print(sep)
    print(f"  Figure written to:")
    print(f"  {output_path}")
    print(sep)


if __name__ == "__main__":
    main()
