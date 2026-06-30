#!/usr/bin/env python3
"""
Plot gffcompare evaluation results produced by eval_training.nf.

Two modes:
  --mode val   : transcript-level F1 over epochs, one line per validation
                 species (single chart with legend).
  --mode test  : transcript-level precision vs sensitivity scatter, one dot
                 per test species, labelled.

Two ways to call it:

  1. Auto-derive everything from the same YAML used by the pipeline:
       plot_eval.py --config_yaml config.yaml --run train1 --mode val
       plot_eval.py --config_yaml config.yaml --run train1 --mode test

  2. Pass paths explicitly (used by the Nextflow PLOT_EVAL process):
       plot_eval.py --stats_dir <dir> --species sp1,sp2 --mode val \\
                    --out my_plot.png --title train1
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXON_RE = re.compile(r"Exon level:\s*([\d.]+)\s*\|\s*([\d.]+)")
TX_RE = re.compile(r"Transcript level:\s*([\d.]+)\s*\|\s*([\d.]+)")
EPOCH_FROM_NAME_RE = re.compile(r"(?P<idx>[^.]+)\.epoch_(?P<epoch>\d+)\.")


def parse_stats(path):
    text = Path(path).read_text()
    m_tx = TX_RE.search(text)
    m_ex = EXON_RE.search(text)
    if not m_tx or not m_ex:
        raise RuntimeError(f"Could not parse Tx/Exon lines in {path}")
    return {
        "Transcript": {
            "sensitivity": float(m_tx.group(1)),
            "precision": float(m_tx.group(2)),
        },
        "Exon": {
            "sensitivity": float(m_ex.group(1)),
            "precision": float(m_ex.group(2)),
        },
    }


def f1(sens, prec):
    return 2 * sens * prec / (sens + prec) if (sens + prec) > 0 else 0.0


def plot_val(stats_dir, species, title, out_path):
    fig, ax = plt.subplots(figsize=(10, 6))
    plotted_any = False
    for sp in species:
        files = sorted(Path(stats_dir).glob(f"*.epoch_*.{sp}.stats"))
        points = []
        for f in files:
            m = EPOCH_FROM_NAME_RE.match(f.name)
            if not m:
                continue
            epoch = int(m.group("epoch"))
            tx = parse_stats(f)["Transcript"]
            points.append((epoch, f1(tx["sensitivity"], tx["precision"])))
        if not points:
            print(f"warn: no validation stats for {sp} in {stats_dir}", file=sys.stderr)
            continue
        points.sort()
        epochs = [p[0] for p in points]
        f1s = [p[1] for p in points]
        ax.plot(epochs, f1s, marker="o", label=sp)
        plotted_any = True

    if not plotted_any:
        raise RuntimeError(f"no stats matched any species in {stats_dir}")

    ax.set_xlabel("Epoch")
    ax.set_ylabel("Transcript F1 [%]")
    ax.set_title(f"Validation accuracy" + (f" -- {title}" if title else ""))
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    print(f"wrote {out_path}")


def plot_test(stats_dir, species, title, out_path):
    fig, ax = plt.subplots(figsize=(8, 8))
    plotted_any = False
    for sp in species:
        # Test-mode filenames look like "test_<i>.epoch_<n>.<sp>.stats".
        # If multiple test_weights are listed, plot the first.
        files = sorted(Path(stats_dir).glob(f"test_*.{sp}.stats"))
        if not files:
            print(f"warn: no test stats for {sp} in {stats_dir}", file=sys.stderr)
            continue
        tx = parse_stats(files[0])["Transcript"]
        ax.scatter(tx["precision"], tx["sensitivity"], s=140, marker="x", linewidth=3)
        ax.annotate(
            sp,
            (tx["precision"], tx["sensitivity"]),
            textcoords="offset points",
            xytext=(8, -4),
            fontsize=11,
        )
        plotted_any = True

    if not plotted_any:
        raise RuntimeError(f"no test stats matched any species in {stats_dir}")

    ax.set_xlabel("Precision [%]")
    ax.set_ylabel("Sensitivity [%]")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_title(f"Test accuracy (transcript level)" + (f" -- {title}" if title else ""))
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    print(f"wrote {out_path}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--config_yaml", help="Same YAML the pipeline uses.")
    p.add_argument("--run", help="Training run name from config.yaml (e.g. train1).")
    p.add_argument("--stats_dir", help="Directory of .stats files. Overrides config.")
    p.add_argument("--species", help="Comma-separated species list. Overrides config.")
    p.add_argument("--mode", choices=["val", "test"], required=True)
    p.add_argument("--out", help="Output PNG. Defaults to <stats_dir>/<run>_<mode>_acc.png.")
    p.add_argument("--title", default="", help="Plot title suffix (e.g. run name).")
    args = p.parse_args()

    # Auto-derive missing args from config when possible.
    if args.config_yaml:
        import yaml

        cfg = yaml.safe_load(Path(args.config_yaml).read_text())
        if not args.run:
            p.error("--run is required when --config_yaml is given")
        run_cfg = cfg.get("training", {}).get(args.run)
        if not run_cfg:
            p.error(f"training.{args.run} not found in {args.config_yaml}")
        if not args.stats_dir:
            args.stats_dir = run_cfg["eval_dir"]
        if not args.species:
            args.species = ",".join(cfg["species_split"][args.mode])
        if not args.title:
            args.title = args.run

    missing = [k for k in ("stats_dir", "species") if not getattr(args, k)]
    if missing:
        p.error(f"missing: {', '.join('--' + m for m in missing)} "
                "(or pass --config_yaml --run)")

    if not args.out:
        run_part = args.run if args.run else "run"
        args.out = str(Path(args.stats_dir) / f"{run_part}_{args.mode}_acc.png")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    species_list = [s.strip() for s in args.species.split(",") if s.strip()]

    if args.mode == "val":
        plot_val(args.stats_dir, species_list, args.title, args.out)
    else:
        plot_test(args.stats_dir, species_list, args.title, args.out)


if __name__ == "__main__":
    main()
