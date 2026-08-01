#!/usr/bin/env python3
"""Visualise the convergence and benchmark CSVs written by option_pricer.

Reads results/convergence.csv and results/benchmark.csv (paths overridable)
and renders three figures:

  convergence.png    log-log |error| vs iterations per method, with O(1/N)
                     and O(1/sqrt(N)) slope guides - the convergence-rate
                     claims in the README, drawn instead of tabulated
  cost_accuracy.png  log-log |error| vs computation time per method: the
                     honest comparison, accuracy per unit of compute
  throughput.png     benchmark throughput per engine (log scale)

Dependencies: matplotlib only (stdlib csv otherwise). Usage:

  python3 scripts/plot_convergence.py [--results-dir results] [--out-dir results/plots]
"""

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # File output only; no display server required.
import matplotlib.pyplot as plt

# Errors at machine precision would drag a log axis to 1e-16; clip instead.
ERROR_FLOOR = 1e-12


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def group_by_method(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["method"]].append(row)
    return grouped


def plot_convergence(rows, out_path):
    grouped = group_by_method(rows)
    fig, ax = plt.subplots(figsize=(8, 6))

    for method, points in grouped.items():
        points.sort(key=lambda r: float(r["iterations"]))
        xs = [float(r["iterations"]) for r in points]
        ys = [max(abs(float(r["error"])), ERROR_FLOOR) for r in points]
        ax.loglog(xs, ys, marker="o", label=method)

    # Slope guides anchored at the first series' starting point.
    if grouped:
        first = min(grouped.values(), key=lambda pts: float(pts[0]["iterations"]))
        x0 = float(first[0]["iterations"])
        y0 = max(abs(float(first[0]["error"])), ERROR_FLOOR)
        xs = [x0 * (2**k) for k in range(12)]
        ax.loglog(xs, [y0 * x0 / x for x in xs], "k--", alpha=0.4, label="O(1/N)")
        ax.loglog(xs, [y0 * math.sqrt(x0 / x) for x in xs], "k:", alpha=0.4,
                  label="O(1/√N)")

    ax.set_xlabel("Iterations (steps or paths)")
    ax.set_ylabel("|price - reference|")
    ax.set_title("Convergence to the analytic reference")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_cost_accuracy(rows, out_path):
    grouped = group_by_method(rows)
    fig, ax = plt.subplots(figsize=(8, 6))

    for method, points in grouped.items():
        points.sort(key=lambda r: float(r["computation_time_ms"]))
        xs = [max(float(r["computation_time_ms"]), 1e-6) for r in points]
        ys = [max(abs(float(r["error"])), ERROR_FLOOR) for r in points]
        ax.loglog(xs, ys, marker="o", label=method)

    ax.set_xlabel("Computation time (ms)")
    ax.set_ylabel("|price - reference|")
    ax.set_title("Cost-accuracy frontier: error per unit of compute")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_throughput(rows, out_path):
    rows = sorted(rows, key=lambda r: float(r["throughput_per_s"]))
    methods = [r["method"] for r in rows]
    throughput = [float(r["throughput_per_s"]) for r in rows]

    fig, ax = plt.subplots(figsize=(8, 0.6 * len(rows) + 2))
    ax.barh(methods, throughput)
    ax.set_xscale("log")
    ax.set_xlabel("Valuations per second")
    ax.set_title("Benchmark throughput per engine")
    ax.grid(True, axis="x", which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results-dir", default="results", type=Path,
                        help="Directory holding convergence.csv and benchmark.csv")
    parser.add_argument("--out-dir", default=None, type=Path,
                        help="Output directory for PNGs (default: <results-dir>/plots)")
    args = parser.parse_args()

    out_dir = args.out_dir if args.out_dir is not None else args.results_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    rendered = []
    convergence_csv = args.results_dir / "convergence.csv"
    if convergence_csv.exists():
        rows = read_rows(convergence_csv)
        plot_convergence(rows, out_dir / "convergence.png")
        plot_cost_accuracy(rows, out_dir / "cost_accuracy.png")
        rendered += ["convergence.png", "cost_accuracy.png"]
    else:
        print(f"skip: {convergence_csv} not found (run option_pricer first)")

    benchmark_csv = args.results_dir / "benchmark.csv"
    if benchmark_csv.exists():
        plot_throughput(read_rows(benchmark_csv), out_dir / "throughput.png")
        rendered.append("throughput.png")
    else:
        print(f"skip: {benchmark_csv} not found (run option_pricer first)")

    if not rendered:
        return 1
    for name in rendered:
        print(f"wrote {out_dir / name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
