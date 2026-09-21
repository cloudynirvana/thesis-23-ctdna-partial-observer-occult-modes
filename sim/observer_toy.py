#!/usr/bin/env python3
"""Hybrid occult-mode toy and liquid-biopsy-style observation maps.

The mode path is an exogenous switch. The latent count is a scalar ODE in
each mode. Two families of maps are applied to that pair:

- a burden map, which reads the latent count and does not read turnover;
- a two-channel map, whose channels are turnover and net growth, scaled by
  the same latent count, then corrupted by additive Gaussian noise.

Synthetic outputs only. Not an assay. Seed 20260921.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SEED = 20260921
N0 = 500_000.0
KAPPA = 1.0e-3
T_END = 240.0
# (t0, t1, mode, birth, death). Switches are exogenous.
PIECES = (
    (0.0, 60.0, "Q", 0.001, 0.001),
    (60.0, 120.0, "A", 0.090, 0.090),
    (120.0, 126.0, "P", 0.120, 0.040),
    (126.0, 150.0, "I", 0.015, 0.075),
    (150.0, 240.0, "Q", 0.001, 0.001),
)
SWITCHES = (60.0, 120.0, 126.0, 150.0)
MODES = ("Q", "A", "P", "I")
MODE_COLOR = {
    "Q": "#0072B2",
    "A": "#E69F00",
    "P": "#D55E00",
    "I": "#009E73",
}

# Burden map: log-normal multiplier, median N. Logarithmic sd.
BURDEN_LOG_SD = 0.02
# Two-channel additive sd on the primary run. Sweep is separate.
SIGMA_PRIMARY = 3.0
SIGMA_SWEEP = (1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 32.0, 48.0, 64.0, 96.0)
# Amplitude gate between quiescent shedding and the next mode's shedding.
TAU = 6.0
HOLD = 2
MATCH_TOL = 4.0
N_REP = 400
DENSE_DT = 2.0
SPARSE_DT = 28.0
# Slope gates for the burden decoder, midway from 0 to the nonzero nets.
SLOPE_POS = 0.04
SLOPE_NEG = -0.03

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"


def piece_at(t: float) -> tuple:
    if t < 0 or t > T_END:
        raise ValueError(t)
    if t == T_END:
        return PIECES[-1]
    for piece in PIECES:
        if piece[0] <= t < piece[1]:
            return piece
    raise ValueError(t)


def mode_at(t: float) -> str:
    return piece_at(t)[2]


def latent_count(t: float) -> float:
    if t < 0 or t > T_END:
        raise ValueError(t)
    n = N0
    for t0, t1, _mode, birth, death in PIECES:
        if t <= t0:
            break
        dt = min(t, t1) - t0
        n *= math.exp((birth - death) * dt)
        if t <= t1:
            break
    return n


def output_mean(t: float) -> tuple[float, float, float, str]:
    """Return N, turnover channel, net-growth channel, mode. Noise-free."""
    _t0, _t1, mode, birth, death = piece_at(t)
    n = latent_count(t)
    y1 = KAPPA * (birth + death) * n
    y2 = KAPPA * (birth - death) * n
    return n, y1, y2, mode


def grid(dt: float, offset: float = 0.0) -> np.ndarray:
    if offset == 0.0:
        return np.arange(0.0, T_END + 1e-9, dt)
    times = np.arange(offset, T_END - 1e-9, dt)
    return times


def rho_targets() -> dict[str, float]:
    targets = {}
    for _t0, _t1, mode, birth, death in PIECES:
        turn = birth + death
        targets[mode] = (birth - death) / turn
    return targets


def decode_two_channel(y: np.ndarray, targets: dict[str, float]) -> list[str]:
    labels = []
    for y1, y2 in y:
        if y1 < TAU:
            labels.append("Q")
            continue
        rho = y2 / y1
        best = min(("P", "A", "I"), key=lambda m: abs(rho - targets[m]))
        labels.append(best)
    return labels


def decode_burden(y: np.ndarray, times: np.ndarray) -> list[str]:
    """Tie-break: a flat slope is labelled Q. A is never emitted."""
    labels = ["Q"]
    log_y = np.log(np.maximum(y, 1e-12))
    for i in range(1, len(times)):
        slope = (log_y[i] - log_y[i - 1]) / (times[i] - times[i - 1])
        if slope > SLOPE_POS:
            labels.append("P")
        elif slope < SLOPE_NEG:
            labels.append("I")
        else:
            labels.append("Q")
    return labels


def confirmed_switches(labels: list[str], times: np.ndarray, hold: int) -> tuple[list[float], list[str]]:
    current = labels[0]
    confirmed = [current]
    run_label = labels[0]
    run_start = 0
    run_len = 1
    events: list[float] = []
    for i in range(1, len(labels)):
        if labels[i] == run_label:
            run_len += 1
        else:
            run_label = labels[i]
            run_start = i
            run_len = 1
        new_current = run_label if run_len >= hold else current
        if new_current != current:
            events.append(float(times[run_start]))
            current = new_current
        confirmed.append(current)
    return events, confirmed


def match_switches(estimated: list[float], truth: tuple[float, ...] = SWITCHES, tol: float = MATCH_TOL):
    used: set[int] = set()
    errors: list[float | None] = []
    for tr in truth:
        best_j = None
        best_d = None
        for j, est in enumerate(estimated):
            if j in used:
                continue
            dist = abs(est - tr)
            if dist <= tol and (best_d is None or dist < best_d):
                best_j = j
                best_d = dist
        if best_j is None:
            errors.append(None)
        else:
            used.add(best_j)
            errors.append(float(best_d))
    false_n = len(estimated) - len(used)
    return errors, false_n


def summarise_matches(rows: list[tuple[list[float | None], int]]) -> dict:
    n = len(rows)
    per = []
    for s in range(len(SWITCHES)):
        hits = [row[0][s] for row in rows if row[0][s] is not None]
        per.append(
            {
                "time": SWITCHES[s],
                "recall": len(hits) / n,
                "median_abs_error": float(np.median(hits)) if hits else None,
            }
        )
    all_four = sum(1 for errs, _f in rows if all(e is not None for e in errs)) / n
    matched = [e for errs, _f in rows for e in errs if e is not None]
    false_mean = float(np.mean([f for _e, f in rows]))
    return {
        "n": n,
        "per_switch": per,
        "recall_mean": float(np.mean([p["recall"] for p in per])),
        "fraction_all_four": all_four,
        "median_abs_error_matched": float(np.median(matched)) if matched else None,
        "false_switches_per_rep": false_mean,
    }


def confusion(truth: list[str], decoded: list[str]) -> np.ndarray:
    index = {m: i for i, m in enumerate(MODES)}
    mat = np.zeros((len(MODES), len(MODES)), dtype=float)
    for a, b in zip(truth, decoded):
        mat[index[a], index[b]] += 1
    return mat


def row_normalise(mat: np.ndarray) -> np.ndarray:
    out = mat.copy()
    for i in range(out.shape[0]):
        s = out[i].sum()
        if s > 0:
            out[i] /= s
    return out


def shade(ax, ymin_mode: bool = False) -> None:
    for t0, t1, mode, _b, _d in PIECES:
        ax.axvspan(t0, t1, color=MODE_COLOR[mode], alpha=0.16, lw=0, zorder=0)


def style_ax(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8.5)
    ax.set_xlim(0, T_END)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 8,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
        }
    )

    targets = rho_targets()
    dense_t = grid(DENSE_DT)
    sparse_t = grid(SPARSE_DT)
    offset_t = grid(DENSE_DT, offset=1.0)

    dense_pack = [output_mean(float(t)) for t in dense_t]
    n_dense = np.array([p[0] for p in dense_pack])
    y1_dense = np.array([p[1] for p in dense_pack])
    y2_dense = np.array([p[2] for p in dense_pack])
    mode_dense = [p[3] for p in dense_pack]
    n_sparse = np.array([latent_count(float(t)) for t in sparse_t])
    mode_sparse = [mode_at(float(t)) for t in sparse_t]

    # Noise-free signals used to place TAU.
    q_mask = np.array([m == "Q" for m in mode_dense])
    a_mask = np.array([m == "A" for m in mode_dense])
    non_q = ~q_mask
    max_q = float(y1_dense[q_mask].max())
    min_non_q = float(y1_dense[non_q].min())
    if not (max_q < TAU < min_non_q):
        raise SystemExit(f"TAU={TAU} is not between Q shedding {max_q} and other modes {min_non_q}")

    short = [t for t in sparse_t if 120.0 <= t < 126.0]
    if short:
        raise SystemExit(f"sparse grid hits the 6-day proliferative dwell: {short}")
    if abs(latent_count(30.0) - latent_count(90.0)) > 1e-6:
        raise SystemExit("Q and A plateaus do not share a latent count")
    if abs(latent_count(30.0) - N0) > 1e-6:
        raise SystemExit("plateau drifted from N0")

    # Interior windows for the equal-law check. Same latent count.
    q_win = np.array([8.0 <= t <= 52.0 for t in dense_t])
    a_win = np.array([68.0 <= t <= 112.0 for t in dense_t])
    if abs(n_dense[q_win].mean() - n_dense[a_win].mean()) > 1e-6:
        raise SystemExit("window means of N differ")

    noise_free_labels = decode_two_channel(np.column_stack([y1_dense, y2_dense]), targets)
    if noise_free_labels != mode_dense:
        bad = [(float(t), a, b) for t, a, b in zip(dense_t, mode_dense, noise_free_labels) if a != b]
        raise SystemExit(f"noise-free two-channel decoder missed {bad[:8]}")
    nf_events, _nf_conf = confirmed_switches(noise_free_labels, dense_t, HOLD)
    nf_err, nf_false = match_switches(nf_events)
    if any(e != 0.0 for e in nf_err) or nf_false:
        raise SystemExit(f"noise-free aligned recovery failed: {nf_err} false={nf_false} events={nf_events}")

    # Offset grid, noise-free: switches sit between samples.
    off_y = np.array([[output_mean(float(t))[1], output_mean(float(t))[2]] for t in offset_t])
    off_labels = decode_two_channel(off_y, targets)
    off_events, _ = confirmed_switches(off_labels, offset_t, HOLD)
    off_err, off_false = match_switches(off_events)
    if any(e is None or e > 1.0 + 1e-9 for e in off_err):
        raise SystemExit(f"offset noise-free recovery failed: {off_err} events={off_events}")

    rng = np.random.default_rng(SEED)

    def burden_stats(times: np.ndarray, n_rep: int) -> tuple[dict, np.ndarray, list[str], np.ndarray]:
        rows = []
        mat = np.zeros((4, 4))
        example_y = None
        example_dec = None
        # Equal-law z-tests on dense windows only when the grid is the dense one.
        z_reject = []
        mean_gaps = []
        for rep in range(n_rep):
            z = rng.normal(size=len(times))
            y = np.array([latent_count(float(t)) for t in times]) * np.exp(BURDEN_LOG_SD * z)
            labels = decode_burden(y, times)
            events, confirmed = confirmed_switches(labels, times, HOLD)
            rows.append(match_switches(events))
            truth = [mode_at(float(t)) for t in times]
            mat += confusion(truth, labels)
            if rep == 0:
                example_y = y
                example_dec = labels
            if len(times) == len(dense_t) and np.allclose(times, dense_t):
                a = y[q_win]
                b = y[a_win]
                # Normal approximation to a two-sample test of equal means.
                se = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
                zstat = abs(a.mean() - b.mean()) / se
                z_reject.append(zstat > 1.96)
                mean_gaps.append(abs(a.mean() - b.mean()) / N0)
        summary = summarise_matches(rows)
        summary["confusion_counts"] = mat.tolist()
        summary["confusion_row_fraction"] = row_normalise(mat).tolist()
        if mean_gaps:
            summary["equal_law"] = {
                "window_Q": [8.0, 52.0],
                "window_A": [68.0, 112.0],
                "n_per_window": int(q_win.sum()),
                "median_abs_mean_gap_over_N0": float(np.median(mean_gaps)),
                "p95_abs_mean_gap_over_N0": float(np.quantile(mean_gaps, 0.95)),
                "normal_rejection_rate": float(np.mean(z_reject)),
            }
        return summary, example_y, example_dec, mat

    burden_dense, ex_b_y, ex_b_dec, mat_b = burden_stats(dense_t, N_REP)
    burden_sparse, ex_s_y, ex_s_dec, _mat_bs = burden_stats(sparse_t, N_REP)

    # Two-channel Monte Carlo, including a sweep over sigma.
    sweep_out = []
    example = None
    mat_primary = np.zeros((4, 4))
    offset_rows = []
    p_called = []
    for sigma in SIGMA_SWEEP:
        rows = []
        mat = np.zeros((4, 4))
        for rep in range(N_REP):
            noise = rng.normal(size=(len(dense_t), 2))
            y = np.column_stack([y1_dense, y2_dense]) + sigma * noise
            labels = decode_two_channel(y, targets)
            events, confirmed = confirmed_switches(labels, dense_t, HOLD)
            rows.append(match_switches(events))
            mat += confusion(mode_dense, labels)
            if sigma == SIGMA_PRIMARY and rep == 0:
                example = {
                    "y": y,
                    "raw_labels": labels,
                    "confirmed": confirmed,
                    "events": events,
                }
            if sigma == SIGMA_PRIMARY:
                # Sparse two-channel on its own draws, same sigma, counted once.
                pass
        summary = summarise_matches(rows)
        summary["sigma"] = sigma
        summary["confusion_counts"] = mat.tolist()
        summary["confusion_row_fraction"] = row_normalise(mat).tolist()
        sweep_out.append(summary)
        if sigma == SIGMA_PRIMARY:
            mat_primary = mat

    # Sparse two-channel and offset dense, at the primary sigma.
    sparse_channel_rows = []
    sparse_p_hits = 0
    for _rep in range(N_REP):
        noise = rng.normal(size=(len(sparse_t), 2))
        means = np.array([[output_mean(float(t))[1], output_mean(float(t))[2]] for t in sparse_t])
        y = means + SIGMA_PRIMARY * noise
        labels = decode_two_channel(y, targets)
        # Hold of 1: a single sparse sample is already the observation.
        events, confirmed = confirmed_switches(labels, sparse_t, hold=1)
        # Place the event at the midpoint between the samples that differ.
        # confirmed_switches with hold 1 already uses the new sample's time.
        # Recompute midpoints explicitly so the clock is not given the later stamp.
        mid_events = []
        for i in range(1, len(labels)):
            if labels[i] != labels[i - 1]:
                mid_events.append(float(0.5 * (sparse_t[i] + sparse_t[i - 1])))
        sparse_channel_rows.append(match_switches(mid_events))
        if "P" in labels:
            sparse_p_hits += 1
    sparse_channel = summarise_matches(sparse_channel_rows)
    sparse_channel["fraction_reps_emitting_P"] = sparse_p_hits / N_REP
    sparse_channel["sigma"] = SIGMA_PRIMARY

    for _rep in range(N_REP):
        noise = rng.normal(size=(len(offset_t), 2))
        means = np.array([[output_mean(float(t))[1], output_mean(float(t))[2]] for t in offset_t])
        y = means + SIGMA_PRIMARY * noise
        labels = decode_two_channel(y, targets)
        events, _confirmed = confirmed_switches(labels, offset_t, HOLD)
        offset_rows.append(match_switches(events))
    offset_summary = summarise_matches(offset_rows)
    offset_summary["noise_free_errors"] = off_err
    offset_summary["noise_free_false"] = off_false
    offset_summary["sigma"] = SIGMA_PRIMARY

    primary = next(s for s in sweep_out if s["sigma"] == SIGMA_PRIMARY)

    # Distances at representative times, noise-free.
    landmarks = {}
    for name, t in (("mid_Q", 30.0), ("mid_A", 90.0), ("P_entry", 120.0), ("I_entry", 126.0), ("I_exit", 150.0), ("late_Q", 200.0)):
        n, y1, y2, mode = output_mean(t)
        landmarks[name] = {
            "t": t,
            "mode": mode,
            "N": n,
            "y1": y1,
            "y2": y2,
            "rho": None if y1 == 0 else y2 / y1,
        }

    # Claims that the manuscript is allowed to make. Fail closed.
    el = burden_dense["equal_law"]
    if not (0.01 < el["normal_rejection_rate"] < 0.12):
        raise SystemExit(f"burden equal-law rejection rate out of null range: {el}")
    if el["median_abs_mean_gap_over_N0"] > 0.01:
        raise SystemExit(f"burden windows separated: {el}")
    qa_switch = burden_dense["per_switch"][0]
    if qa_switch["recall"] > 0.05:
        raise SystemExit(f"burden map recovered Q→A: {qa_switch}")
    if primary["fraction_all_four"] < 0.95:
        raise SystemExit(f"dense two-channel did not recover switches at sigma={SIGMA_PRIMARY}: {primary['fraction_all_four']}")
    if primary["median_abs_error_matched"] > 2.0:
        raise SystemExit(f"timing error too large: {primary['median_abs_error_matched']}")
    if "P" in mode_sparse:
        raise SystemExit(f"sparse grid landed inside P: {list(zip(sparse_t.tolist(), mode_sparse))}")
    # Recovery: every true switch matched, and false alarms rare.
    # Failure: either the matches drop, or false alarms overwhelm the four true times.
    recovering = [
        s for s in sweep_out
        if s["fraction_all_four"] >= 0.90 and s["false_switches_per_rep"] <= 0.25
    ]
    failing = [
        s for s in sweep_out
        if s["fraction_all_four"] < 0.50 or s["false_switches_per_rep"] >= 2.0
    ]
    if not recovering or not failing:
        brief = [
            (s["sigma"], round(s["fraction_all_four"], 3), round(s["false_switches_per_rep"], 3), round(s["recall_mean"], 3))
            for s in sweep_out
        ]
        raise SystemExit(f"sweep did not show both recovery and failure: {brief}")

    results = {
        "seed": SEED,
        "N0": N0,
        "kappa": KAPPA,
        "T_end": T_END,
        "pieces": [
            {"t0": a, "t1": b, "mode": m, "birth": birth, "death": death, "net": birth - death, "turnover": birth + death}
            for a, b, m, birth, death in PIECES
        ],
        "switches": list(SWITCHES),
        "tau": TAU,
        "hold": HOLD,
        "match_tol_days": MATCH_TOL,
        "burden_log_sd": BURDEN_LOG_SD,
        "slope_pos": SLOPE_POS,
        "slope_neg": SLOPE_NEG,
        "rho_targets": targets,
        "n_replicates": N_REP,
        "dense_dt": DENSE_DT,
        "sparse_dt": SPARSE_DT,
        "sparse_times": sparse_t.tolist(),
        "max_Q_shedding": max_q,
        "min_non_Q_shedding": min_non_q,
        "landmarks": landmarks,
        "noise_free_two_channel": {
            "events": nf_events,
            "errors": nf_err,
            "false": nf_false,
            "label_mismatches": 0,
        },
        "burden_dense": {k: v for k, v in burden_dense.items() if k not in ()},
        "burden_sparse": burden_sparse,
        "two_channel_sweep": [
            {k: v for k, v in s.items() if k not in ("confusion_counts",)}
            for s in sweep_out
        ],
        "two_channel_primary_confusion_counts": mat_primary.tolist(),
        "two_channel_primary_confusion_row_fraction": row_normalise(mat_primary).tolist(),
        "two_channel_sparse": sparse_channel,
        "two_channel_offset": offset_summary,
        "modes_order": list(MODES),
    }
    # Keep the full primary confusion; sweep entries already drop raw counts except primary stored above.
    # Restore primary counts inside the sweep entry for the sigma that is cited.
    for s, full in zip(results["two_channel_sweep"], sweep_out):
        if full["sigma"] == SIGMA_PRIMARY:
            s["confusion_counts"] = full["confusion_counts"]

    def clean(obj):
        if isinstance(obj, dict):
            return {str(k): clean(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [clean(v) for v in obj]
        if isinstance(obj, (np.floating, float)):
            value = float(obj)
            if math.isnan(value):
                return None
            return round(value, 6)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return clean(obj.tolist())
        return obj

    payload = clean(results)
    (ROOT / "results.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    # ----- figures -----
    fine_t = np.linspace(0.0, T_END, 961)
    fine_n = np.array([latent_count(float(t)) for t in fine_t])
    fine_y1 = np.array([output_mean(float(t))[1] for t in fine_t])
    fine_y2 = np.array([output_mean(float(t))[2] for t in fine_t])

    fig, axes = plt.subplots(2, 1, figsize=(7.15, 5.4), sharex=True, gridspec_kw={"height_ratios": [1.15, 1]})
    shade(axes[0])
    axes[0].plot(fine_t, fine_n, color="#222", lw=1.4)
    axes[0].set_ylabel("Latent count N")
    axes[0].set_yscale("log")
    axes[0].set_title("Latent count under an exogenous mode path")
    style_ax(axes[0])
    shade(axes[1])
    axes[1].plot(fine_t, fine_y1, color="#222", lw=1.3, label="Turnover channel")
    axes[1].plot(fine_t, fine_y2, color="#333", lw=1.1, ls="--", label="Net-growth channel")
    axes[1].axhline(TAU, color="#666", lw=0.6, ls=":")
    axes[1].set_ylabel("Noise-free output")
    axes[1].set_xlabel("Toy time (day)")
    axes[1].legend(frameon=False, loc="upper right")
    style_ax(axes[1])
    fig.tight_layout()
    fig.savefig(FIG / "latent_and_noise_free_outputs.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(7.15, 5.6), sharex=True)
    shade(axes[0])
    axes[0].plot(fine_t, fine_n, color="#222", lw=1.2, label="N(t)")
    axes[0].scatter(sparse_t, ex_s_y, s=28, color="#000", zorder=3, label="Sparse samples")
    axes[0].plot(dense_t, ex_b_y, color="#0072B2", lw=0.8, alpha=0.9, label="Dense replicate")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Constant-coefficient output")
    axes[0].legend(frameon=False, loc="lower left")
    style_ax(axes[0])
    shade(axes[1])
    y_ex = example["y"]
    axes[1].plot(fine_t, fine_y1, color="#222", lw=1.0, label="Turnover, noise-free")
    axes[1].plot(dense_t, y_ex[:, 0], color="#D55E00", lw=0.9, label="Turnover, one noisy replicate")
    axes[1].plot(dense_t, y_ex[:, 1], color="#0072B2", lw=0.8, alpha=0.85, label="Net-growth, same replicate")
    for sw in SWITCHES:
        axes[1].axvline(sw, color="#444", lw=0.5, ls="--")
    axes[1].set_ylabel("Two-channel output")
    axes[1].set_xlabel("Toy time (day)")
    axes[1].legend(frameon=False, loc="upper right")
    style_ax(axes[1])
    fig.tight_layout()
    fig.savefig(FIG / "confusing_map_versus_two_channel.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(7.15, 4.8), sharex=True)
    for ax, series, title in (
        (axes[0], ex_b_dec, "Constant-coefficient labels, before confirmation"),
        (axes[1], example["raw_labels"], "Dense two-channel labels, before confirmation"),
    ):
        shade(ax)
        code = {"Q": 0, "A": 1, "P": 2, "I": 3}
        truth_code = [code[m] for m in mode_dense]
        dec_code = [code[m] for m in series]
        ax.step(dense_t, truth_code, where="post", color="#222", lw=1.3, label="True mode")
        ax.step(dense_t, np.array(dec_code) - 0.08, where="post", color="#D55E00", lw=1.0, label="Decoded")
        ax.set_yticks([0, 1, 2, 3], ["Q", "A", "P", "I"])
        ax.set_ylabel("Mode")
        ax.set_title(title)
        ax.legend(frameon=False, loc="upper right")
        style_ax(ax)
    axes[1].set_xlabel("Toy time (day)")
    fig.tight_layout()
    fig.savefig(FIG / "decoded_modes.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.15, 3.6), layout="constrained")
    for ax, mat, title in (
        (axes[0], mat_b, "Constant-coefficient map"),
        (axes[1], mat_primary, "Dense two-channel map"),
    ):
        frac = row_normalise(mat)
        im = ax.imshow(frac, vmin=0, vmax=1, cmap="Blues")
        ax.set_xticks(range(4), MODES)
        ax.set_yticks(range(4), MODES)
        ax.set_xlabel("Decoded")
        ax.set_ylabel("True mode")
        ax.set_title(title)
        for i in range(4):
            for j in range(4):
                val = frac[i, j]
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color="#111" if val < 0.7 else "#fff")
    fig.colorbar(im, ax=axes, fraction=0.046, pad=0.04, label="Row fraction")
    fig.savefig(FIG / "confusion_matrices.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(7.15, 5.3), sharex=True)
    sigmas = [s["sigma"] for s in sweep_out]
    ax = axes[0]
    for idx, sw in enumerate(SWITCHES):
        recalls = [s["per_switch"][idx]["recall"] for s in sweep_out]
        ax.plot(sigmas, recalls, marker="o", ms=3.5, lw=1.2, label=f"Day {int(sw)}")
    ax.plot(
        sigmas,
        [s["fraction_all_four"] for s in sweep_out],
        color="#222",
        lw=1.6,
        marker="s",
        ms=3.5,
        label="All four",
    )
    ax.axhline(qa_switch["recall"], color="#0072B2", lw=1.0, ls="--", label="Constant-coefficient, Q to A")
    ax.set_ylabel("Recall within 4 days")
    ax.set_ylim(-0.02, 1.05)
    ax.legend(frameon=False, ncol=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax = axes[1]
    ax.plot(sigmas, [s["false_switches_per_rep"] for s in sweep_out], color="#222", marker="o", ms=3.5, lw=1.3)
    ax.set_xlabel("Additive noise sd on each channel")
    ax.set_ylabel("False switches per replicate")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "recall_versus_noise.png", dpi=160)
    plt.close(fig)

    print(json.dumps({
        "equal_law": payload["burden_dense"]["equal_law"],
        "burden_dense_switches": payload["burden_dense"]["per_switch"],
        "burden_dense_false": payload["burden_dense"]["false_switches_per_rep"],
        "primary_fraction_all_four": primary["fraction_all_four"],
        "primary": {
            "fraction_all_four": primary["fraction_all_four"],
            "recall_mean": primary["recall_mean"],
            "median_abs_error_matched": primary["median_abs_error_matched"],
            "false": primary["false_switches_per_rep"],
            "per_switch": primary["per_switch"],
        },
        "sparse_channel_P": sparse_channel["fraction_reps_emitting_P"],
        "sparse_channel_recall": sparse_channel["per_switch"],
        "offset": {
            "noise_free": off_err,
            "fraction_all_four": offset_summary["fraction_all_four"],
            "median": offset_summary["median_abs_error_matched"],
            "per_switch": offset_summary["per_switch"],
        },
        "sweep": [
            {
                "sigma": s["sigma"],
                "recall_mean": s["recall_mean"],
                "all_four": s["fraction_all_four"],
                "false": s["false_switches_per_rep"],
                "median": s["median_abs_error_matched"],
            }
            for s in sweep_out
        ],
        "landmarks": payload["landmarks"],
        "max_Q": max_q,
        "min_non_Q": min_non_q,
        "burden_row_fraction": payload["burden_dense"]["confusion_row_fraction"],
        "channel_row_fraction": payload["two_channel_primary_confusion_row_fraction"],
        "burden_sparse_switches": payload["burden_sparse"]["per_switch"],
    }, indent=2))


if __name__ == "__main__":
    main()
