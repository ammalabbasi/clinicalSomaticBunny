#!/usr/bin/env python3
"""Sample-identity checks on `somalier relate` output (Clinical SomaticBunny).

Somalier sample names are <patient>_<sample>_<status> (set by SOMALIER_EXTRACT), so every digest maps back to one
sample-sheet row. Flags written to --flags (TSV, header always present, one row per problem):

  swap               two samples of the SAME patient (tumour-normal, or two tumours) with relatedness < --swap-max
  duplicate          samples of DIFFERENT patients with relatedness >= --dup-min (same individual twice, or a mislabel)
  insufficient_sites two samples of the same patient compared at < --min-sites sites: identity not confirmed
  sex_mismatch       inferred sex differs from the sheet's optional `sex` column (see infer_sex for tumour rules)

Tumour copy-number changes and LOH lower tumour-normal relatedness (a heavily rearranged tumour vs its own normal can
sit near 0.6-0.8, against ~1.0 for germline duplicates), so the swap threshold is loose. Pairs compared at fewer than
--min-sites sites are never called duplicates.
Uses only the standard library.
"""
import argparse
import csv
import sys

FLAG_COLS = ["flag", "sample_a", "sample_b", "patient_a", "patient_b", "relatedness", "n_sites", "threshold", "detail"]
SEX_COLS = ["sample", "patient", "status", "declared_sex", "inferred_sex", "x_call", "y_call",
            "x_depth_ratio", "y_depth_ratio", "x_het_fraction", "x_informative_sites"]


def norm_sex(v):
    v = (v or "").strip().lower()
    if v in ("male", "m", "xy", "1"):
        return "male"
    if v in ("female", "f", "xx", "2"):
        return "female"
    return "unknown"


def fnum(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return float("nan")
    return x


def read_sheet(path):
    rows = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if not r.get("patient"):
                continue
            sid = "{}_{}_{}".format(r["patient"], r["sample"], r["status"])
            rows[sid] = {"patient": r["patient"], "status": r["status"], "sex": norm_sex(r.get("sex"))}
    return rows


def read_tsv(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if not lines:
        return []
    header = lines[0].lstrip("#").split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


def infer_sex(s, run_has_y, min_x):
    gt = fnum(s.get("gt_depth_mean"))
    xd = fnum(s.get("X_depth_mean"))
    yd = fnum(s.get("Y_depth_mean"))
    het = fnum(s.get("X_het"))
    homalt = fnum(s.get("X_hom_alt"))
    x_ratio = xd / gt if gt > 0 else float("nan")
    y_ratio = yd / gt if gt > 0 else float("nan")
    informative = (het if het == het else 0) + (homalt if homalt == homalt else 0)
    het_frac = het / informative if informative > 0 else float("nan")

    x_call = "unknown"
    if informative >= min_x:
        if het_frac >= 0.2:
            x_call = "female"
        elif het_frac <= 0.05 and (x_ratio != x_ratio or x_ratio <= 0.75):
            x_call = "male"
    y_call = "unknown"
    if run_has_y and y_ratio == y_ratio:
        if y_ratio >= 0.2:
            y_call = "male"
        elif y_ratio <= 0.05:
            y_call = "female"
    calls = {c for c in (x_call, y_call) if c != "unknown"}
    inferred = calls.pop() if len(calls) == 1 else ("ambiguous" if len(calls) > 1 else "unknown")
    return inferred, x_call, y_call, x_ratio, y_ratio, het_frac, int(informative)


def fmt(x):
    return "NA" if x != x else "{:.3f}".format(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--swap-max", type=float, default=0.6)
    ap.add_argument("--dup-min", type=float, default=0.7)
    ap.add_argument("--min-sites", type=int, default=50)
    ap.add_argument("--min-x-sites", type=int, default=10)
    ap.add_argument("--flags", required=True)
    ap.add_argument("--sex-out", required=True)
    a = ap.parse_args()

    sheet = read_sheet(a.sheet)
    pairs = read_tsv(a.pairs)
    samples = read_tsv(a.samples)
    flags = []

    def meta(sid):
        return sheet.get(sid, {"patient": "NA", "status": "NA", "sex": "unknown"})

    unknown = sorted({p[k] for p in pairs for k in ("sample_a", "sample_b") if p[k] not in sheet}
                     | {s["sample_id"] for s in samples if s["sample_id"] not in sheet})
    if unknown:
        print("WARNING: somalier samples not in the sheet (not checked): " + ", ".join(unknown), file=sys.stderr)

    n_same = n_diff = n_low = 0
    for p in pairs:
        a_id, b_id = p["sample_a"], p["sample_b"]
        if a_id not in sheet or b_id not in sheet:
            continue
        pa, pb = sheet[a_id]["patient"], sheet[b_id]["patient"]
        rel = fnum(p.get("relatedness"))
        n_raw = fnum(p.get("n"))
        n = int(n_raw) if n_raw == n_raw else 0
        row = {"sample_a": a_id, "sample_b": b_id, "patient_a": pa, "patient_b": pb,
               "relatedness": fmt(rel), "n_sites": n}
        if pa == pb:
            n_same += 1
            if n < a.min_sites:
                flags.append(dict(row, flag="insufficient_sites", threshold=a.min_sites,
                                  detail="same patient compared at {} sites (< {}): identity not confirmed".format(n, a.min_sites)))
            elif not (rel >= a.swap_max):
                flags.append(dict(row, flag="swap", threshold=a.swap_max,
                                  detail="same patient but relatedness {} < {}: sample swap or wrong patient".format(fmt(rel), a.swap_max)))
        else:
            n_diff += 1
            if n < a.min_sites:
                n_low += 1
            elif rel >= a.dup_min:
                flags.append(dict(row, flag="duplicate", threshold=a.dup_min,
                                  detail="different patients but relatedness {} >= {}: same individual (duplicate or mislabel)".format(fmt(rel), a.dup_min)))

    run_has_y = any(fnum(s.get("gt_depth_mean")) > 0 and fnum(s.get("Y_depth_mean")) / fnum(s.get("gt_depth_mean")) >= 0.2
                    for s in samples)
    sex_rows = []
    n_sex_checked = 0
    for s in samples:
        sid = s["sample_id"]
        m = meta(sid)
        inferred, x_call, y_call, xr, yr, hf, inf = infer_sex(s, run_has_y, a.min_x_sites)
        sex_rows.append({"sample": sid, "patient": m["patient"], "status": m["status"], "declared_sex": m["sex"],
                         "inferred_sex": inferred, "x_call": x_call, "y_call": y_call, "x_depth_ratio": fmt(xr),
                         "y_depth_ratio": fmt(yr), "x_het_fraction": fmt(hf), "x_informative_sites": inf})
        declared = m["sex"]
        if declared == "unknown" or sid not in sheet:
            continue
        mismatch = None
        if m["status"] == "normal":
            if inferred in ("male", "female"):
                n_sex_checked += 1
                if inferred != declared:
                    mismatch = "normal inferred {} (X {} / Y {})".format(inferred, x_call, y_call)
        else:
            # Tumours lose Y and become X-hemizygous (LOH) often, so only evidence that cannot come from a loss
            # counts: Y reads in a declared female, X heterozygosity in a declared male.
            if x_call != "unknown" or y_call != "unknown":
                n_sex_checked += 1
            if declared == "female" and y_call == "male":
                mismatch = "tumour has Y coverage (Y/autosome depth {})".format(fmt(yr))
            elif declared == "male" and x_call == "female":
                mismatch = "tumour has X heterozygosity (het fraction {})".format(fmt(hf))
        if mismatch:
            flags.append({"flag": "sex_mismatch", "sample_a": sid, "sample_b": "", "patient_a": m["patient"],
                          "patient_b": "", "relatedness": "", "n_sites": inf, "threshold": "",
                          "detail": "declared {}; {}".format(declared, mismatch)})

    with open(a.flags, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FLAG_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(flags)
    with open(a.sex_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SEX_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(sex_rows)

    counts = {}
    for f in flags:
        counts[f["flag"]] = counts.get(f["flag"], 0) + 1
    print("somalier identity check: {} samples, {} same-patient pairs, {} cross-patient pairs ({} below {} sites), "
          "sex checked for {} samples; flags: {}".format(
              len(samples), n_same, n_diff, n_low, a.min_sites, n_sex_checked,
              ", ".join("{} {}".format(v, k) for k, v in sorted(counts.items())) or "none"))


if __name__ == "__main__":
    main()
