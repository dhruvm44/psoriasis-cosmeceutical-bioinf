#!/usr/bin/env python
import argparse, os, pandas as pd, numpy as np

def score_connectivity(up_list, down_list, ref_signature):
    # Very simplified "connectivity": negative correlation between ranks
    # ref_signature: DataFrame with 'gene','rank' (+ up if high in compound)
    ref = ref_signature.set_index('gene')['rank']
    s_up = ref.reindex(up_list).dropna().mean() if up_list else np.nan
    s_dn = ref.reindex(down_list).dropna().mean() if down_list else np.nan
    # We want disease-up genes to be low in compound (negative) and disease-down to be high (positive)
    return -s_up + s_dn

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--up', required=True, help="file with disease up genes, one per line")
    ap.add_argument('--down', required=True, help="file with disease down genes")
    ap.add_argument('--out', default='results')
    args = ap.parse_args()
    os.makedirs(f"{args.out}/tables", exist_ok=True)

    up = [g.strip() for g in open(args.up) if g.strip()]
    dn = [g.strip() for g in open(args.down) if g.strip()]

    # Placeholder compound/mechanism reference signature (rank: smaller is better for 'upregulation')
    ref = pd.DataFrame({
        'gene': ["IL17A","TNF","KRT6A","CXCL8","PPARG","CNR2"],
        'rank':  [  10,     12,     8,        14,     90,     85]
    })

    score = score_connectivity(up, dn, ref)
    pd.DataFrame([{'score': score}]).to_csv(f"{args.out}/tables/connectivity_placeholder.csv", index=False)
    print("Wrote placeholder connectivity score. Replace ref signature with LINCS/perturbation proxy.")

if __name__ == "__main__":
    main()