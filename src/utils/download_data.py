#!/usr/bin/env python
import argparse, os, sys, pathlib

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--accessions', nargs='+', required=True)
    p.add_argument('--outdir', default='data/raw')
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    # NOTE: This is a placeholder. Use GEOquery in R or GEOparse in Python.
    for acc in args.accessions:
        path = pathlib.Path(args.outdir)/f"{acc}.placeholder.txt"
        path.write_text(f"Download {acc} here. Replace this with real downloader.")
        print(f"[stub] wrote {path}")

if __name__ == "__main__":
    main()