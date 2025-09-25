#!/usr/bin/env python
import argparse, os, pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

def calc_props(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {}
    mw = Descriptors.MolWt(m)
    logp = Crippen.MolLogP(m)
    tpsa = rdMolDescriptors.CalcTPSA(m)
    hbd = rdMolDescriptors.CalcNumHBD(m)
    hba = rdMolDescriptors.CalcNumHBA(m)
    rb = rdMolDescriptors.CalcNumRotatableBonds(m)
    return dict(MW=mw, logP=logp, TPSA=tpsa, HBD=hbd, HBA=hba, RB=rb)

def skin_rule(row):
    # Very rough topical heuristics (for screening only)
    return (row['MW'] < 500) and (row['logP'] >= 0 and row['logP'] <= 5) and (row['TPSA'] < 90)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--compounds', required=True, help="CSV with columns: name,smiles")
    ap.add_argument('--out', default='results')
    args = ap.parse_args()
    os.makedirs(f"{args.out}/tables", exist_ok=True)

    df = pd.read_csv(args.compounds)
    props = []
    for _, r in df.iterrows():
        p = calc_props(r['smiles'])
        p.update({'name': r['name']})
        props.append(p)
    props = pd.DataFrame(props)
    props['skin_ok'] = props.apply(skin_rule, axis=1)
    props.to_csv(f"{args.out}/tables/admet_skin.csv", index=False)
    print("Wrote ADMET/skin properties.")

if __name__ == "__main__":
    main()