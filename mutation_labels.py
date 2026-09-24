"""Mutation column -> protein label ("ESR1 p.S463P") from the scIMPACT classifier workbook.

Masterlist mutation columns are SYMBOL___chr___pos___ref___alt[___id]; the
workbook has one sheet per individual with a header row ("Source",
"theVariant", "SYMBOL", "AAchange", ...) where theVariant = chr___pos___ref___alt___id.
Keys are matched on SYMBOL + chr/pos/ref/alt (the trailing id is optional in
some column names). Output: mutation_labels.tsv
"""
import re, pandas as pd

import os
DATA = os.environ.get("EZ_DATA", "data")  # run from the repo root
XLSX = os.environ.get("SCIMPACT_XLSX", "../scIMPACT Mutations_Classifier-AllelicCounts_20250116ZCxlsx.xlsx")  # workbook not in the repo
rows = []
for sheet, df in pd.read_excel(XLSX, sheet_name=None, header=None).items():
    hdr = df.index[df.iloc[:, 1].astype(str).eq("theVariant")][0]
    t = df.iloc[hdr + 1:, :].copy(); t.columns = df.iloc[hdr]
    t = t[t["theVariant"].notna() & t["SYMBOL"].notna()]
    for _, r in t.iterrows():
        v = str(r["theVariant"]).strip()
        rows.append(dict(sheet=sheet, symbol=str(r["SYMBOL"]).strip(), variant=v,
                         key5=str(r["SYMBOL"]).strip() + "___" + "___".join(v.split("___")[:4]),
                         aachange=str(r["AAchange"]).strip(), consequence=r.get("Consequence")))
lab = pd.DataFrame(rows)
lab["label"] = lab.symbol + " p." + lab.aachange
# one label per variant key; flag keys whose label differs between sheets
agg = lab.groupby("key5").agg(label=("label", lambda s: "|".join(sorted(set(s)))),
                              sheets=("sheet", lambda s: ",".join(sorted(set(s)))),
                              consequence=("consequence", "first")).reset_index()
agg["conflict"] = agg.label.str.contains("|", regex=False)
# Corrections to the workbook, with the reason recorded in the output.
OVERRIDES = {
    # NRAS is on the minus strand; 114713909 G>T is Q61K (c.181C>A, in the workbook),
    # so 114713908 is c.182: T>C = c.182A>G, CAA->CGA = Q61R. The workbook's Q61P
    # would need T>G. Q61R is also the label used in the 07-2026 MM-127 figures.
    "NRAS___1___114713908___T___C": ("NRAS p.Q61R", "workbook says Q61P; codon T>C at c.182 gives Q61R"),
}
agg["note"] = ""
for k, (l, why) in OVERRIDES.items():
    agg.loc[agg.key5 == k, ["label", "note"]] = [l, why]
agg.to_csv(os.path.join(DATA, "mutation_labels.tsv"), sep="\t", index=False)
print(len(lab), "rows,", len(agg), "variants,", agg.conflict.sum(), "with conflicting labels")
print(agg[agg.label.str.contains(r"p\.(?:nan|0|\.)$|\|") | agg.note.ne("")].head(20).to_string())
