import sys
import pandas as pd

# usage: python merge_genes.py <out.tsv> <abs.tsv> <relab.tsv> <category>
out, abs_file, relab_file, category = sys.argv[1:5]

abs_df = pd.read_csv(abs_file, sep="\t", index_col=0)
relab_df = pd.read_csv(relab_file, sep="\t", index_col=0)

if list(abs_df.columns) != list(relab_df.columns):
    sys.exit(f"ERROR: sample columns differ between {abs_file} and {relab_file}")

suffix = f"_{category}"
columns = {}
for col in abs_df.columns:
    sample = col[:-len(suffix)] if col.endswith(suffix) else col
    columns[f"{sample}_abs"] = abs_df[col]
    columns[f"{sample}_relab"] = relab_df[col]

merged = pd.DataFrame(columns).fillna(0)     # outer join on gene name
merged.index.name = "Gene"
merged.to_csv(out, sep="\t")
print(f"Wrote {out} ({merged.shape[0]} genes, {len(abs_df.columns)} samples)")