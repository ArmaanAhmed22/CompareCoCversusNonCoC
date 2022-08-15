import pandas as pd
from os import walk
import seaborn as sns
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
COC = "ANALYTIC abstract from facility WITH CoC accreditation"
nCOC = "Abstract from facility WITHOUT CoC accreditation"


divisions_pre_df = {"cancer": [], "ratio": []}
for i, (dirpath, dirnames, filenames) in tqdm(enumerate(walk(snakemake.input[0]))):
    if len(filenames) == 2 and "Abstract from facility WITHOUT CoC accreditation.csv" in filenames and "ANALYTIC abstract from facility WITH CoC accreditation.csv":
        temp_df_dict = {}
        for filename in filenames:
            no_ext = filename.split(".")[0]
            cur_df = pd.read_csv(f"{dirpath}/{filename}")
            temp_df_dict[no_ext] = cur_df[no_ext]
        min_size_index = min([cur_df.shape[0] for cur_df in temp_df_dict.values()]) - 1
        divisions_pre_df["ratio"].append(temp_df_dict[nCOC].iloc[min_size_index] / temp_df_dict[COC].iloc[min_size_index])
        divisions_pre_df["cancer"].append(f"(month={min_size_index}) {dirpath.split('/')[-1]}")


out_df = pd.DataFrame(divisions_pre_df).sort_values(by="ratio", ascending=True)

fig, ax = plt.subplots(figsize=(20,6))
sns.barplot(x="cancer", y="ratio", data=out_df, color="blue", ax=ax)
sns.despine(ax=ax)
plt.xticks(rotation=90)
plt.savefig(snakemake.output[0], dpi=600, bbox_inches="tight")