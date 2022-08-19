import pandas as pd
from os import walk
import seaborn as sns
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from lifelines.statistics import survival_difference_at_fixed_point_in_time_test
from lifelines.fitters.kaplan_meier_fitter import KaplanMeierFitter
COC = "ANALYTIC abstract from facility WITH CoC accreditation"
nCOC = "Abstract from facility WITHOUT CoC accreditation"


divisions_pre_df = {"cancer": [], "ratio": []}
p_values_pre_df = {"cancer": [], "p-value": []}
for i, (dirpath, dirnames, filenames) in tqdm(enumerate(walk(snakemake.input[0]))):
    if len(filenames) == 2 and "Abstract from facility WITHOUT CoC accreditation.csv" in filenames and "ANALYTIC abstract from facility WITH CoC accreditation.csv":
        temp_df_dict = {}
        models = {}
        for filename in filenames:
            no_ext = filename.split(".")[0]
            cur_df = pd.read_csv(f"{dirpath}/{filename}")
            temp_df_dict[no_ext] = cur_df[no_ext]
            
            model_df = pd.read_csv(f"{snakemake.input[1]}/{dirpath.split('/')[-1]}/{filename}")
            models[no_ext] = model_df
        min_size_index = min([cur_df.shape[0] for cur_df in temp_df_dict.values()]) - 1
        divisions_pre_df["ratio"].append(temp_df_dict[nCOC].iloc[min_size_index] / temp_df_dict[COC].iloc[min_size_index])
        divisions_pre_df["cancer"].append(f"(month={min_size_index}) {dirpath.split('/')[-1]}")

        COC_model = models[COC]
        COC_kmf = KaplanMeierFitter()
        COC_kmf.fit(COC_model["T"], COC_model["E"])

        nCOC_model = models[nCOC]
        nCOC_kmf = KaplanMeierFitter()
        nCOC_kmf.fit(nCOC_model["T"], nCOC_model["E"])
        results = survival_difference_at_fixed_point_in_time_test(min_size_index, COC_kmf, nCOC_kmf)
        p_values_pre_df["cancer"].append(dirpath.split('/')[-1])
        p_values_pre_df["p-value"].append(results.p_value)


out_df = pd.DataFrame(divisions_pre_df).sort_values(by="ratio", ascending=True)

fig, ax = plt.subplots(figsize=(20,6))
sns.barplot(x="cancer", y="ratio", data=out_df, color="blue", ax=ax)
sns.despine(ax=ax)
plt.xticks(rotation=90)
plt.savefig(snakemake.output[0], dpi=600, bbox_inches="tight")

pd.DataFrame(p_values_pre_df).to_csv(snakemake.output[1], index=False)