from os import walk
import pandas as pd
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import os

from lifelines import KaplanMeierFitter

from lifelines import statistics as stats


def str_percentage_to_float(percentage):
    try:
        return float(percentage[:-1]) / 100
    except:
        return -1

def get_lifelines_format(df):
    indices = df["Index"]
    died = df["Died"]
    censored = df["Lost to Follow-up"]
    final_index = -1 if len(indices) == 0 else indices.iloc[-1]
    
    T = []
    E = []
    
    for ((_, index), (_, cur_death), (_, cur_censored)) in zip(indices.iteritems(), died.iteritems(), censored.iteritems()):
        
        if index == final_index:
            cur_censored += df["Start"].iloc[-1]
        
        T.extend([index] * (cur_death + cur_censored))
        E.extend([1] * cur_death + [0] * cur_censored)
    
    return pd.DataFrame({"T": T, "E": E})

for _, _, filenames in walk(snakemake.input[0]):
    filenames = filenames
    
pd_dict = {}
pd_dict_non_kmf = {}

for filename in tqdm(filenames):
    if filename[-3:] != "csv":
        continue
    splitted = filename.split("~")
    flag = splitted[0]
    cancer = "~".join(splitted[1:-1])
    
    if cancer not in pd_dict:
        pd_dict[cancer] = {}
        pd_dict_non_kmf[cancer] = {}
    temp_df = pd.read_csv(f"{snakemake.input[0]}/{filename}", converters={"Observed Cum": str_percentage_to_float}, thousands=",")
    temp_df = temp_df[temp_df["Observed Cum"] != -1]

    temp_df = temp_df.astype({"Died": "int32", "Lost to Follow-up": "int32", "Start": "int32"})
    lifelines_format = get_lifelines_format(temp_df)
    
    pd_dict[cancer][flag] = lifelines_format
    pd_dict_non_kmf[cancer][flag] = temp_df
    pd_dict_non_kmf[cancer][flag]["FLAG"] = flag


last_survival = {}
os.makedirs(snakemake.output[0], exist_ok=True)
for cancer, v in tqdm(pd_dict.items()):
    v_non_kmf = pd_dict_non_kmf[cancer]
    try:
        cur_df = pd.concat([cur_df_part for flag, cur_df_part in v_non_kmf.items() if "NON-ANALYTIC" not in flag], ignore_index=True)
        cur_dfs = {}
        for flag, cur_df_part in v.items():
            if "NON-ANALYTIC" in flag:
                continue
            cur_dfs[flag] = cur_df_part
    except ValueError:
        continue
    """plt.clf()
    sns.lineplot(x="Index", y="Observed Cum", data=cur_df, hue="FLAG", drawstyle='steps-pre')
    plt.ylim([0,1])
    plt.savefig(f"outputs/observed/{cancer}.png", dpi=600, bbox_inches="tight")"""
    plt.clf()
    kmf = KaplanMeierFitter()
    fig, ax = plt.subplots()
    
    for flag, cur_df in cur_dfs.items():
        if cur_df.shape[0] == 0:
            continue
        kmf.fit(cur_df["T"], event_observed=cur_df["E"], label=flag)
        kmf.plot_survival_function(ax=ax)
        if cancer not in last_survival:
            last_survival[cancer] = {}
        last_survival[cancer][flag] = pd.concat([kmf.survival_function_, kmf.confidence_interval_], axis = 1)
    
    """for cancer, sub_flag in cur_last_survival.items():
        for flag, surv_df in sub_flag.items():"""
            
    plt.savefig(f"{snakemake.output[0]}/{cancer}.png", dpi=600, bbox_inches="tight")

os.makedirs(snakemake.output[1], exist_ok=True)
for cur_cancer, flag_dict in last_survival.items():
    os.makedirs(f"{snakemake.output[1]}/{cur_cancer}", exist_ok=True)
    for cur_flag, cur_df in flag_dict.items():
        cur_df.to_csv(f"{snakemake.output[1]}/{cur_cancer}/{cur_flag}.csv", index=False)