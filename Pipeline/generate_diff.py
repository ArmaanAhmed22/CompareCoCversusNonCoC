import pandas as pd
from os import walk
import seaborn as sns
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
COC = "ANALYTIC abstract from facility WITH CoC accreditation"
nCOC = "Abstract from facility WITHOUT CoC accreditation"

os.makedirs(snakemake.output[0], exist_ok=True)
for i, (dirpath, dirnames, filenames) in tqdm(enumerate(walk(snakemake.input[0]))):
    if len(filenames) == 2 and "Abstract from facility WITHOUT CoC accreditation.csv" in filenames and "ANALYTIC abstract from facility WITH CoC accreditation.csv":
        dfs = {}
        for filename in filenames:
            name = filename.split(".")[0]
            dfs[name] = pd.read_csv(f"{dirpath}/{filename}")
        
        min_size = min(dfs[COC].shape[0], dfs[nCOC].shape[0])
        dfs[COC] = dfs[COC].head(n=min_size)
        dfs[nCOC] = dfs[nCOC].head(n=min_size)
        
        value_out = dfs[COC][COC] - dfs[nCOC][nCOC] 
        stdCOC = dfs[COC][f"{COC}_upper_0.95"] - dfs[COC][COC]
        stdnCOC = dfs[nCOC][f"{nCOC}_upper_0.95"] - dfs[nCOC][nCOC]
        varCOC = stdCOC.apply(lambda x: x**2)
        varnCOC = stdnCOC.apply(lambda x: x**2)
        std_out = (varCOC + varnCOC).apply(lambda x: x**0.5)

        time_range = range(min_size)
        diff_df = pd.DataFrame({"time": time_range, "data": value_out, "error": std_out})
        upper = value_out + std_out
        lower = value_out - std_out
        plt.clf()
        fig, ax = plt.subplots()

        sns.lineplot(x="time", y="data", data=diff_df, ax=ax, drawstyle="steps-pre")
        sns.lineplot(x=time_range, y=upper, ax=ax, alpha=0, drawstyle="steps-pre")
        sns.lineplot(x=time_range, y=lower, ax=ax, alpha=0, drawstyle="steps-pre")
        ax.fill_between(time_range, lower, upper, alpha=0.2, step="pre")
        plt.savefig(f"{snakemake.output[0]}/{dirpath.split('/')[-1]}.png", bbox_inches="tight")