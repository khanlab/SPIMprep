import pickle

with open(snakemake.input.reg_result_pkl, "rb") as f:
    reg_result = pickle.load(f)


fig = reg_result["pairwise_registration"]["summary_plot"][0]

# Resize the figure before saving
fig.set_size_inches(20, 30)
fig.savefig(snakemake.output.pairwise_plot, bbox_inches="tight", dpi=300)


fig = reg_result["groupwise_resolution"]["summary_plot"][0]

# Resize the figure before saving
fig.set_size_inches(20, 30)
fig.savefig(snakemake.output.groupwise_plot, bbox_inches="tight", dpi=300)
