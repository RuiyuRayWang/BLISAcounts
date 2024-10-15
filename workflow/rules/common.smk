import pandas as pd
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

# # Debug
# import yaml
# with open ('config/config.yaml') as f:
#     config = yaml.safe_load(f)

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")

def get_demultiplex_inputs():
    if config['staggered']:
        demultiplex_inputs = "workflow/data/{user}/{project}/{library}/outs/{sample}_staggered_multiplex",
        # demultiplex_inputs = expand(
        #     "workflow/data/{samples.user}/{samples.project}/{samples.library}/outs/{samples.sample}_staggered_multiplex",
        #     samples = samples.itertuples()
        # )
    else:
        demultiplex_inputs = "workflow/data/{user}/{project}/{library}/outs/{sample}_multiplex",
        # demultiplex_inputs = expand(
        #     "workflow/data/{samples.user}/{samples.project}/{samples.library}/outs/{samples.sample}_multiplex",
        #     samples = samples.itertuples()
        # )
    return demultiplex_inputs

def get_final_outputs():
    final_outputs = expand(
        "workflow/data/{samples.user}/{samples.project}/{samples.library}/outs/{samples.sample}_final_count",
        samples = samples.itertuples()
    )
    return final_outputs
