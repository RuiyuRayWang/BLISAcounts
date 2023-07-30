import pandas as pd

def wash_whitelist(whitelist_path,
                   well_ground_truth_path,
                   ab_ground_truth_path):

    whitelist = pd.read_csv(whitelist_path, sep='\t', names=["Barcode","Corrected","Count_BC","Count_Corrected"], 
                            header=None)
    well_ground_truth = pd.read_csv(well_ground_truth_path)
    ab_ground_truth = pd.read_csv(ab_ground_truth_path)

    whitelist['Well_BC'] = whitelist.Barcode.str.slice(0,6)
    whitelist['Ab_BC'] = whitelist.Barcode.str.slice(6,12)

    # Keep entries whose Well BC and Ab BC are in their respective ground truth lists
    whitelist_washed = whitelist.loc[whitelist['Well_BC'].isin(well_ground_truth['Well_BC']) & \
        whitelist['Ab_BC'].isin(ab_ground_truth['Ab_BC'])]
    
    whitelist_washed = whitelist_washed.drop(columns=['Well_BC','Ab_BC'])

    return(whitelist_washed)

whitelist_washed = wash_whitelist(whitelist_path = snakemake.input[0],
                                  well_ground_truth_path = snakemake.input[1],
                                  ab_ground_truth_path = snakemake.input[2])

whitelist_washed.to_csv(snakemake.output[0], index=False, sep="\t", header=False)
