import gzip
import regex
from Bio import SeqIO
from Levenshtein import distance as lev_dist
import pandas as pd
import numpy as np
from datetime import datetime


def parse_matching_pattern(anchor_1_pattern, anchor_2_pattern, tail_pattern):
    pattern = (r'[ATCGN]*'
               '(?e)(?P<anchor1>' + anchor_1_pattern + '){i<=1,d<=1,e<=3,2i+2d+1s<=5}'
               '(?P<well_bc>[ATCGN]{6})'  # must be 6, otherwise anchor2 contains insertions
               '(?e)(?P<anchor2>' + anchor_2_pattern + '){i<=1,d<=1,e<=3,2i+2d+1s<=5}'
               '(?P<ab_bc>[ATCGN]{5,7})'
               '(?P<umi>[ATCGN]{15})'
               '(?e)(?P<tail>' + tail_pattern + '){s<=1}')
    return pattern


def correct_bc(bc, ground_truth_list, count_correct, count_fail, err_dist = 2):
    #
    d = [lev_dist(bc, bc_GT) for bc_GT in ground_truth_list]
    
    if min(d) > err_dist:  # If Levenshtein distance is greater than set err_dist, throw failure
        count_fail += 1
        success = False
    elif sorted(d)[0] == sorted(d)[1]:  # In case of tied Levenshtein distance, throw failure
        count_fail += 1
        success = False
    else:
        count_correct += 1
        bc = ground_truth_list[np.argmin(d)]
        success = True
        
    return bc, count_correct, count_fail, success


def approx_match(input_file,
                path_well_bc_GT,  # Path to well barcode ground truth list
                path_ab_bc_GT,
                anchor_1_pattern = "GACAAGTGGCCACAAACCACCAG",
                anchor_2_pattern = "CTTGTGGAAAGGACGAAACA", 
                tail_pattern = "GTCTGGAGCATGCG"):
    
    # Read well barcode ground truth
    with open(path_well_bc_GT) as file:
        well_bc_GT = [line.rstrip('\n') for line in file]
        
    # Read antibody barcode ground truth
    with open(path_ab_bc_GT) as file:
        ab_bc_GT = [line.rstrip('\n') for line in file]
    
    # Parse regex pattern
    pattern = parse_matching_pattern(anchor_1_pattern, 
                                     anchor_2_pattern, 
                                     tail_pattern)
    
    lst = []
    
    count_regex_unmatch = 0
    count_correct_fail = 0
    count_well_bc_correction = 0
    count_well_bc_fail = 0
    count_ab_bc_correction = 0
    count_ab_bc_fail = 0
    
    with gzip.open(input_file, "rt") as in_handle:
        for i, record in enumerate(SeqIO.parse(in_handle, 'fastq')):
            success_well = False
            success_ab = False
            
            if i % 1000000 == 0:
                print(datetime.now().strftime("%D %H:%M:%S") + ' Parsed ' + str(i) + ' reads. ' +
                      'Discarded ' + str(count_regex_unmatch) + ' reads due to unmatching regex.\n' +
                      '\t\t  Corrected ' + str(count_well_bc_correction) + ' well barcodes. ' +
                      str(count_well_bc_fail) + ' well barcodes failed to correct.\n' +
                      '\t\t  Corrected ' + str(count_ab_bc_correction) + ' antibody barcodes. ' +
                      str(count_ab_bc_fail) + ' antibody barcodes failed to correct.\n' +
                      '\t\t  Discarded ' + str(count_correct_fail) + ' reads due to failed correction.\n'
                      )

            # m = regex.match(pattern, str(record.seq), overlapped=False) # As of 21.10.15, overlapped in regex.match is no longer supported
            m = regex.match(pattern, str(record.seq))
            
            if m is not None:
                # Parse well barcodes
                well_bc = m.group('well_bc')
                
                if well_bc in well_bc_GT:
                    success_well = True
                else:
                    well_bc, count_well_bc_correction, count_well_bc_fail, success_well = correct_bc(well_bc, well_bc_GT, count_well_bc_correction, count_well_bc_fail)
                
                # Parse antibody barcodes
                ab_bc = m.group('ab_bc')
                
                if ab_bc in ab_bc_GT:
                    success_ab = True
                else:
                    ab_bc, count_ab_bc_correction, count_ab_bc_fail, success_ab = correct_bc(ab_bc, ab_bc_GT, count_ab_bc_correction, count_ab_bc_fail)
                
                umi = m.group('umi')
                
                # Actually make the correction
                if success_well and success_ab:
                    lst.append([well_bc, ab_bc, umi])
                else:
                    count_correct_fail += 1
                
            else:
                count_regex_unmatch += 1
    
    cols = ['well_bc', 'ab_bc', 'umi']
    df = pd.DataFrame(lst, columns = cols)
    
    return df

df = approx_match(snakemake.input[0], path_well_bc_GT=snakemake.input[1], path_ab_bc_GT=snakemake.input[2])

df.to_csv(snakemake.output[0], header=False, index=False, sep="\t")
