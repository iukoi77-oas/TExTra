import pandas as pd
import numpy as np
import os
import subprocess as sp

def compute_key(df, bin_size=10000):
    # Compute a "rounded" column based on the given bin size and generate a unique key for each row
    df.loc[:, 'rounded'] = (df[1] / bin_size).astype(int)
    df.loc[:, 'key'] = df[3].astype(str) + "/" + df[0].astype(str) + "/" + df['rounded'].astype(str) + "/" + df[5].astype(str)
    return df

def match_reads_python(R1_path, R2_path, strandedness, matched_file, unmatched_file1, unmatched_file2, debug=False):
    # Load the read files
    R1 = pd.read_csv(R1_path, sep='\t', header=None)
    R2 = pd.read_csv(R2_path, sep='\t', header=None)

    # First round of matching: group reads by bin size of 10,000
    R1 = compute_key(R1, bin_size=10000)
    R2 = compute_key(R2, bin_size=10000)

    # Sort by key
    R1_sorted = R1.sort_values(by='key')
    R2_sorted = R2.sort_values(by='key')

    # Merge on key
    merged = pd.merge(R1_sorted, R2_sorted, on='key', suffixes=('_1', '_2'))

    # Filter based on strandedness rules
    def filter_conditions(row):
        s1, e1, start1, end1, strand1, flag1 = row[1], row[2], row[1], row[2], row[5], row[4]
        s2, e2, start2, end2, strand2, flag2 = row[13], row[14], row[13], row[14], row[17], row[15]
        insert1 = e1 - s2
        insert2 = e2 - s1

        polyA1 = (flag1 == 1000)
        polyA2 = (flag2 == 1000)

        # Strand matching rules
        if strandedness == 1:
            return (
                (insert1 <= 500 and insert1 >= 0 and start1 >= s2 and strand1 == '+' and not polyA1 and not polyA2) or
                (insert2 <= 500 and insert2 >= 0 and start2 >= s1 and strand1 == '-' and not polyA1 and not polyA2) or
                (strand1 == '+' and (
                    (insert1 <= 500 and insert1 >= 0 and start1 >= s2 and polyA1 and not polyA2) or
                    (insert2 <= 500 and insert2 >= 0 and start2 >= s1 and not polyA1 and polyA2) or
                    (insert1 <= 500 and insert1 >= 0 and start1 >= s2 and polyA1 and polyA2)
                )) or
                (strand1 == '-' and (
                    (insert1 <= 500 and insert1 >= 0 and start1 >= s2 and not polyA1 and polyA2) or
                    (insert2 <= 500 and insert2 >= 0 and start2 >= s1 and polyA1 and not polyA2) or
                    (insert2 <= 500 and insert2 >= 0 and start2 >= s1 and polyA1 and polyA2)
                ))
            )
        elif strandedness == 2:
            return (
                (insert2 <= 500 and insert2 >= 0 and start2 >= s1 and strand1 == '+' and not polyA1) or
                (insert1 <= 500 and insert1 >= 0 and start1 >= s2 and strand1 == '-' and not polyA1) or
                (strand1 == '+' and (
                    (insert1 <= 500 and insert1 >= 0 and start1 >= s2 and not polyA1 and polyA2) or
                    (insert2 <= 500 and insert2 >= 0 and start2 >= s1 and polyA1 and not polyA2) or
                    (insert2 <= 500 and insert2 >= 0 and start2 >= s1 and polyA1 and polyA2)
                )) or
                (strand1 == '-' and (
                    (insert1 <= 500 and insert1 >= 0 and start1 >= s2 and polyA1 and not polyA2) or
                    (insert2 <= 500 and insert2 >= 0 and start2 >= s1 and not polyA1 and polyA2) or
                    (insert2 <= 500 and insert2 >= 0 and start2 >= s1 and polyA1 and polyA2)
                ))
            )
        elif strandedness == 0:
            return (
                (insert2 <= 500 and insert2 >= 0 and start2 >= s1) or
                (insert1 <= 500 and insert1 >= 0 and start1 >= s2)
            )
        else:
            return False

    filtered = merged[merged.apply(filter_conditions, axis=1)]

    # Save matched reads
    # filtered.to_csv(matched_file, sep='\t', header=False, index=False)

    # Get unmatched reads
    matched_ids = set(filtered.iloc[:, 3])  # original read ID
    unmatched_R1 = R1[~R1[3].isin(matched_ids)].drop(columns=['key', 'rounded'])
    unmatched_R2 = R2[~R2[3].isin(matched_ids)].drop(columns=['key', 'rounded'])

    unmatched_R1.to_csv(unmatched_file1, sep='\t', header=False, index=False)
    unmatched_R2.to_csv(unmatched_file2, sep='\t', header=False, index=False)

    if debug:
        print("Matched:", len(filtered))
        print("Unmatched R1:", len(unmatched_R1))
        print("Unmatched R2:", len(unmatched_R2))

    # Second round of matching: increase bin size to 100,000 for unmatched reads
    unmatched_R1 = compute_key(unmatched_R1, bin_size=100000)
    unmatched_R2 = compute_key(unmatched_R2, bin_size=100000)

    # Sort and merge unmatched reads for second round
    R1_sorted = unmatched_R1.sort_values(by='key')
    R2_sorted = unmatched_R2.sort_values(by='key')

    merged_second_round = pd.merge(R1_sorted, R2_sorted, on='key', suffixes=('_1', '_2'))

    filtered_second_round = merged_second_round[merged_second_round.apply(filter_conditions, axis=1)]

    # Save second round matched reads
    # filtered_second_round.to_csv(matched_file, sep='\t', header=False, index=False)

    # Combine first and second round matches
    final_matched = pd.concat([filtered, filtered_second_round], ignore_index=True).drop(columns=['key', 'rounded_1', 'rounded_2'])

    # Save final matched reads
    final_matched.to_csv(matched_file, sep='\t', header=False, index=False)

    os.unlink(R1_path)
    os.unlink(R2_path)

    return matched_file, unmatched_file1, unmatched_file2


# match_reads_python(
#     R1_path='group1_rep1_paired_ulabeled_1.tmp305llqmd',
#     R2_path='group1_rep1_paired_ulabeled_2.tmp0p87dhpt',
#     strandedness=1,
#     matched_file='group1_rep1_paired_matched.tmp7jfc37yx',
#     unmatched_file1='group1_rep1_paired_matched.tmp7jfc37yx',
#     unmatched_file2='group1_rep1_paired_unmatched_2.tmpyshlru98'
# )