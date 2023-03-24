import sys
import os
import pandas as pd
import statistics

paf_file = sys.argv[1]
file_template = paf_file.split("fastq", 1)[0]
previous_read = ""

# df to store arm based evaluation including parental 
map_eval = pd.DataFrame(columns=['chr','Q','error'])
map_eval = map_eval.astype({'chr': str, 'Q': float, 'error': float})

with open(paf_file) as file:
    for alignment in file:
        line = alignment.strip().split("\t")
        read_name = line[0]
        
        if read_name == previous_read:
            continue


        # Extract sequence alignment identity
        mapped_chromosome = line[5]
        mapped_location = line[7]

        try:
            mapped_location = float(mapped_location)
        except ValueError:
            continue

        if mapped_location < 10000000:
            mapped_arm = "p"
        else:
            mapped_arm = "q"
        
        mapped_chr_arm = mapped_chromosome + mapped_arm

        map_Q = line[11]


        try:
            map_Q = float(map_Q)
        except ValueError:
            continue

        
        mapped_chr_arm = mapped_chromosome + mapped_arm

        # Extract true sequence identity
        title = read_name.split("!")
        true_chromosome = title[1]
        true_chromosome = true_chromosome.split("|")
        true_chromosome = true_chromosome[1]
        true_arm = true_chromosome.split(":")[1]
        true_arm = true_arm[0]
        true_chromosome = true_chromosome.split(":")[0]
        
        true_chr_arm = true_chromosome + true_arm


        # Update accuracy df 
        if mapped_chr_arm == true_chr_arm:
            error = 0
        else:
            error = 1

        new_row = {'chr': true_chr_arm,'Q': map_Q,'error':error}
        map_eval = map_eval.append(new_row, ignore_index=True)

        previous_read = read_name




all_arms = list(set(map_eval['chr']))
chr_numbers = [string[3:] for string in all_arms]
chr_order = sorted(range(len(chr_numbers)),key=lambda i: chr_numbers[i])
all_arms = [all_arms[i] for i in chr_order]

# Calculate overall alignment accuracy
overall_acc = pd.DataFrame(columns=['Q','reads','error_rate'])
overall_acc = overall_acc.astype({'Q': float, 'reads': float,'error_rate': float})

for q in [60, 15, 1, 0]:
    errors = map_eval.loc[map_eval['Q'] >= q, 'error']
    read_count = len(errors)
    
    if read_count > 0:
        rate = statistics.mean(errors)
    else:
            rate = 0
    new_row = {'Q': q,'reads': read_count,'error_rate':rate}
    overall_acc = overall_acc.append(new_row, ignore_index=True)

# Calculate arm by arm alignment accuracy
arm_acc = pd.DataFrame(columns=['Q','arm','reads','error_rate'])
arm_acc = arm_acc.astype({'Q': float,'arm': str, 'reads': float,'error_rate': float})

for q in [60, 15, 1, 0]:    
    for current_arm in all_arms:
        subset = map_eval[map_eval['chr'] == current_arm]
        errors = subset.loc[subset['Q'] >= q, 'error']
        read_count = len(errors)

        if read_count > 0:
            rate = statistics.mean(errors)
        else:
            rate = 0
        new_row = {'Q': q,'arm': current_arm,'reads': read_count,'error_rate':rate}
        arm_acc = arm_acc.append(new_row, ignore_index=True)


file_template = file_template + "map_eval.txt"


with open(file_template, 'a') as f:
    f.write("Overall Mapping Accuracy\n")
    overall_acc.to_csv(f, header=True, index=False, sep='\t')

    f.write("\n\nChromosomal Arm Mapping Accuracy\n")
    arm_acc.to_csv(f, header=True, index=False, sep='\t')



