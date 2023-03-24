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
# df to store arm based evaluation excluding parental
map_eval_hap = pd.DataFrame(columns=['chr','Q','error'])
map_eval_hap = map_eval_hap.astype({'chr': str,'Q': float, 'error': float})

with open(paf_file) as file:
    for alignment in file:
        line = alignment.strip().split("\t")
        read_name = line[0]
        
        if read_name == previous_read:
            continue

        # Extract sequence alignment identity
        mapped_chromosome = line[5]
        mapped_location = line[7]
        map_Q = line[11]


        try:
            map_Q = float(map_Q)
        except ValueError:
            continue

        try:
            mapped_location = float(mapped_location)
        except ValueError:
            continue

        if mapped_location < 10000000:
            mapped_arm = "p"
        else:
            mapped_arm = "q"
        
        mapped_chr_arm = mapped_chromosome + mapped_arm
        mapped_chr_arm_hap = mapped_chromosome + mapped_arm

        # Check is parental identity in alignment
        if "PATERNAL" in mapped_chromosome or "MATERNAL" in mapped_chromosome:
            mapped_chr_hap = mapped_chromosome.split("_", 1)[0]
            mapped_chr_arm_hap = mapped_chr_hap + mapped_arm


        # Extract true sequence identity
        title = read_name.split("!")
        true_chromosome = title[1]
        true_location = title[2]
        

        try:
            true_location = float(true_location)
        except ValueError:
            continue

        if true_location < 10000000:
            true_arm = "p"
        else:
            true_arm = "q"
        
        true_chr_arm = true_chromosome + true_arm
        true_chr_arm_hap = true_chromosome + true_arm

        # Check is parental identity in alignment
        if "PATERNAL" in true_chromosome or "MATERNAL" in true_chromosome:
            true_chr_hap = true_chromosome.split("_", 1)[0]
            true_chr_arm_hap = true_chr_hap + true_arm


        # Update diploid accuracy df 
        if mapped_chr_arm == true_chr_arm:
            error = 0
        else:
            error = 1

        new_row = {'chr': true_chr_arm,'Q': map_Q,'error':error}
        map_eval = map_eval.append(new_row, ignore_index=True)

        # Update haploid accuracy df
        if mapped_chr_arm_hap == true_chr_arm_hap:
            error = 0
        else:
            error = 1

        new_row = {'chr': true_chr_arm_hap,'Q': map_Q,'error':error}
        map_eval_hap = map_eval_hap.append(new_row, ignore_index=True)


        previous_read = read_name

all_arms = list(set(map_eval_hap['chr']))
chr_numbers = [string[3:] for string in all_arms]
chr_order = sorted(range(len(chr_numbers)),key=lambda i: chr_numbers[i])
all_arms = [all_arms[i] for i in chr_order]

# Calculate overall alignment accuracy neglecting parental identity
overall_acc_hap = pd.DataFrame(columns=['Q','reads','error_rate'])
overall_acc_hap = overall_acc_hap.astype({'Q': float, 'reads': float,'error_rate': float})

for q in [60, 15, 1, 0]:
    errors = map_eval_hap.loc[map_eval_hap['Q'] >= q, 'error']
    
    read_count = len(errors)
    if read_count > 0:
        rate = statistics.mean(errors)
    else:
            rate = 0
    new_row = {'Q': q,'reads': read_count,'error_rate':rate}
    overall_acc_hap = overall_acc_hap.append(new_row, ignore_index=True)

# Calculate arm by arm alignment accuracy neglecting parental identity
arm_acc_hap = pd.DataFrame(columns=['Q','arm','reads','error_rate'])
arm_acc_hap = arm_acc_hap.astype({'Q': float,'arm': str, 'reads': float,'error_rate': float})

for q in [60, 15, 1, 0]:
    for current_arm in all_arms:
        subset = map_eval_hap[map_eval_hap['chr'] == current_arm]
        errors = subset.loc[subset['Q'] >= q, 'error']

        read_count = len(errors)
        if read_count > 0:
            rate = statistics.mean(errors)
        else:
            rate = 0
        new_row = {'Q': q,'arm': current_arm,'reads': read_count,'error_rate':rate}
        arm_acc_hap = arm_acc_hap.append(new_row, ignore_index=True)

# Calculate overall alignment accuracy including parental idenitity
overall_acc_dip = pd.DataFrame(columns=['Q','reads','error_rate'])
overall_acc_dip = overall_acc_dip.astype({'Q': float, 'reads': float,'error_rate': float})

for q in [60, 15, 1, 0]:
    errors = map_eval.loc[map_eval['Q'] >= q, 'error']
    read_count = len(errors)
    
    if read_count > 0:
        rate = statistics.mean(errors)
    else:
            rate = 0
    new_row = {'Q': q,'reads': read_count,'error_rate':rate}
    overall_acc_dip = overall_acc_dip.append(new_row, ignore_index=True)

# Calculate arm by arm alignment accuracy including parental identity
arm_acc_dip = pd.DataFrame(columns=['Q','arm','reads','error_rate'])
arm_acc_dip = arm_acc_dip.astype({'Q': float,'arm': str, 'reads': float,'error_rate': float})

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
        arm_acc_dip = arm_acc_dip.append(new_row, ignore_index=True)


#standard_template = file_template + "map_eval"
file_template = file_template + "map_eval.txt"
#command = 'nohup bash -c "cat {paf} | paftools.js mapeval -m 1 - > {template}".format(paf=paf_file,template=standard_template) &'
#os.system(command)


with open(file_template, 'a') as f:
    f.write("Overall Mapping Accuracy (Haploid)\n")
    overall_acc_hap.to_csv(f, header=True, index=False, sep='\t')

    f.write("\n\nOverall Mapping Accuracy Requiring Correct Parental Assignment (Diploid)\n")
    overall_acc_dip.to_csv(f, header=True, index=False, sep='\t')

    f.write("\n\nChromosomal Arm Mapping Accuracy (Haploid)\n")
    arm_acc_hap.to_csv(f, header=True, index=False, sep='\t')

    f.write("\n\nChromosomal Arm Mapping Accuracy (Diploid)\n")
    arm_acc_dip.to_csv(f, header=True, index=False, sep='\t')






        
        


