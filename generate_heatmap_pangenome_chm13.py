import sys
import os
import pandas as pd
import statistics
import seaborn as sns
import matplotlib.pyplot as plt

paf_file = sys.argv[1]
file_template = paf_file.split("fastq", 1)[0]
previous_read = ""
chromosomes = ["chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q","chr6p","chr6q","chr7p","chr7q","chr8p","chr8q","chr9p","chr9q","chr10p","chr10q","chr11p","chr11q","chr12p","chr12q","chr13p","chr13q","chr14p","chr14q","chr15p","chr15q","chr16p","chr16q","chr17p","chr17q","chr18p","chr18q","chr19p","chr19q","chr20p","chr20q","chr21p","chr21q","chr22p","chr22q","chrXp","chrXq","chrYp","chrYq"]

df60 = pd.DataFrame(0, index=range(len(chromosomes)), columns=chromosomes)
df15 = pd.DataFrame(0, index=range(len(chromosomes)), columns=chromosomes)
df1 = pd.DataFrame(0, index=range(len(chromosomes)), columns=chromosomes)
df0 = pd.DataFrame(0, index=range(len(chromosomes)), columns=chromosomes)

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

        # Extract true sequence identity
        title = read_name.split("!")
        true_chromosome = title[1]
        true_chromosome = true_chromosome.split("|")
        true_chromosome = true_chromosome[1]
        true_arm = true_chromosome.split(":")[1]
        true_arm = true_arm[0]
        true_chromosome = true_chromosome.split(":")[0]
        true_chr_arm = true_chromosome + true_arm

        # Update haploid accuracy df
        try:
            index = chromosomes.index(true_chr_arm)
        except ValueError:
            continue

        if mapped_chr_arm not in chromosomes:
            continue

        if map_Q >= 60:
            df60.loc[index,mapped_chr_arm] += 1
        elif map_Q >= 15:
            df15.loc[index,mapped_chr_arm] += 1
        elif map_Q >= 1:
            df1.loc[index,mapped_chr_arm] += 1
        else:
            df0.loc[index,mapped_chr_arm] += 1

        previous_read = read_name



labels = [s.split("chr")[1] for s in chromosomes]

df60 = df60.div(df60.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(14, 9))
heatmap = sns.heatmap(df60, cmap='coolwarm', xticklabels=labels, yticklabels=labels)
ax.tick_params(labelsize=12)
heatmap_name = file_template + "Q60.heatmap.png"
heatmap.figure.savefig(heatmap_name, format='png')

df15 = df15.div(df15.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(14, 9))
heatmap = sns.heatmap(df15, cmap='coolwarm', xticklabels=labels, yticklabels=labels)
ax.tick_params(labelsize=12)
heatmap_name = file_template + "Q15.heatmap.png"
heatmap.figure.savefig(heatmap_name, format='png')

df1 = df1.div(df1.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(14, 9))
heatmap = sns.heatmap(df1, cmap='coolwarm', xticklabels=labels, yticklabels=labels)
ax.tick_params(labelsize=12)
heatmap_name = file_template + "Q1.heatmap.png"
heatmap.figure.savefig(heatmap_name, format='png')

df0 = df0.div(df0.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(14, 9))
heatmap = sns.heatmap(df0, cmap='coolwarm', xticklabels=labels, yticklabels=labels)
ax.tick_params(labelsize=12)
heatmap_name = file_template + "Q0.heatmap.png"
heatmap.figure.savefig(heatmap_name, format='png')