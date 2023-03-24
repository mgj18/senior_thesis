import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

read_lengths = sys.argv[1]

dir_path = "../haplotype_merge_fasta/"
file_paths = [os.path.join(dir_path, filename) for filename in os.listdir(dir_path) if filename.endswith('.fasta')]

merged_ref = "./all_pangenome_ref_genomes.merged.fasta"

with open(merged_ref, 'w') as output_file:
    for file_path in file_paths:
        with open(file_path, 'r') as input_file:
            output_file.write(input_file.read())

cmd = "samtools faidx all_pangenome_ref_genomes.merged.fasta"
os.system(cmd)

cmd = "minimap2 -x map-ont -I 8G -t 80 /meyersonlab/maxgj2/for_max/ref_genomes/chm13.draft_v1.0.fasta U2OS.nanopore.telomere_fixed.telomeric.sorted.fasta > U2OS.nanopore.telomere_fixed.telomeric.sorted.fasta.chm13.paf" 
os.system(cmd)

cmd = "minimap2 -x map-ont -I 8G -t 80 {merged_ref} U2OS.nanopore.telomere_fixed.telomeric.sorted.fasta > U2OS.nanopore.telomere_fixed.telomeric.sorted.fasta.pangenome.paf".format(merged_ref=merged_ref)
os.system(cmd)

chromosomes = ["1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13p","13q","14p","14q","15p","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21p","21q","22p","22q","Xp","Xq","Yp","Yq"]

reads_chm13 = pd.DataFrame(columns=['read_name','chr','length'])
reads_pan = pd.DataFrame(columns=['read_name','chr','length'])
previous_read = ""

# Generate chm13 boxplot
with open("U2OS.nanopore.telomere_fixed.telomeric.sorted.fasta.chm13.paf") as file:
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
        mapped_chr_arm = mapped_chr_arm.split("chr")[1]

        if mapped_chr_arm not in chromosomes:
        	continue

        new_row = {'read_name': read_name,'chr': mapped_chr_arm,'length': 0}
        reads_chm13 = reads_chm13.append(new_row, ignore_index=True)

        previous_read = read_name

# Generate chm13 boxplot
with open("U2OS.nanopore.telomere_fixed.telomeric.sorted.fasta.pangenome.paf") as file:
    for alignment in file:
        line = alignment.strip().split("\t")
        read_name = line[0]
        
        if read_name == previous_read:
            continue

        # Extract sequence alignment identity
        mapped_chromosome = line[5]
        mapped_chromosome = mapped_chromosome.split("|")
        mapped_chromosome = mapped_chromosome[1]
        mapped_arm = mapped_chromosome.split(":")[1]
        mapped_arm = mapped_arm[0]
        mapped_chromosome = mapped_chromosome.split(":")[0]
        mapped_chr_arm = mapped_chromosome + mapped_arm
        mapped_chr_arm = mapped_chr_arm.split("chr")[1]

        if mapped_chr_arm not in chromosomes:
        	continue

        new_row = {'read_name': read_name,'chr': mapped_chr_arm,'length': 0}
        reads_pan = reads_pan.append(new_row, ignore_index=True)

        previous_read = read_name


with open(read_lengths) as file:
	for line in file:
		strings = ['left_telomeric','right_telomeric']
		line = line.strip().split("\t")
		read_name = line[0]
		tel_length = line[3]
		classification = line[8]

		try:
			tel_length = float(tel_length)
		except ValueError:
			continue

		if classification not in strings:
			continue

		reads_chm13.loc[reads_chm13['read_name'] == read_name, 'length'] = tel_length
		reads_pan.loc[reads_pan['read_name'] == read_name, 'length'] = tel_length


reads_chm13 = reads_chm13[reads_chm13['length'] != 0]
reads_pan = reads_pan[reads_pan['length'] != 0]


# Generate boxplot for chm13 alignment
fig, ax = plt.subplots(figsize=(20, 8))
ax.grid(True)
ax = sns.boxplot(x='chr', y='length', data=reads_chm13, order=chromosomes,width=0.7)

# rotate the x-axis labels for better readability
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right',fontsize=14)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=16)
ax.set_xlabel("Chromosome Arm", fontsize=16)
ax.set_ylabel("Telomere Length (bp)", fontsize=16)
ax.set_title("Telomere Length Distribution for Chm13 Alignment", fontsize=20)


plt.savefig('tel_length_distribution.chm13.png', dpi=300)


# Generate boxplot for pangenome alignment
fig, ax = plt.subplots(figsize=(20, 8))
ax.grid(True)
ax = sns.boxplot(x='chr', y='length', data=reads_pan, order=chromosomes,width=0.7)

# rotate the x-axis labels for better readability
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right',fontsize=14)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=16)
ax.set_xlabel("Chromosome Arm", fontsize=16)
ax.set_ylabel("Telomere Length (bp)", fontsize=16)
ax.set_title("Telomere Length Distribution for Pangenome Alignment", fontsize=20)


plt.savefig('tel_length_distribution.pangenome.png', dpi=300)




# Generate Read Count Plot for chm13
counts = reads_chm13['chr'].value_counts()
df = pd.DataFrame({'chromosome': counts.index, 'read_count': counts.values})

fig, ax = plt.subplots(figsize=(20, 6))
sns.barplot(x='chromosome', y='read_count', data=df, order=chromosomes, color='blue', ax=ax)

# slant and increase font size of x-axis labels
plt.xticks(rotation=45, ha='right', fontsize=12)

# set axis labels and title
ax.set_xlabel('Chromosomal Arm', fontsize=16)
ax.set_ylabel('Read Count', fontsize=16)
ax.set_title('Chr Arm Read Count for Chm13 Alignment', fontsize=20)

plt.savefig('read_count_distribution.chm13.png')




# Generate Read Count Plot for pangenome
counts = reads_pan['chr'].value_counts()
df = pd.DataFrame({'chromosome': counts.index, 'read_count': counts.values})

fig, ax = plt.subplots(figsize=(20, 6))
sns.barplot(x='chromosome', y='read_count', data=df, order=chromosomes, color='blue', ax=ax)

# slant and increase font size of x-axis labels
plt.xticks(rotation=45, ha='right', fontsize=12)

# set axis labels and title
ax.set_xlabel('Chromosomal Arm', fontsize=16)
ax.set_ylabel('Read Count', fontsize=16)
ax.set_title('Chr Arm Read Count for Pangenome Alignment', fontsize=20)

plt.savefig('read_count_distribution.pangenome.png')











