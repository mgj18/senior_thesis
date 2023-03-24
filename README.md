# Repository Description
These python scripts were written and used for analysis in the creation of Max Garrity-Janger's Harvard Applied Math senior thesis 3/24/2023

## Script descriptions and Utilization Instructions:

### extract_exact_sequence.py
This script was used to simulated so called "HiFi" or 0 error reads. In order to execute the script takes the path to a fastq file containing PBSIM2 simulated reads as command line input as well as the reference from which the reads were simulated. It then parses through the fastq file read by read, extracting the true exact sequence of each read from the reference genome using samtools faidx and writing these exact sequences to a new fastq file. Thus, after running this script you are left with the same reads as in the initial fastq file but with any simulated sequencing error removed. 

example usage: "python {fastq file} {reference genome} >> {desired name for new exact read sequence fastq}



### extract_exact_sequence.py
This script was used to simulated so called "HiFi" or 0 error reads. In order to execute the script takes the path to a fastq file containing PBSIM2 simulated reads as command line input as well as the reference from which the reads were simulated. It then parses through the fastq file read by read, extracting the true exact sequence of each read from the reference genome using samtools faidx and writing these exact sequences to a new fastq file. Thus, after running this script you are left with the same reads as in the initial fastq file but with any simulated sequencing error removed. 

example usage: "python extract_exact_sequence.py {fastq file} {reference genome} >> {desired name for new exact read sequence fastq}



### evaluate_mapping_complex.py
This script was used to evaluate the mapping accuracy after the alignment of PBSIM2 simulated reads. Specifically, it was used to evaluate mappings from any combination of CHM13 and HG002 reference genomes as simulation or alignment references. It generates a file with overall and arm specific mapping error rates. It assumes reference names of the form "chr1", "chr2", or "chr1_MATERNAL", "chr1_PATERNAL" etc.

example usage: "python evaluate_mapping_complex.py {.paf mapping file to be analyzed}"




### evaluate_mapping_pangenome.py
This script was used to evaluate the mapping accuracy after the alignment of PBSIM2 simulated reads. Specifically, it was used to evaluate mappings from one of the pangenome references to the entire pangenome. It generates a file with overall and arm specific mapping error rates. It assumes the reference name form of the specific pangenome used

example usage: "python evaluate_mapping_pangenome.py {.paf mapping file to be analyzed}"




### evaluate_mapping_pangenome_chm13.py
This script was used to evaluate the mapping accuracy after the alignment of PBSIM2 simulated reads. Specifically, it was used to evaluate mappings from one of the pangenome references to chm13. It generates a file with overall and arm specific mapping error rates. It assumes the reference name form of the specific pangenome for the true arm origin and the form "chr1", "chr2", etc. for the mapped arm (consistent with Chm13).

example usage: "python evaluate_mapping_pangenome_chm13.py {.paf mapping file to be analyzed}"




### analyze_pangenome.py
This script runs a single reference from the full pangenome reference through a complete analysis. It simulates reads using PBSIM2, collects these reads into a single .fastq file using pbsim2fq from paftools.js, then extracts reads <10,100,or 1000kb from the chr end and maps these reads to both chm13 and a merged pangenome of all other references in the pangenome besides the current simulation reference. It then analyzes the accuracy of each mapping using evaluate_mapping_pangenome.py and evaluate_mapping_pangenome_chm13.py

example usage: "python analyze_pangenome.py {path to single reference from pangenome for simulation}"



### generate_alignment_heatmap.py
This script makes a chr arm specific heat map of alignment from a .paf mapping file and saves it as .png file in the current directory assuming an HG002 or Chm13 formatted reference naming system e.g. "chr1", "chr2", etc.

example usage: "python generate_alignment_heatmap {path to .paf mapping file to be analyzed}"


### generate_heatmap_pangenome_pangenome.py
Same as generate_alignment_heatmap.py but assumes different reference format consistent with pangenome naming conventions. To be run on a .paf file generated from mapping pangenome reference simulated reads (single reference) to pangenome

example usage: "python generate_pangenome_pangenome {path to .paf mapping file to be analyzed}"



### generate_heatmap_pangenome_chm13.py
Same as generate_alignment_heatmap.py but assumes different reference format consistent with pangenome naming conventions for simulation reference. To be run on a .paf file generated from mapping single reference from pangenome simulated reads to chm13

example usage: "python generate_heatmap_pangenome_chm13 {path to .paf mapping file to be analyzed}"




### plot_chr_arm.py
Generates chr arm specific error rate plot comparing two different mapping accuracy files generated with evaluate_mapping_complex.py. Orders side by side bar plot in ascending error rate order of first of two files 

example usage: "python plot_chr_arm.py {mapping evaluation file 1} {mapping evaluation file 2}"



### U20S.analyze.py
This script was used to map the U20S reads to Chm13 and the pangenome, then parse through a file containing the telomeric length content of each reads and tabulate the telomeric read length distribution by mapped chr. Then it generates a side-by-side box plot of mapped telomeric read length distribution for each chr (for both alignment to pangenome and chm13) as well as a barplot of mapped read count by chr arm (for both alignment to pangenome and chm13).

example usage: "python U20S.analyze.py {file containing estimated telomeric content of each U20S read}"
