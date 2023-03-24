import sys
import os

reference = sys.argv[1]
print(reference)

name_ref = reference.split("/", 1)[1]
name_ref = name_ref.split(".", 1)[0]

dir_path = "./haplotype_merge_fasta/"
file_paths = [os.path.join(dir_path, filename) for filename in os.listdir(dir_path) if filename.endswith('.fasta')]
file_paths = [x for x in file_paths if x != reference]

cmd = "mkdir {name_ref}_simulation".format(name_ref=name_ref)
os.system(cmd)

merged_ref = "./{name_ref}_simulation/{name_ref}.ref_genomes.merged.fasta".format(name_ref=name_ref)

with open(merged_ref, 'w') as output_file:
    for file_path in file_paths:
        with open(file_path, 'r') as input_file:
            output_file.write(input_file.read())

work_dir = "./{name_ref}_simulation".format(name_ref=name_ref)
os.chdir(work_dir)

print("in working directory")

cmd = "samtools faidx {name_ref}.ref_genomes.merged.fasta".format(name_ref=name_ref)
os.system(cmd)

print("merged fasta indexed")


# Simulate reads from given reference genome
cmd = "/meyersonlab/maxgj2/for_max/pbsim --depth 100 --length-min 5000 --length-mean 10000 --accuracy-mean 0.999 --hmm_model /meyersonlab/maxgj2/for_max/data/R94.model /meyersonlab/maxgj2/pangenome/{reference}".format(reference=reference)
os.system(cmd)

print("reads simulated")

cmd = "paftools.js pbsim2fq /meyersonlab/maxgj2/pangenome/{reference}.fai *.maf | pigz > ./{name_ref}.simulate.R94.depth100.accuracy_0_999.fastq.gz".format(reference=reference,name_ref=name_ref)
os.system(cmd)

print("reads compressed")

os.system("rm *.maf")
os.system("rm *.fastq")
os.system("rm *.ref")




#cmd = "perl /meyersonlab/maxgj2/for_max/extract_terminal_simulated_reads.v2.pl /meyersonlab/maxgj2/pangenome/{reference}.fai 1000000 ./{name_ref}.simulate.R94.depth100.accuracy_0_999.fastq.gz | pigz > ./{name_ref}.simulate.R94.depth100.accuracy_0_999.term_1000kb.fastq.gz".format(reference=reference,name_ref=name_ref)
#os.system(cmd)
cmd = "perl /meyersonlab/maxgj2/for_max/extract_terminal_simulated_reads.v2.pl /meyersonlab/maxgj2/pangenome/{reference}.fai 100000 ./{name_ref}.simulate.R94.depth100.accuracy_0_999.fastq.gz | pigz > ./{name_ref}.simulate.R94.depth100.accuracy_0_999.term_100kb.fastq.gz".format(reference=reference,name_ref=name_ref)
os.system(cmd)
cmd = "perl /meyersonlab/maxgj2/for_max/extract_terminal_simulated_reads.v2.pl /meyersonlab/maxgj2/pangenome/{reference}.fai 10000 ./{name_ref}.simulate.R94.depth100.accuracy_0_999.fastq.gz | pigz > ./{name_ref}.simulate.R94.depth100.accuracy_0_999.term_10kb.fastq.gz".format(reference=reference,name_ref=name_ref)
os.system(cmd)

# Align to chm13
cmd = "minimap2 -x map-ont -I 8G -t 20 /meyersonlab/maxgj2/for_max/ref_genomes/chm13.draft_v1.0.fasta {name_ref}.simulate.R94.depth100.accuracy_0_999.term_10kb.fastq.gz > ./{name_ref}.chm13.simulate.R94.depth100.accuracy_0_999.term_10kb.fastq.gz.paf".format(name_ref=name_ref)
os.system(cmd)
cmd = "minimap2 -x map-ont -I 8G -t 20 /meyersonlab/maxgj2/for_max/ref_genomes/chm13.draft_v1.0.fasta {name_ref}.simulate.R94.depth100.accuracy_0_999.term_100kb.fastq.gz > ./{name_ref}.chm13.simulate.R94.depth100.accuracy_0_999.term_100kb.fastq.gz.paf".format(name_ref=name_ref)
os.system(cmd)
cmd = "minimap2 -x map-ont -I 8G -t 20 /meyersonlab/maxgj2/for_max/ref_genomes/chm13.draft_v1.0.fasta {name_ref}.simulate.R94.depth100.accuracy_0_999.fastq.gz > ./{name_ref}.chm13.simulate.R94.depth100.accuracy_0_999.term_1000kb.fastq.gz.paf".format(name_ref=name_ref)
os.system(cmd)


# Align to all other pangenome
cmd = "minimap2 -x map-ont -I 8G -t 20 ./{name_ref}.ref_genomes.merged.fasta {name_ref}.simulate.R94.depth100.accuracy_0_999.term_10kb.fastq.gz > ./{name_ref}.pangenome_ref.simulate.R94.depth100.accuracy_0_999.term_10kb.fastq.gz.paf".format(name_ref=name_ref)
os.system(cmd)
cmd = "minimap2 -x map-ont -I 8G -t 20 ./{name_ref}.ref_genomes.merged.fasta {name_ref}.simulate.R94.depth100.accuracy_0_999.term_100kb.fastq.gz > ./{name_ref}.pangenome_ref.simulate.R94.depth100.accuracy_0_999.term_100kb.fastq.gz.paf".format(name_ref=name_ref)
os.system(cmd)
cmd = "minimap2 -x map-ont -I 8G -t 20 ./{name_ref}.ref_genomes.merged.fasta {name_ref}.simulate.R94.depth100.accuracy_0_999.fastq.gz > ./{name_ref}.pangenome_ref.simulate.R94.depth100.accuracy_0_999.term_1000kb.fastq.gz.paf".format(name_ref=name_ref)
os.system(cmd)



# Analyze mapping to chm13
cmd = "python /meyersonlab/maxgj2/for_max/evaluate_mapping_pangenome_chm13.py ./{name_ref}.chm13.simulate.R94.depth100.accuracy_0_999.term_10kb.fastq.gz.paf".format(name_ref=name_ref)
os.system(cmd)

cmd = "python /meyersonlab/maxgj2/for_max/evaluate_mapping_pangenome_chm13.py ./{name_ref}.chm13.simulate.R94.depth100.accuracy_0_999.term_100kb.fastq.gz.paf".format(name_ref=name_ref)
os.system(cmd)

cmd = "python /meyersonlab/maxgj2/for_max/evaluate_mapping_pangenome_chm13.py ./{name_ref}.chm13.simulate.R94.depth100.accuracy_0_999.term_1000kb.fastq.gz.paf &".format(name_ref=name_ref)
os.system(cmd)



# Analyze mapping to pangenome
cmd = "python /meyersonlab/maxgj2/for_max/evaluate_mapping_pangenome.py ./{name_ref}.pangenome_ref.simulate.R94.depth100.accuracy_0_999.term_10kb.fastq.gz.paf".format(name_ref=name_ref)
os.system(cmd)

cmd = "python /meyersonlab/maxgj2/for_max/evaluate_mapping_pangenome.py ./{name_ref}.pangenome_ref.simulate.R94.depth100.accuracy_0_999.term_100kb.fastq.gz.paf".format(name_ref=name_ref)
os.system(cmd)

cmd = "python /meyersonlab/maxgj2/for_max/evaluate_mapping_pangenome.py ./{name_ref}.pangenome_ref.simulate.R94.depth100.accuracy_0_999.term_1000kb.fastq.gz.paf &".format(name_ref=name_ref)
os.system(cmd)



