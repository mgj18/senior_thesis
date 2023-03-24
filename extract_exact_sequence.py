import sys
import os
import gzip
import subprocess


reads = sys.argv[1]
ref = sys.argv[2]

file_template = reads.split("depth100", 1)
file_template[1] = file_template[1].split("term",1)[1]
write_file = file_template[0] + "depth100.exact_accuracy_1_0.term" + file_template[1]
#os.system("touch %s" %(write_file))


with gzip.open(reads, "rt") as file:
    for header, sequence in zip(*[iter(file)]*2):
        header = header.strip()
        title = header.split("!")
        chromosome = title[1]
        start = title[2]
        end = title[3]
        full_loc = chromosome + ":" + start + "-" + end
        orientation = title[4]

        print(header)

        if orientation == "-":
            cmd = "samtools faidx -i %s %s" %(ref, full_loc)

            # execute the command and capture the output
            result = subprocess.check_output(cmd, shell=True)

            # decode the output from bytes to string
            exact_seq = result.decode('utf-8')
            exact_seq = exact_seq.split('\n',1)[1]
            exact_seq = exact_seq.replace(' ', '').replace('\t', '').replace('\n', '')
            print(exact_seq)


        else:
            cmd = "samtools faidx %s %s" %(ref, full_loc)

            # execute the command and capture the output
            result = subprocess.check_output(cmd, shell=True)

            # decode the output from bytes to string
            exact_seq = result.decode('utf-8')
            exact_seq = exact_seq.split('\n',1)[1]
            exact_seq = exact_seq.replace(' ', '').replace('\t', '').replace('\n', '')
            print(exact_seq)

