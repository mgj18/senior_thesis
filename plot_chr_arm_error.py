import sys
import os
import subprocess
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

file1 = sys.argv[1]
file2 = sys.argv[2]

ref1 = file1.split("txt")[0]
ref2 = file2.split(".",1)[0]
distance = ref1.split(".map")[0]
distance = distance.split("term_")[1]
file_name = ref1 + ref2 +  ".chr_arm_error.png"

cmd = "head -n 210 {file1} | tail -n 48 | cut -f4".format(file1=file1)
result = subprocess.check_output(cmd, shell=True)
errors = [float(x) for x in result.decode('utf-8').strip().split()]

cmd = "head -n 210  {file1} | tail -n 48 | cut -f2".format(file1=file1)
result = subprocess.check_output(cmd, shell=True)
chromosomes = [x for x in result.decode('utf-8').strip().split()]
chromosomes = [s.lstrip("chr") for s in chromosomes]

df1 = pd.DataFrame({'chr':chromosomes,'error':errors})
df1 = df1[df1["chr"]!="Mp"]
df1 = df1[df1["chr"]!="Mq"]

cmd = "head -n 210 {file2} | tail -n 48 | cut -f4".format(file2=file2)
result = subprocess.check_output(cmd, shell=True)
errors = [float(x) for x in result.decode('utf-8').strip().split()]

cmd = "head -n 210  {file2} | tail -n 48 | cut -f2".format(file2=file2)
result = subprocess.check_output(cmd, shell=True)
chromosomes = [x for x in result.decode('utf-8').strip().split()]
chromosomes = [s.lstrip("chr") for s in chromosomes]

df2 = pd.DataFrame({'chr':chromosomes,'error':errors})
df2 = df2[df2["chr"]!="Mp"]
df2 = df2[df2["chr"]!="Mq"]

df1 = df1.sort_values(by='error')
df2 = df2.reindex(df1.index)





ind = np.arange(len(df1["error"]))  # the x locations for the groups
width = 0.35  # the width of the bars

fig, ax = plt.subplots(figsize=(20, 8))
rects1 = ax.bar(ind - width/2, df1["error"], width,
                label='HG00621 Aligned to Chm13')
rects2 = ax.bar(ind + width/2, df2["error"], width,
                label='HG00621 Aligned to Pangenome Reference')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Error Rate',fontsize=16)
ax.set_xlabel('Chromosomal Arm', fontsize=16)
simulated_ref = ref2.upper()
title_name = '{reference} R94 Simulated Reads (<{distance}) Error Rate by Chromosome Arm Identity'.format(reference = simulated_ref,distance=distance)
ax.set_title(title_name,fontsize=20)
ax.set_xticks(ind)
ax.set_xticklabels(df1["chr"])
plt.xticks(rotation=45, ha='right', fontsize=16)
ax.legend(fontsize=20)


fig.tight_layout()


plt.savefig(file_name, dpi=300)


