
import sys
import os
from collections import defaultdict

#program to process trimmed fastq files from Crispr experuments
#command line files inputs: trimmed fastq file and sgRNA librart
#command line string inputs: string for sample alignment statistics file (for output), 
#	string for sample counts file (for output),
#	string for trimmed counts file (for output)
#	string for sequence counts file (for output)

#create a dictionary with: key is sequence, value is array with guide name, sequence, and gene name
#file name is the comma delimited file with guidename, sequence, gene name
#length_dict is empty as input then populated
#dupe_seq is empty as input then populated
def create_library_dict(filename=None,length_dict=None,dupe_seq=None):
  d={}
  seq_array=[]
  with open(filename) as f:
    for line in f:
      sgrna = line.split(",")
      #add sgrna sequence
      if sgrna[1] in seq_array:
        if sgrna[1] not in dup_seq:
          dupe_seq.append(sgrna[1])
      else:
        seq_array.append(sgrna[1])
      #add length of sequence to sequence length dict
      seq_length=len(sgrna[1])
      length_dict.update({seq_length:0})
      d[sgrna[1]]= sgrna
  return d


#input is a trimmed fastq file
#input_dict key is sequence of sgRNA, value is a list with guide  name, sequence, gene target
#output is list with #unmapped reads, #reads, #multimapped reads, #reads with wrong length
def process_fastq(filename=None,input_dict=None,n=0,seq_length_counter=None,dup_seq=None):
  unmappedreads=0
  numreads=0
  nummulti=0
  num_wronglength=0
  with open(filename, 'r') as fh:
    lines = []
    for line in fh:
      lines.append(line.rstrip())
      if len(lines) == n:
        numreads+=1
        sequence = lines[1]
        if sequence in dup_seq:
          nummulti+=1
        elif sequence in input_dict:
          input_dict[sequence]+=1
        else:
          unmappedreads+=1
        if len(sequence) in seq_length_counter:
          seq_length_counter[len(sequence)]+=1
        else:
          num_wronglength+=1
        lines=[]
  return [unmappedreads,numreads,nummulti,num_wronglength]

#setting and checking for command line arguments
try:
    fastq = sys.argv[1]
except IndexError as ie:
    raise SystemError("Error: Specify fastq file name\n")

if not os.path.exists(fastq):
    raise SystemError("Error: fastq file does not exist\n")

try:
    lib_file = sys.argv[2]
except IndexError as ie:
    raise SystemError("Error: Specify library file name\n")

if not os.path.exists(lib_file):
    raise SystemError("Error: lib_file file does not exist\n")

try:
    statistics_outfile = sys.argv[3]
except IndexError as ie:
    raise SystemError("Error: Specify statistics out file name\n")


try:
    count_outfile = sys.argv[4]
except IndexError as ie:
    raise SystemError("Error: Specify count out file name\n")



#number of lines for a fastq file
record_lines = 4
reads_length_dict={}
dupe_sgrna_seq=[]
sgrna_sim_dict={}
#create a dictionary with all the guide sequences based on input file lib_file
#populate reads_length_dict, dupe_sgrna_seq, sgrna_sim_dict
sgrna_dict=create_library_dict(lib_file,reads_length_dict,dupe_sgrna_seq)
sgrna_seq = list(sgrna_dict.keys())
sample_dict = {key: 0 for key in sgrna_seq}

reads_summary=process_fastq(fastq,sample_dict,record_lines,reads_length_dict,dupe_sgrna_seq)

#output_filename=fastq[0:fastq.index('.')]

#write out alignment statistics
with open(statistics_outfile, 'w') as f:
    f.write('Number of processed reads\tNumber of unmapped reads\tNumber of multimapped reads')
    seq_lens=reads_length_dict.keys()
    for i in seq_lens:
        f.write('\tNumber of reads with length ' + str(i))
    f.write('\tNumber reads with different length\n')
    f.write(str(reads_summary[1]) + '\t' + str(reads_summary[0]) + '\t' + str(reads_summary[2]))
    for i in seq_lens:
        f.write('\t'+str(reads_length_dict[i]))
    f.write('\t'+str(reads_summary[3])+'\n')
f.close()

ds = [sgrna_dict, sample_dict]
sample_counter = defaultdict(list)

for d in (sgrna_dict,sample_dict): 
    for key, value in d.items():
        sample_counter[key].append(value)

df=open(count_outfile,"w")
for i in sgrna_seq:
  df.write(sample_counter[i][0][0] + "\t" + str(sample_counter[i][1]) + "\n")
df.close()


