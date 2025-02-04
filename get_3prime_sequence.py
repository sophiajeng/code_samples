import re 
import pandas as pd
 
#read in fasta
hcmv_fasta = ""
with open('/Users/jengs/OneDrive - Oregon Health & Science University/From_Box/Pre-Share/RNA190611JN/TB40E-GFP parental virus FASTA.fasta') as f:
    	hcmv_fasta_list= f.readlines()
hcmv_fasta_list = hcmv_fasta_list[1:]
for line in hcmv_fasta_list:
	hcmv_fasta += line

#forward strand stop codons
stop_codon = "TAA"
stop_codon2 = "TAG"
stop_codon3 = "TGA"
forward_stop_codons=(stop_codon,stop_codon2,stop_codon3)
#forward strand polya
polyA = "AATAAA"

#reverse strand stop codons
reverse_codon = "TTA"
reverse_codon2 = "CTA"
reverse_codon3 = "TCA"
reverse_stop_codons=(reverse_codon,reverse_codon2,reverse_codon3)
#reverse strand polyA
reverse_polyA = "TTTATT"

#read in gtf file
hcmv_gtf=pd.read_csv('/Users/jengs/OneDrive - Oregon Health & Science University/From_Box/Pre-Share/PPG/CEL220308JN/python/tb40e_threeprime_refined2.gtf',sep="\t",header=None)

#for each entry in gtf, subset fasta and search for stop codon and polyA
try:
	with open('hcmv_three_prim_sequences.txt','w') as output:
		output.write("HCMV_Transcript\tStrand\tStop_Codon\tStatus\tTranscript_Threeprime_Start\tTranscript_Threeprime_End\tSequence\n")
		for idx,gene_entry in hcmv_gtf.iterrows():
			strand=gene_entry[6]
			gene_fasta=hcmv_fasta[gene_entry[3]:gene_entry[4]]
			#positive strand then look for stop codon then polyA
			if strand=="+":
				for _, codon in enumerate(forward_stop_codons):
					m=re.search(f"{codon}(.+?){polyA}",gene_fasta)
					if m:
						output.write(gene_entry[8]+"\t")
						output.write(gene_entry[6]+"\t")
						output.write(codon+"\t")
						output.write("Found Stop codon and polyA\t")
						output.write(str(m.span()[0])+"\t")
						output.write(str(m.span()[1])+"\t")
						output.write(gene_fasta[m.span()[0]:m.span()[1]]+"\n")
					else:
						m2=re.search(f"{codon}",gene_fasta)
						if m2:
							output.write(gene_entry[8]+"\t")
							output.write(gene_entry[6]+"\t")
							output.write(codon+"\t")
							output.write("Found Stop codon not polyA\t")
							output.write(str(m2.span()[0])+"\t")
							output.write(str(m2.span()[0]+500)+"\t")
							output.write(gene_fasta[m2.span()[0]:m2.span()[0]+500]+"\n")
						else:
							output.write(gene_entry[8]+"\t")
							output.write(gene_entry[6]+"\t")
							output.write(codon+"\t")
							output.write("No stop codon"+"\t\t\t\n")
			#reverse strand look for polyA then stop codon
			else:
				for _, rev_codon in enumerate(reverse_stop_codons):
					m=re.search(f"{reverse_polyA}(.+?){rev_codon}",gene_fasta)
					if m:
						output.write(gene_entry[8]+"\t")
						output.write(gene_entry[6]+"\t")
						output.write(rev_codon+"\t")
						output.write("Found Stop codon and polyA\t")
						output.write(str(m.span()[0])+"\t")
						output.write(str(m.span()[1])+"\t")
						output.write(gene_fasta[m.span()[0]:m.span()[1]]+"\n")
					else:
						m2=re.search(f"{rev_codon}",gene_fasta)
						if m2:
							output.write(gene_entry[8]+"\t")
							output.write(gene_entry[6]+"\t")
							output.write(rev_codon+"\t")
							output.write("Found Stop codon not polyA\t")
							#output.write(str(m2.span()[0])+"\t")
							if (m2.span()[0]-500)>=0:
								output.write(str(m2.span()[0]-500)+"\t")
								output.write(str(m2.span()[1])+"\t")
								output.write(gene_fasta[m2.span()[0]-500:m2.span()[1]]+"\n")
							else:
								output.write(str(0)+"\t")
								output.write(str(m2.span()[1])+"\t")
								output.write(gene_fasta[0:m2.span()[1]]+"\n")
						else:
							output.write(gene_entry[8]+"\t")
							output.write(gene_entry[6]+"\t")
							output.write(rev_codon+"\t")
							output.write("No stop codon"+"\t\t\t\n")
	output.close()
except IOError:
	print("Error")



 
