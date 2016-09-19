CODON_TO_AA = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
               'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
               'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
               'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
               'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
               'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
               'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
               'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
               'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
               'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
               'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
               'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
               'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*'}

STOP_CODONS = set(['TAA', 'TAG', 'TGA'])  # set with stop codons

NT_TO_NUM = {'A' : 0, 'G' : 1, 'C' : 2, 'T' : 3}
NUM_TO_NT = {0 : 'A', 1 : 'G', 2 : 'C', 3 : 'T'}



#ATG is methionine
#TAA, TAG, TGA are stop codons
#A = 0, G = 1, C = 2, T = 3

def seq_to_num(input):
	newList = list(input)
	return newList

#input: genome as string
#output: ???

#input: sequence as numeric list
#output: reverse complement as numeric list
#method: reverse and 3-x
#rc = [3 - x for x in seq]

#length
#gc content (G's + C's)/total
#codon bias: # of times codon used/total # of amino acid
#protein sequence analysis
#molecular weight, hydrophobicity, aromaticity, basicity/polarity
#quality?