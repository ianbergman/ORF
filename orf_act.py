from sys import argv

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

NT_TO_COMP = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}

script, input_file = argv


# store genome as single string
# input: fasta file
# output: genome as string
def gen_to_string(file):
     with open(file, "rU") as f:
          seq = ""
          for line in f:
            line = line.rstrip("\n")  # remove newline from end of line
            if not line.startswith(">"):  # new sequence header
                seq += line  # append line to end of sequence

          print("I HAVE CHOMPED THE SEQUENCE")
          return seq  # yield last sequence in file


# Return reverse complement of sequence
# input: sequence as string
# output: reverse complement as string
def rev_comp(input):
     newlist = reversed(list(input))
     chars = [NT_TO_COMP[x] for x in newlist]
     return ''.join(chars)


# Extract open reading frames from a given DNA Sequence
# Input: DNA sequence as string
# Output: List of lists of start:int, stop:int, orf:string
def extract(seq):
     out = []
     start = len(seq) + 1  # running in reverse, so start at end
     stop = 0  # edge case: found methionine before stop codon
     for i in range(len(seq) - 1, 1, -3):
          if seq[i - 2:i + 1] in STOP_CODONS:
               stop = i + 1
          elif seq[i - 2:i + 1] == 'ATG':
               start = i - 2
               if start < stop:
                    # print(start, stop)
                    # print[seq[start:stop]]
                    out.append((start, stop - 1, seq[start:stop]))
     return out[::-1]  # put back in order


# Compute GC content of ORF
# Input: ORF as string sequence
# Output: Percentage of GC bases as string fraction and float
def gc_content(seq):
     chars = list(seq[2])
     gc = 0
     for char in chars:
          if char == 'G' or char == 'C':
               gc += 1
     out = "GC Content of ORF " + str(seq[0]) + "," + str(seq[1])
     out += " is " + str(gc) + "/" + str(len(seq[2]))
     out += " " + str(float(gc) / float(len(seq[2])) * 100)[:4] + "%"
     return out


# Translate ORFs to protein seqs
# Input: ORF as string
# Output: Protein sequence as string
def orf_to_prot(seq):
     codons = [seq[2][i:i + 3] for i in range(0, len(seq[2]), 3)]
     aas = [CODON_TO_AA[x] for x in codons]
     return (seq[0], seq[1], ''.join(aas))

seq = gen_to_string(input_file)  # store the sequence
comp_seq = rev_comp(seq)  # reverse complement the sequence
orfs = []  # no orfs yet
for i in range(0, 3):  # check all offset reading frames
     orfs.append(extract(seq[:(len(seq) - i)]))

for i in range(0, 3):  # check all offset reading frames
     orfs.append(extract(comp_seq[:(len(seq) - i)]))
gc = [[gc_content(x) for x in y] for y in orfs]  # calculate gc's
prots = [[orf_to_prot(x) for x in y] for y in orfs]

# print(prots)
# print(gc[:])
# print(list('steps'))

# gc content (G's + C's)/total
# codon bias: # of times codon used/total # of amino acid
# protein sequence analysis
# molecular weight, hydrophobicity, aromaticity, basicity/polarity
# quality
