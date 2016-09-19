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

NT_TO_COMP = {'A' : 'T', 'G' : 'C', 'C' : 'G', 'T' : 'A'}

script, input_file = argv

#ATG is methionine
#TAA, TAG, TGA are stop codons
#A = 0, G = 1, C = 2, T = 3

#store genome as single string
#input: fasta file
#output: genome as string
def gen_to_string(file):
     with open(file, "rU") as f:
          seq = ""
          for line in f:
            line = line.rstrip("\n")  # remove newline from end of line
            if not line.startswith(">"):  # new sequence header
                seq += line  # append line to end of sequence

          print("I HAVE CHOMPED THE SEQUENCE")
          return seq  # yield last sequence in file



#Return reverse complement of sequence
#input: sequence as string
#output: reverse complement as string
def rev_comp(input):
     newlist = reversed(list(input))
     chars = [NT_TO_COMP[x] for x in newlist]
     return ''.join(chars)

#Extract start codon indices from list of codons
#input: list of codons
#output: list of indices associated with start codons
def find_m(codons):
     out = []
     for k,v in enumerate(codons):
          if v == 'M':
               out.append(k)
     return out

#Extract stop codons indices from list of codons
#input: list of codons
#output: list of indices associated with stop codons
def find_stop(codons):
     out = []
     for k,v in enumerate(codons):
          if v == '*':
               out.append(k)
     return out

#return an orf
#input: sequence
#output: orf
def find_orfs(seq):
     orfs = []
     starts = []
     stops = []
     rev_seq = rev_comp(seq)
     codons = [[seq[i:i+3] for i in range(0,len(seq)-2,3)]]
     print("FIRST SET OF CODONS CALCULATED, COMMANDER")
     codons.append([seq[i:i+3] for i in range(1,len(seq)-2,3)])
     print("SECOND SET OF CODONS CALCULATED, COMMANDER")
     codons.append([seq[i:i+3] for i in range(2,len(seq)-2,3)])
     print("THIRD SET OF CODONS CALCULATED, COMMANDER")
     codons.append([rev_seq[i:i+3] for i in range(0,len(seq)-2,3)])
     print("FOURTH SET OF CODONS CALCULATED, COMMANDER")
     codons.append([rev_seq[i:i+3] for i in range(1,len(seq)-2,3)])
     print("FIFTH SET OF CODONS CALCULATED, COMMANDER")
     codons.append([rev_seq[i:i+3] for i in range(2,len(seq)-2,3)])
     print("FINAL SET OF CODONS CALCULATED, COMMANDER")
     aas = [[CODON_TO_AA[x] for x in codons[0]]]
     print("FIRST SET OF AAS CALCULATED, COMMANDER")
     aas.append([CODON_TO_AA[x] for x in codons[1]])
     print("SECOND SET OF AAS CALCULATED, COMMANDER")
     aas.append([CODON_TO_AA[x] for x in codons[2]])
     print("THIRD SET OF AAS CALCULATED, COMMANDER")
     aas.append([CODON_TO_AA[x] for x in codons[3]])
     print("FOURTH SET OF AAS CALCULATED, COMMANDER")
     aas.append([CODON_TO_AA[x] for x in codons[4]])
     print("FIFTH SET OF AAS CALCULATED, COMMANDER")
     aas.append([CODON_TO_AA[x] for x in codons[5]])
     print("FINAL SET OF AAS CALCULATED, COMMANDER")

     for step,seqs in enumerate(codons):
          orfs.append([])
          print("ASSEMBLING READING FRAMES FOR STEP " + str(step) + ", COMMANDER")
          for i,codon in enumerate(seqs):
               open_frame = ""
               if codon == 'ATG':
                    print(i)
                    open_frame += codon
                    for j,cdn in enumerate(seqs[i+1:]):
                         open_frame += cdn
                         if cdn in STOP_CODONS:
                              orfs[step].append((open_frame,(step%3+3*i,step%3+3*(i+j+1))))
                              break
     print(orfs)
     #print(codons)

#find_orfs("ATGTGA")
#find_orfs("TCACAT")

find_orfs("ATGTACCGTATGCAGTAGCAGAGATTTCCAGATATTTGCCCTAA")
find_orfs("TTAGGGCAAATATCTGGAAATCTCTGCTACTGCATACGGTACAT")
find_orfs(gen_to_string(input_file))
#length
#gc content (G's + C's)/total
#codon bias: # of times codon used/total # of amino acid
#protein sequence analysis
#molecular weight, hydrophobicity, aromaticity, basicity/polarity
#quality?