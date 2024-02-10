#Python script to get convert the DNA sequences in RefSeq to RNA Sequences.

import re

#We want to input the sequence of letters, and get rid of all the spaces and numbers

test_sequence = input("Copy and Paste the full DNA sequence from RefSeq Database: ")

def DNAtoRNA(DNA):
        DNA = DNA.replace(' ','')
        #Now, we want to remove all numbers from the sequence
        pattern = r'[0-9]'
        DNA = re.sub(pattern, '', DNA)
        DNA = DNA.upper()
        DNA = DNA.replace('T','U')
        return str(DNA)


#Now we want to pick out the index of the first start codon (AUG), and delete it and everything before$

RNA_seq = DNAtoRNA(test_sequence)
index = RNA_seq.index('AUG')
index = int(index)

RNA_seq_nostart = RNA_seq[index+3:]
#print(RNA_seq_nostart)

j = 0

#creating a for loop that scans every third index for one of the stop codons
while j<len(RNA_seq_nostart):
    if RNA_seq_nostart[j:j+3] != 'UAA' and RNA_seq_nostart[j:j+3] != 'UAG' and RNA_seq_nostart[j:j+3] !='UGA':
        j+=3
    elif RNA_seq_nostart[j:j+3] == 'UAA':
        print('Final Sequence: ', RNA_seq_nostart[:j])
        break
    elif RNA_seq_nostart[j:j+3] == 'UAG':
        print('Final Sequence: ', RNA_seq_nostart[:j])
        break
    elif RNA_seq_nostart[j:j+3] == 'UGA':
        print('Final Sequence: ', RNA_seq_nostart[:j])
        break



