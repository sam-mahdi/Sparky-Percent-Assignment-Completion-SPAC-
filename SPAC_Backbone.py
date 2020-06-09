##This program calculates the percantage of the backbone amides completed using your sequence
##It does exclude prolines

import re

dict={'D':1,'T':1,'S':1,'E':1,'G':1,'A':1,'C':1,'V':1,'M':1,'I':1,'L':1,'Y':1,'F':1,'H':1,'K':1,'R':1,'W':1,'Q':1,'N':1}


sequence_file=input('type in the name of your file that contains your proteins sequence (e.g. 1XXH.fasta.txt): ')
assignment_file=input('type in the name of your peaklist file: ')

backbone_count=0
with open(sequence_file) as f:
    sequence = f.read().strip().upper()
    if sequence.startswith('>') or ('\n>' in sequence):
        sequence_file = sequence[sequence.find('>'):].split('\n', 1)[1]
    else:
        sequence_file = sequence
seq_file = "".join(sequence_file)
for amino_acid in seq_file:
    backbone = dict.get(amino_acid)
    if backbone is not None:
        backbone_count+=backbone
backbone_assigned=0
with open(assignment_file) as file:
  for letter in file:
      A=letter.strip()
      if re.findall(r'\BN-HN',A):
        backbone_assigned+=1
percent_backbone_assigned=int((backbone_assigned/backbone_count)*100)
print(f' Percentage of backbone assigned: {percent_backbone_assigned}%')
