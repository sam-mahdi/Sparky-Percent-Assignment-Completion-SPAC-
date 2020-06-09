##Program calculates the number of alphas, betas, gammas, deltas, and epsilons you have assigned based on your proteins sequence
import re

alphas_betas={'D':2,'T':2,'S':2,'E':2,'P':2,'G':1,'A':2,'C':2,'V':2,'M':2,'I':2,'L':2,'Y':2,'F':2,'H':2,'K':2,'R':2,'W':2,'Q':2,'N':2}
gammas={'T':1,'E':1,'P':1,'V':1,'M':1,'I':2,'L':1,'K':1,'R':1,'Q':1}
deltas={'P':1,'M':1,'I':1,'L':1,'K':1,'R':1}
epsilons={'M':1,'K':1}

sequence_file=input('type in the name of your file that contains your proteins sequence (e.g. 1XXH.fasta.txt): ')
assignment_file=input('type in the name of your peaklist file: ')

alphas_betas_count=0
gammas_count=0
deltas_count=0
epsilons_count=0
#loads up sequence file, removes anything above the > seen in fafsta files
with open(sequence_file) as f:
    sequence = f.read().strip().upper()
    if sequence.startswith('>') or ('\n>' in sequence):
        sequence_file = sequence[sequence.find('>'):].split('\n', 1)[1]
    else:
        sequence_file = sequence
seq_file = "".join(sequence_file)
#using the above dictionary, these calculate how many alphas, betas, gammas, detlas, and epsilons you have based off
#your sequence
for amino_acid in seq_file:
    alpha_beta = alphas_betas.get(amino_acid)
    gamma = gammas.get(amino_acid)
    delta = deltas.get(amino_acid)
    epsilon = epsilons.get(amino_acid)
    if alpha_beta is not None:
        alphas_betas_count+=alpha_beta
    if gamma is not None:
        gammas_count+=gamma
    if delta is not None:
        deltas_count+=delta
    if epsilon is not None:
        epsilons_count+=epsilon

alphas_assigned=0
betas_assigned=0
gammas_assigned=0
deltas_assigned=0
epsilons_assigned=0
#counds how many alphas,betas...you have assigned in your peaklist
#all secondary assignments are removed so each CD2, CG2, is counted only as 1 assignment
with open(assignment_file) as file:
  for letter in file:
      A=letter.strip()
      betas=re.sub('HB2','',A)
      alphas=re.sub('HA2','',A)
      gamma_carbons=re.sub('CG2','',A)
      gamma_hydrogens=re.sub('HG2','',gamma_carbons)
      delta_carbons=re.sub('CD2','',A)
      delta_hydrogens=re.sub('HD2','',delta_carbons)
      epsilon_hydrogens=re.sub('HE2','',A)
      if re.findall(r'\BCE-HE',epsilon_hydrogens):
        epsilons_assigned+=1
      if re.findall(r'\BCD-HD', delta_hydrogens):
          deltas_assigned+=1
      if re.findall(r'\BCG-HG', gamma_hydrogens):
          gammas_assigned+=1
      if re.findall(r'\BCB-HB', betas):
        betas_assigned+=1
      if re.findall(r'\BCA-HA', alphas):
        alphas_assigned+=1

percent_alphas_betas_assigned=int(((alphas_assigned+betas_assigned)/alphas_betas_count)*100)
percent_gammas_assigned=int((gammas_assigned/gammas_count)*100)
percent_deltas_assigned=int((deltas_assigned/deltas_count)*100)
percent_epsilons_assigned=int((epsilons_assigned/epsilons_count)*100)
print(f' Percentage of alphas and betas assigned: {percent_alphas_betas_assigned}% \n',
f'Percentage of gammas assigned: {percent_gammas_assigned}% \n',f'Percentage of deltas assigned:  {percent_deltas_assigned}% \n',
f'Percentage of epsilons assigned: {percent_epsilons_assigned}%')
