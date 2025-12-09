#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import pairwise2
from itertools import combinations_with_replacement
from itertools import permutations
import pandas as pd

# Levure Motif dégénéré T(G)2-3(TG)1-6
# motif = ['TGGGTGTG',
#            'TGGTGTGT',
#            'TGGTGTGG',
#            'TGTGTGTG']

name='TTTGGGG'
motif=['TTTGGGG']

k = len(motif[0])

min_score = k/2

Perfect = list(set([motif[x][i:] + motif[x][:i] for x in range(len(motif)) for i in range(len(motif[x])) for x in range(len(motif))]))


A=['A','T','G','C']

Tab_score=pd.DataFrame(columns=['Mer','score'])

for mer in list(combinations_with_replacement(A,k)):
    for p_mer in set(list(permutations(mer))):
        p_mer=''.join(p_mer)
        if p_mer not in list(Tab_score['Mer']):
            s_max=0
            for tel_mer in Perfect:
                s_ali=pairwise2.align.globalxx(tel_mer,p_mer,score_only=True)
                if s_ali>s_max:
                    s_max=s_ali
            i=len(Tab_score)
            
            Tab_score.at[i,'Mer']=p_mer
            Tab_score.at[i,'score']=max(s_max,min_score)
    

Tab_score.to_csv('Score_for_'+name+'.tab',sep='\t',index=False)
