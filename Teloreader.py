#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd 
import numpy as np
import argparse
import time
import os
import concurrent.futures
import subprocess
import re
from pathlib import Path

start_time = time.time()

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Teloreader.py <strain> <fastafile> <motif>')
    parser.add_argument(
        "strain",
        help="strain name",
    )

    parser.add_argument(
        "fastafile",
        help="fasta file of reads",
    )
    parser.add_argument(
        "motif",
        help="telomeric motif (possible motifs = Yeast / TTAGGG / TTTAGGG / TTAGG / TTGGG / TTTGGG / TTTGGGG / )",
    )
    parser.add_argument(
        "-threads",
        "-t",
        default=4,
        type=int,
        help="Number of threads. default=4",
    )
    parser.add_argument(
        "-max_split_seq",
        "-ms",
        default=1000,
        type=int,
        help="Maximal number of sequences by spliting files. default=1000",
    )
    parser.add_argument(
        "-size_window",
        "-w",
        default=15,
        type=int,
        help="Size of the sliding window for the calculation of the score average (must be smaller than max_size_gap) default=15",
    )

    parser.add_argument(
        "-max_size_gap",
        "-g",
        default=20,
        type=int,
        help="Maximum size of a non telomeric gap in telomere default=20",
    )
    parser.add_argument(
        "-min_mean_window",
        "-mw",
        default=0,
        type=float,
        help="Minimum average score of a telomere window default=mer_size-1",
    )
    parser.add_argument(
        "-min_mean_telo",
        "-mint",
        default=0,
        type=float,
        help="Minimum average score of a telomere default=mer_size*0.8",
    )
    parser.add_argument(
        "-min_len",
        "-minl",
        default=16,
        type=int,
        help="Minimum length of a telomere default=16",
    )
    parser.add_argument(
        "-max_dist_term",
        "-maxt",
        default=50,
        type=int,
        help="Distance from the extremity to be considered as terminal",
    )
    parser.add_argument(
        "-out_pwd",
        "-o",
        default='None',
        help="Directory of output result",
    )

    return parser.parse_args()

def Make_reverse(seq):
    '''
    Takes sequence and returns the reverse complement 

    :param seq: a nucleotide sequence
    :return: the reverse sequence
    '''
    rev={'A':'T','T':'A','C':'G','G':'C'}
    rseq=''
    for i in range(len(seq)-1,-1,-1):
        rseq+=rev[seq[i]]
    return(rseq)


def Fasta_to_Kmertype(fasta,fileKmer,typ,mer_size,Tab_score):
    '''
    Converts a fasta file to an K-mer score file.

    :param fasta: path to input fasta file
    :param fileKmer: path to output K-mer score file
    :param typ: type of scoring system ('C' for C-rich or 'G' for G-rich)
    :param mer_size: size of telomere motif
    :param Tab_score: table of K-mer score
    '''
    f=open(fasta,'r')
    fK=open(fileKmer,'w')
    rest=''
    for l in f :
        if l[0]=='>':
            fK.write(l)
            rest=''
            name = l
        else:
            fK.write(l)
            if l[-1]!='\n':
                fK.write('\n')
            if rest=='':
                fK.write(' '*(mer_size-1))
            if len(rest)>1:
                seq=rest+l.strip('\n')
            else:
                seq=l.strip('\n')
            for i in range(len(seq)-(mer_size-1)):
                mer=seq[i:i+mer_size].upper()
                if re.fullmatch("[ATGC]+", mer):
                    if typ=='C' or typ=='c': # if C rich scoring, reverse complement the K-mer
                        try:    
                            mer=Make_reverse(mer)
                        except KeyError as e:
                            print(f"Error processing file: {fasta}, read: {name}, Error: {str(e)}")
                    fK.write(str(np.base_repr(int(Tab_score.at[mer,'score']),36)))
                else :
                    fK.write(str(0))
            rest=seq[-(mer_size-1):]
            fK.write('\n')
    f.close()
    fK.close()


def Add_Telo(dfTelo,nam,deb,fin,typ,len_read,mean_score,strain,max_dist_term):
    '''
    Adds a telomere to a telomere DataFrame.

    :param dfTelo: DataFrame to add telomere to
    :param nam: name of telomere
    :param deb: start position of telomere
    :param fin: end position of telomere
    :param typ: type of telomere ('C' for C-rich or 'G' for G-rich)
    :param len_read: length of read
    :param mean_score: mean 8-mer score of telomere
    :param strain: strain / label name
    :param max_dist_term: maximal distance from the extremity to be considered as terminal
    :return: updated DataFrame with new telomere
    '''
    ind=len(dfTelo)
    dfTelo.at[ind,'strain']=strain
    dfTelo.at[ind,'name']=nam
    dfTelo.at[ind,'type']=typ
    dfTelo.at[ind,'len']=fin-deb
    dfTelo.at[ind,'start']=deb+1
    dfTelo.at[ind,'end']=fin
    if deb<max_dist_term or fin+max_dist_term > len_read : 
        dfTelo.at[ind,'Loc']='term'
    else:
        dfTelo.at[ind,'Loc']='intern'
    dfTelo.at[ind,'Score_Kmer']=mean_score
    dfTelo.at[ind,'reads_len']=len_read
    return(dfTelo)

def Cal_all_mean(score,size_window,mer_size,min_mean_window):
    '''
    Calculates mean scores for all windows of given size in a list of scores,
    and returns a list of decisions based on whether the mean score in each
    window is above a given threshold.

    :param score: list of integer scores
    :param size_window: size of the sliding window for the calculation of the score average
    :param mer_size: size of telomere motif
    :param min_mean_window: Minimum average score of a telomere window
    :return: list of decisions (0 if mean score is below threshold, 1 if above)
    '''
    mean_score=[]
    for i in range(len(score)-size_window+mer_size):
        mean_score.append(np.mean([int(k,36) for k in score[i:i+size_window-mer_size+1]]))
    decision=[]
    for ms in mean_score:
        if float(ms)<float(min_mean_window):
            decision.append(0)
        else:
            decision.append(1)
    return(decision)

def Eval_One_Read(dfTelo,nam,seq,typ,score,mer_size,size_window,min_mean_window,min_mean_telo,strain,max_dist_term,max_size_gap,min_len):
    '''
    Given a score of a sequencing read, identifies and evaluates telomeric regions.
    Adds the telomeric regions to the data frame dfTelo, with the corresponding read name and type.

    :param dfTelo: Pandas DataFrame to which the telomeric regions will be added
    :param nam: read name
    :param seq: DNA sequence of the read
    :param typ: telomere type, 'G' Grich and 'C' for Crich
    :param score: list of Kmer score (in string format) of the sequence
    :param mer_size: size of telomere motif
    :param size_window: size of the sliding window for the calculation of the score average
    :param min_mean_window: minimum average score of a telomere window
    :param min_mean_telo: minimum average score of a telomere
    :param strain: strain / label name
    :param max_dist_term: maximal distance from the extremity to be considered as terminal
    :param max_size_gap: maximum size of a non telomeric gap in telomere
    :param min_len : minimum length of a telomere
    :returns dfTelo: the updated Pandas DataFrame with the telomeric regions added
    '''
    #At least one perfect telomeric motif. If this is not the case, the function does not execute the rest of the code.
    if str(np.base_repr(mer_size,36)) in score :
        telo=0 
        i=0
        start_gap=''
        decision=Cal_all_mean(score,size_window,mer_size,min_mean_window)
        if typ == 'G' : # if we search a G-rich telomere, we go through the sequence from the end
            ind_seq=range(len(seq)-1,-1,-1)
            ind_score=range(len(score)-1,-1,-1)
            ind_d=range(len(decision)-1,-1,-1)
        else: # if we search a C-rich telomere, we go through the sequence from the start
            ind_seq=range(len(seq))
            ind_score=range(len(score))
            ind_d=range(len(decision))
        while i < len(score): # if we reach the end of the read
            if telo == 0 :
                if score[ind_score[i]] in [str(np.base_repr(mer_size,36))] :
                    telo=1
                    start=ind_score[i]
                i+=1
            else : # if in telomeric stretch
                if i >= len(decision): # if we are close to the end of the sequence
                    if typ == 'G':
                        score_telo = score[:start+1]
                        while score_telo[0] not in  [str(np.base_repr(mer_size,36))] or np.mean([int(k,36) for k in score_telo]) < min_mean_telo: #str(mer_size-1),
                            score_telo=score_telo[1:]
                        # removes the first or/and the last nucleotide at the boundaries, when transitioning from a mer_size to mer_size-1 score
                        if len(score_telo)>1 and  score_telo[-1]==str(np.base_repr(mer_size-1,36)) and score_telo[-2]==str(np.base_repr(mer_size,36)):
                            start-=1
                            score_telo=score_telo[:-1]
                        end=start+mer_size
                        if len(score_telo)>1 and score_telo[0]==str(np.base_repr(mer_size-1,36))  and score_telo[1]==str(np.base_repr(mer_size,36)):
                            score_telo=score_telo[1:]
                        start=max(0,end - len(score_telo)- mer_size + 1)
                    else : 
                        score_telo = score[start:]
                        while score_telo[-1] not in [str(np.base_repr(mer_size-1,36)),str(np.base_repr(mer_size,36))] or np.mean([int(k,36) for k in score_telo]) < min_mean_telo :
                            score_telo=score_telo[:-1]
                        # removes the first or/and the last nucleotide at the boundaries, when transitioning from a mer_size to mer_size-1 score
                        if len(score_telo)>1 and score_telo[0] ==str(np.base_repr(mer_size-1,36)) and score_telo[1]==str(np.base_repr(mer_size,36)):
                            start+=1
                            score_telo=score_telo[1:]
                        end=min(start+len(score_telo)-1,len(score)-1)+mer_size
                        if len(score_telo)>1 and score_telo[-1]==str(np.base_repr(mer_size-1,36)) and score_telo[-2]==str(np.base_repr(mer_size,36)):
                            score_telo=score_telo[:-1]
                    if len(score_telo)+mer_size> min_len and score_telo.count(str(np.base_repr(mer_size,36)))>=mer_size:
                        dfTelo = Add_Telo(dfTelo,nam,start,end,typ,len(seq),np.mean([int(k,36) for k in score_telo]),strain,max_dist_term)
                    i=len(score) 
                elif decision[ind_d[i]] == 1: # if decision is 1, telomeric stretch continues
                    i+=1
                else : # enter in a new gap 
                    if start_gap != '': # if there is already a gap in the telomeric stretch
                        # Verifies if the previous gap is accepted, depending on the mean of telomeric stretch up to the current position before the new gap.
                        if typ == 'G':
                            if np.mean([int(k,36) for k in score[max(ind_score[i]-size_window+mer_size,0):start+1]]) < min_mean_telo :
                                telo=0
                        else :
                            if np.mean([int(k,36) for k in score[start:min(ind_score[i]+size_window-mer_size,len(score))+1]]) < min_mean_telo :
                                telo=0
                    if telo==1 : # if the previous gap is accepted, verifies if the new gap is tested
                        start_gap=ind_score[i]
                        # if the gap is too long, stops the telomere stretch.
                        if 1 not in decision[min(ind_d[i],ind_d[min(i+max_size_gap+mer_size,len(ind_d)-1)]):max(ind_d[i],ind_d[min(i+max_size_gap+mer_size,len(ind_d)-1)])+1] : 
                            telo=0
                        else: # telomeric stretch continues
                            while decision[ind_d[i]]!=1:
                                i+=1
                if telo == 0 :
                    if typ == 'G':
                        score_telo = score[max(0,start_gap-size_window+mer_size):start+1]
                        while score_telo[0] not in  [str(np.base_repr(mer_size-1,36)),str(np.base_repr(mer_size,36))] or np.mean([int(k,36) for k in score_telo]) < min_mean_telo: 
                            score_telo=score_telo[1:]
                        # removes the first or/and the last nucleotide at the boundaries, when transitioning from a 8 to 7 score 
                        if len(score_telo)>1 and  score_telo[-1]==str(np.base_repr(mer_size-1,36)) and score_telo[-2]==str(np.base_repr(mer_size,36)):
                            start-=1
                            score_telo=score_telo[:-1]
                        if len(score_telo)>1 and score_telo[0]==str(np.base_repr(mer_size-1,36)) and score_telo[1]==str(np.base_repr(mer_size,36)):
                            score_telo=score_telo[1:]
                        end=start+mer_size
                        start=max(0,end - len(score_telo)- mer_size + 1)
                        
                        i=ind_seq.index(max(start-1,0))

                    else : 
                        score_telo = score[start:min(start_gap+size_window-mer_size,len(score))+1]
                        j=min(len(score_telo)-(start_gap-start),len(score_telo))+1
                        while score_telo[-1] not in  [str(np.base_repr(mer_size-1,36)),str(np.base_repr(mer_size,36))] or np.mean([int(k,36) for k in score_telo]) < min_mean_telo :
                            score_telo=score_telo[:-1]


                        if len(score_telo)>1 and score_telo[0] ==str(np.base_repr(mer_size-1,36)) and score_telo[1]==str(np.base_repr(mer_size,36)):
                            start+=1
                            score_telo=score_telo[1:]
                        if len(score_telo)>1 and score_telo[-1]==str(np.base_repr(mer_size-1,36)) and score_telo[-2]==str(np.base_repr(mer_size,36)):
                            score_telo=score_telo[:-1]
                        end=min(start+len(score_telo),len(score)-1)+mer_size-1
                        i=ind_seq.index(end)
                    start_gap=''
                    # Check if there are at least three scores of 8 in the telomere
                    if len(score_telo)+mer_size> min_len and score_telo.count(str(np.base_repr(mer_size,36)))>=mer_size :
                        dfTelo = Add_Telo(dfTelo,nam,start,end,typ,len(seq),np.mean([int(k,36) for k in score_telo]),strain,max_dist_term)
    return(dfTelo)

def Add_N(dfTelo):
    '''
    Adds number of telomere by read in the dfTelo DataFrame

    :param dfTelo: DataFrame to add number of telomere by read
    :return: updated DataFrame with number of telomere by read
    '''
    for i in dfTelo.index:
        dfTelo.at[i,'N']=len(dfTelo[dfTelo['name']==dfTelo.at[i,'name']])
    return(dfTelo)

def Make_Fasta_Telo_Complet_dfTelo(csv,fasta,original_fasta):
    '''
    Creates a fasta file with all telomeric regions and adds telomeric lengths in csv. 

    :param csv: DataFrame (dfTelo) with all telomeric region found
    :param fasta: path of the new created fasta file with all telomeric regions
    :param original_fasta: path of original fasta file with all sequences
    :return :updated DataFrame with read lengths 
    '''
    foriginal=open(original_fasta,'r')
    f=open(fasta,'w')
    seq=''
    ok=0
    for l in foriginal:
        if l[0]=='>':
            if ok==1:
                for i in ind:
                    s=int(csv.at[i,'start'])-1
                    e=int(csv.at[i,'end'])
                    csv.at[i,'reads_len']=len(seq)
                    f.write('>'+nam+' Loc : '+csv.at[i,'Loc']+' Telo start: '+str(csv.at[i,'start'])+' Telo end: '+str(csv.at[i,'end'])+'\n')
                    f.write(seq[s:e]+'\n')
            nam=l[1:-1]
            if nam not in list(csv['name']):
                ok=0
            else:
                ok=1
                ind=list(csv[csv['name']==nam].index)
            seq=''
        elif ok==1:
            seq+=(l[:-1])
            
    if ok==1:
        for i in ind:
            s=int(csv.at[i,'start'])-1
            e=int(csv.at[i,'end'])
            csv.at[i,'reads_len']=len(seq)
            f.write('>'+nam+' Loc : '+csv.at[i,'Loc']+' Telo start: '+str(csv.at[i,'start'])+' Telo end: '+str(csv.at[i,'end'])+'\n')
            f.write(seq[s:e]+'\n')
    f.close()
    foriginal.close()
    return(csv)

def Find_Telo_on_Kmer(original_fasta,out_fasta,fileCKmer,fileGKmer,out_csv,mer_size,size_window,min_mean_window,min_mean_telo,strain,max_dist_term,max_size_gap,min_len):
    '''
    This function initializes dfTelo DataFrame, and executes the different functions to obtain output results.

    :param original_fasta: path of original fasta file with all sequences
    :param out_fasta: path of the new created fasta file with all telomeric regions
    :param fileCKmer: path of the file created with the function Fasta_to_Kmertype with C-rich Kmer scores of given sequences.
    :param fileGKmer: path of the file created with the function Fasta_to_Kmertype with G-rich Kmer scores of given sequences.
    :param out_csv: path of the out Dataframe summarizing all telomeric regions found.
    :param mer_size: size of telomere motif
    :param size_window: size of the sliding window for the calculation of the score average
    :param min_mean_window: minimum average score of a telomere window
    :param min_mean_telo: minimum average score of a telomere
    :param strain: strain / label name
    :param max_dist_term: maximal distance from the extremity to be considered as terminal
    :param max_size_gap: maximum size of a non telomeric gap in telomere
    :param min_len : minimum length of a telomere
    '''
    dfTelo=pd.DataFrame(columns=['strain','name','N','type','len','start','end','Loc','Score_Kmer','reads_len'])
    for file,typ in ((fileCKmer,'C'),(fileGKmer,'G')):
        f = open(file,'r')
        seq=''
        score=''
        for l in f :
            if l[0]=='>':
                if seq != '':
                    dfTelo=Eval_One_Read(dfTelo,nam,seq,typ,score,mer_size,size_window,min_mean_window,min_mean_telo,strain,max_dist_term,max_size_gap,min_len)
                nam = l[1:-1]
                seq=''
                score=''
                s=0
            elif s==1:
                score+=l.split(' ')[-1].strip('\n')
                s=0
            else : 
                seq+=l[:-1]
                s=1
        dfTelo=Eval_One_Read(dfTelo,nam,seq,typ,score,mer_size,size_window,min_mean_window,min_mean_telo,strain,max_dist_term,max_size_gap,min_len)
    dfTelo=Add_N(dfTelo)
    dfTelo.index.name = 'index'
    dfTelo.to_csv(out_csv,sep='\t',index=False)
    Make_Fasta_Telo_Complet_dfTelo(dfTelo,out_fasta,original_fasta)

def process_fasta(original_fasta, fileCKmer, fileGKmer, out_fasta, out_csv, mer_size,Tab_score,size_window,min_mean_window,min_mean_telo,strain,max_dist_term,max_size_gap,min_len):
    '''
    This function executes the different functions on a fasta file to obtain output results.

    :param original_fasta: path of original fasta file with all sequences
    :param fileCKmer: path of the file created with the function Fasta_to_Kmertype with C-rich Kmer scores of given sequences.
    :param fileGKmer: path of the file created with the function Fasta_to_Kmertype with G-rich Kmer scores of given sequences.
    :param out_fasta: path of the new created fasta file with all telomeric regions
    :param out_csv: path of the out Dataframe summarizing all telomeric regions found.
    :param mer_size: size of telomere motif
    :param Tab_score: table of K-mer score
    :param size_window: size of the sliding window for the calculation of the score average
    :param min_mean_window: minimum average score of a telomere window
    :param min_mean_telo: minimum average score of a telomere
    :param strain: strain / label name
    :param max_dist_term: maximal distance from the extremity to be considered as terminal
    :param max_size_gap: maximum size of a non telomeric gap in telomere
    :param min_len : minimum length of a telomere
    '''
    print('Run '+original_fasta)
    Fasta_to_Kmertype(original_fasta, fileCKmer, 'C', mer_size,Tab_score)
    Fasta_to_Kmertype(original_fasta, fileGKmer, 'G', mer_size,Tab_score)
    Find_Telo_on_Kmer(original_fasta, out_fasta, fileCKmer, fileGKmer, out_csv,mer_size,size_window,min_mean_window,min_mean_telo,strain,max_dist_term,max_size_gap,min_len)
    os.remove(original_fasta)
    os.remove(fileCKmer)
    os.remove(fileGKmer)
    print('End '+original_fasta)

if __name__ == "__main__":
    args = parse_arguments()

    mer_size = len(args.motif)

    strain=args.strain
    fastafile=args.fastafile.split('/')[-1]
    pwd=os.path.dirname(args.fastafile)+'/'
    motif=args.motif
    if motif == 'Yeast':
        mer_size=8
    else:
        mer_size = len(args.motif)
    threads=args.threads
    max_split_seq=args.max_split_seq
    size_window = args.size_window 
    max_size_gap = args.max_size_gap 

    if args.min_mean_window == 0:
        min_mean_window = mer_size-1
    else :
        min_mean_window = args.min_mean_window 

    if args.min_mean_telo == 0:
        min_mean_telo = mer_size*0.8
    else :
        min_mean_telo = args.min_mean_telo 

    max_dist_term = args.max_dist_term 
    min_len=args.min_len
    
    if args.out_pwd == 'None':
        out_pwd=pwd+'out_TeloReader/'
    else:
        out_pwd=args.out_pwd
        if len(out_pwd)==0:
            out_pwd='./'
        if out_pwd[-1]!='/':
            out_pwd+='/'

    print("strain : ",strain, 
    "\n fastafile : ",fastafile,
    "\n motif: ",motif,
    "\n mer size: ",mer_size,
    "\n threads: ",threads,
    "\n max_split_seq: ",max_split_seq,
    "\n size_window : ",size_window,
    "\n max_size_gap : ",max_size_gap,
    "\n min_mean_window : ",min_mean_window,
    "\n min_mean_telo : ",min_mean_telo,
    "\n max_dist_term : ",max_dist_term,
    "\n min_len : ",min_len,
    "\n out_pwd : ",out_pwd,'\n')


    os.makedirs(out_pwd, exist_ok=True)
    temp_files = out_pwd+"temp_files/"
    os.makedirs(temp_files, exist_ok=True)


    Tab_score=pd.read_csv(f"{Path(__file__).resolve().parent}/Motif_Score_Table/Score_for_{args.motif}.tab",sep='\t',index_col='Mer')

    nb_seq = int(subprocess.check_output(f"grep -c '^>' {os.path.join(pwd, fastafile)}", shell=True))
    if nb_seq <= max_split_seq:
        split_arg = f"-p {threads}"
    else:
        split_arg = f"-s {max_split_seq}"

    print(split_arg)
    split_location = out_pwd+'split/'
    os.makedirs(split_location, exist_ok=True)
    os.system(f"seqkit split2 {os.path.join(pwd, fastafile)} {split_arg} -O {split_location}")


    split_files = os.listdir(split_location)


    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = []
        for fasta_file in split_files :
            fileCKmer = f"{temp_files}{fasta_file.rsplit('.',1)[0]+'_C.kmer.txt'}"
            fileGKmer = f"{temp_files}{fasta_file.rsplit('.',1)[0]+'_G.kmer.txt'}"
            out_fasta = f"{temp_files}{fasta_file}_out.fasta"
            out_csv = f"{temp_files}{fasta_file.rsplit('.',1)[0]+'_out.csv'}"
            futures.append(executor.submit(process_fasta, split_location+fasta_file, fileCKmer, fileGKmer, out_fasta, out_csv, mer_size, Tab_score,size_window,min_mean_window,min_mean_telo,strain,max_dist_term,max_size_gap,min_len))
        # Wait for all tasks to finish
        concurrent.futures.wait(futures)


    # Now that all Teloreader processes have finished, merge the CSV files together

    output_csv = os.path.join(out_pwd, f"merged_output_{strain}_{os.path.splitext(fastafile)[0]}.csv")
    awk_command = f"awk '(NR == 1) || (FNR > 1)' {os.path.join(temp_files, '*.csv')} > {output_csv}"
    subprocess.run(awk_command, shell=True)
    for file in os.listdir(temp_files):
        if file.endswith(".csv"):
            os.remove(os.path.join(temp_files, file))

    output_fasta = os.path.join(out_pwd, f"merged_output_{strain}_{os.path.splitext(fastafile)[0]}.fasta")
    cat_command = f"cat {os.path.join(temp_files, '*.fasta')} >> {output_fasta}"
    subprocess.run(cat_command, shell=True)
    for file in os.listdir(temp_files):
        if file.endswith(".fasta"):
            os.remove(os.path.join(temp_files, file))

    print("All Teloreader processes have finished.")

    os.rmdir(split_location)
    os.rmdir(temp_files)

    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"Elapsed time for {args.fastafile}: {elapsed_time} seconds")
