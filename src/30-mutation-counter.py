import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd # this library is used for advanced data management
import pysam # this is the python version of samtools
import seaborn as sns # this is used for improving plots
import argparse # this library is used for parsing multidimensional complex data
import tqdm # this library is used to monitor progress of the job
import glob
from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO
import sys


parser = argparse.ArgumentParser(description='Finds nucleotide changes from input bam file.')
parser.add_argument('--bam', help='Input bam file path', required=True)
parser.add_argument('--fasta', help='Input reference fasta file path', required=True)
parser.add_argument('--contig', help='Input reference contig', required=True)
parser.add_argument('--start-codon-pos',
    help='First nucleotide of the coding sequence', type=int, required=True)
parser.add_argument('--stop-codon-pos',
    help='First nucleotide after the coding sequence', type=int,required=True)
parser.add_argument('--base-quality-filter',
    help='Minimum accepted phred score for each mutation position [Default=30]', type=int, default=30)
parser.add_argument('--output', help='Specify the output folder', required=True)

args = vars(parser.parse_args())

bam = pysam.AlignmentFile(args["bam"], "rb")
ref = pysam.FastaFile(args["fasta"])
# fetch the reference seq and convert it to uppercase, sequences are case sensitive
ref_seq = ref.fetch(args["contig"]).upper() 
# create a mutations list
mutations=[]
filtered_mutation_count = 0
# start fetching reads from the bam file 
for read in tqdm.tqdm(bam.fetch(args["contig"]), mininterval=30):
    # get all of the reference positions for the read, exact reference positions used
    ref_positions = read.get_reference_positions(full_length = True)
    # in read, go over every base
    for base_pos, query_base, base_quality in zip(ref_positions, read.query_sequence, read.query_qualities):
        # if base is 'none', skip to the next base
        if base_pos is None:
            continue
        # if there is a mutation
        if ref_seq[base_pos] != query_base:
            # if base quality is too low filter out
            if base_quality < args["base_quality_filter"]:
                filtered_mutation_count+=1
                continue
            # pull info for the mutation
            mutation_record = {'reference_pos': base_pos,
                               'reference_base': ref_seq[base_pos],
                               'mutated_base':query_base,
                               'base_quality': base_quality,
                               'read_id': read.query_name}
            # append the mutation to the mutations list
            mutations.append(mutation_record)
# organize the file by creating a data frame
df = pd.DataFrame(mutations)
print("[1/3] Mutations table created.")
print(f"Number of mutations found:\t{df.shape[0]}")
print(f"Percent low quality filtered:\t%{100*filtered_mutation_count/df.shape[0]}")

# calculcate reference amino acid positions wrt the input start codon position
df['reference_aa_pos'] = (df['reference_pos'] - args['start_codon_pos'])//3
# calculate reference codon for each mutation record
ref_codons=[]
obs_codons=[]
for i in df.index:
    p = df.loc[i, 'reference_pos']
    codon_start = p - (p-args['start_codon_pos'])%3
    ref_codons.append("".join(ref_seq[codon_start:codon_start+3]).upper())
# insert reference codon and reference aa fiels into dataframe
df['reference_codon'] = ref_codons
# translate reference codons
df['reference_aa'] = df['reference_codon'].apply(lambda codon: Seq.translate(codon)[0])
# index mutations dataframe by read id and reference aa positions to calculate
# observed codons that are modified more than one time
df_indexed = df.set_index(['read_id', 'reference_aa_pos'])

print("[2/3] Reference codons found and indexed.")

# this kernel func takes dataframe rows indexed by read id and reference_aa_pos
# and records observed codons into the indexed mutation dataframe
def calculate_observed_codon(row):
    ref_pos_frame = (np.array(row['reference_pos'])-args['start_codon_pos'])%3
    observed_codon = list(list(row['reference_codon'].values)[0])
    for mutation_id in range(len(ref_pos_frame)):
        observed_codon[ref_pos_frame[mutation_id]] = row['mutated_base'].values[mutation_id]
    df_indexed.loc[(row['read_id'], row['reference_aa_pos']), 'observed_codon'] = "".join(observed_codon)

# run the above function over the mutation dataframe
df.groupby(['read_id','reference_aa_pos']).apply(calculate_observed_codon)
# translate observed codons
df_indexed['observed_aa'] = df_indexed['observed_codon'].apply(lambda codon: Seq.translate(codon)[0])
# sort the dataframe based on reference position
df_indexed = df_indexed.sort_values(by=['reference_pos', 'mutated_base'])
# write the output to a tsv file
df_indexed.to_csv(args["output"]+"_all_mutations.tsv", sep='\t') 

print("[3/3] Observed codons calculated and printed.")
