# "A","T","C","G","N"

import pyfaidx
import pysam
from operator import add
from pathlib import Path
import argparse


def make_mutation_matrix(chrom,tbx_path,genome_path):
    tbx = pysam.TabixFile(tbx_path)
    genome = pyfaidx.Fasta(genome_path)

    summ_type = {'A': [0,0,0,0,0],
             'T': [0,0,0,0,0],
             'C': [0,0,0,0,0],
             'G': [0,0,0,0,0],
             'N': [0,0,0,0,0]}

    i = 0
    for row in tbx.fetch(chrom, parser=pysam.asTuple()):
        i = i+1
        v = row[1]
        k = genome[chrom][int(v)-1]

        rowcounts = row[2:7]
        rowcounts = [int(v) for v in rowcounts]

        if i % 1_000_000 == 0:
            print(f'{i} reads processed in chromosome {chrom}')

        summ_type[str(k)] = list( map(add, summ_type[str(k)], rowcounts) )

    return(summ_type)

def write_mutation_matrix(chrom,summ_type,output,tbx):
    base_name_tbx  = Path(tbx).stem
    # Open a file to write
    file_name = f"{output}{base_name_tbx}_{chrom}_mutation_matrix.tsv"

    with open(file_name, 'w') as f:
        # Write the header
        header = "\t" + "\t".join(summ_type.keys())
        f.write(header + "\n")

        # Write each row
        for key in summ_type:
            row = [str(key)] + [str(v) for v in summ_type[key]] 

            f.write("\t".join(row) + "\n")
        # print(f"{str(k)}", row[nt_to_idx[str(k)]])

    print(f"chrom {chrom} written to {file_name}")

def main():

    #     genome = pyfaidx.Fasta("/SAN/vyplab/vyplab_reference_genomes/sequence/human/gencode/GRCh38.primary_assembly.genome.fa")

    #/SAN/vyplab/alb_projects/data/dartseq/BULLSEYE/mutation_parse/
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tbx")
    parser.add_argument("-f", "--fasta")
    parser.add_argument("-o", "--outputdir")

    args = parser.parse_args()

    tabix_file = args.tbx
    fasta_path = args.fasta
    output = args.outputdir

    #TODO this shouldn't be hardcoded - add in a function to check if the bam used chr or no
    STANDARD_CHROMSOMES = ["chr" + str(i) for i in range(1, 23)] + ['chrX', 'chrY']

    for c in STANDARD_CHROMSOMES:
        print(f"processing chrom: {c}")
        mtrx = make_mutation_matrix(c,tabix_file,fasta_path)
        print(mtrx)
        write_mutation_matrix(c,mtrx,output,tabix_file)

if __name__ == "__main__":
    main()