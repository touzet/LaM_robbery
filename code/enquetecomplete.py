from dataclasses import dataclass
import sys
import dna

def align_position(read, gene, start):
    i=0
    while  i<len(read) and read[i]==gene[i+start] :
        i=i+1
    return i==len(read)

def align(read, gene):
    for position in range(0, len(gene)-len(read)+1):
        if substitution(read, gene, position):
            return position
    return -1
    
def substitution(read, gene, start):
    '''
    La fonction ne gère que les cas ayant exactement zero ou une substitution
    '''
    i=0
    nb_substitutions=0
    while i<len(read) and nb_substitutions<2 :
        if read[i]!=gene[i+start]:
            nb_substitutions=nb_substitutions+1
        i=i+1
    return  nb_substitutions<2

def align_all(all_reads, all_genes):
    for gene in all_genes:
        print(gene.name)
        list_of_aligned_reads=[]
        gene_seq=gene.dna
        for read in all_reads:
            read_seq=read.dna
            start=align(read_seq, gene_seq)
            if start >-1:
                dna.print_alignment(read, gene, start)
                list_of_aligned_reads.append((len(read_seq),start))
        print("\nNombre de reads :", len(list_of_aligned_reads))
        print("---------------------")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} genes reads")
        sys.exit(1)
    all_genes = dna.parse_sequences(sys.argv[1])
    all_reads = dna.parse_sequences(sys.argv[2])
    align_all(all_reads, all_genes)
    
