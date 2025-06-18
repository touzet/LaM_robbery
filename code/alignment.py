from dataclasses import dataclass
import sys

@dataclass
class Sequence:
    name: str = ""
    sequence: str = ""

def parse_fasta(filename):
    """Parse un fichier FASTA et retourne une liste de Sequence."""
    sequences = []
    current_sequence = Sequence()
    with open(filename, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    current_sequence.sequence = ''.join(seq_lines)
                    sequences.append(current_sequence)
                header = line[1:]
                current_sequence = Sequence(name=header)
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            current_sequence.sequence = ''.join(seq_lines)
            sequences.append(current_sequence)
    return sequences

def match(read, gene, start):
    i=0
    while  i<len(read) and read[i]==gene[i+start] :
        i=i+1
    return i==len(read)

def print_match(read, gene, start):
    print(read)
    print('|'*len(read))
    print(gene[start:start+len(read)])


def match_substitution(read, gene, start):
    i=0
    nb_substitutions=0
    pos_substitution=0
    while i<len(read) and nb_substitutions<2 :
        if read[i]!=gene[i+start]:
            nb_substitutions=nb_substitutions+1
            pos_substitution=i
        i=i+1
    if nb_substitutions==0:
        return True, "match", 0
    elif nb_substitutions==1:
        return True, "substitution",  pos_substitution
    else:
        return False, "", 0


def print_match_substitution(read,gene,start,pos):
    print(read)
    print("|"*pos+"."+"|"*(len(read)-pos-1))
    print(gene[start:start+len(read)] )


def match_error(read,gene,start):
    i=0
    while i<len(read) and read[i]==gene[i+start] :
        i=i+1
    if i==len(read):
        return True, "match", 0
    # substitution
    substitution = match(read[i+1:len(read)], gene, start+i+1)
    if substitution:
        return True, "substitution", i
    # insertion in the read
    insertion = match(read[i+1:len(read)], gene, start+i)
    if insertion:
        return True, "insertion", i
    #deletion in the read
    if i>0 :
        deletion = match(read[i:len(read)], gene, start+i+1)
        if deletion:
            return True, "deletion", i
    return False, "", 0

def print_match_insertion(read, gene, start, pos):
    print(read)
    print("|"*pos+" "+"|"*(len(read)-pos-1))
    print(gene[start:start+pos]+"-"+gene[start+pos:start+len(read)-1])

def print_match_deletion(read, gene, start, pos):
    print(pos)
    print(read[0:pos]+"-"+read[pos:len(read)])
    print("|"*pos+" "+"|"*(len(read)-pos))
    print(gene[start:start+len(read)+1]) 
        
def print_match_error(read, gene, start, edit, pos):
    if edit=="match":
        print_match(read, gene, start)
    elif edit=="substitution":
        print_match_substitution(read, gene, start, pos)
    elif edit=="insertion":
        print_match_insertion(read, gene, start, pos)
    elif edit=="deletion":
        print_match_deletion(read, gene, start, pos)
        
def taux_de_couverture(list_of_aligned_reads, gene):
    list_of_aligned_reads.sort(key=lambda x:x[1])
    couverture=0
    gene_position=0
    for (length, start) in list_of_aligned_reads:
        end=start+length-1
        if start >= gene_position:
            couverture = couverture + length
            gene_position = end + 1
        elif end >= gene_position:
            couverture = couverture + end - gene_position + 1
            gene_position = end + 1
    print("Taux de couverture: "+str(couverture/len(gene))+"\n")


def align_match(all_reads, all_genes):
    print("SANS ERREURS, SANS SUBSTITUTIONS")
    for gene in all_genes:
        print(gene)
        list_of_aligned_reads=[]
        gene_seq=gene.sequence
        for read in all_reads:
            read_seq=read.sequence
            for start in range(0,len(gene_seq)-len(read_seq)+1):
                found = match(read_seq, gene_seq, start)
                if found:
                    print(read)
                    print_match(read_seq,gene_seq, start)
                    print()
                    list_of_aligned_reads.append((len(read_seq),start))
                    break
        taux_de_couverture(list_of_aligned_reads, gene_seq)
        print("--------------------")


def align_substitution(all_reads, all_genes):
    print("SUBSTITUTIONS AUTORISEES")
    for gene in all_genes:
        print(gene)
        list_of_aligned_reads=[]
        gene_seq=gene.sequence
        for read in all_reads:
            read_seq=read.sequence
            for start in range(0,len(gene_seq)-len(read_seq)+1):
                found, edit, pos = match_substitution(read_seq, gene_seq, start)
                if found:
                    print(read)
                    print_match_error(read_seq, gene_seq, start, edit, pos)
                    list_of_aligned_reads.append((len(read_seq),start))
                    print()
                    break
        taux_de_couverture(list_of_aligned_reads, gene_seq)
        print("--------------------")        
        
def align_error(all_reads, all_genes):
    print("TOUTES ERREURS AUTORISEES")
    for gene in all_genes:
        print(gene)
        list_of_aligned_reads=[]
        gene_seq=gene.sequence
        for read in all_reads:
            read_seq=read.sequence
            for start in range(0,len(gene_seq)-len(read_seq)+1):
                found, edit, pos = match_error(read_seq,gene_seq,start)
                if found:
                    print(read.name + " " + edit)
                    print_match_error(read_seq,gene_seq, start, edit, pos)
                    print()
                    if edit=="match" or edit=="substitution":
                        list_of_aligned_reads.append((len(read_seq),start))
                    elif edit=="insertion":
                        list_of_aligned_reads.append((len(read_seq)-1,start))
                    else: # deletion
                        list_of_aligned_reads.append((len(read_seq)+1,start))
                    break
        taux_de_couverture(list_of_aligned_reads, gene_seq)
        print("--------------------")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} genes reads")
        sys.exit(1)
    all_genes = parse_fasta(sys.argv[1])
    all_reads = parse_fasta(sys.argv[2])
    align_match(all_reads, all_genes)
    align_substitution(all_reads, all_genes)
    align_error(all_reads, all_genes)
