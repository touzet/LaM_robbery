from dataclasses import dataclass
import sys

@dataclass
class Sequence:
    name: str = ""
    dna: str = ""

def parse_sequences(filename):
    """Lit un fichier de reads ou de genes et retourne une liste de Sequence."""
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
                    current_sequence.dna = ''.join(seq_lines)
                    sequences.append(current_sequence)
                header = line[1:]
                current_sequence = Sequence(name=header)
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            current_sequence.dna = ''.join(seq_lines)
            sequences.append(current_sequence)
    return sequences

def match_read(read, gene, start):
    i=0
    while  i<len(read) and read[i]==gene[i+start] :
        i=i+1
    return i==len(read)

def print_match(read, gene, start):
    read_length=len(read.dna)
    print(read.name+ " " +str(start)+"-"+str(start+read_length-1))
    print(read.dna)
    print('|'*read_length)
    print(gene[start:start+read_length])

def match_substitution(read, gene, start):
    '''
    La fonction ne gère que les cas ayant exactement une substitution
    '''
    i=0
    nb_substitutions=0
    pos_substitution=0
    while i<len(read) and nb_substitutions<2 :
        if read[i]!=gene[i+start]:
            nb_substitutions=nb_substitutions+1
            pos_substitution=i
        i=i+1
    if nb_substitutions==1:
        return True, pos_substitution
    else:
        return False, 0

def print_match_substitution(read,gene,start,pos):
    read_length=len(read.dna)
    print(read.name+ " "+str(start)+"-"+str(start+read_length-1))
    print(read.dna)
    print("|"*pos+"."+"|"*(read_length-pos-1))
    print(gene[start:start+read_length] )

def match_insertion(read,gene,start):
    '''
    La fonction ne gère que les cas ayant exactement une insertion
    '''
    i=0
    while i<len(read) and read[i]==gene[i+start] :
        i=i+1
    insertion = match_read(read[i+1:len(read)], gene, start+i)
    return insertion, i
    
def match_deletion(read,gene,start):
    '''
    La fonction ne gère que les cas ayant exactement une délétion
    '''
    i=0
    while i<len(read) and read[i]==gene[i+start] :
        i=i+1
    deletion = match_read(read[i:len(read)], gene, start+i+1)
    return deletion, i

def print_match_insertion(read, gene, start, pos):
    read_length=len(read.dna)
    print(read.name+" "+ str(start)+"-"+str(start+read_length-2))
    print(read.dna)
    print("|"*pos+" "+"|"*(read_length-pos-1))
    print(gene[start:start+pos]+"-"+gene[start+pos:start+read_length-1])

def print_match_deletion(read, gene, start, pos):
    read_length=len(read.dna)
    print(read.name+ " " + str(start)+"-"+str(start+read_length))
    print(read.dna[0:pos]+"-"+read.dna[pos:read_length])
    print("|"*pos+" "+"|"*(read_length-pos))
    print(gene[start:start+read_length+1]) 
          
        
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


def align(all_reads, all_genes):
    for gene in all_genes:
        print(gene.name)
        list_of_aligned_reads=[]
        gene_seq=gene.dna
        for read in all_reads:
            read_seq=read.dna
            found = False
            for start in range(0,len(gene_seq)-len(read_seq)+1):
                found = match_read(read_seq, gene_seq, start)
                if found:
                    print_match(read,gene_seq, start)
                    list_of_aligned_reads.append((len(read_seq),start))
                    break
            if not found:
                for start in range(0,len(gene_seq)-len(read_seq)+1):
                    found, pos = match_substitution(read_seq, gene_seq, start)
                    if found:
                        print_match_substitution(read,gene_seq, start, pos)
                        list_of_aligned_reads.append((len(read_seq),start))
                        break
            if not found:
                for start in range(0,len(gene_seq)-len(read_seq)+2):
                    found, pos = match_insertion(read_seq, gene_seq, start)
                    if found:
                        print_match_insertion(read,gene_seq, start, pos)
                        list_of_aligned_reads.append((len(read_seq)-1,start))
                        break
            if not found:
                for start in range(0,len(gene_seq)-len(read_seq)):
                    found, pos = match_deletion(read_seq, gene_seq, start)
                    if found:
                        print_match_deletion(read,gene_seq, start, pos)
                        list_of_aligned_reads.append((len(read_seq)+1,start))
                        break
        
        taux_de_couverture(list_of_aligned_reads, gene_seq)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} genes reads")
        sys.exit(1)
    all_genes = parse_sequences(sys.argv[1])
    all_reads = parse_sequences(sys.argv[2])
    align(all_reads, all_genes)
    
