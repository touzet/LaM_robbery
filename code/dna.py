from dataclasses import dataclass

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

def print_alignment(read,gene,start):
    read_seq = read.dna
    read_length = len(read_seq)
    gene_seq = gene.dna
    print(read.name+ " "+str(start)+"-"+str(start+read_length-1))
    print(read_seq)
    for i in range(0, read_length):
        if gene_seq[i+start]==read_seq[i]:
            print("|", end="")
        else:
            print(".", end="")
    print("")
    print(gene_seq[start:start+read_length] )

        
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




