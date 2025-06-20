import argparse
import os

def parse_fasta(filename):
    """Parse un fichier FASTA et retourne un dictionnaire {header: sequence}."""
    sequences = {}
    with open(filename, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequences[header] = "".join(seq_lines)
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            sequences[header] = "".join(seq_lines)
    return sequences


def latex_header(output_filehandle):
    latex_content = """\\documentclass[a4paper]{article}
\\usepackage[utf8]{inputenc}
\\usepackage{xcolor}
\\usepackage{tikz}
\\usetikzlibrary{chains,matrix}
\\usepackage[margin=3mm,landscape]{geometry}

% Définition des couleurs pour chaque base
\\definecolor{Acolor}{HTML}{ACD700}
\\definecolor{Tcolor}{HTML}{C90075}
\\definecolor{Ccolor}{HTML}{317CC5}
\\definecolor{Gcolor}{HTML}{FF9F0E}

% Commande pour la taille de police
\\newcommand{\\dnafont}{\\fontsize{100pt}{70pt}\\selectfont\\sf}

% Commande de coloration
\\newcommand{\\highlightDNA}[2][]{%
    \\node[start chain=2 going right] {};
    \\foreach \\base in {#2} {%
        \\node[transform shape, on chain=2,fill=\\base color, text=white] {\\textbf{\\base}};
    }%
    \\end{tikzpicture}
}

\\begin{document}
%\\begin{center}
    \\dnafont
\\setlength{\\parskip}{-1.1em}
  \\begin{tikzpicture}
    \\matrix[matrix of nodes, column sep=0pt, row sep=4pt,
    nodes={text width=.9em,inner sep=8pt, outer sep=0pt, node distance=0pt, font=\\sf, scale=1.2, text=white}]
    {"""
    print(latex_content, file=output_filehandle,end="")

def latex_footer(output_filehandle):
    print("};\n\\end{tikzpicture}\n%\\end{center}\n\\end{document}", file=output_filehandle)

def save_split_genes_to_latex(split_genes, output_filehandle):
    """
    Enregistre les segments d'ADN dans un fichier LaTeX avec coloration des bases.
    """
    latex_content = ''
    for seq in split_genes:
        latex_content+= "\n"
        for i, c in enumerate(seq):
            latex_content += f"|[fill={c}color]| \\textbf{{{c}}} "
            if i < len(seq)-1:
                latex_content += "&"
        latex_content += "\\\\"

    print(latex_content, file=output_filehandle, end='')


def main():
    parser = argparse.ArgumentParser(
        description="Convertit un fichier FASTA en un fichier LaTeX avec coloration de séquence.")
    parser.add_argument("input", nargs='+', help="Fichier FASTA en entrée")
    parser.add_argument("-o", "--output", help="Fichier LaTeX en sortie", required=True)

    args = parser.parse_args()
    f = open(args.output, 'w')

    latex_header(f)

    for inputfile in args.input:
        name = os.path.splitext(os.path.basename(inputfile))[0].replace('_', '~')
        all_sequences = parse_fasta(inputfile)
        print(f"\n|[rotate=90,anchor=east,yshift=1.3em,xshift=-3em,overlay,text=black]| {name}\\\\[5pt]", file=f, end='')
        reads = [seq for header, seq in all_sequences.items()]
        #print("\\textbf{\Large "+name+"}\\\\", file=f)
        save_split_genes_to_latex(reads, f)
    latex_footer(f)
    f.close()


if __name__ == "__main__":
    main()
