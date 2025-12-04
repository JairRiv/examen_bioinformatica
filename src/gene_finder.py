# Este programa toma un archivo FASTA y un archivo GFF, extrae las secuencias de genes anotados en el GFF y las guarda en un archivo FASTA de salida.
#
# Uso:
# python gene_finder.py --gff genes.gff --fasta genome.fasta --output genes.fna --min-length 300
# python gene_finder.py --gff genes.gff --fasta genome.fasta --output genes.fna

import argparse
import os

def parse_fasta(fasta_path):
    """
    Lee un archivo FASTA y retorna la secuencia genómica como un string.
    Omite las líneas de encabezado (que comienzan con '>').
    Args:
        fasta_path (str): Ruta al archivo FASTA.
    Returns:
        str: Secuencia genómica concatenada.
    """
    seq = ""
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq


def parse_gff(gff_path, min_length=None):
    """
    Lee un archivo GFF y retorna una lista de diccionarios con información de los genes anotados.
    Solo considera las líneas con la característica 'gene'.
    Args:
        gff_path (str): Ruta al archivo GFF.
        min_length (int, opcional): Longitud mínima de los genes a extraer.
    Returns:
        list: Lista de diccionarios con id, start, end, strand y length de cada gen.
    """
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature.lower() != "gene":
                continue
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]
            length = end - start + 1
            if min_length and length < min_length:
                continue
            gene_id = "unknown"
            for attr in attributes.split(";"):
                if attr.startswith("ID="):
                    gene_id = attr.split("=")[1]
                    break
            genes.append({
                "id": gene_id,
                "start": start,
                "end": end,
                "strand": strand,
                "length": length
            })
    return genes


def reverse_complement(seq):
    """
    Obtiene la secuencia reversa complementaria de una secuencia de ADN.
    Args:
        seq (str): Secuencia de ADN.
    Returns:
        str: Secuencia reversa complementaria.
    """
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]


def extract_genes(seq, genes):
    """
    Extrae las secuencias de los genes anotados según las coordenadas y el strand.
    Args:
        seq (str): Secuencia genómica.
        genes (list): Lista de diccionarios con información de los genes.
    Returns:
        list: Lista de diccionarios con id, secuencia y longitud de cada gen extraído.
    """
    extracted = []
    for gene in genes:
        gene_seq = seq[gene["start"]-1:gene["end"]]
        if gene["strand"] == "-":
            gene_seq = reverse_complement(gene_seq)
        extracted.append({
            "id": gene["id"],
            "seq": gene_seq,
            "length": gene["length"]
        })
    return extracted


def write_fasta(genes, output_path):
    """
    Guarda las secuencias extraídas en un archivo en formato FASTA.
    Args:
        genes (list): Lista de diccionarios con id, secuencia y longitud de cada gen.
        output_path (str): Ruta al archivo de salida FASTA.
    """
    with open(output_path, "w") as f:
        for gene in genes:
            f.write(f">{gene['id']} length={gene['length']}\n")
            # Escribir en líneas de 60 caracteres
            for i in range(0, len(gene['seq']), 60):
                f.write(gene['seq'][i:i+60] + "\n")


def main():
    """
    Función principal del programa. Procesa los argumentos, valida archivos, extrae genes y guarda el resultado.
    """
    parser = argparse.ArgumentParser(
        description="Extrae las secuencias de genes anotados en un archivo GFF desde una secuencia genómica en FASTA y las guarda en un archivo FASTA de salida."
    )
    parser.add_argument("--gff", required=True, help="Archivo GFF con anotaciones de genes.")
    parser.add_argument("--fasta", required=True, help="Archivo FASTA con la secuencia genómica.")
    parser.add_argument("--output", required=True, help="Archivo de salida en formato FASTA.")
    parser.add_argument("--min-length", type=int, default=None, help="Longitud mínima de los genes a extraer.")
    args = parser.parse_args()

    # Validar archivos de entrada
    for path in [args.gff, args.fasta]:
        if not os.path.isfile(path):
            print(f"Error: El archivo {path} no existe.")
            exit(1)
    try:
        seq = parse_fasta(args.fasta)
        genes = parse_gff(args.gff, args.min_length)
        if not genes:
            print("No se encontraron genes que cumplan los criterios.")
            exit(1)
        extracted = extract_genes(seq, genes)
        write_fasta(extracted, args.output)
        print(f"Se extrajeron {len(extracted)} genes en {args.output}")
    except Exception as e:
        print(f"Error: {e}")
        exit(1)

if __name__ == "__main__":
    main()

