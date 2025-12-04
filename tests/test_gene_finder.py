import os
import filecmp
import pytest
from src import gene_finder

test_fasta = "tests/test_assets/test_genome.fasta"
test_gff = "tests/test_assets/test_genes.gff"
expected_fna = "tests/test_assets/expected_genes.fna"
out_fna = "tests/test_assets/out_genes.fna"

def test_extract_genes_basic():
    # Ejecuta el script principal con los archivos de prueba
    gene_finder.write_fasta(
        gene_finder.extract_genes(
            gene_finder.parse_fasta(test_fasta),
            gene_finder.parse_gff(test_gff)
        ),
        out_fna
    )
    assert filecmp.cmp(out_fna, expected_fna, shallow=False)
    os.remove(out_fna)

def test_min_length():
    # Solo debe extraer araC si min_length=24
    gene_finder.write_fasta(
        gene_finder.extract_genes(
            gene_finder.parse_fasta(test_fasta),
            gene_finder.parse_gff(test_gff, min_length=24)
        ),
        out_fna
    )
    with open(out_fna) as f:
        content = f.read()
    assert ">araC" in content
    assert ">crp" not in content
    os.remove(out_fna)

def test_reverse_complement():
    # Verifica que la función reverse_complement funciona correctamente
    seq = "ATGC"
    rc = gene_finder.reverse_complement(seq)
    assert rc == "GCAT"
