# Ejercicio 1 — Extraer secuencias de genes (GFF + FASTA)

## **Parte A — Implementación base**

Escribe un programa llamado `extract_genes.py` que:

### **Uso**
```
python extract_genes.py --gff genes.gff --fasta genome.fasta --output genes.fna
```

### **Objetivo**
- Leer un archivo FASTA con el genoma completo.
- Leer un archivo GFF.
- Extraer las secuencias DNA correspondientes a las features `gene`.
- Guardarlas en un archivo FASTA de salida.

### **Ejemplo de salida esperada (genes.fna)**

```
>araC  gene_coords=3456-41020 strand=+
ATGCGTAGCTAGCTAGCTAGCTAA
>crp  gene_coords=3456-41020 strand=-
ATTTGCGCGGCGCGCGTTAG
```

### **Parte B — Extensión**
Agregar:

```
--min-length N
```

Ejemplo:

```
python extract_genes.py --gff genes.gff --fasta genome.fasta --output genes.fna --min-length 300
```

Solo genes con longitud ≥ 300 serán exportados.

---

## **Requisitos técnicos**
- `argparse` obligatorio.
- Funciones: `load_fasta()`, `parse_gff()`, `extract_gene_seqs()`.
- Manejo de errores con excepciones.
- Docstrings y PEP8.
- Pruebas con asserts + documento de pruebas, bien con pytest.
