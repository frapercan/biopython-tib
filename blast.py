from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

# Leer la secuencia de ADN desde un archivo FASTA
with open("SecuenciaADN.txt") as file:
    sequence = SeqIO.read(file, "fasta")

# Imprimir la secuencia de ADN
print(sequence.seq)

# Realizar una búsqueda BLAST en la base de datos refseq_rna utilizando blastn
result_handle = NCBIWWW.qblast("blastn", "refseq_rna", sequence.seq)

# Leer los resultados BLAST en formato XML
blast_record = NCBIXML.read(result_handle)


print(blast_record)


# Variables para almacenar el nombre del gen y si es homólogo
gene_name = None
is_homologous = False

# Recorrer los alineamientos obtenidos en los resultados BLAST
for alignment in blast_record.alignments:
    # Recorrer los segmentos de alto puntaje (HSP) dentro de cada alineamiento
    for hsp in alignment.hsps:
        # Verificar si el alineamiento pertenece a un gen con acceso que comienza con "NM"
        if alignment.accession.startswith("NM"):
            # Extraer el nombre del gen de la descripción del alineamiento
            gene_name = alignment.title.split("[")[0].strip()
            is_homologous = True
            break
    if is_homologous:
        break

# Imprimir el nombre del gen encontrado y si es homólogo o no
print("Nombre del gen encontrado:", gene_name)
if is_homologous:
    print("El gen puede considerarse homólogo.")
else:
    print("No se encontraron genes homólogos.")
