from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline

# Leer las secuencias desde un archivo FASTA
secuencias = SeqIO.parse("ABCB1_refseq_transcript.fasta", "fasta")

# Guardar las secuencias en un archivo temporal
temp_file = "temp.fasta"
SeqIO.write(secuencias, temp_file, "fasta")

# Ejecutar Clustal Omega para realizar el alineamiento
clustalomega_cline = ClustalOmegaCommandline(infile=temp_file, outfile="alineamiento.fasta", verbose=True, auto=True)
stdout, stderr = clustalomega_cline()

# Leer el alineamiento desde el archivo de salida
alineamiento = AlignIO.read("alineamiento.fasta", "fasta")

# Imprimir el alineamiento
print(alineamiento)
