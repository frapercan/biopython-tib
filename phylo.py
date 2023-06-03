from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Cargar el alineamiento múltiple en formato FASTA
alineamiento = AlignIO.read("alineamiento.fasta", "fasta")

# Calcular las distancias entre las secuencias
calculator = DistanceCalculator('identity')

# Construir el árbol filogenético.
constructor = DistanceTreeConstructor(calculator)
arbol = constructor.build_tree(alineamiento)

# Obtener un diccionario de las descripciones de las secuencias
descripcion_secuencias = {registro.id: registro.description for registro in alineamiento}

# Recorrer las hojas del árbol y cambiar las etiquetas por las descripciones
for clade in arbol.get_terminals():
    etiqueta = clade.name
    nueva_etiqueta = descripcion_secuencias.get(etiqueta, etiqueta)
    clade.name = nueva_etiqueta

# Dibujar y mostrar el árbol con las nuevas etiquetas
Phylo.draw(arbol)