from Bio import Entrez
from Bio import SeqIO

# Establecer la dirección de correo electrónico para identificarse en las solicitudes a Entrez
Entrez.email = "tu_correo_electronico@example.com"

# Realizar una búsqueda en la base de datos "nucleotide" utilizando un término.
handle = Entrez.esearch(db="nucleotide", term="lacy")
record = Entrez.read(handle)
handle.close()

# Obtener una lista de los identificadores de secuencia obtenidos en la búsqueda
id_list = record["IdList"]


print('Busqueda sobre operon lactosa en la base de nucleotidos')
# Recorrer los identificadores y obtener la información de cada secuencia
for seq_id in id_list:
    # Obtener la secuencia en formato FASTA utilizando el identificador de secuencia
    handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
    sequence = SeqIO.read(handle, "fasta")
    handle.close()

    # Imprimir información sobre la secuencia
    print("ID:", sequence.id)
    print("Secuencia:", sequence.seq[:5],"...")
    print("Descripción:", sequence.description)
    print("-----------")




# Obtener información sobre las bases de datos disponibles en Entrez
handle = Entrez.einfo()
record = Entrez.read(handle)
handle.close()

# Recorrer la lista de bases de datos y mostrar sus nombres
print('#####################################')
print('Bases de datos disponibles en Entrez:')
for db in record["DbList"]:
    print(db)
