import logging
import os
import sys
from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Blast import NCBIWWW
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Configure logging to stdout
logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

class SequenceAnalysis:
    def __init__(self, email):
        self.email = email
        Entrez.email = self.email

    def search_sequence(self, term, retmax=5):
        handle = Entrez.esearch(db="nucleotide", term=term, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]
        logging.info(f"Found {len(id_list)} sequence IDs.")
        return id_list

    def fetch_sequences(self, seqs_id):
        # Fetch and print sequences
        sequences = []
        for i,seq_id in enumerate(seqs_id):
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
            sequence = SeqIO.read(handle, "fasta")
            sequences.append(sequence)
            handle.close()
        return sequences

    def choose_sequence(self, sequences):
        # Print the available sequences
        for i, sequence in enumerate(sequences):
            print(f"[{i + 1}] Sequence ID: {sequence.id}")
            print(f"    Description: {sequence.description}")
            print(f"    Sequence: {sequence.seq[:10]}\n")

        # Prompt the user to choose a sequence
        choice = input("Enter the number of the sequence you want to select: ")

        # Validate the user's input
        while not choice.isdigit() or int(choice) < 1 or int(choice) > len(sequences):
            print("Invalid choice. Please enter a valid number.")
            choice = input("Enter the number of the sequence you want to select: ")

        # Convert the user's choice to an index
        index = int(choice) - 1

        # Return the selected sequence
        return sequences[index]

    def blast_sequence(self, sequence):
        logging.info("Performing BLAST...")
        result_handle = NCBIWWW.qblast("blastn", "nr", sequence.seq)
        logging.info("BLAST Finished...")
        file_path = os.path.join(f"./results/blast/{sequence.id}.xml")
        print(file_path)
        with open(file_path, "w") as output_file:
            output_file.write(result_handle.read())

        logging.info(f"BLAST result saved as: {output_file}")

        result_handle.close()


    def save_sequence_to_fasta(self, sequence):
        logging.info("Saving sequences to FASTA file...")
        outfile = f"./results/searchz/f{sequence.id}.fasta"

        with open(outfile, "w") as file:
            SeqIO.write(sequence, file, "fasta")

    def perform_alignment(self, sequences_dir):
        outfile = f"./results/aligned/f{sequences_dir}"
        logging.info("Performing sequence alignment...")
        clustalomega_cline = ClustalOmegaCommandline(infile=f"./samples/multiple/{sequences_dir}", outfile=outfile, verbose=True, auto=True, force=True)
        stdout, stderr = clustalomega_cline()
        print(stdout)
        alineamiento = AlignIO.read(outfile, "fasta")
        logging.info("Alignment completed.")
        return alineamiento

    def generate_phylogenetic_tree(self, alignment):
        # Calcular las distancias entre las secuencias
        calculator = DistanceCalculator('identity')

        # Construir el árbol filogenético.
        constructor = DistanceTreeConstructor(calculator)
        arbol = constructor.build_tree(alignment)

        # Obtener un diccionario de las descripciones de las secuencias
        descripcion_secuencias = {registro.id: registro.description for registro in alignment}

        # Recorrer las hojas del árbol y cambiar las etiquetas por las descripciones
        for clade in arbol.get_terminals():
            etiqueta = clade.name
            nueva_etiqueta = descripcion_secuencias.get(etiqueta, etiqueta)
            clade.name = nueva_etiqueta

        # Dibujar y mostrar el árbol con las nuevas etiquetas
        Phylo.draw(arbol)



# Uso del programa
if __name__ == "__main__":
    email = "franciscoperez@gmail.com"
    analyzer = SequenceAnalysis(email)

    # Búsqueda en Entrez
    term = 'ABCB1[All Fields] AND "Homo sapiens"[Organism]'
    output_folder = "results/"
    seq_ids = analyzer.search_sequence(term, retmax=5)
    seqs = analyzer.fetch_sequences(seq_ids)
    seq = analyzer.choose_sequence(seqs)
    analyzer.save_sequence_to_fasta(seq)

    # # # # BLAST
    blast_sequences = analyzer.blast_sequence(seq)
    #
    #
    # # # Realizar alinemiento múltiple de Secuencias
    input_filename = "ABCB1_refseq_transcript.fasta"
    alignment = analyzer.perform_alignment(input_filename)
    # # # Generar árbol filogenético
    tree = analyzer.generate_phylogenetic_tree(alignment)

