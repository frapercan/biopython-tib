from Bio import SeqIO
from Bio import SeqRecord
from Bio.HMM import MarkovModel
from Bio.HMM import DynamicProgramming
from Bio.Seq import Seq
import random

alphabet = ["A","C","G","T"]
emmision_alphabet = ["H1","H2"]
hmm = MarkovModel.MarkovModelBuilder(alphabet,emmision_alphabet)

# Configura las probabilidades del modelo
hmm.set_random_initial_probabilities()
hmm.set_random_emission_probabilities()
hmm.set_random_initial_probabilities()

# Genera una secuencia de aminoácidos
seq_length = 10
random_sequence = Seq("".join(random.choice(alphabet) for _ in range(seq_length)))

# Print the random sequence
print("Random sequence:", random_sequence)

# Convierte la secuencia en un objeto SeqRecord
seq_record = SeqRecord.SeqRecord(random_sequence)


