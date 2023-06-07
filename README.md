# Experiencia práctica del manejo de secuencias con BioPython

Exploramos diversas funcionalidades relacionadas con el análisis de secuencias biológicas. Algunas de las funciones
clave incluyen la búsqueda de secuencias en una base de datos, la recuperación de las secuencias encontradas, la
selección de una secuencia específica, el uso de BLAST para comparar una secuencia con una base de datos, el guardado de
secuencias en formato FASTA, el alineamiento de secuencias y la generación de un árbol filogenético. El código utiliza
la biblioteca Biopython y las APIs de Entrez y NCBI para acceder a bases de datos y ejecutar las operaciones
correspondientes. 


El ejemplo proporcionado consta de dos partes:
*  Se realiza una búsqueda de secuencias con un término específico, se
selecciona una secuencia de interés, se guarda en un archivo FASTA y se realiza una comparación BLAST con esa secuencia.

* Tras haber seleccionado de forma arbitraria un conjunto de secuencias sobre el gen ABCB1 y haberlos almacenado en un único fichero fasta, procedemos a utilizar clustalOmega para alinear estas secuencias y generar un arbol filogenético para estudiar el parentezco entre especies.

## Resultados
* Ficheros fasta de la busqueda en Entrez.
* Resultado en .XML de Blast.
* Resultado de la alineación multiple de secuencias.
* Arbol filogenético.

![](/home/xaxi/PycharmProjects/biopython-tib/results/trees/ABCB1_refseq_transcript.png)

# Logs de ejecución

```
/home/xaxi/PycharmProjects/biopython-tib/environment/bin/python /home/xaxi/PycharmProjects/biopython-tib/example.py 
2023-06-08 00:37:59,202 INFO Found 5 sequence IDs.
[1] Sequence ID: NG_034070.1
    Description: NG_034070.1 Homo sapiens calcium voltage-gated channel auxiliary subunit alpha2delta 2 (CACNA2D2), RefSeqGene on chromosome 3
    Sequence: TCCCTGTGCT

[2] Sequence ID: NM_001348989.2
    Description: NM_001348989.2 Homo sapiens ATP binding cassette subfamily G member 2 (Junior blood group) (ABCG2), transcript variant 7, mRNA
    Sequence: GTGACTGGGC

[3] Sequence ID: NM_001348945.2
    Description: NM_001348945.2 Homo sapiens ATP binding cassette subfamily B member 1 (ABCB1), transcript variant 1, mRNA
    Sequence: ATCATTACTC

[4] Sequence ID: NM_000927.5
    Description: NM_000927.5 Homo sapiens ATP binding cassette subfamily B member 1 (ABCB1), transcript variant 3, mRNA
    Sequence: ATCATTACTC

[5] Sequence ID: NM_001348944.2
    Description: NM_001348944.2 Homo sapiens ATP binding cassette subfamily B member 1 (ABCB1), transcript variant 2, mRNA
    Sequence: ATCATTACTC

Enter the number of the sequence you want to select: 3
2023-06-08 00:39:02,959 INFO Saving sequences to FASTA file...
2023-06-08 00:39:02,960 INFO Performing BLAST...
2023-06-08 00:40:07,002 INFO BLAST Finished...
./results/blast/NM_001348945.2.xml
2023-06-08 00:40:07,004 INFO BLAST result saved as: <_io.TextIOWrapper name='./results/blast/NM_001348945.2.xml' mode='w' encoding='UTF-8'>
2023-06-08 00:40:07,004 INFO Performing sequence alignment...
Using 12 threads
Read 6 sequences (type: DNA) from ./samples/multiple/ABCB1_refseq_transcript.fasta
not more sequences (6) than cluster-size (100), turn off mBed
Setting options automatically based on input sequence characteristics (might overwrite some of your options).
Auto settings: Enabling mBed.
Auto settings: Setting iteration to 1.
Using 5 seeds (chosen with constant stride from length sorted seqs) for mBed (from a total of 6 sequences)
Calculating pairwise ktuple-distances...
Ktuple-distance calculation progress: 0 % (0 out of 20)
Ktuple-distance calculation progress: 5 % (1 out of 20)
Ktuple-distance calculation progress: 10 % (2 out of 20)
Ktuple-distance calculation progress: 15 % (3 out of 20)
Ktuple-distance calculation progress: 20 % (4 out of 20)
Ktuple-distance calculation progress: 94 % (19 out of 20)
Ktuple-distance calculation progress: 114 % (23 out of 20)
Ktuple-distance calculation progress: 129 % (26 out of 20)
Ktuple-distance calculation progress: 139 % (28 out of 20)
Ktuple-distance calculation progress done. CPU time: 0.63u 0.00s 00:00:00.63 Elapsed: 00:00:00
mBed created 1 cluster/s (with a minimum of 1 and a soft maximum of 100 sequences each)
Distance calculation within sub-clusters: 0 % (0 out of 1)
Distance calculation within sub-clusters done. CPU time: 0.51u 0.00s 00:00:00.51 Elapsed: 00:00:00
Guide-tree computation (mBed) done.
Progressive alignment progress: 20 % (1 out of 5)
Progressive alignment progress: 40 % (2 out of 5)
Progressive alignment progress: 60 % (3 out of 5)
Progressive alignment progress: 80 % (4 out of 5)
Progressive alignment progress: 100 % (5 out of 5)
Progressive alignment progress done. CPU time: 6.82u 1.28s 00:00:08.10 Elapsed: 00:00:05
Iteration step 1 out of 1
Computing new guide tree (iteration step 0)
Calculating pairwise aligned identity distances...
Pairwise distance calculation progress: 0 % (0 out of 21)
Pairwise distance calculation progress: 0 % (0 out of 21)
Pairwise distance calculation progress: 0 % (0 out of 21)
Pairwise distance calculation progress: 4 % (1 out of 21)
Pairwise distance calculation progress: 4 % (1 out of 21)
Pairwise distance calculation progress: 85 % (18 out of 21)
Pairwise distance calculation progress: 109 % (23 out of 21)
Pairwise distance calculation progress: 123 % (26 out of 21)
Pairwise distance calculation progress: 133 % (28 out of 21)
Pairwise distance calculation progress: 138 % (29 out of 21)
Pairwise identity calculation progress done. CPU time: 0.01u 0.00s 00:00:00.01 Elapsed: 00:00:00
Guide-tree computation done.
Computing HMM from alignment
Progressive alignment progress: 20 % (1 out of 5)
Progressive alignment progress: 40 % (2 out of 5)
Progressive alignment progress: 60 % (3 out of 5)
Progressive alignment progress: 80 % (4 out of 5)
Progressive alignment progress: 100 % (5 out of 5)
Progressive alignment progress done. CPU time: 16.31u 3.61s 00:00:19.91 Elapsed: 00:00:17
Alignment written to ./results/aligned/fABCB1_refseq_transcript.fasta

2023-06-08 00:40:29,902 INFO Alignment completed.
```