#Getting ready:
```
wsl
git clone https://github.com/mjlopezramal/coding.git
cd coding
mkdir coding
touch coding/Session_1_exercise.md
nano coding/Session_1_exercise.md
docker pull csicunam/bioinformatics_iamz
cd coding
docker pull csicunam/bioinformatics_iamz
```
Set persistent folder for results files to avoid data loss when docker is turned off
For instance, for the VEP class you could create one named 'vep_class'
```mkdir $HOME/vep_class 
chmod a+w $HOME/vep_class
```
Launch docker binding the persistent folder to internal folder (/data)
```
sudo docker run -t -i -v $HOME/vep_class:/data -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY csicunam/bioinformatics_iamz:latest
```

# 4.1 Annotating coding sequences by alignment to other sequences

## 4.1.1 Formatting a sequence collection for BLAST
The first task is to format our test sequence set, which was obtained from UniProt:
Since I do not have write permissions in some directories, I work inside the mounted persistent folder (/data). From inside the container, I copy the file uniprot_Atha.fasta.gz into /data and decompress it there.
```
vep@487087b874ff:/home/vep$ ls
bioinformatics  environment.yml  get_homologues  test_data  variant_data
vep@487087b874ff:/home/vep$ cd test_data
vep@487087b874ff:/home/vep/test_data$ ls
vep@487087b874ff:/home/vep/test_data$ cd /data
vep@487087b874ff:/data$ ls
vep@487087b874ff:/data$ cp /home/vep/test_data/uniprot_Atha.fasta.gz .
vep@487087b874ff:/data$ gunzip uniprot_Atha.fasta.gz
```
**Exe1)** How many sequences have been formatted and how does this affect the E-value of BLAST searches?
Each sequence in a FASTA file starts with a line beginning with >. Therefore, counting these lines gives the number of sequences.
```
vep@487087b874ff:/data$ grep -c "^>" uniprot_Atha.fasta
15719
```
Result: 15,719 sequences

The E-value (expect value) measures the number of hits one can expect to see by chance when searching a database of a particular size.

A larger database increases the E-value for the same alignment score, making matches appear less statistically significant. Conversely, a smaller database reduces the E-value for the same score, making matches appear more significant.

Thus, having 15,719 sequences in the database affects the statistical significance of BLAST search results.

## 4.1.2 Querying the collection with a sample coding sequence
I encountered issues because BLAST was not installed inside Docker, while it was available in a Conda environment on my system. Therefore, I copied the necessary files from the Docker container to my local WSL environment.
```
(base) mjlopez@DESKTOP-MBGBGC4:/mnt/c/Users/mjlop$ conda activate blast_env
(blast_env) mjlopez@DESKTOP-MBGBGC4:/mnt/c/Users/mjlop$ which blastp

(blast_env) mjlopez@DESKTOP-MBGBGC4:/mnt/c/Users/mjlop$ cd /mnt/c/Users/mjlop/coding/coding
(blast_env) mjlopez@DESKTOP-MBGBGC4:/mnt/c/Users/mjlop/coding/coding$ docker ps
```
Copiar todo lo necesario a la carpeta WSL
```docker cp 487087b874ff:/home/vep/test_data/uniprot_Atha.fasta.gz /mnt/c/Users/mjlop/coding/coding/
docker cp 487087b874ff:/home/vep/test_data/test.faa /mnt/c/Users/mjlop/coding/coding/
docker cp 487087b874ff:/home/vep/test_data/test.fna /mnt/c/Users/mjlop/coding/coding/
```
Check the files:
```(blast_env) mjlopez@DESKTOP-MBGBGC4:/mnt/c/Users/mjlop/coding/coding$ ls -lh /mnt/c/Users/mjlop/coding/coding
(blast_env) mjlopez@DESKTOP-MBGBGC4:/mnt/c/Users/mjlop/coding/coding$ cd /mnt/c/Users/mjlop/coding/coding
```

At this point, I will work with protein ARF6. Its protein and transcript sequence is in files test.faa and test.fna, respectively. And I shall now look for similar sequences in my collection.
**Exe2)*** Can you redirect the output to separate files called test.faa.blast and test.fna.blast?
```
# Unzip FASTA if you din’t before
gunzip uniprot_Atha.fasta.gz

# Create BLAST data base (prefijo: uniprot_Atha_db)
makeblastdb -dbtype prot -in uniprot_Atha.fasta -out uniprot_Atha_db

# BLASTP (protein vs protein)
blastp -db uniprot_Atha_db -query test.faa -outfmt 6 > test.faa.blast

# BLASTX (transcript vs protein)
blastx -db uniprot_Atha_db -query test.fna -outfmt 6 > test.fna.blast
```

**Exe3)** What is the default alignment format, can you show an example?
When using the `-outfmt 6` option in BLAST, the output is in tabular format (tab-delimited), also called BLAST tabular format.
Each line represents a hit (a sequence similar to the query), with columns separated by tabs.

**Columns in format 6:**

1. `qseqid` → Query sequence ID  
2. `sseqid` → Subject sequence ID (database)  
3. `pident` → Percentage of identical matches  
4. `length` → Alignment length  
5. `mismatch` → Number of mismatches  
6. `gapopen` → Number of gap openings  
7. `qstart` → Start of alignment in query  
8. `qend` → End of alignment in query  
9. `sstart` → Start of alignment in subject  
10. `send` → End of alignment in subject  
11. `evalue` → Expect value of the hit  
12. `bitscore` → Bit score of the alignment  

This is one example from my data: 
AT1G30330.2     sp|Q9ZTX8|ARFF_ARATH    100.000 935     0       0       1  935      1       935     0.0     1915

**Exe4**) Are there differences in the results retrieved in both searches?
- `blastp` compares the protein sequence (`test.faa`) directly with the protein database.  
- `blastx` translates the nucleotide sequence (`test.fna`) in all 6 reading frames and compares with the protein database.

Observed differences:
1. Number of hits  
   - `blastp` usually finds more precise and complete hits  
   - `blastx` may find fewer or partial hits depending on reading frame  
2. E-values and bit scores  
   - Hits from `blastp` typically have lower E-values (more significant)  
   - Hits from `blastx` may have higher E-values or fragmentary alignments  
3. Alignment positions  
   - `blastx` may align only fragments of the protein  
   - `blastp` usually aligns the full protein  

## 4.1.3 Producing a sequence profile
Here I’m making a sequence profile out of similar sequences matched with three iterations of BLASTP, using PSI-BLAST:
1.	Make sure I’m in the propor directory:
```
cd /mnt/c/Users/mjlop/coding/coding
conda activate blast_env
```
2.	Create the database if you haven’t done before:
```
makeblastdb -dbtype prot -in uniprot_Atha.fasta -out uniprot_Atha_db
```
3.	Run PSI-BLAST with three iterations and generate a PSSM output:
```
psiblast -db uniprot_Atha_db -query test.faa -num_iterations 3 -out_ascii_pssm profile.out
```
4.	Verify that the file has been generated:
```
ls -lh profile.out
head profile.out
```
The file `profile.out` is generated by PSI-BLAST when using the `-out_ascii_pssm` option.  It contains the Position-Specific Scoring Matrix (PSSM) in ASCII format, representing how conserved each amino acid position is across the sequences found in the BLAST search.

**Exe5)** Can you explain the contents of the output file `profile.out`?
1. **Header information**  
   - Query sequence name and description  
   - Number of residues, number of sequences used in the profile, and BLAST iteration number  
2. **Matrix data (PSSM)**  
   - Each row corresponds to a position in the query protein  
   - Each column corresponds to one of the 20 amino acids  
   - The values indicate **log-odds scores**, reflecting the likelihood of observing each amino acid at that position compared to random expectation  
3. **Additional statistics**  
   - Relative weight of the sequences  
   - Conservation scores per position  
   - Alignment information for each iteration  

**Example snippet from `profile.out`:
Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts
            A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    1 M    -1  -1  -2  -3  -2  -1  -2  -3  -2   1   2  -1   6   0  -3  -2  -1  -2  -1   1    0   0   0   0   0   0   0   0   0   0   0   0 100   0   0   0   0   0   0   0  0.39 0.01
    2 R    -1   3   0  -1  -3   1   1  -2  -1  -3  -2   4  -1  -3  -1   0  -1  -3  -2  -2    0  25   0   0   0   0   0   0   0   0   0  75   0   0   0   0   0   0   0   0  0.53 0.01
    3 L     2  -2  -2  -3  -1  -1  -2  -2  -2   1   2  -2   3   0  -2  -1  -1  -2  -1   0   35   0   0   0   0   0   0   0   0   0  35   0  30   0   0   0   0   0   0   0  0.20 0.01
    4 S     0  -1   0   0  -2   0   1  -1  -1  -2  -3   0  -2  -3   3   3   2  -3  -2  -2    0   0   0   0   0   0  14   0   0   0   0   0   0   0  20  54  12   0   0   0  0.47 0.02
    5 S     0  -1   2  -1  -2  -1  -1   4  -1  -3  -3  -1  -2  -3  -2   2   1  -3  -2  -2    0   0  14   0   0   0   0  45   0   0   0   0   0   0   0  30  11   0   0   0  0.46 0.01
    6 A     0  -1   5   1  -2   0   0   0   0  -3  -3   0  -2  -3  -2   1   0  -4  -2  -3   10   0  79   0   0   0   0   0   0   0   0   0   0   0   0  11   0   0   0   0  0.60 0.03
    7 G     0  -3  -1  -2  -2  -2  -2   4  -3   1  -1  -2  -1  -2  -2  -1  -1  -3  -2   2    0   0   0   0   0   0   0  55   0  12   0   0   0   0   0   0   0   0   0  32  0.40 0.02

## 4.1.4 Making a Hidden Markov Model (HMM) with aligned sequences
This task actually comprises four steps:
1.	Create a FASTA file with the complete protein sequences of the matches of your protein search with bit score > 200. You might find one-liners useful for this.
The hit IDs must be filtered using bit score > 200, which corresponds to column 12:
-	$12 > 200 → filters hits with bit score > 200
-	{print $2} → extracts the ID of the matched protein
-	sort | uniq → removes duplicate entries
-	Save the IDs into high_bitscore_ids.txt
```
(blast_env) mjlopez@DESKTOP-MBGBGC4:/mnt/c/Users/mjlop/coding/coding$ awk '$12 > 200 {print $2}' test.faa.blast | sort | uniq > high_bitscore_ids.txt
```
Y luego extraer las sencuencias completas desde UniProt FASTA:
-	Busca los IDs en high_bitscore_ids.txt
-	-A1 → copia también la línea de secuencia
-	Crea el FASTA hits_bs200.fasta con las proteínas completas
```
(blast_env) mjlopez@DESKTOP-MBGBGC4:/mnt/c/Users/mjlop/coding/coding$ grep -A1 -F -f high_bitscore_ids.txt uniprot_Atha.fasta > hits_bs200.fasta
```
2.	Compute a multiple alignment of these sequences with Clustal Omega. Check the available output formats.
```
clustalo -i hits_bs200.fasta -o hits_bs200.aln --force --outfmt=clu
```
Where:
-	-i → FASTA de entrada
-	-o → archivo de alineamiento
-	--outfmt=clu → formato Clustal
-	--force → sobrescribe si el archivo existe
Otros formatos posibles:
-	--outfmt=fasta
-	--outfmt=msf
-	--outfmt=phylip

3.	Build a HMM out of these aligned sequences with hmmbuild
Tomamos el alineamiento multiple y construimos el HMM, lo guarda en arf6.hmm
```
(blast_env) mjlopez@DESKTOP-MBGBGC4:/mnt/c/Users/mjlop/coding/coding$ hmmbuild arf6.hmm hits_bs200.aln
```
4.	Scan the HMM against your sequence collection with hmmerscan and write a short report on the results. This should be deliverable **Exe6**.
PROBLEM. I started working on the exercises several days after attending the classes, and I could not remember where HMMER was installed. 
The availability of HMMER tools was checked using `which`. Although `hmmerscan` was not present, other HMMER binaries such as `hmmbuild` and `hmmsearch` were available. I leave here the steps I had to follow for future reference.
HMMER is be installed, at least partially, because the following command worked:

```
hmmbuild arf6.hmm hits_bs200.aln 
```
Initial attempt to use `hmmerscan`:
```
hmmerscan --tblout hmm_results.tbl arf6.hmm uniprot_Atha.fasta > hmm_results.txt
```
Result:
```hmmerscan: command not found
```
The `hmmerscan` command is not available on the system. Therefore, I systematically checked the HMMER installation.
I tried to check whether hmmerscan is in the PATH:
```which hmmerscan```
Resultado:
```(no output)```

`hmmerscan` is not installed or not accessible.
Next, I checked for other HMMER binaries:
```which hmmbuild
which hmmsearch
which hmmpress
```
Result:
```/usr/bin/hmmbuild
/usr/bin/hmmsearch
/usr/bin/hmmpress
```

HMMER is installed at the system level, but the installation is partial and does not include `hmmerscan`.
According to the HMMER documentation:
- `hmmsearch` compares an HMM profile against a sequence database
- `hmmerscan` compares a sequence against a database of HMM profiles
For this exercise, `hmmsearch` was used to scan the HMM profile against the protein sequence collection:
```
hmmsearch --tblout hmm_results.tbl arf6.hmm uniprot_Atha.fasta > hmm_results.txt
```
Finally, the results were verified:
```
ls -lh hmm_results.tbl hmm_results.txt
head hmm_results.tbl
```


## 4.1.5 Annotating function by predicted structural similarity
Search for homologous sequences and define protein domains of AT1G30330.2 using HHPred.  
Note that it is also possible to search for proteins with AlphaFold structural predictions using EBI BLAST against AlphaFoldDB.
**Exe 7)** Produce a table with: i) domains defined by the boundaries of matched entries from the Protein Data Bank and Pfam  ii) similar sequences in AlphaFoldDB
On the website https://toolkit.tuebingen.mpg.de/tools/hhpred, the protein sequence obtained from:
```grep -A 1 "AT1G30330.2" test.faa```
was submitted and compared against PDB and Pfam databases. The results were downloaded, and the summary table is included below.
```
No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
  1 4LDU_A Auxin response factor 5  99.9 1.2E-22 1.2E-27  146.2   5.1   65   16-80     45-109 (397)
  2 4LDV_A Auxin response factor 1  99.7 6.4E-18 6.6E-23  116.1   3.6   63   16-79     13-75  (362)
  3 8OJ2_A Auxin response factor;   99.7 4.5E-17 4.6E-22  111.7   4.0   64   15-79     17-80  (366)
  4 PF07202.18 ; Tcp10_C ; T-compl  54.4      24 0.00025   18.3   1.7   23   32-54     12-34  (35)
  5 2MX4_A Eukaryotic translation   36.6      54 0.00056   20.0   1.4   15   32-46     23-37  (45)
  6 PF11995.13 ; DUF3490 ; Domain   36.2      43 0.00044   23.9   1.1   11   25-35     11-21  (155)
  7 PF02408.26 ; CUB_2 ; CUB-like   29.9      67  0.0007   21.1   1.2   40   30-69      6-51  (121)
  8 PF08961.15 ; NRBF2 ; Nuclear r  27.5 1.4E+02  0.0015   23.4   2.7   48   17-73    132-179 (199)
  9 PF10105.14 ; DUF2344 ; Unchara  26.7      67 0.00069   22.0   0.8    8    5-12     38-45  (183)
 10 1SSZ_A Pulmonary surfactant-as  25.6      46 0.00048   20.3  -0.1   19   61-79      9-27  (34)
```
## 4.1.6 Annotating function on orthology grounds
Search for functional annotations for protein AT1G30330.2 with help from eggNOG-mapper. Make sure you set one-to-one orthologues only.
**Exe 8)** What are the GO terms and eggNOG orthology groups of this protein?
No puedo usar eggNOG mapper vía web en este momento porque el servicio público está fuera de servicio o devolviendo errores.

## 4.1.7 Annotating function with Gene Ontology (GO) terms
In this exercise I will use the GO-web browser QuickGo as it is an intuitive and weekly updated resource. Although I won’t be using them in this session, there are many other related tools out there, such as UniProt or Ensembl Plants.
The goal of this task is for you to learn how to: Search for GO annotations using different inputs (protein IDs, gene IDs and GO terms), search for all products annotated to specific GO term or GO ID and customize the research to better fit the user preferences.
Let’s practice then:
1.	Search for the GO terms and the functional categories of the following GO IDs GO:0009414, GO:0035618, GO:0016491. Tip : For multiple search GO IDs needs to be separated by a space.
2.	What are the GO ID and the functional category corresponding to photosynthesis?
3.	What are the immediate parent(s) and children of the photosynthesis GO term?
4.	Search for the GO annotation terms of the following protein A0A068LKP4,A0A097PR28, A0A059Q6N8? What do you observe?
5.	How many gene products are involved in leaf development? Give the GO ID corresponding to this term.
6.	How many proteins of Arabidopsis thaliana, Prunus perisca and Zea mays are assigned to the leaf development GO term. Tip : Zea mays. Taxonomy ID=4577
7.	Check the total number of BP annotations and proteins supported by the experimental evidence codes in both Arabidopsis thaliana and Prunus persica. (see the evidence codes) Tip : check the ‘Statistics’ box.
For arabidopsis: Your current result set contains 28,553 annotations to 8,929 distinct gene products.
In prunus pérsica: Your current result set contains 6 annotations to 5 distinct gene products.

**Exe 9)** Summarize your results in a table.
|Busqueda|	GO ID |	Term |Name	| Ontology |
|--------|------------|------|----------|----------|
|GO IDs	|GO:0009414	|response to water deprivation	|BP |
|	|GO:0035618	|regulation of auxin transport	|BP|
|	|GO:0016491	|oxidoreductase activity	|MF|
|Photosyntesis term	|GO:0015979|	photosynthesis	|BP|
|Photosyntesis parents	|GO:0008150 ||		BP|
|	|GO:0009765|	photosynthesis, light reaction||
|	
|	|GO:0008152	|Metabolic process	||
|Photosyntesis childrens|	GO:0009767	|photosynthetic electron transport in photosystem II ||	
|	|GO:0009768|
|	photosynthetic electron transport in photosystem I ||
|	|GO:0009769|	photosynthetic electron transport in photosystem II, oxygen evolving complex||
|A0A068LKP4|	-|||		
|A0A097PR28|	GO:0008270|||
||GO:0046872|||
||GO:0006281|||
||GO:0006289|||
||GO:0006351|||
||GO:0006974|||
||GO:0005634|||
||GO:0000439|||
|A 0A059Q6N8|	GO:0015979|||
||GO:0019684|||
||GO:0005737|||
||GO:0009535|||
||GO:0016020|||
||GO:0009523|||
|Leaf product|	GO:0048366	1018 gen products |||		
|Leaf development 	|GO:0048366 entradas en uniprot-->	|Arabidopsis thaliana:637;Prunus perisca:59;		Zea mays:178 ||
	

##4.1.8 Predicting 3D structure
Choose from the options here and model the structure of one of the sequences obtained in **Exe5**.
**Exe 10)** Save a screen capture of your model and a table with associated quality scores.
Done with ARF6 protein.
![Diagrama](coding/Ejercicio_10.png)

 

