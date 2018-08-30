---
layout: default
title:  'Exercise: Data Quality Assessment'
---

# Exercise: Data Quality Assessment

## Exercises:

1. Use `md5sum` to calculate the checksums of the data files in the folder `/sw/courses/assembly/QC_Data/`.
Redirect (`>` operator) the output into a file called `checksums.txt` in your workspace.
	<details>
	<summary> Solution: (click here) </summary>
	
	The correct solution requires you to change to the directory the data is in. The
	file paths must reflect the directory structure you want to check.

	```bash
	cd /sw/courses/assembly/QC_Data/
	md5sum */* > /proj/g2017025/nobackup/$USER/checksums.txt
	cd /proj/g2017025/nobackup/$USER
	cat checksums.txt
	  b3acf269e095e93120cf7ee96126e3e4  Bacteria/bacteria_R1.fastq.gz
	  93a06cff64ecf32c4f56ebe98cf24d82  Bacteria/bacteria_R2.fastq.gz
	  9a5a3305130fe3e53785ac4fbf0e0584  Ecoli/E01_1_135x.fastq.gz
	```

	Incorrect Solution:
	
	```bash
	md5sum /sw/courses/assembly/QC_Data/*/* > checksums.txt
	cat checksums.txt
	  b3acf269e095e93120cf7ee96126e3e4  /sw/courses/assembly/QC_Data/Bacteria/bacteria_R1.fastq.gz
	  93a06cff64ecf32c4f56ebe98cf24d82  /sw/courses/assembly/QC_Data/Bacteria/bacteria_R2.fastq.gz
	  9a5a3305130fe3e53785ac4fbf0e0584  /sw/courses/assembly/QC_Data/Ecoli/E01_1_135x.fastq.gz
	```
	When these checksums are checked, they check the files in `/sw/courses/assembly/QC_Data/` rather
	than the files in your directory.
	</details>

2. Make a copy of the data in your workspace (note the `.` at the end of the command):

	```bash
	cp -vr /sw/courses/assembly/QC_Data/* .
	```

	Use `md5sum -c` to check the checksums are complete.

	**solution:**
	
	When you check the checksums of the transferred files, make sure the files you check are correct.
	```bash
	cp -vr /sw/courses/assembly/QC_Data/* .
	chmod u+w -R Bacteria Ecoli
	md5sum -c checksums.txt
	  Bacteria/bacteria_R1.fastq.gz: OK
	  Bacteria/bacteria_R2.fastq.gz: OK
	  Ecoli/E01_1_135x.fastq.gz: OK
	```


3. Use `file` to get the properties of the data files. In which format are they compressed?

	**solution:**
	
	```bash
	file */*
	  Bacteria/bacteria_R1.fastq.gz: gzip compressed data, was "bacteria_R1.fastq", from Unix, last modified: Thu Nov 10 15:16:50 2016
	  Bacteria/bacteria_R2.fastq.gz: gzip compressed data, was "bacteria_R2.fastq", from Unix, last modified: Thu Nov 10 15:16:52 2016
	  Ecoli/E01_1_135x.fastq.gz:     gzip compressed data, was "E01_1_135x.fastq", from Unix, last modified: Wed Oct  5 14:15:05 2016
	```
	This tells you all the files are gzip compressed.


4. Use `zcat` and `less` to inspect the contents of the data files. From which sequencing technology are the following files and do you notice anything else?

	a. `Bacteria/bacteria_R{1,2}.fastq.gz`

	b. `Ecoli/E01_1_135x.fastq.gz`

	**solution:**
	
	```bash
	zcat Bacteria/bacteria_R1.fastq.gz | less
	  - result omitted -
	zcat Bacteria/bacteria_R2.fastq.gz | less
	  - result omitted -
	zcat Ecoli/E01_1_135x.fastq.gz | less -S
	  -- result omitted -
	```

	a. Inspecting the headers of `Bacteria/bacteria_R{1,2}.fastq.gz` indicates the sequences
	come from the Illumina platform. One can also see the sequences are of different lengths,
	indicating the files have been filtered in some way already.

	b. Inspecting the headers of `Ecoli/E01_1_135x.fastq.gz` indicates the sequences come from
	the PacBio platform. This is further supported by the long sequences in the file.


5. Identify the different parts of the Illumina header information:
	```
	@HWI-ST486:212:D0C8BACXX:6:1101:2365:1998 1:N:0:ATTCCT
	```

	**solution:**
	
	Each part is:
	* Machine number
	* Run number
	* Flowcell ID
	* Lane
	* Tile
	* X-coordinate on the Tile
	* Y-coordinate on the Tile

	* First/Second in pair
	* If it has failed Illumina QC (Chastity)
	* Control bit
	* Barcode index sequence/ID


6. Identify the different parts of the Pacific Biosciences header information:
	```
	@m151121_235646_42237_c100926872550000001823210705121647_s1_p0/81/22917_25263
	```

	**solution:**
	
	Each part is:
	* m = movie
	* time of run (yymmdd_hhmmss)
	* Instrument Serial Number
	* SMRT Cell Barcode
	* Set Number (deprecated)
	* Part Number
	* ZMW hole number
	* Subread region (start _ end)


7. What does each tool in this compound command do, and what is the purpose of this command?
	```bash
	zcat *.fastq.gz | seqtk seq -A - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m
	```
	
	**solution:**
	
	```bash
	zcat *.fastq.gz      # concatenates compressed files to one output stream
	seqtk seq -A -       # seqtk is a toolkit for manipulating sequence data. The -A converts input to fasta output.
	grep -v "^>"         # grep searches for lines beginning (^) with the string > and excludes them (-v).
	tr -dc "ACGTNacgtn"  # tr translates characters from one set to another. The -dc deletes characters not in the "ACGTNacgtn" set.
	wc -m                # wc is the word count tool. wc -m counts characters.
	```
	The purpose of the command is to count the total number of bases in your files.


8. How many bases are in:

	a. `Bacteria/bacteria_R{1,2}.fastq.gz`?

	b. `Ecoli/E01_1_135x.fastq.gz`?

	**solution:**
	
	a. The number of bases in `Bacteria/bacteria_R{1,2}.fastq.gz` is:
	
	```bash
	zcat Bacteria/bacteria_R{1,2}.fastq.gz | seqtk seq -A - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m
	  225890464
	```

	b. The number of bases in `Ecoli/E01_1_135x.fastq.gz` is:

	```bash
	zcat Ecoli/E01_1_135x.fastq.gz | seqtk seq -A - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m
	  748508257
	```


9. In the data set `Ecoli/E01_1_135x.fastq.gz`, how many bases are in reads of size 10kb or longer?
	
	**solution:**
	
	```bash
	zcat Ecoli/E01_1_135x.fastq.gz | seqtk seq -A -L 10000 - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m
	  510546313
	```


10. Run FastQC on the data sets. How many sequences are in each file?
	
	**solution:**
	
	```bash
	fastqc -t 6 */*.fastq.gz
	firefox */*.html
	```
	
	* There are 766616 sequences in `Bacteria/bacteria_R{1,2}.fastq.gz`.
	* There are 87217 sequences in `Ecoli/E01_1_135x.fastq.gz`.

11. What is the average GC% in each data set?
	
	**solution:**
	
	* Bacteria/bacteria_R{1,2}.fastq.gz: 40%
	* E01/E01_1_135x.fastq.gz: 49%


12. Which quality score encoding is used?
	
	**solution:**
	
	All files use the Sanger / Illumina 1.9 quality score encoding.


13. What does a quality score of 20 mean?
	
	**solution:**
	
	An expectation of 1 error in 100bp.


14. What does a quality score of 40 mean?
	
	**solution:**
	
	An expectation of 1 error in 10000bp.


15. Which distribution should the **per base sequence** plot be similar to in the FastQC output for Illumina data?
	
	**solution:**
	
	A Uniform distribution.


16. Which distribution should the **per sequence GC** plot be similar to in the FastQC output for Illumina data?
	
	**solution:**
	
	A Normal/Gaussian distribution.


17. Which value should the **per sequence GC** distribution be centered on?
	
	**solution:**
	
	The average GC content.


18. How much duplication is present in `Bacteria/bacteria_R{1,2}.fastq.gz`?
	
	**solution:**
	
	24% in R1 and 15% in R2.


19. What is adapter read-through?
	
	**solution:**
	
	When the read sequence continues past the end of the DNA insert/fragment into the adapter sequence on the other end.


20. Use `trimmomatic` to trim adapters from the data set `Bacteria/bacteria_R{1,2}.fastq.gz`. The `trimmomatic` jar file
	can be found in `$TRIMMOMATIC_HOME`, and the adapter files can be found in `$TRIMMOMATIC_HOME/adapters/`.

	a. Trim only the adapters. How much is filtered out?

	b. Quality trim the reads as well. How much is filtered out?

	**solution:**
	
	a.
	```bash
	cd Bacteria
	java -jar $TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE bacteria_R1.fastq.gz bacteria_R2.fastq.gz \
	bacteria_R1.clean_pair.fastq.gz bacteria_R1.clean_unpair.fastq.gz \
	bacteria_R2.clean_pair.fastq.gz bacteria_R2.clean_unpair.fastq.gz \
	ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-PE.fa:2:30:10
	  TrimmomaticPE: Started with arguments:
 	  bacteria_R1.fastq.gz bacteria_R2.fastq.gz bacteria_R1.clean_pair.fastq.gz bacteria_R1.clean_unpair.fastq.gz bacteria_R2.clean_pair.fastq.gz bacteria_R2.clean_unpair.fastq.gz ILLUMINACLIP:/sw/apps/bioinfo/trimmomatic/0.36/milou/adapters/TruSeq3-PE.fa:2:30:10
	  Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
	  ILLUMINACLIP: Using 1 prefix pairs, 0 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
	  Quality encoding detected as phred33
	  Input Read Pairs: 766616 Both Surviving: 764633 (99.74%) Forward Only Surviving: 1983 (0.26%) Reverse Only Surviving: 0 (0.00%) Dropped: 0 (0.00%)
	  TrimmomaticPE: Completed successfully
	```

	b.
	```bash
	java -jar $TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE bacteria_R1.fastq.gz bacteria_R2.fastq.gz \
	bacteria_R1.clean_qc_pair.fastq.gz bacteria_R1.clean_qc_unpair.fastq.gz \
	bacteria_R2.clean_qc_pair.fastq.gz bacteria_R2.clean_qc_unpair.fastq.gz \
	ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 \
	TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	  TrimmomaticPE: Started with arguments:
 	  bacteria_R1.fastq.gz bacteria_R2.fastq.gz bacteria_R1.clean_qc_pair.fastq.gz bacteria_R1.clean_qc_unpair.fastq.gz bacteria_R2.clean_qc_pair.fastq.gz bacteria_R2.clean_qc_unpair.fastq.gz ILLUMINACLIP:/sw/apps/bioinfo/trimmomatic/0.36/milou/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	  Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
	  ILLUMINACLIP: Using 1 prefix pairs, 0 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
	  Quality encoding detected as phred33
	  Input Read Pairs: 766616 Both Surviving: 573596 (74.82%) Forward Only Surviving: 170485 (22.24%) Reverse Only Surviving: 5253 (0.69%) Dropped: 17282 (2.25%)
	  TrimmomaticPE: Completed successfully
	```
