---
title: De novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 22nd November 2018
---

# Data Quality Assessment

1. Run FastQC on your data. Then open the results using `firefox` or `fastqc`.
	```bash
	# Run FastQC with 4 threads on all the files in all the folders in this folder
	fastqc -t 4 */*.fastq.gz
	# Open all the FastQC reports in a browser
	firefox */*_fastqc.html
	```

	How many sequences are in each file?

	<details>
	<summary> Solution - click to expand </summary>

	* Bacteria/bacteria_Str1_R{1,2}.fastq.gz: 766616 sequences
	* Bacteria/bacteria_Str2_R{1,2}.fastq.gz: 957939 sequences
	* SRR492065/SRR492065_{1,2}.fastq.gz: 5354356 sequences
	* E01/E01_1_135x.fastq.gz: 87217 sequences

	</details>

2. What is the average GC% in each data set?

	<details>
	<summary> Solution - click to expand </summary>

	* Bacteria/bacteria_Str1_R{1,2}.fastq.gz: 40%
	* Bacteria/bacteria_Str2_R{1,2}.fastq.gz: 40%
	* SRR492065/SRR492065_{1,2}.fastq.gz: 40%
	* E01/E01_1_135x.fastq.gz: 49%

	</details>

3. Which quality score encoding is used?

	<details>
	<summary> Solution - click to expand </summary>

	Sanger / Illumina 1.9

	</details>

4. What does a quality score of 20 (Q20) mean?

	<details>
	<summary> Solution - click to expand </summary>

	An expectation of 1 error in 100bp.

	</details>

5. What does a quality score of 40 (Q40) mean?

	<details>
	<summary> Solution - click to expand </summary>

	An expectation of 1 error in 10000bp.

	</details>

6. For **Bacteria/bacteria_Str1_R{1,2}.fastq.gz**, what percentage should the G and C lines be at in the per base sequence plot, and why?

	<details>
	<summary> Solution - click to expand </summary>

	20, because the mean GC is 40% and G and C should be in equal proportions and therefore half of the mean GC%.

	</details>

7. For **Bacteria/bacteria_Str1_R{1,2}.fastq.gz**, what percentage should the A and T lines be at in the per base sequence plot, and why?

	<details>
	<summary> Solution - click to expand </summary>

	30, because the mean AT is 60% and A and T should be in equal proportions and therefore half of the mean AT%.

	</details>

8. What distribution should the per base sequence plot follow?

	<details>
	<summary> Solution - click to expand </summary>

	A Uniform distribution.

	</details>

9. What distribution should the per base GC plot follow?

	<details>
	<summary> Solution - click to expand </summary>

	A Gaussian/Normal distribution.

	</details>

10. What value should the per base GC distribution be centered on?

	<details>
	<summary> Solution - click to expand </summary>

	Average GC content.

	</details>

11. How much duplication is present in **Bacteria/bacteria_Str1_R{1,2}.fastq.gz**?

	<details>
	<summary> Solution - click to expand </summary>

	24% (R1) and 15% (R2).

	</details>

12. What is adapter read through?

	<details>
	<summary> Solution - click to expand </summary>

	When the sequence reads past the insert into the adapter sequence on the other end.

	</details>

13. Let's look at the adapter sequence in the **SRR492065/SRR492065_{1,2}.fastq.gz** fastq files. Illumina uses different adapters
for different libraries. It is important to know which adapter sequence it is. Since this is public data, it is sometimes difficult to
find out what the adapters were. Use `bbmerge` to discover the adapter sequence.
	```bash
	bbmerge.sh in=<read1.fastq.gz> in2=<read2.fastq.gz> outa=adapters.fa
	more adapters.fa
	```

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	bbmerge.sh in=SRR492065/SRR492065_1.fastq.gz in2=SRR492065/SRR492065_2.fastq.gz outa=adapters.fa
	more adapters.fa
	```

	```
	>Read1_adapter
	AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATNTCGTATGCCGTCTTNTGNTT
	>Read2_adapter
	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
	```
	</details>

14. Use the command below to view the reads that have matching adapter sequence in your files.
	```bash
	paste <( zcat SRR492065/SRR492065_1.fastq.gz ) <( zcat SRR492065/SRR492065_1.fastq.gz ) \
	| grep -A2 -B1 --colour=always "AGATCGGAAGAGC" | less -SR
	```

	In the next step we will use Trimmomatic to trim adapters. It needs the correct adapter file. Use
	`grep` to identify the necessary adapter file to use. Trimmomatic's adapter files can be found in
	`/opt/byod/byod/Trimmomatic-0.36/adapters/`.
	```bash
	grep -B1 --colour=always "AGATCGGAAGAGCACACGTCTGAACTCC" /opt/byod/byod/Trimmomatic-0.36/adapters/*PE*.fa
	grep -B1 --colour=always "AGATCGGAAGAGCGTCGTGTAGGGAAAG" /opt/byod/byod/Trimmomatic-0.36/adapters/*PE*.fa
	```

	Which adapter file should be used?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	$ grep -B1 --colour=always "AGATCGGAAGAGCACACGTCTGAACTCC" /opt/byod/byod/Trimmomatic-0.36/adapters/*PE*.fa
	/opt/byod/byod/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa->PE2_rc
	/opt/byod/byod/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
	$ grep -B1 --colour=always "AGATCGGAAGAGCGTCGTGTAGGGAAAG" /opt/byod/byod/Trimmomatic-0.36/adapters/*PE*.fa
	/opt/byod/byod/Trimmomatic-0.36/adapters/TruSeq2-PE.fa->PCR_Primer1_rc
	/opt/byod/byod/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
	--
	/opt/byod/byod/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa->PE1_rc
	/opt/byod/byod/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
	```

	Therefore the file to use is:

	```
	/opt/byod/byod/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa
	```

	Since we have paired end data, we only use the files with PE in their name. Also as we are looking to remove adapter
	read-through, we are searching for the reverse compliment of the adapter `*_rc`. These two sequences are only common
	to one file, so we will use that one.

	</details>

15. Run Trimmomatic on **SRR492065/SRR492065_{1,2}.fastq.gz** to only remove adapters. How many reads were trimmed for adapters?
	```bash
	TRIMMOMATIC=/opt/byod/byod/Trimmomatic-0.36/trimmomatic-0.36.jar
	ADAPTER_DB=/opt/byod/byod/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa
	java -jar "$TRIMMOMATIC" PE SRR492065/SRR492065_1.fastq.gz SRR492065/SRR492065_2.fastq.gz \
	SRR492065/SRR492065_1.clean_pair.fastq.gz SRR492065/SRR492065_1.clean_unpair.fastq.gz \
	SRR492065/SRR492065_2.clean_pair.fastq.gz SRR492065/SRR492065_2.clean_unpair.fastq.gz \
	ILLUMINACLIP:"$ADAPTER_DB":2:30:10
	```

	<details>
	<summary> Solution - click to expand </summary>

	```
	Input Read Pairs: 5354356 Both Surviving: 5339516 (99.72%) Forward Only Surviving: 11947 (0.22%) Reverse Only Surviving: 1888 (0.04%) Dropped: 1005 (0.02%)
	```

	</details>

16. Run Kraken on **SRR492065/SRR492065_{1,2}.fastq**. What do you see? Is this consistent with what you see from `kat gcp`?
	```bash
	KRAKEN_DB=/home/data/byod/minikraken_20141208
	kraken --threads 4 --db "$KRAKEN_DB" --fastq-input --gzip-compressed --paired <read_{1,2}.fastq.gz> > <kraken.tsv>
	kraken-report --db "$KRAKEN_DB" <kraken.tsv> > <kraken.rpt>
	ktImportTaxonomy <( cut -f2,3 <kraken.tsv> ) -o <krona.html>
	```

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	KRAKEN_DB=/home/data/byod/minikraken_20141208
	kraken --threads 4 --db "$KRAKEN_DB" --fastq-input --gzip-compressed --paired SRR492065/SRR492065_{1,2}.fastq.gz > SRR492065.kraken.tsv
	kraken-report --db "$KRAKEN_DB" SRR492065.kraken.tsv > SRR492065.kraken.rpt
	ktImportTaxonomy <( cut -f2,3 SRR492065.kraken.tsv ) -o SRR492065.krona.html
	firefox SRR492065.krona.html
	```

	To make an image open the html file, and click on the snapshot button. Then save the resulting image to `SRR492065.krona.svg`.

	![A Krona plot of the Kraken analysis of SRR492065.](Data_QC_Exercises_17-10-23/SRR492065.krona.svg)

	The Kraken analysis is consistent with the KAT gcp analysis since it shows at least three organisms in the sample; Enterococcus, Staphylococcus, and Cutibacterium. Enterococcus also shows a higher abundance than both Staphylococcus and Cutibacterium, which are both in similar proportions.

	</details>

17. Run Kraken on **Bacteria/bacteria_Str1_R{1,2}.fastq.gz**. This bacteria is a strain of Moritella viscosa. Why is the bacteria not identified here?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	kraken --threads 4 --db "$KRAKEN_DB" --fastq-input --gzip-compressed --paired Bacteria/bacteria_Str1_R{1,2}.fastq.gz > bacteria_Str1.kraken.tsv
	kraken-report --db "$KRAKEN_DB" bacteria_Str1.kraken.tsv > bacteria_Str1.kraken.rpt
	ktImportTaxonomy <( cut -f2,3 bacteria_Str1.kraken.tsv ) -o bacteria_Str1.krona.html
	firefox bacteria_Str1.krona.html
	```

	To make an image open the html file, and click on the snapshot button. Then save the resulting image to `bacteria_Str1.krona.svg`.

	![A Krona plot of the Kraken analysis of SRR492065.](Data_QC_Exercises_17-10-23/bacteria_Str1.krona.svg)

	No species is identified here. This is because the database is not comprehensive enough to classify the bacteria in the sample.
	Although a small portion is classified, this may be a close relative, or a false positive classification.

	</details>

