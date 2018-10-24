---
title: De novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 23rd November 2018
---
# Assembly Validation Exercises.

### K-mer Analysis Toolkit

KAT is useful tool for high accuracy sequence data. The spectra-cn (copy number spectra) graph shows a
decomposition of k-mers in the assembly vs k-mers in the reads.
The black portion are k-mers not present in the assembly, the red portion is found once in the assembly, and so on.
This shows the completeness of an assembly, i.e. are all the reads assembled into contigs representative of the sequence data.

4. Use KAT to compare the reads to the assemblies.
	```bash
	mkdir KAT
	cd KAT
	ln -s ../*.fasta .
	apply_KAT () {
		ASSEMBLY="$1" # The assembly is the first parameter to this function
		READ1="$2" # The first read pair is the second parameter to this function
		READ2="$3" # The second read pair is the third parameter to this function
		COMBINED_DATA=combined_data.fastq # The file combining the read data
		mkfifo "$COMBINED_DATA" && zcat "$READ1" "$READ2" > "$COMBINED_DATA" & # Make a named pipe and combine reads
		kat comp -t "${CPUS:-4}" -o "${ASSEMBLY/.fasta/_vs_reads.cmp}" "$COMBINED_DATA" "$ASSEMBLY" # Compare Reads to Assembly
		rm -f "$COMBINED_DATA" # Clean up
	}
	```

	* How complete are the assemblies?
	* What would an incomplete assembly show?
	* How do the assemblies compare?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	for FASTA in *.fasta; do
		apply_KAT "$FASTA" ../bacteria_R1.fastq.gz ../bacteria_R2.fastq.gz
	done
	```

	![ CLC contigs K-mer spectra-cn plot ](images/kat/CLC_contigs_vs_reads.cmp-main.mx.spectra-cn.png)

	![ Mira large contigs K-mer spectra-cn plot ](images/kat/Mira_large_contigs_vs_reads.cmp-main.mx.spectra-cn.png)

	![ Mira trimmed contigs K-mer spectra-cn plot ](images/kat/Mira_trimmed_data_contigs_vs_reads.cmp-main.mx.spectra-cn.png)

	![ Spades contigs K-mer spectra-cn plot ](images/kat/SPAdes_contigs_vs_reads.cmp-main.mx.spectra-cn.png)

	![ Spades trimmed contigs K-mer spectra-cn plot ](images/kat/SPAdes_trimmed_data_contigs_vs_reads.cmp-main.mx.spectra-cn.png)

	![ Velvet K-mer spectra-cn plot ](images/kat/Velvet_contigs_vs_reads.cmp-main.mx.spectra-cn.png)

	All the assemblies are nearly complete, although it can be seen that the CLC and Mira assemblies are more complete
	than the Velvet assembly. The velvet assembler has failed to assemble some low coverage portions of the genome. This is seen
	by the black portion in the left hand side of the distribution with high frequency.

	</details>

### Blobtools

During the sequence quality assessment stage we tried to discern whether contamination was present. Sometimes this is
not feasible at the read level. By plotting Contig GC content vs Contig Read Coverage we can look for clusters of contigs that
share similar coverage. The appearance of multiple clusters can indicate multiple organisms. Occasionally, contigs can also be
taxonomically classified, providing further evidence for contaminants.

5. Run Blobtools on each assembly. Blobtools requires both a BAM file as input and blast output for the classification step.
	```bash
	mkdir Blast
	cd Blast
	ln -s ../*.fasta .
	apply_Blast () {
		ASSEMBLY="$1" # The assembly is the first parameter to this function. The file must end in .fasta
		BLAST_DB=/home/data/opt-byod/nt-database/nt # The location of the blast database nt
		echo "Blast: $ASSEMBLY"
		blastn -task megablast -query "$ASSEMBLY" -db "$BLAST_DB" -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
                -culling_limit 5 -num_threads "${CPUS:-4}" -evalue 1e-25 -out "${ASSEMBLY/.fasta/_blast_alignment.tsv}"
	}
	# Run blast for each assembly. Each blast takes approximately 5 minutes (Use the time to solve another exercise).
	cd ..
	mkdir Blobtools
	cd Blobtools
	ln -s ../*.fasta .
	ln -s ../BWA/*.bam .
	ln -s ../Blast/*.tsv .
	apply_Blobtools () {
		ASSEMBLY="$1" # The assembly is the first parameter to this function. The file must end in .fasta
		BAM="$2" # The BAM file is the second parameter to this function
		BLAST="$3" # The BLAST file is the third parameter to this function
		BLOB_DB=/opt/blobtools/data/nodesDB.txt
		blobtools create -i "$ASSEMBLY" -b "$BAM" -t "$BLAST" -o "${ASSEMBLY/.fasta/_blobtools}" --db "$BLOB_DB"
		blobtools blobplot -i "${ASSEMBLY/.fasta/_blobtools}.blobDB.json" -o "${ASSEMBLY/.fasta/_blobtools}"
	}
	```

	* Is there contamination in the assembly?
	* Do any assemblies show strange clustering?
	* Why might coverage vary across contigs within an assembly?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	for FASTA in *.fasta; do
		apply_Blobtools "$FASTA" "${FASTA/.fasta/_bwa_alignment.bam}" "${FASTA/.fasta/_blast_alignment.tsv}"
	done
	```

	![ CLC contigs Blobplot ](images/blobtools/CLC_contigs_blobtools.CLC_contigs_blobtools.blobDB.json.bestsum.phylum.p7.span.100.blobplot.bam0.png)

	![ Mira large contigs Blobplot ](images/blobtools/Mira_large_contigs_blobtools.Mira_large_contigs_blobtools.blobDB.json.bestsum.phylum.p7.span.100.blobplot.bam0.png)

	![ Mira trimmed contigs Blobplot ](images/blobtools/Mira_trimmed_data_contigs_blobtools.Mira_trimmed_data_contigs_blobtools.blobDB.json.bestsum.phylum.p7.span.100.blobplot.bam0.png)

	![ Spades contigs Blobplot ](images/blobtools/SPAdes_contigs_blobtools.SPAdes_contigs_blobtools.blobDB.json.bestsum.phylum.p7.span.100.blobplot.bam0.png)

	![ Spades trimmed contigs Blobplot ](images/blobtools/SPAdes_trimmed_data_contigs_blobtools.SPAdes_trimmed_data_contigs_blobtools.blobDB.json.bestsum.phylum.p7.span.100.blobplot.bam0.png)

	![ Velvet contigs Blobplot ](images/blobtools/Velvet_contigs_blobtools.Velvet_contigs_blobtools.blobDB.json.bestsum.phylum.p7.span.100.blobplot.bam0.png)

	The Blobplots all indicate a single cluster. Some contigs show fairly high coverage in comparison to the rest of the genome which could be repetitive elements in the genome.

	</details>

### Kraken Taxonomic Classification

As we pointed out before, occasionally classification might not be informative at the read level. By applying Kraken to the longer contigs,
we can get a better idea of what is in the assembly as long as the classification database contains that information.

6. Run Kraken on each assembly.
 	```bash
	mkdir Kraken
	cd Kraken
	ln -s ../*.fasta .
	apply_Kraken () {
		ASSEMBLY="$1" # The assembly is the first parameter to this function. This file must end in .fasta
		KRAKEN_DB=/home/data/byod/minikraken_20141208 # The location of the kraken database
		echo "Running Kraken: $ASSEMBLY"
		kraken --threads "${CPUS:-4}" --db "$KRAKEN_DB" --fasta-input "$ASSEMBLY" > "${ASSEMBLY/.fasta/.kraken.tsv}"
		kraken-report --db "$KRAKEN_DB" "${ASSEMBLY/.fasta/.kraken.tsv}" > "${ASSEMBLY/.fasta/.kraken.rpt}"
		ktImportTaxonomy <( cut -f2,3 "${ASSEMBLY/.fasta/.kraken.tsv}" ) -o "${ASSEMBLY/.fasta/.krona.html}"
	}
	```

	* What is identifiable here?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	for FASTA in *.fasta; do
		apply_Kraken "$FASTA"
	done
	```

	![CLC contigs Kraken output](Assembly_Validation_Exercises_17-10-23/Kraken/CLC_contigs.krona.svg)

	![Mira large contigs Kraken output](Assembly_Validation_Exercises_17-10-23/Kraken/Mira_large_contigs.krona.svg)

	![Mira trimmed contigs Kraken output](Assembly_Validation_Exercises_17-10-23/Kraken/Mira_trimmed_data_contigs.krona.svg)

	![Velvet contigs Kraken output](Assembly_Validation_Exercises_17-10-23/Kraken/Velvet_contigs.krona.svg)

	Here we see Kraken classifying around half the contigs, in contrast to the reads where it struggled to find a match. There are a high
	number of identified species, which is in contrast to all other analysis that indicate a single organism. This indicates that the
	sample is more likley closely related to the identified organisms, also supported by the point that many contigs are unclassified as well.
	The real organism is not present in the database used with Kraken.

	</details>

