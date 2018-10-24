---
title: De novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 23rd November 2018
---
# Assembly Validation Exercises.

## Introduction

The following exercises are intended to introduce you to the tools involved in assembly validation.

Many commands below are wrapped up in functions to make it easier for you to run.

Functions can be copy-and-pasted into your terminal window or you can write them in a file
and `source` the file (e.g. `source functions.sh`).

A function looks like the following:
```bash
run_many_commands () {
	ARG_1="$1"  # This is the first parameter
	ARG_2="$2"  # This is the second parameter
	ARG_3="$3"  # This is the third parameter
	# Then follows a complex series of commands to get a result
	command_1 "$ARG_1" "$ARG_2" > new.file
	command_2 new.file "$ARG_3" > result.file
}
```

Once this command is copied into your terminal, you can use the function to run the complex series of commands:
```bash
run_many_commands File1.ext File2.ext File3.ext
# Your results will now be in result.file as written in the function above
```

## Exercises

### QUAST

QUAST is a good starting point to help evaluate the quality of assemblies. It provides many helpful contiguity statistics.

1. Run QUAST on all the assemblies at once and generate a report.
	```bash
	quast.py -t 4 --est-ref-size <int> <draft_assembly1.fasta> <draft_assembly2.fasta> [ <draft_assembly3.fasta ... ]
	# -t 4 : use 4 threads
	# --est-ref-size <int> : Estimated reference size (for computing NGx metrics without a reference)
	# <draft_assembly1.fasta> <draft_assembly2.fasta> [ <draft_assembly3.fasta ... ] : All the draft assemblies you want to compare.
	```

	* What is the ideal shape of the Nx graph for a bacteria?
	* What can you see in the Nx graph?
	* Does the GC content graph indicate contamination?
	* Which assembly has the best contiguity metrics?

	<details>
	<summary> Solution - click to expand </summary>

	First run Quast on all the assemblies.
	```bash
	quast.py -t 4 --est-ref-size 5000000 *.fasta
	```

	The ideal shape of the Nx graph would be vertical line if the bacteria is a single contig. As more contigs are present, the cumulative
	length graph shifts to the right slightly as contigs get shorter.

	The Nx graph shows that the CLC assembly is the making the longest and least contigs.

	The GC content graph does not indicate contamination.

	The CLC assembly has the best contiguity metrics.

	![Quast Cumulative Length Plot](images/quast/cumulative_plot.png)

	![Quast NGx Plot](images/quast/NGx_plot.png)

	![Quast GC Plot](images/quast/GC_content_plot.png)

	</details>

### Read alignment statistics

Read congruency is an important measure in determining assembly accuracy. Clusters of read pairs that align incorrectly are
strong indicators of mis-assembly.

2. How well do the reads align back to the draft assemblies? Use `bwa` and `samtools` to assess the basic alignment statistics.

	Make a folder for your results.
	```bash
	mkdir BWA
	cd BWA
	ln -s ../*.fasta . # link the fasta files in this directory
	```

	Then copy this function into your terminal.

	```bash
	align_reads () {
		ASSEMBLY="$1" # The assembly is the first parameter to this function. Must end in .fasta
		READ1="$2" # The first read pair is the second parameter to this function
		READ2="$3" # The second read pair is the third parameter to this function
		bwa index "$ASSEMBLY" # Index the assembly prior to alignment
		bwa mem -t "${CPUS:-4}" "$ASSEMBLY" "$READ1" "$READ2" | samtools sort -@ 4 -T "${ASSEMBLY/.fasta/}" -O BAM -o "${ASSEMBLY/.fasta/_bwa_alignment.bam}" -
		samtools index "${ASSEMBLY/.fasta/_bwa_alignment.bam}"
		# bwa mem : Align reads to the assembly
		# samtools sort : Sort the output by coordinate
		#    -O BAM : save the output as a BAM file
		#    -@ <int> : use <int> cores
		#    -T <temp> : Write temporary files to <temp>.nnnn.bam
		# samtools index : index the BAM file
		samtools flagstat "${ASSEMBLY/.fasta/_bwa_alignment.bam}" > "${ASSEMBLY/.fasta/_bwa_alignment.bam.stats}"
		# samtools flagstat : basic alignment statistics
	}
	```

	To run the function above, copy the function into your terminal window and use in the following way:
	```bash
	align_reads SPAdes_assembly.fasta read1.fastq.gz read2.fastq.gz
	```
	This will then run `bwa index`, `bwa mem`, `samtools sort`, `samtools index`, and `samtools flagstat` in the correct
	order and with the correct parameters.

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	for FASTA in *.fasta; do
		align_reads "$FASTA" ../bacteria_R1.fastq.gz ../bacteria_R2.fastq.gz
	done
	```

	```shell
	# CLC contigs alignment rate
	1510883 + 0 mapped (98.08% : N/A)
	1451556 + 0 properly paired (94.67% : N/A)

	# Mira large contigs alignment rate
	1508358 + 0 mapped (98.35% : N/A)
	1462312 + 0 properly paired (95.37% : N/A)

	# Mira trimmed contigs alignment rate
	1506432 + 0 mapped (98.19% : N/A)
	1438640 + 0 properly paired (93.83% : N/A)

	# Spades contigs alignment rate
	1508412 + 0 mapped (98.35% : N/A)
	1445892 + 0 properly paired (94.30% : N/A)

	# Spades trimmed contigs alignment rate
	1509790 + 0 mapped (98.36% : N/A)
	1451312 + 0 properly paired (94.66% : N/A)

	# Velvet contigs alignment rate
	1483277 + 0 mapped (96.29% : N/A)
	1404476 + 0 properly paired (91.60% : N/A)
	```

	In contrast to the Quast results, the alignment rate is showing greater alignment statistics for the Mira contigs.

	</details>

3. Next, run FRCBAM on the alignments.
	```bash
	mkdir FRC
	cd FRC
	ln -s ../BWA/*.bam .
	apply_FRC () {
		ALIGNMENT="$1" # The BAM alignment file is the first parameter to this function. The filename must end in .bam.
		GENOME_SIZE="$2" # The estimated genome size is the second parameter to this function
		FRC --pe-sam "$ALIGNMENT" --genome-size "$GENOME_SIZE" --output "${ALIGNMENT/.bam/}"
	}
	```

	Plot the FRC curves (`<output>_FRC.txt`) together in a plot using Gnuplot. Which assembly has the best feature curve?
	```bash
	gnuplot << EOF
	set terminal png size 800,600
	set output 'FRC_Curve_all_assemblies.png'
	set title "FRC Curve" font ",14"
	set key right bottom font ",8"
	set autoscale
	set ylabel "Approximate Coverage (%)"
	set xlabel "Feature Threshold"
	files = system('find -name "*alignment_FRC.txt"')
	plot for [data in files] data using 1:2 with lines title data
	EOF
	```

	How do these results compare to the Quast results?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	for BAM in *.bam; do
		apply_FRC "$BAM" 5000000
	done
	gnuplot << EOF
	set terminal png size 800,600
	set output 'FRC_Curve_all_assemblies.png'
	set title "FRC Curve" font ",14"
	set key right bottom font ",8"
	set autoscale
	set ylabel "Approximate Coverage (%)"
	set xlabel "Feature Threshold"
	files = system('find -name "*alignment_FRC.txt"')
	plot for [data in files] data using 1:2 with lines title data
	EOF
	```

	![FRC curve](images/frc/FRC_Curve_all_assemblies.png)

	The FRC curve indicates that the Velvet contigs assembly have the least features (misassembly signals), i.e. is the most correct.

	</details>

