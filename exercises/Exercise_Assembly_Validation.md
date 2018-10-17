---
layout: default
title:  'Exercise: Assembly Validation'
---
# Assembly Validation.

## Introduction

The following exercises are intended to introduce you to the tools involved in assembly validation.

## Unix Guidance

Many commands below are wrapped up in functions to make it easier for you to run.

Functions can be copy-and-pasted into your terminal window or you can write them in a file
and `source` the file (e.g. `source functions.sh`).

A function looks like the following:
```bash
function run_many_commands {
	ARG_1=$1  # This is the first parameter
	ARG_2=$2  # This is the second parameter
	ARG_3=$3  # This is the third parameter
	# Then follows a complex series of commands to get a result
	command_1 $ARG_1 $ARG_2 > new.file
	command_2 new.file $ARG_3 > result.file
}
```

Once this command is copied into your terminal, you can use the function to run the complex series of commands:
```bash
run_many_commands File1.ext File2.ext File3.ext
# Your results will now be in result.file as written in the function above
```

## Exercises

For these exercises we will use the assemblies generated during your Illumina exercise.

### QUAST

```bash
module load bioinfo-tools quast/4.5.4
```

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

### Read alignment statistics

```bash
module load bwa/0.7.17 samtools/1.5 gnuplot
export PATH="/proj/g2017025/tools/FRC_align/bin:$PATH"
```

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
	function align_reads {
		ASSEMBLY=$1 # The assembly is the first parameter to this function
		READ1=$2 # The first read pair is the second parameter to this function
		READ2=$3 # The second read pair is the third parameter to this function
		bwa index $ASSEMBLY # Index the assembly prior to alignment
		bwa mem -t 4 $ASSEMBLY $READ1 $READ2 | samtools sort -@ 4 -T ${ASSEMBLY/.fasta/} -O BAM -o ${ASSEMBLY/.fasta/_bwa_alignment.bam} -
		samtools index ${ASSEMBLY/.fasta/_bwa_alignment.bam}
		# bwa mem : Align reads to the assembly
		# samtools sort : Sort the output by coordinate
		#    -O BAM : save the output as a BAM file
		#    -@ <int> : use <int> cores
		#    -T <temp> : Write temporary files to <temp>.nnnn.bam
		# samtools index : index the BAM file
		samtools flagstat ${ASSEMBLY/.fasta/_bwa_alignment.bam} > ${ASSEMBLY/.fasta/_bwa_alignment.bam.stats}
		# samtools flagstat : basic alignment statistics
	}
	```

	To run the function above, copy the function into your terminal window and use in the following way:
	```bash
	align_reads SPAdes_assembly.fasta read1.fastq.gz read2.fastq.gz
	```
	This will then run `bwa index`, `bwa mem`, `samtools sort`, `samtools index`, and `samtools flagstat` in the correct
	order and with the correct parameters.

3. Next, run FRCBAM on the alignments.
	```bash
	mkdir FRC
	cd FRC
	ln -s ../BWA/*.bam .
	function apply_FRC {
		ALIGNMENT=$1 # The BAM alignment file is the first parameter to this function
		GENOME_SIZE=$2 # The estimated genome size is the second parameter to this function
		FRC --pe-sam $ALIGNMENT --genome-size $GENOME_SIZE --output ${ALIGNMENT/.bam/}
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

	* How do these results compare to the Quast results?

### K-mer Analysis Toolkit

```bash
module load KAT/2.3.4
```

KAT is useful tool for high accuracy sequence data. The spectra-cn (copy number spectra) graph shows a
decomposition of k-mers in the assembly vs k-mers in the reads.
The black portion are k-mers not present in the assembly, the red portion is found once in the assembly, and so on.
This shows the completeness of an assembly, i.e. are all the reads assembled into contigs representative of the sequence data.

4. Use KAT to compare the reads to the assemblies.
	```bash
	mkdir KAT
	cd KAT
	ln -s ../*.fasta .
	function apply_KAT {
		ASSEMBLY=$1 # The assembly is the first parameter to this function
		READ1=$2 # The first read pair is the second parameter to this function
		READ2=$3 # The second read pair is the third parameter to this function
		COMBINED_DATA=combined_data.fastq # The file combining the read data
		mkfifo $COMBINED_DATA && zcat $READ1 $READ2 > $COMBINED_DATA & # Make a named pipe and combine reads
		kat comp -t 4 -o ${ASSEMBLY/.fasta/_vs_reads.cmp} $COMBINED_DATA $ASSEMBLY # Compare Reads to Assembly
		rm -f $COMBINED_DATA # Clean up
	}
	```

	* How complete are the assemblies?
	* What would an incomplete assembly show?
	* How do the assemblies compare?

### Blobtools

```bash
module load blast/2.6.0+ blobtools/0.9.17 
```

During the sequence quality assessment stage we tried to discern whether contamination was present. Sometimes this is
not feasible at the read level. By plotting Contig GC content vs Contig Read Coverage we can look for clusters of contigs that
share similar coverage. The appearance of multiple clusters can indicate multiple organisms. Occasionally, contigs can also be
taxonomically classified, providing further evidence for contaminants.

5. Run Blast and Blobtools on each assembly. Blobtools requires both a BAM file as input and blast output for the classification step.
	```bash
	mkdir Blast
	cd Blast
	ln -s ../*.fasta
	function apply_Blast {
		ASSEMBLY=$1 # The assembly is the first parameter to this function
		echo "Blast: $ASSEMBLY"
		blastn -task megablast -query $ASSEMBLY -db $BLASTDB/nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
                -culling_limit 5 -num_threads 4 -evalue 1e-25 -out ${ASSEMBLY/.fasta/_blast_alignment.blast.out}
	}
	cd ..
	mkdir Blobtools
	cd Blobtools
	ln -s ../*.fasta .
	ln -s ../BWA/*.bam .
	function apply_Blobtools {
		ASSEMBLY=$1 # The assembly is the first parameter to this function
		BAM=$2 # The BAM file is the second parameter to this function
		BLAST=$3 # The BLAST file is the third parameter to this function
		NAMES_DB=/opt/byod/byod/taxdump/names.dmp # The location of the names database
		NODES_DB=/opt/byod/byod/taxdump/nodes.dmp # The location of the nodes database
		blobtools create -i $ASSEMBLY -b $BAM -t $BLAST -o ${ASSEMBLY/.fasta/_blobtools} --names $NAMES_DB --nodes $NODES_DB
		blobtools blobplot -i ${ASSEMBLY/.fasta/_blobtools}.blobDB.json -o ${ASSEMBLY/.fasta/_blobtools}
	}
	```

	* Is there contamination in the assembly?
	* Do any assemblies show strange clustering?
	* Why might coverage vary across contigs within an assembly?

### Kraken Taxonomic Classification

```bash
module load Kraken/1.0 Krona/2.7
```

As we pointed out before, occasionally classification might not be informative at the read level. By applying Kraken to the longer contigs,
we can get a better idea of what is in the assembly as long as the classification database contains that information.

6. Run Kraken on each assembly.

	```bash
	mkdir Kraken
	cd Kraken
	ln -s ../*.fasta .
	function apply_Kraken {
		ASSEMBLY=$1 # The assembly is the first parameter to this function
		KRAKEN_DB=/sw/courses/assembly/minikraken_20141208 # The location of the kraken database
		echo "Running Kraken: $ASSEMBLY"
		kraken --threads 4 --db $KRAKEN_DB --fasta-input $ASSEMBLY > ${ASSEMBLY/.fasta/.kraken.out}
		kraken-report --db $KRAKEN_DB ${ASSEMBLY/.fasta/.kraken.out} > ${ASSEMBLY/.fasta/.kraken.rpt}
		cut -f2,3 ${ASSEMBLY/.fasta/.kraken.out} > ${ASSEMBLY/.fasta/.krona.in}
		ktImportTaxonomy ${ASSEMBLY/.fasta/.krona.in} -o ${ASSEMBLY/.fasta/.krona.html}
	}
	```

	* What is identifiable here?

### BUSCO

```bash
module load BUSCO/2.0.1
```

Assessing gene space is a core aspect of knowing whether or not you have a good assembly.

Benchmarking Universal Single-Copy Orthologs
(BUSCO) sets are collections of orthologous groups with near-universally-distributed single-copy genes in each species, selected from
OrthoDB root-level orthology delineations across arthropods, vertebrates, metazoans, fungi, and eukaryotes. BUSCO groups were selected
from each major radiation of the species phylogeny requiring genes to be present as single-copy orthologs in at least 90% of the species;
in others they may be lost or duplicated, and to ensure broad phyletic distribution they cannot all be missing from one sub-clade.
The species that define each major radiation were selected to include the majority of OrthoDB species, excluding only those with unusually
high numbers of missing or duplicated orthologs, while retaining representation from all major sub-clades. Their widespread presence means
that any BUSCO can therefore be expected to be found as a single-copy ortholog in any newly-sequenced genome from the appropriate
phylogenetic clade. [Busco Supplementary Online Material.]

7. Run Busco on each assembly.
	```bash
	mkdir BUSCO
	cd BUSCO
	ln -s ../*.fasta .
	# run the following line to setup the busco environment on milou
	source $BUSCO_SETUP
	function apply_BUSCO {
		ASSEMBLY=$1 #The assembly is the first parameter to this function
		LINEAGE=$BUSCO_LINEAGE_SETS/bacteria_odb9
		BUSCO -i $ASSEMBLY -l $LINEAGE -c 4 -m genome -o ${ASSEMBLY/.fasta/_busco}
	}
	```

	* What does the gene space look like for each assembly?

### Comparative Alignment

```bash
module load mauve/2015-02-13
```

Comparative alignment is a useful tool to see how assemblies compare to each other. This can be useful to compare assemblies to a reference, or
to see if assemblies have large structural differences.

8. Use Mauve to compare all the assemblies. First use Mauve contig mover (`Tools > Move contigs`) to reorder contigs with respect to your chosen first assembly. Then use `align with ProgressiveMauve` to align the reordered contigs.
	```bash
	Mauve
	```

9. Use Gepard to see how assemblies compare to each other.
	```bash
	mkdir Gepard
	cd Gepard
	ln -s ../*.fasta .
	function apply_Gepard {
		ASSEMBLY1=$1 # An assembly is the first parameter to this function
		ASSEMBLY2=$2 # An assembly is the second parameter to this function
		GEPARD_HOME=/sw/courses/assembly/Tools
		SUBSTITUTION_MATRIX=$GEPARD_HOME/matrices/edna.mat
		java -Djava.awt.headless=true -Xmx1024m -cp $GEPARD_HOME/Gepard-1.40.jar org.gepard.client.cmdline.CommandLine \
			-seq1 $ASSEMBLY1 -seq2 $ASSEMBLY2 -matrix $SUBSTITUTION_MATRIX -outfile ${ASSEMBLY1/.fasta/}_vs_${ASSEMBLY2/.fasta/}_dotplot.png
	}
	```
	This bash snippet will generate every pairwise comparison of the set.
	```bash
	set -- *.fasta
	for a; do
		shift;
		for b; do
			<function> $a $b
		done
	done
	```

	* Do you see any large scale differences between the assemblies?
