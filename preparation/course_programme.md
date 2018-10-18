# Course Program

 Time | Activity
:-------------|:----------------
Day 1                 |
09.00 - 09.15         | Introduction
09.15 - 10.00         | [**Session 1**](lectures/session_1.md): Getting started in Genome Assembly
10.00 - 10.15         | Coffee break / fika
10.15 - 10.45         | Session 1 Exercises/Discussion
10.45 - 11.15         | [**Session 2**](lectures/session_2.md): Sequence Data Quality Assessment - Corruption, Coverage, and Subsampling
11.15 - 12.00         | Session 2 Exercises/Discussion
12.00 - 13.00         | Lunch
13.00 - 13.45         | [**Session 3**](lectures/session_3.md): Sequence Data Quality Assessment - FastQC, Contamination, Trimming, Duplication, and filtering
13.45 - 15.00         | Session 3 Exercises/Discussion
15.00 - 15.15         | Coffee break / fika
15.15 - 16.00         | [**Session 4**](lectures/session_4.md): Sequence Data Quality Assessment - Kmer Analysis using KAT
16.00 - 17.00         | Session 4 Exercises/Discussion
Day 2                 |
09.00 - 09.30         | Discussion of previous days exercises
09.30 - 10.00         | [**Session 5**](lectures/session_5.md): Overview of de novo assembler principles
10.00 - 10.15         | Coffee break / fika
10.15 - 10.45         | [**Session 6**](lectures/session_6.md): Assembly Validation - Assembly Statistics, and misassembly detection
10.45 - 12.00         | Session 6 Exercises/Discussion
12.00 - 13.00         | Lunch
13.00 - 13.30         | [**Session 7**](lectures/session_7.md): Assembly Validation - Assembly Graphs and contamination detection
13.30 - 15.00         | Session 7 Exercises/Discussion 
15.00 - 15.15         | Coffee break / fika
15.15 - 15.45         | [**Session 8**](lectures/session_8.md): Assembly Validation - Gene space statistics and comparitive alignment (and wrap)
15.45 - 17.00         | Session 8 Exercises/Discussion
Day 3                 |
09.00 - 9.30          | Discussion of previous days exercises
09.30 - 10.00         | [**Session 9**](lectures/session_9.md): Reference guided Assembly, Small genome Assembly with Illumina and PacBio
10.00 - 10.15         | Coffee break / fika
10.15 - 12.00         | Session 9 Exercises/Discussion
12.00 - 13.00         | Lunch
13.00 - 13.30         | [**Session 10**](lectures/session_10.md): Scaffolding, Gap Filling, Assembly Reconcilliation
13.30 - 15.00         | Session 10 Exercises/Discussion
15.00 - 15.15         | Coffee break / fika
15.15 - 15.45         | [**Session 11**](lectures/session_11.md): Challenges in Eukaryote assembly
15.45 - 17.00         | Session 11 Exercises/Discussion

## Session 1: Getting started in Genome Assembly

Time: 45 m + 30 m

Learning Objectives:

* Investigating properties of the genome of interest
* Extract high quality DNA
* Choice of sequencing technology - overview of properties like length, data quantity, starting material, price, output
* Estimation of necessary computing resources
* Scales of difficulty in genome assembly - Chart of what is easy to achieve -> what is hard to achieve (explain contig, scaffold, consensus, phased, etc here)
* When to make a pilot study - Sequence an Illumina library to estimate properties, etc.

Learning Outcomes:

* Students will know which properties of the genome should be investigated before starting a sequencing project.
* Students will know what might hinder DNA extraction.
* Students will know how to select appropriate sequencing technology.
* Students will know how to estimate computational resources.
* Students will know how to estimate project difficultly.
* Students can explain what contig, scaffold, consensus, phased means.

Exercises:

* **Test Connection to UPPMAX** and sort out problems during following lecture.

* Can they describe properties of their genome of interest?
* What might hinder DNA extraction for their organism?
* Which sequencing technologies are suitable for assembly?
* What kind of computational resources might they need?
* How difficult is what they're trying to achieve?

## Session 2: Sequence Data Quality Assessment - Corruption, Coverage, and Subsampling

Time: 30 m + 45 m

Learning objectives:

* Understand what data corruption is and how it occurs and how to remedy it. - How to use `md5sum` and `file`.
* Understand the concepts of Depth of Coverage and Breadth of Coverage.
* Understand what are suitable ranges of these values for each platform.
* Understand how to calculate these values from their data.
* Understand how to remedy sequence data outside of the suitable ranges. 

Learning outcomes:

* The student will be able to tell if a file is corrupt.
* The student will be able to tell if they have enough coverage in their sequence data.
* The student will be able to modify their data to have suitable coverage.

Exercises:

* Run `md5sum`
* Run `file`
* Run data quantity calculation
* Run `seqtk`
* Run `bbnorm`?

## Session 3: Sequence Data Quality Assessment - FastQC, Contamination, Trimming, Duplication, and filtering

Time: 45 m + 1 hr 15 m

Learning objectives:
* Understand how to use FastQC
* Understand how to interpret specific plots
* Understand the purpose of trimming
* Understand how to remove adapters
* Understand how to remove duplication
* Understand when this is applicable for PacBio and Nanopore data
* Understand how to detect contamination
* Understand how to screen for contamination - including spike-in
* Understand how to filter for contamination
* Understand the impact of contamination databases
* Understand how to use FastQC screen, Kraken, screen for spike-in 

Learning outcomes:

* The student will be able to use and interpret FastQC analyses
* The student will be able to use and interpret Kraken analyses
* The student will be able to use and interpret FastQ screen analyses
* The student will be able to use appropriate tools to filter data.

Exercises:

* Run `fastqc`
* Run `kraken`
* Run `fastqscreen`
* Run `bbmap`
* Run `trimmomatic`
* Run deduplication
* Run bwa, join, filtering pipeline 

## Session 4: Sequence Data Quality Assesment - Kmer Analysis using KAT

Time: 45 m + 1 hr

Learning objectives:
* Understand what distinct and unique k-mers are
* Understand how to use KAT
* Understand how to interpret a k-mer histogram
* Understand how to estimate genome properties from a k-mer histogram - genome size, ploidy, heterozygosity, repeat content.
* Understand how to detect GC bias or contamination
* Understand how to detect biases between datasets

Exercises:

* Run `kat hist`, `kat gcp`, `kat comp`

## Session 5: An overview of de novo assembler principles

Time: 30 mins

Learning objectives:

* Understand that all assemblers look for overlaps between reads
* Understand that a graph of some form is made - that the graph can not always be visualised (note common assemblers that make fastg,gfa ?)
* Understand that a path through the graph needs to be found
* Understand what heterozygosity does
* Understand what repeats do
* General stages of assembly - Read Error correction, assembly, consensus
* Polishing

Excercises:

* none

## Session 6: Assembly Validation - Assembly Statistics, and misassembly detection

Time 30 m + 1 hr 15 m

Learning objectives:
* Assembly statistics - Number of contigs, N50, Nx50, etc.
* Quast 
* FRC, and TigMint
* IGV

Exercises

* Run `quast`
* Run `bwa`
* Run `FRC`
* Run `TigMint`
* Run `IGV`

## Session 7: Assembly Validation - Assembly Graphs and contamination detection

Time 30 m + 1 hr 30

Learning objectives:
* Detecting contamination 
* Kraken, Blast
* Bandage
* Blobtools, KAT cold?
* Filtering contamination

Exercises:

* Run `bandage`
* Run `kraken`
* Run `blast`
* Run `blobtools`
* Run `KAT cold`
* Make a filter table
* Filter contigs with `samtools`

## Session 8: Gene space statistics and comparitive alignment

Time 30 m + 1hr 15 m

Learning objectives:
* Assessing Gene space statistics: Busco
* Comparative assessment: Make and interpret Dot plots
	- Mauve, Mummer, Gepard, D-Genies


## Session 9: Reference guided assembly and Small genome assembly with Illumina and PacBio

Time 30 m + 1 hr 45 m 

Learning objectives:

* Practical example of how to assemble a prokaryote
* Assemble with Spades, MaSuRCA?, Canu, Miniasm
* Effect of polishing
* Recircularizing the assembly (extra)
* Learn how reference guided assembly differs from de novo assembly
* Learn how to make an assembly guided by the reference

Exercises:

* Run `bwa` + `bcftools` etc.
* Run `spades`
* Run `masurca`
* Run `pilon` + `bwa`
* Run `canu`
* Run `miniasm`
* Run `racon`
* Run `quiver`
* Run `minimap2`
* Eval gene space before and after polishing

Extra:

* Break and recircularize the assembly

## Session 10: Scaffolding, Gap Filling, Assembly Reconcilliation, 

Time: 30 m + 1 hr 30 m

Learning objectives:

* How to scaffold - long reads, 10x?
* How to gap fill using other data
* How to do Assembly reconcilliation

Exercises:

* ?
* Run Arcs, Links?
* Run Hi-C validation?
* Run PbJelly?

## Session 12: Challenges in eukaryote assembly

Time 30 m + 1 hr 15 m

Learning objectives:

* Current challenges in assembly and how to solve them
	- Repeats in long reads -> Marvel
	- Phasing -> Falcon unzip
	- Too little DNA -> WGA methods
	- Parameter optimization, e.g. choosing optimal K, optimal overlap length, optimal error rate

Exercises:

* Choice of exercise? Repetitive genome, phasing a genome, WGA sample?

