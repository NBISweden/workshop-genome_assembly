---
layout: default
title:  'Exercise: Assembly Assessment'
---

# Exercise: Assembly Assessment

**Questions**

If you have questions about the lab after the course, you are welcome to contact me: martin.norling@nbis.se

## Introduction

In this exercise we will look further into assembly assessment. We will focus on illumina assemblies, as the tools for these assemblies are more mature, but many things can be used for pacbio as well.

Start with one assembly and go through all the steps, then continue with the other assemblies if you have the time.

## Read Mapping

As the  assemblers we used for illumina assembly didn't provide read mappings, we need to map the reads back on the assemblies. We will do this with BWA, The Burrows-Wheeler Aligner. We will need to do this for each assembly. These parts can be stored in a script, like in the illumina assembly exercise, or run directly on the command line - the choice is yours!

To run BWA, we need to load the `bioinfo-tools`, `bwa/0.7.13`, and `samtools/1.3` modules.

The first thing you need to do to get a mapping is to index the target sequence, in our case this is the assembly.

```
$ bwa index <assembly.fasta>
```

We can then map the raw reads to the assembly using the "mem" algorithm. Run `$ bwa mem` to get some help on running BWA, then map the the reads to the assembly like:

```
$ bwa mem -t 8 assembly.fasta reads.fastq >assembly_lib_type.sam
```

This will produce a SAM file. These are bulky plain-text files, so we will sort and convert it into a BAM file, which are more efficient. You can convert the mapping using samtools like:

```
$ samtools view -bS assembly_lib_type.sam -o assembly_lib_type.bam
```

We can now sort the BAM file to get a file that can be more efficiently parsed:

```
$ samtools sort assembly_lib_type.bam -o assembly_lib_type.sorted.bam
```

And finally we can index the sorted BAM file:

```
$ samtools index assembly_lib_type.sorted.bam
```

Now that we have a sorted, indexed BAM file, we can also remove the SAM file, and the unsorted BAM file to save some space and keep the project neat! Now - what do we do with a BAM file? We will use it a little bit more later - but for now, we can use picard tools to get some insert statistics to see whether our reads seem to be in the correct place. Picard tools is a great set of tools, but they suffer from being a tad user hostile. To spare you some time, just use the following commands: (just replace <assembly name> and <sorted-bam-filename>

```
module load bioinfo-tools picard/2.0.1
java -Xmx16g -XX:PermSize=8g -jar $PICARD_HOME/picard.jar CollectInsertSizeMetrics MINIMUM_PCT=0 HISTOGRAM_FILE=<assembly name>.pdf  INPUT=<sorted-bam-filename>  OUTPUT=<assembly name>.sorted.collectInseSize HISTOGRAM_WIDTH=500
```

This will produce a pdf report of insert sizes. Download this file using `scp`, and have a look.

## KMER analysis

A good continuation of the validation is looking at the kmer content of the assembly. Before the assembly we looked at the kmer content of the read set, and tried to determine if this information was suitable for assembly. Now we need to check if that information is actually present in the assembly as well as in the read set.

### KAT

Using KAT again (You will need the modules: `KAT/2.1.1` and `gnuplot/4.6.5`) – we can plot the kmer content of the assembly compared to the kmer content of the read set. The first thing we need to do is to combine the reads into a single file, for gzipped files, this can be done with `zcat`, or for unzipped files `cat`.
 
Ex.
```
$ cat reads_R1.fastq >> combined.fastq
$ zcat reads_R2.fastq.gz >> combined.fastq
```

We will now use `kat comp` to create a kmer content comparison.
Use `kat comp --help` to get help for the program, then create a comparison between the combined reads and the assembly. Make sure that you use the flags for **canonical hashes** for both sequence 1 and 2, as well as **8 threads**.
Finally, clean up you working directory by removing the combined fasta file, and re-zipping any unzipped files. Then download the output files to you computer using `scp` and look at the png file that was produced.

 - Does the kmer content look good to you?
 - How much of the kmer “noise” is part of the final assembly?
 - What do you think contamination would look like in the kmer plot?
 - What can the kmer graph tell you about the ploidy of the organism?
 
### REAPR

[Reapr](http://www.sanger.ac.uk/science/tools/reapr) is a tool trying to find explicit errors in the assembly based on incongruently mapped reads. It is heavily based on too low span coverage, or reads mapping too far or too close to each other. The program will also break up contigs/scaffolds at spurious sites to form smaller (but hopefully correct) contigs. Reapr runs pretty slowly, sadly, 

Reapr is a bit fuzzy with contig names, but luckily it's given us a tool to check if things are ok before we proceed!
The command `reapr facheck <assembly.fasta>` will tell you if everything's ok! in this case, no output is good output, since the only output from the command is the potential problems with the contig names.
If you run into any problems, run `reapr facheck <assembly.fasta> <renamed_assembly.fasta>`, and you will get an assembly file with renamed contigs.

Once the names are ok, we continue:

The first thing we reapr needs, is a list of all "perfect" reads. This is reads that have a perfect map to the reference. Reapr is finicky though, and can't use libraries with different read lengths, so you'll have to use assemblies based on the raw data for this. Run the command `reapr perfectmap` to get information on how to create a perfect mapping file, and create a perfect mapping called `<assembler>_perfect`. This should take about a minute.

The next tool we need is `reapr smaltmap` which creates a bam file of read-pair mappings. Do the same thing you did with `perfectmap` and create an output file called `<assembler>_smalt.bam`. This should take about twenty minutes.

Finally we can use the smalt mapping, and the perfect mapping to run the reapr pipeline. Run `reapr pipeline` to get help on how to run, and then run the pipeline. Store the results in `reapr_<assembler>`.  This should take about ten minutes.

There are several checks you can do after running Reapr (detailed [here](ftp://ftp.sanger.ac.uk/pub/resources/software/reapr/Reapr_1.0.18.manual.pdf)) but for now we'll stick to looking at the split output file, called `04.break.broken_assembly.fa`. Use this file together with the original assembly to generate a quast report. How does the results look after reapr?

