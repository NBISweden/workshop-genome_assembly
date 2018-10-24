---
layout: default
title:  'Exercise: Illumina Assembly'
---

## Exercise: Illumina Assembly

In this exercise you will assemble genomes de novo using commonly used assembly software. You will work with Illumina data of _Rhodobacter sphaerioides_, data that was used in the GAGE-B comparison of assemblers. All settings used for the different programs are the ones used by the GAGE-B project. Due to time-restrictions we are forced to stick to assemblers that run quickly, but if you have an assembly project of your own you are encouraged to try other assemblers too. Remember that no assembler is best for all projects; if possible you need to try several.
One of the assemblers, MaSuRCA, will take a bit longer to run (about an hour) so be prepared to start this to run over lunch!

### Setup

Make a copy of the folder which includes the read data using `cp -vr /sw/courses/assembly/illumina_assembly/data .`

### Rhodobacter HiSeq data

Use the files `Rhodo_Hiseq_trimmed_read{1,2}.fastq` for the assemblies. 

#### Spades

In your working directory make a folder called spades, and load the spades module:

```bash
mkdir spades
module load bioinfo-tools
module load spades/3.10.1
```

Spades is very easy to run in the basic configuration. Usually you want to run with the "--careful" flag, but this will take too long for this exercise. It will take around 16 mins anyway, grab a coffee.

```bash
spades.py -t 8 --pe1-1 ../data/Rhodo_Hiseq_trimmed_read1.fastq --pe1-2 ../data/Rhodo_Hiseq_trimmed_read2.fastq -o Spades_Rhodobacter
```

In Spades_Rhodobacter you will now have a number of files, including contigs.fasta and scaffolds.fasta.

Take a look at the files using `less`. Can you see any regions where contigs have been scaffolded together?

#### Abyss

Make a new folder for the Abyss assembly and load the necessary modules:

```bash
mkdir abyss
module load abyss/2.0.2
module load bowtie
```

Start Abyss using the abyss-pe script:

```bash
abyss-pe k=31 l=1 n=5 s=100 np=8 name=asm lib='reads' reads='../data/Rhodo_Hiseq_trimmed_read1.fastq ../data/Rhodo_Hiseq_trimmed_read2.fastq' aligner=bowtie
```

Once done you will have two files called asm-contigs.fa and asm.scaffolds.fa.

#### SOAP denovo

Soap de novo is one of very few program that can assemble large genomes using high coverage Illumina data.

Make a directory for the Soap de novo assembly and load the module:

```bash
mkdir soapdenovo
module load soapdenovo/2.04-r240
```

SoapDeNovo requires a config file that you need to create yourself using a text editor like nano:

```bash
nano soap.config
```

Enter this information:

```
[LIB]

avg_ins=220

reverse_seq=0

asm_flags=3

rank=1

q1=../data/Rhodo_Hiseq_trimmed_read1.fastq

q2=../data/Rhodo_Hiseq_trimmed_read2.fastq
```

Exit and save the file by ctrl-x (if using nano) and answer yes when asked to save.

Start SoapDeNovo by:

```bash
SOAPdenovo-63mer all -s soap.config -o asm -F -R -E -w -u -K 55 -p 8 | tee SOAPdenovo.log
# tee saves the output to a file as well as printing to screen
```

Examine the contigs properties.

SoapDeNovo also comes with a GapCloser utlity that tries to improve the assemblies by closing gaps in the scaffolds. Try it out using:

```bash
GapCloser -b soap.config -a asm.scafSeq -o asm.new.scafSeq -t 8 | tee -a SOAPdenovo.log
```

Are there any improvements?

#### MaSuRCA

MaSuRCA takes quite some time to run, so start it over lunch! (...or if you're not into waiting,
there are finished results in `/sw/courses/assembly/illumina_assembly/assemblies`).
Make an directory called masurca and load the necessary modules:

```bash
mkdir masurca
module load MaSuRCA/3.2.1
```

Like SOAP, MaSuRCA also needs a configuration file. Luckily, it can make a template for you! Run

```bash
masurca -g masurca_config
```

to make a template, and open it in nano for editing.

Find the line called `PE` and add the FULL paths to your `Rhodo_Hiseq_trimmed_read1.fastq` and `Rhodo_Hiseq_trimmed_read2.fastq` files. Generally, you should **NOT** give pre-trimmed data to MaSuRCA, as this will deteriorate the assembly (it does it's own pre-processing), but in this case we do it to get a slightly shorter run time. Next change the prefix to `pe 220 20` to set the proper type, insert length, and standard deviation. Also find the JUMP, PACBIO, and OTHER lines, and add a `#` to remove them from the assembly.
Finally change the number of threads to 8, save the file, and exit nano.

Now run the command:

```bash
masurca masurca_config
```

To generate an assembly script, and FINALLY run:

```bash
./assemble.sh
```

to start the assembly!

Once done you will have two files called CA/10-gapclose/genome.ctg.fasta and CA/10-gapclose/genome.scf.fasta.

#### Quast

```bash
module load bioinfo-tools quast/4.5.4
```

Now load the masurca files into Quast together with the earlier spades, abyss, and soap-denovo contigs. Are there any major differences? How do the assemblies compare to each other?

Since the Rhodobacter genome is available, we will cheat and use it for comparison, by including it as a reference.

```bash
quast.py -R R_sphaeroides.fasta -o assembly_metrics -l label1[,label2,label3,...] -t 1 <draft_genome1.fasta> [<draft_genome2.fasta> <draft_genome3.fasta> ...]
```

#### Pilon

```bash
module load Pilon/1.22
module load bwa/0.7.17
module load samtools/1.5
```

Use Pilon to polish the assemblies. Do they improve? Does running Pilon again on the polished assemblies make any improvements?

```bash
bwa index <assembly>
bwa mem -t 4 <assembly> <read1> <read2> | samtools sort -@ 4 -T <temp_file_name> -O BAM -o <output_bwa_alignment.bam> -
samtools index <output_bwa_alignment.bam>
java -jar $PILON_HOME/pilon.jar --genome <assembly> --frags <output_bwa_alignment.bam> --threads 4 --outdir <polished_assembly_dir> --output <polished_assembly_prefix> --changes
```

#### Finding the best k

Use KmerGenie to find the optimal k.  

```bash
module load KmerGenie/1.7039
kmergenie list_of_readfiles.fofn
```

Also examine the assembly graphs from different values of k to see which gives the best layout. Do they agree?

The tool Bandage can be used to make an image of the graph.
```bash
module load Bandage/0.8.0
Bandage load assembly.gfa --draw
```
