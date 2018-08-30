---
layout: default
title:  'Exercise: Illumina Assembly'
---

## Exercise: Illumina Assembly

In this exercise you will assemble genomes de novo using commonly used assembly software. You will work with Illumina data of _Rhodobacter sphaerioides_, data that was used in the GAGE-B comparison of assemblers. All settings used for the different programs are the ones used by the GAGE-B project. Due to time-restrictions we are forced to stick to assemblers that run quickly, but if you have an assembly project of your own you are encouraged to try other assemblers too. Remember that no assembler is best for all projects; if possible you need to try several.
One of the assemblers, MaSuRCA, will take a bit longer to run (about an hour) so be prepared to start this to run over lunch! 

**Questions**

If you have questions about the lab after the course, you are welcome to contact me: martin.norling@nbis.se

### Setup

If you're not already logged in, start by logging in to UPPMAX using the earlier supplied instructions and be sure to log in to your reserved node.

Before you start working with the assemblers you need to setup a project structure to have easy access to the data, results, and run scripts.

To make certain you are in your home folder, type: `cd`.

Make a directory for all of today’s exercises using `mkdir illumina_assembly`, and enter it using `cd illumina_assembly`. Finally, make a directory called work, `mkdir work`. The directory "work" is where you can keep the files you work with and run things. Do your work in this directory, and then copy results back to the main project directory to keep them structured. (For real projects though - remember though that uppmax is **not** a backup! Data can, and has, been lost!)

Anyway - make a copy of the folder which includes the read data using `cp -r /proj/g2016024/nobackup/illumina_assembly/data .`
When you are done, type `ls -l` and you should see something like this:

```bash
student@milou1:~/illumina_assembly$ ls -l
total 64
drwxrwxr-x 2 student student 2048 Nov 15 11:14 data
lrwxrwxrwx 1 student student   36 Nov 15 11:08 work -> /home/student/glob/illumina_assembly
```

You are now ready to proceed to working with the first assembly program, Spades. You will start by using the HiSeq data for all exercises, and later if time allows, redo with the MiSeq data and compare results. 

### Part1, HiSeq data

#### Spades

Start by making a folder called spades in the illumina_assembly/work directory, and enter it with `cd spades`.

Now load the spades module:

```
module load bioinfo-tools
module load spades
```

Spades is very easy to run in the basic configuration. Usually you want to run with the "--careful" flag, but this will take too long for this exercise. It will take around 16 mins anyway, grab a coffee.

```
spades.py -t 8 --pe1-1 ../data/Rhodo_Hiseq_trimmed_read1.fastq --pe1-2 ../data/Rhodo_Hiseq_trimmed_read2.fastq -o SpadesOut
```

In SpadesOut you will now have a number of files, including contigs.fasta and scaffolds.fasta.

Take a look at the files using `less`. Can you see any regions where contigs have been scaffolded together?

We then calculate some statistics and generate plots using a program called "Quast". Go back to your project directory (it is useful to keep all reports in the same place), and use this command, but fill int the complete paths to your spades output files:

```
module load quast/3.2
quast.py -o quast_spades -l Spades_scaffolds,Spades_contigs -t 1 scaffolds.fasta contigs.fasta
```

Download the whole quast result-folder (spades) to your own computer using `scp` and click on the reports.html file. Any big differences between the scaffolds and contigs files?

(OBS! You can also supply a reference genome to Quast that it will compare your assemblies with. You can find a reference genome at /proj/g2016024/nobackup/illumina_assembly/data/reference/R_sphaeroides.fasta. This was copied into your project along with the rest of your files, so it's already available to you! Again, make sure that all the paths are correct and run quast again with the reference:

```
quast.py -R R_sphaeroides.fasta -o spades -l Spades_scaffolds,Spades_contigs -t 1 scaffolds.fasta contigs.fasta
```
Does it tell you anything about misassemblies?

Now go on the next assembly program:

#### Abyss

First go to your work directory and make a new folder using `mkdir Abyss`, and enter it.

Now load the necessary modules:

```
module load abyss/1.3.7
module load bowtie
```

If you have paired end data you can start Abyss using the abyss-pe script:

```
abyss-pe k=31 l=1 n=5 s=100 np=8 name=asm lib='reads' reads=' ../../data/Rhodo_Hiseq_trimmed_read1.fastq ../../data/Rhodo_Hiseq_trimmed_read2.fastq' aligner=bowtie
```

Once done you will have two files called asm-contigs.fa and asm.scaffolds.fa. Now load these files into Quast together with the earlier Spades contigs. Can you based on these numbers say which assembler does the best job? (Note that this is a trick question!)

The next assembler we'll try is

#### SOAP denovo

SoapDeNovo is one of very few program that can assemble large genomes using high coverage Illumina data.

First, load the module:

```
module load soapdenovo/2.04-r240
```

Now make and enter a folder called 'soap'

SoapDeNovo requires a config file that you need to create yourself using a text editor like nano:

```
nano soap.config
```

Enter this information:

```
[LIB]

avg_ins=220

reverse_seq=0

asm_flags=3

rank=1

q1=../../data/Rhodo_Hiseq_trimmed_read1.fastq

q2=../../data/Rhodo_Hiseq_trimmed_read2.fastq
```

Exit and save the file by ctrl-x (if using nano) and answer yes when asked to save.

Start SoapDeNovo by:

```
SOAPdenovo-63mer all -s soap.config -o asm -F -R -E -w -u -K 55 -p 8 >>SOAPdenovo.log
```

Check the result-files asm.contig and asm.scafSeq for N50 size and number of contigs/scaffolds and compare with earlier results using Quast.

SoapDeNovo also comes with a GapCloser utlity that tries to improve the assemblies by closing gaps in the scaffolds. Try it out using:

```
GapCloser -b soap.config -a asm.scafSeq -o asm.new.scafSeq -t 8 >> SOAPdenovo.log
```

Any improvements? 

#### MaSuRCA

MaSuRCA takes quite some time to run, so start it over lunch! (...or if you're not into waiting, there are finished results in `/proj/g2016024/nobackup/illumina_assembly/assemblies/`)
Once again - start by making an output directory, this time called "masurca", and load the necessary modules:

```
module load MaSuRCA/3.2.1
```

Like SOAP, MaSuRCA also needs a configuration file. Luckily, it can make a template for you! Run

```
masurca -g masurca_config
```

to make a template, and open it in nano for editing.

Find the line called `PE` and add the FULL paths to your `Rhodo_Hiseq_trimmed_read1.fastq` and `Rhodo_Hiseq_trimmed_read2.fastq` files. Generally, you should **NOT** give pre-trimmed data to MaSuRCA, as this will deteriorate the assembly (it does it's own pre-processing), but in this case we do it to get a slightly shorter run time. Next change the prefix to `pe 220 20` to set the proper type, insert length, and standard deviation. Also find the JUMP, PACBIO, and OTHER lines, and add a `#` to remove them from the assembly.
Finally change the number of threads to 8, save the file, and exit nano.

Now run the command:

```
masurca masurca_config
```

To generate an assembly script, and FINALLY run:

```
./assemble.sh
```

to start the assembly!

Once done you will have two files called CA/10-gapclose/genome.ctg.fasta and CA/10-gapclose/genome.scf.fasta. Now load these files into Quast together with the earlier Spades, abyss, and soap-denovo contigs. Are there any major differences? Was running for an hour instead of a few minutes worth it in this particular case? 

### Part 2, MiSeq data

There is also MiSeq data for the same organism. Links are in the data folder. You should now try running the three programs you already tried using MiSeq data to see if you get any improvements. If you feel adventurous, you can also try to run the programs with both HiSeq and MiSeq at the same time.

*IMPORTANT! Remember to work in different directories and/or change the name of output directories so that you do not overwrite your old data!*

The MiSeq data have longer reads, and you therefore need to change the following parameters (and remember to use the MiSeq files this time):

- Abyss-pe - change k to 49
- SoapDeNovo - change in the config file avg_ins to 540, use on the command line a K value of 79. **Note**: when using a kmer size greater than 63, you'll need to use `SOAPdenovo-127mer` instead of `SOAPdenovo-63mer`.
- MaSuRCA - change insert size from 220 to 540, and standard deviation to 50.

Compare with your HiSeq results. Differences? 

### Part 3, Not had enough?

This is an optional exercise. Here you should try to use Mira on the HiSeq data. Load the module using: `module load mira/4.0rc4`

Here you will receive no help. Try to figure out how to use the program by googling, in particular try to find the manual. Once you get it to run, kill it. No, I am serious. Press ctrl-c to kill the process. It simply takes to long to run, but now you have probably learned a lot by trying to get it to work. smile Then check the result-files in `/proj/g2015027/private/assembly_workshop/mira/`

Was the longer running-time worth it?

By **Henrik Lantz** and **Martin Norling**, NBIS/SciLife/Uppsala University 
