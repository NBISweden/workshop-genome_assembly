---
layout: default
title:  'Exercise: Data Quality Assessment'
---

# Exercise: Data Quality Assessment

## Loading modules on Milou / Rackham:

To use bioinformatic tools on Milou / Rackham, first the library of tools must be made available using the command:

```bash
module load bioinfo-tools
```

Then specific tools can be loaded in a similar fashion. If a particular version is needed, it can be appended to the end.

```bash
module load FastQC/0.11.5
module load seqtk/1.2-r101
module load trimmomatic/0.36
```

If you have trouble finding a tool, use the `module spider` function to search.

```bash
module spider fastqc
```

## Exercises:

1. Use `md5sum` to calculate the checksums of the data files in the folder `/sw/courses/assembly/QC_Data/`.
Redirect (`>` operator) the output into a file called `checksums.txt` in your workspace.

2. Make a copy of the data in your workspace (note the `.` at the end of the command):

	```bash
	cp -vr /sw/courses/assembly/QC_Data/* .
	```

	Use `md5sum -c` to check the checksums are complete.

3. Use `file` to get the properties of the data files. In which format are they compressed?

4. Use `zcat` and `less` to inspect the contents of the data files. From which sequencing technology are the following files and do you notice anything else?

	a. `Bacteria/bacteria_R{1,2}.fastq.gz`

	b. `Ecoli/E01_1_135x.fastq.gz`

5. Identify the different parts of the Illumina header information:
```
@HWI-ST486:212:D0C8BACXX:6:1101:2365:1998 1:N:0:ATTCCT
```

6. Identify the different parts of the Pacific Biosciences header information:
```
@m151121_235646_42237_c100926872550000001823210705121647_s1_p0/81/22917_25263
```

7. What does each tool in this compound command do, and what is the purpose of this command?
```bash
zcat *.fastq.gz | seqtk seq -A - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m
```

8. How many bases are in:

	a. `Bacteria/bacteria_R{1,2}.fastq.gz`?

	b. `Ecoli/E01_1_135x.fastq.gz`?

9. In the data set `Ecoli/E01_1_135x.fastq.gz`, how many bases are in reads of size 10kb or longer?

10. Run FastQC on the data sets. How many sequences are in each file?

11. What is the average GC% in each data set?

12. Which quality score encoding is used?

13. What does a quality score of 20 mean?

14. What does a quality score of 40 mean?

15. Which distribution should the **per base sequence** plot be similar to in the FastQC output for Illumina data?

16. Which distribution should the **per sequence GC** plot be similar to in the FastQC output for Illumina data?

17. Which value should the **per sequence GC** distribution be centered on?

18. How much duplication is present in `Bacteria/bacteria_R{1,2}.fastq.gz`?

19. What is adapter read-through?

20. Use `trimmomatic` to trim adapters from the data set `Bacteria/bacteria_R{1,2}.fastq.gz`. The `trimmomatic` jar file
can be found in `$TRIMMOMATIC_HOME`, and the adapter files can be found in `$TRIMMOMATIC_HOME/adapters/`.

	a. Trim only the adapters. How much is filtered out?

	b. Quality trim the reads as well. How much is filtered out?
