---
layout: default
title:  'Exercise: K-mer and contamination analysis'
---

# Exercise: K-mer and contamination analysis

## Loading modules on Milou / Rackham:

To use bioinformatic tools on Milou / Rackham, first the library of tools must be made available using the command:

```bash
module load bioinfo-tools
```

Then specific tools can be loaded in a similar fashion. If a particular version is needed, it can be appended to the end.

```bash
module load KAT/2.3.4
module load Kraken/1.0
module load Krona/2.7
```

If you have trouble finding a tool, use the `module spider` function to search.

```bash
module spider fastqc
```

## Unix help:

Many tools can work with compressed files. In the rare case a tool cannot, it may be able to read data from a pipe `|` by
passing information through STDIN `-`. If a tool needs a filename and cannot read from STDIN, a "named pipe" can help you
process compressed data without decompressing it.
```bash
# Make a named pipe called sequence.fastq so it mimics a file name.
$ mkfifo sequence.fastq
# Read data into the named pipe and put the process in the background (&)
$ zcat read1.fastq.gz read2.fastq.gz > sequence.fastq &
# Now run the command with the named pipe
$ command sequence.fastq
# A named pipe can only be used once. Remove it afterwards.
$ rm sequence.fastq
```

## Exercises:

1. What is a k-mer?

2. Use the following set of commands to extract the list of k-mers in `Bacteria/bacteria_R{1,2}.fastq.gz`.
	```bash
	mkfifo <named_pipe.fastq> && zcat <reads.fastq.gz> > <named_pipe.fastq> & # Make a named pipe and run in the background
	kat hist -t 4 -d -o <output.hist> <named_pipe.fastq> # Run KAT reading from the named pipe
	kat_jellyfish dump <output.hist>-hash.jf27 > <kmer.lst> # Use Jellyfish to print out a human readable list
	rm <named_pipe.fastq> # named pipes can only be used once, and so are removed after use.
	```
	The **<kmer.lst>** file has the following format.
	```
	>frequency
	kmer_sequence
	```
	How many distinct k-mers were found? Use the line count command `wc -l` to find out.

3. How many k-mers have a frequency of 1?  Use the following command to find out.
	```bash
	paste - - < kmer.lst | cut -c2- | awk '$1 == 1 { sum++ } END { print sum+0 }'
	# paste - - : reads two consecutive lines onto the same line.
	# cut -c2- : prints from the second character up to the last character in a line.
	# awk '$1 == 1 { sum++ } END { print sum+0 }' : if column 1 has a frequency of 1, increase the variable "sum". Print the value of "sum" at the end.
	```

4. How many k-mers have a frequency greater than 5?

5. A k-mer histogram was plotted using `kat hist` in a file `*.hist.png`. Open the image using `display` and estimate the mean k-mer frequency.

6. The following command prints the frequency of each k-mer frequency between 5 and 45. What is the mean k-mer frequency?
	```bash
	paste - - < kmer.lst | cut -c2- | awk '$1 > 5 && $1 < 45 {sum[$1]++ } END { for (freq in sum) {print freq" "sum[freq]} }' | sort -k1,1n
	# paste - - : reads two consecutive lines onto the same line.
	# cut -c2- : prints from the second character up to the last character in a line.
	# awk '$1 > 5 && $1 < 45 {sum[$1]++ } END { for (freq in sum) {print freq" "sum[freq]} }' :
	# 	if column 1 has a frequency greater than 5 and less than 45, increase the value of the array "sum[frequency]" by 1.
	# 	Then for each frequency in sum print the value of sum[frequency] at the end.
	# sort -k1,1n : Perform a numerical sort on the data sorted only by column 1
	```

7. Use `kat gcp` to plot the gc content vs k-mer frequency.
	```bash
	mkfifo <named_pipe.fastq> && zcat <reads.fastq.gz> > <named_pipe.fastq> & # Make a named pipe and run in the background
	kat gcp -t 4 -o <output.gcp> <named_pipe.fastq> # Run KAT reading from the named pipe
	rm <named_pipe.fastq> # named pipes can only be used once, and so are removed after use.
	```
	Open the plot of GC vs coverage. On what scale is the GC content measured and how is this converted to GC%?

8. Use `kat comp` to compare `Bacteria/bacteria_R{1,2}.fastq.gz`.
	```bash
	mkfifo <named_pipe_read1.fastq> && zcat <read1.fastq.gz> > <named_pipe_read1.fastq> & # Make a named pipe for read 1 and run in background
	mkfifo <named_pipe_read2.fastq> && zcat <read2.fastq.gz> > <named_pipe_read2.fastq> & # Make a named pipe for read 2 and run in background
	kat comp -t 4 -o <output.cmp> --density_plot <named_pipe_read1.fastq> <named_pipe_read2.fastq> # run KAT on the named pipes and print a density plot
	kat plot spectra-mx -x 50 -y 500000 -i -o <output.cmp>-main.mx.spectra-mx.png <output.cmp>-main.mx # Make a spectra-mx plot
	rm <named_pipe_read1.fastq> <named_pipe_read2.fastq> # names pipes can only be used once, and so are removed after use
	```

	Why is there a difference in the distribution means between the two datasets?

9. Run Kraken on `Bacteria/bacteria_R{1,2}.fastq.gz`. What is the reason for this result? Can one do better?
	```bash
	KRAKEN_DB=/sw/courses/assembly/minikraken_20141208
	kraken --threads 4 --db $KRAKEN_DB --fastq-input --gzip-compressed --paired <read_{1,2}.fastq.gz> > <kraken.out>
	kraken-report --db $KRAKEN_DB <kraken.out> > <kraken.rpt>
	cut -f2,3 <kraken.out> > <krona.in>
	ktImportTaxonomy <krona.in> -o <krona.html>
	```
	
	**note**: `ktImportTaxonomy` is now a broken link. Use this file instead:
	```bash
	/sw/apps/bioinfo/Krona/2.7/src/KronaTools-2.7/scripts/ImportTaxonomy.pl
	```

10. Run Kraken on `Ecoli/E01_1_135x.fastq.gz`. What do you find here and how does the error rate influence this finding?
