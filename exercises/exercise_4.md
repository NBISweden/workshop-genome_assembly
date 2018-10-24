---
title: De novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 22nd November 2018
---

# Data Quality Assessment

23. What is a k-mer?

	<details>
	<summary> Solution - click to expand </summary>

	A sequence of characters of length k.

	</details>

24. How many 4-mers are in the following sequence?
	```
	ACGTTTATCCTATACGGTAATAC
	```

	<details>
	<summary> Solution - click to expand </summary>

	20   (L-k+1)=(23-4+1)

	</details>

25. What are the frequencies of the 4-mers in the sequence above.

	<details>
	<summary> Solution - click to expand </summary>

	```
	ACGT 1	CGTT 1	GTTT 1	TTTA 1	TTAT 1
	TATC 1	ATCC 1	TCCT 1	CCTA 1	CTAT 1
	TATA 1	ATAC 2	TACG 1	ACGG 1	CGGT 1
	GGTA 1	GTAA 1	TAAT 1	AATA 1
	```

	</details>

26. Use the following commands to get a list of all the k-mers in **Bacteria/bacteria_Str1_R{1,2}.fastq.gz**.
	```bash
	# Make a named pipe and run in the background
	mkfifo <named_pipe.fastq> && zcat <reads.fastq.gz> > <named_pipe.fastq> &
	# Run KAT reading from the named pipe
	kat hist -t 4 -d -o <output.hist> <named_pipe.fastq>
	# Use Jellyfish to print out a human readable list
	kat_jellyfish dump <output.hist>-hash.jf27 > <kmer.lst>
	# named pipes can only be used once, so remove them after use.
	rm <named_pipe.fastq>
	```
	The **<kmer.lst>** file has the following format.
	```
	>frequency
	kmer_sequence
	```
	How many distinct k-mers were found? Use the line count command `wc -l` to find out.

 	<details>
	<summary> Solution - click to expand </summary>

	```bash
	mkfifo bacteria_Str1.fastq && zcat Bacteria/bacteria_Str1_R{1,2}.fastq.gz > bacteria_Str1.fastq &
	kat hist -t 4 -d -o bacteria_Str1.hist bacteria_Str1.fastq
	kat_jellyfish dump bacteria_Str1.hist-hash.jf27 > bacteria_Str1.kmer.lst
	rm bacteria_Str1.fastq
	# Inspect the output
	less bacteria_Str1.kmer.lst
	# Use paste to put 2 lines on 1 line, and then count lines
	paste - - < bacteria_Str1.kmer.lst | wc -l
	```

	There are 41531545 distinct k-mers.

	</details>

27. How many k-mers have a frequency of 1?  Use the following command to find out.
	```bash
	paste - - < kmer.lst | cut -c2- | awk '$1 == 1 { sum++ } END { print sum+0 }'
	# paste - - : reads two consecutive lines onto the same line.
	# cut -c2- : prints from the second character up to the last character in a line.
	# awk '$1 == 1 { sum++ } END { print sum+0 }' : if column 1 has a frequency of 1, increase the variable "sum". Print the value of "sum" at the end.
	```

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	paste - - < bacteria_Str1.kmer.lst | cut -c2- | awk '$1 == 1 { sum++ } END { print sum+0 }'
	```
	35701246 distinct k-mers have a frequency of 1.

	</details>

28. How many k-mers have a frequency greater than 5?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	paste - - < bacteria_Str1.kmer.lst | cut -c2- | awk '$1 > 5 { sum++ } END { print sum+0 }'
	```
	4969071 distinct k-mers have a frequency greater than 5.

	</details>

29. `kat hist` plotted a histogram in **<output.hist>.png**. Open this image using `display`. What is the estimated mean k-mer frequency (k-mer coverage)?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	display bacteria_Str1.hist.png
	```

	![A k-mer histogram of bacteria Str1.](Data_QC_Exercises_17-10-23/bacteria_Str1.hist.png)

	The peak is around 25 on the x-axis (k-mer frequency/coverage).

	</details>

30. The following command prints the frequency of each k-mer frequency between 5 and 45.  What is the mean k-mer frequency?
	```bash
	paste - - < kmer.lst | cut -c2- | awk '$1 > 5 && $1 < 45 {sum[$1]++ } END { for (freq in sum) {print freq" "sum[freq]} }' | sort -k1,1n | column -t
	# paste - - : reads two consecutive lines onto the same line.
	# cut -c2- : prints from the second character up to the last character in a line.
	# awk '$1 > 5 && $1 < 45 {sum[$1]++ } END { for (freq in sum) {print freq" "sum[freq]} }' :
	# 	if column 1 has a frequency greater than 5 and less than 45, increase the value of the array "sum[frequency]" by 1.
	# 	Then for each frequency in sum print the value of sum[frequency] at the end.
	# sort -k1,1n : Perform a numerical sort on the data sorted only by column 1
	# column -t : Format tabular data nicely on screen
	```

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	paste - - < bacteria_Str1.kmer.lst | cut -c2- | awk '$1 > 5 && $1 < 45 {sum[$1]++ } END { for (freq in sum) {print freq" "sum[freq]} }' | sort -k1,1n | column -t
	```

	The mean k-mer frequency is 26 with a count of 222001 distinct k-mers having a frequency of 26.

	</details>

31. Use the following command to plot the k-mer frequency vs gc content.
	```bash
	mkfifo <named_pipe.fastq> && zcat <reads.fastq.gz> > <named_pipe.fastq> & # Make a named pipe and run in the background
	kat gcp -t 4 -o <output.gcp> <named_pipe.fastq> # Run KAT reading from the named pipe
	rm <named_pipe.fastq> # named pipes can only be used once, and so are removed after use.
	```
	Open the plot of k-mer frequency (X-axis) vs GC (Y-axis). Why is the Y axis range from 0 to 27?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	mkfifo bacteria_Str1.fastq && zcat Bacteria/bacteria_Str1_R{1,2}.fastq.gz > bacteria_Str1.fastq &
	kat gcp -t 4 -o bacteria_Str1.gcp bacteria_Str1.fastq
	rm bacteria_Str1.fastq
	display bacteria_Str1.gcp.mx.png
	```

	![A GC content plot of bacteria Str1.](Data_QC_Exercises_17-10-23/bacteria_Str1.gcp.mx.png)

	The GC content scale (Y-axis) is the absolute GC count per k-mer. The default k-mer size is 27 and therefore the y-axis is from 0 to 27.

	</details>

32. Use `kat comp` to compare **Bacteria/bacteria_Str1_R{1,2}.fastq.gz**.
	```bash
	mkfifo <named_pipe_read1.fastq> && zcat <read1.fastq.gz> > <named_pipe_read1.fastq> & # Make a named pipe for read 1 and run in background
	mkfifo <named_pipe_read2.fastq> && zcat <read2.fastq.gz> > <named_pipe_read2.fastq> & # Make a named pipe for read 2 and run in background
	kat comp -t 4 -o <output.cmp> --density_plot <named_pipe_read1.fastq> <named_pipe_read2.fastq> # run KAT on the named pipes and print a density plot
	kat plot spectra-mx -x 50 -y 500000 -n -o <output.cmp>-main.mx.spectra-mx.png <output.cmp>-main.mx # Make a spectra-mx plot
	rm <named_pipe_read1.fastq> <named_pipe_read2.fastq> # names pipes can only be used once, and so are removed after use
	```

	Why is there a difference in the distribution means between the two datasets?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	# I need to make 2 named pipes; 1 for each input file to kat.
	mkfifo bacteria_Str1.R1.fastq && zcat Bacteria/bacteria_Str1_R1.fastq.gz > bacteria_Str1.R1.fastq &
	mkfifo bacteria_Str1.R2.fastq && zcat Bacteria/bacteria_Str1_R2.fastq.gz > bacteria_Str1.R2.fastq &
	kat comp -t 4 -o bacteria_Str1.R1vR2.cmp --density_plot bacteria_Str1.R1.fastq bacteria_Str1.R2.fastq
	kat plot spectra-mx -x 50 -y 500000 -n -o bacteria_Str1.R1vR2.cmp-main.mx.spectra-mx.png bacteria_Str1.R1vR2.cmp-main.mx
	rm bacteria_Str1.R1.fastq bacteria_Str1.R2.fastq
	display bacteria_Str1.R1vR2.cmp-main.mx.density.png
	display bacteria_Str1.R1vR2.cmp-main.mx.spectra-mx.png
	```

	![A density plot of R1 vs R2 in bacteria Str1.](Data_QC_Exercises_17-10-23/bacteria_Str1.R1vR2.cmp-main.mx.density.png)

	![A spectra-mx plot of R1 vs R2 in  bacteria Str1.](Data_QC_Exercises_17-10-23/bacteria_Str1.R1vR2.cmp-main.mx.spectra-mx.png)

	The spectra-mx plot shows the shared content for dataset 2 (R2) is shifted to the left and slightly higher than
	the shared content of dataset 1 (R1). From the previous FastQC analyses of these files we saw that R2 has lower read qualities than R1.
	Lower quality reads mean more errors and ambiguous bases, resulting in lower k-mer frequency counts (in R2) which shifts the mean
	k-mer frequency to the left.

	</details>

33.	**Bacteria_Str1** and **Bacteria_Str2** are from two different strains of the same bacteria. Use `kat comp` to compare k-mers from both reads to each other. What can be seen in the density plot?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	mkfifo bacteria_Str1.fastq && zcat Bacteria/bacteria_Str1*.fastq.gz > bacteria_Str1.fastq &
	mkfifo bacteria_Str2.fastq && zcat Bacteria/bacteria_Str2*.fastq.gz > bacteria_Str2.fastq &
	kat comp -t 4 -o bacteria_S1vS2.cmp --density_plot bacteria_Str1.fastq bacteria_Str2.fastq
	rm bacteria_Str1.fastq bacteria_Str2.fastq
	display bacteria_S1vS2.cmp-main.mx.density.png
	```

	![A density plot of bacterira Str1 vs Str2.](Data_QC_Exercises_17-10-23/bacteria_S1vS2.cmp-main.mx.density.png)

	The density plot shows a high frequency of k-mers along the x=0 axis and y=0 axis indicating there are significant quantities of k-mers
	unique to both libraries.

	</details>

34. Use `kat plot spectra-mx` to plot the shared and exclusive k-mers between the two datasets. How is what you see in the density plot shown?
	```bash
	kat plot spectra-mx -x 100 -y 300000 -n -o <Str1_vs_Str2.cmp>-main.mx.spectra-mx.png <Str1_vs_Str2.cmp>-main.mx
	```

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	kat plot spectra-mx -x 100 -y 300000 -n -o bacteria_S1vS2.cmp-main.mx.spectra-mx.png bacteria_S1vS2.cmp-main.mx
	display bacteria_S1vS2.cmp-main.mx.spectra-mx.png
	```

	![A spectra-mx plot of bacterira Str1 vs Str2.](Data_QC_Exercises_17-10-23/bacteria_S1vS2.cmp-main.mx.spectra-mx.png)

	The dataset 1 exclusive content and dataset 2 exclusive content have shallow peaks indicating sequence content unique to each dataset.

	</details>

35. Use `kat hist` and `kat gcp` on **SRR492065/SRR492065_{1,2}.fastq**. What can be inferred from these plots?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	mkfifo SRR492065.fastq && zcat SRR492065/SRR492065_*.fastq.gz > SRR492065.fastq &
	kat hist -t 4 -o SRR492065.hist SRR492065.fastq
	rm SRR492065.fastq
	mkfifo SRR492065.fastq && zcat SRR492065/SRR492065_*.fastq.gz > SRR492065.fastq &
	kat gcp -t 4 -o SRR492065.gcp SRR492065.fastq
	rm SRR492065.fastq
	display SRR492065.hist.png
	display SRR492065.gcp.mx.png
	```

	![A k-mer histogram of SRR492065.](Data_QC_Exercises_17-10-23/SRR492065.hist.png)

	![A GC content plot of SRR492065.](Data_QC_Exercises_17-10-23/SRR492065.gcp.mx.png)

	The histogram shows a very high trough around 50x k-mer coverage indicating a significant proportion of low frequency k-mers.

	The GC content plot shows three fairly distinct clusters with differing GC content indicating there are at least 3 different
	organisms in this sample. One organism is in much higher abundance than the other two.

	</details>

