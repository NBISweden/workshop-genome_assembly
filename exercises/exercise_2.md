---
title: De novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 22nd November 2018
---

# Data Quality Assessment

## Exercises

1. Change directory to the exercise data folder (`/home/data/byod/Sequence_Assessment`). Use `md5sum` to calculate the checksum of all the data files in the exercise folder. Redirect the checksum values to a file called **checksums.md5** in your working directory. Change directory to your working directory.
   
   <details>
   <summary> Solution - click to expand </summary>
   Simple solution:
   
   {% highlight bash %}
# Change to working directory
cd /home/data/byod/Sequence_Assessment/
# Check contents of folder
ls -R
# Use the **relative** paths of the files with md5sum and redirect STDOUT to a file in your working directory
md5sum */* > "$WORKDIR/checksums.md5"
   {% endhighlight %}
   
   Advanced solution (this is a more generally applicable solution):
   
   {% highlight bash %}
cd /home/data/byod/Sequence_Assessment/
# Get the **relative** paths of the files and use md5sum. Redirect the output to a file and screen
find -type "f" -exec md5sum {} \; | tee $WORKDIR/checksums.md5 # redirected output
   {% endhighlight %}
   
   </details>

2. Copy the exercise files to your working directory, but interrupt transfer with `ctrl + c`.

	{% highlight bash %}
cd $WORKDIR # Change directory to my working directory
cp -vr /home/data/byod/Sequence_Assessment/* . # Copy all the files to my working directory
	{% endhighlight %}

	Use the `-c` option of `md5sum` to check the files are complete.

	<details>
	<summary> Solution - click to expand </summary>

	{% highlight bash %}
md5sum -c checksums.md5
	{% endhighlight %}

	</details>

	Transfer the files again, this time making sure the files are complete.

3. Use `file` to check the fastq files. In what format is the data compressed?
	{% highlight bash %}
	file */*.fastq.gz
	{% endhighlight %}

	<details>
	<summary> Solution - click to expand </summary>

	All the fastq files are gzip compressed except one, which is in ascii text.

	</details>

4. One file is not compressed. Fix the problematic file by renaming the file and compressing it.

	<details>
	<summary> Solution - click to expand </summary>

	Bacteria/bacteria_Str2_R1.fastq.gz: ASCII text is the problematic file.

	{% highlight bash %}
	mv Bacteria/bacteria_Str2_R1.fastq.gz Bacteria/bacteria_Str2_R1.fastq  # Correct the filename to indicate it is not compressed
	gzip Bacteria/bacteria_Str2_R1.fastq                                   # Recompress the file.
	{% endhighlight %}

	</details>

5. Use `zcat` and `head` to view the first 20 lines of **Bacteria/bacteria_Str1_R1.fastq.gz**. What property of this data is
unusual for raw Illumina data?

	<details>
	<summary> Solution - click to expand </summary>

	{% highlight bash %}
	zcat Bacteria/bacteria_Str1_R1.fastq.gz | head -n 20
	{% endhighlight %}

	The property that is unusual is the read length. Unprocessed Illumina data contain fixed-length reads, i.e.
	every read is the same length. This data has been preprocessed.

	</details>

6. What does each tool in this command do?
	{% highlight bash %}
	zcat <fastq.gz> | seqtk seq -A - | grep -v “^>” | tr -dc “ACGTNacgtn” | wc -m
	{% endhighlight %}

	<details>
	<summary> Solution - click to expand </summary>

	{% highlight bash %}
	zcat <fastq.gz >     # concatenates compressed files to one output stream
	seqtk seq -A -       # seqtk is a toolkit for manipulating sequence data. The -A converts input to fasta output.
	grep -v "^>"         # grep searches for lines beginning (^) with the string > and excludes them (-v).
	tr -dc "ACGTNacgtn"  # tr translates characters from one set to another. The -dc deletes characters not in the "ACGTNacgtn" set.
	wc -m                # wc is the word count tool. wc -m counts characters.
	{% endhighlight %}

	</details>

7. How many bases in total are in these files?

    a. **Bacteria/bacteria_Str1_R{1,2}.fastq.gz**

    b. **SRR492065/SRR492065_{1,2}.fastq.gz**

    c. **E01/E01_1_135x.fastq.gz**

	<details>
	<summary> Solution - click to expand </summary>
	
	a. **Bacteria/bacteria_Str1_R{1,2}.fastq.gz**
	
	{% highlight bash %}
	zcat Bacteria/bacteria_Str1_R{1,2}.fastq.gz | seqtk seq -A - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m
	{% endhighlight %}
	
	225890464 (nucleotides)
	
	b. **SRR492065/SRR492065_{1,2}.fastq.gz**
	
	{% highlight bash %}
	zcat SRR492065/SRR492065_{1,2}.fastq.gz | seqtk seq -A - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m
	{% endhighlight %}
	
	1070871200 (nucleotides)
	
	c. **E01/E01_1_135x.fastq.gz**
	
	{% highlight bash %}
	zcat E01/E01_1_135x.fastq.gz | seqtk seq -A - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m
	{% endhighlight %}

	748508257 (nucleotides)

	</details>

10. How many bases in **E01/E01_1_135x.fastq.gz** are contained in reads 10kb or longer?

	<details>
	<summary> Solution - click to expand </summary>

	The `-L <int>` option in `seqtk` drops sequences smaller than `<int>` bases.

	{% highlight bash %}
	zcat E01/E01_1_135x.fastq.gz | seqtk seq -A -L 10000 - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m
	{% endhighlight %}

	510546313 (nucleotides)

	</details>

