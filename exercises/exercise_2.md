---
title: De novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 22nd November 2018
---

# Data Quality Assessment

## Exercises

### Task 1.

Change directory to the exercise data folder (`/proj/sllstore2017027/workshop-GA2018/data/QC_files`). Use `md5sum` to calculate the checksum of all the data files in the exercise folder. Redirect the checksum values to a file called **checksums.md5** in your working directory. Change directory to your working directory.

<details>
<summary> Solution - click to expand </summary>
Simple solution:

{% highlight bash %}
# Change to working directory
cd /proj/sllstore2017027/workshop-GA2018/data/QC_files
# Check contents of folder
ls -R
# Use the **relative** paths of the files with md5sum and redirect STDOUT to a file in your working directory
md5sum */* > "$WORKDIR/checksums.md5"
{% endhighlight %}

Advanced solution (this is a more generally applicable solution):

{% highlight bash %}
cd /proj/sllstore2017027/workshop-GA2018/data/QC_files
# Get the **relative** paths of the files and use md5sum. Redirect the output to a file and screen
find -type "f" -exec md5sum {} \; | tee $WORKDIR/checksums.md5 # redirected output
{% endhighlight %}

</details>

### Task 2.

Copy the exercise files to your working directory, but interrupt transfer with `ctrl + c`.

{% highlight bash %}
cd $WORKDIR # Change directory to my exercise directory
cp -vlr /proj/sllstore2017027/workshop-GA2018/data/QC_files/* . # Copy all the files to my exercise directory
{% endhighlight %}

Use the `-c` option of `md5sum` to check the files are complete.

<details>
<summary> Solution - click to expand </summary>

{% highlight bash %}
md5sum -c checksums.md5
{% endhighlight %}

</details>


Transfer the files again, this time making sure the files are complete.

### Task 3.

The PacBio data has been converted from RSII platforms' hd5 files to the Sequel platforms' unaligned BAM format using
the `bax2bam` tool in SMRT command line package. Use SMRT tools to extract the fastq from the BAM file.

{% highlight bash %}
module load bioinfo-tools SMRT/5.0.1
bam2fastq --help
{% endhighlight %}

<details>
<summary> Solution - click to expand </summary>

Only the subreads BAM file needs to be given as an argument. The scraps file contains poor quality sequence and adapters.

{% highlight bash %}
module load SMRT/5.0.1
bam2fastq -o Ecoli_pacbio Ecoli_pb.subreads.bam
{% endhighlight %}

</details>

### Task 4.

What does each tool in this command do?
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

### Task 5.

Load `seqtk` using:

{% highlight bash %}
module load bioinfo-tools seqtk/1.2-r101
{% endhighlight %}

How many bases in total are in these files?

  a. **Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz**:

  b. **Escherichia_coli/ERR022075_{1,2}.fastq.gz**:

  c. **Escherichia_coli/Ecoli_pacbio.fastq.gz**:

  d. **Escherichia_coli/Ecoli_nanopore.fasta**:

<details>
<summary> Solution - click to expand </summary>
<div markdown="1">

**Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz**

{% highlight bash %}
zcat Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz | seqtk seq -A - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m
{% endhighlight %}

1070871200 (nucleotides)

**Escherichia_coli/ERR022075_{1,2}.fastq.gz**

{% highlight bash %}
zcat Escherichia_coli/ERR022075_{1,2}.fastq.gz | seqtk seq -A - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m
{% endhighlight %}

4589460200 (nucleotides)

**Escherichia_coli/Ecoli_pacbio.fastq.gz**

{% highlight bash %}
zcat Escherichia_coli/Ecoli_pacbio.fastq.gz | seqtk seq -A - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m
{% endhighlight %}

748508361 (nucleotides)

**Escherichia_coli/Ecoli_nanopore.fasta**

{% highlight bash %}
grep -v "^>" Escherichia_coli/Ecoli_nanopore.fasta | tr -dc "ACGTNacgtn" | wc -m
{% endhighlight %}

410782292 (nucleotides)

</div>
</details>

### Task 6.

How many bases in **Escherichia_coli/Ecoli_pacbio.fastq.gz** are contained in reads 10kb or longer?

<details>
<summary> Solution - click to expand </summary>
<div markdown="1">

The `-L <int>` option in `seqtk` drops sequences smaller than `<int>` bases.

{% highlight bash %}
zcat Escherichia_coli/Ecoli_pacbio.fastq.gz | seqtk seq -A -L 10000 - | grep -v "^>" | tr -dc "ACGTNacgtn" | wc -m
{% endhighlight %}

510546352 (nucleotides)

</div>
</details>

### Task 7.

What is the depth of coverage of these datasets?

a. **Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz**:

b. **Escherichia_coli/ERR022075_{1,2}.fastq.gz**:

c. **Escherichia_coli/Ecoli_pacbio.fastq.gz**:

d. **Escherichia_coli/Ecoli_nanopore.fasta**:

<details>
<summary> Solution - click to expand </summary>
<div markdown="1">

**Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz**

Searching for the Enterococcus faecalis genome size gives and approximate value of 3.22 Mb.

{% highlight bash %}
echo "1070871200 / 3220000" | bc -l
{% endhighlight %}

Approximately 332x depth of coverage

**Escherichia_coli/ERR022075_{1,2}.fastq.gz**

Searching for the Escherichia coli genome size gives and approximate value of 4.6 Mb.

{% highlight bash %}
echo "4589460200 / 4600000" | bc -l
{% endhighlight %}

Approximately 998x depth of coverage.

**Escherichia_coli/Ecoli_pacbio.fastq.gz**

{% highlight bash %}
echo "748508361 / 4600000" | bc -l
{% endhighlight %}

Approximately 163x depth of coverage.

**Escherichia_coli/Ecoli_nanopore.fasta**

{% highlight bash %}
echo "410782292 / 4600000" | bc -l
{% endhighlight %}

Approximately 89x depth of coverage.

</div>
</details>

### Task 8.

Use `seqtk` to subsample **Escherichia_coli/ERR022075_{1,2}.fastq.gz** to approximately 100x coverage.

<details>
<summary> Solution - click to expand </summary>

Since we want approximately 10% of the reads, we use a value of 0.1 as the fraction of reads to sample.

{% highlight bash %}
seqtk sample -s100 Escherichia_coli/ERR022075_1.fastq.gz 0.1 > Escherichia_coli/ERR022075_100x_1.fastq.gz
seqtk sample -s100 Escherichia_coli/ERR022075_2.fastq.gz 0.1 > Escherichia_coli/ERR022075_100x_2.fastq.gz
{% endhighlight %}

</details>


### Task 9.

Use `bbmap` to normalize the **Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz** data.

{% highlight bash %}
module load bioinfo-tools bbmap/38.08
bbnorm.sh --help
{% endhighlight %}

<details>
<summary> Solution - click to expand </summary>

Since we want approximately 10% of the reads, we use a value of 0.1 as the fraction of reads to sample.

{% highlight bash %}
bbnorm.sh in=Enterococcus_faecalis/SRR492065_1.fastq.gz in2=Enterococcus_faecalis/SRR492065_2.fastq.gz \
 out=Enterococcus_faecalis/SRR492065_normalized_1.fastq.gz out2=Enterococcus_faecalis/SRR492065_normalized_2.fastq.gz \
 target=100 min=5
{% endhighlight %}

</details>
