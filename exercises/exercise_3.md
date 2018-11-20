---
title: De novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 22nd November 2018
---

# Data Quality Assessment

## Exercises

Nanoplot and FastQC produce html files as output. These can be opened on the login node using `firefox`,
however running graphical applications across a network is slow. Alternatively you can download the
html files to your computer using either `scp` or `rsync`. You can then open the downloaded files with your
own html browser.

In a terminal, on your computer (not the login node) type the following:

{% highlight bash %}
scp username@rackham.uppmax.uu.se:/path/to/files/filename /path/to/download
{% endhighlight %}

The command is similar using `rsync`:

{% highlight bash %}
rsync -av username@rackham.uppmax.uu.se:/path/to/files/filename /path/to/download
{% endhighlight %}

### Task 1.

Run NanoPlot on your PacBio data. The results can be opened using `firefox`.
What are the average and median length of the long reads.

{% highlight bash %}
source /proj/sllstore2017027/workshop-GA2018/tools/anaconda/miniconda2/etc/profile.d/conda.sh
conda activate GA2018
NanoPlot --help
{% endhighlight %}

<details>
<summary> Solution - click to expand </summary>

{% highlight bash %}
NanoPlot --fastq Ecoli_pacbio.fastq.gz
{% endhighlight %}

The average read length of the PacBio data is 8.6kb, and the median read length is 6.7kb.

</details>

### Task 2.

Run FastQC on your fastq data. Then open the results using `firefox` or `fastqc`.

{% highlight bash %}
module load bioinfo-tools FastQC/0.11.5
fastqc --help
{% endhighlight %}

How many sequences are in each fastq file?

<details>
<summary> Solution - click to expand </summary>
<div markdown="1">

{% highlight bash %}
fastqc -t 6 */*.fastq.gz
{% endhighlight %}

**Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz**: 5354356 each

**Escherichia_coli/ERR022075_{1,2}.fastq.gz**: 22720100 each

**Escherichia_coli/Ecoli_pacbio.fastq.gz**: 87225

</div>
</details>

### Task 3.

What is the average GC% in each data set?

<details>
<summary> Solution - click to expand </summary>
<div markdown="1">

**Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz**: 40%

**Escherichia_coli/ERR022075_{1,2}.fastq.gz**: 49%

**Escherichia_coli/Ecoli_pacbio.fastq.gz**: 49%

</div>
</details>

### Task 4.

Which quality score encoding is used?

<details>
<summary> Solution - click to expand </summary>

Sanger / Illumina 1.9

</details>

### Task 5.

What does a quality score of 20 (Q20) mean?

<details>
<summary> Solution - click to expand </summary>

An expectation of 1 error in 100bp.

</details>

### Task 6.

What does a quality score of 40 (Q40) mean?

<details>
<summary> Solution - click to expand </summary>

An expectation of 1 error in 10000bp.

</details>

### Task 7.

What distribution should the per base sequence plot follow?

<details>
<summary> Solution - click to expand </summary>

A Uniform distribution.

</details>

### Task 8.

What value should the per base GC distribution be centered on?

<details>
<summary> Solution - click to expand </summary>

Average GC content.

</details>

### Task 9.

How much duplication is present in each fastq file?

<details>
<summary> Solution - click to expand </summary>
<div markdown="1">

**Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz**: 29.4% and 17.24%

**Escherichia_coli/ERR022075_{1,2}.fastq.gz**: 61.71% and 27.87%

**Escherichia_coli/Ecoli_pacbio.fastq.gz**: 0.12% but this value is uninformative for pacbio due to the error rate.

</div>
</details>

### Task 10.

What is adapter read through?

<details>
<summary> Solution - click to expand </summary>

When the sequence reads past the insert into the adapter sequence on the other end.

</details>

### Task 11.

Let's look at the adapter sequence in the **Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz** fastq files. Illumina uses different adapters
for different libraries. It is important to know which adapter sequence it is. Since this is public data, it is sometimes difficult to
find out what the adapters were. Use `bbmerge` to discover the adapter sequence.

{% highlight bash %}
module load bioinfo-tools bbmap/38.08
bbmerge.sh --help
{% endhighlight %}

<details>
<summary> Solution - click to expand </summary>

{% highlight bash %}
bbmerge.sh in=Enterococcus_faecalis/SRR492065_1.fastq.gz in2=Enterococcus_faecalis/SRR492065_2.fastq.gz outa=adapters.fa
more adapters.fa
{% endhighlight %}

{% highlight bash %}
>Read1_adapter
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATNTCGTATGCCGTCTTNTGNTT
>Read2_adapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
{% endhighlight %}

</details>

### Task 12.

Use the command below to view the reads that have matching adapter sequence in your files.

{% highlight bash %}
paste <( zcat Enterococcus_faecalis/SRR492065_1.fastq.gz ) <( zcat Enterococcus_faecalis/SRR492065_2.fastq.gz ) \
| grep -A2 -B1 --colour=always "AGATCGGAAGAGC" | less -SR
{% endhighlight %}

In the next step we will use Trimmomatic to trim adapters. It needs the correct adapter file. Use
`grep` to identify the necessary adapter file to use. Trimmomatic's adapter files can be found in
`$TRIMMOMATIC_HOME/adapters/`.

{% highlight bash %}
module load bioinfo-tools trimmomatic/0.36
ls $TRIMMOMATIC_HOME
ls $TRIMMOMATIC_HOME/adapters
grep -B1 --colour=always "AGATCGGAAGAGCACACGTCTGAACTCC" $TRIMMOMATIC_HOME/adapters/*PE*.fa
grep -B1 --colour=always "AGATCGGAAGAGCGTCGTGTAGGGAAAG" $TRIMMOMATIC_HOME/adapters/*PE*.fa
{% endhighlight %}

Which adapter file should be used?

<details>
<summary> Solution - click to expand </summary>

{% highlight bash %}
$ grep -B1 --colour=always "AGATCGGAAGAGCACACGTCTGAACTCC" $TRIMMOMATIC_HOME/adapters/*PE*.fa
/sw/apps/bioinfo/trimmomatic/0.36/rackham/adapters/TruSeq3-PE-2.fa->PE2_rc
/sw/apps/bioinfo/trimmomatic/0.36/rackham/adapters/TruSeq3-PE-2.fa:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
$ grep -B1 --colour=always "AGATCGGAAGAGCGTCGTGTAGGGAAAG" $TRIMMOMATIC_HOME/adapters/*PE*.fa
/sw/apps/bioinfo/trimmomatic/0.36/rackham/adapters/TruSeq2-PE.fa->PCR_Primer1_rc
/sw/apps/bioinfo/trimmomatic/0.36/rackham/adapters/TruSeq2-PE.fa:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
--
/sw/apps/bioinfo/trimmomatic/0.36/rackham/adapters/TruSeq3-PE-2.fa->PE1_rc
/sw/apps/bioinfo/trimmomatic/0.36/rackham/adapters/TruSeq3-PE-2.fa:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
{% endhighlight %}

Therefore the file to use is: `/sw/apps/bioinfo/trimmomatic/0.36/rackham/adapters/TruSeq3-PE-2.fa`.

Since we have paired end data, we only use the files with PE in their name. Also as we are looking to remove adapter
read-through, we are searching for the reverse compliment of the adapter `*_rc`. These two sequences are only common
to one file, so we will use that one.

</details>

### Task 13.

Run Trimmomatic on **Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz** to only remove adapters. How many reads were trimmed for adapters?

{% highlight bash %}
java -jar $TRIMMOMATIC_HOME/trimmomatic-0.36.jar PE Enterococcus_faecalis/SRR492065_1.fastq.gz Enterococcus_faecalis/SRR492065_2.fastq.gz \
Enterococcus_faecalis/SRR492065_clean_paired_1.fastq.gz Enterococcus_faecalis/SRR492065_clean_unparied_1.fastq.gz \
Enterococcus_faecalis/SRR492065_clean_paired_2.fastq.gz Enterococcus_faecalis/SRR492065_clean_unparied_2.fastq.gz \
ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-PE-2.fa:2:30:10
{% endhighlight %}

<details>
<summary> Solution - click to expand </summary>

{% highlight bash %}
Input Read Pairs: 5354356 Both Surviving: 5339516 (99.72%) Forward Only Surviving: 11947 (0.22%) Reverse Only Surviving: 1888 (0.04%) Dropped: 1005 (0.02%)
{% endhighlight %}

</details>

### Task 14.

Using the bacterial database, run Kraken on **Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz**.
How many sequences are classified? What are they classified as?

{% highlight bash %}
module load bioinfo-tools Kraken2/2.0.7-beta-bc14b13 Krona/2.7
KRAKEN2DB=$SNIC_TMP/kraken_bacterial_db
rsync -av "/proj/sllstore2017027/workshop-GA2018/databases/kraken_bacterial_db/" "$KRAKEN2DB"
kraken2 --threads "$CPUS" --db "$KRAKEN2DB" --report "${PREFIX}_kraken.rpt" --gzip-compressed --paired "$READ1" "$READ2" > "${PREFIX}_kraken.tsv"
ktImportTaxonomy <( cut -f2,3 "${PREFIX}_kraken.tsv" ) -o "${PREFIX}_kraken_krona.html"
{% endhighlight %}

It takes a very long time to build the database, so we have provided a build for you above.

<details>
<summary> Click here to see how the database was built </summary>

This is how we built the database for you. See the Kraken2 homepage for how to build more comprehensive databases.

{% highlight bash %}
CPUS=${SLURM_NPROCS:-10}
TMPDB=$(mktemp -u)
kraken2-build --download-taxonomy --db "$SNIC_TMP/$TMPDB"
kraken2-build --download-library bacteria --db "$SNIC_TMP/$TMPDB"
kraken2-build --build --threads "$CPUS" --db "$SNIC_TMP/$TMPDB"
rsync -av "$SNIC_TMP/$TMPDB/" kraken_bacterial_db
{% endhighlight %}

</details>

<details>
<summary> Solution - click to expand </summary>
<div markdown="1">

{% highlight bash %}
CPUS=${SLURM_NPROCS:-10}
READ1=Enterococcus_faecalis/SRR492065_1.fastq.gz
READ2=Enterococcus_faecalis/SRR492065_2.fastq.gz
PREFIX=$(basename "$READ1" _1.fastq.gz )
KRAKEN2DB=$SNIC_TMP/kraken_bacterial_db
rsync -av "/proj/sllstore2017027/workshop-GA2018/databases/kraken_bacterial_db/" "$KRAKEN2DB"
kraken2 --threads "$CPUS" --db "$KRAKEN2DB" --report "${PREFIX}_kraken.rpt" --gzip-compressed --paired "$READ1" "$READ2" > "${PREFIX}_kraken.tsv"
ktImportTaxonomy <( cut -f2,3 "${PREFIX}_kraken.tsv" ) -o "${PREFIX}_kraken_krona.html"
{% endhighlight %}

{% highlight bash %}
5354356 sequences (1070.87 Mbp) processed in 42.233s (7606.9 Kseq/m, 1521.39 Mbp/m).
  4466191 sequences classified (83.41%)
  888165 sequences unclassified (16.59%)
{% endhighlight %}

To make an image, open the html file and click on the snapshot button. Then save the resulting image to `SRR492065_kraken_krona.svg`.

![A Krona plot of the Kraken analysis of SRR492065.]({{ site.url }}/workshop-genome_assembly/exercises/img/SRR492065_kraken_krona.svg)

The Kraken analysis shows at least three organisms in the sample; Enterococcus, Staphylococcus, and Cutibacterium. Enterococcus also shows a higher abundance than both Staphylococcus and Cutibacterium, which are both in similar proportions.

</div>
</details>

### Task 15.

Using the references, filter the reads that align to Staphylococcus and Cutibacterium.

{% highlight bash %}
/proj/sllstore2017027/workshop-GA2018/data/Illumina_SRR492065/Staphylococcus_aureus.fasta
/proj/sllstore2017027/workshop-GA2018/data/Illumina_SRR492065/Cutibacterium_avidum.fasta
{% endhighlight %}

{% highlight bash %}
# Load modules
module load bioinfo-tools bwa/0.7.17 samtools/1.8
# Align to the references.
PREFIX=$( basename "$REFERENCE" .fasta )    # Make a PREFIX from the reference name without .fasta on the end
bwa index "$REFERENCE"
bwa mem -t "$CPUS" "$REFERENCE" "$READ1" "$READ2" | samtools sort -@ "$CPUS" -T "$SNIC_TMP/$PREFIX" -O BAM -o "${PREFIX}_bwa_alignment.bam" -
samtools index "${PREFIX}_bwa_alignment.bam"
samtools flagstat "${PREFIX}_bwa_alignment.bam" > "${PREFIX}_bwa_alignment.bam.stats"
# Extract read names that align to reference
samtools view -@ "$CPUS" -F 4 "${PREFIX}_bwa_alignment.bam" | cut -f1 | sort -u -o "${PREFIX}_aligned_reads.tsv"
# Extract all read names
samtools view -@ "$CPUS" "${PREFIX}_bwa_alignment.bam" | cut -f1 | sort -u -o "${PREFIX}_all_reads.tsv"
# Extract read names unique to list1 between list1 (*.tsv) and list2 (*.tsv)
sort --parallel="$CPUS" "${LIST_1_TSV}" "${LIST_2_TSV}" "${LIST_2_TSV}" | uniq -u > "${LIST_3_TSV}"
# Use the list to extract the reads unique to list3
join -t ' ' <( zcat "$READ1" | paste - - - - | sort -k1,1 ) <( sed 's/^/@/' "${LIST_3_TSV}" ) | tr '\t' '\n' | pigz -c > "${LIST_3_TSV%%_*}_cleaned_R1.fastq.gz"
join -t ' ' <( zcat "$READ2" | paste - - - - | sort -k1,1 ) <( sed 's/^/@/' "${LIST_3_TSV}" ) | tr '\t' '\n' | pigz -c > "${LIST_3_TSV%%_*}_cleaned_R2.fastq.gz"
{% endhighlight %}

<details>
<summary> Solution - click to expand </summary>

{% highlight bash %}
cp /proj/sllstore2017027/workshop-GA2018/data/Illumina_SRR492065/Staphylococcus_aureus.fasta .
cp /proj/sllstore2017027/workshop-GA2018/data/Illumina_SRR492065/Cutibacterium_avidum.fasta .
CPUS=${SLURM_NPROCS:-10}
READ1=Enterococcus_faecalis/SRR492065_1.fastq.gz
READ2=Enterococcus_faecalis/SRR492065_2.fastq.gz
# Reads aligned to Staphylococcus
REFERENCE=Staphylococcus_aureus.fasta
PREFIX=$( basename "$REFERENCE" .fasta )
bwa index "$REFERENCE"
bwa mem -t "$CPUS" "$REFERENCE" "$READ1" "$READ2" | samtools sort -@ "$CPUS" -T "$SNIC_TMP/$PREFIX" -O BAM -o "${PREFIX}_bwa_alignment.bam" -
samtools index "${PREFIX}_bwa_alignment.bam"
samtools flagstat "${PREFIX}_bwa_alignment.bam" > "${PREFIX}_bwa_alignment.bam.stats"
# List of reads aligned to Staphylococcus
samtools view -@ "$CPUS" -F 4 "${PREFIX}_bwa_alignment.bam" | cut -f1 | sort -u -o "${PREFIX}_aligned_reads.tsv"
# List of all reads (not just aligned to Staphylococcus)
samtools view -@ "$CPUS" "${PREFIX}_bwa_alignment.bam" | cut -f1 | sort -u -o "${PREFIX}_all_reads.tsv"
# Reads aligned to Cutibacterium
REFERENCE=Cutibacterium_avidum.fasta
PREFIX=$( basename "$REFERENCE" .fasta )
bwa index "$REFERENCE"
bwa mem -t "$CPUS" "$REFERENCE" "$READ1" "$READ2" | samtools sort -@ "$CPUS" -T "$SNIC_TMP/$PREFIX" -O BAM -o "${PREFIX}_bwa_alignment.bam" -
samtools index "${PREFIX}_bwa_alignment.bam"
samtools flagstat "${PREFIX}_bwa_alignment.bam" > "${PREFIX}_bwa_alignment.bam.stats"
# List of reads aligned to Cutibacterium
samtools view -@ "$CPUS" -F 4 "${PREFIX}_bwa_alignment.bam" | cut -f1 | sort -u -o "${PREFIX}_aligned_reads.tsv"
# List of read names not mapped to Staphylococcus
LIST_1_TSV=Staphylococcus_aureus_all_reads.tsv
LIST_2_TSV=Staphylococcus_aureus_aligned_reads.tsv
sort --parallel="$CPUS" "${LIST_1_TSV}" "${LIST_2_TSV}" "${LIST_2_TSV}" | uniq -u > "Staphylococcus_aureus_removed_reads.tsv"
# List of read names not mapped to Cutibacterium
LIST_1_TSV=Staphylococcus_aureus_removed_reads.tsv
LIST_2_TSV=Cutibacterium_avidum_aligned_reads.tsv
sort --parallel="$CPUS" "${LIST_1_TSV}" "${LIST_2_TSV}" "${LIST_2_TSV}" | uniq -u > "SRR492065_Enterococcus_only_reads.tsv"
# Get filtered fastqs
LIST_3_TSV=SRR492065_Enterococcus_only_reads.tsv
join -t ' ' <( zcat "$READ1" | paste - - - - | sort -k1,1 ) <( sed 's/^/@/' "${LIST_3_TSV}" ) | tr '\t' '\n' | pigz -c > "${LIST_3_TSV%%_*}_cleaned_R1.fastq.gz"
join -t ' ' <( zcat "$READ2" | paste - - - - | sort -k1,1 ) <( sed 's/^/@/' "${LIST_3_TSV}" ) | tr '\t' '\n' | pigz -c > "${LIST_3_TSV%%_*}_cleaned_R2.fastq.gz"
{% endhighlight %}

</details>
