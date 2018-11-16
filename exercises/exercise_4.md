---
title: De novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 22nd November 2018
---

# Data Quality Assessment

## Exercises

### Task 1.

What is a k-mer?

<details>
<summary> Solution - click to expand </summary>

A sequence of characters of length k.

</details>

### Task 2.

How many 4-mers are in the following sequence?

```
ACGTTTATCCTATACGGTAATAC
```

<details>
<summary> Solution - click to expand </summary>

20   (L-k+1)=(23-4+1)

</details>

### Task 3.

What are the frequencies of the 4-mers in the sequence above.

<details>
<summary> Solution - click to expand </summary>

```
ACGT 1	CGTT 1	GTTT 1	TTTA 1	TTAT 1
TATC 1	ATCC 1	TCCT 1	CCTA 1	CTAT 1
TATA 1	ATAC 2	TACG 1	ACGG 1	CGGT 1
GGTA 1	GTAA 1	TAAT 1	AATA 1
```

</details>

### Task 4.

Use `kat hist` to obtain a histogram of data **Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz**.
What is the approximate k-mer coverage? Is anything wrong with this histogram?

{% highlight bash %}
# Make a named pipe and run in the background
TMPFILE=$(mktemp -u --suffix ".fastq")
mkfifo "$TMPFILE" && zcat "$READ1" "$READ2" > "$TMPFILE" &
# Run KAT reading from the named pipe
kat hist -t "$CPUS" -o "${PREFIX}.hist" "$TMPFILE"
rm "$TMPFILE"
{% endhighlight %}

<details>
<summary> Solution - click to expand </summary>

{% highlight bash %}
READ1=Enterococcus_faecalis/SRR492065_1.fastq.gz
READ2=Enterococcus_faecalis/SRR492065_2.fastq.gz
PREFIX=$(basename ${READ1%_1.*} )
TMPFILE=$(mktemp -u --suffix ".fastq")
mkfifo "$TMPFILE" && zcat "$READ1" "$READ2" > "$TMPFILE" &
kat hist -t "$CPUS" -o "${PREFIX}.hist" "$TMPFILE"
rm "$TMPFILE"
{% endhighlight %}

![A k-mer histogram of the SRR492065 data.]()

The homozygous peak in the histogram appears to be at 50x k-mer coverage.

The histogram is unusual because there is a higher than expected frequency of low frequency k-mers.

</details>

### Task 5.

Use the histogram to estimate the genome size from the homozygous peak.
Does this correspond to the estimated genome size of Enterococcus faecalis?

<details>
<summary> Solution - click to expand </summary>

The peak of the histogram is at around 50x coverage.

</details>

### Task 6.

Use `kat gcp` to plot the k-mer frequency vs gc content.
Why is the Y-axis from 0-27, and how can you estimate the average GC% of the blobs from this scale?

{% highlight bash %}
# Make a named pipe and run in the background
TMPFILE=$(mktemp -u --suffix ".fastq")
mkfifo "$TMPFILE" && zcat "$READ1" "$READ2" > "$TMPFILE" &
# Run KAT reading from the named pipe
kat gcp -t "$CPUS" -o "${PREFIX}.gcp" "$TMPFILE"
rm "$TMPFILE"
{% endhighlight %}

<details>
<summary> Solution - click to expand </summary>

{% highlight %}
READ1=Enterococcus_faecalis/SRR492065_1.fastq.gz
READ2=Enterococcus_faecalis/SRR492065_2.fastq.gz
PREFIX=$(basename ${READ1%_1.*} )
TMPFILE=$(mktemp -u --suffix ".fastq")
mkfifo "$TMPFILE" && zcat "$READ1" "$READ2" > "$TMPFILE" &
kat gcp -t "$CPUS" -o "${PREFIX}.gcp" "$TMPFILE"
rm "$TMPFILE"
{% endhighlight %}

![A GC content plot of SRR492065 data.]()

The GC content scale (Y-axis) is the absolute GC count per k-mer. The default k-mer size is 27 and therefore the y-axis is from 0 to 27.
One can estimate the average GC% of a blob from this scale by multiplying the approximate GC content value by 4 (since 4 * 27 == 108 =~ 100).

</details>

### Task 7.

Use `kat comp` to compare **Enterococcus_faecalis/SRR492065_{1,2}.fastq.gz**.

{% highlight bash %}
TMPFILE1=$(mktemp -u --suffix ".fastq")
TMPFILE2=$(mktemp -u --suffix ".fastq")
mkfifo "$TMPFILE1" && zcat "$READ1" > "$TMPFILE1" &
mkfifo "$TMPFILE2" && zcat "$READ2" > "$TMPFILE2" &
kat comp -t "$CPUS" -o "${PREFIX}_r1vr2.cmp" --density_plot "$TMPFILE1" "$TMPFILE2"
kat plot spectra-mx -x 50 -y 500000 -n -o "${PREFIX}_r1vr2.cmp-main.mx.spectra-mx.png" "${PREFIX}_r1vr2.cmp-main.mx"
rm "$TMPFILE1" "$TMPFILE2"
{% endhighlight %}

Why is there a difference in the distribution means between the two datasets?

<details>
<summary> Solution - click to expand </summary>

{% highlight bash %}
READ1=Enterococcus_faecalis/SRR492065_1.fastq.gz
READ2=Enterococcus_faecalis/SRR492065_2.fastq.gz
PREFIX=$(basename ${READ1%_1.*} )
TMPFILE1=$(mktemp -u --suffix ".fastq")
TMPFILE2=$(mktemp -u --suffix ".fastq")
mkfifo "$TMPFILE1" && zcat "$READ1" > "$TMPFILE1" &
mkfifo "$TMPFILE2" && zcat "$READ2" > "$TMPFILE2" &
kat comp -t "$CPUS" -o "${PREFIX}_r1vr2.cmp" --density_plot "$TMPFILE1" "$TMPFILE2"
kat plot spectra-mx -x 50 -y 500000 -n -o "${PREFIX}_r1vr2.cmp-main.mx.spectra-mx.png" "${PREFIX}_r1vr2.cmp-main.mx"
rm "$TMPFILE1" "$TMPFILE2"
{% endhighlight %}

![A density plot of READ1 vs READ2 in the SRR492065 data.]()

![A spectra-mx plot of READ1 vs READ2 in the SRR492065 data.]()

The spectra-mx plot shows the shared content for dataset 2 (READ2) is shifted to the left and slightly higher than
the shared content of dataset 1 (READ1). From the previous FastQC analyses of these files we saw that READ2 has lower read qualities than READ1.
Lower quality reads mean more errors and ambiguous bases, resulting in lower k-mer frequency counts (in READ2) which shifts the mean
k-mer frequency to the left.

</details>
