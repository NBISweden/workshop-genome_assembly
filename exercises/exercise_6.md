---
title: De novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 23rd November 2018
---
# Assembly Validation Exercises.

## Introduction

The following exercises are intended to introduce you to the tools involved in assembly validation.

Many commands below are wrapped up in functions to make it easier for you to run.

Functions can be copy-and-pasted into your terminal window or you can write them in a file
and `source` the file (e.g. `source functions.sh`).

A function looks like the following:

{% highlight bash %}
run_many_commands () {
	ARG_1="$1"  # This is the first parameter
	ARG_2="$2"  # This is the second parameter
	ARG_3="$3"  # This is the third parameter
	# Then follows a complex series of commands to get a result
	command_1 "$ARG_1" "$ARG_2" > new.file
	command_2 new.file "$ARG_3" > result.file
}
{% endhighlight %}

Once this command is copied into your terminal, you can use the function to run the complex series of commands:

{% highlight bash %}
run_many_commands File1.ext File2.ext File3.ext
# Your results will now be in result.file as written in the function above
{% endhighlight %}

## Exercises

Some assemblies have been provided for you to assess. They can be found here:

{% highlight bash %}
{% endhighlight %}

This is how the assemblies were made:

Assembly | Description
:--- | :---
spades_k21-55_full.fasta | spades using k=21,33,55 and SRR492065_{1,2}.fastq.gz
spades_k21-127_full.fasta | spades using k=21,33,55,77,99,127 and SRR492065_{1,2}.fastq.gz
spades_k21-55_normalized.fasta | spades using k=21,33,55 and SRR492065_normalized_{1,2}.fastq.gz
spades_k21-127_normalized.fasta | spades using k=21,33,55,77,99,127 and SRR492065_normalized_{1,2}.fastq.gz
spades_k21-55_cleaned.fasta | spades using k=21,33,55 and SRR492065_cleaned_{1,2}.fastq.gz
spades_k21-127_cleaned.fasta | spades using k=21,33,55,77,99,127 and SRR492065_cleaned_{1,2}.fastq.gz
shovill_full_spades.fasta | shovill using assembler=spades and SRR492065_{1,2}.fastq.gz
shovill_full_megahit.fasta | shovill using assembler=megahit and SRR492065_{1,2}.fastq.gz
masurca_cleaned.fasta | masurca using SRR492065_cleaned_{1,2}.fastq.gz
abyss_k35_cleaned.fasta | abyss using k=35 and SRR492065_cleaned_{1,2}.fastq.gz

### Task 1.

QUAST is a good starting point to help evaluate the quality of assemblies. It provides many helpful contiguity statistics.

Run QUAST on all the assemblies at once and generate a report.

{% highlight bash %}
module load bioinfo-tools quast/4.5.4
quast.py --help
{% endhighlight %}

Which assembly looks the best?

<details>
<summary> Solution - click to expand </summary>

First run Quast on all the assemblies.

{% highlight bash %}
quast.py -t "${CPUS:-10}" --est-ref-size 3200000 *.fasta
{% endhighlight %}

![Quast Cumulative Length Plot](images/quast/cumulative_plot.png)

![Quast NGx Plot](images/quast/NGx_plot.png)

![Quast GC Plot](images/quast/GC_content_plot.png)

</details>

### Task 2.

Read congruency is an important measure in determining assembly accuracy. Clusters of read pairs that align incorrectly are
strong indicators of mis-assembly.

How well do the reads align back to the draft assemblies? Use `bwa` and `samtools` to assess the basic alignment statistics.

Make a folder for your results.

{% highlight bash %}
mkdir BWA
cd BWA
ln -s ../*.fasta . # link the fasta files in this directory
{% endhighlight %}

Then copy this function into your terminal.

{% highlight bash %}
align_reads () {
	ASSEMBLY="$1" # The assembly is the first parameter to this function. Must end in .fasta
	READ1="$2" # The first read pair is the second parameter to this function
	READ2="$3" # The second read pair is the third parameter to this function
	PREFIX=$(basename "$ASSEMBLY" .fasta )
	bwa index "$ASSEMBLY" # Index the assembly prior to alignment
	bwa mem -t "${CPUS:-10}" "$ASSEMBLY" "$READ1" "$READ2" | samtools sort -@ "${CPUS:-10}" -T "${PREFIX}" -O BAM -o "${PREFIX}_bwa_alignment.bam" -
	samtools index "${PREFIX}_bwa_alignment.bam"
}
{% endhighlight %}

To run the function above, copy the function into your terminal window and use in the following way:

{% highlight bash %}
align_reads assembly.fasta read1.fastq.gz read2.fastq.gz
{% endhighlight %}

This will then run `bwa index`, `bwa mem`, `samtools sort`, `samtools index`, and `samtools flagstat` in the correct
order and with the correct parameters.

<details>
<summary> Solution - click to expand </summary>

{% highlight bash %}
for FASTA in *.fasta; do
	align_reads "$FASTA" ../bacteria_R1.fastq.gz ../bacteria_R2.fastq.gz
done
{% endhighlight %}

{% highlight bash %}

{% endhighlight %}

</details>

### Task 3.

Use the alignments to search for signals of misassembly.

{% highlight bash %}
mkdir FRC
cd FRC
ln -s ../BWA/*.bam .
apply_FRC () {
	ALIGNMENT="$1" # The BAM alignment file is the first parameter to this function. The filename must end in .bam.
	GENOME_SIZE="$2" # The estimated genome size is the second parameter to this function
	PREFIX=$( basename "$ALIGNMENT" .bam )
	FRC --pe-sam "$ALIGNMENT" --genome-size "$GENOME_SIZE" --output "${PREFIX}"
}
{% endhighlight %}

Plot the FRC curves (`${PREFIX}_FRC.txt`) together in a plot using Gnuplot. Which assembly has the best feature response curve?

{% highlight bash %}
gnuplot << EOF
set terminal png size 800,600
set output 'FRC_Curve_all_assemblies.png'
set title "FRC Curve" font ",14"
set key right bottom font ",8"
set autoscale
set ylabel "Approximate Coverage (%)"
set xlabel "Feature Threshold"
files = system('find -name "*alignment_FRC.txt"')
plot for [data in files] data using 1:2 with lines title data
EOF
{% endhighlight %}

How do these results compare to the Quast results?

<details>
<summary> Solution - click to expand </summary>

{% highlight bash %}
for BAM in *.bam; do
	apply_FRC "$BAM" 5000000
done
gnuplot << EOF
set terminal png size 800,600
set output 'FRC_Curve_all_assemblies.png'
set title "FRC Curve" font ",14"
set key right bottom font ",8"
set autoscale
set ylabel "Approximate Coverage (%)"
set xlabel "Feature Threshold"
files = system('find -name "*alignment_FRC.txt"')
plot for [data in files] data using 1:2 with lines title data
EOF
{% endhighlight %}

![An FRC curve comparison of the assemblies](images/frc/FRC_Curve_all_assemblies.png)

</details>

### Task 4.

Use IGV to load an assembly, BAM file, and the corresponding FRC bed file.

Note that running a graphical application over the network is a very slow process. If you have
IGV installed on your computer, please download the files locally and view them there. Otherwise
IGV can be started on your node using:

{% highlight bash %}
module load bioinfo-tools IGV/2.4.2
igv.sh
{% endhighlight %}

### Task 5.

KAT is useful tool for high accuracy sequence data. The spectra-cn (copy number spectra) graph shows a
decomposition of k-mers in the assembly vs k-mers in the reads.
The black portion are k-mers not present in the assembly, the red portion is found once in the assembly, and so on.
This shows the completeness of an assembly, i.e. are all the reads assembled into contigs representative of the sequence data.

Use KAT to compare the assembly against the reads.

{% highlight bash %}
source /proj/sllstore2017027/workshop-GA2018/tools/anaconda/miniconda2/etc/profile.d/conda.sh
conda activate GA2018

apply_katcomp () {
	ASSEMBLY="$1"		# The assembly is the first parameter to this function. It must end in .fasta
	READ1="$2" 			# The first read pair is the second parameter to this function
	READ2="$3" 			# The second read pair is the third parameter to this function
	PREFIX=$( basename "${ASSEMBLY}" .fasta)
	TMP_FASTQ=$(mktemp -u --suffix ".fastq")
	mkfifo "${TMP_FASTQ}" && zcat "$READ1" "$READ2" > "${TMP_FASTQ}" &		# Make a named pipe and combine reads
	sleep 5																	# Give a little time for the pipe to be made
	kat comp -H 800000000 -t "${CPUS:-10}" -o "${PREFIX}_vs_reads.cmp" "${TMP_FASTQ}" "$ASSEMBLY" 	# Compare Reads to Assembly
	rm "${TMP_FASTQ}"
}
{% endhighlight %}

### Task 6.
