---
title: De Novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 23rd November 2018
---
# Assembly Validation Exercises.

### Task 1.

Run BUSCO on your assemblies.

{% highlight bash %}
module load bioinfo-tools BUSCO/3.0.2b
mkdir BUSCO
cd BUSCO
ln -s ../*.fasta .
# Use UPPMAX's script to setup the Busco environment.
source $BUSCO_SETUP

apply_BUSCO () {
	ASSEMBLY="$1" # The assembly is the first parameter to this function. The file must end in .fasta
	LINEAGE="$BUSCO_LINEAGE_SETS/bacteria_odb9"
	PREFIX=$( basename "$ASSEMBLY" .fasta )
	run_BUSCO.py -i "$ASSEMBLY" -l "$LINEAGE" -c "${CPUS:-10}" -m genome -o "${PREFIX}_busco"
}
{% endhighlight %}

What does the gene space look like for each assembly?

<details>
<summary> Solution - click to expand </summary>

{% highlight bash %}
::::::::::::::
run_abyss_k35_cleaned_busco/short_summary_abyss_k35_cleaned_busco.txt
::::::::::::::
	C:97.3%[S:48.0%,D:49.3%],F:0.0%,M:2.7%,n:148

	144	Complete BUSCOs (C)
	71	Complete and single-copy BUSCOs (S)
	73	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	4	Missing BUSCOs (M)
	148	Total BUSCO groups searched
::::::::::::::
run_masurca_cleaned_busco/short_summary_masurca_cleaned_busco.txt
::::::::::::::
	C:93.9%[S:50.0%,D:43.9%],F:0.0%,M:6.1%,n:148

	139	Complete BUSCOs (C)
	74	Complete and single-copy BUSCOs (S)
	65	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	9	Missing BUSCOs (M)
	148	Total BUSCO groups searched
::::::::::::::
run_shovill_full_megahit_busco/short_summary_shovill_full_megahit_busco.txt
::::::::::::::
	C:0.0%[S:0.0%,D:0.0%],F:0.0%,M:100.0%,n:148

	0	Complete BUSCOs (C)
	0	Complete and single-copy BUSCOs (S)
	0	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	148	Missing BUSCOs (M)
	148	Total BUSCO groups searched
::::::::::::::
run_spades_k21-127_cleaned_busco/short_summary_spades_k21-127_cleaned_busco.txt
::::::::::::::
	C:98.6%[S:35.1%,D:63.5%],F:0.0%,M:1.4%,n:148

	146	Complete BUSCOs (C)
	52	Complete and single-copy BUSCOs (S)
	94	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	2	Missing BUSCOs (M)
	148	Total BUSCO groups searched
::::::::::::::
run_spades_k21-127_full_busco/short_summary_spades_k21-127_full_busco.txt
::::::::::::::
	C:99.3%[S:13.5%,D:85.8%],F:0.0%,M:0.7%,n:148

	147	Complete BUSCOs (C)
	20	Complete and single-copy BUSCOs (S)
	127	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	1	Missing BUSCOs (M)
	148	Total BUSCO groups searched
::::::::::::::
run_spades_k21-127_normalized_busco/short_summary_spades_k21-127_normalized_busco.txt
::::::::::::::
	C:99.3%[S:14.2%,D:85.1%],F:0.0%,M:0.7%,n:148

	147	Complete BUSCOs (C)
	21	Complete and single-copy BUSCOs (S)
	126	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	1	Missing BUSCOs (M)
	148	Total BUSCO groups searched
::::::::::::::
run_spades_k21-55_cleaned_busco/short_summary_spades_k21-55_cleaned_busco.txt
::::::::::::::
	C:98.7%[S:36.5%,D:62.2%],F:0.0%,M:1.3%,n:148

	146	Complete BUSCOs (C)
	54	Complete and single-copy BUSCOs (S)
	92	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	2	Missing BUSCOs (M)
	148	Total BUSCO groups searched
::::::::::::::
run_spades_k21-55_full_busco/short_summary_spades_k21-55_full_busco.txt
::::::::::::::
	C:99.4%[S:12.2%,D:87.2%],F:0.0%,M:0.6%,n:148

	147	Complete BUSCOs (C)
	18	Complete and single-copy BUSCOs (S)
	129	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	1	Missing BUSCOs (M)
	148	Total BUSCO groups searched
::::::::::::::
run_spades_k21-55_normalized_busco/short_summary_spades_k21-55_normalized_busco.txt
::::::::::::::
	C:99.3%[S:14.2%,D:85.1%],F:0.0%,M:0.7%,n:148

	147	Complete BUSCOs (C)
	21	Complete and single-copy BUSCOs (S)
	126	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	1	Missing BUSCOs (M)
	148	Total BUSCO groups searched

{% endhighlight %}

</details>

### Task 2.

Use Mauve to compare the assemblies to the reference.

{% highlight bash %}
module load bioinfo-tools mauve/2015-02-13
Mauve
{% endhighlight %}

The reference is here:
{% highlight bash %}
/proj/sllstore2017027/workshop-GA2018/data/Illumina_SRR492065/Enterococcus_faecalis.fasta
{% endhighlight %}

Reorder the assemblies with respect to the reference (`Tools > Move contigs`), and then make an alignment (`align with ProgressiveMauve`).

Hint: Reordered are in an alignment folder. You can use `grep ">" assembly.fasta | less -S` to see if the contigs have been reordered.

### Task 3.

Let's change datasets now and re-circularize a different bacterial assembly.

The assembly is made from PacBio RSII data, and assembled using PacBio's HGAP assembler. The result is in fastq format.
Use `seqtk` to convert it to fasta in your directory.

{% highlight bash %}
module load bioinfo-tools seqtk/1.2-r101
seqtk seq -A /proj/sllstore2017027/workshop-GA2018/data/PacBio_P6C4_20kb_Ecoli/polished_assembly.fastq.gz > Ecoli_polished_assembly.fasta
{% endhighlight %}

Circular assemblies are written out as linear contigs with an overlap at the end to piece at the beginning.
Use Mummer to find the coordinates of the overlap.

{% highlight bash %}
conda activate GA2018
nucmer --maxmatch --nosimplify -p assembly Ecoli_polished_assembly.fasta Ecoli_polished_assembly.fasta
show-coords -lrcT assembly.delta | less -S
{% endhighlight %}

### Task 4.

We would like to start this assembly somewhere at the origin of replication, between the genes `rpmH` and `dnaA`.
Use Prokka to annotate the assembly and find a point to break the assembly near the origin of replication.

{% highlight bash %}
conda activate GA2018_prokka
prokka --cpus "${CPUS:-10}" --outdir Prokka_annotation Ecoli_polished_assembly.fasta
grep "rpmH\|dnaA" "$GFF_FILE"
{% endhighlight %}

Use `Bedtools` to check there are no other genes between this region. Use the start of the rpmH gene as `$START`
and the end of the dnaA gene as `$END`. The GFF file written by Prokka is not
strictly formatted as GFF and contains other data. The awk command retains only the needed lines of the file.
The output is all the genes in this region.

{% highlight bash %}
printf "unitig_0_quiver\t%d\t%d" $START $END | tee replication_origin_region.bed
bedtools intersect -wa -a <( awk 'NF > 3 || $0 ~ /^#/' "$GFF_FILE" ) -b replication_origin_region.bed
{% endhighlight %}

### Task 5.

Select a point between the genes (not within a gene) to use as a break point.
Use samtools to break the polished assembly at this point. Redirect the output
of both commands to a file called `Ecoli_broken.fasta`. Modify the commands
as necessary to get a section to overlap at the ends of the contigs.

{% highlight bash %}
# This selects all the sequence on 'unitig_0|quiver' starting from $BREAK until the end of the sequence
# The output is redirected ( > ) to Ecoli_broken.fasta
samtools faidx Ecoli_polished_assembly.fasta 'unitig_0|quiver':"$BREAK" > Ecoli_broken.fasta
# This selects all the sequence on 'unitig_0|quiver' starting from 1 until $BREAK (inclusive)
# The output is appended ( >> ) to Ecoli_broken.fasta
samtools faidx Ecoli_polished_assembly.fasta 'unitig_0|quiver':1-"$BREAK" >> Ecoli_broken.fasta
{% endhighlight %}

Merge the broken pieces again using AMOS.

{% highlight bash %}
conda activate GA2018
toAmos -s Ecoli_broken.fasta -o Ecoli_fixed.afg
minimus2 Ecoli_fixed
{% endhighlight %}

<details>
<summary> Solution - click to expand </summary>

The overlap shown in the previous task was near the end (4642500-4660550), but not up to it (4681865).
In order to make a successful reassembly on the overlap we need to trim out the part on the end that
does not overlap, by not including it in the selection.

{% highlight bash %}
samtools faidx Ecoli_polished_assembly.fasta 'unitig_0|quiver':1985200-4660550 > Ecoli_broken.fasta
samtools faidx Ecoli_polished_assembly.fasta 'unitig_0|quiver':1-1985200 >> Ecoli_broken.fasta
toAmos -s Ecoli_broken.fasta -o Ecoli_fixed.afg
minimus2 Ecoli_fixed
{% endhighlight %}

</details>

### Task 6.

Polish the assembly again to reduce errors in the overlap region. Use the reads in your Ecoli_pb.subreads.bam file.

Using the entire dataset takes a long time to run. For this example, subsample the reads using samtools to make the remaining
tools run quicker.

{% highlight bash %}
module load bioinfo-tools SMRT/5.0.1 samtools/1.8

# Subsample in the dataset in the interests of time (use the full dataset for a normal run)
samtools view -b -@ "${CPUS:-10}" -s 0.10 -o "${PREFIX}.subsampled.subreads.bam" "${PREFIX}.subreads.bam"

# Index the assembly for arrow
samtools faidx "$ASSEMBLY"

# Input files are "${PREFIX}.subsampled.subreads.bam" and "$ASSEMBLY"
# Output files are "${PREFIX}_pbalign_alignment.bam"
pbalign --nproc "${CPUS:-10}" --tmpDir "$SNIC_TMP" "${PREFIX}.subsampled.subreads.bam" "$ASSEMBLY" "${PREFIX}_pbalign_alignment.bam"

# Assembly is polished using Arrow.
arrow --numWorkers "${CPUS:-10}" "${PREFIX}_pbalign_alignment.bam" -r "$ASSEMBLY" -o "${PREFIX}_polished.fasta" -o "${PREFIX}_polished.fastq" -o "${PREFIX}_polished_changes.gff"
{% endhighlight %}
