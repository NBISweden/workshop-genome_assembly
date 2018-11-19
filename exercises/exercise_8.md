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

### Task 2.

Use Mauve to compare the assemblies to the reference.

The reference is here:
{% highlight bash %}
/proj/sllstore2017027/workshop-GA2018/data/Illumina_SRR492065/Enterococcus_faecalis.fasta
{% endhighlight %}

Reorder the assemblies with respect to the reference (`Tools > Move contigs`), and then make an alignment (`align with ProgressiveMauve`).

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

Use `Bedtools` to check there are no other genes between this region. The GFF file written by Prokka is not
strictly formatted as GFF and contains other data. The awk command retains only the needed lines of the file.

{% highlight bash %}
printf "unitig_0_quiver\t%d\t%d" $START $END > replication_origin_region.bed
bedtools intersect -wa -a <( awk 'NF > 3 || $0 ~ /^#/' "$GFF_FILE" ) -b replication_origin_region.bed
{% endhighlight %}

### Task 5.

Select a point between the genes to use as a break point. Use samtools to break the polished assembly at this point.

{% highlight bash %}
samtools faidx 'unitig_0|quiver':"$BREAK" > Ecoli_broken.fasta
samtools faidx 'unitig_0|quiver':1-"$BREAK" >> Ecoli_broken.fasta
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

{% highlight bash %}
module load bioinfo-tools SMRT/5.0.1
# Input files are "${PREFIX}.subreads.bam" and "$ASSEMBLY"
# Output files are "${PREFIX}_pbalign_alignment.bam"
pbalign --nproc "${CPUS:-10}" --tmpDir "$SNIC_TMP" "${PREFIX}.subreads.bam" "$ASSEMBLY" "${PREFIX}_pbalign_alignment.bam"

# Assembly is polished using Arrow.
arrow --numWorkers "${CPUS:-10}" "${PREFIX}_pbalign_alignment.bam" -r "$ASSEMBLY" -o "${PREFIX}_polished.fasta" -o "${PREFIX}_polished.fastq" -o "${PREFIX}_polished.gff"
{% endhighlight %}

Optional:
Long read correction can still leave errors. Try using Pilon to see if that reduces the error further.
