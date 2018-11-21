---
title: De novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 23rd November 2018
---
# Assembly Validation Exercises.

### Task 1.

Contamination is not always easy to detect at the read level. Let's assess the assemblies again using Kraken2, but this time
at the assembly level.

{% highlight bash %}
module load bioinfo-tools Kraken2/2.0.7-beta-bc14b13 Krona/2.7
mkdir Kraken
cd Kraken
ln -s ../*.fasta .
KRAKEN2DB=$SNIC_TMP/kraken_bacterial_db
rsync -av "/proj/sllstore2017027/workshop-GA2018/databases/kraken_bacterial_db/" "$KRAKEN2DB"

apply_kraken () {
   ASSEMBLY="$1" # The assembly is the first parameter to this function. This file must end in .fasta
   PREFIX=$( basename "$ASSEMBLY" .fasta )
   echo "Running Kraken2: $ASSEMBLY"
   kraken2 --threads "${CPUS:-10}" --db "$KRAKEN2DB" --report "${PREFIX}_kraken.rpt" "$ASSEMBLY" > "${PREFIX}_kraken.tsv"
   ktImportTaxonomy <( cut -f2,3 "${PREFIX}_kraken.tsv" ) -o "${PREFIX}_kraken_krona.html"
}
{% endhighlight %}

### Task 2.

One could also classify the sequences using Blast, however this is a time consuming process.
Running a single assembly on the resources allocated to you can take up to an hour.

Copy the results we have already provided and inspect them.

{% highlight bash %}
/proj/sllstore2017027/workshop-GA2018/data/SRR492065_blast/
{% endhighlight %}

The results were generated in the following way:

{% highlight bash %}
module load bioinfo-tools blast/2.7.1+ Krona/2.7
mkdir Blast
cd Blast
ln -s ../*.fasta .
set_difference () {
	sort "$1" "$2" "$2" | uniq -u
}
apply_Blast () {
    ASSEMBLY="$1" # The assembly is the first parameter to this function. The file must end in .fasta
    PREFIX=$( basename "$ASSEMBLY" .fasta )
    echo "Blast: $ASSEMBLY"
    blastn -task megablast -query "$ASSEMBLY" -db "${BLASTDB:-/sw/data/uppnex/blast_databases}/nt" \
        -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
        -culling_limit 5 -num_threads "${CPUS:-10}" -evalue 1e-25 -out "${PREFIX}_blast_alignment.tsv"
    # Blast tabular output does not include `no hit`. The set difference is used to include remaining unclassified sequence.
    ktImportTaxonomy <( cat <( cut -f1,2 "${PREFIX}_blast_alignment.tsv" | sort -u ) <( set_difference <( grep ">" "$ASSEMBLY" | cut -c2- ) <(cut -f1 "${PREFIX}_blast_alignment.tsv" ) )) -o "${PREFIX}_blast_krona.html"
}
{% endhighlight %}

### Task 3.  

Run Blobtools on each assembly. Blobtools requires both a BAM file as input and blast output for the classification step.
What do these plots show?

{% highlight bash %}
module load bioinfo-tools blobtools/1.0-af55ced
mkdir Blobtools
cd Blobtools
ln -s ../*.fasta .
ln -s ../BWA/*.bam .
ln -s ../Blast/*.tsv .
apply_blobtools () {
    ASSEMBLY="$1" # The assembly is the first parameter to this function. The file must end in .fasta
    BAM="$2" # The BAM file is the second parameter to this function
    BLAST="$3" # The BLAST file is the third parameter to this function
    PREFIX=$( basename "$ASSEMBLY" .fasta )
    blobtools create -i "$ASSEMBLY" -b "$BAM" -t "$BLAST" -o "${PREFIX}_blobtools"
    blobtools blobplot -i "${PREFIX}_blobtools.blobDB.json" -o "${PREFIX}_blobtools"
}
{% endhighlight %}

### Task 4.

Bandage is a great tool to visualise assembly graphs. Load and draw some of the assembly graphs.

The assembly graphs for some of the assemblies can be found here:
{% highlight bash %}
/proj/sllstore2017027/workshop-GA2018/data/SRR492065_graphs
{% endhighlight %}

Bandage is loaded through the conda environment `GA2018`.

{% highlight bash %}
conda activate GA2018
Bandage --help
{% endhighlight %}

### Task 5.

Use the Kraken and Blast results to create a label csv to load into bandage and identify scaffolds.

Look at the workshop wiki for a brief description of the GFA file format, and the Bandage webpage for
information on how to construct the csv.

Hint: Use the unix command line tools such as `grep`, `cut`, and `join` to manipulate the data into csv.
