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

Look at the workshop wiki for a brief description of the GFA file format, and the
[Bandage webpage](https://github.com/rrwick/Bandage/wiki/CSV-labels) for information on how to construct the csv.

Hint: Use the unix command line tools such as `grep`, sort, `cut`, and `join` to manipulate the data into csv.

<details>
<summary> Solution - click to expand </summary>

Bandage displays the segment lines (S) of a GFA file. The second column of a segment line is the node name.
Spades contigs correspond to the path lines (P).
This means we need to know two things. First, what are the spades contigs annotated as, and second, which
paths correspond to which segments.

Let's use the data for spades_k21-55_full as the example.

The aim is to combine the data from the GFA and Blast into one file based on a common feature, the contig names.
This means this is a perfect task for `join`. However all the data is not in a format friendly for join, so let's
go through how to manipulate the data into two files that join can work with.

Starting with the Blast file, first we need column 1 which are the contig names.
Then let's use column 15 of the blast output to annotate the contigs, and use sort to remove duplicate entries.

{% highlight bash %}
cut -f1,15 spades_k21-55_full_blast_alignment.tsv | sort -u > spades_k21-55_full_blast_annotation.tsv
head spades_k21-55_full_blast_annotation.tsv
{% endhighlight %}

These are the first 10 lines of that output, to show you what the data should look like.

{% highlight bash %}
NODE_1000_length_306_cov_2.509960	Staphylococcus epidermidis
NODE_1000_length_306_cov_2.509960	Staphylococcus epidermidis ATCC 12228
NODE_1000_length_306_cov_2.509960	Staphylococcus epidermidis RP62A
NODE_1001_length_305_cov_4.088000	Staphylococcus epidermidis
NODE_1001_length_305_cov_4.088000	Staphylococcus epidermidis ATCC 12228
NODE_1001_length_305_cov_4.088000	Staphylococcus epidermidis RP62A
NODE_1002_length_305_cov_3.296000	Staphylococcus epidermidis
NODE_1002_length_305_cov_3.296000	Staphylococcus epidermidis ATCC 12228
NODE_1002_length_305_cov_3.296000	Staphylococcus epidermidis PM221
NODE_1003_length_305_cov_2.700000	Staphylococcus hominis
{% endhighlight %}

We can see already that even though we sorted and removed duplicate lines, there are still multiple entries of
a species for certain contigs. So, let's change the command above to only keep the first instance of every contig name.

{% highlight bash %}
cut -f1,15 spades_k21-55_full_blast_alignment.tsv | sort -u -k1,1 > spades_k21-55_full_blast_annotation.tsv
{% endhighlight %}

This now leaves us with one label per contig name.

{% highlight bash %}
NODE_1000_length_306_cov_2.509960	Staphylococcus epidermidis
NODE_1001_length_305_cov_4.088000	Staphylococcus epidermidis
NODE_1002_length_305_cov_3.296000	Staphylococcus epidermidis
NODE_1003_length_305_cov_2.700000	Staphylococcus hominis
NODE_1004_length_304_cov_39.807229	Paenibacillus sp. FSL R7-0331
NODE_1005_length_304_cov_2.815261	Staphylococcus aureus
NODE_1006_length_304_cov_2.112450	Staphylococcus hominis
NODE_1007_length_303_cov_3.193548	Staphylococcus epidermidis
NODE_1008_length_303_cov_2.286290	Staphylococcus hominis
NODE_1009_length_303_cov_1.326613	Staphylococcus hominis
{% endhighlight %}

Next we want the path lines of the GFA file, and specifically the information from column 2 and 3, which are the
contig names, and the segment(s) linked to that name.

{% highlight bash %}
grep "^P" spades_k21-55_full.gfa | cut -f2,3 | head
{% endhighlight %}

In order to be able to merge this dataset with the one above, the contig names need to be identical, however all the names now
end with an underscore followed by a number. Furthermore, the contig name belongs to many segment names.

{% highlight bash %}
NODE_1_length_1448318_cov_80.178312_1	10311124+,10031624+,10311308+,9911828-,10258069-,10057802+,10272111+,10118480+,9985602+,10086712+,1392052-,9996864-,10171502-,10320089-,10273179-,10319496-,10273179-,10319498-,9886408-,10074530-,10294167-,10135462-,9886408-,10294169-,10294167-,10178112-,10128536+,10309890+,10128536+,10317320-,723032-,10317322-,723032-,10306893-,10226805-,10245575-,10226805-,10319919-,10318590+,10320057+,9955186-,10319875-,9711424-,10320055+,9955186-,10319933-,10318590+,10315156-,9711424-,10319640-,10302373+,10317464+,10314624+,10314630+,10302021+,10135030+,10310524+,10310194+,10310388+,235774+,10313776+,10300185+,10314676+,10314684+,10297987+,10313836+,10316434-,10017942-,10317812-,10314202-,1548974+,10309092+,10314132+,10314140+,10313694+,10314684-,10314676-,10318362-,10314186+,10314194+,10314202+,10317812+,10312142-,10218043+,10313772-,10310880+,9796758-,10311458+,10171996-,10270995+,666778+,10307503+,1936534-,10312728+,9906056+,10312142+,10317812-,10314202-,1548974+,10313566+,10287551+,10313866+,10297503+,10318937+,10085614+,10311126+,10085614+,10310820-,10200554-,10115886-,10085372-,10313984-,10314700-,10314692-,10314684-,10314676-,10300185-,10313832+,9946622+,10210134+,9946622+,10315032-,10311484+,10311344-,2957484+,10317710+,2957484+,10317708+,10171502+,10248517-,1392052+,10184644+,9985602-,10186332-,10272111-,10204768+,10258069+,10068726+,10311310-,10031624-,10320101-,10200554-,10186842+,10085372-,10287176-,10319791-,10313030-,10048488-,10110418-,9804040+,10309472-
NODE_2_length_336381_cov_10.156140_1	10285711+,10243057+,10268313+,10243057+,10287871-
NODE_2_length_336381_cov_10.156140_2	10320938+
NODE_3_length_320122_cov_11.770314_1	10247199-,10281667+
NODE_3_length_320122_cov_11.770314_2	10289913-
NODE_3_length_320122_cov_11.770314_3	10287919-
NODE_3_length_320122_cov_11.770314_4	10299737+,1484104+,10291809+
NODE_3_length_320122_cov_11.770314_5	10321250+,328964-,2720874-
NODE_4_length_310123_cov_12.870461_1	10291757+
NODE_4_length_310123_cov_12.870461_2	10300875+
{% endhighlight %}

To break this problem down, first let's take only the first segment name as the node to annotate. We can
use the `,` as a delimiter, which means everything before the first column is contained in column 1.

{% highlight bash %}
grep "^P" spades_k21-55_full.gfa | cut -f2,3 | cut -f1 -d, | head
{% endhighlight %}

This leaves us with data of this form:

{% highlight bash %}
NODE_1_length_1448318_cov_80.178312_1	10311124+
NODE_2_length_336381_cov_10.156140_1	10285711+
NODE_2_length_336381_cov_10.156140_2	10320938+
NODE_3_length_320122_cov_11.770314_1	10247199-
NODE_3_length_320122_cov_11.770314_2	10289913-
NODE_3_length_320122_cov_11.770314_3	10287919-
NODE_3_length_320122_cov_11.770314_4	10299737+
NODE_3_length_320122_cov_11.770314_5	10321250+
NODE_4_length_310123_cov_12.870461_1	10291757+
NODE_4_length_310123_cov_12.870461_2	10300875+
{% endhighlight %}

If we use `_` as the cut delimiter, rather than `\t` (tab) we can cut off that last underscore and number. The problem
is then that we also cut off the node names. Swapping the positions of the node names and contig names would prevent
that by making the node names part of column 1 instead of column 7.

{% highlight bash %}
grep "^P" spades_k21-55_full.gfa | cut -f2,3 | cut -f1 -d, | awk '{ print $2 "\t" $1 }' | cut -f1-6 -d"_" | head
{% endhighlight %}

The data now looks like this:

{% highlight bash %}
10311124+	NODE_1_length_1448318_cov_80.178312
10285711+	NODE_2_length_336381_cov_10.156140
10320938+	NODE_2_length_336381_cov_10.156140
10247199-	NODE_3_length_320122_cov_11.770314
10289913-	NODE_3_length_320122_cov_11.770314
10287919-	NODE_3_length_320122_cov_11.770314
10299737+	NODE_3_length_320122_cov_11.770314
10321250+	NODE_3_length_320122_cov_11.770314
10291757+	NODE_4_length_310123_cov_12.870461
10300875+	NODE_4_length_310123_cov_12.870461
{% endhighlight %}

Let's save the full output to a file:

{% highlight bash %}
grep "^P" spades_k21-55_full.gfa | cut -f2,3 | cut -f1 -d, | awk '{ print $2 "\t" $1 }' | cut -f1-6 -d"_" > spades_k21-55_full_node_contig_names.tsv
{% endhighlight %}

So now we can try and merge the files using join. However, join needs data to be sorted on the column you want to merge on.
Also we need to tell it that we want it to merge on column 1 of the filtered down blast file, and column 2 of the path segment
relationship file. I also tell join that I want it to separate only on the `\t` (tab) character, and that I want the output
in the order of column 1 from file 1 (node names), column 2 of file 1 (contig names), followed by column 2 of file 2 (blast annotation).

{% highlight bash %}
join -1 2 -2 1 <( sort -k2,2 spades_k21-55_full_node_contig_names.tsv ) <( sort -k1,1 spades_k21-55_full_blast_annotation.tsv) -t $'\t' -o 1.1,1.2,2.2 | head
{% endhighlight %}

The output now looks like this:

{% highlight bash %}
1511768+	NODE_1000_length_306_cov_2.509960	Staphylococcus epidermidis
9815238+	NODE_1001_length_305_cov_4.088000	Staphylococcus epidermidis
10091228+	NODE_1002_length_305_cov_3.296000	Staphylococcus epidermidis
1099986+	NODE_1003_length_305_cov_2.700000	Staphylococcus hominis
10320035+	NODE_1004_length_304_cov_39.807229	Paenibacillus sp. FSL R7-0331
10115802+	NODE_1005_length_304_cov_2.815261	Staphylococcus aureus
9844216+	NODE_1006_length_304_cov_2.112450	Staphylococcus hominis
10271579+	NODE_1007_length_303_cov_3.193548	Staphylococcus epidermidis
10031180+	NODE_1008_length_303_cov_2.286290	Staphylococcus hominis
9908708+	NODE_1009_length_303_cov_1.326613	Staphylococcus hominis
{% endhighlight %}

This now needs to be converted to a CSV file, so the `\t` (tab) characters need to be changed to `,`.

{% highlight bash %}
join -1 2 -2 1 <( sort -k2,2 spades_k21-55_full_node_contig_names.tsv ) <( sort -k1,1 spades_k21-55_full_blast_annotation.tsv) -t $'\t' -o 1.1,1.2,2.2 | tr "\t" "," > spades_k21-55_full_bandage_labels.csv
head spades_k21-55_full_bandage_labels.csv
{% endhighlight %}

Now we have a CSV file that can be loaded into bandage.

{% highlight bash %}
1511768+,NODE_1000_length_306_cov_2.509960,Staphylococcus epidermidis
9815238+,NODE_1001_length_305_cov_4.088000,Staphylococcus epidermidis
10091228+,NODE_1002_length_305_cov_3.296000,Staphylococcus epidermidis
1099986+,NODE_1003_length_305_cov_2.700000,Staphylococcus hominis
10320035+,NODE_1004_length_304_cov_39.807229,Paenibacillus sp. FSL R7-0331
10115802+,NODE_1005_length_304_cov_2.815261,Staphylococcus aureus
9844216+,NODE_1006_length_304_cov_2.112450,Staphylococcus hominis
10271579+,NODE_1007_length_303_cov_3.193548,Staphylococcus epidermidis
10031180+,NODE_1008_length_303_cov_2.286290,Staphylococcus hominis
9908708+,NODE_1009_length_303_cov_1.326613,Staphylococcus hominis
{% endhighlight %}

After loading the network into bandage, we can see that the largest network is a mix of Enterococcus, Staphylococcus, and Cutibacterium.

{% highlight bash %}
Bandage load spades_k21-55_full.gfa --draw
{% endhighlight %}

</details>
