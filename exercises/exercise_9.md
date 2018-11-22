---
title: De novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 24th November 2018
---
# Ecoli genome assembly project

## Exercise

This is a team project, so split up the workload as you see fit.

You have three datasets from the Ecoli K12 substrain MG1655, sequenced using Illumina, PacBio, and Nanopore.

{% highlight bash %}
/proj/sllstore2017027/workshop-GA2018/data/QC_files/Escherichia_coli
{% endhighlight %}

The aim is for you to try and explore different assemblers and see what you get. It will be impossible to
do evaluate every combination, so choose your tasks wisely. Document your commands and share them with
each other.

Working directories have been created for each team:
{% highlight bash %}
/proj/sllstore2017027/nobackup_GA2018/team_Turtle
/proj/sllstore2017027/nobackup_GA2018/team_Wolf
/proj/sllstore2017027/nobackup_GA2018/team_Rhino
/proj/sllstore2017027/nobackup_GA2018/team_Rooster
{% endhighlight %}

### Illumina data.

* Subsample the data to 50x, 150x, and 250x coverage.
* Produce summaries of the subsampled data.
* Assemble the data using Spades, Abyss, and MaSuRCA.
* Polish using Pilon or Racon.
* Evaluate the assemblies with Quast, Busco, KAT, Bandage and FRC

Running spades example:

{% highlight bash %}
spades.py -k 21,33,55 --careful --pe1-1 "$READ1" --pe1-2 "$READ2" -o "${PREFIX}-spades_assembly"
{% endhighlight %}

Running abyss example:

{% highlight bash %}
abyss-pe name=abyss_k35_cleaned k=35 in='SRR492065_cleaned_R1.fastq.gz SRR492065_cleaned_R2.fastq.gz'
{% endhighlight %}

Running MaSuRCA:

{% highlight bash %}
# MaSuRCA needs a config file - You can use nano instead of cat
cat <<-EOF > "${PREFIX}_masurca.cfg"
DATA
PE= pe 500 50 $READ1 $READ2
END

PARAMETERS
GRAPH_KMER_SIZE = auto
END
EOF
masurca "${PREFIX}_masurca.cfg"
bash assemble.sh
{% endhighlight %}

### PacBio data.

* Subsample the data to 10x, 30x, and 70x coverage.
* Produce summaries of the subsampled data.
* Assemble the data using Canu, Miniasm, and wtdbg2.
* Polish with Arrow or Racon.
* Evaluate the assemblies with Quast, Busco, and Bandage.

{% highlight bash %}
canu -p ecoli -d ecoli-pacbio useGrid=false genomeSize=4.6m -pacbio-raw pacbio.fastq
{% endhighlight %}

{% highlight bash %}
# Overlap for PacBio reads (or use "-x ava-ont" for nanopore read overlapping)
minimap2/minimap2 -x ava-pb -t8 pb-reads.fq pb-reads.fq | gzip -1 > reads.paf.gz
# Layout
miniasm/miniasm -f reads.fq reads.paf.gz > reads.gfa
{% endhightlight %}

{% highlight bash %}
# assemble long reads
./wtdbg2 -t 16 -i reads.fa.gz -fo prefix -L 5000
# derive consensus
./wtpoa-cns -t 16 -i prefix.ctg.lay -fo prefix.ctg.lay.fa
{% endhighlight %}

### Nanopore data.

* Subsample the data to 10x, 30x, and 70x coverage.
* Produce summaries of the subsampled data.
* Assemble the data using Canu, Miniasm, and wtdbg2.
* Polish with Medaka or Racon.
* Evaluate the assemblies with Quast, Busco, and Bandage.

{% highlight bash %}
canu -p ecoli -d ecoli-oxford useGrid=false genomeSize=4.8m -nanopore-raw oxford.fasta
{% endhighlight %}

### How to load the tools.

You have already used most of the tools needed for this task. Here is how to
load the tools you have not encountered already.

Spades:
{% highlight bash %}
module load bioinfo-tools spades/3.12.0
{% endhighlight %}

Abyss:
{% highlight bash %}
module load bioinfo-tools abyss/2.0.2
{% endhighlight %}

MaSuRCA:
{% highlight bash %}
module load bioinfo-tools MaSuRCA/3.2.3
{% endhighlight %}

Pilon:
{% highlight bash %}
module load bioinfo-tools Pilon/1.22
{% endhighlight %}

Racon:
{% highlight bash %}
conda activate GA2018
{% endhighlight %}

Minimap2:
{% highlight bash %}
conda activate GA2018
{% endhighlight %}

Miniasm:
{% highlight bash %}
conda activate GA2018
{% endhighlight %}

Wtdbg2:
{% highlight bash %}
export PATH="$PATH:/proj/sllstore2017027/workshop-GA2018/tools/wtdbg2"
{% endhighlight %}

Medaka:
{% highlight bash %}
source /proj/sllstore2017027/workshop-GA2018/tools/medaka/venv/bin/activate
# to unload the virtual environment use
deactivate
{% endhighlight %}
