---
title: De novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 24th November 2018
---
# Ecoli genome assembly project

## Exercise

This is a team project, so split up the workload as you see fit.

You have three datasets from the Ecoli K12 substrain MG1655, sequenced using Illumina, PacBio, and Nanopore.

The aim is for you to try and explore different assemblers and see what you get. It will be impossible to
do evaluate every combination, so choose your tasks wisely.

### Illumina data.

a. Subsample the data to 50x, 150x, and 250x coverage.
b. Produce summaries of the subsampled data.
c. Assemble the data using Spades, Abyss, and MaSuRCA.
d. Polish using Pilon or Racon.
e. Evaluate the assemblies with Quast, Busco, KAT, Bandage and FRC

### PacBio data.

a. Subsample the data to 10x, 30x, and 70x coverage.
b. Produce summaries of the subsampled data.
c. Assemble the data using Canu, Miniasm, and wtdbg2.
d. Polish with Arrow or Racon.
e. Evaluate the assemblies with Quast, Busco, and Bandage.

### Nanopore data.

a. Subsample the data to 10x, 30x, and 70x coverage.
b. Produce summaries of the subsampled data.
c. Assemble the data using Canu, Miniasm, and wtdbg2.
d. Polish with Medaka or Racon.
e. Evaluate the assemblies with Quast, Busco, and Bandage.

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
{% endhighlight}

Miniasm:
{% highlight bash %}
conda activate GA2018
{% endhighlight}

Wtdbg2
{% highlight bash %}
export PATH="$PATH:/proj/sllstore2017027/workshop-GA2018/tools/wtdbg2"
{% endhighlight}

Medaka:
{% highlight bash %}
source /proj/sllstore2017027/workshop-GA2018/tools/medaka/venv/bin/activate
# to unload the virtual environment use
deactivate
{% endhighlight}
