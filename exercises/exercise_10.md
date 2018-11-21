---
title: De novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 24th November 2018
---
# Challenges in Assembly

## Exercise

### Task 1

Scaffolding genomes is usually only performed on large or complex genomes. As a result, this makes
it difficult to provide a dataset that runs in a short time on the given resources. For this task,
we provide you with links to familiarise yourself with so you can try these tools on your own resources
at a later date (this may be a good test to determine if your own resources are sufficient, develop expectations,
and perhaps even write a workflow).

Firstly, datasets can be hard to find. Here are some datasets that are publicly available.

[NBIS Genome Assembly Workshop Wiki - Datasets](https://github.com/NBISweden/workshop-genome_assembly/wiki/Datasets)

Next, tools to process that data are required. These should give you a good starting base:

[Long read scaffolding with Links](https://github.com/bcgsc/LINKS)

[10X Genomics scaffolding with Arcs](https://github.com/bcgsc/arcs)

[Hi-C Scaffolding with SALSA](https://github.com/machinegun/SALSA)

### Task 2.

10X Genomics is a fairly recent platform designed for assembling large genomes, but also with phased output.

Use the `Tiny` dataset and run the `supernova` pipeline.

{% highlight bash %}
module load bioinfo-tools bcl2fastq/2.20.0
PATH="$PATH:/proj/sllstore2017027/workshop-GA2018/tools/supernova-2.1.1"
{% endhighlight %}

Make a copy of the Tiny dataset files in your own working directory:

{% highlight bash %}
/proj/sllstore2017027/workshop-GA2018/data/10x_tiny
{% endhighlight %}

Extract the tar archive using the `tar` utility.

{% highlight bash %}
tar xvf tiny-bcl-2.0.0.tar.gz
{% endhighlight %}

Supernova needs to covert the raw Illumina data to fastq. Use the `supernova mkfastq` module to do this.

{% highlight bash %}
supernova mkfastq --run tiny-bcl-2.0.0 --id=tiny-bcl --samplesheet=tiny-bcl-samplesheet-2.1.0.csv
{% endhighlight %}

Now you can run Supernova on the fastq you just called.
{% highlight bash %}
supernova run --id=tiny --maxreads=all --fastqs tiny-bcl/outs/fastq_path/
{% endhighlight %}

Supernova offers various styles of output of the assembly. Take a look at the
[Supernova Support Page](https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/output/generating)
to understand the different kinds of output, and generate one of them.

### Task 3.

Write a slurm script for calculating the data quantity in a fastq file. Incorporate the `SLURM_ARRAY_TASK_ID` variable.

As we are working with limited resources, run this script in the terminal on your reserved node, instead of submitting to the slurm queue.
Normally Slurm will set `$SLURM_ARRAY_TASK_ID` for you when you use the `-a` option, but since we are running in the
terminal, we should set it ourselves and use the `export` command to make sure the script see's it.

{% highlight bash %}
export SLURM_ARRAY_TASK_ID=0
# Make your executable.
chmod 755 script.sh
# run the script by using the folder location ( ./ ) and the name of the script
./script.sh
{% endhighlight %}

Further reading:

[Uppmax SLURM guide](https://www.uppmax.uu.se/support/user-guides/slurm-user-guide/)

[Uppmax Disk Usage guide](https://www.uppmax.uu.se/support/user-guides/disk-storage-guide/)

[Uppmax Rackham User guide](https://www.uppmax.uu.se/support/user-guides/rackham-user-guide/)
