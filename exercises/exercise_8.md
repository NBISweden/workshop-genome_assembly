---
title: De Novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 23rd November 2018
---
# Assembly Validation Exercises.

### Task 1.

Run BUSCO on your assemblies.

{% highlight bash %}
mkdir BUSCO
cd BUSCO
ln -s ../*.fasta .
# run the following two lines to setup the busco environment
rsync -r /opt/miniconda3/pkgs/augustus-3.2.2-boost1.61_3/config/{species,model} .
export AUGUSTUS_CONFIG_PATH=$PWD

apply_BUSCO () {
	ASSEMBLY="$1" #The assembly is the first parameter to this function. The file must end in fasta
	LINEAGE=/home/data/opt-byod/busco/lineages/bacteria_odb9
	busco -i "$ASSEMBLY" -l "$LINEAGE" -c "${CPUS:-4}" -m genome -o "${ASSEMBLY/.fasta/_busco}"
}
{% endhighlight %}

What does the gene space look like for each assembly?

### Task 2.

Use Mauve to compare the assemblies to the reference.

The reference is here:
{% highlight bash %}
{% endhighlight %}

Reorder the assemblies with respect to the reference (`Tools > Move contigs`), and then make an alignment (`align with ProgressiveMauve`).

### Task 3.

Let's recircularize a circular assembly.

Use Mummer to find the coordinates of the overlap.

{% highlight bash %}
{% endhighlight %}

### Task 4.

Use Prokka to annotate the assembly and find a point to break the assembly near the origin of replication.

{% highlight bash %}
{% endhighlight %}

### Task 5.

Break the assembly overlap the ends using AMOS.

{% highlight bash %}
{% endhighlight %}

### Task 6.

Polish the assembly again to reduce errors in the overlap region.

{% highlight bash %}
{% endhighlight %}
