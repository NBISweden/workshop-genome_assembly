# Getting started in Genome Assembly

## Exercises:

### Test: Test your connection to uppmax.

Follow the instructions to log into the [UPPMAX cluster Rackham](../uppmax_login.md).

Test if graphical applications are available to you by typing the following:

{% highlight bash %}
module load bioinfo-tools FastQC
fastqc
module unload FastQC
{% endhighlight %}

Some tools are not available using the `module` system. We've provided these through a package manager called
`conda`. Please check you can use these environments.

{% highlight bash %}
source /proj/sllstore2017027/workshop-GA2018/tools/anaconda/miniconda2/etc/profile.d/conda.sh
conda activate GA2018
kat --help
conda deactivate
{% endhighlight %}

To give you enough disk space to work, please use the following commands to make a directory for yourself:

{% highlight bash %}
WORKDIR=/proj/sllstore2017027/nobackup_GA2018/$USER
mkdir -p "$WORKDIR"
cd "$WORKDIR"
{% endhighlight %}

### Task 1: Describe your own project

* Describe the properties of your genome of interest. If you can not think of a good example, use the insect ladybird, with an estimated genome size of 200 Mbp. For this example, you want to assemble as much as possible of the gene space.
* What might hinder DNA extraction for your organism?
* Which sequencing technologies are suitable for your assembly?
* What kind of computational resources might you need?
* How difficult is the level of assembly you're trying to achieve?
