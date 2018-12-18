# Getting started in Genome Assembly

## Exercises:

### Task 1: Test your connection to uppmax.

1. Follow the instructions to log into the [UPPMAX cluster Rackham](../uppmax_login.md).

2. The tools needed for this course will be made available using a package manager called
`conda`. Check you can use these tools by doing the following:

{% highlight bash %}
source /proj/sllstore2017027/workshop-GA2019/tools/anaconda/miniconda2/etc/profile.d/conda.sh
conda activate GA2018
kat --help
conda deactivate
{% endhighlight %}

3. Test if graphical applications are available to you by typing the following:

{% highlight bash %}
conda activate GA2018
fastqc
conda deactivate
{% endhighlight %}

4. Please make a working directory for yourself in the workshop exercise directory using the
following commands (Note: `$USER` is a preset environment variable with your username in it):

{% highlight bash %}
mkdir -p /proj/sllstore2017027/nobackup_GA2019/$USER
cd /proj/sllstore2017027/nobackup_GA2019/$USER
{% endhighlight %}

### Task 2: Describe your own project

Describe your own assembly project. If you do not have a project of your own, use the ladybird, 
an insect with an estimated genome size of 200 Mbp. For this example, you want to assemble as 
much as possible of the gene space.

* Describe the properties of your genome of interest.
  * Estimated genome size:
  * Estimated repeat content:
  * Ploidy, and estimated heterozygosity if relevent:
  * Estimated GC content:
* How difficult is the level of assembly you're trying to achieve?
  * Are you aiming to assemble a genomic consensus?
  * Do you want phased alleles?
* What might hinder DNA extraction for your organism?
  * Can you get enough DNA from a single individual?
  * What are common contaminants from your species?
* Which sequencing technologies are suitable for your assembly?
* What kind of computational resources might you need?

