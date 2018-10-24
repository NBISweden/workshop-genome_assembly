---
layout: default
title:  'Exercise: Assembly Validation'
---

## Exercise: Assembly Validation

**Questions**

If you have questions about the lab after the course, you are welcome to contact me: martin.norling@nbis.se

## Gene Space

Finding genes and other genomic regions is generally the domain of annotation, not assembly – we can cheat a bit though!
Making a rough estimate of the gene content compared to what is expected from a complete assembly can be a hint at the state of the assembly (at least of the assembled gene space).

### BUSCO

Busco, Benchmarking Universal Single-Copy Orthologs, is a program that attempts to estimate how well an assembly covers the gene space of a particular linage. Load the BUSCO module with:

```
module load bioinfo-tools BUSCO/1.22
```

There will be some information printed when you load the BUSCO module - this information is really important so it's a good idea to read it! Run `BUSCO –h` to get information on how to run. You will want to run Busco using `-l $BUSCO_LINEAGE_SETS/bacteria` to look for bacteria, and `-c 8` to use 8 threads. 

When the program finishes, look through the output directory for the files `full_table_*`, `missing_buscos_list_*`, and `short_summary_*`.

 - Does the assembly seem to include the full gene set based on these results?

## Whole Genome Alignment

When a reference is available, we can use whole genome alignment to compare the assemblies and look for differences and similarities.

### MUMmer

The program we will use to look do whole genome alignment is called `nucmer` and is from the MUMmer package. To use nucmer, load the `MUMmer/3.23` module. Running nucmer is fairly straight forward, type `nucmer -h` to get some help, and then run with your assembly, the provided reference from yesterday (from the illumina assembly part) and use the `-p` flag to set output name.

This will produce a mapping between the reference and the assembly, but it won't be sorted, so the output will be very noisy. You can use the `show-tiling` command to sort the delta file like this:
```
show-tiling [filename].delta > [filename].tiling
```

#### Dotplots

There is a built-in tool called `mummerplot` that can help us plot the output. It will default to trying to render directly to X11 though, so run with the `--png` flag to write an image file. As with `nucmer` you should also use the `-p` flag to set output name.
Run mummerplot on the delta-file that nucmer produced, as well as the tiling file. Download both files and have a look!

 - Does your assembly match well to the reference?
 - Did the plots look like you thought they would?

The mummerplot of a tiling dataset isn't the same kind of plot as the original dotplot. To remedy this, we've also supplied a simple python script to make dotplots with. The script is called `dotplot.py` and it can read .coord and .tiling files. It is available in `/proj/g2016024/nobackup/illumina_assembly/scripts/dotplot.py`, and you need to load the `biopython` module for it to run.

### Mauve

For the final part we will use a program called Mauve, that you'll have to run on your own machine! Start by downloading Mauve for your operating system from the [Mauve website](http://darlinglab.org/mauve/).

You will also need to download the assemblies that you wish to look at, as well as the reference to you local machine using `scp`.

Start Mauve, and select `File` -> `Align with ProgressiveMauve`, and select the sequences that you wish to align. You can give it multiple sequences, but I'd recommend doing pairwise alignments. Once the alignment is done you can run `Tools` -> `Move Contigs` to have mauve recursively move the contigs around to match better each other better.

