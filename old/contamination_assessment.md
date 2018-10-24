---
layout: default
title:  'Exercise: Contamination Assessment'
---

## Exercise: Contamination Assessment

**Questions**

If you have questions about the lab after the course, you are welcome to contact me: martin.norling@nbis.se

## Running BLAST

Many of you may have used BLAST before - either of the web or on the command line - but it's quite useful so we'll go through it again!
BLAST can be quite slow though, especially when running towards one of the more general databases! Therefore we offer a choice of two databases: either use `/sw/data/uppnex/blast_databases/nt` to get general, "real", results - or use `/proj/g2016024/nobackup/illumina_assembly/blastdb/bacterial` to run against a small bacterial database. The expected times to query an assembly against the two databases are about 5 minutes for nt and about 2 seconds for the small bacterial database.

To run blast you need to load the `module load blast/2.4.0+` database, then run `blastn -h` to get help on how to run nucleotide blast. Try different values of the argument `-outfmt` to get different kinds of output. The options are:
 - 0 = pairwise,
 - 1 = query-anchored showing identities,
 - 2 = query-anchored no identities,
 - 3 = flat query-anchored, show identities,
 - 4 = flat query-anchored, no identities,
 - 5 = XML Blast output,
 - 6 = tabular,
 - 7 = tabular with comment lines,
 - 8 = Text ASN.1,
 - 9 = Binary ASN.1
 - 10 = Comma-separated values
 - 11 = BLAST archive format (ASN.1)

With more detail available [here](https://www.ncbi.nlm.nih.gov/books/NBK279675/).

## Running BlobTools

In the read preparation you attempted to verify that only reads from the intended organism went into the assembly. We continue this verification here, by classifying the contigs and plotting them on a GC vs. coverage blob plot.
To do the classification – we will need taxonomy information! This can be downloaded (using the command `wget`) from NCBI, at this location: ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz. This file needs to be unpacked, so use the command `$ tar xzf taxdump.tar.gz` to unpack the contents.

We also need a BAM file of mapped reads. You already made on of these in the [assessment tutorial](assembly_assessment) though, so you're set! If not, just go back there for instructions on how to make one.

The third requirement is to actually have a BLAST classification of the contigs. We will use `blastn`, towards either of the databases mentioned in the BLAST part. Even though we just ran BLAST, in this case we will need a quite specific command to get a file that’s useful for blobtools, so use this command:

```
blastn -task megablast \
    -query [your assembly file].fasta \
    -db path-to-db \
    -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
    -culling_limit 5 \
    -num_threads 8 \
    -evalue 1e-25 \
    -out [output_name].blast.out
```

Now we can finally run blobtools! The uppmax module for blobtools is called `blobtools/0.9.17`. Sadly, the uppmax module is suffering from a small, but devastating, permission problem right now though, so it won't work. We've put a second version that should work at `/proj/g2016024/nobackup/illumina_assembly/blobtools/blobtools` though! You'll still need to load the module to get the dependencies right, but then use the full path to this version of the tool instead of the module version.
Blobtools works by using a creating a blob database - then you can query the database to get information and plots. But first - run `blobtools –help` to get more information on how to run blobtools.

Continue by creating a database using `blobtools create`. Use `blobtools create –help` to get help on the syntax. Now supply the **assembly**, **BWA mapping**, and **blast result**, as well as the nodes.dmp and names.dmp from the downloaded taxdump. Name the output based on the assembly you used.

You can now run `blobtools view` to create an output table, and `blobtools blobplot` to create a blobplot. 



    

