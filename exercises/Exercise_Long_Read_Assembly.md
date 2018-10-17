---
layout: default
title:  'Exercise: Long Read Assembly'
---

## Exercise: Long Read Assembly

For this exercise we will focus on assembling an E. coli dataset sequenced using the Pacific Biosciences platform.

### Setup

Copy the data to your working directory.

```bash
cp -vr /sw/courses/assembly/PacBio_Assembly/Ecoli_Data .
```

### Miniasm

```bash
module load bioinfo-tools
module load miniasm minimap
export PATH="/proj/g2017025/tools/racon/bin:$PATH"
```

Use miniasm to assemble the **p6_25x_ecoli.fastq.gz**.

```bash
minimap -Sw5 -L100 -m0 -t4 p6_25x_ecoli.fastq.gz p6_25x_ecoli.fastq.gz | gzip -1 > minimap_ecoli_overlaps.paf.gz
miniasm -f p6_25x_ecoli.fastq.gz minimap_ecoli_overlaps.paf.gz > miniasm_ecoli_contigs.gfa
awk '/^S/ {print ">"seq++"\n"$3}' miniasm_ecoli_contigs.gfa > contigs.fasta
minimap -t 4 contigs.fasta p6_25x_ecoli.fastq.gz | racon -t 4 p6_25x_ecoli.fastq.gz - contigs.fasta consensus_contigs.fasta
```

### Canu

```bash
module load canu/1.5
```

Use Canu to assemble the **p6_25x_ecoli.fastq.gz**. Assembly takes approximately 20 mins, so please start writing the config file
for the next exercise while you wait.

Canu can auto-detect cluster settings but this is not recommended for running on milou. Please use the options
`useGrid=false` and `maxThreads=8` to limit the resources used when running on a node.

```bash
canu -p ecoli_canu -d ecoli_canu_assembly genomeSize=46m useGrid=false maxThreads=8 -pacbio-raw p6_25x_ecoli.fastq.gz
```

### Falcon

```bash
module load FALCON-integrate/20161113
module load seqtk
```

Create a Falcon config file and then run Falcon.

```bash
seqtk seq -l 5000 -A p6_25x_ecoli.fastq.gz > p6_25x_ecoli.fasta
echo $PWD/p6_25x_ecoli.fasta > input.fofn
fc_run.py falcon.cfg
```

This is what an example Falcon config might look like, including an explanation of the parameters. 
```
[General]
job_type = local
input_fofn = input.fofn
input_type = raw

length_cutoff = 5000

length_cutoff_pr = 5000

# grid settings for...
# daligner step of raw reads
sge_option_da = -pe smp 8
# las-merging of raw reads
sge_option_la = -pe smp 8
# consensus calling for preads
sge_option_cns = -pe smp 8
# daligner on preads
sge_option_pda = -pe smp 8
# las-merging on preads
sge_option_pla = -pe smp 8
# final overlap/assembly 
sge_option_fc = -pe smp 8

# concurrency settings
pa_concurrent_jobs = 8
ovlp_concurrent_jobs = 8
cns_concurrent_jobs = 8

# overlapping options for Daligner 
pa_HPCdaligner_option = -B4 -t16 -e.70 -l1000 -s1000 
ovlp_HPCdaligner_option = -B4 -t32 -h60 -e.96 -l500 -s1000

# How the database is split up for making comparison blocks
pa_DBsplit_option = -x1000 -s50 -a
ovlp_DBsplit_option = -x1000 -s50 -a

# error correction consensus option 
falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --max_n_read 200 --n_core 6

# overlap filtering options 
overlap_filtering_setting = --max_diff 100 --max_cov 100 --min_cov 20 --bestn 10

```

A more up to date explanation can be found in the [Falcon Documentation](http://pb-falcon.readthedocs.io/en/latest/), but I try to explain the parameters below.
```
[General]
jobtype = local             # other values sge, slurm
input_fofn = input.fofn
input_type = raw            # uncorrected reads
#input_type = preads        # falcon corrected reads

# The length cutoff used for seed reads used in initial mapping - these make the corrected reads
length_cutoff = 12000       # use longest 30X coverage

# The length cutoff used for seed reads used for pre-assembly - the min length of corrected reads
length_cutoff_pr = 12000    # 0-5000 lower than above

# concurrency settings
pa_concurrent_jobs = 8 # pre-assembly
ovlp_concurrent_jobs = 8 # overlap
cns_concurrent_jobs = 8 # consensus

# overlapping options for Daligner 
pa_HPCdaligner_option = -dal4 -t16 -e.70 -l1000 -s1000 
ovlp_HPCdaligner_option = -dal4 -t32 -h60 -e.96 -l500 -s1000

# -B <int>, -dal <int> 
# blocks to compare => higher = less but longer jobs

# -e <int>   
# average correlation rate (def 70%)

# -v # turns on verbose
# -l <int>
# the length in base pairs of the minimum local alignment (def. 1000)

# -s <int>
# how frequently trace alignments measured in bases are recorded (def. 100)

# -b 
# daligner assumes the data has a strong compositional bias (e.g. >65% AT rich).

# -t <int>,-M <int>       # Limits the effects of repeats
# Invariably, some k-mers are significantly over-represented (e.g. homopolymer runs). These k-mers create an excessive number of matching k-mer pairs and left unaddressed would cause daligner to overflow the available physical memory.  One way to deal with this is to explicitly set the -t parameter which suppresses the use of any k-mer that occurs more than t times in either the subject or target block.  However, a better way to handle the situation is to let the program automatically select a value of t that meets a given memory usage limit specified (in Gb) by the -M parameter.  By default daligner will use the amount of physical memory as the choice for -M.  If you want to use less, say only 8Gb on a 24Gb HPC cluster node because you want to run 3 daligner jobs on the node, then specify -M8.  Specifying -M0 basically indicates that you do not want daligner to self adjust k-mer suppression to fit within a given amount of memory.

# -H <int> 
# By default daligner compares all overlaps between reads in the database that are greater than the minimum cutoff set when the DB or DBs were split, typically 1 or 2 Kbp.  However, the HGAP assembly pipeline only wants to correct large reads, say 8Kbp or over, and so needs only the overlaps where the a-read is one of the large reads.  By setting the -H parameter to say N, one alters daligner so that it only reports overlaps where the a-read is over N base-pairs long.
# Essentially limits making alignments of reads of any size only to reads longer than <int>

# -k <int>, -h <int>, -w <int>
# The options -k, -h, and -w control the initial filtration search for possible matches between reads.  Specifically, the daligner search code looks for a pair of diagonal bands of width 2^w (default 2^6 = 64) that contain a collection of exact matching k-mers (default 14) between the two reads, such that the total number of bases covered by the k-mer hits is h (default 35). k cannot be larger than 32 in the current implementation.

# How the database is split up for making comparison blocks
pa_DBspliAt_option = -x1000 -s50 -a
ovlp_DBsplit_option = -x1000 -s50 -a

# -x <int>
# Ignore reads lower than length

# -s <int>
# specifies number of mb in each DB chunk - larger numbers makes smaller numbers of longer jobs (should be 400 mb or so for large genomes)

# -a # ignore secondary reads from the same well

# error correction consensus option 
falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --local_match_count_threshold 2 --max_n_read 200 --n_core 6

# --min_cov <int>
# break/trim seed read lower than <int>

# --max_n_read <int>
# max reads used for error correction - reduce value for highly repetitive genomes

# overlap filtering options 
overlap_filtering_setting = --max_diff 100 --max_cov 100 --min_cov 20 --bestn 10

# --bestn <int>
# Use the <int> best overlaps to simplify transitive edges in the graph

# --max_cov <int>, --min_cov <int>
# filter overlaps that are too high or too low (e.g. reads ending in repeats, or many sequencing errors)

# --max_diff <int>
# Max difference of coverage between 5’ and 3’ ends
```

### Quast

Use quast to compare the assemblies with the E. coli reference.

### Gepard

Compare the draft assemblies and the reference to each other using Gepard.

```bash
GEPARD_HOME=/sw/courses/assembly/Tools
java -cp $GEPARD_HOME/Gepard-1.40.jar org.gepard.client.cmdline.CommandLine -seq1 <fasta1> -seq2 <fasta2> -outfile <comparison.png> -matrix $GEPARD_HOME/matrices/edna.mat
```

### Polishing assemblies using Quiver.

```bash
module load SMRT/2.3.0
```

Polishing is a necessary step to remove indel and SNP errors. Use Quiver from the SMRT portal to polish a draft assembly. 
Note that polishing requires the pulsefield information that is encoded in the raw \*.h5 files. This means that the fastq files are insufficient.

The full PacBio data for the Ecoli dataset can be found here:
```bash
/proj/g2017025/exercises/PacBio_Assembly/Ecoli_raw_data
```

Modify the environment variables to the appropriate settings for your data (ASSEMBLY).

```bash
source $SMRT_SETUP_SCRIPT

PACBIO_DATA_DIR=/proj/g2017025/exercises/PacBio_Assembly/Ecoli_raw_data
PROTOCOL_XML=quiver_and_bridge_mapping.xml

ASSEMBLY=draft_assembly.fasta

referenceUploader -c -p $PWD -n E_coli -f $ASSEMBLY --saw='sawriter -welter' --samIdx='/sw/apps/bioinfo/samtools/1.5/milou/bin/samtools faidx'
SMRT_REFERENCE=$PWD/Ecoli
DIPLOID=False

# The following protocol is the BridgeMapper Protocol, and includes Quiver polishing as part of it. This was protocol was taken from a settings file from a dummy Bridge Mapping Protocol job.
cat > $PROTOCOL_XML <<EOF
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<smrtpipeSettings>
    <protocol version="2.3.0" id="RS_BridgeMapper.1" editable="false">
        <param name="name" label="Protocol Name">
            <value>RS_BridgeMapper</value>
            <input type="text"/>
            <rule required="true"/>
        </param>
        <param name="description">
            <value>Troubleshoot de novo assemblies,variants, indels, and so on.</value>
            <textarea></textarea>
        </param>
        <param name="version" hidden="true">
            <value>1</value>
            <input type="text"/>
            <rule type="digits" required="true" min="1.0"/>
        </param>
        <param name="state" hidden="true">
            <value>active</value>
            <input value="active" type="radio"/>
            <input value="inactive" type="radio"/>
        </param>
        <param name="reference" hidden="true">
            <value>$SMRT_REFERENCE</value>
        </param>
        <param name="control" hidden="true"/>
        <param name="fetch" hidden="true">
            <value>common/protocols/preprocessing/Fetch.1.xml</value>
        </param>
        <param name="filtering">
            <value>common/protocols/filtering/SFilter.1.xml</value>
            <select multiple="true">
                <import extension="xml" contentType="text/directory">common/protocols/filtering</import>
            </select>
		</param>
        <param name="spikeinControl" hidden="true">
            <value>common/protocols/control/SControl.1.xml</value>
        </param>
        <param name="mapping">
            <value>common/protocols/mapping/BLASR_Resequencing.1.xml</value>
            <select multiple="true">
                <import extension="xml" contentType="text/directory">common/protocols/mapping</import>
            </select>
        </param>
        <param name="bridgemapper">
            <value>common/protocols/postprocessing/BridgeMapper.1.xml</value>
            <select multiple="true">
                <import extension="xml" contentType="text/directory">common/protocols/postprocessing</import>
            </select>
        </param>
        <param name="consensus">
            <value>common/protocols/consensus/Quiver.1.xml</value>
            <select multiple="true">
                <import extension="xml" contentType="text/directory">common/protocols/consensus</import>
            </select>
        </param>
    </protocol>
    <moduleStage name="fetch" editable="true">
        <module label="Fetch v1" id="P_Fetch" editableInJob="true">
            <description>Sets up inputs</description>
        </module>
    </moduleStage>
    <moduleStage name="filtering" editable="true">
        <module label="SFilter v1" id="P_Filter" editableInJob="true">
            <description>This module filters reads based on a minimum subread length, polymerase read quality and polymerase read length.</description>
            <param name="minSubReadLength" label="Minimum Subread Length">
                <value>50</value>
                <title>Subreads shorter than this value (in base pairs) are filtered out and excluded from analysis.</title>
                <input type="text" size="3"/>
                <rule type="number" min="0.0" message="Value must be a positive integer"/>
            </param>
			<param name="readScore" label="Minimum Polymerase Read Quality">
                <value>75</value>
                <title>Polymerase reads with lower quality than this value are filtered out and excluded from analysis.</title>
                <input type="text"/>
                <rule type="number" min="0.0" message="Value must be between 0 and 100" max="100.0"/>
            </param>
            <param name="minLength" label="Minimum Polymerase Read Length">
                <value>50</value>
                <title>Polymerase reads shorter than this value (in base pairs) are filtered out and excluded from analysis.</title>
                <input type="text" size="3"/>
                <rule type="number" min="0.0" message="Value must be a positive integer"/>
            </param>
        </module>
        <module label="SFilter Reports v1" id="P_FilterReports" editableInJob="false"/>
    </moduleStage>
    <moduleStage name="spikeinControl" editable="true">
        <module label="SControl v1" id="P_Control" editableInJob="true">
            <param name="pbalign_opts" hidden="true">
                <value>--maxHits=1 --minAccuracy=0.75 --minLength=50 --algorithmOptions="-useQuality"</value>
            </param>
        </module>
        <module label="SControl Reports v1" id="P_ControlReports" editableInJob="false"/>
    </moduleStage>
    <moduleStage name="mapping" editable="true">
        <module label="BLASR v1" id="P_Mapping" editableInJob="true">
            <description>
BLASR maps reads to genomes by finding the highest scoring local alignment or set of local alignments between the read and the genome. The first set of alignments is found by querying an index of the reference genome, and then refining until only high scoring alignments are retained.  Additional pulse metrics are loaded into the resulting cmp.h5 file to enable downstream use of the Quiver algorithm.
            </description>
            <param name="maxHits" label="Maximum number of hits per read" hidden="true">
                <value>10</value>
                <title>
                The maximum number of matches of each read to the reference
                sequence that will be evaluated. maxHits should be greater
                than the expected number of repeats if you want to spread hits
                out on the genome.
                </title>
                <input type="text"/>
                <rule type="digits" message="Value must be an integer between 0 and 1000"/>
            </param>
			<param name="maxDivergence" label="Maximum divergence (%)">
                <value>30</value>
                <title>The maximum allowed divergence (in %) of a read from the reference sequence.</title>
                <input type="text"/>
                <rule type="digits" message="Value must be an integer between 0 and 100"/>
            </param>
            <param name="minAnchorSize" label="Minimum anchor size">
                <value>12</value>
                <title>The minimum size of the read (in base pairs) that must match against the reference.</title>
                <input type="text"/>
                <rule type="digits" message="Value must be an integer between 8 and 30"/>
            </param>
            <param name="samBam" label="Write output as a BAM file">
                <value>True</value>
                <title>Specify whether or not to output a BAM representation of the cmp.h5 file.</title>
                <input type="checkbox"/>
            </param>
            <param name="gff2Bed" label="Write BED coverage file">
                <value>True</value>
                <title>Specify whether or not to output a BED representation of the depth of coverage summary.</title>
                <input type="checkbox"/>
            </param>
            <param name="placeRepeatsRandomly" label="Place repeats randomly">
                <value>True</value>
                <title>Specify that if BLASR maps a read to more than one location with equal probability, then it randomly selects which location it chooses as the best location. If not set, defaults to the first on the list of matches.</title>
                <input type="checkbox"/>
            </param>
            <param name="pbalign_opts" hidden="true">
                <value>--seed=1 --minAccuracy=0.75 --minLength=50 --concordant --algorithmOptions="-useQuality"</value>
            </param>
            <param name="pbalign_advanced_opts" label="Advanced pbalign options">
                <value> </value>
                <title>Specify additional Pbalign options. For advanced users only.</title>
                <input type="text"/>
            </param>
            <param name="pulseMetrics" hidden="true">
                <value>DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag</value>
            </param>
            <param name="loadPulsesOpts" hidden="true">
                <value>bymetric</value>
                <title>The default option of loadPulses is 'byread'. Option 'bymetric'
               is desined to sacrifice memory for increased speed, especially
               for jobs of which the number of reference contigs is large. </title>
            </param>
        </module>
        <module label="BLASR Reports v1" id="P_MappingReports" editableInJob="false"/>
    </moduleStage>
	<moduleStage name="consensus" editable="true">
        <module label="Quiver v1" id="P_GenomicConsensus" editableInJob="true">
            <description>Quiver identifies haploid SNPs and indels by performing a local realignment of reads using the full range of sequence quality metrics.</description>
            <param name="algorithm" label="Consensus algorithm">
                <value>quiver</value>
                <title>Specify the consensus/variant algorithm to use in the analysis.</title>
                <input value="plurality" type="radio"/>
                <input value="quiver" type="radio"/>
            </param>
            <param name="outputConsensus" label="Output consensus FASTA and FASTQ files">
                <value>True</value>
                <title>Specify whether or not to output FASTA and FASTQ representations of the consensus sequence.</title>
                <input type="checkbox"/>
            </param>
            <param name="makeVcf" label="Write SNPs/Variants as VCF file">
                <value>True</value>
                <title>Specify whether or not to output a VCF representation of the variants.</title>
                <input type="checkbox"/>
            </param>
            <param name="makeBed" label="Write SNPs/Variants as BED file">
                <value>True</value>
                <title>Specify whether or not to output a BED representation of the variants.</title>
                <input type="checkbox"/>
            </param>
            <param name="enableMapQVFilter" label="Use only unambiguously mapped reads">
                <value>True</value>
                <title>Specify whether or not to filter out reads where Map QV is less than 10. Reduces coverage in repeat regions that are shorter than the read length.</title>
                <input type="checkbox"/>
            </param>
            <param name="minConfidence" label="Minimum variant confidence" hidden="true">
                <value>40</value>
                <title>Minimum variant confidence</title>
            </param>
            <param name="minCoverage" label="Minimum coverage requirement" hidden="true">
                <value>5</value>
                <title>Minimum variant coverage</title>
            </param>
            <param name="diploidMode" label="Diploid analysis">
                <value>$DIPLOID</value>
                <title>Specify whether or not Quiver operates in diploid mode and calls variants with the assumption that there are two copies of the genome in the sample; the mapping specificity should be higher.</title>
                <input type="checkbox"/>
            </param>
        </module>
        <module label="Genomic Consensus Reports v1" id="P_ConsensusReports" editableInJob="false"/>
    </moduleStage>
	<moduleStage name="bridgemapper" editable="true">
        <module label="BridgeMapper" id="P_BridgeMapper" editableInJob="true">
            <description>BridgeMapper returns split alignments of PacBio reads using blasr.</description>
            <param name="minAffixLength" label="Min Affix Length">
                <value>50</value>
                <title>Only split alignments with secondary alignments longer than this value (in base pairs) are reported.</title>
                <input type="text"/>
                <rule type="digits" min="0.0" message="Value must be between than 0 and 1000 inclusive" max="1000.0"/>
            </param>
            <param name="minRootLength" label="Minimum Root Length">
                <value>250</value>
                <title>The minimum length (in base pairs) of a subread's primary alignment for its affixes to be considered.</title>
                <input type="text"/>
                <rule type="digits" min="0.0" message="Value must be between 0 and 5000 inclusive" max="5000.0"/>
            </param>
        </module>
    </moduleStage>
    <fileName>RS_BridgeMapper.1.xml</fileName>
</smrtpipeSettings>
EOF

find $PACBIO_DATA_DIR -name "*.bax.h5" > input.fofn
fofnToSmrtpipeInput.py input.fofn > input.xml
smrtpipe.py --params=$PROTOCOL_XML xml:input.xml
rm input.fofn
```
