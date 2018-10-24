---
title: De Novo Assembly Workshop
author: Mahesh Binzer-Panchal
date: 23rd November 2018
---
# Assembly Validation Exercises.

### BUSCO

Assessing gene space is a core aspect of knowing whether or not you have a good assembly.

Benchmarking Universal Single-Copy Orthologs
(BUSCO) sets are collections of orthologous groups with near-universally-distributed single-copy genes in each species, selected from
OrthoDB root-level orthology delineations across arthropods, vertebrates, metazoans, fungi, and eukaryotes. BUSCO groups were selected
from each major radiation of the species phylogeny requiring genes to be present as single-copy orthologs in at least 90% of the species;
in others they may be lost or duplicated, and to ensure broad phyletic distribution they cannot all be missing from one sub-clade.
The species that define each major radiation were selected to include the majority of OrthoDB species, excluding only those with unusually
high numbers of missing or duplicated orthologs, while retaining representation from all major sub-clades. Their widespread presence means
that any BUSCO can therefore be expected to be found as a single-copy ortholog in any newly-sequenced genome from the appropriate
phylogenetic clade. [Busco Supplementary Online Material.]

7. Run Busco on each assembly.
	```bash
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
	```

	* What does the gene space look like for each assembly?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	for FASTA in *.fasta; do
		apply_BUSCO "$FASTA"
	done
	```

	```
	# CLC contigs
	C:98.0%[S:98.0%,D:0.0%],F:0.0%,M:2.0%,n:148

	145	Complete BUSCOs (C)
	145	Complete and single-copy BUSCOs (S)
	0	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	3	Missing BUSCOs (M)
	148	Total BUSCO groups searched

	# Mira large contigs
	C:97.3%[S:97.3%,D:0.0%],F:0.0%,M:2.7%,n:148

	144	Complete BUSCOs (C)
	144	Complete and single-copy BUSCOs (S)
	0	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	4	Missing BUSCOs (M)
	148	Total BUSCO groups searched

	# Mira trimmed contigs
	C:96.6%[S:96.6%,D:0.0%],F:0.7%,M:2.7%,n:148

	143	Complete BUSCOs (C)
	143	Complete and single-copy BUSCOs (S)
	0	Complete and duplicated BUSCOs (D)
	1	Fragmented BUSCOs (F)
	4	Missing BUSCOs (M)
	148	Total BUSCO groups searched

	# Spades contigs
	C:98.0%[S:98.0%,D:0.0%],F:0.0%,M:2.0%,n:148

	145	Complete BUSCOs (C)
	145	Complete and single-copy BUSCOs (S)
	0	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	3	Missing BUSCOs (M)
	148	Total BUSCO groups searched

	# Spades trimmed contigs
	C:98.0%[S:98.0%,D:0.0%],F:0.0%,M:2.0%,n:148

	145	Complete BUSCOs (C)
	145	Complete and single-copy BUSCOs (S)
	0	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	3	Missing BUSCOs (M)
	148	Total BUSCO groups searched

	# Velvet contigs
	C:96.6%[S:96.6%,D:0.0%],F:0.7%,M:2.7%,n:148

	143	Complete BUSCOs (C)
	143	Complete and single-copy BUSCOs (S)
	0	Complete and duplicated BUSCOs (D)
	1	Fragmented BUSCOs (F)
	4	Missing BUSCOs (M)
	148	Total BUSCO groups searched
	```

	These results show that the CLC and Spades contigs have found the most bacterial BUSCO's. These are therefore likely
	to be more complete assemblies than the others, although they are still fairly complete assemblies.

	</details>

### Comparative Alignment

Comparative alignment is a useful tool to see how assemblies compare to each other. This can be useful to compare assemblies to a reference, or
to see if assemblies have large structural differences.

8. Use Mauve to compare all the assemblies. First use Mauve contig mover (`Tools > Move contigs`) to reorder contigs with respect to your chosen first assembly. Then use `align with ProgressiveMauve` to align the reordered contigs.
	```bash
	Mauve
	```

	<details>
	<summary> Solution - click to expand </summary>
	```bash
	Mauve
	```

	First select `Tools` from the menu, and then `Move contigs`. Then click `add sequence` and add the CLC contigs, followed by one of
	the other assemblies. This will re-order contigs with respect to the CLC assembly. At the same time, Mauve will also show the resulting
	conservation between the two assemblies. The final alignment is also the folder name that the reordered contigs will be in.

	Use grep to make proper fasta file that can be used with other programs
	```bash
	grep -v "^#" alignment5/alignment5 > Mira_large_reordered_contigs.fasta
	```

	</details>

9. Use Gepard to see how assemblies compare to each other.
	```bash
	mkdir Gepard
	cd Gepard
	ln -s ../*.fasta .
	apply_Gepard () {
		ASSEMBLY1="$1" # An assembly is the first parameter to this function. The file must end in .fasta
		ASSEMBLY2="$2" # An assembly is the second parameter to this function. The file must end in .fasta
		SUBSTITUTION_MATRIX=/home/data/opt-byod/gepard/resources/matrices/edna.mat
		java -Djava.awt.headless=true -Xmx1024m -cp /home/data/opt-byod/gepard/dist/Gepard-1.40.jar org.gepard.client.cmdline.CommandLine \
			-seq1 "$ASSEMBLY1" -seq2 "$ASSEMBLY2" -matrix "$SUBSTITUTION_MATRIX" -outfile "${ASSEMBLY1/.fasta/}_vs_${ASSEMBLY2/.fasta/}_dotplot.png"
	}
	```
	This bash snippet will generate every pairwise comparison of the set.
	```bash
	set -- *.fasta
	for a; do
		shift;
		for b; do
			<function> "$a" "$b"
		done
	done
	```

	* Do you see any large scale differences between the assemblies?

	<details>
	<summary> Solution - click to expand </summary>

	```bash
	set -- *.fasta
	for a; do
		shift;
		for b; do
			apply_Gepard "$a" "$b"
		done
	done
	```

	If you compare assemblies that have not had their contigs reordered with respect to another, the Gepard plots do not show a linear correspondence between contigs.

	![CLC vs Unordered contigs](images/gepard/CLC_contigs_vs_Mira_large_contigs_dotplot.png)

	![CLC vs Ordered contigs](images/gepard/CLC_contigs_vs_Mira_large_reordered_contigs_dotplot.png)

	</details>
