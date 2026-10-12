# nucleosome_methylation

Nucleosome positioning and DNA methylation analysis establishing that methylation directs nucleosome positioning in RNA-mediated transcriptional silencing in *Arabidopsis thaliana*.

Published as Rothi MH\*, **Sethuraman S**\*, Dolata J, Boyle AP, Wierzbicki AT. *DNA methylation directs nucleosome positioning in RNA-mediated transcriptional silencing.* bioRxiv. [doi:10.1101/2020.10.29.359794](https://doi.org/10.1101/2020.10.29.359794) (\*co-first authors)

## The question

DNA methylation and nucleosome positioning are tightly correlated at silenced loci, and either could be upstream of the other. Methylation might determine where nucleosomes sit, or nucleosome occupancy might determine where the methylation machinery can reach. Correlation cannot distinguish the two, so the question has to be approached by removing one and asking what happens to the other.

## Approach

The design works in both directions, using mutants that disrupt one process at a time.

**Removing methylation.** Differentially methylated regions were identified in the methyltransferase mutants *drm2-2*, *met1-3* and *cmt3-11*. Where methylation was lost, nucleosome signal dropped with it, which places methylation upstream of positioning.

**Removing positioning.** In *swi3b*/+, a heterozygote for a subunit of the SWI/SNF chromatin remodelling complex (null plants are embryo lethal), nucleosome signal fell at Pol V-stabilised nucleosomes while CHH methylation on SWI3B-stabilised nucleosomes was unchanged. Remodeller-driven positioning therefore does not determine methylation.

Together these give the direction rather than the association.

**A limitation the paper states explicitly**, and worth repeating here: these mutants lose methylation globally, so the effects observed may act through both cis- and trans-acting factors. Fully separating them would need tools that remove methylation at specific loci rather than genome-wide.

## Methods

**Nucleosome positioning** — paired-end MNase H3 ChIP-seq, two biological replicates, aligned to TAIR10 with bowtie2. Duplicates removed with PICARD. Reads retained at fragment lengths of 120–170 bp and MAPQ ≥ 2. Differential nucleosomes called with DANPOS2 at greater than two-fold change and FDR < 0.05, then tested across replicates with a negative binomial test (NBPseq, same thresholds). Well-positioned nucleosomes selected as main-peak nucleosomes using iNPS. Signal RPM-normalised and counted at dyads with BEDTools v2.15.0. SWI3B-stabilised nucleosomes required higher counts in Col-0 than in *swi3b* at FDR < 0.05, followed by an NBPseq test at FDR < 0.01.

**DNA methylation** — WGBS reads mapped to TAIR10 with Bismark allowing no mismatches. Methylation level calculated as C/(C+T) for cytosines covered by at least five reads. DMRs identified with methylKit using 100 bp bins at a 50 bp step, a minimum of 10 bases per tile, a methylation difference of at least 10%, and FDR 0.01. MNase H3 ChIP-seq signal across these DMRs was plotted as profiles.

**Genotypes** — Col-0 (wild type), *nrpe1* (*nrpd1b-11*), *ago4-1*, *idn2-1*, *swi3b*/+, *drm2-2*, *met1-3*, *cmt3-11*.

## Repository layout

| Path | Contents |
|---|---|
| `nucleosome_calling/nucleosome_calling/` | Nucleosome identification, differential calling and positioning analysis from MNase H3 ChIP-seq |
| `nucleosome_calling/DMR_calling/` | Differential methylation analysis from WGBS, and the comparison of MNase signal across DMRs in the methyltransferase mutants |

Scripts expect to be run from their own directory, with shared helper scripts at the repository root and intermediate files one level above the analysis folders. The relative paths inside the scripts follow that arrangement.

## Scripts not included

Two Perl helper scripts are called by the analysis scripts but are not part of this repository:

| Script | Called as | What it does |
|---|---|---|
| `convert_2d_2_MNase_2.pl` | `perl ../../convert_2d_2_MNase_2.pl <counts>.bed 6 0 > <counts>.csv` | Reshapes per-region MNase count BED files into the CSV matrices used for plotting |
| `convert_2d.pl` | `perl convert_2d.pl` | Equivalent reshaping step in the nucleosome analysis |

Both were small format-conversion utilities that sat outside the analysis directories at the time and were not preserved with the rest of the code. The descriptions above are inferred from how they are called, since the scripts themselves are no longer available.

This repository is published as a record of the analysis methods behind the paper rather than as a turnkey pipeline. The statistical and analytical steps are fully specified in the scripts and in the methods section above; reproducing the figures end to end would require reimplementing these two conversion utilities.

## Data

Sequencing data generated for this study are in GEO under accession [GSE148173](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148173).

Previously published datasets reanalysed here: [GSE70290](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70290), [GSE39901](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39901), and SRA054962.

## Software

bowtie2 · PICARD · DANPOS2 · iNPS · NBPseq (R) · Bismark · methylKit (R) · BEDtools v2.15.0

Reference genome: TAIR10.
