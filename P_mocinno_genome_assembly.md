# Genome Assembly and Reference-Guided Scaffolding of *Pharomachrus mocinno*

## 1. Primary Assembly Using Oxford Nanopore Technologies (ONT)

The primary genome assembly of *Pharomachrus mocinno* was generated using Oxford Nanopore long-read sequencing data. Reads were trimmed, merged, and purged to remove redundant haplotigs, resulting in a high-contiguity draft assembly composed exclusively of contigs (i.e., no scaffolding at this stage).

Assembly quality and contiguity statistics were assessed using **QUAST**, and the main metrics are summarized in Table 1.

### Table 1. QUAST statistics for the primary ONT assembly

| Metric | Value |
|------|------|
| Assembly name | Pharomachrus_mocinno_ONT_trim_2000bp_merged.purged_FILT |
| Number of contigs (≥ 0 bp) | 91 |
| Number of contigs (≥ 1 kb) | 91 |
| Number of contigs (≥ 5 kb) | 91 |
| Number of contigs (≥ 10 kb) | 91 |
| Number of contigs (≥ 25 kb) | 90 |
| Number of contigs (≥ 50 kb) | 90 |
| Total assembly length (bp) | 1,137,361,244 |
| Largest contig (bp) | 88,601,189 |
| GC content (%) | 42.65 |
| N50 (bp) | 31,929,196 |
| N90 (bp) | 9,028,569 |
| auN (bp) | 42,851,402.8 |
| L50 | 10 |
| L90 | 34 |
| Ns per 100 kbp | 0.00 |

The assembly shows very high contiguity for a primary long-read assembly, with only 91 contigs spanning approximately 1.14 Gb and no ambiguous bases (Ns), indicating a high-quality ONT-based reconstruction prior to scaffolding.

---

## 2. Reference Genome Used for Scaffolding

To improve chromosomal organization, the assembly was scaffolded using **RagTag** with the genome of *Trogon surrucura* as a reference.

**Reference genome details**
- Accession: GCA_020746105.1  
- Assembly name: bTroSur1.pri.cur  
- Assembly level: Chromosome  

| Feature | Value |
|------|------|
| Genome size | ~1.2 Gb |
| Total ungapped length | ~1.2 Gb |
| Number of chromosomes | 42 |
| Number of organelles | 1 |
| Number of scaffolds | 181 |
| Scaffold N50 | 82.6 Mb |
| Scaffold L50 | 5 |
| Number of contigs | 434 |
| Contig N50 | 16.2 Mb |
| Contig L50 | 23 |
| GC content (%) | 43 |
| Genome coverage | 78× |

---

## 3. Reference-Guided Scaffolding Strategy

Because scaffolding was performed using a reference genome from a different species (*Trogon surrucura*), specific precautions were taken to minimize biologically implausible joins and sex-chromosome–related artifacts. Two key decisions guided the scaffolding strategy: (i) removal of the W chromosome from the reference genome, and (ii) the use of a stringent confidence threshold during RagTag scaffolding.

---

### 3.1 Removal of the W chromosome from the reference genome

The *T. surrucura* reference genome includes a W chromosome. In birds, females are heterogametic (ZW) and males are homogametic (ZZ), and the presence, structure, and completeness of the W chromosome can vary substantially among species and assemblies. Moreover, the sex of the sequenced *Pharomachrus mocinno* individual was unknown (or potentially male), making the inclusion of a W chromosome in the reference problematic.

To avoid erroneous placement of autosomal or Z-linked contigs onto the W chromosome, and to prevent spurious scaffolding driven by repetitive or sex-specific sequences, the W chromosome was removed from the reference genome prior to scaffolding.

```bash
awk '
BEGIN{skip=0}
/^>/{ skip = ($0 ~ /chromosome W/ || $0 ~ /SUPER_W/ || $0 ~ / chrW\\b/ || $0 ~ /chrW/ || $0 ~ / W, /) }
skip==0 {print}
' GCA_020746105.1_bTroSur1.pri.cur_genomic.fna \
> GCA_020746105.1_bTroSur1.noW.fna
````
---

### 3.2 RagTag scaffolding using the no-W reference with high-confidence alignments (-i 0.95) 

After removing the W chromosome, scaffolding was performed using **RagTag** with a stringent confidence threshold (-i 0.95). This ensures RagTag only uses highly reliable alignments to guide joins and ordering.

This is especially important for cross-species scaffolding, where:  
- lineage-specific rearrangements (inversions/translocations) may exist,
- repetitive regions can generate spurious or ambiguous mappings,
- over-aggressive scaffolding can introduce incorrect joins that are hard to detect later.

By enforcing **-i 0.95**, scaffolding becomes more conservative and prioritizes conserved syntenic blocks, reducing the risk of mis-scaffolding while still improving chromosome-scale organization.

```bash
ragtag.py scaffold \
  -t 48 \
  -u \
  -i 0.95 \
  --aligner minimap2 \
  -o ragtag_pharomachrus_vs_trogon_noW_u_i0.95 \
  GCA_020746105.1_bTroSur1.noW.fna \
  Pharomachrus_mocinno_ONT_trim_2000bp_merged.purged_FILT.fa
```
---
## 4. Post-scaffolding Assembly Evaluation

After reference-guided scaffolding with RagTag using a stringent confidence threshold (`-i 0.95`), the resulting assembly was evaluated using RagTag summary statistics, QUAST, and BUSCO to assess placement success, contiguity gains, and gene-space completeness.

---

### 4.1 RagTag scaffolding statistics

RagTag provides a summary of how many sequences were successfully placed onto reference-guided scaffolds, how many remained unplaced, and how many gaps were introduced during scaffolding. The results obtained from `ragtag.scaffold.stats` are summarized below.

| Metric | Value |
|------|------|
| Placed sequences | 59 |
| Placed bases (bp) | 1,053,587,136 |
| Unplaced sequences | 32 |
| Unplaced bases (bp) | 83,774,108 |
| Gap bases (bp) | 2,800 |
| Gap sequences | 28 |

The majority of the genome (~92.6% of the total assembly length) was successfully placed onto reference-guided scaffolds. Unplaced sequences likely correspond to lineage-specific regions, unresolved repetitive elements, or contigs lacking sufficiently high-confidence alignments to the reference genome, consistent with the conservative scaffolding strategy.

---

### 4.2 Assembly contiguity assessment with QUAST

Assembly contiguity and structural statistics of the scaffolded genome were evaluated using **QUAST**. Key metrics are reported below.

| Metric | Value |
|------|------|
| Assembly name | ragtag.scaffold |
| Number of contigs (≥ 0 bp) | 63 |
| Number of contigs (≥ 1 kb) | 63 |
| Number of contigs (≥ 5 kb) | 63 |
| Number of contigs (≥ 10 kb) | 63 |
| Number of contigs (≥ 25 kb) | 62 |
| Number of contigs (≥ 50 kb) | 62 |
| Total assembly length (bp) | 1,137,364,044 |
| Largest contig (bp) | 217,862,911 |
| GC content (%) | 42.65 |
| N50 (bp) | 81,037,523 |
| N90 (bp) | 15,303,621 |
| auN (bp) | 93,079,038.5 |
| L50 | 5 |
| L90 | 19 |
| Ns per 100 kbp | 0.25 |

Compared to the primary contig-level assembly, reference-guided scaffolding resulted in a substantial increase in contiguity, with the N50 increasing from ~31.9 Mb to ~81.0 Mb and the L50 decreasing from 10 to 5. The small number of gap bases reflects scaffold joins introduced by RagTag and remains negligible at the genome scale.

---

### 4.3 BUSCO assessment of the scaffolded assembly

Genome completeness was assessed using **BUSCO v5** with the *aves_odb10* lineage dataset:

```bash
busco -i ragtag.scaffold.fasta \
  -l aves_odb10 \
  -m genome \
  -o busco_ragtag \
  -c 24
```

#### BUSCO results
        C:99.4%[S:98.9%,D:0.5%],F:0.0%,M:0.6%,n:8338,E:10.5%
        8287    Complete BUSCOs (C)     (of which 868 contain internal stop codons)
        8248    Complete and single-copy BUSCOs (S)
        39      Complete and duplicated BUSCOs (D)
        2       Fragmented BUSCOs (F)
        49      Missing BUSCOs (M)
        8338    Total BUSCO groups searched

---

## 5. Mapping-based validation of the scaffolded assembly

To validate the structural correctness of the RagTag-scaffolded assembly, raw Oxford Nanopore reads were mapped back to the scaffolded genome. This analysis provides an independent assessment of global assembly accuracy and local validation of scaffold joins introduced during reference-guided scaffolding.

---

### 5.1 Mapping of raw ONT reads to the scaffolded genome

Raw ONT reads were aligned to the RagTag-scaffolded assembly using **minimap2** with parameters optimized for long-read mapping, followed by sorting and indexing with **samtools**.

```bash
ASM=/mnt/data/OVHlord/annotation/scaffolding/ragtag_pharomachrus_vs_trogon_noW_u/ragtag.scaffold.fasta
READS=/mnt/data/OVHlord/assemblies/projects/Pharomachrus_mocinno/Pharomachrus_mocinno_ONT_rawreads_trim_80_percent.fq.gz
T=24

minimap2 -t $T -ax map-ont $ASM $READS | \
  samtools sort -@ 8 -o ont_vs_ragtag.bam

samtools index ont_vs_ragtag.bam
```

Mapping statistics were obtained using `samtools flagstat`:

```bash
samtools flagstat ont_vs_ragtag.bam
```

**Mapping summary:**

| Metric | Value |
|------|------|
| Total reads | 57,710,002 |
| Mapped reads | 57,332,000 (99.34%) |
| Primary reads | 38,765,918 |
| Primary mapped reads | 38,387,916 (99.02%) |
| Secondary alignments | 6,997,009 |
| Supplementary alignments | 11,947,075 |
| Duplicates | 0 |

The very high overall and primary mapping rates (>99%) indicate excellent global consistency between the raw ONT reads and the scaffolded genome, supporting the correctness of the assembly at the genome-wide scale.

---

### 5.2 Identification of RagTag scaffold joins

To specifically evaluate the validity of scaffold joins introduced by RagTag, gap coordinates were extracted from the RagTag AGP file (`ragtag.scaffold.agp`). These gap regions represent scaffold joins where contigs were linked based on reference guidance.

```bash
cd /mnt/data/OVHlord/annotation/scaffolding/ragtag_pharomachrus_vs_trogon_noW_u
ls *.agp

AGP=ragtag.scaffold.agp

awk '$5=="N" || $5=="U" {print $1"\t"$2"\t"$3"\tJOIN_"NR"\t0\t+"}' \
  $AGP > ragtag_joins.bed

wc -l ragtag_joins.bed
head ragtag_joins.bed
```

The resulting BED file contains the genomic coordinates of all scaffold joins introduced during RagTag scaffolding.

---

### 5.3 Read-depth validation at scaffold joins

To assess whether scaffold joins are supported by sequencing data, read depth and read counts were evaluated in windows flanking each join. For each join, coverage statistics were calculated for regions immediately upstream (left) and downstream (right) of the gap.

```text
chr         join_start   join_end   gap_len  meanDepth_L  meanDepth_R  reads_L  reads_R
CM036617.1  16614648     16614748   100      59.078       64.358       1511     988
CM036617.1  89059618     89059718   100      61.279       95.946       421      892
CM036617.1  130622133    130622233  100      73.947       176.560      3359     6349
CM036617.1  132309650    132309750  100      554.594      61.590       22083    824
CM036617.1  133071808    133071908  100      216.056      78.694       3535     1593
CM036618.1  30107252     30107352   100      33.319       50.924       428      569
CM036618.1  116169687    116169787  100      254.658      90.928       3869     989
CM036618.1  125198356    125198456  100      70.367       356.192      485      3823
CM036619.1  88601189     88601289   100      87.662       121.283      1889     3009
CM036620.1  510519       510619     100      69.330       282.751      2147     2434
CM036620.1  18368978     18369078   100      255.008      172.931      5443     3324
CM036620.1  36139586     36139686   100      34.600       64.114       974      981
CM036621.1  13059471     13059571   100      143.717      408.188      5420     9947
CM036621.1  39402018     39402118   100      168.609      1046.483     2459     17223
CM036622.1  363430       363530     100      73.011       169.610      600      7752
CM036623.1  5140422      5140522    100      549.028      519.832      10135    8052
CM036625.1  16957475     16957575   100      31.152       13.111       324      86
CM036625.1  17410534     17410634   100      385.698      1091.843     6679     12167
CM036630.1  1033044      1033144    100      988.029      345.193      9981     4667
CM036631.1  1730660      1730760    100      157.453      80.605       2837     2886
CM036641.1  2411182      2411282    100      651.205      368.263      9442     6412
CM036644.1  3756572      3756672    100      34.243       571.453      572      11142
CM036658.1  7820692      7820792    100      223.599      201.928      2273     3504
CM036658.1  14271196     14271296   100      1239.190     80.465       31754    1587
CM036658.1  46200492     46200592   100      130.361      779.620      1753     13158
CM036658.1  46429913     46430013   100      46.552       69.250       1299     1260
CM036658.1  49702159     49702259   100      110.020      124.826      2254     2343
CM036658.1  74925696     74925796   100      56.246       45.242       804      845
```

All scaffold joins show substantial read depth and read support on both sides of the gap. Although coverage asymmetry is observed for some joins, such variation is expected in long-read datasets and may reflect local repeat structure, copy number variation, or differences in mappability.

Crucially, no join exhibits a complete loss of coverage on either flank, which would be indicative of unsupported or erroneous joins. These results support the structural validity of the scaffold joins introduced by RagTag under a stringent confidence threshold.

---

## 7. Repeat annotation and genome masking

Repetitive elements in the *Pharomachrus mocinno* genome were annotated using **EDTA (Extensive de novo TE Annotator)**. EDTA integrates multiple tools to identify, classify, and annotate transposable elements (TEs) in a unified framework and is particularly well suited for long-read, chromosome-scale assemblies.

---

### 7.1 Preparation of the input genome

Repeat annotation was performed on a renamed and corrected version of the RagTag-scaffolded genome. By default, RagTag produces scaffold headers that reflect the reference-guided scaffolding paths and may incorporate reference-derived identifiers. While informative, these composite headers are not ideal for downstream annotation workflows.

Prior to running EDTA, scaffold headers were therefore **normalized** to produce a clean, assembly-native naming scheme. Scaffolded sequences were assigned sequential identifiers (`Pmoc_scafXXXX`), while unplaced contigs were labeled separately (`Pmoc_unplacedXXXXXX`). This resulted in a consistent, species-specific identifier system that clearly distinguishes placed scaffolds from unplaced sequences and is independent of the reference genome used during scaffolding.

Example headers from the corrected genome file (`pharomachrus.ragtag.renamed_corrected.fa`):

```text
>Pmoc_scaf0001
>Pmoc_scaf0002
>Pmoc_scaf0003
...
>Pmoc_scaf0031
>Pmoc_unplaced000001
>Pmoc_unplaced000002
>Pmoc_unplaced000003
...
>Pmoc_unplaced000015
```

This normalization step improves clarity, reproducibility, and downstream interpretability of the assembly, particularly for repeat annotation and gene prediction.

---

### 7.2 EDTA repeat annotation

Repeat annotation was performed using EDTA in sensitive mode with annotation enabled. The full script used to run EDTA is shown below.

```bash
#!/bin/bash
set -euo pipefail

# Increase LMDB map size for BLAST/RepeatMasker
export BLASTDB_LMDB_MAP_SIZE=10000000

# ===== Paths =====
GENOME="/mnt/data/OVHlord/annotation/scaffolding/syri/pharomachrus.ragtag.renamed_corrected.fa"
WORKDIR="/mnt/data/OVHlord/annotation/repetitive/Quetzal_EDTA"
LOG="Quetzal_EDTA.log"

# ===== Resources =====
THREADS=24

# ===== Run =====
mkdir -p "$WORKDIR"
cd "$WORKDIR"

echo "[INFO] Copying genome..."
cp -f "$GENOME" Quetzal.raw.fa

echo "[INFO] Cleaning FASTA headers..."
awk '
  BEGIN{OFS=""}
  /^>/{
    h=$0
    sub(/^>/,"",h)
    split(h,a," ")
    gsub(/[^A-Za-z0-9_.:-]/,"",a[1])
    print ">",a[1]
    next
  }
  {print}
' Quetzal.raw.fa > Quetzal.fa

rm -f Quetzal.raw.fa

echo "[INFO] Starting EDTA..."
EDTA.pl \
  --threads "$THREADS" \
  --genome Quetzal.fa \
  --sensitive 1 \
  --force 1 \
  --anno 1 \
  2>&1 | tee "$LOG"

echo "[DONE] EDTA finished. Log saved to $WORKDIR/$LOG"
```

The `--sensitive` option was used to maximize repeat detection, which is appropriate for a first-pass annotation of a non-model avian genome. Annotation mode (`--anno 1`) was enabled to produce RepeatMasker-compatible outputs and summary statistics.

---

### 7.3 Generation of the soft-masked genome

After EDTA completed, a soft-masked version of the genome was generated using the EDTA utility `make_masked.pl`. This step applies the TE annotations produced by EDTA to the genome sequence, converting repeat regions to lowercase while retaining the full sequence length and coordinates.

Soft masking was preferred over hard masking to preserve sequence information for downstream gene prediction and annotation pipelines.

```bash
perl /home/cauz/miniconda3/envs/edta/share/EDTA/bin/make_masked.pl \
  -genome Quetzal.fa \
  -rmout ./Quetzal.fa.mod.EDTA.anno/Quetzal.fa.mod.EDTA.TEanno.out \
  -hardmask 0 \
  -minlen 80 \
  -t 32
```

The resulting soft-masked genome file was `Pharomachrus_mocinno_softmasked.fa`.

---

### 7.4 Summary statistics of repetitive elements

Repeat annotation statistics for the soft-masked genome are summarized below.

```text
File name: Quetzal.fa.mod
Sequences:            63
Total length:         1,137,364,044 bp (1,137,361,244 bp excluding N/X runs)
GC content:           42.65 %

Bases masked:         107,527,468 bp (9.45 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
Retroelements       202,045     95,654,049 bp    8.41 %
  SINEs               2,330        281,587 bp    0.02 %
  LINEs             175,581     81,088,781 bp    7.13 %
    L2/CR1/Rex       174,138     80,082,162 bp    7.04 %
    RTE/Bov-B          1,443      1,006,619 bp    0.09 %
  LTR elements       24,134     14,283,681 bp    1.26 %
    Ty1/Copia             98         32,665 bp    0.00 %
    Gypsy/DIRS1        1,255        636,445 bp    0.06 %
    Retroviral        1,412        450,834 bp    0.04 %

DNA transposons        6,422        909,215 bp    0.08 %
  hobo-Activator        474         30,824 bp    0.00 %
  Tc1-IS630-Pogo      1,341        498,261 bp    0.04 %
  Tourist/Harbinger     123         10,285 bp    0.00 %

Rolling-circles            0              0 bp    0.00 %

Unclassified           31,277     10,939,012 bp    0.96 %

Total interspersed repeats: 107,502,276 bp    9.45 %

Small RNA                 345         48,585 bp    0.00 %
Satellites                  0              0 bp    0.00 %
Simple repeats              0              0 bp    0.00 %
Low complexity               0              0 bp    0.00 %
==================================================
```

Approximately 9.45% of the *P. mocinno* genome was annotated as interspersed repetitive elements. As expected for avian genomes, repeats are dominated by **LINE elements**, particularly the **L2/CR1/Rex** superfamily, which alone accounts for ~7.0% of the genome. LTR retrotransposons contribute a smaller fraction (~1.3%), while DNA transposons are rare.

The relatively low overall repeat content is consistent with patterns observed in other bird genomes and supports the high contiguity and mappability of the assembly. 

---

## 8. Gene Annotation 

Gene annotation of the resplendent quetzal (*Pharomachrus mocinno*) was generated using **BRAKER3** under two complementary evidence strategies:

1. **Protein-only annotation**, based exclusively on conserved vertebrate protein evidence.
2. **Protein + RNA-seq annotation**, integrating protein evidence with heterospecific RNA-seq data from *Trogon surrucura*.

---

### 8.1 Protein-only annotation

Protein-based gene prediction was performed using **Vertebrata OrthoDB v11** protein sequences (`Vertebrata.fa`) as external evidence.

OrthoDB provides a curated set of orthologous proteins sampled across vertebrates and is commonly used to guide ab initio gene prediction in newly assembled genomes. These proteins were aligned to the *P. mocinno* soft-masked genome using the BRAKER3 protein evidence pipeline, which leverages spliced protein alignments to infer exon–intron structures.

While the protein-only annotation successfully identified coding regions across the genome, downstream structural statistics revealed:

- A strongly inflated number of predicted genes.
- Short gene spans and coding sequences.
- Reduced exon counts per gene.

These patterns are consistent with **gene fragmentation**, a known limitation of protein-only annotation approaches in vertebrate genomes, particularly when protein evidence is phylogenetically broad. As a result, the protein-only annotation was retained for comparison purposes but not selected as the final annotation.

---

### 8.2 RNA-seq evidence for annotation refinement

To improve gene model completeness and exon–intron structure reconstruction, RNA-seq data from a closely related species were incorporated as transcriptional evidence.

#### RNA-seq dataset

RNA-seq data were obtained from the following SRA accession:

- **Species:** *Trogon surrucura*
- **Tissue:** Skin
- **Sex:** Unknown
- **Accession:** SRR10852885 (from SRX7523463)
- **Study:** Transcriptome sequence-based phylogeny of birds  
- **BioProject:** PRJNA599522  
- **Instrument:** Illumina HiSeq 2500  
- **Library layout:** Paired-end  
- **Read length:** 50 bp  
- **Sequencing depth:**  
  - 41.8 million spots  
  - 4.2 Gb of sequence data  

The dataset was generated by the Max Planck Institute for Ornithology–Seewiesen.

---

### 8.3 RNA-seq mapping strategy and validation

Because the RNA-seq data originate from a **different species**, a **precision-first and conservative mapping strategy** was adopted. The goal was to retain only high-confidence transcriptional signals suitable for **splice site validation and exon boundary support**, while minimizing spurious cross-species alignments.

#### Mapping approach

RNA-seq reads were aligned to the *P. mocinno* soft-masked genome using **STAR v2.7.11b** in a **two-pass alignment framework**:

**First pass:** 
- Moderate stringency parameters were used to detect candidate splice junctions.
- Detected junctions were collected to inform the second pass.

**Second pass:**
The second-pass alignment was performed with increased stringency to prioritize accuracy over sensitivity. Key parameters included:

- Restriction to **uniquely mapped reads only** (`outFilterMultimapNmax = 1`)
- Strict mismatch threshold (maximum mismatch rate ≤ 4%)
- Minimum alignment and match length fractions ≥ 66% of read length
- Elevated splice junction overhang requirements
- Exclusion of ambiguous and low-quality alignments

#### Mapping results:

- **Unique mapping rate (second pass):** 73.6%
- **Multi-mappers retained:** 0 (by design)
- **Unmapped reads:** 24.1%, primarily classified as *too short*
- **Reads unmapped due to mismatches:** 0
- **Per-base mismatch rate:** 2.07%

The observed mismatch rate is well within expectations for cross-species RNA-seq mapping between closely related avian taxa and indicates high sequence similarity between *T. surrucura* transcripts and the *P. mocinno* genome.

Splice junction detection was robust, with over **5.7 million splice events** identified. All detected junctions corresponded exclusively to **canonical splice motifs** (GT/AG, GC/AG, or AT/AC), and no non-canonical junctions were observed.

#### Justification for use in genome annotation:

Despite originating from a different species and a single tissue, the RNA-seq dataset provides reliable transcriptional evidence for *P. mocinno* genome annotation due to:

- High unique mapping rates under strict, uniqueness-enforcing conditions.
- Low mismatch rates and absence of mismatch-driven read exclusion.
- Exclusive detection of canonical splice junctions, supporting biologically plausible exon–intron structures.
- Conservative filtering that limits RNA-seq influence to **structural validation** rather than de novo gene discovery in poorly conserved regions.

Accordingly, RNA-seq evidence was used to refine gene models by improving exon connectivity, intron placement, and coding sequence completeness.

---

### 8.4 Final annotation strategy

The final gene annotation integrates:

- **Protein evidence** from Vertebrata OrthoDB v11, providing broad evolutionary constraints.
- **RNA-seq evidence** from *Trogon surrucura*, providing high-confidence splice and exon boundary support.

Comparative structural analyses demonstrated that the combined protein + RNA-seq annotation yielded:

- A biologically realistic gene count for an avian genome.
- Substantially increased exon counts per gene.
- Longer gene spans and coding sequences.
- Markedly improved BUSCO completeness relative to protein-only annotation.


**Limitations:** RNA-seq evidence was derived from a single tissue and a different species. As such, the annotation should be regarded as a **high-quality structural annotation**, suitable for genome description and comparative genomics, but not an exhaustive representation of tissue-specific or developmental gene expression.

---
### 8.5 Comparison of genome annotations

Structural statistics and BUSCO completeness metrics were used to assess annotation quality and guide selection of the final annotation.

#### Structural annotation statistics

| Metric | Protein-only | Protein + RNA |
|------|-------------|---------------|
| Predicted genes | 34,441 | **15,054** |
| Predicted transcripts | 36,180 | **20,797** |
| Total exons | 179,972 | 173,476 |
| Total CDS features | 179,972 | 173,476 |
| Transcripts with exons | 36,180 | 18,816 |
| Single-exon transcripts | 5,124 | 3,079 |
| Multi-exon transcripts | 31,056 | 15,737 |
| Single-exon transcripts (%) | 14.16 | 16.36 |
| Mean exons per transcript | 4.97 | **9.22** |
| Median exons per transcript | 3 | **7** |
| Mean exons per gene | 5.23 | **11.52** |
| Median exons per gene | 3 | **6** |
| Mean exon length (bp) | 192.84 | 170.22 |
| Median exon length (bp) | 129 | 126 |
| Mean gene span (bp) | 4,820 | **16,231** |
| Median gene span (bp) | 2,498 | **6,463** |
| Mean CDS length per gene (bp) | 1,008 | **1,962** |
| Median CDS length per gene (bp) | 675 | **1,179** |

**OBS**: The gene annotation for *Trogon surrucura* from the Vertebrate Genome Project (VGP) contains 14,763 predicted genes (Info: https://projects.ensembl.org/vgp/).

#### BUSCO completeness (aves_odb10)

BUSCO analysis was performed using the **aves_odb10** dataset (n = 8,338 conserved avian genes).

| BUSCO category | Protein-only | Protein + RNA |
|---------------|-------------|---------------|
| Complete BUSCOs (%) | 24.0 | **81.1** |
| Complete, single-copy (%) | 21.0 | **65.5** |
| Complete, duplicated (%) | 3.0 | 15.7 |
| Fragmented (%) | 5.1 | **2.1** |
| Missing (%) | 71.0 | **16.8** |
| Total BUSCOs | 8,338 | 8,338 |


Although the protein-only annotation identified a large number of coding regions, it resulted in a strongly inflated gene count, short gene spans, and reduced coding sequence lengths, consistent with extensive gene fragmentation. This is further reflected in low BUSCO completeness (24.0%) and a high proportion of missing BUSCO genes (71.0%).

In contrast, integration of RNA-seq evidence substantially improved gene model reconstruction. The protein + RNA annotation yielded a biologically realistic number of genes for an avian genome, with increased exon counts per gene, longer gene spans, and nearly doubled CDS lengths. BUSCO completeness increased to 81.1%, with low fragmentation (2.1%), indicating a structurally coherent annotation.

Based on these results, the protein + RNA annotation was selected as the final gene annotation for the *Pharomachrus mocinno* genome.

---
##  9. Functional annotation

Protein-coding genes were functionally annotated using a combination of domain-based and orthology-based approaches. Conserved protein domains and Gene Ontology (GO) terms were identified using InterProScan v5.77, while orthologous groups, functional descriptions, and pathway annotations were inferred using eggNOG-mapper v2.1.13 with DIAMOND searches restricted to Aves (taxon ID 8782).

### 9.1 Results
```bash
cat quetzal_annotation_summary_stats.tsv
metric  value
interpro_proteins_with_hit      14251
interpro_proteins_with_go       14251
eggnog_proteins_with_record     12480
eggnog_proteins_with_go 1330
```

InterProScan identified conserved domains in 14,251 out of 15,054 predicted proteins (94.7%), all of which were associated with at least one Gene Ontology term. Orthology-based annotation using eggNOG-mapper assigned functional descriptions to 12,480 proteins (82.9%), with 1,330 proteins linked to curated GO terms. These results indicate a high level of completeness and functional characterization of the predicted gene set.

### 9.2 Technical validation: genomic organization of multi-copy gene families

As an internal quality control of the genome annotation, we evaluated the genomic organization of large gene families. In vertebrate genomes, several well-characterized gene families (for example, keratins, protocadherins, immunoglobulins, histones, and zinc-finger proteins) are typically arranged as local tandem arrays. Recovery of such patterns provides independent support for the structural validity of predicted gene models.

Protein-coding genes predicted by BRAKER were assigned to orthologous groups using eggNOG-mapper (Aves taxonomic scope). Multi-copy ortholog groups were mapped back to genomic coordinates using the BRAKER GTF file. For each group, we quantified family size, scaffold distribution, and the largest local gene cluster.

#### Results
Many of the largest orthologous groups occur as local tandem clusters, often with all or most family members located on a single scaffold. These families predominantly correspond to gene categories known to form tandem arrays in avian genomes.

| rank | Aves OG ID | genes | scaffolds | max genes on one scaffold | scaffold of max cluster | clustered fraction | top description |
|------|------------|-------|-----------|---------------------------|------------------------|--------------------|-----------------|
| 1 | 4GWI7 | 23 | 1 | 23 | Pmoc_scaf0031 | 1.000 | Ring finger |
| 2 | 4GVZ7 | 8 | 2 | 6 | Pmoc_scaf0004 | 0.750 | Inositol 1,4,5-trisphosphate receptor-interacting protein-like |
| 3 | 4GVIK | 8 | 1 | 8 | Pmoc_scaf0028 | 1.000 | Keratin |
| 4 | 4GVG4 | 8 | 3 | 6 | Pmoc_scaf0016 | 0.750 | Proline-rich protein |
| 5 | 4GV08 | 7 | 1 | 7 | Pmoc_scaf0028 | 1.000 | Scale keratin-like |
| 6 | 4GTGM | 7 | 1 | 7 | Pmoc_scaf0019 | 1.000 | Myosin family protein |
| 7 | 4GP9D | 6 | 1 | 6 | Pmoc_scaf0014 | 1.000 | Protocadherin |
| 8 | 4GIZJ | 6 | 1 | 6 | Pmoc_scaf0014 | 1.000 | Protocadherin |
| 9 | 4GR9J | 6 | 1 | 6 | Pmoc_scaf0018 | 1.000 | Immunoglobulin V-type |
|10 | 4GVMJ | 5 | 1 | 5 | Pmoc_unplaced000009 | 1.000 | Immunoglobulin V-type |

The complete results, including gene identifiers and additional families, are available in: `functional_annotation/merged/aves_og_genomic_clustering.tsv`

The observed genomic organization is consistent with known duplication mechanisms and gene family architectures in birds, and no atypical expansions of transposable element-associated or uncharacterized gene families were detected. Together, these results support the structural validity of the gene annotation and indicate that the predicted gene models capture biologically plausible gene family organization.

##  10. BUSCO TEST - Construction of protein evidence datasets for BRAKER3 tests

Two protein evidence datasets were generated to test whether alternative protein inputs could improve gene prediction completeness when combined with the same RNA-seq evidence and genome assembly.

All steps below were performed identically for both datasets where applicable.

### 10.1 OrthoDB v12 – Aves Protein Dataset Construction

**Aves-specific protein FASTA** was generated from **OrthoDB v12**.

Metadata tables required to map proteins to taxonomic clades were downloaded from the official OrthoDB v12 distribution:

```bash
wget https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_levels.tab
wget https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_level2species.tab
wget https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_species.tab
```

These files provide: 
- odb12v2_levels.tab – taxonomic hierarchy and clade names
- odb12v2_level2species.tab – mapping of species to taxonomic levels
- odb12v2_species.tab – species metadata

The Aves clade was identified as:

```bash
grep -w -i "Aves" odb12v2_levels.tab
# TaxID: 8782
```

#### Download OrthoDB v12 protein sequences

The complete amino acid FASTA for OrthoDB v12 was downloaded and decompressed:

```bash
wget https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_aa_fasta.gz
gunzip odb12v2_aa_fasta.gz
```

#### Extract Aves-specific proteins

Protein sequences belonging to the clade Aves were extracted using the selectClade.py script from the orthodb-clades repository, which filters FASTA entries based on OrthoDB taxonomic metadata:

```bash
python3 selectClade.py \
  odb12v2_aa_fasta \
  odb12v2_levels.tab \
  odb12v2_level2species.tab \
  Aves \
  > OrthoDB12_Aves.fa
```

Final dataset summary was observed

```bash
grep -c '^>' OrthoDB12_Aves.fa
```
#### 6,079,582 protein sequences

---

### 10.2 Input protein sources

The following protein FASTA files were used as inputs:

- **OrthoDB v12 (Aves clade)**: `OrthoDB12_Aves.fa`
- **Avian reference proteomes** (FASTA):
  - *Amazona aestiva*
  - *Anas platyrhynchos*
  - *Aptenodytes forsteri*
  - *Calypte anna*
  - *Corvus brachyrhynchos*
  - *Dryobates pubescens*
  - *Falco tinnunculus*
  - *Ficedula albicollis*
  - *Gallus gallus*
  - *Meleagris gallopavo*
  - *Taeniopygia guttata*
- **Closest relative proteome**: *Trogon surrucura*
- **Curated avian proteins**: UniProt reviewed Aves entries only

---

### 10.3 Dataset 1: Broad Aves protein evidence (PROT_COMBO)

This dataset was designed to maximize protein coverage across birds.

#### Concatenation
```bash
cat \
  OrthoDB12_Aves.fa \
  Amazona_aestiva.fasta \
  Anas_platyrhynchos.fasta \
  Aptenodytes_forsteri.fasta \
  Calypte_anna.fasta \
  Corvus_brachyrhynchos.fasta \
  Dryobates_pubescens.fasta \
  Falco_tinnunculus.fasta \
  Ficedula_albicollis.fasta \
  Gallus_gallus.fasta \
  Meleagris_gallopavo.fasta \
  Taeniopygia_guttata.fasta \
  Trogonsur_protein.fa \
  uniprotkb_aves_AND_reviewed_true_2026_01_28.fasta \
  > PROT_COMBO.raw.faa
```

#### Filtering
```bash
# Removal of exact duplicate sequences
seqkit rmdup -s PROT_COMBO.raw.faa > PROT_COMBO.rmdup.faa

# Removal of short proteins (< 50 aa)
seqkit seq -m 50 PROT_COMBO.rmdup.faa > PROT_COMBO.rmdup.min50aa.faa

#Protein count checks
grep -c '^>' PROT_COMBO.raw.faa
grep -c '^>' PROT_COMBO.rmdup.faa
grep -c '^>' PROT_COMBO.rmdup.min50aa.faa
```

Final dataset size: ~5.2 million proteins

---

### 10.4 Dataset 2: Curated avian proteomes (PROT_BIRDS)

This dataset was designed to reduce noise while retaining phylogenetically relevant, high-quality protein evidence.

#### Concatenation
```bash
cat \
  Amazona_aestiva.fasta \
  Anas_platyrhynchos.fasta \
  Aptenodytes_forsteri.fasta \
  Calypte_anna.fasta \
  Corvus_brachyrhynchos.fasta \
  Dryobates_pubescens.fasta \
  Falco_tinnunculus.fasta \
  Ficedula_albicollis.fasta \
  Gallus_gallus.fasta \
  Meleagris_gallopavo.fasta \
  Taeniopygia_guttata.fasta \
  Trogonsur_protein.fa \
  uniprotkb_aves_AND_reviewed_true_2026_01_28.fasta \
  > PROT_BIRDS.raw.faa
```

#### Filtering
```bash
# Removal of exact duplicate sequences
seqkit rmdup -s PROT_BIRDS.raw.faa > PROT_BIRDS.rmdup.faa

# Removal of short proteins (< 50 aa)
seqkit seq -m 50 PROT_BIRDS.rmdup.faa > PROT_BIRDS.rmdup.min50aa.faa

# Protein count checks
grep -c '^>' PROT_BIRDS.raw.faa
grep -c '^>' PROT_BIRDS.rmdup.faa
grep -c '^>' PROT_BIRDS.rmdup.min50aa.faa
```

Final dataset size: ~194,000 proteins

---

### 10.5 BRAKER Protein Evidence Comparison

#### Annotation Feature Statistics

| Feature        | PROT_BIRDS | PROT_COMBO |
|----------------|------------|------------|
| Genes          | 14,087     | 15,094     |
| Transcripts    | 17,606     | 18,986     |
| mRNA           | 1,718      | 2,013      |
| CDS            | 168,207    | 175,921    |
| Exons          | 168,207    | 175,921    |
| Introns        | 150,601    | 156,935    |
| Start codons   | 17,589     | 18,954    |
| Stop codons    | 17,602     | 18,980    |


#### BUSCO Results (aves_odb10, n = 8,338)

| Metric                         | PROT_BIRDS | PROT_COMBO |
|--------------------------------|------------|------------|
| Complete BUSCOs (C)            | 6,752 (81.0%) | 6,805 (81.6%) |
| ├─ Single-copy (S)             | 5,416 (65.0%) | 5,460 (65.5%) |
| ├─ Duplicated (D)              | 1,336 (16.0%) | 1,345 (16.1%) |
| Fragmented BUSCOs (F)          | 121 (1.5%)   | 161 (1.9%)   |
| Missing BUSCOs (M)             | 1,465 (17.6%) | 1,372 (16.5%) |
| Total BUSCO groups searched    | 8,338        | 8,338        |

**PROT_COMBO** yields a modest increase in gene and transcript counts, and BUSCO completeness improves slightly with **PROT_COMBO**, mainly due to fewer missing BUSCOs.

---

## 11. Technical validation: synteny (MCScanX) between resplendent quetzal and trogon surrucura

This workflow generates a pairwise synteny/collinearity analysis strictly for **technical validation** of assembly/annotation structure. We infer putative 1:1 ortholog candidates using **reciprocal best hits (RBH)** and then detect collinear blocks with **MCScanX**.

#### Inputs
**Quetzal (Pharomachrus mocinno)**
- Proteins (canonical, stop-codon cleaned): `functional_annotation/inputs/quetzal.canonical.clean.faa`
- Gene coordinates: `braker.gtf` (gene rows contain the gene ID as field 9, e.g. `g12723`)

**Trogon surrucura**
- Proteins (canonical): `/path/to/trogon.canonical.faa`
- Gene coordinates: `/path/to/trogon.gtf` or `/path/to/trogon.gff3`

#### Output
- `synteny_mcscanx/QTTS.collinearity` (MCScanX collinear blocks)
- `synteny_mcscanx/QTTS.summary.tsv` (counts for technical validation)
- Intermediate: RBH pairs + formatted `.gff`/`.blast`

#### Results

![Comparative karyotype alignment](Figures/karyotype_alllabels_page-0001.jpg)

The figure illustrates the comparative karyotype alignment between *P. mocinno* and *T. surrucura*, based on the gene annotations of both genome assemblies. Chromosomes are ordered and labeled sequentially for each species, highlighting conserved syntenic blocks and revealing large-scale structural relationships between homologous genomic regions. The connecting ribbons represent correspondences inferred from annotated gene positions, allowing the visualization of chromosomal collinearity as well as potential rearrangements between the two genomes. The figure also highlights the presence of unplaced scaffolds in *P. mocinno* and the identification of the Z chromosome in *T. surrucura*, providing additional insights into assembly structure and sex chromosome differentiation. Overall, this annotation-based karyotype comparison provides a genome-wide perspective of chromosomal conservation and structural divergence between the two species.


---

### 11.1 Cluster script (build inputs → RBH → MCScanX → summary)

Save as: `run_synteny_mcscanx.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail

# =========================
# USER SETTINGS (EDIT THESE)
# =========================
WORKDIR="/lisc/data/scratch/botany/cauz/Genome_Annotation/Braker3/test/Quetzal/Annot_out_RNA_vertebrata/braker"
OUTDIR="${WORKDIR}/synteny_mcscanx"
PREFIX="QTTS"   # prefix for MCScanX input/output files

# Quetzal
Q_FAA="${WORKDIR}/functional_annotation/inputs/quetzal.canonical.clean.faa"
Q_GTF="${WORKDIR}/braker.gtf"

# Trogon (EDIT)
T_FAA="/path/to/trogon.canonical.faa"
T_GTF="/path/to/trogon.gtf"   # can be GTF or GFF3 (see notes below)

# Compute resources
THREADS=16
TOPHITS=5        # keep top hits per query for RBH filtering

# =========================
# MODULES (ADAPT TO CLUSTER)
# =========================
# Example (edit based on module avail):
# module purge
# module load DIAMOND/2.x
# module load MCScanX/...
# If MCScanX is not available as a module, ask your HPC admin or install locally.

# =========================
# SETUP
# =========================
mkdir -p "$OUTDIR"/{01_inputs,02_diamond,03_mcscanx,logs,tmp}

echo "[INFO] OUTDIR: $OUTDIR"

# =========================
# 2) Prepare gene coordinate files in MCScanX .gff-like format
#    MCScanX expects tab-delimited: <chr/scaffold>  <gene_id>  <start>  <end>
#    We'll prefix scaffold names to avoid collisions between species.
# =========================

echo "[INFO] Extracting Quetzal gene coordinates from GTF..."
awk -F'\t' 'BEGIN{OFS="\t"}
  $3=="gene" {
    scaf="Q|" $1; gid=$9; gsub(/^[[:space:]]+|[[:space:]]+$/,"",gid);
    print scaf, gid, $4, $5
  }' "$Q_GTF" \
  > "$OUTDIR/01_inputs/Q.gff"

echo "[INFO] Extracting Trogon gene coordinates..."
# Case A: trogon is GTF with gene id as field 9 (like your BRAKER GTF)
# If your trogon file is standard GTF/GFF3 with attributes, see Notes below.
awk -F'\t' 'BEGIN{OFS="\t"}
  $3=="gene" {
    scaf="T|" $1; gid=$9; gsub(/^[[:space:]]+|[[:space:]]+$/,"",gid);
    print scaf, gid, $4, $5
  }' "$T_GTF" \
  > "$OUTDIR/01_inputs/T.gff"

cat "$OUTDIR/01_inputs/Q.gff" "$OUTDIR/01_inputs/T.gff" \
  > "$OUTDIR/03_mcscanx/${PREFIX}.gff"

# =========================
# 3) Prepare protein FASTA headers so they map to gene IDs (not transcript IDs)
#    Quetzal proteins are like g12723.t1; MCScanX gene IDs are g12723
#    We'll strip ".tX" to create gene-level IDs for both species.
# =========================

echo "[INFO] Normalizing FASTA headers to gene IDs..."
python3 - <<'PY'
from pathlib import Path
import re

def norm_faa(in_faa, out_faa, species_prefix):
    out = Path(out_faa)
    with open(in_faa) as f, out.open("w") as w:
        for line in f:
            if line.startswith(">"):
                h = line[1:].strip().split()[0]
                # strip transcript suffix like .t1
                h = re.sub(r"\.t\d+$", "", h)
                w.write(f">{h}\n")
            else:
                w.write(line)
    # prepend species prefix later by scaffold prefixing; gene IDs must match .gff

WORKDIR = Path(r"""'"$OUTDIR"'/01_inputs""")
norm_faa(r"""'"$Q_FAA"'""", str(WORKDIR/"Q.geneid.faa"), "Q")
norm_faa(r"""'"$T_FAA"'""", str(WORKDIR/"T.geneid.faa"), "T")
PY

# =========================
# 4) DIAMOND searches (Q vs T and T vs Q), then compute RBH
# =========================

echo "[INFO] Building DIAMOND db for Trogon..."
diamond makedb --in "$OUTDIR/01_inputs/T.geneid.faa" -d "$OUTDIR/02_diamond/T.dmnd" >/dev/null

echo "[INFO] Running DIAMOND Q->T..."
diamond blastp \
  -q "$OUTDIR/01_inputs/Q.geneid.faa" \
  -d "$OUTDIR/02_diamond/T.dmnd" \
  -o "$OUTDIR/02_diamond/Q_vs_T.tsv" \
  --outfmt 6 qseqid sseqid evalue bitscore \
  --threads "$THREADS" \
  --max-target-seqs "$TOPHITS" \
  --evalue 1e-5

echo "[INFO] Building DIAMOND db for Quetzal..."
diamond makedb --in "$OUTDIR/01_inputs/Q.geneid.faa" -d "$OUTDIR/02_diamond/Q.dmnd" >/dev/null

echo "[INFO] Running DIAMOND T->Q..."
diamond blastp \
  -q "$OUTDIR/01_inputs/T.geneid.faa" \
  -d "$OUTDIR/02_diamond/Q.dmnd" \
  -o "$OUTDIR/02_diamond/T_vs_Q.tsv" \
  --outfmt 6 qseqid sseqid evalue bitscore \
  --threads "$THREADS" \
  --max-target-seqs "$TOPHITS" \
  --evalue 1e-5

echo "[INFO] Computing reciprocal best hits (RBH)..."
python3 - <<'PY'
from pathlib import Path
import math

qvt = Path(r"""'"$OUTDIR"'/02_diamond/Q_vs_T.tsv""")
tvq = Path(r"""'"$OUTDIR"'/02_diamond/T_vs_Q.tsv""")
out = Path(r"""'"$OUTDIR"'/03_mcscanx/QTTS.rbh.tsv""")

def best_hits(path):
    # best by bitscore, tie-break by lowest evalue
    best = {}
    with path.open() as f:
        for line in f:
            q,s,e,b = line.rstrip("\n").split("\t")
            e = float(e)
            b = float(b)
            if q not in best:
                best[q] = (s,e,b)
            else:
                s0,e0,b0 = best[q]
                if (b > b0) or (b == b0 and e < e0):
                    best[q] = (s,e,b)
    return best

A = best_hits(qvt)  # Q -> best T
B = best_hits(tvq)  # T -> best Q

rbh = []
for q,(t,e1,b1) in A.items():
    if t in B and B[t][0] == q:
        # output format: q  t  evalue  bitscore
        rbh.append((q,t,e1,b1))

with out.open("w") as w:
    for q,t,e,b in rbh:
        w.write(f"{q}\t{t}\t{e}\t{b}\n")

print(f"RBH pairs: {len(rbh)}")
print(f"Wrote: {out}")
PY

# =========================
# 5) Create MCScanX .blast file
#    For MCScanX: tab-delimited: gene1 gene2 evalue score
#    Use RBH file as the similarity graph (strict 1:1 ortholog candidates).
# =========================

cp "$OUTDIR/03_mcscanx/QTTS.rbh.tsv" "$OUTDIR/03_mcscanx/${PREFIX}.blast"

# =========================
# 6) Run MCScanX
# =========================
echo "[INFO] Running MCScanX..."
cd "$OUTDIR/03_mcscanx"
MCScanX "$PREFIX"

# =========================
# 7) Summarize collinearity for technical validation
# =========================
echo "[INFO] Summarizing MCScanX collinearity..."
python3 - <<'PY'
from pathlib import Path
import re

col = Path(r"""'"$OUTDIR"'/03_mcscanx/'"$PREFIX"'.collinearity""")
out = Path(r"""'"$OUTDIR"'/03_mcscanx/'"$PREFIX"'.summary.tsv""")

blocks = 0
pairs = 0

# MCScanX collinearity file contains blocks separated by "## Alignment"
# and gene pairs as lines with two gene IDs.
pair_re = re.compile(r"^\s*\d+\-\s*\d+:\s+(\S+)\s+(\S+)\s+")

with col.open() as f:
    for line in f:
        if line.startswith("## Alignment"):
            blocks += 1
        m = pair_re.search(line)
        if m:
            pairs += 1

with out.open("w") as w:
    w.write("metric\tvalue\n")
    w.write(f"collinear_blocks\t{blocks}\n")
    w.write(f"collinear_gene_pairs\t{pairs}\n")

print(f"Blocks: {blocks}")
print(f"Gene pairs: {pairs}")
print(f"Wrote: {out}")
PY

echo "[DONE] Outputs:"
echo "  - $OUTDIR/03_mcscanx/${PREFIX}.collinearity"
echo "  - $OUTDIR/03_mcscanx/${PREFIX}.summary.tsv"
```
--- 

### 11.2 From MCScanX outputs → JCVI karyotype figure 

From the MCScanX run (prefix `QT_TS`), we assume the following exist:

-   `QT_TS.gff` (MCScanX gene coordinates; scaffold gene start end)
-   `QT_TS.blast` (RBH/BLAST-like pairs; gene1 gene2 evalue bitscore)
-   `QT_TS.collinearity` (MCScanX collinearity blocks)

Additionally, plotting uses BED files with original genome scaffold
names:

-   `quetzal.bed`
-   `trogon.bed`


#### Convert MCScanX `.collinearity` → JCVI anchors

``` bash
python3 - <<'PY'
from pathlib import Path
import re

col = Path("QT_TS.collinearity")
out = Path("quetzal.trogon.anchors")

pair_re = re.compile(r"^\s*\d+\-\s*\d+:\s+(\S+)\s+(\S+)\s+")

blocks = 0
pairs = 0

with col.open() as f, out.open("w") as w:
    for line in f:
        if line.startswith("## Alignment"):
            if blocks > 0:
                w.write("###\n")
            blocks += 1
            continue

        m = pair_re.search(line)
        if m:
            g1, g2 = m.group(1), m.group(2)
            g1 = g1.replace("Q|", "").replace("T|", "")
            g2 = g2.replace("Q|", "").replace("T|", "")
            w.write(f"{g1}\t{g2}\n")
            pairs += 1

print(f"[OK] Wrote {out}")
print(f"[OK] Blocks: {blocks}")
print(f"[OK] Pairs:  {pairs}")
PY
```


#### Create numeric-label BEDs

**Create `trogon.num.bed`** 

``` bash
python3 - <<'PY'
from pathlib import Path

inp = Path("trogon.bed")
out = Path("trogon.num.bed")

allowed = {str(i) for i in range(1, 31)} | {"Z"}

with inp.open() as f, out.open("w") as w:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        seqid, start, end, gid = line.rstrip("\n").split("\t")[:4]
        if seqid not in allowed:
            continue
        w.write("\t".join([seqid, start, end, gid]) + "\n")
PY
```


#### Identify Quetzal unplaced scaffold matching Trogon chr24

``` bash
python3 - <<'PY'
from pathlib import Path
from collections import Counter

qbed = Path("quetzal.bed")
tbed = Path("trogon.bed")
anchors = Path("quetzal.trogon.anchors")

q_gene2seq = {}
with qbed.open() as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        seqid, start, end, gid = line.rstrip("\n").split("\t")[:4]
        q_gene2seq[gid] = seqid

t_gene2seq = {}
with tbed.open() as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        seqid, start, end, gid = line.rstrip("\n").split("\t")[:4]
        t_gene2seq[gid] = seqid

def is_q(g):
    return g.startswith("g") and g[1:].isdigit()

counts = Counter()
unplaced_counts = Counter()

with anchors.open() as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        a, b = line.split()[:2]

        if is_q(a) and not is_q(b):
            qg, tg = a, b
        elif is_q(b) and not is_q(a):
            qg, tg = b, a
        else:
            continue

        if t_gene2seq.get(tg) != "24":
            continue

        qseq = q_gene2seq.get(qg)
        if not qseq:
            continue

        counts[qseq] += 1
        if "unplaced" in qseq.lower():
            unplaced_counts[qseq] += 1

print(unplaced_counts.most_common(1))
PY
```

**Build `quetzal.num.bed`**  

``` bash
python3 - <<'PY'
from pathlib import Path
import re

inp = Path("quetzal.bed")
out = Path("quetzal.num.bed")

U1 = "Pmoc_unplaced000001"
U2 = "Pmoc_unplaced000002"

def map_seqid(s):
    m = re.fullmatch(r"Pmoc_scaf(\d{4})", s)
    if m:
        n = int(m.group(1))
        if 1 <= n <= 31 and n != 30:
            return str(n)
    if s == U1:
        return "unplaced1"
    if s == U2:
        return "unplaced2"
    return None

with inp.open() as f, out.open("w") as w:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        seqid, start, end, gid = line.rstrip("\n").split("\t")[:4]
        new = map_seqid(seqid)
        if new:
            w.write("\t".join([new, start, end, gid]) + "\n")
PY
```

#### Write `seqids`

``` bash
q=$(python3 - <<'PY'
nums = [str(i) for i in range(1,32) if i != 30]
print(",".join(nums + ["unplaced1","unplaced2"]))
PY
)

t=$(python3 - <<'PY'
print(",".join([str(i) for i in range(1,31)] + ["Z"]))
PY
)

printf "%s\n%s\n" "$q" "$t" > seqids
```

#### Edit Layout file

    .70, .05, .95, 0, , P. mocinno, top, quetzal.num.bed, center
    .30, .05, .95, 0, , T. surrucura, bottom, trogon.num.bed, center

    e, 0, 1, quetzal.trogon.blocks.simple


#### Rebuild `blocks.simple`

``` bash
python3 rebuild_blocks.py
```

(See repository script for full implementation.)


#### Run JCVI

``` bash
python -m jcvi.graphics.karyotype \
  --keep-chrlabels \
  --chrstyle=rect \
  --figsize=30x14 \
  --dpi=600 \
  -o karyotype_alllabels.pdf \
  seqids layout
```

#### JCVI Source Patch

**Disable label downsampling:** 

``` python
if (not keep_chrlabels) and nseqids >= 2 * MaxSeqids and (i + 1) % step != 0:
    continue
```

**Increase chromosome label size:** 

    TextCircle(... size=18 ...)

**Fix species label clipping:** 

    ha="right"

---
## 12. Structural validation of Quetzal contig_1125 using minimap2 alignments and ONT read support

This workflow was used to assess whether **Quetzal contig_1125** (`contig_1125_LN:i:28940123_RC:i:700486_XC:f:1.000000_1`) represents a chimeric misassembly
or a biologically plausible breakpoint when compared to the *Trogon surrucura* genome. All analyses are based on nucleotide-level whole-genome alignments (minimap2) and
ONT read mappings.

### 12.1 Input data

#### Assemblies
```bash
# Quetzal assembly
QFA=/mnt/data/OVHlord/annotation/scaffolding/Quetzal/syri/pharomachrus.ragtag.renamed.fasta

# Trogon assembly
TFA=/mnt/data/OVHlord/annotation/scaffolding/Quetzal/GCA_020746105.1_bTroSur1.noW.fna

# Contig under investigation
CTG='contig_1125_LN:i:28940123_RC:i:700486_XC:f:1.000000_1'

Indexing
samtools faidx "$QFA"
samtools faidx "$TFA"
```

### 12.2 Whole-contig alignment to Trogon genome (minimap2)
```bash
samtools faidx "$QFA" "$CTG" > contig_1125.fa
minimap2 -t 16 -x asm5 "$TFA" contig_1125.fa > c1125_vs_trogon.paf
```

Extract long alignment blocks (≥50 kb) and sort by query coordinates:
```bash
awk 'BEGIN{OFS="\t"}
($11>=50000){
  print $3,$4,$6,$5
}' c1125_vs_trogon.paf | sort -k1,1n > c1125.blocks.tsv
```
Result: two dominant blocks on different Trogon scaffolds:
- CM036622.1 (Quetzal positions up to ~32.64 Mb)
- CM036624.1 (Quetzal positions from ~32.65 Mb onward)

### 12.3 Breakpoint definition
```bash
awk '
$6=="CM036622.1"{ if($4>max1) max1=$4 }
$6=="CM036624.1"{ if(min2==0 || $3<min2) min2=$3 }
END{
  print "max_qend_CM036622.1:", max1;
  print "min_qstart_CM036624.1:", min2;
  print "gap:", min2-max1;
  print "QB_midpoint:", int((max1+min2)/2);
}' c1125_vs_trogon.paf
```

Final values:
- End of CM036622.1 block: 32,643,015
- Start of CM036624.1 block: 32,647,979
- **Gap size: ~4,964 bp**
- Breakpoint midpoint (QB): 32,645,497

### 12.4 ONT read mapping statistics

ONT reads mapped to the Quetzal assembly:

```bash
BAM=/mnt/data/OVHlord/annotation/scaffolding/Quetzal/mapping/ont_vs_ragtag.bam

Read length distribution in the region:

samtools view region.bam | awk '{print length($10)}' | sort -n | awk '
{a[NR]=$1}
END{
  print "N_reads:",NR;
  print "Median_len:", a[int(NR*0.5)];
  print "P90_len:", a[int(NR*0.9)];
}'
```

Observed:
- Median read length ≈ 812 bp
- P90 ≈ 4.1 kb

### 12.5 Read spanning across the breakpoint

Given realistic read lengths, spanning was evaluated at ±200, ±500, and ±1000 bp.

```bash 
QB=32645497
WIN=20000
REG="${CTG}:$((QB-WIN))-$((QB+WIN))"


for SPAN in 200 500 1000; do
  samtools view "$BAM" "$REG" \
  | awk -v bp=$QB -v s=$SPAN '
    {
      pos=$4; cigar=$6;
      len=0; tmp=cigar;
      while (match(tmp, /([0-9]+)([MDN=X])/, a)) {
        len+=a[1]; tmp=substr(tmp, RSTART+RLENGTH)
      }
      end=pos+len;
      if (pos < bp-s && end > bp+s) c++
    }
    END{print "Reads_spanning_±"s":", c+0}
  '
done
```

Results:
±200 bp: 68 reads
±500 bp: 58 reads
±1000 bp: 39 reads

A nearby control position (32,700,000) showed comparable spanning support,
**indicating no local loss of continuity at the breakpoint**.

### 12.6 Soft-clipping profile
```bash
samtools view "$BAM" "$REG" \
| awk '
function clip(cg,a){
  l=0;r=0;
  if(match(cg,/^([0-9]+)S/,a)) l=a[1];
  if(match(cg,/([0-9]+)S$/,a)) r=a[1];
  return (l>=50 || r>=50);
}
clip($6){
  bin=int($4/10)*10;
  c[bin]++
}
END{for(b in c) print b, c[b]}' \
| sort -k2,2nr | head -20
```

Soft-clipped reads are distributed across the region rather than forming a
single sharp breakpoint, inconsistent with a classical chimeric misjoin.

### 12.7 Extraction of the non-aligning gap sequence
```bash
G1=32643016
G2=32647978
samtools faidx "$QFA" "${CTG}:${G1}-${G2}" > gap_4964bp.fa

Composition check:

tail -n +2 gap_4964bp.fa | tr -d '\n' | awk '
{
  len=length($0);
  n=gsub(/[Nn]/,"");
  print "length:",len,"Ns:",n,"N_fraction:",(len? n/len:0)
}'
```

Result:
- Length ≈ 4963 bp
- Ns = 0 (real assembled sequence)

### 12.8 Mapping the gap sequence to the Trogon genome
```bash
minimap2 -t 8 -x asm10 "$TFA" gap_4964bp.fa > gap_vs_trogon.paf
```

Representative alignments:
```bash
head gap_vs_trogon.paf
```

Observation:
- Multiple short, partial alignments
- Hits to multiple Trogon contigs/scaffolds
- No long unique alignment bridging CM036622.1 and CM036624.1

**This pattern is consistent with repeat-rich or low-complexity sequence.**  

**Conclusion:** ONT reads span the inferred breakpoint at levels comparable to nearby regions. No sharp soft-clipping wall is observed. 
The ~5 kb gap contains real sequence (no Ns). The gap aligns to the Trogon genome as multiple short, dispersed fragments, consistent with repetitive sequence. 
There is no evidence supporting a chimeric misassembly of
Quetzal contig_1125. The observed whole-genome alignment split between
CM036622.1 and CM036624.1 is best explained by repeat-rich sequence and/or
structural or scaffolding differences between the Quetzal and Trogon assemblies.
Accordingly, contig_1125 was retained intact in the final assembly, and synteny
across this locus is interpreted cautiously.

---
