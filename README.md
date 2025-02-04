Analysis of COL5A1 Gene Data

Overview

This project involves retrieving, filtering, and analyzing genetic variation data for the COL5A1 gene from the 1000 Genomes Project. The analysis includes sample extraction, Hardy-Weinberg equilibrium testing, linkage disequilibrium (LD) calculations, phylogenetic tree construction, and selective sweep detection.

Data Retrieval

Downloading the COL5A1 Gene Data

tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr9.filtered.shapeit2-duohmm-phased.vcf.gz chr9:134641803-134844843 > /projectnb/anth333/rwhelan3/COL5A1/COL5A1_all_regions.vcf

Explanation:

Extracts genetic variation data for the COL5A1 gene region (chromosome 9, positions 134641803-134844843).

Stores the data in COL5A1_all_regions.vcf.

Downloading the Sample Panel File

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

Explanation:

Retrieves the sample metadata file containing population information for all individuals in the dataset.

Extracting Population-Specific Samples

Creating Subject Lists

grep CHS /projectnb/anth333/rwhelan3/COL5A1/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > /projectnb/anth333/rwhelan3/COL5A1/CHS.samples.list
grep CHB /projectnb/anth333/rwhelan3/COL5A1/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > /projectnb/anth333/rwhelan3/COL5A1/CHB.samples.list
grep JPT /projectnb/anth333/rwhelan3/COL5A1/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > /projectnb/anth333/rwhelan3/COL5A1/JPT.samples.list

Explanation:

Extracts sample IDs for three East Asian populations: CHS (Southern Han Chinese), CHB (Han Chinese in Beijing), JPT (Japanese in Tokyo).

Saves the sample IDs to separate files for filtering.

Filtering Population-Specific Variants

vcftools --vcf COL5A1_all_regions.vcf --keep CHS.samples.list --recode --out COL5A1_CHS.vcf --non-ref-ac-any 1
vcftools --vcf COL5A1_all_regions.vcf --keep CHB.samples.list --recode --out COL5A1_CHB.vcf --non-ref-ac-any 1
vcftools --vcf COL5A1_all_regions.vcf --keep JPT.samples.list --recode --out COL5A1_JPT.vcf --non-ref-ac-any 1

Explanation:

Filters COL5A1 variants for CHS, CHB, and JPT populations.

Keeps only variants with at least one non-reference allele.

Analyzing Ancestral Alleles

tabix -h https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr9.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz 9:134641803-134844843 | bgzip -c > CHS_ANCESTRAL.vcf.gz
vcftools --gzvcf CHS_ANCESTRAL.vcf.gz --keep CHS.samples.list --maf 0.05 --remove-indels --recode-INFO AA --recode --stdout | bgzip -c > CHS_ANCESTRAL_FILTERED.vcf.gz

Explanation:

Downloads ancestral allele data and extracts COL5A1 region.

Filters for minor allele frequency (MAF) ≥ 0.05 and removes indels.

Hardy-Weinberg Equilibrium (HWE) Test

pegas.CHS <- vcfR2genind(COL5A1vcf, sep = "[|/]")
HWE <- hw.test(pegas.CHS, B = 0)
HWE.sig <- HWE %>% as_tibble(HWE) %>% mutate(snp = rownames(HWE)) %>% rename(chi2 = 'chi^2', pval = 'Pr(chi^2 >)') %>% select(snp,chi2,df,pval) %>% filter(pval <= 0.05)

Explanation:

Tests variants for deviation from Hardy-Weinberg equilibrium.

Extracts statistically significant deviations (p ≤ 0.05).

Linkage Disequilibrium (LD) Analysis

COL5A1vcf <- read.vcfR("COL5A1_CHS.vcf")
COL5A1matrix <- vcfR2SnpMatrix(COL5A1vcf)
COL5A1_LD <- snpStats::ld(COL5A1matrix$data, depth = 1616, stats = "R.squared")

Explanation:

Computes pairwise LD using R² for CHS.

Phylogenetic Tree Construction

CHSdna <- vcfR2DNAbin(COL5A1vcf)
dist <- dist.dna(CHSdna, model = "K80")
CHStree <- nj(dist)
plot(CHStree, cex=0.5)

Explanation:

Converts VCF data to DNA sequences.

Builds a neighbor-joining phylogenetic tree.

Multitree Analysis

CHB <- read.vcfR("COL5A1_CHB/COL5A1_CHB.vcf")
JPT <- read.vcfR("COL5A1_JPT/COL5A1_JPT.vcf")
CHBdna <- vcfR2DNAbin(CHB)
JPTdna <- vcfR2DNAbin(JPT)
CHBtree <- nj(dist.dna(CHBdna, model = "K80"))
JPTtree <- nj(dist.dna(JPTdna, model = "K80"))
densitree <- densiTree(list(CHStree, CHBtree, JPTtree), type="phylogram", col=c("red", "green", "blue"))

Explanation:

Constructs and overlays phylogenetic trees for CHS, CHB, and JPT.

Selective Sweep Detection

hap <- data2haplohh(hap_file="CHS_ANCESTRAL_FILTERED.vcf.gz", polarize_vcf = FALSE)
ehh <- scan_hh(hap, limhaplo=2, limehh=0.05, limehhs=0.05)
ihs <- ihh2ihs(ehh)
cr.se <- calc_candidate_regions(ihs, threshold=2, window_size=3000, overlap=300)
manhattanplot(ihs, pval=TRUE, threshold=2, main="iHS")

Explanation:

Detects regions under selection using integrated haplotype score (iHS).

Visualizes significant sweeps with a Manhattan plot.

Conclusion

This analysis investigates genetic variation in COL5A1 across East Asian populations, highlighting HWE deviations, LD structure, phylogenetic relationships, and selective sweeps.
