# Terminal Code

## COL5A1 gene data:
```bash
tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr9.filtered.shapeit2-duohmm-phased.vcf.gz chr9:134641803-134844843 > /projectnb/anth333/rwhelan3/COL5A1/COL5A1_all_regions.vcf
```
## Call Samples
```bash
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```
## Subject Lists
```bash
grep CHS /projectnb/anth333/rwhelan3/COL5A1/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > /projectnb/anth333/rwhelan3/COL5A1/CHS.samples.list
grep CHB /projectnb/anth333/rwhelan3/COL5A1/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > /projectnb/anth333/rwhelan3/COL5A1/CHB.samples.list
grep JPT /projectnb/anth333/rwhelan3/COL5A1/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > /projectnb/anth333/rwhelan3/COL5A1/JPT.samples.list
```
## Filter 
```bash
vcftools --vcf /projectnb/anth333/rwhelan3/COL5A1//COL5A1_all_regions.vcf --keep /projectnb/anth333/rwhelan3/COL5A1/COL5A1_JPT/JPT.samples.list --recode --out /projectnb/anth333/rwhelan3/COL5A1/COL5A1_JPT/COL5A1_JPT.vcf --non-ref-ac-any 1
vcftools --vcf /projectnb/anth333/rwhelan3/COL5A1/COL5A1_all_regions.vcf --keep /projectnb/anth333/rwhelan3/COL5A1/COL5A1_CHS/CHS.samples.list --recode --out /projectnb/anth333/rwhelan3/COL5A1/COL5A1_CHS/COL5A1_CHS.vcf --non-ref-ac-any 1
vcftools --vcf /projectnb/anth333/rwhelan3/COL5A1/COL5A1_all_regions.vcf --keep /projectnb/anth333/rwhelan3/COL5A1/COL5A1_CHB/CHB.samples.list --recode --out /projectnb/anth333/rwhelan3/COL5A1/COL5A1_CHB/COL5A1_CHB.vcf --non-ref-ac-any 1
```
## Ancestral Alleles
```bash
tabix -h https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr9.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz 9:134641803-134844843 | bgzip -c > CHS_ANCESTRAL.vcf.gz
vcftools --gzvcf CHS_ANCESTRAL.vcf.gz --keep CHS.samples.list --maf 0.05 --remove-indels --recode-INFO AA --recode --stdout | bgzip -c > CHS_ANCESTRAL_FILTERED.vcf.gz
```
# R Studio Code
## HWE
```R
pegas.CHS <- vcfR2genind(COL5A1vcf, sep = "[|/]") 
HWE <- hw.test(pegas.CHS, B = 0)
library(tidyverse)
HWE.sig <- HWE %>% 
  as_tibble(HWE) %>% 
  mutate(snp = rownames(HWE)) %>% 
  rename(chi2 = 'chi^2', pval = 'Pr(chi^2 >)') %>% 
  select(snp,chi2,df,pval) %>% 
  filter(pval <= 0.05)
HWE.sig
```
## LD1
```R
COL5A1vcf <- read.vcfR("COL5A1_CHS.vcf")
COL5A1matrix <- vcfR2SnpMatrix(COL5A1vcf)
COL5A1_LD <- snpStats::ld(COL5A1matrix$data, depth = 1616, stats = "R.squared")
R2heatmapCOL5A1 <- LDheatmap(COL5A1matrix$data,
                             genetic.distances=COL5A1matrix$genetic.distances,
                             distances="physical",
                             LDmeasure="r",
                             title="Pairwise LD in COL5A1 for CHS using R^2",
                             add.map=TRUE, add.key=TRUE,
                             geneMapLocation=0.15,
                             SNP.name=c("9:134651496:G:A","9:134654639:C:A"),
                             color=NULL, newpage=TRUE,
                             name="COL5A1 LD Heatmap (R^2)")
```
## Important LD Segments
```R
col5a1_res = BigLD(genofile = "COL5A1_CHS.vcf", LD="r2")
```
## LD2
```R
rgb.palette <- colorRampPalette(rev(c("grey85", "green", "darkgreen")), space = "rgb"))

png(filename="COL5A1_CHS_LD2.png", units="in", res=600, width=6, height=5)

R2heatmapCOL5A1 <- LDheatmap(COL5A1matrix$data,
                             genetic.distances=COL5A1matrix$genetic.distances,
                             distances="physical",
                             LDmeasure="r",
                             title="Pairwise LD in COL5A1 for CHS using R^2",
                             add.map=TRUE, add.key=TRUE,
                             geneMapLocation=0.15,
                             SNP.name=c("9:134651496:G:A","9:134654639:C:A"),
                             color=NULL, newpage=TRUE,
                             name="COL5A1 LD Heatmap (R^2)")

LDheatmap.highlight(R2heatmapCOL5A1, i = 5, j = 375, fill = "NA", col = "gold", lwd = 2, lty = 1, flipOutline=FALSE, crissCross = FALSE)
LDheatmap.highlight(R2heatmapCOL5A1, i = 402, j = 419, fill = "NA", col = "gold3", lwd = 2, lty = 1, flipOutline=FALSE, crissCross = FALSE)
LDheatmap.highlight(R2heatmapCOL5A1, i = 438, j = 761, fill = "NA", col = "gold4", lwd = 2, lty = 1, flipOutline=FALSE, crissCross = FALSE)
LDheatmap.highlight(R2heatmapCOL5A1, i = 1089, j = 1302, fill = "NA", col = "goldenrod", lwd = 2, lty = 1, flipOutline=FALSE, crissCross = FALSE)

dev.off()
```
## Phylogenetic tree
```R
CHSdna <- vcfR2DNAbin(COL5A1vcf)
View the variance as a table:
par(mar=c(5,8,4,2))
ape::image.DNAbin(CHSdna[,ape::seg.sites(CHSdna)], cex.lab=0.5)
Now as a tree:
dist <- dist.dna(CHSdna, model = "K80")
heatmap <- as.data.frame(as.matrix(dist))
table.paint(heatmap, cleg=0, clabel.row=.5, clabel.col=.5)
CHStree <- nj(dist)
plot(CHStree, cex=0.5)
To get the pretty tree:
ggtree(CHStree, layout="daylight") + geom_tiplab(aes(angle=angle), color='blue')
```
## Multitree
```R
CHB <- read.vcfR("COL5A1_CHB/COL5A1_CHB.vcf")
JPT <- read.vcfR("COL5A1_JPT/COL5A1_JPT.vcf")
CHBdna <- vcfR2DNAbin(CHB)
JPTdna <- vcfR2DNAbin(JPT)

D <- dist.dna(CHBdna, model = "K80")
CHBtree <- nj(D)
D <- dist.dna(JPTdna, model = "K80")
JPTtree <- nj(D)

multitree <- c(CHStree, CHBtree, JPTtree)
densitree <- densiTree(multitree, type="phylogram", col=c("red", "green", "blue"), tip.color=c("red", "green", "blue"), cex = 0.5, width=2, jitter=list(amount=.3, random=FALSE), alpha=1)
```
## Tajima's D
```R
tajima <- tajima.test(CHSdna)
```
## Selective Sweeps
```R
hap <- data2haplohh(hap_file="CHS_ANCESTRAL_FILTERED.vcf.gz", polarize_vcf = FALSE)
ehh <- scan_hh(hap, limhaplo=2, limehh=0.05, limehhs=0.05, maxgap=NA, threads=1)
ihs <- ihh2ihs(ehh, freqbin=0.025)
cr.se <- calc_candidate_regions(ihs, threshold=2, pval=TRUE, window_size=3000, overlap=300, min_n_extr_mrk=1)
print(cr.se, max.print=T)
manhattanplot(ihs, pval=TRUE, threshold=2, pch=16, chr.name="4", cr=cr.se, cr.col="skyblue", main="iHS")
```
