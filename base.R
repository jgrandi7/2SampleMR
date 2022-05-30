library(TwoSampleMR)
library(plyr)
library(dplyr)
library(data.table)
library(biomaRt)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)


# List available GWASs
ao <- available_outcomes()

genome <- BSgenome.Hsapiens.UCSC.hg19
all_snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
seqlevelsStyle(genome) <- "NCBI"

system("zcat /home/jacob/continuous-1160-both_sexes.tsv.gz | sort -gk 8 | head -n 10000 > /home/jacob/ukbbtop1000.tsv")

sleep_file <- fread("/home/jacob/continuous-1160-both_sexes.tsv.gz")
no_na <- na.omit(sleep_file) 

some_na <- sleep_file %>%
  filter_at(vars(chr, pos, se_meta, beta_meta), all_vars(!is.na(.)))

top5220 <- head(some_na[order(some_na$pval_meta)], 5220)


top5220 <- rename(top5220, c("beta" = "beta_meta"))
top5220 <- rename(top5220, c("se" = "se_meta"))
top5220 <- rename(top5220, c("pval" = "pval_meta"))
top5220 <- rename(top5220, c("effect_allele" = "ref"))
top5220 <- rename(top5220, c("other_allele" = "alt"))
top5220 <- rename(top5220, c("eaf" = "af_meta"))
sleep_clip <- top5220[,c(1,2)] 
positions <- GPos(seqnames = sleep_clip$chr, pos = sleep_clip$pos)
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)

rsid_mod <- data.frame(pos = start(my_snps), mcols(my_snps)[,c("RefSNP_id", "ref_allele", "alt_alleles")])

alltogether5220 <- left_join(top5220, rsid_mod, by = "pos") 
alltogether5220 <- rename(alltogether5220, c("SNP" = "RefSNP_id"))
alltogether5220e <- format_data(alltogether5220, type="exposure")
exposure_dat <- extract_instruments(alltogether5220e)
sleep_exp_dat <- clump_data(alltogether5220e)
outcome_dat <- extract_outcome_data(snps=sleep_exp_dat$SNP, outcomes="ieu-a-9")
dat <- harmonise_data(sleep_exp_dat, outcome_dat)
res <- mr(dat)
#scatter plot
p1 <- mr_scatter_plot(res, dat)
p1[[1]]
#funnel plot
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p4[[1]]
#mr steiger directionality test
out <- directionality_test(dat)
#reverse reverse
exposure_dat2 <- extract_instruments(outcomes="ieu-a-9")
expog <- format_data(alltogether5220, type="outcome")
exposure_dat2_clump <- clump_data(exposure_dat2)
outcome_dat2 <- extract_outcome_data(snps=exposure_dat2_clump$SNP, outcomes=expog)
dat2 <- harmonise_data(exposure_dat2, outcome_dat2)
res2 <- mr(dat2)
#scatter plot 2
p2 <- mr_scatter_plot(res2, dat2)
p2[[1]]
