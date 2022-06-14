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

#read in sleep study file
sleep_file <- fread("/home/jacob/continuous-1160-both_sexes.tsv.gz")

#removing snps where any of the 4 required categories for 2sample MR are N/A
some_na <- sleep_file %>%
  filter_at(vars(chr, pos, se_meta, beta_meta), all_vars(!is.na(.)))

#sorting significant snps with pval threshold (5228 smaller than 5E-8)
significant <- some_na[pval_meta <= 5E-8]

#renaming column names to match 2sample requirements
significant <- dplyr::rename(significant, beta = beta_meta)
significant <- dplyr::rename(significant, se = se_meta)
significant <- dplyr::rename(significant, pval = pval_meta)
significant <- dplyr::rename(significant, effect_allele = ref)
significant <- dplyr::rename(significant, other_allele = alt)
significant <- dplyr::rename(significant, eaf = af_meta)

# taking chr and pos information from top5220 so that we can convert to rsid (which is required by 2sampleMR)
sleep_clip <- significant[,c(1,2)] 
positions <- GPos(seqnames = sleep_clip$chr, pos = sleep_clip$pos)
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)

#adding snp rsid column back into sleep clip dataframe
rsid_mod <- data.frame(pos = start(my_snps), mcols(my_snps)[,c("RefSNP_id", "ref_allele", "alt_alleles")])

#merging and renaming rsid info into original significant modified dataframe with all variables
alltogethersig <- left_join(significant, rsid_mod, by = "pos") 
#alltogether5220 <- rename(alltogether5220, c("SNP" = "RefSNP_id"))
alltogethersig <- dplyr::rename(alltogethersig, SNP = RefSNP_id)

#2 sample MR steps can now be run
#formatting UKBB sleep study as exposure
alltogethersigE <- format_data(alltogethersig, type="exposure")
exposure_dat <- extract_instruments(alltogethersigE)
sleep_exp_dat <- clump_data(alltogethersigE)

#formatting ieugwas heart study as outcome
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

# reverse reverse
# Now using CHD ieugwas as exposure and PanUKBB sleep study as outcome
# THIS AREA IS UNDER CONSTRUCTION - code is non functional and errors out, still trying to figure it out

exposure_dat_CHD <- extract_instruments(outcomes = "ieu-a-9")

outcome_dat_sleep <- format_data(alltogethersig, type="outcome")
outcome_dat_sleep_c <- extract_outcome_data(snps=exposure_dat_CHD$SNP, outcomes=outcome_dat_sleep)

dat2 <- harmonise_data(exposure_dat_CHD, outcome_dat_sleep)
res2 <- mr(dat2)

#scatter plot 2
p2 <- mr_scatter_plot(res2, dat2)
p2[[1]]
