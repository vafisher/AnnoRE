#' Calculate random effect variances from annotations and partion of heritability.
#' 
#'  @param gwas_locus data frame containing GWAS summmary statistics for fine-mapping locus.
#'  Columns SNP (variant ID, any character string format), b (effect estimate), 
#'  se (standard error), N (sample size), Freq1 (alternate allele frequency) are required.
#'  All SNPs are assumed to be bi-allelic.
#'  
#'  @param annotation data frame with functional annotations including the fine mapping locus.
#'  Must include columns CHR, BP, SNP, CM, base, and additional columns named for the annotation categories
#'  (i.e. format of ENCODE annotations at https://alkesgroup.broadinstitute.org/LDSCORE/)
#'  
#'  @param h2partition data frame with estimated heritability partition
#'  Must include column "Coefficient" with estimates, and "[annotation name?]" with names of annotation categories.
#'  (i.e. format of ldsc output with --h2 option)
#'  
#'  @returns data frame with all columns in [input summary] plus RE variance Vbeta and its variance varVbeta.
#'  
#'  @export
#'  
annotateREvariance <- function(gwas_locus, annotation, h2partition) {
  #TK: check all the inputs have the right columns
  
  anno_names <- annotation %>%
    select(-CHR, -BP, -SNP, -CM) %>%
    names()
  
  anno_locus <- annotation %>%
    filter(SNP %in% gwas_locus$SNP)
  
  # set negative heritability estimates to zero
  coef_pos <- pmax(h2partition$Coefficient, 0)
  
  # calculate RE variances as inner product of annotations and heritabilities
  
  # add these to GWAS summary file 
}
