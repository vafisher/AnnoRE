#' Calculate random effect variances from annotations and partion of heritability.
#' 
#' @param gwas_locus data frame containing GWAS summmary statistics for fine-mapping locus.
#'  Columns SNP (variant ID, any character string format), b (effect estimate), 
#'  se (standard error), N (sample size), Freq1 (alternate allele frequency) are required.
#'  All SNPs are assumed to be bi-allelic.
#' @param annotation data frame with functional annotations including the fine mapping locus.
#'  Must include columns CHR, BP, SNP, CM, and additional columns named for the annotation categories
#'  (i.e. format of ENCODE annotations at https://alkesgroup.broadinstitute.org/LDSCORE/)
#' @param h2partition data frame with estimated heritability partition
#'  Must include column "Coefficient" with estimates, and "Category" with names of annotation categories.
#'  (i.e. format of ldsc output with --h2 option)
#' @return data frame with all columns in gwas_locus plus RE variance Vbeta and its variance varVbeta.
#'  Rows represent SNPs present in both gwas_locus and annotation.
#'  
#' @export
#'  
annotateREvariance <- function(gwas_locus, annotation, h2partition) {
  #check all the inputs have the right columns
  required_column_names <- c("CHR", "BP", "SNP", "CM")
  if (sum(required_column_names %in% names(annotation)) < length(required_column_names)) {
    missing_cols <- required_column_names[which(!(required_column_names %in%
                                                    names(annotation)))]
    stop(paste("input annotation data frame does not contain required variable(s)",
               paste(missing_cols, collapse = ", ")))
  }
  
  anno_names <- names(annotation)
  anno_names <- anno_names[which(!(anno_names %in% required_column_names))]
  
  anno_locus <- annotation[which(annotation$SNP %in% gwas_locus$SNP),]
  

  
  
  # calculate RE variances as inner product of annotations and heritabilities
  anno_out <- anno_locus[, required_column_names]
  
  # set negative heritability estimates to zero
  coef_pos <- pmax(h2partition$Coefficient, 0)
  names(coef_pos) <- gsub("_0$", "", h2partition$Category)
  coef_pos <- coef_pos[anno_names]
  anno_out$Vbeta <- as.matrix(anno_locus[, anno_names]) %*% coef_pos

  # heritability variances
  coef_var <- h2partition$Coefficient_std_error^2*as.numeric(h2partition$Coefficient <= 0)
  names(coef_var) <- gsub("_0$", "", h2partition$Category)
  coef_var <- coef_var[anno_names]
  anno_out$varVbeta <- as.matrix(anno_locus[, anno_names]) %*% coef_var
  
  if (F) {
    
    varbj <- function(k) 
      sum(loc[k,5:57] * coef_pos)
    loc$Vbeta<-sapply(1:nrow(loc),varbj)
    
    coef_var<-h2partition$Coefficient_std_error^2*as.numeric(coef_pos==0)
    var_var<-function(k) sum(loc[k,5:57] * coef_var)
    loc$varVbeta<-sapply(1:nrow(loc),var_var)
  }
  
  # add these to GWAS summary file 
  gwas_anno_out <- merge(gwas_locus, anno_out)
  gwas_anno_out <- gwas_anno_out[order(gwas_anno_out$BP), ]
  return(gwas_anno_out)
}
