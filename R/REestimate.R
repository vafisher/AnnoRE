#' Estimate best linear unbiased predictor of SNP random effects in a GWAS locus
#'
#' @param loc_summary_df A data frame containing summary statistics for the locus.
#' Columns SNP, CHR, BP, b, se, N, Freq1, Vbeta, varVbeta are required.
#' @param trait_h2 A number between 0 and 1 representing total trait heritability
#' @param nominal_p_filter A number between 0 and 1 representing p-value threshold for consideration in fine mapping. 
#' @param loc_ld_matrix A numeric matrix containing linkage disequilibrium $R^2$ for all SNPs in loc_summary_df.
#' @param loc_ld_file A file path for the VCFtools output containing LD for all SNPs in loc_summary_df.
#' Only used if loc_ld_matrix is NULL.
#' @return A data frame matching loc_summary_df with columns containing AnnoRE results appended.
#' @export
REestimate <- function(loc_summary_df, trait_h2, nominal_p_filter = 0.05, 
                       loc_ld_matrix = NULL, loc_ld_file = NULL) {

  # check dimensions of input data
  input_vars <- c("SNP", "CHR", "BP", "b", "se", "N", "Freq1", "Vbeta", "varVbeta")
  if (sum(input_vars %in% names(loc_summary_df)) < length(input_vars)) {
    missing_vars <- input_vars[which(! input_vars %in% names(loc_summary_df))]
    return(paste0("input data frame must include column ", missing_vars))
  }
  #check and process pairwise LD
  if (!is.null(loc_ld_matrix)) {
    if (!is.matrix(loc_ld_matrix)) {
      return("loc_ld_matrix must be a numeric matrix")
    } else if (!(dim(loc_ld_matrix)[1] == dim(loc_ld_matrix)[2] &
                 dim(loc_ld_matrix)[1] == nrow(loc_summary_df))) {
      return("size of ld matrix does not match summary data frame")
    }
  } else { #read LD from vcftools output file and convert to matrix; loc_ld_file = "inst/data/rs3736485.1000G.loc.hap.ld"
    ld_long <- read.table(loc_ld_file, header=T, as.is=T)
    ld_long <- ld_long[which(ld_long$POS1 %in% loc_summary_df$BP &
                               ld_long$POS2 %in% loc_summary_df$BP),]
    loc_ld_matrix <- reshape(ld_long, direction="wide",
                   v.names=c("R.2"),
                   idvar=c("POS1"),
                   timevar=c("POS2"),
                   drop=c("CHR","D","Dprime"))
    loc_ld_matrix <- rbind(cbind(NA,loc_ld_matrix[,-(1:2)]),NA)
    loc_ld_matrix[lower.tri(loc_ld_matrix)]<-t(loc_ld_matrix)[lower.tri(loc_ld_matrix)]
    diag(loc_ld_matrix)<-1
    loc_ld_matrix <- as.matrix(loc_ld_matrix)
  }

  loc_df <- loc_summary_df[which(!is.na(loc_summary_df$Freq1) & 
                                   loc_summary_df$p < nominal_p_filter), ]
  loc_df$Z <- loc_df$b/loc_df$se
  loc_df$D <- 2*loc_df$N*loc_df$Freq1*(1-loc_df$Freq1)
  loc_df$sig2R <- loc_df$D*(loc_df$varVbeta*(loc_df$N-1) + loc_df$b^2)
  loc_df$h2 <- (2*loc_df$b^2*loc_df$Freq1*(1-loc_df$Freq1))/median(loc_df$sig2R, na.rm = T)
  loc_df$D_std = sqrt((2*loc_df$N*loc_df$Freq1*(1-loc_df$Freq1))/loc_df$sig2R)
  
  loc_snp_idx <- which(loc_df$SNP %in% loc_summary_df$SNP)
  loc_ld_matrix <- loc_ld_matrix[loc_snp_idx, loc_snp_idx]

  gwa_h2<-trait_h2/1e6
  locus_factor <- mean(loc_df$h2, na.rm = T)/gwa_h2

  G.vec <- as.numeric(loc_df$Vbeta*locus_factor)
  G.inv <- diag(1/G.vec)
  D.mat <- diag(as.numeric(loc_df$D_std))

  B <- (D.mat %*% loc_ld_matrix %*% D.mat) + G.inv
  B.inv <- MASS::ginv(B)
  
  b.tilde <- B.inv %*% D.mat %*% as.numeric(loc_df$Z)
  var.b.tilde <- diag(B.inv)
  
  Wald.scores <- b.tilde^2/var.b.tilde
  Wald.p <- pchisq(Wald.scores, 1, lower.tail=F)
  
  loc_df <- loc_df[, names(loc_summary_df)]
  loc_df$b_AnnoRE <- b.tilde
  loc_df$chisq_AnnoRE <- Wald.scores
  loc_df$p_AnnoRE <- Wald.p
  
  return(merge(loc_summary_df, loc_df, all.x = TRUE))
}
