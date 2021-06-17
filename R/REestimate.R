#' Estimate best linear unbiased predictor of SNP random effects in a GWAS locus
#'
#' @param loc_summary_df A data frame containing summary statistics for the locus.
#' Columns b, se, N, Freq1, Vbeta, varVbeta are required.
#' @param trait_h0 A number between 0 and 1 representing total trait heritability
#' @param loc_ld_matrix A numeric matrix containing linkage disequilibrium $R^2$ for all SNPs in loc_summary_df
#' @param loc_ld_file A file path for the VCFtools output containing LD for all SNPs in loc_summary_df.
#' Only used if loc_ld_matrix is NULL.
#' @return A data frame matching loc_summary_df with columns containing AnnoRE results appended.
REestimate <- function(loc_summary_df, trait_h2, loc_ld_matrix = NULL, loc_ld_file = NULL) {

  # check dimensions of input data
  input_vars <- c("SNP", "b", "se", "N", "Freq1", "Vbeta", "varVbeta")
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
    ld_long<-read.table(loc_ld_file, header=T, as.is=T)
    loc_ld_matrix<-reshape(ld_long, direction="wide",
                   v.names=c("R.2"),idvar=c("POS1"),
                   timevar=c("POS2"),
                   drop=c("CHR","D","Dprime"))
    loc_ld_matrix <- rbind(cbind(NA,loc_ld_matrix[,-(1:2)]),NA)
    loc_ld_matrix[lower.tri(loc_ld_matrix)]<-t(loc_ld_matrix)[lower.tri(loc_ld_matrix)]
    diag(loc_ld_matrix)<-1
    loc_ld_matrix <- as.matrix(loc_ld_matrix)
  }

  loc_df <- loc_summary_df %>%
    filter(!is.na(Freq1)) %>%
    select(all_of(input_vars)) %>%
    mutate(Z = b/se,
           D = 2*N*Freq1*(1-Freq1),
           sig2R = D*(varVbeta*(N-1) + b^2),
           h2 = (2*b^2*Freq1*(1-Freq1))/median(sig2R, na.rm = T),
           D_std = sqrt((2*N*Freq1*(1-Freq1))/sig2R))
  
  loc_snp_idx <- which(loc_df$SNP %in% loc_summary_df$SNP)
  loc_ld_matrix <- loc_ld_matrix[loc_snp_idx, loc_snp_idx]

  gwa_h2<-trait_h2/1e6
  loc_fac<-mean(loc_df$h2, na.rm = T)/gwa_h2

  G.vec<-loc_df$Vbeta*loc_fac
  G.inv<-diag(1/G.vec)
  D.mat<-diag(loc_df$D_std)

  B <- (D.mat %*% loc_ld_matrix %*% D.mat) + G.inv
  B.inv <- MASS::ginv(B)
  b.tilde <- B.inv %*% D.mat %*% loc_df$Z

  var.b.tilde<-diag(B.inv)
  Wald.scores <- b.tilde^2/var.b.tilde

  Wald.p <- pchisq(Wald.scores, 1, lower.tail=F)

  loc_df <- loc_df %>%
    mutate(b_AnnoRE = b.tilde,
           chisq_AnnoRE = Wald.scores,
           p_AnnoRE = Wald.p)

  loc_summary_df %>%
    left_join(loc_df)
}
