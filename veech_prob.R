veech_prob <- function(asv_df, conf.level = 0.05) {
  library(lubridate)
  # Observed number of sites having both species
  occ_df <- crossprod((asv_df))
  calc_df <- NULL
  
  # Total number of sites
  N <- nrow(asv_df)
  
  print("Calculating Co-occurrences/Co-exclusions")
  
  pb <- txtProgressBar(min = 0, max = nrow(occ_df), style = 3, width = nrow(occ_df), char = "=") 
  
  # Calculate Pgl and Plt according to Veech (2013) for each pair of species
  for (r in 1:nrow(occ_df)) {
    for (c in 1:ncol(occ_df)) {
      if (c >= r) {
        next
      }
      
      # Number of sites occupied for species 1
      N1 <- occ_df[r,r]
      # Number of sites occupied for species 2
      N2 <- occ_df[c,c]
      
      # Probability that by chance the observed co-occurrence is exactly equal to the observed co-occurrence
      j.o <- occ_df[r,c]
      Pobs <- (choose(N, j.o) * choose(N-j.o, N2-j.o) * choose(N-N2, N1-j.o)) / (choose(N, N2) * choose(N, N1))
      
      # P-values testing if two species co-occur significantly more often that expected
      Pgt <- 0
      for (j.g in occ_df[r,c]+1:N) {
        Pj <- (choose(N, j.g) * choose(N-j.g, N2-j.g) * choose(N-N2, N1-j.g)) / (choose(N, N2) * choose(N, N1))
        Pgt <- Pgt + Pj
      }
      p.val_Pgt <- Pgt + Pobs 
      
      # P-values testing if two species co-occur significantly less often that expected
      Plt <- 0
      for (j.l in 0:occ_df[r,c]-1) {
        Pj <- (choose(N, j.l) * choose(N-j.l, N2-j.l) * choose(N-N2, N1-j.l)) / (choose(N, N2) * choose(N, N1))
        Plt <- Plt + Pj
      }
      p.val_Plt <- Plt + Pobs
      
      if (p.val_Pgt <= conf.level) {
        # Generate a data.frame with the results
        calc_df <- rbind(calc_df,
                         cbind.data.frame(sp1 = row.names(occ_df)[r], sp2 = row.names(occ_df)[c],
                                          p.value = p.val_Pgt,
                                          type = "Co-occurrence"))
        
      }
      if (p.val_Plt <= conf.level) {
        calc_df <- rbind(calc_df,
                         cbind.data.frame(sp1 = row.names(occ_df)[r], sp2 = row.names(occ_df)[c],
                                          p.value = p.val_Plt,
                                          type = "Co-exclusion"))
      }
    }
    
  setTxtProgressBar(pb, r)

  }
  return(calc_df)
}
