# Code for calculating the adjusted Chao estimator for species richness (according to Chiu and Chao, 2016)

# Function to calculate Chao adjusted estimator from simulated incidence data
adjChao <- function(file1,nrep){
  # read the data
  dat <- as.matrix(read.delim(file1,header = F))
  # number of observed species
  nsp <- as.numeric(dat[3,1])
  # number of days sampled
  nd <- as.numeric(dat[3,2])
  # define a matrix that will contain the values of the Chao estimators
  # and a binary value indicating if the adjusted, or the classic, estimator was calculated
  chaos <- matrix(NA,nrow=nrep,ncol=2)
  #
  for(i in 1:nrep){
    # first row of the i-th matrix
    ini <- 4 + (i-1)*(nsp+2)
    # last row of the i-th matrix
    fin <- 1 + i*(nsp+2)
    # variable to store the i-th matrix
    subM <- matrix(NA,nrow=nsp,ncol=nd)
    # index used to fill the i-th matrix
    k <- 1
    # index used to indicate if the adjusted estimator is calculated
    adj <- 0 # 0 means the classic estimator is calculated
    for(j in ini:fin){
      # make sure that the new row to be added is numeric
      subM[k,] <- as.numeric(dat[j,])
      # go to the next row of the sub-matrix
      k <- k+1
    }
    # once the sub-matrix is defined, let's calculate the adjusted Chao2 estimator
    # first, we need to calculate the sample frequencies and store them in a matrix
    sp.freqs <- rowSums(subM)
    # calculate the number of uniques, duplicates, etc.
    smpl.freqs <- tabulate(sp.freqs)
    # actually, we only need Q1, Q2, Q3 and 
    Q1 <- smpl.freqs[1]
    Q2 <- smpl.freqs[2]
    Q3 <- smpl.freqs[3]
    Q4 <- smpl.freqs[4]
    # with the frequencies vector, we calculate the number of observed species
    Sobs <- sum(smpl.freqs)
    if(Q3>0 && Q4>0){
      # second, calculate the estimator of the true number of uniques in terms of Q2,Q3,Q4
      Q1.hat <- 2*Q2 * ((2*Q2/(3*Q3)) - (Q3/(4*Q4)))
      adj <- 1
      # third, calculate the adjusted Chao estimator
      if(Q2>0){
        ChaoEst <- Sobs - Q1 + Q1.hat + ((nd-1)*(Q1.hat^2))/(2*Q2*nd)
      } else {
        ChaoEst <- Sobs - Q1 + Q1.hat + (Q1.hat*(Q1.hat-1))/2
      }
    } else {
      # when Q3 and/or Q4 are zeros, calculate the classic estimator
      if(Q2>0){
        ChaoEst <- Sobs + ((nd-1)*(Q1^2))/(2*Q2*nd)
      } else {
        ChaoEst <- Sobs + (Q1*(Q1-1)*(nd-1))/(2*nd)
      }
    }
    chaos[i,] <- c(adj,ChaoEst)
  }
  return(chaos)
}

# MAIN ------------
#library(breakaway)
home <- ""
setwd(home)
for (file in list.files(path = home, pattern = "*.txt")){
  # number of replicates per case
  nrep <- 100
  # Calculate the Chao estimators
  results <- adjChao(file,nrep)
  write.table(results, file = paste0(paste0("./output/",file,"_output"), ".txt"))
}

### Change files name----------------
outputs <- ""
setwd(outputs)
for (txtfile in list.files(path = outputs, pattern = "*.txt")){
  file.rename(from=txtfile,to=sub(pattern=".txt_output",replacement="_output",txtfile))
}

### First version created by Laura Jimenez on July 2018