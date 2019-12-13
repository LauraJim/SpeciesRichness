# First version: December, 2018

# Code for simulating tables of incidence data for different samples sizes (number of days sampled)
# We used a Poisson mixed model with either a Gamma or a Log-normal distribution for the mean abundances
# of the species in a community, and added error species to explore the performance of the non-parametric
# estimators of species richnes, classic Chao and adjusted Chao estimates.

# Generate sampling days: 100 real species, 10 error species, and 100 replicates 
# number of days sampled
binsize=c(5,7,10,15,20,25,50,75,100,125,150,175,200,300,400,500,600,700,800,900,1000)
# number of species in the community
nsp <- 100
# number of error species
nerr <- 10
# number of replicates
nsim <- 100

# Select the parameters of the gamma/lognormal distribution
par1 <- 0.8
par2 <- 0.5
# generate mean abundances
lam <- rgamma(nsp,shape=par1,scale=par2)
# select a level of error
per.gam <- 0.000001
# generate mean abundances of error species
lam.error <- qgamma(per.gam,shape=par1,scale=par2)
# create vector with all mean abundances
lams <- c(lam,rep(lam.error,nerr))

# Loop over each binsize (sampling days) to creat eand export incidence tables
loop=0
for(k in binsize){
  input_file <- file(paste0(paste(nsp+nerr,"sp",per.gam, "err", nsim,"sim",k,"days",sep="_"),".txt"), open = "at")
  loop=loop+1
  for(j in 1:nsim){
    repmatrix=matrix(nrow=(nsp+nerr),ncol=k)
    for(i in 1:(nsp+nerr))
      # abundance matrix is filled with simulated data from Poisson distributions with given mean abundances
      repmatrix[i,]=rpois(k,lams[i])
      # convert to incidence matrix
      repmatrix <- (repmatrix>0)*1
      if(j==1){
      # use the format of input files for EstimateS
      cat("*MultipleSampleSets*\t", nsim, "\t", paste0(k, "_days"), "\n", file = input_file, sep = "")}
      cat(paste("rep", j, sep = " "), "\t*SampleSet*\t", "1\t", "0\t", "0\n", file = input_file, sep = "")
      cat(nsp+nerr, "\t", k, "\n", file = input_file, sep = "")
      write.table(repmatrix,file=input_file, sep = "\t",row.names=F, col.names=F, append = T)
  }
  flush(input_file)
  close(input_file)
}

# END