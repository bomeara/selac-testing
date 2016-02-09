rm(list=ls())
library(selac)

phy <- read.tree('yeast.tre')

set.size <- 10

gtr.params <- runif(6, 0.1, 2)
acgt.freq <- runif(4)
acgt.freq <- (acgt.freq / sum(acgt.freq))[-4]

for(g in 1:set.size)
{
  gene.length <- round(runif(1, 100, 1000))
  optim.aa <- sample(x = c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"), 
                     size = gene.length, replace = T)
  
  #They are ordered as follows: C.q.phi, alpha, beta, Ne, base.freqs for A C G, and the rates for the nucleotide model.
  model.param <- c(runif(1, 10^-10, 10^-5), runif(2, 0.01, 1), 57852168, acgt.freq, gtr.params)
  names(model.param) <- c("qphic", "alpha", "beta", "Ne", "A.freq", "C.freq", "G.freq", 
                          "GTR.alpha", "GTR.beta", "GTR.gamma", "GTR.delta", "GTR.epsilon", "GTR.eta")
  root.codon.freq <- vector(mode = "numeric", length = 64*21)
  for(i in 1:64)
  {
    x <- runif(21)
    x <- x / sum(x)
    root.codon.freq[((i-1)*21+1):(i*21)] <- x
  }
  
  nucl.data <- SelacSimulator(phy, pars = model.param, aa.optim_array = optim.aa, 
                 root.codon.frequencies = root.codon.freq, numcode=1, aa.properties=NULL, 
                 nuc.model = "GTR", k.levels=0, diploid=TRUE)
  
  sim.name <- paste("gene", g, sep = "")
  write.dna(x = nucl.data, file = paste(sim.name, "fasta", sep = "."), format = "fasta", nbcol = 1, colw = 60)
  save(list = c("model.param", "optim.aa", "root.codon.freq", "phy"), file = paste(sim.name, "Rda", sep = "."))
}
