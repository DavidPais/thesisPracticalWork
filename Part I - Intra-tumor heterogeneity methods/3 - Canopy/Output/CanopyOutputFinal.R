#' Created on Sun Jun 19 23:45:02 2021

#' @author: david

#  Canopy heterogeneity estimation algorithm output processing


# Instructions: If first time running, run the code until storing the canopy executions object, so that next time the data is needed
#               it is available for the plots output processing that follows, not needing to run all the executions again.
#               The next time this allows the loading of the executions object, starting with the samples' heterogeneity index calculation 
#               based on the clones frequency.


#------------------------------------------------------------------------------------------------------------------------------------------



#####  Methods  #####  


# Calculate samples' Shannon's Index
shannonIndex <- function(clonal_freq) { # for each clonal frequency value 
  if(clonal_freq == 0) {
    return(0)
  } else {
    return(clonal_freq * log(clonal_freq))
  }
}



#------------------------------------------------------------------------------------------------------------------------------------------



#######  Load the output file from the MCMC sampling  #######


# On linux
# setwd("~/Desktop/Dados/CNVS/canopyImages")

# On windows
setwd("C:/Users/david/Desktop/Dados/CNVS/canopyImages")


# R Canopy package
library(Canopy)


# The different Canopy executions (with the 4 mice separately or with all mice together)
miceExecs <- c("M49", "M55", "M61", "M62", "AllMice")
# The projects name list
pnames <- list()


for(i in 1 : as.numeric(length(miceExecs))) {

  if(miceExecs[i] == "AllMice") {
    pnames[[i]] <- "MDA231ctrlsunit"
  }


  if(miceExecs[i] == "M49" |  miceExecs[i] == "M55") {
    pnames[[i]] <- paste("MDA231Ctrl", miceExecs[i], sep = "")
  }


  if(miceExecs[i] == "M61" |  miceExecs[i] == "M62") {
    pnames[[i]] <- paste("MDA231Sunit", miceExecs[i], sep = "")
  }


}



# To increase the memory used by RStudio in order to be able to run the executions
if(.Platform$OS.type == "windows") withAutoprint({
  memory.size()
  memory.size(TRUE)
  memory.limit()
})
memory.limit(size=80000)




# using m variable because i was used in the other program and is load with a value different from the current being used here
for (m in 1 : as.numeric(length(pnames))) {
  
  projectname <- pnames[[m]]
  
  # Load the canopy.sample results for each execution configuration
  load(paste(projectname, '_postmcmc_image.rda', sep=''))
  
  
  if(as.numeric(m) == 5) { # since the input variable is not named canopyInputAllMice, assign the canopyInput for all mice content to this new variable name
    canopyInputAllMice <- canopyInput
  } else {
    canopyInputAllMice <- "" # so that the rm of canopyInputAllMice does not give any warning
  }

  
  # The matrix of the intersection between the CNA regions and CNAs will be needed
  assign(paste("C", miceExecs[m], sep = ""), get(paste("canopyInput", miceExecs[m], sep = ""),  envir = .GlobalEnv)[[3]])
  
  
  # Delete environment variables to save space
  rm(cMatrix, cMatrixM49, cMatrixM55, cMatrixM61, cMatrixM62, majorCopyMatrix, majorCopyMatrixM49, majorCopyMatrixM55, majorCopyMatrixM61, 
     majorCopyMatrixM62, majorEpsilon, majorEpsilonM49, majorEpsilonM55, majorEpsilonM61, majorEpsilonM62, minorCopyMatrix, minorCopyMatrixM49,
     minorCopyMatrixM55, minorCopyMatrixM61, minorCopyMatrixM62, minorEpsilon, minorEpsilonM49, minorEpsilonM55, minorEpsilonM61, minorEpsilonM62,
     rMatrix, rMatrixM49, rMatrixM55, rMatrixM61, rMatrixM62, xMatrix, xMatrixM49, xMatrixM55, xMatrixM61, xMatrixM62, yMatrix, yMatrixM49, yMatrixM55,
     yMatrixM61, yMatrixM62, allCopyNumber, allSamples, allSamplesCNA, cnasInfo, copyNumberInfo, genomeInfo, miceCnasInfo, miceCopyNumberInfo, miceGenomeInfo,
     C, epsilonm, epsilonM, R, Wm, WM, X, Y, mice, canopyInput, canopyInputM49, canopyInputM55, canopyInputM61, canopyInputM62, canopyInputAllMice,
     regions, get_identifier, read_counts, g, groups, h, i, j, k, projectname, projectnameM49, projectnameM55, projectnameM61, projectnameM62)
  
  
  # store the sampchain variable with a different identifier for each execution
  sc <- get("sampchain", envir = .GlobalEnv)
  assign(paste("sampchain", miceExecs[m], sep = ""), sc)
  
  rm(sampchain, sc)
  gc()

}



#------------------------------------------------------------------------------------------------------------------------------------------



# In the canopy package manual (https://cran.r-project.org/web/packages/Canopy/Canopy.pdf), the sections/functions to look at are:
# canopy.sample, canopy.BIC, canopy.post, canopy.output, canopy.plottree 


# On linux
# setwd("~/Desktop/Dados/CNVS/canopyOutputFiles")

# On windows
setwd("C:/Users/david/Desktop/Dados/CNVS/canopyOutputFiles")



#######   BIC to determine number of subclones  #######


K <- get("K", envir = .GlobalEnv) # number of subclones
numchain <- get("numchain", envir = .GlobalEnv) # number of chains with random initiations


# If there is error in the bic and canopy.post step below, make sure burnin and thinning parameters are wisely selected so that there are posterior trees left
# Typically a number of the initial states are attributed to burnin and removed, whilst the remainder of the chain is thinned if compression is also required

# Use burnin = 100 and thin = 5 for M49 and M55 (chains converge faster by looking at the likelihood pdfs), and burnin = 200 and thin = 5 for M61, M62, and all mice together

burnin <- c(100, 100, 200, 200, 200) # burnin of MCMC chains
thin <- 5 # MCMC chains thinning  
              

for (m in 1 : as.numeric(length(miceExecs))) {

  sc <- get(paste("sampchain", miceExecs[m], sep = "")) # list of sampled trees in subtree space with different number of subclone
 
  
  # To get BIC as a model selection criterion from MCMC sampling results
  # Returns: BIC values (vector) for model selection with plot generated (pdf format).
  bic <- canopy.BIC(sampchain = sc, projectname = pnames[[m]], K = K,
                    numchain = numchain, burnin = burnin[[m]], thin = thin, pdf = TRUE)
  
  # optK (optimal number of subclones) including the normal subclone
  optK <- K[which.max(bic)]
  
  # store each of the executions bic and optK
  assign(paste("bic", miceExecs[m], sep = ""), bic)
  assign(paste("optK", miceExecs[m], sep = ""), optK)
  
  
  rm(bic, optK, sc)
  gc()
  
}




#------------------------------------------------------------------------------------------------------------------------------------------




#######  Posterior tree evaluation  #######


for (m in 1 : as.numeric(length(miceExecs))) {

# get needed parameters for each execution  
  
sc <- get(paste("sampchain", miceExecs[m], sep = ""),  envir = .GlobalEnv)

optK <- get(paste("optK", miceExecs[m], sep = ""),  envir = .GlobalEnv)
  
C <- get(paste("C", miceExecs[m], sep = ""),  envir = .GlobalEnv)


# Posterior evaluation of sampled trees by MCMC

# Returns: samptreethin - list of sampled posterior trees (like sampchain, but after burnin and thinning)
#          samptreethin.lik - vector of likelihood of sampled posterior trees
#          config - vector of configuration of sampled posterior trees (integer values)
#          config.summary - summary of configurations of sampled posterior trees
post <- canopy.post(sampchain = sc, projectname = pnames[[m]], K = K,
                    numchain = numchain, burnin = burnin[[m]], thin = thin, 
                    optK = optK, C = C, post.config.cutoff = 0.05)

rm(sc, optK, C)
gc()


assign(paste("post", miceExecs[m], sep = ""), post)


# first column: tree configuration
# second column: posterior configuration probability in the entire tree space
# third column: posterior configuration likelihood in the subtree space
config.summary <- post[[4]] 

# For each iteration in the MCMC, we calculate the likelihood. For trees with the same configurations, 
# we calculate the mean of their likelihood and output Mean_post_lik. Their configurations are the same,
# but their VAFs can be different resulting in slightly different likelihoods. 
print(config.summary)


# Delete environment variables just to save space
# list argument to remove dynamically created variables
rm(list = paste("sampchain", miceExecs[m], sep = ""), post, config.summary)
gc()

}




#------------------------------------------------------------------------------------------------------------------------------------------




#######  Tree output and plot  #######


# One can then use Canopy to output and plot the most likely tree (i.e., tree with the highest posterior likelihood).
# Mutations, clonal frequencies, tree topology, and other information of the tree are obtained from 
# the posterior distributions of subtree space with trees having the same configuration


for (m in 1 : as.numeric(length(miceExecs))) {

  post <- get(paste("post", miceExecs[m], sep = ""),  envir = .GlobalEnv)
  
  # not saving the config.summary for each mouse since it can be obtained from the post variable
  config.summary <- post[[4]]
  
  
  # choose the configuration with the highest posterior likelihood
  # configuration of sub-tree space to be output
  assign(paste("config.i", miceExecs[m], sep = ""), config.summary[which.max(config.summary[,3]),1])
  config.i <- get(paste("config.i", miceExecs[m], sep = ""), envir = .GlobalEnv)
  
  cat('Configuration', config.i, 'has the highest posterior likelihood.\n')


  C <- get(paste("C", miceExecs[m], sep = ""),  envir = .GlobalEnv)

  
  # To generate a posterior tree from the sub-tree space of trees with the same configurations.
  # Returns: posterior tree from the sub-tree space of trees with the same configurations.
  output.tree <- canopy.output(post, config.i, C)

  assign(paste("output.tree", miceExecs[m], sep = ""), output.tree)


  # Delete environment variables just to save space
  rm(list = paste("post", miceExecs[m], sep = ""), C, post, config.summary, config.i)
  gc()


  
  # Note: A separate txt file can be generated (with txt=TRUE and txt.name= txt name) 
  # if the figure legend of mutational profiles (texts below the phylogeny tree) in the plot 
  # is too long to be fitted entirely
  
  
  # clonal frequencies (i.e., the composition of the clones for the tumor samples, which is the P matrix in our notation)
  # output.tree[["P"]]
  pdf.name <- paste(pnames[[m]], '_config_highest_likelihood.pdf', sep='')
  
  # output.tree[["clonalmut"]]
  txt.name <- paste(pnames[[m]], '_config_mutational_profiles.txt', sep ='')
  
  
  # Each sampled tree is modeled as a list by Canopy.
  # The most likely tree is obtained from the posterior distribution in the tree space from the MCMC sampling
  
  # To plot Canopy's reconstructed phylogeny (the most likely tree in the posterior tree space inferred by Canopy). Major plotting function of Canopy
  # Returns: Plot of tree structure, clonal frequency and mutation legends (pdf format)
  canopy.plottree(output.tree, pdf = TRUE, txt = TRUE, pdf.name = pdf.name, txt.name = txt.name)

  
  rm(list = paste("C", miceExecs[m], sep = ""), output.tree, pdf.name, txt.name)
  gc()  
  
}


# save the Canopy executions output.tree matrices so that every time we start the program we do not have to run all the files first again
save(output.treeM49, output.treeM55, output.treeM61, output.treeM62, output.treeAllMice,
     burnin, thin, K, miceExecs, numchain, optKM49, optKM55, optKM61, optKM62, shannonIndex, file = "canopyExecutions.RData")




#------------------------------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------------------------------




#####  Samples' heterogeneity index value based on the clonal frequencies  #####


# When loading the Canopy executions object, start the program here


# On linux
# setwd("~/Desktop/Dados/CNVS/canopyOutputFiles")

# On windows
setwd("C:/Users/david/Desktop/Dados/CNVS/canopyOutputFiles")


load("canopyExecutions.RData")


miceExecs <- get("miceExecs", envir = .GlobalEnv) # number of subclones


for (m in 1 : as.numeric(length(miceExecs))) {

  output.tree <- get(paste("output.tree", miceExecs[m], sep = ""), envir = .GlobalEnv)

  
  # the composition of the clones for the tumor samples, which is the P matrix in our notation
  c_f <- output.tree[["P"]] 
  
  
  # To put it like the clonal frequency matrix of the Canopy output
  
  # transpose the matrix
  c_f <- t(c_f)
  
  # change column names
  i <- 1
  
  for(col_name in colnames(c_f)) {
    
    clone_num <- as.numeric(substr(col_name, nchar(col_name), nchar(col_name)))
    clone_num <- clone_num - 1
    
    colnames(c_f)[i] <- paste("Clone", toString(clone_num), sep = "")
    
    i <- i + 1
  }
  
  colnames(c_f)[1] <- "Normal"
  
  
  
  # Calculate the heterogeneity index for each sample based on the clones frequency 
  
  
  # lists of clonal frequencies values for each sample
  list_cf <- split(c_f, seq(nrow(c_f)))
  
  
  # create the dataframe with the heterogeneity index values based on the clonal frequencies for each sample
  HI <- c()
  
  for (l in list_cf){
    HI <- append(HI, -sum(sapply(l, shannonIndex)))
  }
  
  HI <- data.frame(HI)
  row.names(HI) = row.names(c_f)


  assign(paste("heterogeneityIndex", miceExecs[m], sep = ""), HI)


  rm(c_f, list_cf, HI)
  gc()
}




#------------------------------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------------------------------




#####  Plot mice samples' heterogeneity index value based on the clonal frequencies for all mice together and comparing different mice  #####


# install required packages
# install.packages('dplyr')

# import libraries

# for bind_rows()
library(dplyr)

# for plotting
library(ggplot2)


# Dataframe plot input processing
HI_plots <- list(heterogeneityIndexAllMice, list(heterogeneityIndexM49, heterogeneityIndexM55, heterogeneityIndexM61, heterogeneityIndexM62))
  

for (m in 1 : as.numeric(length(HI_plots))) {

  sample_HI <- HI_plots[[m]]
  
  # bind the rows
  sample_HI <- bind_rows(sample_HI)
  
  samples <- row.names(sample_HI)
  
  # x represents each sample
  mouse <- sapply(samples, function(x) substr(x, 1, 3))
  group <- sapply(mouse, function(x) if(x == "M49" | x == "M55") {"Ctrl"} else {"Sunit"} )
  
  sample_HI$sample <- samples
  sample_HI$mouse <- mouse
  sample_HI$group <- group
  
  
  
  if(as.numeric(m) == 1) {
    subt = "Using all mice together in a file execution"
    figure_name = "Comparison of the control and treated groups samples heterogeneity index variation using all mice together.png"
  } else {
    subt = "Using different mice files executions"
    figure_name = "Comparison of the control and treated groups samples heterogeneity index variation using different mice.png"
  }
  
  
  
  # Plotting heterogeneity index 
  
  
  print( # When in a for loop, you have to explicitly print your resulting ggplot object
  
  sample_HI %>%
    
    ggplot(aes(x = mouse, y = HI, color = group)) +
    
    # https://sphweb.bumc.bu.edu/otlt/MPH-Modules/PH717-QuantCore/PH717-Module6-RandomError/PH717-Module6-RandomError11.html
    
    # Plot with mean value estimator and confidence intervals of 95%
    # CI formula: x_mean +- t-value (since n < 30) * (s_deviation (or s_error) / sqrt (sample_size))
    
    # t-value for 95% = 3.182
    # sample size is 4 (because for each mouse there are 4 samples) 
    
    # Using the t.test function
    stat_summary(fun = function(mouse_points) t.test(mouse_points)[5][[1]][[1]], # mean value 
                 fun.min = function(mouse_points) t.test(mouse_points)[4][[1]][1], # CI lower bound
                 fun.max = function(mouse_points) t.test(mouse_points)[4][[1]][2], # CI higher bound
                 geom = "pointrange") +
    
    labs(
      title = "Comparison of the control and treated groups samples' heterogeneity index variation",
      subtitle = subt,
      x = "Mouse ID",
      y = "Sample TH Index"
    ) + 
    
    # colors in R: http://sape.inf.usi.ch/quick-reference/ggplot2/colour
    scale_colour_manual(values = c("dodgerblue3", "darkorange"), name = "Group", labels = c("Control", "Sunitinib")) +
    
    theme_classic() + 
    
    # all ggplot methods for theme - https://ggplot2.tidyverse.org/reference/theme.html
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, margin=margin(2,0,17,0), size = 11.5),  
          axis.title.x = element_text(margin=margin(7,0,0,0), size = 11),
          axis.title.y = element_text(margin=margin(0,7,0,0), size = 11),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 11.5),
          legend.key.size = unit(1, 'cm'),
          legend.key.width = unit(0.3, 'cm'),
          legend.spacing.y = unit(0.25, 'cm'),
          legend.box.spacing = unit(0.23, 'cm')
    ) 
  
  )
  
  # Save image (will save into the current directory) and keep seeing the plot in R
  # width and height of the plot in Rstudio when priting without measures
  ggsave(figure_name, width = 8.99, height = 6.47, units = "in")


}




