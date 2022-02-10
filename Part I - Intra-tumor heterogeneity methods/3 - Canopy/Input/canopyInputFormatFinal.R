#' Created on Mon Apr 19 22:50:07 2021

#' @author: david


# Install Canopy Heterogeneity Estimation Algorithm 

# install.packages(c("ape", "fields", "pheatmap", "scatterplot3d", "devtools"))
# devtools::install_github("yuchaojiang/Canopy/package")


#------------------------------------------------------------------------------------------------------------------------------------------

#####  Methods  #####  

read_counts <- function(alleleInfo) { # for each line of tumorInfo column 
  refCounts <- substr(alleleInfo, unlist(gregexpr(":", alleleInfo))[1] + 1, unlist(gregexpr(",", alleleInfo))[1] - 1)
  mutCounts <- substr(alleleInfo, unlist(gregexpr(",", alleleInfo))[1] + 1, unlist(gregexpr(":", alleleInfo))[2] - 1)
  counts <- c(refCounts, mutCounts)
  return (counts)
}


get_identifier <- function(identifier_row) { # extract the identifier from a SNV/CNA, not including the sample info
  id <- substr(identifier_row, 1, as.numeric(unlist(gregexpr(':', identifier_row))[3]) - 1)
  return(id)
  
}

#------------------------------------------------------------------------------------------------------------------------------------------



##### Load somatic mutations input (SNVs) #####


# On linux
# setwd("~/Desktop/Dados/variantes_somaticas/excel_format")


# On windows
setwd("C:/Users/david/Desktop/Dados/variantes_somaticas/excel_format")



# for all mice together
genomeInfo <- list()


# for the 4 mice separately
miceGenomeInfo <- list()


for (i in 1 : 4) {
  miceGenomeInfo[[i]] <- list()
}


auxMiceGenome <- list()


library(data.table)


mice <- c("M49", "M55", "M61", "M62")

groups <- c("Ctrl", "Sunit")

regions <- c("R1", "R2", "R3", "R4")



i <- 1

for(file in sort(list.files(getwd()))) {
  
  # On linux
  # data <- read.csv(file, header = TRUE, sep = ";")
  
  # On windows
  data <- read.csv(file, fileEncoding="UTF-8-BOM", header = TRUE, sep = ";") 
  
  data <- data[data$FILTER == "PASS", ]
  
  # sex chromosomes are not enabled by canopy
  data <- data[data$CHROM != 'chrX', ]
  
  data <- data[ , c("CHROM", "POS", "TUMOR", "ALT")]
  
  genomeInfo[[i]] <- data
  
  
  tumorInfo <- data$TUMOR
  
  readCounts <- sapply(tumorInfo, read_counts)
  
  
  genomeInfo[[i]] <- cbind(genomeInfo[[i]], "REF COUNTS" = readCounts[c(1), ])
  genomeInfo[[i]] <- cbind(genomeInfo[[i]], "ALT COUNTS" = readCounts[c(2), ])
  
  genomeInfo[[i]] <- cbind(genomeInfo[[i]], "MAJOR COPY" = -1)
  genomeInfo[[i]] <- cbind(genomeInfo[[i]], "MINOR COPY" = -1)
  
  # remove "TUMOR" column
  genomeInfo[[i]] <- genomeInfo[[i]] [  , -3]
  
  
  # insert mutation identifier
  mutation_id <- paste(paste(
    toupper(genomeInfo[[i]]$CHROM), genomeInfo[[i]]$POS, sep = ":"), genomeInfo[[i]]$ALT, sep = ":")
  
  genomeInfo[[i]] <- cbind(genomeInfo[[i]], "MUTATION ID" = mutation_id)
  
  
  # insert sample identifier
  
  # because here indexes start in 1, and remainder of integer division is %%
  genomeInfo[[i]] <- cbind(genomeInfo[[i]], "SAMPLE ID" = paste(paste(mice[as.numeric(as.integer((i-1)/4)) + 1], 
                                                                      groups[as.numeric(as.integer(i/9)) + 1], sep = ""), regions[((i - 1) %% 4) + 1], sep = ""))
  
  
  # remove chrom, position and alt columns
  genomeInfo[[i]] <- genomeInfo[[i]] [ , c(length(genomeInfo[[i]]) - 1, length(genomeInfo[[i]]), 2, 4:(length(genomeInfo[[i]]) - 2))]
  
  
  # reset index numbers
  row.names(genomeInfo[[i]]) <- NULL
  
  
  #------------------------------------------------------------------------------------------------------------------------------------------
  
  
  ##### For each of the 4 mice #####
  
  
  auxMiceGenome <- c(auxMiceGenome, genomeInfo[i])
  
  
  # Collect the 4 samples of each mouse
  
  if(i %% 4 == 0) {
    
    miceGenomeInfo[[i / 4]] <- rbindlist(auxMiceGenome) 
    
    row.names(miceGenomeInfo[[i / 4]]) <- NULL
    
    auxMiceGenome <- list()
    
  }
  
  i <- i + 1
  
}



#------------------------------------------------------------------------------------------------------------------------------------------



#####  Load copy number alterations input (CNAs)  #####


# on linux
# setwd("~/Desktop/Dados/CNVS/excel_format_germline")


# on windows
setwd("C:/Users/david/Desktop/Dados/CNVS/excel_format_germline")



# for all mice together


# with the information of the copy number of all the regions of all samples (including CNA regions and Non-CNA regions)
copyNumberInfo <- list()

# only with the CNAs (cn copy number != 2) for each sample
cnasInfo <- list()



# for the separate 4 mice
miceCopyNumberInfo <- list()

miceCnasInfo <- list() 

auxMiceCopyNumber <- list()

auxMiceCnas <- list()


for (i in 1 : 4) {
  miceCopyNumberInfo[[i]] <- list()
  miceCnasInfo[[i]] <- list()
  
}


i <- 1


for(file in sort(list.files(getwd()))) {
  
  # On linux
  # data <- read.csv(file, header = TRUE, sep = "\t")
  
  # On windows
  data <- read.csv(file,  fileEncoding="UTF-8-BOM", header = TRUE, sep = ";") 
  
  # Removing last row (Y chromosomes, since it is breast cancer)
  data <- data[ -c(nrow(data)), ]
  
  # sex chromosomes are not enable by canopy
  data <- data[data$chromosome != 'chrX', ]
  
  # eliminate those cnas that might be noise, with a small difference between start and end position
  data <- data[ ( as.numeric(data$end) - as.numeric(data$start) ) > 1048576, ]
  
  # complete.cases() returns the rows that have all columns ! NA, removing the rows that have some column with NA
  data <- data[complete.cases(data), ]
  
  names(data)[names(data) == "chromosome"] <- "CHROM"
  names(data)[names(data) == "start"] <- "START POS"
  names(data)[names(data) == "end"] <- "END POS"
  names(data)[names(data) == "cn"] <- "COPY NUMBER"
  names(data)[names(data) == "cn1"] <- "MAJOR COPY"
  names(data)[names(data) == "cn2"] <- "MINOR COPY"
  names(data)[names(data) == "ci_lo"] <- "CI LOWER"
  names(data)[names(data) == "ci_hi"] <- "CI UPPER"
  names(data)[names(data) == "baf"] <- "BAF"
  names(data)[names(data) == "log2"] <- "LOG2"
  
  copyNumberInfo[[i]] <- data
  
  copy_number_id = paste(paste(paste(
    toupper(copyNumberInfo[[i]]$CHROM), copyNumberInfo[[i]]$'START POS', sep = ":"), copyNumberInfo[[i]]$'END POS', sep = ":"),
    paste(paste(paste(mice[as.numeric(as.integer((i-1)/4)) + 1], groups[as.numeric(as.integer(i/9)) + 1], sep = ""),
                regions[((i - 1) %% 4) + 1], sep = "")), sep = ":")
  
  
  copyNumberInfo[[i]] <- cbind(copyNumberInfo[[i]], "COPY ALTERATION ID" = copy_number_id)
  
  
  # because here indexes start in 1, and remainder of integer division is %%
  copyNumberInfo[[i]] <- cbind(copyNumberInfo[[i]], "SAMPLE ID" = paste(paste(mice[as.numeric(as.integer((i-1)/4)) + 1],
                                                                              groups[as.numeric(as.integer(i/9)) + 1], sep = ""), regions[((i - 1) %% 4) + 1], sep = ""))
  
  
  copyNumberInfo[[i]] <- copyNumberInfo[[i]][ , c("COPY ALTERATION ID", "SAMPLE ID", "CHROM", "START POS", "END POS", "COPY NUMBER", 
                                                  "MAJOR COPY", "MINOR COPY", "CI LOWER", "CI UPPER", "BAF", "LOG2")]
  
  
  # removing CNAs with major copy number = 2 and minor copy number = 0
  cnas_to_remove <- as.numeric(row.names(copyNumberInfo[[i]] [copyNumberInfo[[i]]$`MAJOR COPY` == 2 & copyNumberInfo[[i]]$`MINOR COPY` == 0, ]))
  
  # remove CNAs that may have been incorrectly identified
  copyNumberInfo[[i]] <- copyNumberInfo[[i]] [! (as.numeric(row.names(copyNumberInfo[[i]])) %in% cnas_to_remove), ]

  
  row.names(copyNumberInfo[[i]]) <- NULL
 
  
  
  cnas_to_remove <- as.numeric(row.names(copyNumberInfo[[i]] [copyNumberInfo[[i]]$`MAJOR COPY` == 1 & copyNumberInfo[[i]]$`MINOR COPY` == 1, ]))
  
  # only keep CNAS (with copy number != 2, major copy != 1 and minor copy != 1)
  cnasInfo[[i]] <- copyNumberInfo[[i]] [! (as.numeric(row.names(copyNumberInfo[[i]])) %in% cnas_to_remove), ]

  
  row.names(cnasInfo[[i]]) <- NULL



  #------------------------------------------------------------------------------------------------------------------------------------------



  # For each of the 4 mice separately

  auxMiceCopyNumber <- c(auxMiceCopyNumber, copyNumberInfo[i])

  auxMiceCnas <- c(auxMiceCnas, cnasInfo[i])


  if(i %% 4 == 0) {

    miceCopyNumberInfo[[i / 4]] <- rbindlist(auxMiceCopyNumber)

    row.names(miceCopyNumberInfo[[i / 4]]) <- NULL


    miceCnasInfo[[i / 4]] <- rbindlist(auxMiceCnas)

    row.names(miceCnasInfo[[i / 4]]) <- NULL


    auxMiceCopyNumber <- list()

    auxMiceCnas <- list()

  }
 
  i <- i + 1
   
}


# if the variable is not needed again, use this
rm(auxMiceGenome, auxMiceCnas, auxMiceCopyNumber, data, readCounts, copy_number_id, file, mutation_id, tumorInfo, cnas_to_remove)




#------------------------------------------------------------------------------------------------------------------------------------------



#####  Insert SNVs in the samples that they are not part of  #####


# For all mice

allSamples <- rbindlist(genomeInfo)


# install.packages(c("sqldf"))

library(sqldf)


samples_count <- sqldf("SELECT `MUTATION ID`, COUNT(*) AS 'SAMPLES COUNT' from allSamples 
                        GROUP BY `MUTATION ID`
                        ORDER BY COUNT(*) DESC")

# select the entire row , with all columns
to_introduce <- samples_count[samples_count$`SAMPLES COUNT` != 16, ] 

row.names(to_introduce) <- NULL


samples_id <- unique(allSamples$`SAMPLE ID`)



new_mut_samples <- list()

h <- 1

for (i in 1:nrow(to_introduce)) {
  
  in_samples <- allSamples[allSamples$`MUTATION ID` == to_introduce[i, 'MUTATION ID']]$`SAMPLE ID`
  
  samples_to_include <- setdiff(samples_id, in_samples)
  
  for (j in 1:length(samples_to_include)) {
    
    # ref counts can be approximated with the average sequencing depth
    new_mut_samples[[h]] <- c('MUTATION ID' = to_introduce[i, 'MUTATION ID'], 'SAMPLE ID' = samples_to_include[j], 
                              'POS' = -1, 'REF COUNTS' = 20, 'ALT COUNTS' = 0, 'MAJOR COPY' = -1, 'MINOR COPY' = -1)
    
    h <- h + 1 
    
  }
  
  
}


df <- as.data.frame(do.call(rbind, new_mut_samples))


# All mutations in all samples
allSamples <- rbind(allSamples, df)

allSamples$'VAF' <- as.numeric(allSamples$`ALT COUNTS`) / (as.numeric(allSamples$`REF COUNTS`) + as.numeric(allSamples$`ALT COUNTS`))



#------------------------------------------------------------------------------------------------------------------------------------------



# For each mice separately


for(i in 1 : as.numeric(length(miceGenomeInfo))) {
  
  currentMice <- miceGenomeInfo[[i]]
  
  miceSamplesCount <- sqldf("SELECT `MUTATION ID`, COUNT(*) AS 'SAMPLES COUNT' from currentMice 
                             GROUP BY `MUTATION ID`
                             ORDER BY COUNT(*) DESC")
  
  miceToIntroduce <- miceSamplesCount[miceSamplesCount$`SAMPLES COUNT` != 4, ] 
  
  row.names(miceToIntroduce) <- NULL
  
  
  samples_id <- unique(currentMice$`SAMPLE ID`)
  
  
  miceMutSamples <- list()
  
  h <- 1
  
  
  
  for (j in 1 : as.numeric(nrow(miceToIntroduce))) {
    
    in_samples <- currentMice[currentMice$`MUTATION ID` == miceToIntroduce[j, 'MUTATION ID']]$`SAMPLE ID`
    
    samples_to_include <- setdiff(samples_id, in_samples)
    
    
    for (g in 1 : as.numeric(length(samples_to_include))) {
      
      # ref counts can be approximated with the average sequencing depth
      miceMutSamples[[h]] <- c('MUTATION ID' =  miceToIntroduce[j, 'MUTATION ID'], 'SAMPLE ID' = samples_to_include[g],
                               'POS' = -1, 'REF COUNTS' = 20, 'ALT COUNTS' = 0, 'MAJOR COPY' = -1, 'MINOR COPY' = -1)
      
      h <- h + 1
      
    }
    
    
  }
  
  df <- as.data.frame(do.call(rbind, miceMutSamples))
  
  
  # All mutations in all samples for each mouse 
  allSamplesMice <- rbind(currentMice, df)
  
  allSamplesMice$'VAF' <- as.numeric(allSamplesMice$`ALT COUNTS`) / (as.numeric(allSamplesMice$`REF COUNTS`) + as.numeric(allSamplesMice$`ALT COUNTS`))
  
  
  miceGenomeInfo[[i]] <- allSamplesMice
  
}



rm(in_samples, samples_to_include, currentMice, df, miceMutSamples, miceSamplesCount, miceToIntroduce, new_mut_samples, samples_count,
   to_introduce, allSamplesMice)



#------------------------------------------------------------------------------------------------------------------------------------------



#####  Create R and X matrices  #####


# order mutations according their chrm

# like x is the dataframe column allSamplesCNA$`MUTATION ID`
allSamples$'CHRM NUMBER' <- sapply(allSamples$`MUTATION ID`, function(x) as.numeric(substr(x,
                                                                                           unlist(gregexpr("R", x))[1] + 1,  unlist(gregexpr(":", x))[1] - 1)))

allSamples <- allSamples[order(as.numeric(allSamples$`CHRM NUMBER`)) , ]



mut_identifiers <- unique(allSamples$`MUTATION ID`)

samples_id <- unique(allSamples$`SAMPLE ID`)



rMatrix <- matrix(0, nrow = length(mut_identifiers), ncol = length(samples_id), byrow = TRUE, 
                  dimnames = list(paste("SNA", mut_identifiers, sep = ":"), samples_id))

xMatrix <- matrix(0, nrow = length(mut_identifiers), ncol = length(samples_id), byrow = TRUE, 
                  dimnames = list(paste("SNA", mut_identifiers, sep = ":"), samples_id))




for (i in 1 : length(mut_identifiers)) {
  
  
  for (j in 1 : length(samples_id)){
    
    mut_read_counts <- allSamples[(allSamples$`MUTATION ID` == mut_identifiers[[i]]) & 
                                    (allSamples$`SAMPLE ID` == samples_id[[j]])]$`ALT COUNTS`
    
    ref_read_counts <- allSamples[(allSamples$`MUTATION ID` == mut_identifiers[[i]]) & 
                                    (allSamples$`SAMPLE ID` == samples_id[[j]])]$`REF COUNTS`
    
    # mut counts
    rMatrix[i, j] <- as.numeric(mut_read_counts)
    
    # total read counts
    xMatrix[i, j] <- as.numeric(mut_read_counts) + as.numeric(ref_read_counts)
    
    
    j <- j + 1
    
  }
  
  i <- i + 1
  
}


# Changes matrices from double to integer

mode(rMatrix) <- "integer"

mode(xMatrix) <- "integer"


#------------------------------------------------------------------------------------------------------------------------------------------



# Create R and X matrices for each mice separately


for(i in 1 : as.numeric(length(miceGenomeInfo))) {
  
  currentMice <- miceGenomeInfo[[i]]
  
  
  # Order rows according chrm number
  
  # like x is the dataframe column allSamplesCNA$`MUTATION ID`
  currentMice$'CHRM NUMBER' <- sapply(currentMice$`MUTATION ID`, function(x) as.numeric(substr(x,
                                                                                               unlist(gregexpr("R", x))[1] + 1,  unlist(gregexpr(":", x))[1] - 1)))
  
  
  currentMice <- currentMice[order(as.numeric(currentMice$`CHRM NUMBER`)) , ]
  
  
  miceGenomeInfo[[i]] <- currentMice
  
  
  
  mut_identifiers <- unique(currentMice$`MUTATION ID`)
  
  samples_id <- unique(currentMice$`SAMPLE ID`)
  
  
  
  assign(paste("rMatrix", mice[i], sep = ""), matrix(0, nrow = length(mut_identifiers), ncol = length(samples_id), byrow = TRUE, 
                                                     dimnames = list(paste("SNA", mut_identifiers, sep = ":"), samples_id)))
  
  
  assign(paste("xMatrix", mice[i], sep = ""), matrix(0, nrow = length(mut_identifiers), ncol = length(samples_id), byrow = TRUE, 
                                                     dimnames = list(paste("SNA", mut_identifiers, sep = ":"), samples_id)))
  
  
  
  
  for (j in 1 : as.numeric(length(mut_identifiers))) {
    
    
    for (g in 1 : as.numeric(length(samples_id))) {
      
      mut_read_counts <- currentMice[(currentMice$`MUTATION ID` == mut_identifiers[[j]]) & 
                                       (currentMice$`SAMPLE ID` == samples_id[[g]])]$`ALT COUNTS`
      
      ref_read_counts <- currentMice[(currentMice$`MUTATION ID` == mut_identifiers[[j]]) & 
                                       (currentMice$`SAMPLE ID` == samples_id[[g]])]$`REF COUNTS`
      
      
      
      r_m <- get(paste("rMatrix", mice[i], sep = ""), envir = .GlobalEnv)
      
      r_m[paste("SNA", mut_identifiers[j], sep = ":"), samples_id[g]] <- as.numeric(mut_read_counts)
      
      
      mode(r_m) <- "integer"
      
      assign(paste("rMatrix", mice[i], sep = ""), r_m)
      
      
      
      x_m <- get(paste("xMatrix", mice[i], sep = ""), envir = .GlobalEnv)
      
      x_m[paste("SNA", mut_identifiers[j], sep = ":"), samples_id[g]] <- as.numeric(mut_read_counts) + as.numeric(ref_read_counts)
      
      
      mode(x_m) <- "integer"
      
      assign(paste("xMatrix", mice[i], sep = ""), x_m)
      
      
    }
    
    
  }
  
  
}


rm(currentMice, r_m, x_m, mut_read_counts, ref_read_counts)



#------------------------------------------------------------------------------------------------------------------------------------------



#####  Define CNA regions (per sample, with merges)  #####


# For all mice together


# CNAS, with copy number != 2
allSamplesCNA <- rbindlist(cnasInfo)


chromosomes_identifiers = unique(allSamplesCNA$CHROM)


# needed to reduce the number of chromosomes visited in the process of defining the CNA regions
chrmRegionCNAS <- list()


h <- 1



for (i in 1 : length(chromosomes_identifiers)) {
  
  
  chrm_region <- allSamplesCNA [allSamplesCNA$CHROM == chromosomes_identifiers[i]]
  
  chrm_region$`COPY ALTERATION ID` <- sapply(chrm_region$`COPY ALTERATION ID`, get_identifier)
  
  
  # means it is a whole genome CNA region
  if(length(unique(chrm_region$`COPY ALTERATION ID`)) == 1){
    
    # [1] since they are all equal
    start_pos <- chrm_region[1]$`START POS`
    end_pos <- chrm_region[1]$`END POS`
    
    
    chrmRegionCNAS[[h]] <- c(as.numeric(start_pos), as.numeric(end_pos))
    
    names(chrmRegionCNAS)[[h]] <- toupper(chromosomes_identifiers[i])
    
    
    h <- h + 1
    
  }
  
  
}



# before getting below, we already have those CNAs which chromosome is present in all samples, so it is only one CNA in that chromosome across all samples





#------------------------------------------------------------------------------------------------------------------------------------------





# For all mice together


#####  Merging/union of CNAs per sample  #####


# The regions found per chromosome (CNA and non-CNA regions)
sample_chrms_list <- list()

toMerge <- list()


# just not to get variable removal warning that it does not exist if this case does not happen

equal_copy <- 0

is_next <- NULL

not_nested <- 0


# install.packages("installr")

# for the is.empty function
suppressPackageStartupMessages(library(installr))



for (i in 1 : as.numeric(length(cnasInfo))) {
  
  sample <- cnasInfo[[i]]
  
  
  # remove from the list the chromosomes we know that only have one CNA in it, across all genome, that are CNA regions  
  toRemove <- tolower(names(chrmRegionCNAS))
  
  # because here is removing by line number
  sample <- sample[ ! (sample$CHROM %in% toRemove),  ]
  
  
  
  for (j in 1: as.numeric(length(unique(sample$CHROM)))) {
    
    
    t <- 1
    
    
    toMerge[[j]] <- list()
    
    sample_chrm_region <- sample[sample$CHROM == unique(sample$CHROM)[j], ]
    
    
    sample_chrm_region$`COPY ALTERATION ID` <- sapply(sample_chrm_region$`COPY ALTERATION ID`, get_identifier)
    
    # to put the indexes from 1 to nrow(sample_chrm_region)
    row.names(sample_chrm_region) <- NULL
    
    
    # we found a CNA region
    if(nrow(sample_chrm_region) == 1) {
      
      toMerge[[j]] <- sample_chrm_region[1, ]$`COPY ALTERATION ID`
      
      names(toMerge)[[j]] <- toupper(sample_chrm_region$CHROM)
      
      
      
    } else {
      
      g <- 1
      
      
      # - 1 since we will not see the last one, we already compared it with the other behind
      while(g <= as.numeric(nrow(sample_chrm_region)) - 1) {
        
        
        current_CNA <- sample_chrm_region[g, ]
        
        
        cna_1 <- as.numeric(current_CNA$`MAJOR COPY`)
        cna_2 <- as.numeric(current_CNA$`MINOR COPY`)
        
        
        next_CNA <- sample_chrm_region[g + 1, ]
        
        
        major_copy <- as.numeric(next_CNA$`MAJOR COPY`)
        minor_copy <- as.numeric(next_CNA$`MINOR COPY`)
        
        
        
        if(cna_1  == major_copy & cna_2 == minor_copy) {
          
          
          # They have the same major and minor. But can they be merged? check difference
          
          
          # less than 10 MB
          if(next_CNA$`START POS` - current_CNA$`END POS` <  10485760 ) {
            
            
            if(g + 1 == as.numeric(nrow(sample_chrm_region))) {
              
              
              if (g > 1) {
                
                
                if(as.numeric(length(toMerge[[j]] [[t - 1]]) > 1)) {
                  
                  
                  if (current_CNA$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [length(toMerge[[j]] [[t - 1]] ) ]) {
                    
                    
                    start_CNA <- toMerge[[j]] [[t - 1]] [1]
                    
                    
                    if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                        ( sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                      
                        
                      
                      toMerge [[j]] [[t - 1]] <- NULL  
                      
                      toMerge [[j]] [[t - 1]] <- c(start_CNA)
                      
                      toMerge [[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                      
                      break
                      
                      
                    } else { 
                      
                      toMerge [[j]] [[t - 1]] <- c(  toMerge [[j]] [[t - 1]], next_CNA$`COPY ALTERATION ID` )
                      
                      toMerge [[j]] [[t]] <- NULL
                      
                      break
                    }
                    
                    
                    
                  } else { 
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                    
                    break
                    
                  }
                  
                  
                  
                } else {
                  
                  
                  if (current_CNA$`COPY ALTERATION ID` != toMerge[[j]] [[t - 1]] [1]) {
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                    
                    break
                    
                    
                  } else {
                    
                    toMerge[[j]] [[t - 1]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)  
                    
                    toMerge [[j]] [[t]] <- NULL
                    
                    break
                    
                  } 
                  
                }
                
                
                
              } else { 
                
                toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                
                break
                
              }
              
              
      
              
              
            } else { # the next CNA is not the final one
              
              
              
              if (g > 1) {
                
                
                if(as.numeric(length(toMerge[[j]] [[t - 1]]) > 1)) {
                  
                  
                  if (current_CNA$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [length(toMerge[[j]] [[t - 1]] )]) { # ok
                    
                    
                    start_CNA <- toMerge[[j]] [[t - 1]] [1]
                    
                    previousToMerge <- toMerge[[j]] [[t - 1]]
                    
                    
                    
                    # the next time with the same copy number appears, if till then, all the others that have distinct copy number are all the same, and all close,
                    # it is a nested segment within other segment
                    
                    
                    if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                        ( sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                      
                      
                      equal_copy = g : (as.numeric(nrow(sample_chrm_region)))
                      
                      equal_copy <- equal_copy[ equal_copy != g]
                      
                      
                      is_next <- as.numeric( equal_copy [ (sample_chrm_region[equal_copy, ]$`MAJOR COPY` == 
                                                             sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == start_CNA, ]$`MAJOR COPY`) &
                                                            (sample_chrm_region[equal_copy, ]$`MINOR COPY` == 
                                                               sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == start_CNA, ]$`MINOR COPY`) ] ) 
                      
                      if(length(is_next) > 1) {
                        
                        is_next <- min(is_next)
                      }
                      
                      
                      not_nested <- 0
                      
                      
                      if(!(is.empty(is_next))) {
                        
                        
                        if( is_next - g > 1) {
                          
                          
                          if(as.numeric(g) + 1 == is_next - 1) {
                            
                            toMerge [[j]] [[t - 1]] <- c( toMerge [[j]] [[t - 1]], sample_chrm_region [is_next - 1, ]$`COPY ALTERATION ID`, sample_chrm_region [is_next, ]$`COPY ALTERATION ID`)
                            
                            g <- is_next
                            
                          } else {
                            
                            
                            for (index in (as.numeric(g) + 1) : (is_next - 1)) {
                              
                              if ((sample_chrm_region [index, ]$`MAJOR COPY` == current_CNA$`MAJOR COPY`) &   
                                  (sample_chrm_region [index, ]$`MINOR COPY` == current_CNA$`MINOR COPY`) & 
                                  sample_chrm_region [index, ]$`START POS` - sample_chrm_region [index - 1, ]$`END POS` <  10485760) {
                                
                                
                                toMerge [[j]] [[t - 1]] <- c( toMerge [[j]] [[t - 1]], sample_chrm_region [index, ]$`COPY ALTERATION ID`)
                                
                                
                                
                              } else { 
                                
                                not_nested <- not_nested + 1
                                
                                toMerge [[j]] [[t - 1]] <- previousToMerge
                              }
                              
                              
                            }
                            
                            
                            
                            if(not_nested == 0) { 
                              
                              toMerge [[j]] [[t - 1]] <- c( toMerge [[j]] [[t - 1]], sample_chrm_region [is_next, ]$`COPY ALTERATION ID`)
                              
                              toMerge [[j]] [[t]] <- NULL
                              
                              g <- is_next
                              
                            }
                            
                            
                          }
                          
                        }
                        
                      } else { 
                        
                        toMerge [[j]] [[t - 1]] <- NULL
                        
                        
                        toMerge[[j]] [[t - 1]] <- c(start_CNA)  
                        
                        toMerge [[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID` ) 
                        
                        
                        t <- t + 1
                        
                        g <- g + 1
                        
                      }
                      
                      
                      
                    } else { 
                      
                      toMerge [[j]] [[t - 1]] <- c(toMerge [[j]] [[t - 1]], next_CNA$`COPY ALTERATION ID`)
                      
                      toMerge [[j]] [[t]] <- NULL
                      
                      g <- g + 1
                      
                    }
                    
                    
                    
                    
                  } else { 
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                    
                    t <- t + 1
                    
                    g <- g + 1
                    
                  }
                  
                  
                  
                } else {  
                  
                  
                  
                  if (current_CNA$`COPY ALTERATION ID` != toMerge[[j]] [[t - 1]] [1]) { 
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                    
                    t <- t + 1
                    
                    g <- g + 1
                    
                    
                  } else { 
                    
                    toMerge[[j]] [[t - 1]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)       
                    
                    toMerge [[j]] [[t]] <- NULL
                    
                    # does not increment t here to stay on the actual index (since the previous index is the one changed)
                    
                    g <- g + 1
                    
                  } 
                  
                  
                }
                
                
                
              } else { 
                
                toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                
                t <- t + 1
                
                g <- g + 1
                
              }
              
              
            }
            
            
            
        
            
          # cna_1 == major_copy & cna_2 == minor_copy, with distance > 10MB
            
            
          } else { 
            
            if(g + 1 == as.numeric(nrow(sample_chrm_region))) {
              
              if ( g > 1 ) {
                
                if(as.numeric(length(toMerge[[j]] [[t - 1]]) > 1)) { # ok
                  
                  if (current_CNA$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [length(toMerge[[j]] [[t - 1]] ) ]) {
                    
                    
                    if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                        ( sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                      
                      
                      toMerge [[j]] [[t - 1]] <- NULL  
                      
                      toMerge [[j]] [[t - 1]] <- c(start_CNA)
                      
                      toMerge [[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                      
                      t <- t + 1
                      
                      toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`) 
                      
                      break
                      
                    } else { 
                      
                      toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`) 
                      
                      break
                      
                    }
                    
                    
                  } else { 
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                    
                    t <- t + 1
                    
                    toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                    
                    break
                    
                  }
                  
                  
                  
                } else { 
                  
                  
                  if (current_CNA$`COPY ALTERATION ID` != toMerge[[j]] [[t - 1]] [1]) {
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                    
                    t <- t + 1
                    
                    toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                    
                    break
                    
                    
                  } else {
                    
                    toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)  
                    
                    break
                    
                  } 
                  
                  
                }     
                
                
                
              } else { 
                
                toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                
                t <- t + 1
                
                toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                
                break
                
              }
              
              
              
              
            # else if the end of the array is not in the following position
              
              
            } else {
              
              if ( g > 1 ) {
                
                if(as.numeric(length(toMerge[[j]] [[t - 1]]) > 1)) { # ok
                  
                  if (current_CNA$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [length(toMerge[[j]] [[t - 1]] ) ]) {
                    
                    
                    if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                        ( sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                      
                      
                      toMerge [[j]] [[t - 1]] <- NULL  
                      
                      toMerge [[j]] [[t - 1]] <- c(start_CNA)
                      
                      toMerge [[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                      
                      t <- t + 1
                      
                      toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`) 
                      
                      g <- g + 1
                      
                    } else { 
                      
                      toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`) 
                      
                      t <- t + 1
                      
                      g <- g + 1
                      
                    }
                    
                    
                    
                  } else { 
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)       
                    
                    t <- t + 1
                    
                    g <- g + 1
                    
                  } 
                  
                  
                  
                } else { 
                  
                  if (current_CNA$`COPY ALTERATION ID` != toMerge[[j]] [[t - 1]] [1]) {
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)       
                    
                    t <- t + 1
                    
                    g <- g + 1
                    
                  } else {
                    
                    toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)       
                    
                    t <- t + 1
                    
                    g <- g + 1
                    
                  }      
                  
                }     
                
                
                
              } else { 
                
                toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                
                t <- t + 1
                
                g <- g + 1
                
              }
              
            }
            
          } 
          
          
          
          
          
          
          #------------------------------------------------------------------------------------------------------------------------------------------
          
        
          
          
          
          
        # if ! (cna_1  == major_copy & cna_2 == minor_copy)
          
        } else { 
          
          
          if(next_CNA$`START POS` - current_CNA$`END POS` <  10485760 ) {
            
            
            
            if(g + 1 == as.numeric(nrow(sample_chrm_region))) {
              
              
              if ( g > 1 ) {
                
                
                if(as.numeric(length(toMerge[[j]] [[t - 1]]) > 1)) {
                  
                  
                  if (current_CNA$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [length(toMerge[[j]] [[t - 1]] )]) {
                    
                    start_CNA <- toMerge[[j]] [[t - 1]] [1]
                    
                    
                    
                    if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                        ( sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                      
                      
                      if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MAJOR COPY` != current_CNA$`MAJOR COPY`) |
                          (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MINOR COPY` != current_CNA$`MINOR COPY`) ) {
                        
                        
                        toMerge [[j]] [[t - 1]] <- NULL  
                        
                        toMerge [[j]] [[t - 1]] <- c(start_CNA)
                        
                        toMerge [[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                        
                        t <- t + 1
                        
                        toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`) 
                        
                        break
                        
                        
                      }
                      
                      else { 
                        
                        toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)         
                        
                        break
                        
                      }
                      
                      
                    } else { 
                      
                      toMerge [[j]] [[t - 1]] <- c(toMerge [[j]] [[t - 1]], next_CNA$`COPY ALTERATION ID` )
                      
                      toMerge [[j]] [[t]] <- NULL
                      
                      break
                      
                    }
                    
                    
                    
                  } else {
                    
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                    
                    t <- t + 1
                    
                    toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                    
                    break
                    
                    
                  }
                  
                  
                  
                  
                } else {
                  
                
                  if (toMerge[[j]] [[t - 1]] [1] != current_CNA$`COPY ALTERATION ID`) {
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                    
                    t <- t + 1
                    
                    toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                    
                    break
                    
                  }
                  
                  
                  else {
                    
   
                    toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                    
                    break
                    
                    
                  }
                  
                  
                  
                  
                }
                
                
                
                
                
              } else { # only has two rows
                
                
                toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                
                t <- t + 1
                
                toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                
                break
                
              }
              
              
              
              
            } else { 
              
              
              if ( g > 1 ) {
                
                
                if(as.numeric(length(toMerge[[j]] [[t - 1]]) > 1)) {
                  
                  if (current_CNA$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [length(toMerge[[j]] [[t - 1]] )]) {
                    
                    start_CNA <- toMerge[[j]] [[t - 1]] [1]
                    
                    previousToMerge <- toMerge[[j]] [[t - 1]]
                    
                    
                    
                    # the next time with the same copy number appears, if till then, all the others that have distinct copy number are all the same, and all close,
                    # it is a nested segment within other segment
                    
                    
                    if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                        ( sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                      
                      
                      equal_copy = g : (as.numeric(nrow(sample_chrm_region)))
                      
                      equal_copy <- equal_copy[ equal_copy != g]
                      
                      
                      is_next <- as.numeric( equal_copy [ (sample_chrm_region[equal_copy, ]$`MAJOR COPY` == 
                                                             sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == start_CNA, ]$`MAJOR COPY`) &
                                                            (sample_chrm_region[equal_copy, ]$`MINOR COPY` == 
                                                               sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == start_CNA, ]$`MINOR COPY`) ] ) 
                      
                      if(length(is_next) > 1) {
                        
                        is_next <- min(is_next)
                      }
                      
                      
                      
                      not_nested <- 0
                      
                      
                      if(!(is.empty(is_next))) {
                        
                        
                        if( is_next - g > 1) {
                          
                          
                          if(as.numeric(g) + 1 == is_next - 1) {
                            
                            toMerge [[j]] [[t - 1]] <- c( toMerge [[j]] [[t - 1]], sample_chrm_region [is_next - 1, ]$`COPY ALTERATION ID`, sample_chrm_region [is_next, ]$`COPY ALTERATION ID`)
                            
                            g <- is_next
                            
                            
                          } else {
                            
                            for (index in (as.numeric(g) + 1) : (is_next - 1)) {
                              
                              if ((sample_chrm_region [index, ]$`MAJOR COPY` == current_CNA$`MAJOR COPY`) &   
                                  (sample_chrm_region [index, ]$`MINOR COPY` == current_CNA$`MINOR COPY`) & 
                                  sample_chrm_region [index, ]$`START POS` - sample_chrm_region [index - 1, ]$`END POS` <  10485760) {
                                
                                
                                toMerge [[j]] [[t - 1]] <- c( toMerge [[j]] [[t - 1]], sample_chrm_region [index, ]$`COPY ALTERATION ID`)
                                
                                
                                
                              } else {
                                
                                not_nested <- not_nested + 1
                                
                                toMerge [[j]] [[t - 1]] <- previousToMerge
                              }
                              
                              
                            }
                            
                            
                            
                            if(not_nested == 0) { 
                              
                              toMerge [[j]] [[t - 1]] <- c( toMerge [[j]] [[t - 1]], sample_chrm_region [is_next, ]$`COPY ALTERATION ID`)
                              
                              toMerge [[j]] [[t]] <- NULL 
                              
                              g <- is_next
                              
                            }
                            
                            
                          }
                          
                        }  
                        
                        
                        
                      } else { 
                        
                        toMerge [[j]] [[t - 1]] <- NULL
                        
                        
                        
                        toMerge[[j]] [[t - 1]] <- c(start_CNA)  
                        
                        toMerge [[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID` ) 
                        
                        
                        t <- t + 1
                        
                        g <- g + 1
                        
                      }
                      
                      
                      
                      
                    } else { 
                      
                      toMerge [[j]] [[t - 1]] <- c(toMerge [[j]] [[t - 1]], next_CNA$`COPY ALTERATION ID`)
                      
                      toMerge [[j]] [[t]] <- NULL
                      
                      g <- g + 1
                      
                    }
                    
                    
                    
                    
                  } else { 
                    
                    c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                    
                    t <- t + 1
                    
                    g <- g + 1
                  }
                  
                  
                  
                } else { 
                  
                  if (current_CNA$`COPY ALTERATION ID` != toMerge[[j]] [[t - 1]] [1]) {
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)     
                    
                    t <- t + 1
                    
                    g <- g + 1
                    
                  } else {
                    
                    
                    toMerge [[j]] [[t - 1]] <- NULL
                    
                    
                    toMerge[[j]] [[t - 1]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)        
                    
                    # t <- t + 1 we remain on the current t, only the previous t was updated
                    
                    g <- g + 1
                    
                  }      
                  
                  
                }
                
                
                
                
                
              } else { 
                
                toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                
                t <- t + 1
                
                g <- g + 1
                
              }
              
              
              
              
            }
            
            
            
            
           
          # cna1 != major copy | cna2 != minor_copy, but their distance is bigger than 10MB 
            
          } else { 
            
            
            if(g + 1 == as.numeric(nrow(sample_chrm_region))) {
              
              
              if ( g > 1 ) {
                
                
                if(as.numeric(length(toMerge[[j]] [[t - 1]]) > 1)) {
                  
                  
                  if (current_CNA$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [length(toMerge[[j]] [[t - 1]] ) ]) {
                    
                    
                    start_CNA <- toMerge[[j]] [[t - 1]] [1]
                    
                    
                    
                    if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                        ( sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                      
                      
                      if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MAJOR COPY` != current_CNA$`MAJOR COPY`) |
                          (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[j]] [[t - 1]] [1], ] $ `MINOR COPY` != current_CNA$`MINOR COPY`) ) {
                        
                        
                        toMerge [[j]] [[t - 1]] <- NULL  
                        
                        toMerge [[j]] [[t - 1]] <- c(start_CNA)
                        
                        toMerge [[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                        
                        t <- t + 1
                        
                        toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`) 
                        
                        break
                        
                        
                      }
                      
                      else { 
                        
                        toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)         
                        
                        break
                        
                      }
                      
                      
                      
                      
                      
                    } else { 
                      
                      
                      toMerge [[j]] [[t - 1]] <- NULL
                      
                      
                      
                      toMerge[[j]] [[t - 1]] <- c(start_CNA)
                      
                      toMerge [[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                      
                      
                      t <- t + 1
                      
                      toMerge [[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID` )
                      
                      
                      break
                      
                    }
                    
                    
                    
                    
                  } else {
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                    
                    t <- t + 1
                    
                    toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                    
                    break
                    
                  }
                  
                  
                } else {
                  
                  if (current_CNA$`COPY ALTERATION ID` != toMerge[[j]] [[t - 1]] [1]) {
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                    
                    t <- t + 1
                    
                    toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                    
                    break
                    
                  } else {
                    
                    toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)
                    
                    break
                    
                  }      
                  
                }    
                
                
                
              
                
              } else { 
                
                toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                
                t <- t + 1
                
                toMerge[[j]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                
                break
                
              }
              
              
              
              
              
              
            } else { # not only has two rows 
              
              
              if ( g > 1 ) {
                
                if(as.numeric(length(toMerge[[j]] [[t - 1]]) > 1)) {
                  
                  if (current_CNA$`COPY ALTERATION ID` != toMerge[[j]] [[t - 1]] [length(toMerge[[j]] [[t - 1]] )]) {
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                    
                    t <- t + 1
                    
                    g <- g + 1
                    
                  } else {
                    
                    g <- g + 1
                    
                  } # else here we go to the next cna without adding
                  
                  
                } else {
                  
                  if (current_CNA$`COPY ALTERATION ID` != toMerge[[j]] [[t - 1]] [1]) {
                    
                    toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)     
                    
                    t <- t + 1
                    
                    g <- g + 1
                    
                  } else {
                    
                    g <- g + 1
                    
                  } # else here we go to the next cna without adding     
                  
                }
                
                
                
                
              } else { 
                
                toMerge[[j]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                
                t <- t + 1
                
                g <- g + 1
                
              }
              
              
              
            }
            
            
            
          }   # else of if cna1 != cna2 and start_pos next - end_pos current > 10MB
          
          
        } # else of if cna1 != cna2   
        
        
      } # g
      
      names(toMerge)[[j]] <- toupper(sample_chrm_region[1, ]$CHROM)
      
    } # else if nrow(chrm_region =! 1)
    
    
    
  } # j
  
  
  
  sample_chrms_list[[i]] <- toMerge
  
  names(sample_chrms_list)[[i]] <- sample[1, ]$`SAMPLE ID` 
  
  toMerge <- NULL
  
  
} # i





# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# For each mice separately


#####  Define all regions (per sample, with merges)


for (i in 1 : as.numeric(length(miceCnasInfo))) {
  
  currentMiceCNAS <- miceCnasInfo[[i]]
  
  chromosomes_identifiers = unique( currentMiceCNAS $CHROM)
  
  
  assign(paste("chrmRegionsCNAS", mice[i], sep = ""), list())
  
  
  h <- 1
  
  
  for (j in 1 : length(chromosomes_identifiers)) {
    
    chrm_region <- currentMiceCNAS [ currentMiceCNAS $CHROM == chromosomes_identifiers[j]]
    
    chrm_region$`COPY ALTERATION ID` <- sapply(chrm_region$`COPY ALTERATION ID`, get_identifier)
    
    # means it is a whole genome CNA region
    if(length(unique(chrm_region$`COPY ALTERATION ID`)) == 1){
      
      start_pos <- chrm_region[1]$`START POS`
      end_pos <- chrm_region[1]$`END POS`
      
      
      
      chr_r <- get(paste("chrmRegionsCNAS", mice[i], sep = ""), envir = .GlobalEnv)
      
      chr_r[[h]] <- c(as.numeric(start_pos), as.numeric(end_pos))
      
      names(chr_r)[[h]] <- toupper(chromosomes_identifiers[j])
      
      
      
      assign(paste("chrmRegionsCNAS", mice[i], sep = ""), chr_r)
      
      
      
      h <- h + 1
      
      
    }
    
  }
  
  
  
}





# before getting below, we already have those CNAs which chromosome is present in all samples, so it is only one CNA in that chromosome across all samples





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------




# For the 4 mice separately


#####  Merging/union of CNAs per sample  #####

toMerge <- list()

mice_chrms_list <- list()


# install.packages("installr")

# suppressPackageStartupMessages(library(installr))
# for the is.empty function


for (i in 1 : as.numeric(length(miceCnasInfo))) {
  
  # The regions found per chromosome (CNA and non-CNA regions)
  assign(paste("sampleChrmsList", mice[i], sep = ""), list())
  
  current_mice <- miceCnasInfo[[i]]
  
  mice_samples <- unique(current_mice$`SAMPLE ID`)
  
  
  chrm_r <- get(paste("chrmRegionsCNAS", mice[i], sep = ""), envir = .GlobalEnv)
  
  toRemove <- tolower(names(chrm_r))  
  
  
  
  for (j in 1 : as.numeric(length(mice_samples))) {
    
    sample <- current_mice[ current_mice$`SAMPLE ID` == mice_samples[j],  ]
    
    sample <- sample[ ! (sample$CHROM %in% toRemove),  ]
    
    
    s_chrms_list <- get(paste("sampleChrmsList", mice[i], sep = ""), envir = .GlobalEnv)
    
    s_chrms_list [[j]] <- list()
    
    
    
    
    for (g in 1: as.numeric(length(unique(sample$CHROM)))) {
      
      t <- 1
      
      toMerge[[g]] <- list()
      
      
      sample_chrm_region <- sample[sample$CHROM == unique(sample$CHROM)[g], ]
      
      
      sample_chrm_region$`COPY ALTERATION ID` <- sapply(sample_chrm_region$`COPY ALTERATION ID`, get_identifier)
      
      # to put the indexes from 1 to nrow(sample_chrm_region)
      row.names(sample_chrm_region) <- NULL
      
      
      
      if( nrow(sample_chrm_region) == 1) {
        
        toMerge[[g]] <- sample_chrm_region[1, ]$`COPY ALTERATION ID`
        
        names(toMerge)[[g]] <- toupper(sample_chrm_region$CHROM)
        
        
      } else {
        
        
        h <- 1
        
        
        
        while(h <= as.numeric(nrow(sample_chrm_region)) - 1) {
          
          
          current_CNA <- sample_chrm_region[h, ]
          
          
          cna_1 <- as.numeric(current_CNA$`MAJOR COPY`)
          cna_2 <- as.numeric(current_CNA$`MINOR COPY`)
          
          
          next_CNA <- sample_chrm_region[h + 1, ]
          
          
          major_copy <- as.numeric(next_CNA$`MAJOR COPY`)
          minor_copy <- as.numeric(next_CNA$`MINOR COPY`)
          
          

          if(cna_1  == major_copy & cna_2 == minor_copy) {
            
            
            # they have the same major and minor. But can they be merged? check difference
            
            
            # less than 10 MB
            if(next_CNA$`START POS` - current_CNA$`END POS` <  10485760 ) {
              
              
              if(h + 1 == as.numeric(nrow(sample_chrm_region))) {
                
                
                if (h > 1) {
                  
                  
                  if(as.numeric(length(toMerge[[g]] [[t - 1]]) > 1)) {
                    
                    
                    if (current_CNA$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [length(toMerge[[g]] [[t - 1]] ) ]) {
                      
                      
                      start_CNA <- toMerge[[g]] [[t - 1]] [1]
                      
                      
                      if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                          ( sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                        
                        
                        toMerge [[g]] [[t - 1]] <- NULL  
                        
                        toMerge [[g]] [[t - 1]] <- c(start_CNA)
                        
                        toMerge [[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                        
                        break
                        
                        
                      } else { 
                        
                        toMerge [[g]] [[t - 1]] <- c(  toMerge [[g]] [[t - 1]], next_CNA$`COPY ALTERATION ID` )
                        
                        toMerge [[g]] [[t]] <- NULL
                        
                        break
                      }
                      
                      
                      
                    } else { 
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                      
                      break
                      
                    }
                    
                    
                    
                  } else { 
                    
                    
                    if (current_CNA$`COPY ALTERATION ID` != toMerge[[g]] [[t - 1]] [1]) {
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                      
                      break
                      
                      
                    } else {
                      
                      toMerge[[g]] [[t - 1]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)  
                      
                      toMerge [[g]] [[t]] <- NULL
                      
                      break
                      
                    } 
                    
                  }
                  
                  
                  
                } else { 
                  
                  toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                  
                  break
                  
                }
                
                
                
                
              } else {
                
                
                
                if (h > 1) {
                  
                  
                  if(as.numeric(length(toMerge[[g]] [[t - 1]]) > 1)) {
                    
                    
                    if (current_CNA$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [length(toMerge[[g]] [[t - 1]] )]) { 
                      
                      
                      start_CNA <- toMerge[[g]] [[t - 1]] [1]
                      
                      previousToMerge <- toMerge[[g]] [[t - 1]]
                      
                      
                      # the next time with the same copy number appears, if till then, all the others that have distinct copy number are all the same, and all close,
                      # it is a nested segment within other segment
                      
                      
                      if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                          (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                        
                        
                        equal_copy = h : (as.numeric(nrow(sample_chrm_region)))
                        
                        equal_copy <- equal_copy[ equal_copy != h]
                        
                        
                        is_next <- as.numeric( equal_copy [ (sample_chrm_region[equal_copy, ]$`MAJOR COPY` == 
                                                               sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == start_CNA, ]$`MAJOR COPY`) &
                                                              (sample_chrm_region[equal_copy, ]$`MINOR COPY` == 
                                                                 sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == start_CNA, ]$`MINOR COPY`) ] ) 
                        
                        if(length(is_next) > 1) {
                          
                          is_next <- min(is_next)
                        }
                        
                        
                        not_nested <- 0
                        
                        
                        if(!(is.empty(is_next))) {
                          
                          
                          if( is_next - h > 1) {
                            
                            
                            if(as.numeric(h) + 1 == is_next - 1) {
                              
                              toMerge [[g]] [[t - 1]] <- c( toMerge [[g]] [[t - 1]], sample_chrm_region [is_next - 1, ]$`COPY ALTERATION ID`, sample_chrm_region [is_next, ]$`COPY ALTERATION ID`)
                              
                              h <- is_next
                              
                              
                            } else {
                              
                              
                              for (index in (as.numeric(h) + 1) : (is_next - 1)) {
                                
                                
                                
                                if ((sample_chrm_region [index, ]$`MAJOR COPY` == current_CNA$`MAJOR COPY`) &   
                                    (sample_chrm_region [index, ]$`MINOR COPY` == current_CNA$`MINOR COPY`) & 
                                    sample_chrm_region [index, ]$`START POS` - sample_chrm_region [index - 1, ]$`END POS` <  10485760) {
                                  
                                  
                                  toMerge [[g]] [[t - 1]] <- c( toMerge [[g]] [[t - 1]], sample_chrm_region [index, ]$`COPY ALTERATION ID`)
                                  
                                  
                                  
                                } else { 
                                  
                                  not_nested <- not_nested + 1
                                  
                                  toMerge [[g]] [[t - 1]] <- previousToMerge
                                }
                                
                                
                              }
                              
                              
                              
                              if(not_nested == 0) { 
                                
                                toMerge [[g]] [[t - 1]] <- c( toMerge [[g]] [[t - 1]], sample_chrm_region [is_next, ]$`COPY ALTERATION ID`)
                                
                                toMerge [[g]] [[t]] <- NULL
                                
                                h <- is_next
                                
                              }
                              
                              
                            }
                            
                          }
                          
                        } else { 
                          
                          toMerge [[g]] [[t - 1]] <- NULL
                          
                          
                          toMerge[[g]] [[t - 1]] <- c(start_CNA)  
                          
                          toMerge [[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID` ) 
                          
                          
                          t <- t + 1
                          
                          h <- h + 1
                          
                        }
                        
                        
                        
                        
                        
                      } else { 
                        
                        toMerge [[g]] [[t - 1]] <- c(toMerge [[g]] [[t - 1]], next_CNA$`COPY ALTERATION ID`)
                        
                        toMerge [[g]] [[t]] <- NULL
                        
                        h <- h + 1
                        
                      }
                      
                      
                      
                      
                    } else { 
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                      
                      t <- t + 1
                      
                      h <- h + 1
                      
                    }
                    
                    
                    
                  } else { 
                    
                    
                    
                    if (current_CNA$`COPY ALTERATION ID` != toMerge[[g]] [[t - 1]] [1]) { # ok
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                      
                      t <- t + 1
                      
                      h <- h + 1
                      
                      
                    } else {
                      
                      toMerge[[g]] [[t - 1]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)       
                      
                      toMerge [[g]] [[t]] <- NULL
                      
                      # does not increment t here to stay on the actual index
                      
                      h <- h + 1
                      
                    } 
                    
                    
                  }
                  
                  
                  
                } else { 
                  
                  toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                  
                  t <- t + 1
                  
                  h <- h + 1
                  
                }
                
                
              }
              
              
              
              
  
            # cna_1 == major_copy & cna_2 == minor_copy, with distance > 10MB
              
              
              
            } else { 
              
              if(h + 1 == as.numeric(nrow(sample_chrm_region))) {
                
                if ( h > 1 ) {
                  
                  if(as.numeric(length(toMerge[[g]] [[t - 1]]) > 1)) { 
                    
                    if (current_CNA$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [length(toMerge[[g]] [[t - 1]] ) ]) {
                      
                      
                      if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                          ( sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                        
                        
                        toMerge [[g]] [[t - 1]] <- NULL  
                        
                        toMerge [[g]] [[t - 1]] <- c(start_CNA)
                        
                        toMerge [[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                        
                        t <- t + 1
                        
                        toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`) 
                        
                        break
                        
                      } else { 
                        
                        toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`) 
                        
                        break
                        
                      }
                      
                      
                    } else { 
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                      
                      t <- t + 1
                      
                      toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                      
                      break
                      
                    }
                    
                    
                    
                  } else {
                    
                    
                    if (current_CNA$`COPY ALTERATION ID` != toMerge[[g]] [[t - 1]] [1]) {
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                      
                      t <- t + 1
                      
                      toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                      
                      break
                      
                      
                    } else {
                      
                      toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)  
                      
                      break
                      
                    } 
                    
                    
                  }     
                  
                  
                  
                } else { 
                  
                  toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                  
                  t <- t + 1
                  
                  toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                  
                  break
                  
                }
                
                
                
                
              # else if the end of the array is not in the following position
                
              } else {
                
                if ( h > 1 ) {
                  
                  if(as.numeric(length(toMerge[[g]] [[t - 1]]) > 1)) { 
                    
                    if (current_CNA$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [length(toMerge[[g]] [[t - 1]] ) ]) {
                      
                      
                      if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                          (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                        
                        
                        toMerge [[g]] [[t - 1]] <- NULL  
                        
                        toMerge [[g]] [[t - 1]] <- c(start_CNA)
                        
                        toMerge [[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                        
                        t <- t + 1
                        
                        toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`) 
                        
                        h <- h + 1
                        
                        
                      } else { 
                        
                        toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`) 
                        
                        t <- t + 1
                        
                        h <- h + 1
                        
                      }
                      
                      
                      
                    } else { 
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)       
                      
                      t <- t + 1
                      
                      h <- h + 1
                      
                    } 
                    
                    
                    
                  } else { 
                    
                    if (current_CNA$`COPY ALTERATION ID` != toMerge[[g]] [[t - 1]] [1]) {
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)       
                      
                      t <- t + 1
                      
                      h <- h + 1
                      
                    } else {
                      
                      toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)       
                      
                      t <- t + 1
                      
                      h <- h + 1
                      
                    }      
                    
                  }     
                  
                  
                  
                } else { 
                  
                  toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                  
                  t <- t + 1
                  
                  h <- h + 1
                  
                }
                
              }
              
            } 
            
            
            
            
            
          # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            
     
            
          # if ! (cna_1  == major_copy & cna_2 == minor_copy)
            
            
          } else {
            
            
            if(next_CNA$`START POS` - current_CNA$`END POS` <  10485760 ) {
              
              
              
              if(h + 1 == as.numeric(nrow(sample_chrm_region))) {
                
                
                if ( h > 1 ) {
                  
                  
                  if(as.numeric(length(toMerge[[g]] [[t - 1]]) > 1)) {
                    
                    
                    if (current_CNA$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [length(toMerge[[g]] [[t - 1]] )]) {
                      
                      start_CNA <- toMerge[[g]] [[t - 1]] [1]
                      
                      
                      if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                          ( sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                        
                        
                        if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MAJOR COPY` != current_CNA$`MAJOR COPY`) |
                            (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MINOR COPY` != current_CNA$`MINOR COPY`) ) {
                          
                          
                          toMerge [[g]] [[t - 1]] <- NULL  
                          
                          toMerge [[g]] [[t - 1]] <- c(start_CNA)
                          
                          toMerge [[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                          
                          t <- t + 1
                          
                          toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`) 
                          
                          break
                          
                          
                        }
                        
                        else { 
                          
                          toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)         
                          
                          break
                          
                        }
                        
                        
                      } else { 
                        
                        toMerge [[g]] [[t - 1]] <- c(toMerge [[g]] [[t - 1]], next_CNA$`COPY ALTERATION ID` )
                        
                        toMerge [[g]] [[t]] <- NULL
                        
                        break
                        
                      }
                      
                      
                      
                    } else {
                      
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                      
                      t <- t + 1
                      
                      toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                      
                      break
                      
                      
                    }
                    
                    
                    
                    
                  } else {
                    
                    
                    if (toMerge[[g]] [[t - 1]] [1] != current_CNA$`COPY ALTERATION ID`) {
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                      
                      t <- t + 1
                      
                      toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                      
                      break
                      
                    }
                    
                    
                    else {
                      
                      toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                      
                      break
                      
                      
                    }
                    
                    
                    
                  }
                  
                  
                  
                  
           
                } else { # only has two rows
                  
                  
                  toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                  
                  t <- t + 1
                  
                  toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                  
                  break
                  
                }
                
                
                
                
                
              } else { 
                
                
                if ( h > 1 ) {
                  
                  
                  if(as.numeric(length(toMerge[[g]] [[t - 1]]) > 1)) {
                    
                    if (current_CNA$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [length(toMerge[[g]] [[t - 1]] )]) {
                      
                      
                      start_CNA <- toMerge[[g]] [[t - 1]] [1]
                      
                      previousToMerge <- toMerge[[g]] [[t - 1]]
                      
                      
                      # the next time with the same copy number appears, if till then, all the others that have distinct copy number are all the same, and all close,
                      # it is a nested segment within other segment
                      
                      
                      if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                          ( sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                        
                        
                        equal_copy = h : (as.numeric(nrow(sample_chrm_region)))
                        
                        equal_copy <- equal_copy[ equal_copy != h]
                        
                        
                        is_next <- as.numeric( equal_copy [ (sample_chrm_region[equal_copy, ]$`MAJOR COPY` == 
                                                               sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == start_CNA, ]$`MAJOR COPY`) &
                                                              (sample_chrm_region[equal_copy, ]$`MINOR COPY` == 
                                                                 sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == start_CNA, ]$`MINOR COPY`) ] ) 
                        
                        if(length(is_next) > 1) {
                          
                          is_next <- min(is_next)
                        }
                        
                        
                        
                        not_nested <- 0
                        
                        
                        if(!(is.empty(is_next))) {
                          
                          
                          if( is_next - h > 1) {
                            
                            
                            if(as.numeric(h) + 1 == is_next - 1) {
                              
                              toMerge [[g]] [[t - 1]] <- c( toMerge [[g]] [[t - 1]], sample_chrm_region [is_next - 1, ]$`COPY ALTERATION ID`, sample_chrm_region [is_next, ]$`COPY ALTERATION ID`)
                              
                              h <- is_next
                              
                              
                            } else {
                              
                              
                              for (index in (as.numeric(h) + 1) : (is_next - 1)) {
                                
                                
                                if ((sample_chrm_region [index, ]$`MAJOR COPY` == current_CNA$`MAJOR COPY`) &   
                                    (sample_chrm_region [index, ]$`MINOR COPY` == current_CNA$`MINOR COPY`) & 
                                    sample_chrm_region [index, ]$`START POS` - sample_chrm_region [index - 1, ]$`END POS` <  10485760) {
                                  
                                  
                                  toMerge [[g]] [[t - 1]] <- c( toMerge [[g]] [[t - 1]], sample_chrm_region [index, ]$`COPY ALTERATION ID`)
                                  
                                  
                                  
                                } else {
                                  
                                  not_nested <- not_nested + 1
                                  
                                  toMerge [[g]] [[t - 1]] <- previousToMerge
                                }
                                
                                
                              }
                              
                              
                              
                              if(not_nested == 0) { 
                                
                                toMerge [[g]] [[t - 1]] <- c( toMerge [[g]] [[t - 1]], sample_chrm_region [is_next, ]$`COPY ALTERATION ID`)
                                
                                toMerge [[g]] [[t]] <- NULL 
                                
                                h <- is_next
                                
                              }
                              
                              
                            }
                            
                          }  
                          
                        } else { 
                          
                          toMerge [[g]] [[t - 1]] <- NULL
                          
                          
                          toMerge [[g]] [[t - 1]] <- c(start_CNA)  
                          
                          toMerge [[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID` ) 
                          
                          
                          t <- t + 1
                          
                          h <- h + 1
                          
                        }
                        
                        
                        
                        
                      } else { 
                        
                        toMerge [[g]] [[t - 1]] <- c(toMerge [[g]] [[t - 1]], next_CNA$`COPY ALTERATION ID`)
                        
                        toMerge [[g]] [[t]] <- NULL
                        
                        h <- h + 1
                        
                      }
                      
                      
                      
                      
                    } else { 
                      
                      c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                      
                      t <- t + 1
                      
                      h <- h + 1
                    }
                    
                    
                    
                  } else { 
                    
                    if (current_CNA$`COPY ALTERATION ID` != toMerge[[g]] [[t - 1]] [1]) {
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)     
                      
                      t <- t + 1
                      
                      h <- h + 1
                      
                      
                    } else {
                      
                      
                      toMerge [[g]] [[t - 1]] <- NULL
                      
                      
                      toMerge[[g]] [[t - 1]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)        
                      
                      # t <- t + 1  to stay on the same index, since we only updated the previous index position
                      
                      h <- h + 1
                      
                    }      
                    
                    
                  }
                  
                  
                  
                  
                  
                } else { 
                  
                  toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`, next_CNA$`COPY ALTERATION ID`)
                  
                  t <- t + 1
                  
                  h <- h + 1
                  
                }
                
                
                
              }
              
              
              
              
              
              
            # cna1 != major copy | cna2 != minor_copy, but their distance is bigger than 10MB 
              
            } else { 
              
              
              if(h + 1 == as.numeric(nrow(sample_chrm_region))) {
                
                
                if ( h > 1 ) {
                  
                  
                  if(as.numeric(length(toMerge[[g]] [[t - 1]]) > 1)) {
                    
                    
                    if (current_CNA$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [length(toMerge[[g]] [[t - 1]] ) ]) {
                      
                      
                      start_CNA <- toMerge[[g]] [[t - 1]] [1]
                      
                      
                      if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MAJOR COPY` != next_CNA$`MAJOR COPY`) | 
                          ( sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MINOR COPY` != next_CNA$`MINOR COPY`) ) {
                        
                        
                        if( (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MAJOR COPY` != current_CNA$`MAJOR COPY`) |
                            (sample_chrm_region [ sample_chrm_region$`COPY ALTERATION ID` == toMerge[[g]] [[t - 1]] [1], ] $ `MINOR COPY` != current_CNA$`MINOR COPY`) ) {
                          
                          
                          toMerge [[g]] [[t - 1]] <- NULL  
                          
                          toMerge [[g]] [[t - 1]] <- c(start_CNA)
                          
                          toMerge [[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                          
                          t <- t + 1
                          
                          toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`) 
                          
                          break
                          
                          
                        }
                        
                        else {
                          
                          toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)         
                          
                          break
                          
                        }
                        
                        
                        
                        
                      } else { 
                        
                        
                        toMerge [[g]] [[t - 1]] <- NULL
                        
                        
                        toMerge [[g]] [[t - 1]] <- c(start_CNA)
                        
                        toMerge [[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                        
                        
                        t <- t + 1
                        
                        toMerge [[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID` )
                        
                        
                        break
                        
                      }
                      
                      
                      
                      
                    } else {
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                      
                      t <- t + 1
                      
                      toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                      
                      break
                      
                    }
                    
                    
                    
                    
                    
                  } else {
                    
                    if (current_CNA$`COPY ALTERATION ID` != toMerge[[g]] [[t - 1]] [1]) {
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)
                      
                      t <- t + 1
                      
                      toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                      
                      break
                      
                    } else {
                      
                      toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)
                      
                      break
                      
                    }      
                    
                  }    
                  
                  
                  
                  
                  
                } else { 
                  
                  toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                  
                  t <- t + 1
                  
                  toMerge[[g]] [[t]] <- c(next_CNA$`COPY ALTERATION ID`)              
                  
                  break
                  
                }
                
                
                
                
                
                
              } else { # else not only has two rows 
                
                
                if ( h > 1 ) {
                  
                  if(as.numeric(length(toMerge[[g]] [[t - 1]]) > 1)) {
                    
                    if (current_CNA$`COPY ALTERATION ID` != toMerge[[g]] [[t - 1]] [length(toMerge[[g]] [[t - 1]] )]) {
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                      
                      t <- t + 1
                      
                      h <- h + 1
                      
                    } else {
                      
                      h <- h + 1
                      
                    } # else here we go to the next cna without adding
                    
                    
                  } else {
                    
                    if (current_CNA$`COPY ALTERATION ID` != toMerge[[g]] [[t - 1]] [1]) {
                      
                      toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)     
                      
                      t <- t + 1
                      
                      h <- h + 1
                      
                    } else {
                      
                      h <- h + 1
                      
                    } # else here we go to the next cna without adding     
                    
                  }
                  
                  
                  
                  
                } else { 
                  
                  toMerge[[g]] [[t]] <- c(current_CNA$`COPY ALTERATION ID`)  
                  
                  t <- t + 1
                  
                  h <- h + 1
                  
                }
                
                
                
              }
              
              
              
            }   # else of if cna1 != cna2 and start_pos next - end_pos current > 10MB
            
            
            
          } # else of if cna1 != cna2   
          
          
        } # h
        
        names(toMerge)[[g]] <- toupper(sample_chrm_region[1, ]$CHROM)  
        
        
        
      } # else if nrow(chrm_region =! 1)
      
      
    } # g
    
    
    s_chrms_list [[j]] <- toMerge
    
    names(s_chrms_list) [[j]]  <- sample[1, ]$`SAMPLE ID` 
    
    
    assign(paste("sampleChrmsList", mice[i], sep = ""), s_chrms_list)
    
    
    toMerge <- NULL
    
    
    
  } # j 
  
  
  mice_chrms_list[[i]] <- s_chrms_list
  
  names(mice_chrms_list)[[i]] <- mice[i]
  
  
} # i





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------




#####  Calculate unions of nested CNAs resulting in the new CNA  #####


# For all mice together


chrms_samples <- names(sample_chrms_list)

k <- 1


union_cnas <- list()



for(i in 1 : as.numeric(length(chrms_samples))) {
  
  current_sample <- chrms_samples[i]
  
  chrms_to_iterate <- names(sample_chrms_list[[current_sample]])
  
  new_cnas <- list()
  
  union_cnas[[i]] <- cnasInfo[[i]]
  
  
  
  for (j in 1 : as.numeric(length(chrms_to_iterate))) {
    
    
    current_chrm <- sample_chrms_list[[current_sample]][[chrms_to_iterate[j]]]
    
    
    for(g in 1 : as.numeric(length(current_chrm))) {
      
      
      if(length(current_chrm[[g]]) < 2) {
        
        next
      }
      
      
      # calculate new CNA
      else {
        
        chrm_cnas <- current_chrm[[g]]
        
        
        cna_size <- as.numeric(length(chrm_cnas))
        
        
        start_pos <- as.numeric(substr(chrm_cnas[1], unlist(gregexpr(":", chrm_cnas[1]))[1] + 1,  unlist(gregexpr(":", chrm_cnas[1]))[2] - 1))
        
        end_pos <- as.numeric(substr(chrm_cnas[length(chrm_cnas)],
                                     unlist(gregexpr(":", chrm_cnas[length(chrm_cnas)]))[2] + 1,  nchar(chrm_cnas[length(chrm_cnas)]) ))
        
        
        cna_range <- end_pos - start_pos
        
        
        cna_sample <- cnasInfo[[i]]
        
        
        cna_sample <- cna_sample [cna_sample$`CHROM` == tolower(chrms_to_iterate[j]), ]
        
        
        
        cna_sample$`COPY ALTERATION ID` <- sapply(cna_sample$`COPY ALTERATION ID`, get_identifier)
        
        # to put the indexes from 1 to nrow(sample_chrm_region)
        row.names(cna_sample) <- NULL
        
        
        first_cna <- cna_sample[cna_sample$`COPY ALTERATION ID` == chrm_cnas[1], ]
        
        
        cna_1 <- first_cna$`MAJOR COPY`
        
        cna_2 <- first_cna$`MINOR COPY`
        
        
        
        union_cnas[[i]] <- union_cnas[[i]] [! (sapply(union_cnas[[i]]$`COPY ALTERATION ID`, get_identifier) %in% first_cna), ]
        
        
        
        # does not start in one cause we want to compare the copy number with the next CNA to see if it is nested
        for (h in 2 : cna_size) {
          
          current_cna <- cna_sample[cna_sample$`COPY ALTERATION ID` == chrm_cnas[h], ]
          
          
          if (current_cna$`MAJOR COPY` == cna_1 & current_cna$`MINOR COPY` == cna_2) {
            
            union_cnas[[i]] <- union_cnas[[i]] [! (sapply(union_cnas[[i]]$`COPY ALTERATION ID`, get_identifier) %in% current_cna), ]
            
          }
          
          
          # else, we do not remove, is a nested
          
        }
        
        
        
        new_cnas[[k]] <- c('COPY ALTERATION ID' = paste(paste(paste(chrms_to_iterate[j], start_pos, sep = ":"), end_pos, sep = ":"), current_sample, sep = ":"),
                           'SAMPLE ID' = current_sample , 'CHROM' = tolower(chrms_to_iterate[j]), 'START POS' = start_pos, 'END POS' = end_pos, 
                           # this is just to have the same columns as cnasInfo
                           'COPY NUMBER' = 0, 'MAJOR COPY' = 0, 'MINOR COPY' = 0, 'CI LOWER' = 0, 'CI UPPER' = 0, 'BAF' = 0, 'LOG2' = 0)  
        
        
        k <- k + 1
        
        
      }
      
      
    }
    
    
  }
  
  
  df <- as.data.frame(do.call(rbind, new_cnas))
  
  
  union_cnas[[i]] <- rbind(union_cnas[[i]], df)
  
  
  
  union_cnas[[i]]$'CHRM NUMBER' <- sapply(union_cnas[[i]]$`COPY ALTERATION ID`, function(x) as.numeric(substr(x,
                                                                                                              unlist(gregexpr("R", x))[1] + 1,  unlist(gregexpr(":", x))[1] - 1)))
  
  
  union_cnas[[i]] <- union_cnas[[i]] [order(as.numeric(union_cnas[[i]]$`CHRM NUMBER`), as.numeric(union_cnas[[i]]$`START POS`)) , ]
  
  
  row.names(union_cnas[[i]]) <- NULL
  
  
}





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------



# For all mice together


#####  Identify CNA regions  #####


union_cnas <- rbindlist(union_cnas)

union_cnas$`COPY ALTERATION ID` <- sapply(union_cnas$`COPY ALTERATION ID`, get_identifier)



cna_regions <- c()

chrm_number <- ""



for(i in 1 : as.numeric(length(unique(allSamplesCNA$CHROM)))) { 
  
  current_chrm <- union_cnas [union_cnas$CHROM == unique(allSamplesCNA$CHROM)[i], ]
  
  current_chrm <- current_chrm[!duplicated(current_chrm$`COPY ALTERATION ID`), ]
  
  
  chrm_number <- current_chrm[1, ]$CHROM
  
  
  
  if(length(unique(current_chrm$`COPY ALTERATION ID`)) == 1 ) {
    
    cna_regions <- c(cna_regions, paste("R", unique(current_chrm$`COPY ALTERATION ID`), sep = ":"))
    
    
    
  } else {
    
    max_iterations <- length(current_chrm)
    
    
    for(g in 1 : max_iterations){
      
      region_start <- min(as.numeric(current_chrm$`START POS`))
      
      region_end <- max(as.numeric(current_chrm[current_chrm$`START POS` == region_start, ]$`END POS`))
      
      
      next_intersection <- current_chrm
      
      
      
      while ( nrow(next_intersection) != 0 ) {
        
        min_start_pos <- min(as.numeric(next_intersection$`START POS`))
        
        max_connection <- max(as.numeric(next_intersection[as.numeric(next_intersection$`START POS`) == min_start_pos, ]$`END POS`))
        
        
        current_chrm <- current_chrm [ ! (as.numeric(current_chrm$`START POS`) >= min_start_pos & as.numeric(current_chrm$`END POS`) <= max_connection)  ,  ]
        
        
        next_intersection <- current_chrm[as.numeric(current_chrm$`START POS`) < max_connection, ]
        
      }
      
      
      
      
      if(nrow(next_intersection) == 0) { 
        
        region_end <- max_connection
        
        cna_regions <- c(cna_regions, paste("R", paste(toupper(chrm_number), paste(toString(region_start), toString(region_end), sep = ":"), sep = ":"), sep = ":"))
        
        
        if(nrow(current_chrm) == 0) {
          
          break
          
        }
        
        
      }
      
      
    }
    
    
  }
  
  
}   





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------




#####  Construct C Matrix  #####


# For all mice together



# All the distinct CNAS that exist, a CNA may be present in all samples

# union_cnas has all (even the repeated ones, the cna with the same copy alteration ID in different samples)


unique_cnas <- union_cnas[!duplicated(union_cnas$`COPY ALTERATION ID`), ]



unique_cnas <- unique_cnas[ order(as.numeric(unique_cnas$`CHRM NUMBER`), as.numeric(unique_cnas$`START POS`), as.numeric(unique_cnas$`END POS`)),  ]



cna_regions <- cna_regions [ order( as.numeric(sapply(cna_regions, function(x) as.numeric(substr(x,
                                                                                                 unlist(gregexpr("R", x))[2] + 1,  unlist(gregexpr(":", x))[2] - 1)))),                     
                                    as.numeric(sapply(cna_regions, function(x) as.numeric(substr(x,
                                                                                                 unlist(gregexpr(":", x))[2] + 1,  unlist(gregexpr(":", x))[3] - 1)))),
                                    as.numeric(sapply(cna_regions, function(x) as.numeric(substr(x,
                                                                                                 unlist(gregexpr(":", x))[3] + 1,  nchar(x) - 1)))) ) ]





cMatrix <- matrix(0, nrow = as.numeric(length(cna_regions)), ncol = as.numeric(length(unique_cnas$`COPY ALTERATION ID`)), 
                  byrow = TRUE, dimnames = list(cna_regions, paste("CNA", unique_cnas$`COPY ALTERATION ID`, sep = ":")))





for(i in 1 : as.numeric(length(unique(unique_cnas$CHROM)))) { 
  
  current_chrm <- unique_cnas [unique_cnas$CHROM == unique(unique_cnas$CHROM)[i], ]
  
  row.names(current_chrm) <- NULL
  
  
  chrm_regions <- cna_regions[ sapply(cna_regions, function (x) 
    substr(x, unlist(gregexpr(":", x))[1] + 1, unlist(gregexpr(":", x))[2] - 1) == toupper(unique(unique_cnas$CHROM)[i]))]
  
  
  
  # iterate CNA regions of a chromosome
  for (j in 1 : as.numeric(length(chrm_regions))) {
    
    current_region <- chrm_regions[j]
    
    region_start <- as.numeric(substr(current_region, unlist(gregexpr(":", current_region))[2] + 1, unlist(gregexpr(":", current_region))[3] - 1))
    
    region_end <- as.numeric(substr(current_region, unlist(gregexpr(":", current_region))[3] + 1, nchar(current_region)))
    
    
    
    # iterate CNAs in that chromosome
    for (g in 1 : nrow(current_chrm) ) {
      
      current_cna <- current_chrm[g,  ]
      
      start_pos <- as.numeric(current_cna$`START POS`)
      
      end_pos <- as.numeric(current_cna$`END POS`)
      
      
      # The CNA is out of the region limits, so it does not intersect it
      if ((region_start > end_pos - 0.05 & region_start > start_pos + 0.05) | 
          (region_end < start_pos + 0.05 & region_end < end_pos - 0.05)) {
        
        g <- g + 1
        next
        
        
        
      # When the CNA start and end is the same as the region (CNA stays completely within the bounds of the region, in the middle) 
      } else if ((start_pos + 0.05 > region_start & end_pos - 0.05 < region_end) |
                 
                 # When the CNA intersects a part of the region, and it is to its left or right 
                 (end_pos - 0.05 > region_start & start_pos + 0.05 < region_end)) {
        
        
        cna_name <- paste("CNA", current_cna$`COPY ALTERATION ID`, sep = ":")
        
        cMatrix [current_region, cna_name] <- 1
        
      }
      
      
    }
    
    
  }
  
  
}





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------






#####  Matrices Wm, WM and Epsilon Wm and WM  #####

# Formulas to use: 
# BAF = Wm / (WM + Wm)
# Log2 (Depth.ratio) = (WM + Wm)/2


samples_id <- unique(allSamples$`SAMPLE ID`)




# Wm (then assign minorCopyMatrix to Wm in the algorithm parameters)
minorCopyMatrix <- matrix(0, nrow = as.numeric(length(row.names(cMatrix))), ncol = as.numeric(length(samples_id)),
                          byrow = TRUE, dimnames = list(row.names(cMatrix), samples_id))

# WM (then assign majorCopyMatrix to WM in the algorithm parameters)
majorCopyMatrix <- matrix(0, nrow = as.numeric(length(row.names(cMatrix))), ncol = as.numeric(length(samples_id)),
                          byrow = TRUE, dimnames = list(row.names(cMatrix), samples_id))




# epsilon Wm
minorEpsilon <- matrix(0, nrow = as.numeric(length(row.names(cMatrix))), ncol = as.numeric(length(samples_id)),
                       byrow = TRUE, dimnames = list(row.names(cMatrix), samples_id))

# epsilon WM
majorEpsilon <- matrix(0, nrow = as.numeric(length(row.names(cMatrix))), ncol = as.numeric(length(samples_id)),
                       byrow = TRUE, dimnames = list(row.names(cMatrix), samples_id)) 





# Non-CNAs (copy number == 2) and CNAs (copy number != 2)
allCopyNumber <- rbindlist(copyNumberInfo)


allCopyNumber$`COPY ALTERATION ID` <- sapply(allCopyNumber$`COPY ALTERATION ID`, get_identifier)




for (i in 1 : as.numeric(length(row.names(cMatrix)))) {
  
  
  currentCNARegion <- row.names(cMatrix)[i]
  
  
  region_start <- as.numeric(substr(currentCNARegion, unlist(gregexpr(":", currentCNARegion))[2] + 1, unlist(gregexpr(":", currentCNARegion))[3] - 1 ))
  
  region_end <- as.numeric(substr(currentCNARegion, unlist(gregexpr(":", currentCNARegion))[3] + 1, as.numeric(nchar(currentCNARegion))))
  
  
  currentChr <- substr(currentCNARegion, unlist(gregexpr(":", currentCNARegion))[1] + 1, unlist(gregexpr(":", currentCNARegion))[2] - 1 )
  
  
  currentChrCNAS <- allCopyNumber[allCopyNumber$CHROM == tolower(currentChr), ]
  
  
  
  for (j in 1 : as.numeric(length(samples_id))) {
    
    sample <- samples_id[j]
    
    currentChrSample <- currentChrCNAS[currentChrCNAS$`SAMPLE ID` == sample, ]
    
    
    toCompareWith <- list()
    
    
    
    # When the CNA (or copy number info) start and end is the same as the region (CNA stays completely within the bounds of the region, in the middle) 
    # Or when the CNA (or copy number info) intersects a part of the region, and it is to its left or right 
    toCompareWith <- currentChrSample[ (as.numeric(currentChrSample$`START POS`) + 0.05 > region_start & as.numeric(currentChrSample$`END POS`) - 0.05 < region_end) |
                                         (as.numeric(currentChrSample$`END POS`) - 0.05 > region_start & as.numeric(currentChrSample$`START POS`) + 0.05 < region_end) , ]
    
    
    
    # just for verification purposes, this should not happen
    if(nrow(toCompareWith) == 0){
      j <- j + 1
      next
    }
    
    
    if(nrow(toCompareWith) == 1) {
      
      # BAF = Wm / (WM + Wm)
      # Depth.ratio = (WM + Wm)/2
      
      baf <- toCompareWith$BAF
      
      # depth_ratio = 2 ^ log2 column
      depth_ratio <- 2^(as.numeric(toCompareWith$LOG2))
      
      
      # Wm
      minor_copy <- 2 * baf * depth_ratio
      
      # WM
      major_copy <- 2 * depth_ratio - minor_copy
      
      
      if (baf > 0.5) { # major_copy is going to be < minor_copy
        minorCopyMatrix[currentCNARegion, toCompareWith$`SAMPLE ID`] <- major_copy
        majorCopyMatrix[currentCNARegion, toCompareWith$`SAMPLE ID`] <- minor_copy
        
      } else {
        
        # add minor_copy to position [i, j] of minorCopyMatrix
        # add major_copy to position [i, j] of majorCopyMatrix
        minorCopyMatrix[currentCNARegion, toCompareWith$`SAMPLE ID`] <- minor_copy
        majorCopyMatrix[currentCNARegion, toCompareWith$`SAMPLE ID`] <- major_copy
        
      }
      
      
      
      # SE = (upper limit - lower limit) / 3.92.
      
      # Confidence interval of the segment mean (--ci), estimated by bootstrap (100 resamplings) of the bin-levellog2 ratio values 
      # within the segment.  The upper and lower bounds are output as separate columns ci_lo and ci_hi
      # That is why it is 2 ^ the ci_lo and ci_hi
      
      
      ci_lower <- 2^toCompareWith$`CI LOWER`
      
      ci_upper <- 2^toCompareWith$`CI UPPER`
      
      
      epsilon <- (ci_upper - ci_lower) / 3.92        
      
      
      
      minorEpsilon[currentCNARegion, toCompareWith$`SAMPLE ID`]  <- epsilon
      majorEpsilon[currentCNARegion, toCompareWith$`SAMPLE ID`]  <- epsilon
      
      
      
      
      
    } else { # cna region intersects different CNAs
      
      
      # weighted average 
      
      ci_lower <- ci_upper <- baf <- depth_ratio <- 0 
      
      
      # Calculate which intersection is bigger, with which CNA info
      for (g in 1 : as.numeric(nrow(toCompareWith))) { 
        
        cna_intersecting <- toCompareWith[g]
        
        # we are going use the info that is there, even if the green intersects a part of that 
        # inside a sample none of the cnas in that sample intersects one another
        
        
        # CNA stays in between region:
        if(toCompareWith[g]$`START POS` + 0.05 > region_start & toCompareWith[g]$`END POS` - 0.05 < region_end) { 
          
          cna_range <- as.numeric(cna_intersecting$`END POS`) - as.numeric(cna_intersecting$`START POS`)
          
          # percentage of the area covered by that CNA
          p <- cna_range / (region_end - region_start)
          
        }
        
        
        # When the CNA intersects a part of the region, and it is to its left or right 
        else {
          
          if (toCompareWith[g]$`END POS` - 0.05 > region_start & toCompareWith[g]$`START POS` + 0.05 < region_end) {
            
            
            # Region stays to the left of the CNA (or copy number info)
            if(toCompareWith[g]$`START POS` + 0.05 > region_start) {
              
              cna_range <- region_end - as.numeric(cna_intersecting$`START POS`)
              
              p <- cna_range / (region_end - region_start)
              
            } 
            
            
            # Region stays to the right of the CNA (or copy number info)
            else if (toCompareWith[g]$`END POS` - 0.05 < region_end){
              
              cna_range <- as.numeric(cna_intersecting$`END POS`) - region_start
              
              
              p <- cna_range / (region_end - region_start)
              
            }
            
            
            # CNA is bigger than the region, and so the region stays contained inside it
            else {
              
              # (region_end - region_start) / (region_end - region_start)
              p <- 1
              
            }
            
            
          }
          
          
        }
        
        
        
        baf <- baf + (p * as.numeric(cna_intersecting$`BAF`)) 
        
        depth_ratio <- depth_ratio + (p * (2^as.numeric(cna_intersecting$`LOG2`) ))  
        
        
        ci_lower <- ci_lower + (p * (2^as.numeric(cna_intersecting$`CI LOWER`) ))          
        
        ci_upper <- ci_upper + (p * (2^as.numeric(cna_intersecting$`CI UPPER`) ))     
        
        
      }  
      
      
      
      # weigthed average
      
      baf <- baf / nrow(toCompareWith)
      
      depth_ratio <- depth_ratio / nrow(toCompareWith)
      
      
      ci_lower <- ci_lower  / nrow(toCompareWith)
      
      ci_upper <- ci_upper / nrow(toCompareWith)
      
      
      
      
      # Wm
      minor_copy <- 2 * baf * depth_ratio
      
      # WM
      major_copy <- 2 * depth_ratio - minor_copy
      
      
      
      if (baf > 0.5) { # major_copy is going to be < minor_copy
        minorCopyMatrix[currentCNARegion, cna_intersecting$`SAMPLE ID`] <- major_copy
        majorCopyMatrix[currentCNARegion, cna_intersecting$`SAMPLE ID`] <- minor_copy
        
      } else {
        
        # add minor_copy to position [i, j] of minorCopyMatrix
        # add major_copy to position [i, j] of majorCopyMatrix
        minorCopyMatrix[currentCNARegion, cna_intersecting$`SAMPLE ID`] <- minor_copy
        majorCopyMatrix[currentCNARegion, cna_intersecting$`SAMPLE ID`] <- major_copy
        
      }
      
      
      
      # SE = (upper limit - lower limit) / 3.92.
      
      # Confidence interval of the segment mean (--ci), estimated by bootstrap (100 resamplings) of the bin-levellog2 ratio values 
      # within the segment.  The upper and lower bounds are output as separate columns ci_lo and ci_hi
      # That is why it is 2 ^ the ci_lo and ci_hi
      
      
      # lower_limit <- 2^toCompareWith[toCompareWith$`COPY ALTERATION ID` == cna_intersecting$`COPY ALTERATION ID`, ]$`CI LOWER`
      
      # upper_limit <- 2^toCompareWith[toCompareWith$`COPY ALTERATION ID` == cna_intersecting$`COPY ALTERATION ID`, ]$`CI UPPER`
      
      
      epsilon <- (ci_upper - ci_lower) / 3.92        
      
      
      minorEpsilon[currentCNARegion, cna_intersecting$`SAMPLE ID`]  <- epsilon
      majorEpsilon[currentCNARegion, cna_intersecting$`SAMPLE ID`]  <- epsilon
      
      
    }
    
  }
  
}  





# After this check Wm and WM and see if is better to remove some CNA regions. Then, remove those CNA regions from the C matrix, and use only those CNA regions for Y, 
# but before doing Y also verify if some SNVs must be removed, using a heatmap, which also means the SNVs removed need to be removed from the R and X matrices

toRemove <- c()

# major and minor copy matrix have the same number of rows so it can be any of both nrow
for (i in 1: as.numeric(nrow(minorCopyMatrix))) {
  
  # since both rows with all zeros in minorCopyMatrix and majorCopyMatrix are the same, we just need to calculate for one of them
  r_content <- lapply(minorCopyMatrix[i, ], function(x) 0 %in% x) 
  
  hasZero <- which (r_content == TRUE)
  
  
  # remove region from matrix
  if(length(hasZero) > 0) {
    toRemove <- c(toRemove, row.names(minorCopyMatrix)[i])
  }
  
  i <- i + 1
  
}



# is only necessary one of these 
# one of the matrices, cause the number of rows with 0 is the same for both minor and major since it goes to next iteration without entering the condition

minorCopyMatrix <- minorCopyMatrix[!(row.names(minorCopyMatrix) %in% toRemove), ]
majorCopyMatrix <- majorCopyMatrix[!(row.names(majorCopyMatrix) %in% toRemove), ]



minorEpsilon <- minorEpsilon[!(row.names(minorEpsilon) %in% toRemove), ]
majorEpsilon <- majorEpsilon[!(row.names(majorEpsilon) %in% toRemove), ]



# remove CNAs with Wm ~= WM (these never happens in this case, all Wm and WM are different)

# Delete these CNAs regions with samples without copy number info (and that do not intersect a CNA) from matrix C
cMatrix <- cMatrix[!(row.names(cMatrix) %in% toRemove) , ]




# Delete CNAs that, after regions without copy number in some samples were deleted, do not intersect any region

toRemove <- c()

for (i in 1: as.numeric(ncol(cMatrix))) {
  
  # evaluate if does not have any 1 in that column
  r_content <- lapply(cMatrix[ , i], function(x) 1 %in% x) 
  
  hasZero <- which (r_content == TRUE)
  
  
  # remove CNA from matrix
  if(length(hasZero) == 0) {
    toRemove <- c(toRemove, colnames(cMatrix)[i])
  }
  
  i <- i + 1
  
}


cMatrix <- cMatrix[ , !(colnames(cMatrix) %in% toRemove)]





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------





#####  Calculate unions of nested CNAs resulting in the new CNA  #####


# For each mouse separately 


for (i in 1 : as.numeric(length(miceCnasInfo))) {
  
  
  mice_chrms_samples <- get(paste("sampleChrmsList", mice[i], sep = ""), envir = .GlobalEnv)
  
  t <- 1
  
  
  assign(paste("unionCNAS", mice[i], sep = ""), list())
  
  
  auxSamples <- NULL
  
  
  for(j in 1 : as.numeric(length(names(mice_chrms_samples)))) {
    
    
    current_sample <- names(mice_chrms_samples)[j]
    
    chrms_to_iterate <- names(mice_chrms_samples[[current_sample]])
    
    
    new_cnas <- list()
    
    # mice_union_cnas <- get(paste("unionCNAS", mice[i], sep = ""), envir = .GlobalEnv)
    
    # previous line would get the previous value of the variable, which was the union of the cnas of the previous sample
    # we would make modifications based on that, which would give incorrect results
    
    # because the each position holds the 4 samples for each mouse
    mice_union_cnas <- miceCnasInfo[[i]] [miceCnasInfo[[i]]$`SAMPLE ID` == current_sample,  ]
    
    
    
    for (g in 1 : as.numeric(length(chrms_to_iterate))) {
      
      
      current_chrm <- mice_chrms_samples[[current_sample]][[chrms_to_iterate[g]]]
      
      
      
      for(h in 1 : as.numeric(length(current_chrm))) {
        
        
        if(length(current_chrm[[h]]) < 2) {
          
          next
        }
        
        
        # calculate new CNA
        else {
          
          chrm_cnas <- current_chrm[[h]]
          
          
          # substring of start CNA and of end of CNA
          cna_size <- as.numeric(length(chrm_cnas))
          
          
          start_pos <- as.numeric(substr(chrm_cnas[1], unlist(gregexpr(":", chrm_cnas[1]))[1] + 1,  unlist(gregexpr(":", chrm_cnas[1]))[2] - 1))
          
          end_pos <- as.numeric(substr(chrm_cnas[length(chrm_cnas)],
                                       unlist(gregexpr(":", chrm_cnas[length(chrm_cnas)]))[2] + 1,  nchar(chrm_cnas[length(chrm_cnas)]) ))
          
          
          cna_range <- end_pos - start_pos
          
          
          # because for all mice together each sample is separated, and here all 4 samples of a mouse are together
          cna_sample <- miceCnasInfo[[i]] [miceCnasInfo[[i]]$`SAMPLE ID` == current_sample,    ]
          
          
          cna_sample <- cna_sample [cna_sample$`CHROM` == tolower(chrms_to_iterate[g]), ]
          
          
          
          cna_sample$`COPY ALTERATION ID` <- sapply(cna_sample$`COPY ALTERATION ID`, get_identifier)
          
          
          # to put the indexes from 1 to nrow(sample_chrm_region)
          row.names(cna_sample) <- NULL
          
          
          
          first_cna <- cna_sample[cna_sample$`COPY ALTERATION ID` == chrm_cnas[1], ]
          
          
          cna_1 <- first_cna$`MAJOR COPY`
          
          cna_2 <- first_cna$`MINOR COPY`
          
          
          
          mice_union_cnas <- mice_union_cnas [! (sapply(mice_union_cnas$`COPY ALTERATION ID`, get_identifier) %in% first_cna), ]
          
          
          
          # does not start in one cause we want to compare the copy number to see if it is nested
          for (k in 2 : cna_size) {
            
            current_cna <- cna_sample[cna_sample$`COPY ALTERATION ID` == chrm_cnas[k], ]
            
            
            if (current_cna$`MAJOR COPY` == cna_1 & current_cna$`MINOR COPY` == cna_2) {
              
              mice_union_cnas <-  mice_union_cnas [! (sapply( mice_union_cnas$`COPY ALTERATION ID`, get_identifier) %in% current_cna), ]
              
            }
            
            # else, we do not remove, is a nested
            
          }
          
          
          
          
          new_cnas[[t]] <- c('COPY ALTERATION ID' = paste(paste(paste(chrms_to_iterate[g], start_pos, sep = ":"), end_pos, sep = ":"), current_sample, sep = ":"),
                             'SAMPLE ID' = current_sample , 'CHROM' = tolower(chrms_to_iterate[g]), 'START POS' = start_pos, 'END POS' = end_pos,
                             # this is just to have the same columns as cnasInfo
                             'COPY NUMBER' = 0, 'MAJOR COPY' = 0, 'MINOR COPY' = 0, 'CI LOWER' = 0, 'CI UPPER' = 0, 'BAF' = 0, 'LOG2' = 0)
          
          
          
          t <- t + 1
          
        }
        
        
        
      }
      
      
      
    }
    
    
    df <- as.data.frame(do.call(rbind, new_cnas))
    
    mice_union_cnas <- rbind(mice_union_cnas, df)
    
    
    
    mice_union_cnas$'CHRM NUMBER' <- sapply( mice_union_cnas$`COPY ALTERATION ID`, function(x) as.numeric(substr(x,
                                                                                                                 unlist(gregexpr("R", x))[1] + 1,  unlist(gregexpr(":", x))[1] - 1)))
    
    
    mice_union_cnas <- mice_union_cnas [order( mice_union_cnas$`SAMPLE ID`, as.numeric( mice_union_cnas$`CHRM NUMBER`), as.numeric( mice_union_cnas$`START POS`)) , ]
    
    
    
    row.names(mice_union_cnas) <- NULL
    
    
    auxSamples <- rbind(auxSamples, mice_union_cnas)
    
  }
  
  
 assign(paste("unionCNAS", mice[i], sep = ""), auxSamples)
  
  
}





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------





# For each mouse separately


#####  Identify CNA regions  #####



for (i in 1 : as.numeric(length(miceCnasInfo))) {
  
  
  mice_union_cnas <- get(paste("unionCNAS", mice[i], sep = ""), envir = .GlobalEnv) 
  
  mice_union_cnas$`COPY ALTERATION ID` <- sapply(mice_union_cnas$`COPY ALTERATION ID`, get_identifier)
  
  
  mice_cna_regions <- c()
  
  chrm_number <- ""
  
  
  for(j in 1: as.numeric(length(unique(mice_union_cnas$CHROM)))) {
    
    current_chrm <- mice_union_cnas [mice_union_cnas $CHROM == unique(mice_union_cnas$CHROM)[j],  ]
    
    current_chrm <- current_chrm[!duplicated(current_chrm$`COPY ALTERATION ID`), ]
    
    
    chrm_number <- current_chrm[1, ]$`CHROM`  
    
    
    
    if( length(unique(current_chrm$`COPY ALTERATION ID`)) == 1 ) {
      
      mice_cna_regions <- c(mice_cna_regions, paste("R", unique(current_chrm$`COPY ALTERATION ID`), sep = ":"))
      
      
      
    } else {
      
      max_iterations <- length(current_chrm)
      
      
      
      for(g in 1 : max_iterations) {
        
        region_start <- min(as.numeric(current_chrm$`START POS`))
        
        region_end <- max(as.numeric(current_chrm[current_chrm$`START POS` == region_start, ]$`END POS`))
        
        
        next_intersection <- current_chrm
        
        
        
        while ( nrow(next_intersection) != 0 ) {
          
          min_start_pos <- min(as.numeric(next_intersection$`START POS`))
          
          max_connection <- max(as.numeric(next_intersection[as.numeric(next_intersection$`START POS`) == min_start_pos, ]$`END POS`))
          
          
          current_chrm <- current_chrm [ ! (as.numeric(current_chrm$`START POS`) >= min_start_pos & as.numeric(current_chrm$`END POS`) <= max_connection)  ,  ]
          
          next_intersection <- current_chrm[as.numeric(current_chrm$`START POS`) < max_connection, ]
          
        }
        
        
        
        if(nrow(next_intersection) == 0) { 
          
          region_end <- max_connection
          
          mice_cna_regions <- c(mice_cna_regions, paste("R", paste(toupper(chrm_number), paste(toString(region_start), toString(region_end), sep = ":"), sep = ":"), sep = ":"))
          
          
          if(nrow(current_chrm) == 0) {
            
            break
          }
          
          
        }
        
        
      }
      
      
    }
    
    
  }
  
  assign(paste("unionCNAS", mice[i], sep = ""), mice_union_cnas)
  
  
  assign(paste("miceCNARegions", mice[i], sep = ""), mice_cna_regions)
  
}






# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------






#####  Construct C Matrix  #####



# For each mice separately 


for (i in 1 : as.numeric(length(miceCnasInfo))) {
  
  
  mice_union_cnas <- get(paste("unionCNAS", mice[i], sep = ""), envir = .GlobalEnv) 
  
  
  
  mice_unique_cnas <- mice_union_cnas [!duplicated(mice_union_cnas$`COPY ALTERATION ID`), ]
  
  mice_unique_cnas <- mice_unique_cnas [order(as.numeric(mice_unique_cnas$`CHRM NUMBER`), 
                                              as.numeric(mice_unique_cnas$`START POS`), as.numeric(mice_unique_cnas$`END POS`)),  ]
  
  
  assign(paste("miceUniqueCNAS", mice[i], sep = ""), mice_unique_cnas)
  
  
  
  
  mice_cna_regions <- get(paste("miceCNARegions", mice[i], sep = ""), envir = .GlobalEnv) 
  
  
  
  mice_cna_regions <- mice_cna_regions [ order( as.numeric(sapply(mice_cna_regions, function(x) as.numeric(substr(x,
                                                                                                                  unlist(gregexpr("R", x))[2] + 1,  unlist(gregexpr(":", x))[2] - 1)))),                     
                                                as.numeric(sapply(mice_cna_regions, function(x) as.numeric(substr(x,
                                                                                                                  unlist(gregexpr(":", x))[2] + 1,  unlist(gregexpr(":", x))[3] - 1)))),
                                                as.numeric(sapply(mice_cna_regions, function(x) as.numeric(substr(x,
                                                                                                                  unlist(gregexpr(":", x))[3] + 1,  nchar(x) - 1)))) ) ]
  
  
  assign(paste("miceCNARegions", mice[i], sep = ""), mice_cna_regions)
  
  
  
  
  mice_unique_cnas <- get(paste("miceUniqueCNAS", mice[i], sep = ""), envir = .GlobalEnv) 
  
  mice_cna_regions <- get(paste("miceCNARegions", mice[i], sep = ""), envir = .GlobalEnv) 
  
  
  
  
  
  assign(paste("cMatrix", mice[i], sep = ""), matrix(0, nrow = as.numeric(length(mice_cna_regions)), 
                                                     ncol = as.numeric(length(unique(mice_unique_cnas$`COPY ALTERATION ID`))), byrow = TRUE, 
                                                     dimnames = list(mice_cna_regions, paste("CNA", unique(mice_unique_cnas$`COPY ALTERATION ID`), sep = ":"))))       
  
  
  c_matrix <- get(paste("cMatrix", mice[i], sep = ""), envir = .GlobalEnv) 
  
  


  for(j in 1 : as.numeric(length(unique(mice_unique_cnas$CHROM)))) {

    current_chrm <- mice_unique_cnas [mice_unique_cnas$CHROM == unique(mice_unique_cnas$CHROM)[j], ]

    row.names(current_chrm) <- NULL

    chrm_regions <- mice_cna_regions[ sapply(mice_cna_regions, function (x)
      substr(x, unlist(gregexpr(":", x))[1] + 1, unlist(gregexpr(":", x))[2] - 1) == toupper(unique(mice_unique_cnas$CHROM)[j]))]



    # iterate CNA regions of a chromosome
    for (g in 1 : as.numeric(length(chrm_regions))) {

      current_region <- chrm_regions[g]

      region_start <- as.numeric(substr(current_region, unlist(gregexpr(":", current_region))[2] + 1, unlist(gregexpr(":", current_region))[3] - 1))

      region_end <- as.numeric(substr(current_region, unlist(gregexpr(":", current_region))[3] + 1, nchar(current_region)))



      # iterate CNAs in that chromosome
      for (h in 1 : nrow(current_chrm) ) {

        current_cna <- current_chrm[h,  ]

        start_pos <- as.numeric(current_cna$`START POS`)

        end_pos <- as.numeric(current_cna$`END POS`)


        # The CNA is out of the region limits, so it does not intersect it
        if ((region_start > end_pos - 0.05 & region_start > start_pos + 0.05) |
            (region_end < start_pos + 0.05 & region_end < end_pos - 0.05)) {

          h <- h + 1
          next



        # When the CNA start and end is the same as the region (CNA stays completely within the bounds of the region, in the middle)
        } else if ((start_pos + 0.05 > region_start & end_pos - 0.05 < region_end) |

                   # When the CNA intersects a part of the region, and it is to its left or right
                   (end_pos - 0.05 > region_start & start_pos + 0.05 < region_end)) {



          cna_name <- paste("CNA", current_cna$`COPY ALTERATION ID`, sep = ":")


          c_matrix [current_region, cna_name] <- 1


        }


      }


    }


  }



  assign(paste("cMatrix", mice[i], sep = ""), c_matrix)

  
}





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------





#####  Matrices Wm, WM and Epsilon Wm and WM  #####


# For each mice separately


# Formulas to use: 
# BAF = Wm / (WM + Wm)
# Log2 (Depth.ratio) = (WM + Wm)/2



for (i in 1 : as.numeric(length(miceCnasInfo))) {
  
  mice_snvs <- miceGenomeInfo[[i]]
  
  samples_id <- unique(mice_snvs$`SAMPLE ID`)
  
  
  
  c_matrix <- get(paste("cMatrix", mice[i], sep = ""), envir = .GlobalEnv) 
  
  
  
  
  # Wm (then assign minorCopyMatrix to Wm in the algorithm parameters)
  assign(paste("minorCopyMatrix", mice[i], sep = ""), matrix(0, nrow = as.numeric(length(row.names(c_matrix))), 
                                                             ncol = as.numeric(length(samples_id)), byrow = TRUE, 
                                                             dimnames = list(row.names(c_matrix), samples_id)))    
  
  
  # WM (then assign majorCopyMatrix to WM in the algorithm parameters)
  assign(paste("majorCopyMatrix", mice[i], sep = ""), matrix(0, nrow = as.numeric(length(row.names(c_matrix))), 
                                                             ncol = as.numeric(length(samples_id)), byrow = TRUE, 
                                                             dimnames = list(row.names(c_matrix), samples_id)))    
  
  
  
  # epsilon Wm
  assign(paste("minorEpsilon", mice[i], sep = ""), matrix(0, nrow = as.numeric(length(row.names(c_matrix))), 
                                                          ncol = as.numeric(length(samples_id)), byrow = TRUE, 
                                                          dimnames = list(row.names(c_matrix), samples_id)))    
  
  
  
  # epsilon WM
  assign(paste("majorEpsilon", mice[i], sep = ""), matrix(0, nrow = as.numeric(length(row.names(c_matrix))), 
                                                          ncol = as.numeric(length(samples_id)), byrow = TRUE, 
                                                          dimnames = list(row.names(c_matrix), samples_id)))    
  
  
  
  
  
  # Non-CNAs (copy number == 2) and CNAs (copy number != 2)
  miceCopyNumber <- miceCopyNumberInfo[[i]]
  
  
  miceCopyNumber$`COPY ALTERATION ID` <- sapply(miceCopyNumber$`COPY ALTERATION ID`, get_identifier)
  
  
  
  
  minorCopy <- get(paste("minorCopyMatrix", mice[i], sep = ""), envir = .GlobalEnv) 
  
  majorCopy <- get(paste("majorCopyMatrix", mice[i], sep = ""), envir = .GlobalEnv) 
  
  minorEpsilonMice <- get(paste("minorEpsilon", mice[i], sep = ""), envir = .GlobalEnv) 
  
  majorEpsilonMice <- get(paste("majorEpsilon", mice[i], sep = ""), envir = .GlobalEnv) 
  
  
  
  
  for (j in 1 : as.numeric(length(row.names(c_matrix)))) {
    
    
    currentCNARegion <- row.names(c_matrix)[j]
    
    
    region_start <- as.numeric(substr(currentCNARegion, unlist(gregexpr(":", currentCNARegion))[2] + 1, unlist(gregexpr(":", currentCNARegion))[3] - 1 ))
    
    region_end <- as.numeric(substr(currentCNARegion, unlist(gregexpr(":", currentCNARegion))[3] + 1, as.numeric(nchar(currentCNARegion))))
    
    
    currentChr <- substr(currentCNARegion, unlist(gregexpr(":", currentCNARegion))[1] + 1, unlist(gregexpr(":", currentCNARegion))[2] - 1 )
    
    
    currentChrCNAS <- miceCopyNumber[miceCopyNumber$CHROM == tolower(currentChr), ]
    
    
    
    
    for (g in 1 : as.numeric(length(samples_id))) {
      
      sample <- samples_id[g]
      
      currentChrSample <- currentChrCNAS[currentChrCNAS$`SAMPLE ID` == sample, ]
      
      
      toCompareWith <- list()
      
      
      
      # When the CNA (or copy number info) start and end is the same as the region (CNA stays completely within the bounds of the region, in the middle) 
      # Or when the CNA (or copy number info) intersects a part of the region, and it is to its left or right 
      toCompareWith <- currentChrSample[ (as.numeric(currentChrSample$`START POS`) + 0.05 > region_start & as.numeric(currentChrSample$`END POS`) - 0.05 < region_end) |
                                           (as.numeric(currentChrSample$`END POS`) - 0.05 > region_start & as.numeric(currentChrSample$`START POS`) + 0.05 < region_end) , ]
      
      
      
      
      # just for verification purposes, this should not happen
      if(nrow(toCompareWith) == 0){
        g <- g + 1
        next
      }
      
      
      
      if(nrow(toCompareWith) == 1) {
        
        # BAF = Wm / (WM + Wm)
        # Depth.ratio = (WM + Wm)/2
        
        baf <- toCompareWith$BAF
        
        # depth_ratio = 2 ^ log2 column
        depth_ratio <- 2^(as.numeric(toCompareWith$LOG2))
        
        
        # Wm
        minor_copy <- 2 * baf * depth_ratio
        
        # WM
        major_copy <- 2 * depth_ratio - minor_copy
        
        
        if (baf > 0.5) { # major_copy is going to be < minor_copy
          
          minorCopy [currentCNARegion, toCompareWith$`SAMPLE ID`] <- major_copy 
          majorCopy [currentCNARegion, toCompareWith$`SAMPLE ID`] <- minor_copy
          
          
        } else {
          
          # add minor_copy to position [i, j] of minorCopyMatrix
          # add major_copy to position [i, j] of majorCopyMatrix
          minorCopy [currentCNARegion, toCompareWith$`SAMPLE ID`] <- minor_copy
          majorCopy [currentCNARegion, toCompareWith$`SAMPLE ID`] <- major_copy
          
        }
        
        
        
        # SE = (upper limit - lower limit) / 3.92.
        
        # Confidence interval of the segment mean (--ci), estimated by bootstrap (100 resamplings) of the bin-levellog2 ratio values 
        # within the segment.  The upper and lower bounds are output as separate columns ci_lo and ci_hi
        # That is why it is 2 ^ the ci_lo and ci_hi
        
        
        ci_lower <- 2^toCompareWith$`CI LOWER`
        
        ci_upper <- 2^toCompareWith$`CI UPPER`
        
        
        epsilon <- (ci_upper - ci_lower) / 3.92        
        
        
        
        minorEpsilonMice [currentCNARegion, toCompareWith$`SAMPLE ID`]  <- epsilon
        majorEpsilonMice [currentCNARegion, toCompareWith$`SAMPLE ID`]  <- epsilon
        
        
        
      } else { # cna region intersects different CNAs
        
        
        # weighted average 
        
        ci_lower <- ci_upper <- baf <- depth_ratio <- 0 
        
        
        # Calculate which intersection is bigger, with which CNA info
        for (h in 1 : as.numeric(nrow(toCompareWith))) { 
          
          cna_intersecting <- toCompareWith[h]
          
          # we are going use the info that is there, even if the green intersects a part of that 
          # inside a sample none of the cnas in that sample intersects one another
          
          
          # CNA stays in between region:
          if(toCompareWith[h]$`START POS` + 0.05 > region_start & toCompareWith[h]$`END POS` - 0.05 < region_end) { 
            
            cna_range <- as.numeric(cna_intersecting$`END POS`) - as.numeric(cna_intersecting$`START POS`)
            
            # percentage of the area covered by that CNA
            p <- cna_range / (region_end - region_start)
            
          }
          
          
          # When the CNA intersects a part of the region, and it is to its left or right 
          else {
            
            if (toCompareWith[h]$`END POS` - 0.05 > region_start & toCompareWith[h]$`START POS` + 0.05 < region_end) {
              
              
              # Region stays to the left of the CNA (or copy number info)
              if(toCompareWith[h]$`START POS` + 0.05 > region_start) {
                
                cna_range <- region_end - as.numeric(cna_intersecting$`START POS`)
                
                p <- cna_range / (region_end - region_start)
                
                
              } 
              
              
              # Region stays to the right of the CNA (or copy number info)
              else if (toCompareWith[h]$`END POS` - 0.05 < region_end){
                
                cna_range <- as.numeric(cna_intersecting$`END POS`) - region_start
                
                
                p <- cna_range / (region_end - region_start)
                
                
              }
              
              
              # CNA is bigger than the region, and so the region stays contained inside it
              else {
                
                p <- 1
                
                
              }
              
              
            }
            
            
          }
          
          
          
          baf <- baf + (p * as.numeric(cna_intersecting$`BAF`)) 
          
          depth_ratio <- depth_ratio + (p * (2^as.numeric(cna_intersecting$`LOG2`) ))  
          
          
          ci_lower <- ci_lower + (p * (2^as.numeric(cna_intersecting$`CI LOWER`) ))          
          
          ci_upper <- ci_upper + (p * (2^as.numeric(cna_intersecting$`CI UPPER`) ))     
          
          
        }
        
        
        
        # weigthed average
        
        baf <- baf / nrow(toCompareWith)
        
        depth_ratio <- depth_ratio / nrow(toCompareWith)
        
        
        ci_lower <- ci_lower  / nrow(toCompareWith)
        
        ci_upper <- ci_upper / nrow(toCompareWith)
        
        
        
        # Wm
        minor_copy <- 2 * baf * depth_ratio
        
        # WM
        major_copy <- 2 * depth_ratio - minor_copy
        
        
        
        if (baf > 0.5) { # major_copy is going to be < minor_copy
          
          minorCopy [currentCNARegion, cna_intersecting$`SAMPLE ID`] <- major_copy 
          majorCopy [currentCNARegion, cna_intersecting$`SAMPLE ID`] <- minor_copy
          
          
          
        } else {
          
          # add minor_copy to position [i, j] of minorCopyMatrix
          # add major_copy to position [i, j] of majorCopyMatrix
          minorCopy [currentCNARegion, cna_intersecting$`SAMPLE ID`] <- minor_copy 
          majorCopy [currentCNARegion, cna_intersecting$`SAMPLE ID`] <- major_copy
          
        }
        
        
        
        # SE = (upper limit - lower limit) / 3.92.
        
        # Confidence interval of the segment mean (--ci), estimated by bootstrap (100 resamplings) of the bin-levellog2 ratio values 
        # within the segment.  The upper and lower bounds are output as separate columns ci_lo and ci_hi
        # That is why it is 2 ^ the ci_lo and ci_hi
        
        epsilon <- (ci_upper - ci_lower) / 3.92        
        
        
        minorEpsilonMice[currentCNARegion, cna_intersecting$`SAMPLE ID`]  <- epsilon
        majorEpsilonMice[currentCNARegion, cna_intersecting$`SAMPLE ID`]  <- epsilon
        
        
      }
      
    }
    
  }  
  
  
  assign(paste("minorCopyMatrix", mice[i], sep = ""), minorCopy)
  
  assign(paste("majorCopyMatrix", mice[i], sep = ""), majorCopy)
  
  assign(paste("minorEpsilon", mice[i], sep = ""), minorEpsilonMice)
  
  assign(paste("majorEpsilon", mice[i], sep = ""), majorEpsilonMice)
  
  
  
  
  # After this check Wm and WM and see if is better to remove some CNA regions. Then, remove those CNA regions from the C matrix, and use only those CNA regions for Y, 
  # but before doing Y also verify if some SNVs must be removed, using a heatmap, which also means the SNVs removed need to be removed from the R and X matrices
  
  miceToRemove <- c()
  
  minor_matrix <- get(paste("minorCopyMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  major_matrix <- get(paste("majorCopyMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  
  
  minor_epsilon <- get(paste("minorEpsilon", mice[i], sep = ""), envir = .GlobalEnv)
  major_epsilon <- get(paste("majorEpsilon", mice[i], sep = ""), envir = .GlobalEnv)
  
  
  
  for(k in 1: as.numeric(nrow(minor_matrix))) {
    
    r_content <- lapply(minor_matrix[k, ], function(x) 0 %in% x) 
    
    hasZero <- which (r_content == TRUE)
    
    
    if(length(hasZero) > 0) {
      miceToRemove <- c(miceToRemove, row.names(minor_matrix)[k])
    }
    
  }
  
  
  minor_matrix <- minor_matrix[!(row.names(minor_matrix) %in% miceToRemove), ]
  major_matrix <- major_matrix[!(row.names(major_matrix) %in% miceToRemove), ]
  
  
  assign(paste("minorCopyMatrix", mice[i], sep = ""), minor_matrix)
  assign(paste("majorCopyMatrix", mice[i], sep = ""), major_matrix)
  
  
  
  
  minor_epsilon <- minor_epsilon[!(row.names(minor_epsilon) %in% miceToRemove), ]
  major_epsilon <- major_epsilon[!(row.names(major_epsilon) %in% miceToRemove), ]
  
  
  assign(paste("minorEpsilon", mice[i], sep = ""), minor_epsilon)
  assign(paste("majorEpsilon", mice[i], sep = ""), major_epsilon)
  
  
  
  # Delete these CNAs regions with samples without copy number info (and that do not intersect a CNA) from matrix C
  c_m <- get(paste("cMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  
  c_m <- c_m[!(row.names(c_m) %in% miceToRemove) , ]
  
  assign(paste("cMatrix", mice[i], sep = ""), c_m)
  
  
  
  
  
  # Delete CNAs that, after regions without copy number in some samples, were deleted do not intersect any region
  
  miceToRemove <- c()
  
  c_m <- get(paste("cMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  
  
  for (k in 1: as.numeric(ncol(c_m))) {
    
    # evaluate if does not have any 1 in that column
    r_content <- lapply(c_m[ , k], function(x) 1 %in% x) 
    
    hasZero <- which (r_content == TRUE)
    
    
    # remove CNA from matrix
    if(length(hasZero) == 0) {
      miceToRemove <- c(miceToRemove, colnames(c_m)[k])
    }
    
    k <- k + 1
    
  }
  
  
  c_m <- c_m[ , !(colnames(c_m) %in% miceToRemove)]
  
  
  assign(paste("cMatrix", mice[i], sep = ""), c_m)
  
}


rm(auxSamples, c_m, c_matrix, chr_r, chrm_r, chrm_region, chrmRegionCNAS, chrmRegionsCNASM49, chrmRegionsCNASM55, chrmRegionsCNASM61, chrmRegionsCNASM62,
   cna_intersecting, cna_sample, current_chrm, current_cna, current_CNA, current_mice, currentChrCNAS, currentChrSample, currentMiceCNAS, df, first_cna,
   major_epsilon, major_matrix, majorEpsilonMice, majorCopy, mice_chrms_list, mice_chrms_samples, mice_snvs, miceCopyNumber, minor_epsilon, minor_matrix,
   minorCopy, minorEpsilonMice, new_cnas, next_CNA, next_intersection, r_content, s_chrms_list, sample_chrm_region, sample_chrms_list,
   sampleChrmsListM49, sampleChrmsListM55, sampleChrmsListM61, sampleChrmsListM62, toCompareWith, baf, chrm_cnas, chrm_number, chrm_regions, chrms_samples,
   chrms_to_iterate, chromosomes_identifiers, ci_lower, ci_upper, cna_1, cna_2, cna_name, cna_range, cna_regions, cna_size, current_region, current_sample,
   currentChr, currentCNARegion, depth_ratio, end_pos, epsilon, equal_copy, hasZero, is_next, major_copy, max_connection, max_iterations, mice_samples,
   miceToRemove, min_start_pos, minor_copy, not_nested, p, previousToMerge, region_end, region_start, sample, samples_id, start_CNA, start_pos,
   t, toMerge, toRemove, mut_identifiers, mice_cna_regions, miceCNARegionsM49, miceCNARegionsM55, miceCNARegionsM61, miceCNARegionsM62, unique_cnas,
   union_cnas, unionCNASM49, unionCNASM55, unionCNASM61, unionCNASM62, miceUniqueCNASM49, miceUniqueCNASM55, miceUniqueCNASM61,
   miceUniqueCNASM62, mice_union_cnas, mice_unique_cnas)






# ----------------------------------------------------------------------------------------------------------------------------------------------





#####  Heatmap of mutations quality in each sample  #####

# For all mice together



# remove SNAs with similar VAF, remove them from R and X matrices

# Then, calculate Y with the CNAs regions and SNAs chosen


mut_identifiers <- unique(allSamples$`MUTATION ID`)

samples_id <- unique(allSamples$`SAMPLE ID`)


vafMatrix <- matrix(0, nrow = as.numeric(length(mut_identifiers)), ncol = as.numeric(length(samples_id)), 
                    byrow = TRUE, dimnames = list(mut_identifiers, samples_id))



for (i in 1 : length(mut_identifiers)) {
  
  for (j in 1 : length(samples_id)) {
    
    current_SNA <- allSamples[allSamples$`MUTATION ID` == mut_identifiers[i] & allSamples$`SAMPLE ID` == samples_id[j] , ]
    
    vafMatrix[i, j] <- current_SNA$`VAF`
    
  }
  
  
}



# the pics are saved to the current directory defined


# On linux
#setwd("~/Desktop/Dados/CNVS/canopyImages")


# On windows
setwd("C:/Users/david/Desktop/Dados/CNVS/canopyImages")


# install.packages("pheatmap")

library("pheatmap")


# The default behavior of the function includes the hierarchical clustering of both rows and columns, 
# in which we can observe similar players and stats types in close positions.

# pheatmap(vafMatrix, main = "pheatmap default")


# Not all SNAs seem to have a distinct pattern between samples


# cluster_cols = FALSE would put samples in order in columns, but the most similar samples would not be closer to each other
snvsHeatmap <- pheatmap(vafMatrix, cutree_rows = 6, main = "pheatmap row cut")

dev.copy(jpeg, filename="snvsHeatmap.jpg", width = 2000, height = 10000)

dev.off()



# See in how many samples each mutation has a VAF of 0
samples_count <- sqldf("SELECT `MUTATION ID`, COUNT(`VAF`) AS 'SAMPLES COUNT' from allSamples 
                        WHERE `VAF` > 0
                        GROUP BY `MUTATION ID`
                        ORDER BY COUNT(*) DESC")


# Each mutation needs to be present in at least 6 samples
toRemove <- samples_count[samples_count$`SAMPLES COUNT` < 6 , ]

# Keep all the mutations that are not part of the ones to be removed
allSamples <- allSamples[ !(allSamples$`MUTATION ID` %in% toRemove$`MUTATION ID`), ]



# VAF matrix after removing the SNVs that have a more constant pattern across all samples


# For the case that, when the SNVs were removed, some sample did no longer have any mutation (unique(allSamples$`SAMPLE ID`))
vafMatrix <- matrix(0, nrow = as.numeric(length(unique(allSamples$`MUTATION ID`))), ncol = as.numeric(length(unique(allSamples$`SAMPLE ID`))), 
                    byrow = TRUE, dimnames = list(unique(allSamples$`MUTATION ID`), unique(allSamples$`SAMPLE ID`)))


for (i in 1 : length(unique(allSamples$`MUTATION ID`))) {
  
  for (j in 1 : length(unique(allSamples$`SAMPLE ID`))) {
    
    current_SNA <- allSamples[allSamples$`MUTATION ID` == unique(allSamples$`MUTATION ID`)[i] & allSamples$`SAMPLE ID` == samples_id[j] , ]
    
    vafMatrix[i, j] <- current_SNA$`VAF`
    
  }
  
}


snvsHeatmap <- pheatmap(vafMatrix, cutree_rows = 6, main = "pheatmap row cut")

dev.copy(jpeg, filename="snvsHeatmapAfterRemoving.jpg", width = 2000, height = 10000)

dev.off()



# remove SNAs removed from vaf matrix in R and X matrices:
rMatrix <- rMatrix[row.names(rMatrix) %in% paste("SNA", row.names(vafMatrix), sep = ":") , ]
xMatrix <- xMatrix[row.names(xMatrix) %in% paste("SNA", row.names(vafMatrix), sep = ":") , ]






#####  Matrix Y (CNA regions and SNAs)  #####


# for all mice together


# add non-cna region column to Y matrix
overlap_cols <- c("non-cna_region", row.names(cMatrix))



#####  Intersect SNVs and CNAS genome regions and build Y matrix of overlap of SNVs and CNAS  #####

yMatrix <- matrix(0, nrow = as.numeric(length(row.names(vafMatrix))), ncol = as.numeric(length(overlap_cols)),
                  byrow = TRUE, dimnames = list(paste("SNA", row.names(vafMatrix), sep = ":") , overlap_cols))



for (i in 1 : as.numeric(length(row.names(vafMatrix)))) {
  
  cSNA <- row.names(yMatrix)[i]
  
  snaPos <- as.numeric(substr(cSNA, unlist(gregexpr(":", cSNA))[2] + 1, unlist(gregexpr(":", cSNA))[3] - 1))
  
  
  currentCharSNA <- substr(cSNA, unlist(gregexpr(":", cSNA))[1] + 1, unlist(gregexpr(":", cSNA))[2] - 1)
  
  # get all CNA regions in that chromosome
  chrRegions <- row.names(cMatrix)[ which(lapply(row.names(cMatrix), function(x) identical(substr(x, as.numeric(unlist(gregexpr(":", x))[1]) + 1,
                                                                                                  as.numeric(unlist(gregexpr(":", x))[2]) - 1) , currentCharSNA)) == TRUE) ]
  
  
  # if the SNV is intersected by any CNA region
  if(as.numeric(length(chrRegions > 0))) {
    
    for (j in 1 : length(chrRegions)){
      
      region_start <- as.numeric(substr(chrRegions[j], unlist(gregexpr(":", chrRegions[j]))[2] + 1, unlist(gregexpr(":", chrRegions[j]))[3] - 1))
      
      region_end <- as.numeric(substr(chrRegions[j], unlist(gregexpr(":", chrRegions[j]))[3] + 1, nchar(chrRegions[j])))
      
      
      # The SNV is only one position, and so to intersect the CNA region must be between its position ranges
      
      if(snaPos + 0.05 > region_start & snaPos + 0.05 < region_end){
        
        yMatrix[cSNA, chrRegions[j]] <- 1
        
      }
      
      j <- j + 1
      
    }
    
  }
  
  
  r_content <- lapply(yMatrix[cSNA, ], function(x) identical(yMatrix[cSNA, x], numeric(1)))
  
  hasSNA <- which (r_content == TRUE)
  
  
  # If the SNV is not intersected by any CNA region (does not have any [SNV, CNA region] = 1), indicate that is not overlapped by a CNA region
  if(length(hasSNA) == 0) {
    yMatrix[cSNA, overlap_cols[1]] <- 1
  }
  
  
  i <- i + 1
  
}

mode(yMatrix) <- "integer"




# ----------------------------------------------------------------------------------------------------------------------------------------------





#####  Heatmap of mutations quality in each sample  #####

# For each mouse separately 



for(i in 1 : as.numeric(length(miceGenomeInfo))) {
  
  currentMice <- miceGenomeInfo[[i]]
  
  
  mut_identifiers <- unique(currentMice$`MUTATION ID`)
  
  samples_id <- unique(currentMice$`SAMPLE ID`)
  
  
  assign(paste("vafMatrix", mice[i], sep = ""), matrix(0, nrow = as.numeric(length(mut_identifiers)), ncol = as.numeric(length(samples_id)),
                                                       byrow = TRUE, dimnames = list(mut_identifiers, samples_id)))
  
  
  
  for(j in 1 : as.numeric(length(mut_identifiers))) {
    
    for(g in 1 : as.numeric(length(samples_id))) {
      
      current_SNA <- currentMice[currentMice$`MUTATION ID` == mut_identifiers[j] & currentMice$`SAMPLE ID` == samples_id[g] , ]
      
      
      vaf_matrix <- get(paste("vafMatrix", mice[i], sep = ""), envir = .GlobalEnv)
      
      vaf_matrix[j, g] <- current_SNA$`VAF`
      
      assign(paste("vafMatrix", mice[i], sep = ""), vaf_matrix)
      
    }
    
  }
  
  
  vaf_matrix <- get(paste("vafMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  
  
  miceSNVsHeatmap <- pheatmap(vaf_matrix, cutree_rows = 4, main = "pheatmap row cut")
  
  dev.copy(jpeg, filename = paste(paste("miceSNVsHeatmap", mice[i], sep = ""), ".jpg", sep = ""), width = 1000, height = 5000)
  
  dev.off()
  
  
  # Remove mutations that are only in one sample for each mouse
  
  miceSamplesCount <- sqldf("SELECT `MUTATION ID`, COUNT(`VAF`) AS 'SAMPLES COUNT' from currentMice
                        WHERE `VAF` > 0
                        GROUP BY `MUTATION ID`
                        ORDER BY COUNT(*) DESC")
  
  
  miceToRemove <- miceSamplesCount[miceSamplesCount$`SAMPLES COUNT` < 2, ]
  
  currentMice <- currentMice[ !(currentMice$`MUTATION ID` %in% miceToRemove$`MUTATION ID`), ]
  
  
  miceGenomeInfo[[i]] <- currentMice
  
  
  assign(paste("vafMatrix", mice[i], sep = ""), matrix(0, nrow = as.numeric(length(unique(currentMice$`MUTATION ID`))),
                                                       ncol = as.numeric(length(unique(currentMice$`SAMPLE ID`))), byrow = TRUE,
                                                       dimnames = list(unique(currentMice$`MUTATION ID`), unique(currentMice$`SAMPLE ID`))))
  
  
  # Vaf matrix after removing these mutations
  
  for (j in 1 : as.numeric(length(unique(currentMice$`MUTATION ID`)))) {
    
    for (g in 1 : as.numeric(length(unique(currentMice$`SAMPLE ID`)))) {
      
      current_SNA <- currentMice[currentMice$`MUTATION ID` == unique(currentMice$`MUTATION ID`)[j] &
                                   currentMice$`SAMPLE ID` == unique(currentMice$`SAMPLE ID`)[g] , ]
      
      
      vaf_matrix <- get(paste("vafMatrix", mice[i], sep = ""), envir = .GlobalEnv)
      
      vaf_matrix[j, g] <- current_SNA$`VAF`
      
      assign(paste("vafMatrix", mice[i], sep = ""), vaf_matrix)
      
    }
    
  }
  
  
  
  vaf_matrix <- get(paste("vafMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  
  
  miceSNVsHeatmap <- pheatmap(vaf_matrix, cutree_rows = 4, main = "pheatmap row cut")
  
  dev.copy(jpeg, filename = paste(paste("miceSNVsHeatmap", mice[i], sep = ""), "AfterRemoving.jpg", sep = ""), width = 1000, height = 5000)
  
  dev.off()
  
  
  
  # remove SNAs removed from vaf matrix in R and X matrices:
  
  r_m <- get(paste("rMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  
  r_m <- r_m[row.names(r_m) %in% paste("SNA",row.names(vaf_matrix), sep = ":") , ]
  
  assign(paste("rMatrix", mice[i], sep = ""), r_m)
  
  
  
  x_m <- get(paste("xMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  
  x_m <- x_m[row.names(x_m) %in% paste("SNA",row.names(vaf_matrix), sep = ":") , ]
  
  assign(paste("xMatrix", mice[i], sep = ""), x_m)
  
  
  
  
  
  
  #####  Matrix Y (CNA regions and SNAs)  #####
  
  
  # for each mouse separately
  
  
  # add non-CNA region column to Y matrix
  
  c_m <- get(paste("cMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  
  miceOverlapCols <- c("non-cna_region", row.names(c_m))
  
  
  assign(paste("yMatrix", mice[i], sep = ""), matrix(0, nrow = as.numeric(length(row.names(vaf_matrix))), ncol = as.numeric(length(miceOverlapCols)),
                                                     byrow = TRUE, dimnames = list(paste("SNA", row.names(vaf_matrix), sep = ":") , miceOverlapCols)))
  
  
  y_m <- get(paste("yMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  
  
  
  for (j in 1 : as.numeric(length(row.names(vaf_matrix)))) {
    
    micecSNA <- row.names(y_m)[j]
    
    
    snaPos <- as.numeric(substr(micecSNA, unlist(gregexpr(":", micecSNA))[2] + 1, unlist(gregexpr(":", micecSNA))[3] - 1))
    
    currentCharSNA <- substr(micecSNA, unlist(gregexpr(":", micecSNA))[1] + 1, unlist(gregexpr(":", micecSNA))[2] - 1)
    
    
    chrRegions <- row.names(c_m)[ which(lapply(row.names(c_m), function(x) identical(substr(x, as.numeric(unlist(gregexpr(":", x))[1]) + 1,
                                                                                            as.numeric(unlist(gregexpr(":", x))[2]) - 1) , currentCharSNA)) == TRUE) ]
    
    
    
    # Check if any of the CNA regions in that chromosome intersects the SNV
    if(as.numeric(length(chrRegions)) != 0 ) {
      
      for (g in 1 : as.numeric(length(chrRegions))) {
        
        region_start <- as.numeric(substr(chrRegions[g], unlist(gregexpr(":", chrRegions[g]))[2] + 1, unlist(gregexpr(":", chrRegions[g]))[3] - 1))
        
        region_end <- as.numeric(substr(chrRegions[g], unlist(gregexpr(":", chrRegions[g]))[3] + 1, nchar(chrRegions[g])))
        
        
        if(snaPos + 0.05 > region_start & snaPos + 0.05 < region_end){
          
          y_m[micecSNA, chrRegions[g]] <- 1
          
          assign(paste("yMatrix", mice[i], sep = ""), y_m)
          
        }
        
      }
      
    }
    
    
    y_m <- get(paste("yMatrix", mice[i], sep = ""), envir = .GlobalEnv)
    
    
    r_content <- lapply(y_m[micecSNA, ], function(x) identical(y_m[micecSNA, x], numeric(1)))
    
    hasSNA <- which (r_content == TRUE)
    
    
    # mark SNV has not being overlapped by any CNA region
    if(length(hasSNA) == 0) {
      y_m[micecSNA, overlap_cols[1]] <- 1
      
      assign(paste("yMatrix", mice[i], sep = ""), y_m)
    }
    
    
    y_m <- get(paste("yMatrix", mice[i], sep = ""), envir = .GlobalEnv)
    
    mode(y_m) <- "integer"
    
    assign(paste("yMatrix", mice[i], sep = ""), y_m)
    
  }
  
  
}


rm(x_m, r_m, y_m, c_m, current_SNA, currentMice, miceSamplesCount, miceSNVsHeatmap, r_content,
   samples_count, snvsHeatmap, toRemove, vaf_matrix, vafMatrix, vafMatrixM49, vafMatrixM55, vafMatrixM61,
   miceToRemove, vafMatrixM62, chrRegions, cSNA, currentCharSNA, hasSNA, micecSNA, miceOverlapCols, 
   mut_identifiers, overlap_cols, region_end, region_start, samples_id, snaPos)





# -----------------------------------------------------------------------------------------------------------------------------------------------





#####  Add matrices to final input list  #####


# For all mice together

projectname <- "MDA231ctrlsunit"


canopyInput <- list(rMatrix, xMatrix, cMatrix, majorCopyMatrix, minorCopyMatrix, majorEpsilon, minorEpsilon, yMatrix, projectname)

names(canopyInput) <- c("R", "X", "C", "WM", "Wm", "epsilonM", "epsilonm", "Y", "projectname")





# For each mouse separately 

for(i in 1 : as.numeric(length(miceGenomeInfo))) {
  
  currentMice <- miceGenomeInfo[[i]]
  
  
  # project name for each mouse
  
  if(substr(currentMice[1]$`SAMPLE ID`, 1, 3) == "M49" |  substr(currentMice[1]$`SAMPLE ID`, 1, 3) == "M55") {
    
    assign(paste("projectname", mice[i], sep = ""), paste("MDA231Ctrl", mice[i], sep = ""))
  }
  
  
  if(substr(currentMice[1]$`SAMPLE ID`, 1, 3) == "M61" |  substr(currentMice[1]$`SAMPLE ID`, 1, 3) == "M62") {
    
    assign(paste("projectname", mice[i], sep = ""), paste("MDA231Sunit", mice[i], sep = ""))
  }
  
  
  
  c_m <- get(paste("cMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  r_m <- get(paste("rMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  x_m <- get(paste("xMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  y_m <- get(paste("yMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  
  minor_matrix <- get(paste("minorCopyMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  major_matrix <- get(paste("majorCopyMatrix", mice[i], sep = ""), envir = .GlobalEnv)
  
  em <- get(paste("minorEpsilon", mice[i], sep = ""), envir = .GlobalEnv)
  eM <- get(paste("majorEpsilon", mice[i], sep = ""), envir = .GlobalEnv)
  
  project_name <- get(paste("projectname", mice[i], sep = ""), envir = .GlobalEnv)
  
  
  
  
  assign(paste("canopyInput", mice[i], sep = ""), 
         list(r_m, x_m, c_m, major_matrix, minor_matrix, eM, em, y_m, project_name))
  
  
  canopy_input <- get(paste("canopyInput", mice[i], sep = ""), envir = .GlobalEnv )
  
  
  
  names(canopy_input) <- c(paste("R", mice[i], sep = ""), paste("X", mice[i], sep = ""), paste("C", mice[i], sep = ""),
                           paste("WM", mice[i], sep = ""), paste("Wm", mice[i], sep = ""), paste("epsilonM", mice[i], sep = ""), 
                           paste("epsilonm", mice[i], sep = ""), paste("Y", mice[i], sep = ""), paste("projectname", mice[i], sep = ""))
  
  
  assign(paste("canopyInput", mice[i], sep = ""), canopy_input)
  
  
}


rm(c_m, em, eM, r_m , x_m, y_m, major_matrix, minor_matrix, canopy_input, project_name, currentMice)





# -----------------------------------------------------------------------------------------------------------------------------------------------





# Running Canopy Heterogeneity Estimation Method


library(Canopy)



#####  CNA and SNV input  ##### (below change according the input to be run, all mice or a specific mouse, canopyInputM49/M55/M61/M62)


projectname <- canopyInput$projectname # name of the project

R <- canopyInput$R # mutant allele read depth (for SNVs)

X <- canopyInput$X # total depth (for SNVs)

WM <- canopyInput$WM # observed major copy number (for CNA regions)

Wm <- canopyInput$Wm # observed minor copy number (for CNA regions)

epsilonM <- canopyInput$epsilonM # standard error of WM

epsilonm <- canopyInput$epsilonm # standard error of Wm

C <- canopyInput$C # whether CNA regions harbor specific CNAs (only needed for overlapping CNAs)

Y <- canopyInput$Y # whether SNVs are affected by CNAs






#####  MCMC sampling  ##### 


K <- 3:6 # number of subclones (since pyClone determined 4 - 5, we will put values around that estimate)

numchain <- 50 # number of chains with random initiations


# To sample the posterior trees. Major function of Canopy.
# Returns: List of sampleed trees in subtree space with different number of subclones; plot of posterior likelihoods in each subtree space generated (pdf format).

# writeskip parameter
# To speed up computation in MCMC and to avoid autocorrelation between neighboring trees that are samples consecutively, Canopy saves samples trees with a skip

sampchain <- canopy.sample(R = R, X = X, WM = WM, Wm = Wm, epsilonM = epsilonM,
                          epsilonm = epsilonm, C = C, Y = Y, K = K,
                          numchain = numchain, max.simrun = 50000,
                          min.simrun = 10000, writeskip = 200,
                          projectname = projectname, cell.line = TRUE,
                          plot.likelihood = TRUE)

save.image(file = paste(projectname, '_postmcmc_image.rda',sep=''),
           compress = 'xz')






