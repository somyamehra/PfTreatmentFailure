# Bayesian Approach

######### subroutines

# recode alleles
# genotypedata is genetic data, where first column (name "Sample ID") has the id of the sample, and rest of columns have the format nameoflocus_X, where X is the Xth allele detected
# the first column has to have a unique ID followed by either " Day 0" or " Day Failure"
# alleles_definitions is list of length number of loci
# each entry is a number of alleles by 2 matrix, where the first column is the lower bound, and the last column the upper bound

# output
# [[1]] is list of length number of loci
# each entry in the length is a vector of length number of ids
# each entry in the vector is a string with the following format:
# A-B-C/D-E
# the letters represent the Ath, Bth, Cth etc allele (as defined in alleles_definitions)
# the letters before the "/" represent day 0 alleles, and after day of failure alleles
# Example: 5-4/2 represents an individual that had alleles 5 and 4 on day 0, and then allele 2 on day of failure
# when nothing appears, either no alleles were detected, or alleles fell outside ranges in alleles_definitions 
# row names are the ids
# [[2]] is a number of ids X 2 matrix of multiplicity of infection, where first column is day 0 MOI and second column is day of failure MOI

recodeallele = function(alleles_definitions_subset,proposed) {
  
  ret = which(proposed > alleles_definitions_subset[,1] & proposed <= alleles_definitions_subset[,2])
  if (length(ret) == 0) {
    ret = NA
  }
  ret
}

recode_alleles = function(genotypedata, alleles_definitions) {
  
  ########### generate MOI for each sample
  
  ids = unique(unlist(strsplit(genotypedata$Sample.ID[grepl("Day 0",genotypedata$Sample.ID)]," Day 0")))
  locinames = unique(sapply(colnames(genotypedata)[-1],function(x) strsplit(x,"_")[[1]][1]))
  nids = length(ids)
  nloci = length(locinames)
  
  
  MOI0 = rep(0,nids)
  MOIf = rep(0,nids)
  
  # for each individual,
  # cycle through each locus and count number of separate alleles
  
  for (i in 1:nids) {
    for (j in 1:nloci) {
      locicolumns = grepl(paste(locinames[j],"_",sep=""),colnames(genotypedata))
      nalleles0 = sum(!is.na(genotypedata[grepl(paste(ids[i],"Day 0"),genotypedata$Sample.ID),locicolumns]))
      nallelesf = sum(!is.na(genotypedata[grepl(paste(ids[i],"Day Failure"),genotypedata$Sample.ID),locicolumns]))
      
      MOI0[i] = max(MOI0[i],nalleles0)
      MOIf[i] = max(MOIf[i],nallelesf)
    }
  }
  
  
  
  observeddatamatrix = list()
  for (j in 1:nloci) {
    locus = locinames[j]
    locicolumns = grepl(paste(locus,"_",sep=""),colnames(genotypedata))
    oldalleles = as.vector(genotypedata[,locicolumns])
    if (length(dim(oldalleles)[2]) == 0) {
      oldalleles = matrix(oldalleles,length(oldalleles),1)
    }
    newalleles = oldalleles
    ncolumns = dim(oldalleles)[2]
    for (i in 1:ncolumns) {
      newalleles[,i] = (sapply(1:dim(oldalleles)[1],function (x) recodeallele(alleles_definitions[[j]],oldalleles[x,i])))
    }
    newalleles[is.na(newalleles)] = ""
    
    tempobservedata = c()
    for (i in 1:nids) {
      locicolumns = grepl(paste(locinames[j],"_",sep=""),colnames(genotypedata))
      day0alleles = newalleles[grepl(paste(ids[i],"Day 0"),genotypedata$Sample.ID),]
      day0alleles = day0alleles[day0alleles != ""]
      dayfalleles = newalleles[grepl(paste(ids[i],"Day Failure"),genotypedata$Sample.ID),]
      dayfalleles = dayfalleles[dayfalleles != ""]
      tempobservedata[i] = paste(paste(sort(unique(as.numeric(day0alleles))),collapse="-"),paste(sort(unique(as.numeric(dayfalleles))),collapse="-"),sep="/")
    }
    observeddatamatrix[[j]] = tempobservedata
  }
  MOItemp = cbind(MOI0,MOIf)
  rownames(MOItemp) = ids
  list(observeddatamatrix = observeddatamatrix, MOI = MOItemp)
}

# generate definitions of alleles (ie binning)
# input:
# genotypedata is genetic data, where first column (name "Sample ID") has the id of the sample, and rest of columns have the format nameoflocus_X, where X is the Xth allele detected
# locirepeats is vector of length number of loci with type of locus (dinucleotide, trinucleotide etc repeats)
# maxk is a vector of length length number of loci with the maximum number of alleles for each locus 

# output
# list of length number of loci
# each entry is a number of alleles by 2 matrix, where the first column is the lower bound, and the last column the upper bound


define_alleles = function(genotypedata, locirepeats, maxk) {
  
  ids = genotypedata$Sample.ID
  
  locinames = unique(sapply(colnames(genotypedata)[-1],function(x) strsplit(x,"_")[[1]][1]))
  
  nids = length(ids)
  nloci = length(locinames)
  
  alleles = list()
  observed_data = list()
  for (j in 1:nloci) {
    locicolumns = grepl(paste(locinames[j],"_",sep=""),colnames(genotypedata))
    raw_alleles = c(as.matrix(genotypedata[,locicolumns]))
    raw_alleles = raw_alleles[!is.na(raw_alleles)]
    
    if (diff(range(raw_alleles)) < locirepeats[j]) {
      alleles[[j]] = matrix(c(min(raw_alleles)-locirepeats[j]/2,max(raw_alleles)+locirepeats[j]/2,length(raw_alleles)),1,3)
    } else {
      # remove outliers
      breaks = seq(from = floor(min(raw_alleles))-0.5, to = (max(raw_alleles)+1), by = 1)
      allele_values = round((breaks[2:length(breaks)] + breaks[1:(length(breaks)-1)]) / 2)
      hist_alleles = hist(raw_alleles, breaks = breaks, plot = FALSE)
      
      counts_by_offset = sapply(1:locirepeats[j], function (x) sum(hist_alleles$counts[seq(from = x, to = length(hist_alleles$counts), by = locirepeats[j])]))
      possible_alleles = allele_values[seq(from = which.max(counts_by_offset), to = length(allele_values), by = locirepeats[j])]
      
      if (min(raw_alleles) <= (min(possible_alleles)-locirepeats[j]/2)) {
        possible_alleles = c(min(possible_alleles-locirepeats[j]),possible_alleles)
      }
      if (max(raw_alleles) > (max(possible_alleles)+locirepeats[j]/2)) {
        possible_alleles = c(possible_alleles,max(possible_alleles+locirepeats[j]))
      }
      
      # assign clusters
      clusters = sapply(raw_alleles, function (x) which.min(abs(possible_alleles - x)))
      k = length(unique(clusters))
      
      colv = rep("white",length(possible_alleles))
      colv[1:length(possible_alleles) %in% unique(clusters)] = rainbow(k)

      # find break values (lower and upper)
      lower_break_value = sort(possible_alleles[unique(clusters)] - locirepeats[j]/2)
      upper_break_value = sort(possible_alleles[unique(clusters)] + locirepeats[j]/2)
      counts = sapply(1:length(lower_break_value), function (x) sum(raw_alleles > lower_break_value[x] & raw_alleles <= upper_break_value[x]))
      alleles[[j]] = cbind(lower_break_value, upper_break_value, counts)
    }
  }
  
  #### compress
  # take maxk most frequent alleles
  
  alleles2 = list()
  for (j in 1:nloci) {
    sortedindex = sort.int(alleles[[j]][,3],decreasing = TRUE,index.return = TRUE)$ix[1:maxk[j]]
    if (length(alleles[[j]][,3]) <= maxk[j]) {
      sortedindex = sort.int(alleles[[j]][,3],decreasing = TRUE,index.return = TRUE)$ix
    }
    alleles2[[j]] = cbind(alleles[[j]][sortedindex,1],alleles[[j]][sortedindex,2])
  }
  alleles2
}

# calculate frequencies of alleles
# input:
# genotypedata is genetic data, where first column (name "Sample ID") has the id of the sample, and rest of columns have the format nameoflocus_X, where X is the Xth allele detected
# alleles_definitions is list of length number of loci
# each entry is a number of alleles by 2 matrix, where the first column is the lower bound, and the last column the upper bound
# output:
# list of length number of loci
# each entry contains a vector with frequencies of each allele (might not sum to 1 if allele definitions do not cover all observed fragment lengths)
# output[[3]] is mean SD of within allele length

calculate_frequencies3 = function(genotypedata, alleles_definitions) {
  
  ids = genotypedata$Sample.ID
  locinames = unique(sapply(colnames(genotypedata)[-1],function(x) strsplit(x,"_")[[1]][1]))
  nids = length(ids)
  nloci = length(locinames)
  
  frequencies = list()
  
  variability = c()
  
  for (j in 1:nloci) {
    locicolumns = grepl(paste(locinames[j],"_",sep=""),colnames(genotypedata))
    raw_alleles = c(as.matrix(genotypedata[,locicolumns]))
    raw_alleles = raw_alleles[!is.na(raw_alleles)]
    low = alleles_definitions[[j]][,1]
    high = alleles_definitions[[j]][,2]
    frequencies[[j]] = sapply(1:dim(alleles_definitions[[j]])[1],function (x) sum(raw_alleles > low[x] & raw_alleles <= high[x]))
    meanSD = mean(sapply(1:dim(alleles_definitions[[j]])[1],function (x) sd(raw_alleles[raw_alleles > low[x] & raw_alleles <= high[x]])),na.rm=TRUE)
    if(is.na(meanSD)) {meanSD = 0}
    variability[j] = meanSD
    frequencies[[j]] = frequencies[[j]] / length(raw_alleles)
  }
  freqmatrix = matrix(0,nloci,max(unlist(lapply(frequencies,length))))
  
  for (j in 1:nloci) {
    freqmatrix[j,1:length(frequencies[[j]])] = frequencies[[j]]
  }
  
  ret = list()
  ret[[1]] = unlist(lapply(frequencies,length))
  ret[[2]] = freqmatrix
  ret[[3]] = variability
  ret
}

rowMeans2 = function(x){
  if (length(dim(x)) == 0) {
    ret = mean(x)
  } else {
    ret = rowMeans(x)
  }
  ret
}

# locirepeats: this is the size of the repeat region for each microsatellite (ie. for dimers use 2, for trimers use 3)
# nruns: number of runs for MCMC algorithm
run_Plucinski_classifier <- function(genotype_data, nruns, locirepeats=c(2,2,3,3,3,3,3), maxalleles=30) {
  # calculate burnin (number of runs to discard) and record interval (which n_th iterations should be recorded)
  record_interval = ceiling(nruns / 1000);
  burnin = ceiling(nruns * 0.25);
  
  ### IMPORT
  
  options(java.parameters = "-Xmx4096m")
  
  genotypedata_latefailures <- genotype_data %>% subset(Type=="Paired") %>% select(-Type)
  additional_genotypedata <- genotype_data %>% subset(Type=="Additional") %>% select(-Type)
  
  ### RUN ALGORITHM
  
  ### identify arms (based on site column)
  
  site_names = unique(genotypedata_latefailures$Site)
  
  state_classification_all = c()
  state_parameters_all = c()
  ids_all = c()
  
  summary_statisics <- posterior_recrudesence <- list()
  
  for (site in site_names) {
    jobname = site
    genotypedata_RR = genotypedata_latefailures[genotypedata_latefailures$Site==site,-c(2)]
    additional_neutral = additional_genotypedata[additional_genotypedata$Site==site,-c(2)]
    if (dim(additional_neutral)[1] > 0) { additional_neutral$Sample.ID = paste("Additional_",1:dim(additional_neutral)[1],sep="")}
    
    maxMOI = max(as.numeric(sapply(1:length(colnames(genotypedata_RR)), function (x) strsplit(colnames(genotypedata_RR)[x],"_")[[1]][2])),na.rm=TRUE)
    
    ids = unique(unlist(strsplit(genotypedata_RR$Sample.ID[grepl("Day 0",genotypedata_RR$Sample.ID)]," Day 0")))
    
    locinames = unique(sapply(colnames(genotypedata_RR)[-1],function(x) strsplit(x,"_")[[1]][1]))
    nloci = length(locinames)
    
    nids = length(ids)
    
    k = rep(maxalleles, nloci)
    alleles_definitions_RR  = define_alleles(rbind(genotypedata_RR,additional_neutral),locirepeats,k)
    
    ##### calculate MOI
    MOI0 = rep(0,nids)
    MOIf = rep(0,nids)
    for (i in 1:nids) {
      for (j in 1:nloci) {
        locicolumns = grepl(paste(locinames[j],"_",sep=""),colnames(genotypedata_RR))
        nalleles0 = sum(!is.na(genotypedata_RR[grepl(paste(ids[i],"Day 0"),genotypedata_RR$Sample.ID),locicolumns]))
        nallelesf = sum(!is.na(genotypedata_RR[grepl(paste(ids[i],"Day Failure"),genotypedata_RR$Sample.ID),locicolumns]))
        
        MOI0[i] = max(MOI0[i],nalleles0)
        MOIf[i] = max(MOIf[i],nallelesf)
      }
    }
    #maxMOI = max(MOI0,MOIf)
    
    
    ##### define statevector
    
    alleles0 = matrix(0,nids,maxMOI*nloci)
    recoded0 = matrix(0,nids,maxMOI*nloci)
    hidden0 = matrix(NA,nids,maxMOI*nloci)
    recr0 = matrix(NA,nids,nloci)
    recr_repeats0 = matrix(NA,nids,nloci) # number of times recrudescing allele is repeated on day 0
    recr_repeatsf = matrix(NA,nids,nloci) # number of times recrudescing allele is repeated on day 0
    allelesf = matrix(0,nids,maxMOI*nloci)
    recodedf = matrix(0,nids,maxMOI*nloci)
    hiddenf = matrix(NA,nids,maxMOI*nloci)
    recrf = matrix(NA,nids,nloci)
    if (length(additional_neutral) > 0) { if (dim(additional_neutral)[1] > 0) {
      recoded_additional_neutral = matrix(0,dim(additional_neutral)[1],maxMOI*nloci)
    }}
    mindistance = matrix(0,nids,nloci)
    alldistance = array(NA,c(nids,nloci,maxMOI*maxMOI))
    allrecrf = array(NA,c(nids,nloci,maxMOI*maxMOI))
    classification = rep(0,nids)
    ##### create state 0
    
    for (j in 1:nloci) {
      locus = locinames[j]
      locicolumns = grepl(paste(locus,"_",sep=""),colnames(genotypedata_RR))
      oldalleles = as.vector(genotypedata_RR[,locicolumns]) %>% do.call(cbind, .) # add cbind
      if (length(dim(oldalleles)[2]) == 0) {
        oldalleles = matrix(oldalleles,length(oldalleles),1)
      }
      newalleles = oldalleles
      ncolumns = dim(oldalleles)[2]
      for (i in 1:ncolumns) {
        newalleles[,i] = (sapply(1:dim(oldalleles)[1],function (x) recodeallele(alleles_definitions_RR[[j]],oldalleles[x,i])))
      }
      newalleles = matrix(as.numeric(unlist(c(newalleles))),dim(newalleles)[1],dim(newalleles)[2])
      newalleles[is.na(newalleles)] = 0
      oldalleles = matrix(as.numeric(unlist(c(oldalleles))),dim(oldalleles)[1],dim(oldalleles)[2])
      oldalleles[is.na(oldalleles)] = 0
      
      oldalleles[newalleles == 0] = 0
      alleles0[,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + dim(oldalleles)[2])] = oldalleles[grepl("Day 0",genotypedata_RR$Sample.ID),]
      allelesf[,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + dim(oldalleles)[2])] = oldalleles[grepl("Day Failure",genotypedata_RR$Sample.ID),]
      recoded0[,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + dim(newalleles)[2])] = newalleles[grepl("Day 0",genotypedata_RR$Sample.ID),]
      recodedf[,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + dim(newalleles)[2])] = newalleles[grepl("Day Failure",genotypedata_RR$Sample.ID),]
      
    }
    
    if (length(additional_neutral) > 0) { if (dim(additional_neutral)[1] > 0) {
      recoded_additional_neutral = matrix(0,dim(additional_neutral)[1],maxMOI*nloci)
      ##### recode additional_neutral
      for (j in 1:nloci) {
        locus = locinames[j]
        locicolumns = grepl(paste(locus,"_",sep=""),colnames(genotypedata_RR))
        oldalleles = as.vector(additional_neutral[,locicolumns]) %>% do.call(cbind, .) # add cbind
        if (length(dim(oldalleles)[2]) == 0) {
          oldalleles = matrix(oldalleles,length(oldalleles),1)
        }
        newalleles = oldalleles
        ncolumns = dim(oldalleles)[2]
        for (i in 1:ncolumns) {
          newalleles[,i] = (sapply(1:dim(oldalleles)[1],function (x) recodeallele(alleles_definitions_RR[[j]],oldalleles[x,i])))
        }
        newalleles = matrix(as.numeric(unlist(c(newalleles))),dim(newalleles)[1],dim(newalleles)[2])
        newalleles[is.na(newalleles)] = 0
        oldalleles = matrix(as.numeric(unlist(c(oldalleles))),dim(oldalleles)[1],dim(oldalleles)[2])
        oldalleles[is.na(oldalleles)] = 0
        
        oldalleles[newalleles == 0] = 0
        recoded_additional_neutral[,(maxMOI*(j-1)+1) : (maxMOI*(j-1) + dim(oldalleles)[2])] = newalleles
      }
    } else {
      recoded_additional_neutral = c()
    }}
    
    ## estimate frequencies
    
    frequencies_RR = calculate_frequencies3(rbind(genotypedata_RR,additional_neutral),alleles_definitions_RR)
    
    ## assign random hidden alleles
    for (i in 1:nids) {
      for (j in 1:nloci) {
        nalleles0 = sum(alleles0[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] != 0)
        nmissing0 = MOI0[i] - nalleles0
        whichnotmissing0 = ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(alleles0[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOI0[i])] != 0)]
        whichmissing0 = ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(alleles0[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOI0[i])] == 0)]
        
        if (nalleles0 > 0) {
          hidden0[i,whichnotmissing0] = 0
        }
        if (nmissing0 > 0) {
          newhiddenalleles0 = sample(1:(frequencies_RR[[1]][j]),nmissing0,replace=TRUE,frequencies_RR[[2]][j,1:(frequencies_RR[[1]][j])])
          recoded0[i,whichmissing0] = newhiddenalleles0
          alleles0[i,whichmissing0] = rowMeans(alleles_definitions_RR[[j]])[newhiddenalleles0] # hidden alleles get mean allele length
          hidden0[i,whichmissing0] = 1
        }
        nallelesf = sum(allelesf[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] != 0)
        nmissingf = MOIf[i] - nallelesf
        whichnotmissingf = ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(allelesf[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOIf[i])] != 0)]
        whichmissingf = ((maxMOI*(j-1)+1) : (maxMOI*(j)))[which(allelesf[i,(maxMOI*(j-1)+1) : (maxMOI*(j-1)+MOIf[i])] == 0)]
        
        if (nallelesf > 0) {
          hiddenf[i,whichnotmissingf] = 0
        }
        if (nmissingf > 0) {
          newhiddenallelesf = sample(1:(frequencies_RR[[1]][j]),nmissingf,replace=TRUE,frequencies_RR[[2]][j,1:(frequencies_RR[[1]][j])])
          recodedf[i,whichmissingf] = newhiddenallelesf
          allelesf[i,whichmissingf] = rowMeans(alleles_definitions_RR[[j]])[newhiddenallelesf] # hidden alleles get mean allele length
          hiddenf[i,whichmissingf] = 1
        }
      }
    }
    
    ## initial estimate of q (probability of an allele being missed)
    qq = mean(c(hidden0,hiddenf),na.rm=TRUE)
    
    ## initial estimate of dvect (likelihood of error in analysis)
    dvect = rep(0,1+round(max(sapply(1:nloci,function (x) diff(range(c(alleles_definitions_RR[[x]])))))))
    #dvect[1] = 0.75
    #dvect[2] = 0.2
    #dvect[3] = 0.05
    dvect = dgeom(0:(length(dvect)-1),0.75)
    ## randomly assign recrudescences/reinfections
    for (i in 1:nids) {
      z = runif(1)
      if (z < 0.5) {
        classification[i] = 1
      }
      for (j in 1:nloci) { # determine which alleles are recrudescing (for beginning, choose closest pair)
        allpossiblerecrud = expand.grid(1:MOI0[i],1:MOIf[i])
        closestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (x) abs(alleles0[i,maxMOI*(j-1)+allpossiblerecrud[x,1]] - allelesf[i,maxMOI*(j-1)+allpossiblerecrud[x,2]])))
        mindistance[i,j] = abs(alleles0[i,maxMOI*(j-1)+allpossiblerecrud[closestrecrud,1]] - allelesf[i,maxMOI*(j-1)+allpossiblerecrud[closestrecrud,2]])
        alldistance[i,j,1:dim(allpossiblerecrud)[1]] = sapply(1:dim(allpossiblerecrud)[1], function (x) abs(alleles0[i,maxMOI*(j-1)+allpossiblerecrud[x,1]] - allelesf[i,maxMOI*(j-1)+allpossiblerecrud[x,2]]))
        allrecrf[i,j,1:dim(allpossiblerecrud)[1]] = recodedf[i,maxMOI*(j-1)+allpossiblerecrud[,2]]
        recr0[i,j] = maxMOI*(j-1)+allpossiblerecrud[closestrecrud,1]
        recrf[i,j] = maxMOI*(j-1)+allpossiblerecrud[closestrecrud,2]
        recr_repeats0[i,j] = sum(recoded0[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] == recoded0[i,recr0[i,j]])
        recr_repeatsf[i,j] = sum(recodedf[i,(maxMOI*(j-1)+1) : (maxMOI*(j))] == recodedf[i,recrf[i,j]])
      }
    }
    
    
    #### correction factor (reinfection)
    correction_distance_matrix = list() # for each locus, matrix of distances between each allele
    for (i in 1:nloci) {
      correction_distance_matrix[[i]] = as.matrix(dist(rowMeans(alleles_definitions_RR[[i]])))
    }
    
    # Switch hidden - update hidden observations as part of Gibbs MCMC 
    switch_hidden = function(x) {
      z = runif(1)
      if (sum(hidden0[x,], hiddenf[x,],na.rm=TRUE) > 0) { # if hidden alleles exist
        if (length(which(c(hidden0[x,], hiddenf[x,])==1))>1) {
          chosen = sample(which(c(hidden0[x,], hiddenf[x,])==1),1)
        } else {
          chosen = which(c(hidden0[x,], hiddenf[x,])==1)
        }
        if (classification[x] == 0) { # reinfection
          if (chosen <= nloci*maxMOI) { # day 0 hidden allele
            chosenlocus = ceiling(chosen/maxMOI)
            old = recoded0[x,chosen]
            new = sample(1:frequencies_RR[[1]][chosenlocus],1)
            
            
            oldalleles = recoded0[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hidden0[x,] == 0))]
            repeatedold = qq
            repeatednew = qq
            if (sum(oldalleles == old) >= 1) { # if old allele is a repeat, don't penalize with missing probability
              repeatedold = 1;
            }
            if (sum(oldalleles == new) >= 1) { # if new allele is a repeat, don't penalize with missing probability
              repeatednew = 1;
            }
            alpha = (sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,new]+1]) * repeatednew) / 
              (sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,old]+1]) * repeatedold)
            if (z < alpha) { # switch made
              recoded0[x,chosen] <<- new
              ####### new allele should have some variability (ie what's being see in real data)
              newallele_length = mean(alleles_definitions_RR[[chosenlocus]][new,]) + rnorm(1,mean=0,sd=frequencies_RR[[3]][chosenlocus])
              alleles0[x,chosen] <<- newallele_length
              allpossiblerecrud = expand.grid(1:MOI0[x],1:MOIf[x])
              closestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (y) abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]])))
              mindistance[x,chosenlocus] <<- abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,2]])
              alldistance[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- sapply(1:dim(allpossiblerecrud)[1], function (y) abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]]))
              allrecrf[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- recodedf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[,2]]
              recr0[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,1]
              recrf[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,2]
              recr_repeats0[x,chosenlocus] <<- sum(recoded0[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recoded0[x,recr0[x,chosenlocus]])
              recr_repeatsf[x,chosenlocus] <<- sum(recodedf[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recodedf[x,recrf[x,chosenlocus]])
              # }
            }
          } else { # day f hidden allele
            chosen = chosen - nloci*maxMOI
            chosenlocus = ceiling(chosen/maxMOI)
            old = recodedf[x,chosen]
            new = sample(1:frequencies_RR[[1]][chosenlocus],1)
            oldalleles = recodedf[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hiddenf[x,] == 0))]
            repeatedold = qq
            repeatednew = qq
            if (sum(oldalleles == old) >= 1) { # if old allele is a repeat, don't penalize with missing probability
              repeatedold = 1;
            }
            if (sum(oldalleles == new) >= 1) { # if new allele is a repeat, don't penalize with missing probability
              repeatednew = 1;
            }
            alpha = (sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,new]+1]) * repeatednew) / 
              (sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,old]+1]) * repeatedold)
            if (z < alpha) { # switch made
              recodedf[x,chosen] <<- new
              newallele_length = mean(alleles_definitions_RR[[chosenlocus]][new,]) + rnorm(1,mean=0,sd=frequencies_RR[[3]][chosenlocus])
              
              ####### new allele should have some variability (ie what's being see in real data)
              allelesf[x,chosen] <<- newallele_length
              allpossiblerecrud = expand.grid(1:MOI0[x],1:MOIf[x])
              closestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (y) abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]])))
              mindistance[x,chosenlocus] <<- abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,2]])
              alldistance[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- sapply(1:dim(allpossiblerecrud)[1], function (y) abs(alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]]))
              allrecrf[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- recodedf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[,2]]
              recr0[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,1]
              recrf[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[closestrecrud,2]
              recr_repeats0[x,chosenlocus] <<- sum(recoded0[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recoded0[x,recr0[x,chosenlocus]])
              recr_repeatsf[x,chosenlocus] <<- sum(recodedf[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recodedf[x,recrf[x,chosenlocus]])
              # }
            }
          }
        } else { # recrudescence
          if (chosen <= nloci*maxMOI) { # day 0 hidden allele
            chosenlocus = ceiling(chosen/maxMOI)
            old = recoded0[x,chosen]
            new = sample(1:frequencies_RR[[1]][chosenlocus],1)
            newallele_length = mean(alleles_definitions_RR[[chosenlocus]][new,]) + rnorm(1,mean=0,sd=frequencies_RR[[3]][chosenlocus])
            
            oldalleles = recoded0[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hidden0[x,] == 0))]
            repeatedold = qq
            repeatednew = qq
            if (sum(oldalleles == old) >= 1) { # if old allele is a repeat, don't penalize with missing probability
              repeatedold = 1;
            }
            if (sum(oldalleles == new) >= 1) { # if new allele is a repeat, don't penalize with missing probability
              repeatednew = 1;
            }
            allpossiblerecrud = expand.grid(1:MOI0[x],1:MOIf[x])
            tempalleles = alleles0[x,maxMOI*(chosenlocus-1)+1:maxMOI]
            tempalleles[chosen-(chosenlocus-1)*maxMOI] = newallele_length 
            temprecoded = recoded0[x,maxMOI*(chosenlocus-1)+1:maxMOI]
            temprecoded[chosen-(chosenlocus-1)*maxMOI] = new
            
            newclosestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (y) abs(tempalleles[allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]])))
            newmindistance = abs(tempalleles[allpossiblerecrud[newclosestrecrud,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,2]])
            newalldistance = sapply(1:dim(allpossiblerecrud)[1], function (y) abs(tempalleles[allpossiblerecrud[y,1]] - allelesf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,2]]))
            newallrecrf = recodedf[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[,2]]
            
            # calculate new multiple-comparisons coefficient
            newrecr0 = maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,1]
            newrecrf = maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,2]
            newrecr_repeats0 = sum(temprecoded == temprecoded[allpossiblerecrud[newclosestrecrud,1]],na.rm=TRUE)
            newrecr_repeatsf = sum(recodedf[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recodedf[x,newrecrf])
            
            likelihoodnew = mean(dvect[round(newalldistance)+1]/sapply(1:length(newallrecrf), function (z) sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,newallrecrf[z]]+1])),na.rm=TRUE) * repeatednew
            likelihoodold = mean(dvect[round(alldistance[x,chosenlocus,])+1]/sapply(1:(maxMOI*maxMOI), function (z) sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,allrecrf[x,chosenlocus,z]]+1])),na.rm=TRUE) * repeatedold 
            
            if (likelihoodnew  == likelihoodold) {
              # if both num and denominator are equal (for case when both are 0..., otherwise 0/0 gives NaN)
              alpha = 1							 
            } else {
              alpha = likelihoodnew / likelihoodold 						 
            }
            
            if (z < alpha) { # switch made
              recoded0[x,chosen] <<- new
              alleles0[x,chosen] <<- newallele_length
              mindistance[x,chosenlocus] <<- newmindistance
              alldistance[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- newalldistance
              allrecrf[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- newallrecrf
              recr0[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,1]
              recrf[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,2]
              recr_repeats0[x,chosenlocus] <<- sum(recoded0[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recoded0[x,recr0[x,chosenlocus]])
              recr_repeatsf[x,chosenlocus] <<- sum(recodedf[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recodedf[x,recrf[x,chosenlocus]])
            }
          } else { # day f hidden allele
            chosen = chosen - nloci*maxMOI
            chosenlocus = ceiling(chosen/maxMOI)
            old = recodedf[x,chosen]
            new = sample(1:frequencies_RR[[1]][chosenlocus],1)
            newallele_length = mean(alleles_definitions_RR[[chosenlocus]][new,]) + rnorm(1,mean=0,sd=frequencies_RR[[3]][chosenlocus])
            
            oldalleles = recodedf[x,intersect(((chosenlocus-1)*maxMOI+1):((chosenlocus)*maxMOI),which(hiddenf[x,] == 0))]
            repeatedold = qq
            repeatednew = qq
            if (sum(oldalleles == old) >= 1) { # if old allele is a repeat, don't penalize with missing probability
              repeatedold = 1;
            }
            if (sum(oldalleles == new) >= 1) { # if new allele is a repeat, don't penalize with missing probability
              repeatednew = 1;
            }
            allpossiblerecrud = expand.grid(1:MOI0[x],1:MOIf[x])
            tempalleles = allelesf[x,maxMOI*(chosenlocus-1)+1:maxMOI]
            tempalleles[chosen-(chosenlocus-1)*maxMOI] = newallele_length 
            temprecoded = recodedf[x,maxMOI*(chosenlocus-1)+1:maxMOI]
            temprecoded[chosen-(chosenlocus-1)*maxMOI] = new
            newclosestrecrud = which.min(sapply(1:dim(allpossiblerecrud)[1], function (y) abs(tempalleles[allpossiblerecrud[y,2]] - alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]])))
            newmindistance = abs(tempalleles[allpossiblerecrud[newclosestrecrud,2]] - alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,1]])
            newalldistance = sapply(1:dim(allpossiblerecrud)[1], function (y) abs(tempalleles[allpossiblerecrud[y,2]] - alleles0[x,maxMOI*(chosenlocus-1)+allpossiblerecrud[y,1]]))
            newallrecrf = temprecoded[allpossiblerecrud[,2]]
            
            # calculate new multiple-comparisons coefficient
            newrecr0 = maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,1]
            newrecrf = maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,2]
            newrecr_repeats0 = sum(recoded0[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recoded0[x,newrecr0])
            newrecr_repeatsf = sum(temprecoded == temprecoded[allpossiblerecrud[newclosestrecrud,2]],na.rm=TRUE)
            
            likelihoodnew = mean(dvect[round(newalldistance)+1]/sapply(1:length(newallrecrf), function (z) sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,newallrecrf[z]]+1])),na.rm=TRUE) * repeatednew
            likelihoodold = mean(dvect[round(alldistance[x,chosenlocus,])+1]/sapply(1:(maxMOI*maxMOI), function (z) sum(frequencies_RR[[2]][chosenlocus,1:frequencies_RR[[1]][chosenlocus]]*dvect[correction_distance_matrix[[chosenlocus]][,allrecrf[x,chosenlocus,z]]+1])),na.rm=TRUE) * repeatedold 
            
            if (likelihoodnew  == likelihoodold) {
              # if both num and denominator are equal (for case when both are 0..., otherwise 0/0 gives NaN)
              alpha = 1							 
            } else {
              alpha = likelihoodnew / likelihoodold 						 
            }
            if (z < alpha) { # switch made
              recodedf[x,chosen] <<- new
              allelesf[x,chosen] <<- newallele_length
              mindistance[x,chosenlocus] <<- newmindistance
              alldistance[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- newalldistance
              allrecrf[x,chosenlocus,1:dim(allpossiblerecrud)[1]] <<- newallrecrf
              recr0[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,1]
              recrf[x,chosenlocus] <<- maxMOI*(chosenlocus-1)+allpossiblerecrud[newclosestrecrud,2]
              recr_repeats0[x,chosenlocus] <<- sum(recoded0[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recoded0[x,recr0[x,chosenlocus]])
              recr_repeatsf[x,chosenlocus] <<- sum(recodedf[x,(maxMOI*(chosenlocus-1)+1) : (maxMOI*(chosenlocus))] == recodedf[x,recrf[x,chosenlocus]])
            }
          }
        }
      }
    }
    
    findposteriorfrequencies = function(x,tempdata) {
      data = tempdata[,c(1:maxMOI)+(x-1)*maxMOI];
      nalleles = frequencies_RR[[1]][x]
      freq_prior_alpha = rep(1,nalleles);
      freq_posterior_alpha = freq_prior_alpha + table(factor(c(data),levels=c(1:nalleles)));
      frequencies_RR[[2]][x,1:nalleles] <<- rdirichlet(1, freq_posterior_alpha);
    }
    
    
    state_classification = matrix(NA,nids,(nruns-burnin)/record_interval)
    state_alleles0 = array(NA,c(nids,maxMOI*nloci,(nruns-burnin)/record_interval))
    state_allelesf = array(NA,c(nids,maxMOI*nloci,(nruns-burnin)/record_interval))
    state_parameters = matrix(NA,2+2*nloci,(nruns-burnin)/record_interval)
    
    count = 1
    dposterior = 0.75
    runmcmc = function() {
      # propose new classification
      likelihoodratio = sapply(1:nids, function (x) exp(sum(log(sapply(1:nloci, function (y) mean(dvect[round(alldistance[x,y,])+1]/sapply(1:(maxMOI*maxMOI), function (z) sum(frequencies_RR[[2]][y,1:frequencies_RR[[1]][y]]*dvect[correction_distance_matrix[[y]][,allrecrf[x,y,z]]+1])),na.rm=TRUE))))))
      
      z = runif(nids)
      newclassification = classification
      newclassification[classification == 0 & z < likelihoodratio] = 1
      newclassification[classification == 1 & z < 1/likelihoodratio] = 0
      classification <<- newclassification
      
      # propose new hidden states
      sapply(1:nids, function (x) switch_hidden(x))
      
      # propose q (beta distribution is conjugate distribution for binomial process)
      q_prior_alpha = 0;
      q_prior_beta = 0;
      q_posterior_alpha = q_prior_alpha + sum(c(hidden0,hiddenf) == 1,na.rm=TRUE)
      q_posterior_beta = q_prior_beta + sum(c(hidden0,hiddenf)==0,na.rm=TRUE)
      if (q_posterior_alpha == 0) {
        q_posterior_alpha =1
      }
      qq <<- rbeta(1, q_posterior_alpha , q_posterior_beta)
      
      #  update dvect (approximate using geometric distribution)
      # only if there is at least 1 recrudescent infection
      if (sum(classification==1) >= 1) {
        d_prior_alpha = 0;
        d_prior_beta = 0;
        d_posterior_alpha = d_prior_alpha + length(c(mindistance[classification==1,]))
        d_posterior_beta = d_prior_beta + sum(c(round(mindistance[classification==1,])))
        if (d_posterior_beta == 0) {
          d_posterior_beta = sum(c((mindistance[classification==1,])))
        }
        if (d_posterior_beta == 0) { ## algorithm will get stuck if dposterior is allowed to go to 1
          d_posterior_beta = 1
        }	
        
        
        dposterior <<- rbeta(1, d_posterior_alpha , d_posterior_beta)
        dvect = (1-dposterior) ^ (1:length(dvect)-1) * dposterior
        dvect <<- dvect / (sum(dvect))
      }
      
      # update frequencies
      # remove recrudescing alleles from calculations
      tempdata = recoded0
      sapply(which(classification == 1), function (x) tempdata[x,recr0[x,]] <<- 0)
      tempdata = rbind(tempdata, recodedf)
      sapply(1:nloci, function (x) findposteriorfrequencies(x,rbind(tempdata,recoded_additional_neutral)))
      
      # record state
      if (count > burnin & count %% record_interval == 0) {
        state_classification[,(count-burnin)/record_interval] <<- classification
        state_alleles0[,,(count-burnin)/record_interval] <<- alleles0
        state_allelesf[,,(count-burnin)/record_interval] <<- allelesf
        state_parameters[1,(count-burnin)/record_interval] <<- qq
        state_parameters[2,(count-burnin)/record_interval] <<- dposterior
        state_parameters[3:(3+nloci-1),(count-burnin)/record_interval] <<- apply(frequencies_RR[[2]],1,max)
        state_parameters[(3+nloci):(3+2*nloci-1),(count-burnin)/record_interval] <<- sapply(1:nloci,function (x) sum(frequencies_RR[[2]][x,]^2))
        
      }
      count <<- count + 1
    }
    
    replicate(nruns,runmcmc())
    
    ## make sure no NAs in result matrices
    state_parameters = state_parameters[,!is.na(colSums(state_parameters))]
    state_classification = state_classification[,!is.na(colSums(state_classification))]
    
    
    posterior_recrudesence[[jobname]] <- data.frame(isolate=ids, prob=rowMeans2(state_classification))
    
    summary_statisticsmatrix = cbind(format(rowMeans(state_parameters),digits=2),
                                     apply(format(t(sapply(1:dim(state_parameters)[1], function (x) quantile(state_parameters[x,],c(0.25,0.75)))),digits=2),1, function (x) paste(x,collapse="–")))
    summary_statisticsmatrix = rbind(summary_statisticsmatrix, c(format(mean(state_parameters[(3+nloci):(3+2*nloci-1),]),digits = 2),paste(format(quantile(state_parameters[(3+nloci):(3+2*nloci-1),],c(0.25,0.75)),digits=2),collapse="–")))
    summary_statisticsmatrix = as.matrix(sapply(1:dim(summary_statisticsmatrix)[1], function (x) paste(summary_statisticsmatrix[x,1], " (",summary_statisticsmatrix[x,2],")",sep="")))
    rownames(summary_statisticsmatrix) = c("q","d",locinames,locinames,"Mean diversity")
    summary_statisics[[jobname]] <- summary_statisticsmatrix[c("q", "d", "Mean diversity"),]
  }
  
  run_summary <- list(posterior_recrudesence=bind_rows(posterior_recrudesence, .id="Site"),
                      summary_statisics=bind_rows(summary_statisics, .id="Site"))
  
  return(run_summary)
}

