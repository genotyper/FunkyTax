TaFuR <-
function(dat, met, no.processors = 1, method = "bray", fact = "Genus", func_col = "Gene", samp_col = "Sample", otu_col = "OTU.ID", count_col = "CountContributedByOTU"){
  #the result of this function is a list of length 4 containing 1) per function test statistics, 2) per function ADONIS tables, 3) per function community matrices, and 4) functions private to each factor level
  #dat = data.frame summarizing each OTU's functional contribution for each sample. (e.g. metagenome contributions from PICRUSt (Langille et al. 2013) or other method) 
  #met = data.frame of metadata including factor of interest and sample IDs that correspond to those used in dat
  #no.processors = the number of processors you would like to use for distance matrix construction and ADONIS
  #method = any dissimilarity index method applicable to function vegdist (vegan)
  #fact = factor of interest occuring in met
  #func_col = ID of the column in dat containing the genes/function names
  #samp_col = ID of the column in dat containing the sample names (which correspond to a column in met with the same name)
  #otu_col = ID of the column in dat containing the OTU names
  #count_col = ID of the column in dat containing the gene counts per OTU per sample
  
  
#   library(reshape2)
#   library(vegan)
  print("building a list of community matrices per function")
  temp = data.frame(dat[func_col], dat[samp_col], dat[otu_col],  dat[count_col])
#   list = as.character(met[[samp_col]])
#   temp = subset(temp, temp$Sample %in% list)
  temp = merge(temp, met, by = "Sample")
  mylist = split(temp, temp[func_col], drop = T)
  mylist = lapply(mylist, function(x) {x[func_col] <- NULL; x})
  for (i in 1:length(mylist)) {colnames(mylist[[i]])[2] = names(mylist)[[i]]}
  nl = names(mylist)
    myfunc <- function(varname) {                                                
    formula <- as.formula(paste(samp_col, " + ", fact, " ~",varname, sep=""))
    reshape2::dcast(mylist[[varname]], formula, value.var= count_col, fun.aggregate = sum)
  }
  mylist = lapply(names(mylist), myfunc)
  myfunc2 <- function(x) {row.names(x) = x[,1]
  return(x)}
  mylist = lapply(mylist, myfunc2)
  myfunc3 <- function(x) {x = x[,-1] 
  return(x)}
  mylist = lapply(mylist, myfunc3)
  names(mylist) = nl

  #if a function was only observed once its list position will not be a data.frame, removing these
  amount = (lapply(mylist, dim))
  amount[sapply(amount, is.null)] <- NULL
  mylist = subset(mylist, names(mylist) %in% names(amount))
  
  #keeping those contributed by more than 1 OTU (only multivariate tables)
  list = sapply(mylist, function(i) length(i) == 2)
  list = names(which(list == FALSE))
  mylist = subset(mylist, names(mylist) %in% list)
  
  #keeping those that occur in at least 2 of the factor levels (so adonis can proceed)
  list = sapply(mylist, function(i) length(unique((i[[fact]]))) < 2)
  notprivate = names(which(list == FALSE))
  
  #keeping a list of those that only occur in 1 factor level (private) so as to add to the final results
  private = names(which(list == TRUE))
  private = subset(mylist, names(mylist) %in% private)
  private = sapply(private, function(x) as.character(private[[1]][1,1]))
  private = data.frame(func_col = names(private), Private = unname(private))
  
  mylist = subset(mylist, names(mylist) %in% notprivate)
  library(foreach)
  #building distance matrices and doing ADONIS
  if(no.processors > 1){
    print(paste("building the distance matrices using ", no.processors, " processors", sep = ""))
    cl = snow::makeCluster(no.processors)
    doSNOW::registerDoSNOW(cl)
    mydis = foreach(i = 1:length(mylist)) %dopar% {
      vegan::vegdist(sqrt(mylist[[i]][,-1]/rowSums(mylist[[i]][-1])), method = method)
    }
    print(paste("performing ADONIS for each distance matrix in parallel mode using ", no.processors, " processors", sep = ""))
    ado = vector("list", length(mydis))
    for (i in 1:length(mylist)) {
      print(paste(i, " of ", length(mylist), sep = ""))
      ado[[i]] = suppressMessages(vegan::adonis(as.dist(mydis[[i]]) ~ mylist[[i]][,1]))
    }
    snow::stopCluster(cl)
  }
  
  
  if(no.processors == 1){
    print("building the distance matrices using a single processor")
    mydis = vector("list", length(mylist))
    for (i in 1:length(mylist)){
      mydis[[i]] = vegan::vegdist(sqrt(mylist[[i]][,-1]/rowSums(mylist[[i]][-1])), method = method)
    }
    print("performing ADONIS for each distance matrix using a single processor")
    ado = vector("list", length(mydis))
    for (i in 1:length(mylist)) {
      ado[[i]] = suppressMessages(vegan::adonis(as.dist(mydis[[i]]) ~ mylist[[i]][,1]))
    }
  }

  #extract the F and p-values
  pvals = c()
  for (i in 1:length(ado)) {pvals[i] = ado[[i]]$aov.tab[1,6]}
  fvals = c()
  for (i in 1:length(ado)) {fvals[i] = ado[[i]]$aov.tab[1,4]}
  rsqs = c()
  for (i in 1:length(ado)) {rsqs[i] = ado[[i]]$aov.tab[1,5]}
  nl = names(mylist)
  results = data.frame(temp = nl, Fstatistic = fvals, P.values = pvals, Rsq = rsqs, P.adjusted = stats::p.adjust(pvals, method ="fdr"))
  names(results)[1] = func_col
  #building a OTU x Function matrix
  temp = reshape2::dcast(dat, OTU.ID ~ Gene, sum, value.var = count_col)
  row.names(temp) = temp[,1]
  temp = temp[-1]
  
  #wrapping things up
  final = vector("list", 5)
  final[[1]] = results
  # final[[1]] = rbind(final[[1]], private)
  final[[2]] = ado
  final[[3]] = mydis
  final[[4]] = private
  final[[5]] = temp
  return(final)
}
