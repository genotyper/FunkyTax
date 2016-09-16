CatFun <-
function(tf, dds, func_col = "Gene", alpha = 0.05) {
  
  #categorize functions as significant or not for both tests
  sigs = tf[[1]][,5] < alpha
  sigs = tf[[1]][sigs,]
  nonsigs = tf[[1]][,5] >= alpha
  nonsigs = tf[[1]][nonsigs,]
  dds = dds[complete.cases(dds),]
  dds = dds[dds$padj < alpha,]
  
  #Compare the results of DESeq2 and TaFuR to make the evolutionary classifications
  Dsig = row.names(dds)
  sig_in_both = subset(sigs, sigs[[func_col]] %in% Dsig)
  sig_in_ds = subset(Dsig, !Dsig %in% sigs[[func_col]])
  ns_ds = data.frame(t(tf[[5]]))
  ns_ds = subset(ns_ds, !row.names(ns_ds) %in% Dsig)
  ns_ds = unique(row.names(ns_ds))
  ns_in_both = subset(nonsigs, nonsigs[[func_col]] %in% ns_ds)
  sig_in_kegg = subset(sigs, sigs[[func_col]] %in% ns_ds)
  
  #just organizing a data.frame for the plot and summary
  who_master = data.frame(names(sort(colMeans(tf[[5]]), decreasing = TRUE)))
  row.names(who_master) = who_master[,1]
  names(who_master)[1] = func_col
  sig_in_both$Divergent = rep("Divergent", length(sig_in_both[,1]))
  test = merge(who_master, sig_in_both, by = func_col, all.x = T)
  sig_in_ds = data.frame(sig_in_ds)
  names(sig_in_ds)[1] = func_col
  sig_in_ds$Enhanced = rep("Enhanced", length(sig_in_ds[,1]))
  test = merge(test, sig_in_ds, by = func_col, all.x = T)
  ns_in_both$Conserved = rep("Conserved", length(ns_in_both[,1]))
  test = merge(test, ns_in_both, by = func_col, all.x = T)
  sig_in_kegg$Equivalent = rep("Equivalent", length(sig_in_kegg[,1]))
  test = merge(test, sig_in_kegg, by = func_col, all.x = T)
  groups = test[,c(1,6, 7, 12, 17)]
  groups$Not.Tested = rep(NA, length(groups[,1]))
  groups$Not.Tested = ((is.na(groups$Divergent) == is.na(groups$Conserved)) == (is.na(groups$Enhanced) == is.na(groups$Equivalent)))
  groups$Not.Tested = gsub("TRUE", "Not Tested", groups$Not.Tested)
  groups$Not.Tested = gsub("FALSE", NA, groups$Not.Tested)
  grouper = reshape2::melt(groups, id.vars = func_col)
  grouper = grouper[complete.cases(grouper),]
  grouper = grouper[,-3]
  who = names(sort(colMeans(tf[[5]]), decreasing = TRUE))
  grouper = grouper[match(who, grouper[[func_col]]),]
  grouper$variable = gsub("Not.Tested", "No Test", grouper$variable)
  grouper$variable = ordered(grouper$variable, levels = c("Divergent", "Enhanced", "Conserved", "Equivalent", "No Test"))
  samplelevels = c(as.character(grouper[[func_col]]))
  grouper[func_col] = factor(grouper[[func_col]], ordered=T, levels = samplelevels)
  grouper["func_col"] = grouper[func_col]
  library(ggplot2)
  p = ggplot(grouper, aes(func_col, fill = variable)) + geom_bar() + 
    theme_bw(base_size = 15) +  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) + 
    theme(legend.text = element_text(face = "italic")) + xlab("") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    theme(legend.title = element_blank()) + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) + 
    ylab("") + facet_wrap(~variable, ncol = 1)+ 
    theme(strip.background = element_blank(), strip.text.x = element_blank()) + 
    theme(panel.margin = unit(.02, "lines"))
  final = vector("list", 2)
  final[[1]] = p
  final[[2]] = grouper
  return(final)
}
