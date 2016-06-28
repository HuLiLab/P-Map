setwd('/home/edroaldo/Documents/Projects/PMAP/GitHub')
library(gplots)
library(igraph)
source("iGraph_GexF_Exporter.R")

######### Loading required files ##########
expMat <- get(load('expClean_drug_annotated_removedGenes_May_18_2015.R'))
stTable <-  get(load('sampTable_May_07_2015.R'))
drug2sample <- get(load('drug2sample_mapping_May_08_2015.R'))
lcl.drugs <- get(load('lcl.drugs_EC50_May_08_2015.R'))

ec50.expr.tables <- list()
ec50.samp.tables <- list()
for(drug in colnames(lcl.drugs)){
  for(group in as.vector(unique(stTable$group))){
    samples.tmp <- stTable[which(stTable$group == group),]
    tmpExp <- expMat[,rownames(samples.tmp)]
    tmp <- data.frame()
    for(sample in rownames(samples.tmp)){
      tmp['EC50', sample] <- lcl.drugs[sample, drug]
    }
    tmpExp <- rbind(tmp, tmpExp)
    tmpExp2 <- tmpExp[,order(tmpExp['EC50',])]
    tmpExp2<-tmpExp2[,colSums(is.na(tmpExp2)) == 0]
    columns <- c(1:15, ((ncol(tmpExp2) - 14):ncol(tmpExp2)))
    tmpExp3 <- tmpExp2[2:nrow(tmpExp2), columns]
    
    sensitive <- rep('sensitive', 15)
    resistant <- rep('resistant', 15)
    conditions <- c(sensitive, resistant)
    x <- data.frame(cell_lines=colnames(tmpExp3), conditions=conditions, stringsAsFactors=FALSE)
    
    ec50.expr.tables[[drug]][[group]] <- tmpExp3
    ec50.samp.tables[[drug]][[group]] <- x
    #filename <- paste('Choong_tables/', drug, '_', group, '.csv',sep='')
    #write.csv(tmpExp2, file=filename)
  }
}

iNetwork <- get(load('PPI_network.R'))

######### Antracyclines ##############
geneSigs.ant <- getGeneSigs(lcl.drugs, stTab, ec50.samp.tables, ec50.expr.tables, c('doxorubicin', 'epirubicin'), 'results/')
## strategy 1
df.sens.ant <- createDF(geneSigs.ant, 'sensitive', 2)
df.res.ant <- createDF(geneSigs.ant, 'resistant', 2)
write.csv(df.sens.ant, file='results/anthracycline_sensitive.csv')
write.csv(df.res.ant, file='results/anthracycline_resistant.csv')

##strategy 4
s4.sens.ant <- getNeighs.strategy4(iNetwork, df.sens.ant)
s4.res.ant <- getNeighs.strategy4(iNetwork, df.res.ant)
write.csv(s4.sens.ant$df, file='results/sensitive_for_frequency_calculation_anthracyclines.csv')
write.csv(s4.res.ant$df, file='results/resistant_for_frequency_calculation_anthracyclines.csv')

freq.gene1.sens.ant <- computeFrequency(s4.sens.ant$n1, 'results/frequency_gene1_neighbors_sensitive_anthracyclines.csv')
freq.gene2.sens.ant <- computeFrequency(s4.sens.ant$n2, 'results/frequency_gene2_neighbors_sensitive_anthracyclines.csv')
freq.gene1.res.ant <- computeFrequency(s4.res.ant$n1, 'results/frequency_gene1_neighbors_resistant_anthracyclines.csv')
freq.gene2.res.ant <- computeFrequency(s4.res.ant$n2, 'results/frequency_gene2_neighbors_resistant_anthracyclines.csv')

gene1.tab <- combFreqTables(freq.gene1.sens.ant, freq.gene1.res.ant, 'results/combined_frequency_table_gene1_sens_res_anthracyclines.csv')
gene2.tab <- combFreqTables(freq.gene2.sens.ant, freq.gene2.res.ant, 'results/combined_frequency_table_gene2_sens_res_anthracyclines.csv')

#order by 
gene1.tab <- gene1.tab[order(gene1.tab$difference, decreasing=FALSE),]
gene1.tab <- gene1.tab[-1,] #remove UBC
selected.genes <- gene1.tab[1:31,]

network <- createDRN(iNetwork, rownames(selected.genes), df.res.ant, "Anthracycline_Resistant_N1_for_Network")


######### Taxol ##############
geneSigs.taxol <- getGeneSigs(lcl.drugs, stTab, ec50.samp.tables, ec50.expr.tables, c('docetaxel', 'pactilatxel'), 'results/')
## strategy 1
df.sens.taxol <- createDF(geneSigs.taxol, 'sensitive', 2)
df.res.taxol <- createDF(geneSigs.taxol, 'resistant', 2)
write.csv(df.sens.taxol, file='results/taxol_sensitive.csv')
write.csv(df.res.taxol, file='results/taxol_resistant.csv')

##strategy 4
s4.sens.taxol <- getNeighs.strategy4(iNetwork, df.sens.taxol)
s4.res.taxol <- getNeighs.strategy4(iNetwork, df.res.taxol)
write.csv(s4.sens.taxol$df, file='results/sensitive_for_frequency_calculation_taxol.csv')
write.csv(s4.res.taxol$df, file='results/resistant_for_frequency_calculation_taxol.csv')

freq.gene1.sens.taxol <- computeFrequency(s4.sens.taxol$n1, 'results/frequency_gene1_neighbors_sensitive_taxol.csv')
freq.gene2.sens.taxol <- computeFrequency(s4.sens.taxol$n2, 'results/frequency_gene2_neighbors_sensitive_taxol.csv')
freq.gene1.res.taxol <- computeFrequency(s4.res.taxol$n1, 'results/frequency_gene1_neighbors_resistant_taxol.csv')
freq.gene2.res.taxol <- computeFrequency(s4.res.taxol$n2, 'results/frequency_gene2_neighbors_resistant_taxol.csv')

gene1.tab <- combFreqTables(freq.gene1.sens.taxol, freq.gene1.res.taxol, 'results/combined_frequency_table_gene1_sens_res_taxol.csv')
gene2.tab <- combFreqTables(freq.gene2.sens.taxol, freq.gene2.res.taxol, 'results/combined_frequency_table_gene2_sens_res_taxol.csv')

#order by 
gene1.tab <- gene1.tab[order(gene1.tab$difference, decreasing=TRUE),]
gene1.tab <- gene1.tab[-1,] #remove UBC
selected.genes <- gene1.tab[1:32,]

network <- createDRN(iNetwork, rownames(selected.genes), df.sens.taxol, "Taxane_Sensitive_N1_for_Network")

#functions
createDRN <- function(iNetwork, nodes2map, df, filename){
  g<-simplify(iNetwork, remove.multiple=TRUE, remove.loops=TRUE)
  nodes2map <- intersect(nodes2map, V(g)$name)
  first.neighs <- induced.subgraph(graph=g,vids=unlist(neighborhood(graph=g,order=1,nodes=nodes2map)))
  all.nodes <- V(first.neighs)$name
  first.nodes <- setdiff(all.nodes, nodes2map)
  nodes2remove <- setdiff(first.nodes, rownames(df)) #nodes to remove from first.neighs subgraph...
  clean.subgraph <- delete.vertices(first.neighs, nodes2remove)
  V(clean.subgraph)$label <- V(clean.subgraph)$name
  V(clean.subgraph)$colors <- ''
  V(clean.subgraph)[nodes2map]$colors <-  "#30EC1C" #'darkgreen' N1
  V(clean.subgraph)[intersect(first.nodes, rownames(df))]$colors <-  '#FA4141'#"#DA34DA" #PDEG
  
  V(clean.subgraph)$nodeType <- ''
  V(clean.subgraph)[nodes2map]$nodeType <-  "N1" #'darkgreen' N1
  V(clean.subgraph)[intersect(first.nodes, rownames(df))]$nodeType <- "PDEG" #PDEG
  
  saveAsGEXF(clean.subgraph, paste(filename, '.gexf', sep=''));
  
  df <- as.data.frame(list(genes=V(clean.subgraph)$name, type=V(clean.subgraph)$nodeType, 
                           color=V(clean.subgraph)$colors))
  write.csv(df, paste(filename, '.csv', sep=''))
  
  clean.subgraph
}

combFreqTables <- function(freq.sens, freq.res, filename){
  all.genes <- unique(c(names(freq.sens), names(freq.res)))
  df <- data.frame(matrix(0, nrow=length(all.genes), ncol=3))
  rownames(df) <- all.genes
  colnames(df) <- c('frequency_sensitive', 'frequency_resistant', 'difference')
  for(gene in rownames(df)){
    if(gene %in% names(freq.sens) & gene %in% names(freq.res)){
      df[gene, 'frequency_sensitive'] <- freq.sens[gene]
      df[gene, 'frequency_resistant'] <- freq.res[gene]
      df[gene, 'difference'] <- freq.sens[gene] - freq.res[gene]
    }else if(gene %in% names(freq.sens)){
      df[gene, 'frequency_sensitive'] <- freq.sens[gene]
      df[gene, 'difference'] <- freq.sens[gene]
    }else{
      df[gene, 'frequency_resistant'] <- freq.res[gene]
      df[gene, 'difference'] <- freq.res[gene]
    }
  }
  df <- df[order(df$difference, decreasing=TRUE),]
  
  write.csv(df, file=filename)
  df
}


computeFrequency <- function(x, filename){
  #x <- s4.sens.ant$n2
  x.all <- unique(as.vector(unlist(x)))
  table <- as.data.frame(matrix(0, nrow=length(x.all), ncol=length(names(x))))
  rownames(table) <- x.all
  colnames(table) <- names(x)
  for(g in x.all){
    for(edge in names(x)){
      if(g %in% x[[edge]]){
        table[g, edge] <- 1
      }
    }
  }
  frequency <- rowSums(table)
  frequency <- frequency[order(frequency, decreasing=TRUE)]
  
  write.csv(frequency, file=filename)
  
  frequency
}

getGeneSigs <- function(lcl.drugs, stTab, ec50.samp.tables, ec50.expr.tables, drugs, prefix){
  geneSigs <- list()
  for(drug in colnames(lcl.drugs)){
    #if(grepl('doxorubicin', drug) | grepl('epirubicin', drug)){
    if(grepl(drugs[1], drug) | grepl(drugs[2], drug)){
      for(group in as.vector(unique(stTable$group))){
        cat(drug, '\t', group, '\n')
        x <- ec50.samp.tables[[drug]][[group]]
        xx <- ec50.expr.tables[[drug]][[group]]
        specGenes <- cn_findSpecGenes(xx, x, 0.95, TRUE, 'conditions')
        #specGenes <- cn_findSpecGenes(xx, x, 0.99, TRUE, 'conditions')
        
        aux <- x
        group.colors <- rainbow(length(unique(aux$conditions)))
        aux[which(aux$conditions == "sensitive"),'colors'] <- group.colors[2]
        aux[which(aux$conditions == "resistant"),'colors'] <- group.colors[1]
        colorVector <- aux$colors
        
        tmp <- xx[as.vector(specGenes$sensitive),]
        hr <- hclust(as.dist(1-abs(cor(t(tmp), method="spearman"))), method="ward.D")
        filename <- paste(prefix, 'sig_genes_SENSITIVE_', drug, '_', group, '_', Sys.Date(), '.pdf', sep='')
        #pdf(file=filename, width=6, height=60);
        pdf(file=filename, width=6, height=7);
        bk = unique(c(seq(-2, -1, length=70), seq(-1, 1, length=70), seq(1, 2, length=70)))
        hmcols<- colorRampPalette(c("blue","white","red"))#(length(bk)-1)
        heatmap.2(as.matrix(tmp), breaks=bk, col=hmcols,trace="none",
                  density.info="none", scale="row",margin=c(15,15), key=TRUE, Colv=FALSE, Rowv=FALSE, ylab="",
                  srtCol=60, dendrogram="none", cexCol=0.75, cexRow=0.65, ColSideColors=colorVector)
        dev.off();
        
        tmp <- xx[as.vector(specGenes$resistant),]
        hr <- hclust(as.dist(1-abs(cor(t(tmp), method="spearman"))), method="ward.D")
        filename <- paste(prefix, 'sig_genes_RESISTANT_', drug, '_', group, '_', Sys.Date(), '.pdf', sep='')
        #pdf(file=filename, width=6, height=60);
        pdf(file=filename, width=6, height=7);
        bk = unique(c(seq(-2, -1, length=70), seq(-1, 1, length=70), seq(1, 2, length=70)))
        hmcols<- colorRampPalette(c("blue","white","red"))#(length(bk)-1)
        heatmap.2(as.matrix(tmp), breaks=bk, col=hmcols,trace="none",
                  density.info="none", scale="row",margin=c(15,15), key=TRUE, Colv=FALSE, Rowv=FALSE, ylab="",
                  srtCol=60, dendrogram="none", cexCol=0.75, cexRow=0.65, ColSideColors=colorVector)
        dev.off();
        
        geneSigs[[drug]][[group]] <- specGenes
        
      }
    }
  }
  geneSigs
}

createDF <- function(geneSigs, phenotype, threshold){
  allGenes <- unique(as.vector(unlist(geneSigs)))
  df <- data.frame()
  for(gene in allGenes){
    for(drug in names(geneSigs)){
      for(group in names(geneSigs[[drug]])){
        signature <- geneSigs[[drug]][[group]][[phenotype]]
        cname <- paste(drug, group, sep='_')
        if(gene %in% signature){
          df[gene, cname] <- 1
        }else{
          df[gene, cname] <- 0
        }
      }
    }
  }
  df2 <- df[rowSums(df) >= threshold,]
  df2
}

#PDEGs
getNeighs.strategy4 <- function(iNetwork, df){
  allGenes <- rownames(df)
  cat('number of genes: ', length(allGenes), '\n')
  
  g2<-igraph::simplify(iNetwork, remove.multiple=TRUE, remove.loops=TRUE) 
  nodes <- allGenes
  nodes <- intersect(nodes, V(g2)$name)  #all PDEGs in the PPI
  
  ##strategy 4
  df3 <- data.frame()
  n1 <- list()
  n2 <- list()
  for(gene1 in nodes){ #for each PDEG
    #get neighbors of genes that survived.
    g1_neighs <- induced.subgraph(graph=g2,vids=unlist(neighborhood(graph=g2,order=1,nodes=gene1)))
    neigh.names <- V(g1_neighs)$name #first interacting partners of PDEG gene1
    xOverlap <- intersect(neigh.names, nodes)
    if(length(xOverlap) >= 2){
      for(gene2 in neigh.names){
        g2_neighs <- induced.subgraph(graph=g2,vids=unlist(neighborhood(graph=g2,order=1,nodes=gene2)))
        neigh2.names <- V(g2_neighs)$name
        
        xxOverlap <- intersect(neigh2.names, nodes) #nodes ->PDEGs
        if(length(xxOverlap) >= 2){
          if(!('UBC' %in% neigh2.names)){
            edge <- paste(gene1, gene2, sep='_')
            df3[edge, 'num_neighbors_gene1'] <- length(neigh.names)
            df3[edge, 'num_neighbors_gene2'] <- length(neigh2.names)
            df3[edge, 'neighbors_gene1'] <- paste(neigh.names, collapse=',')
            df3[edge, 'neighbors_gene2'] <- paste(neigh2.names, collapse=',')
            n1[[edge]] <- neigh.names
            n2[[edge]] <- neigh2.names
          }
        }
      }
    }
  }
  list(df=df3, n1=n1, n2=n2)
}


#### The code below was taken from the CellNet source code
## CellNet: Network Biology Applied to Stem Cell  Engineering. Cahan P*, Li H*, Morris SA*, Lummertz da Rocha E, Daley GQ, Collins JJ.
#  Cell. 2014 Aug 14;158(4):903-15. PMID: 25126793

cn_findSpecGenes<-function# find genes that are preferentially expressed in specified samples
(expDat, ### expression matrix
 sampTab, ### sample table
 qtile=0.95, ### quantile
 remove=FALSE,
 dLevel="population_name" #### annotation level to group on
){
  cat("Template matching...\n")
  myPatternG<-cn_sampR_to_pattern(as.vector(sampTab[,dLevel]));
  specificSets<-apply(myPatternG, 1, cn_testPattern, expDat=expDat);
  
  # adaptively extract the best genes per lineage
  cat("First pass identification of specific gene sets...\n")
  cvalT<-vector();
  ctGenes<-list();
  ctNames<-unique(as.vector(sampTab[,dLevel]));
  for(ctName in ctNames){
    x<-specificSets[[ctName]];
    cval<-quantile(x$cval, qtile, na.rm = TRUE);
    tmp<-rownames(x[x$cval>cval,]);
    ctGenes[[ctName]]<-tmp;
    cvalT<-append(cvalT, cval);
  }
  #specificSets[['306G']][1:5,]
  if(remove){
    cat("Prune common genes...\n");
    # now limit to genes exclusive to each list
    specGenes<-list();
    for(ctName in ctNames){
      others<-setdiff(ctNames, ctName);
      x<-setdiff( ctGenes[[ctName]], unlist(ctGenes[others]));
      specGenes[[ctName]]<-x;
    }
    result <- specGenes
  }else {
    result <- ctGenes;
  }
  result
}

cn_sampR_to_pattern<-function# return a pattern for use in cn_testPattern (template matching)
(sampR){
  d_ids<-unique(as.vector(sampR));
  nnnc<-length(sampR);
  ans<-matrix(nrow=length(d_ids), ncol=nnnc);
  for(i in seq(length(d_ids))){
    x<-rep(0,nnnc);
    x[which(sampR==d_ids[i])]<-1;
    ans[i,]<-x;
  }
  colnames(ans)<-as.vector(sampR);
  rownames(ans)<-d_ids;
  ans;
}

cn_testPattern<-function(pattern, expDat){
  pval<-vector();
  cval<-vector();
  geneids<-rownames(expDat);
  llfit<-ls.print(lsfit(pattern, t(expDat)), digits=25, print=FALSE);
  xxx<-matrix( unlist(llfit$coef), ncol=8,byrow=TRUE);
  ccorr<-xxx[,6];
  cval<- sqrt(as.numeric(llfit$summary[,2])) * sign(ccorr);
  pval<-as.numeric(xxx[,8]);
  
  #qval<-qvalue(pval)$qval;
  holm<-p.adjust(pval, method='holm');
  #data.frame(row.names=geneids, pval=pval, cval=cval, qval=qval, holm=holm);
  data.frame(row.names=geneids, pval=pval, cval=cval,holm=holm);
}
