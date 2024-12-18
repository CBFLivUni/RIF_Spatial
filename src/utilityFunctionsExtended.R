#utility function to summary the QC flags for a NanoStringGeoMxSet
QCSummary <- function(object) {
  # Collate QC Results
  QCResults <- protocolData(object)[["QCFlags"]]
  flag_columns <- colnames(QCResults)
  QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                           Warning = colSums(QCResults[, flag_columns]))
  
  QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
  })
  
  QC_Summary["TOTAL FLAGS", ] <-
    c(sum(QCResults[, "QCStatus"] == "PASS"),
      sum(QCResults[, "QCStatus"] == "WARNING"))
  
  
  warn_formatter <- formatter("span", 
                              style = x ~ style( "background-color" = ifelse(x > 0 , "yellow", "white")))
  
  pass_formatter <- formatter("span", 
                              style = x ~ style( "background-color" = ifelse(x > 0 , "lightgreen", "white")))
  
  QCTable <- formattable(QC_Summary, list(Pass=pass_formatter,Warning=warn_formatter),caption = "Summary of QC flags")
  
  return(list(QCResults=QCResults,QCTable=QCTable))
}


# Graphical summaries of QC statistics plot function- taken from the GeoMxWorkflows vignette
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = "segment",
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

#function to filter a geomxdataset to keep samples suitable for running a random effects model
#i.e contains multiple observations per slide
filterDataWithinSlide <- function(dataset){
  
  pheno <- as.data.frame(pData(dataset))%>% rownames_to_column("ID") 
  
  #keep samples where the staining is present on multiple slides per Disease.
  phenoFiltered <- pheno %>% group_by(Disease,slide) %>% filter(n_distinct(segment)>1) %>%
    group_by(Disease,segment) %>% filter(n() > 1)  %>% group_by(Disease) %>%
    mutate(segmentTypes=n_distinct(segment)) %>% group_by(Disease,slide) %>%
    filter(!segmentTypes > n_distinct(segment))
  
  ind <- which(pheno$ID %in% phenoFiltered$ID)
  
  return(dataset[,ind])
  
}

#function to allow group split but keeping the names of the groups in the resulting list
#modified from https://github.com/tidyverse/dplyr/issues/4223#issuecomment-469269857
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!group_keys(grouped), sep = "_"))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}

#function to plot the lmerTest output as a volcano plot
plotVolcano <- function(expData,title,foldChangeColumn="Estimate",FDRColumn="FDR",n=10){
  
  #make a new column indicating if the gene is regarded as differentially expressed or not
  expData$test <- ifelse(abs(expData[,foldChangeColumn])>=log2(1.5) & expData[,FDRColumn]<=0.05,"yes","no")
  expData$logFDR <- -log10(expData[,FDRColumn])
  
  expData.selected <- expData[order(expData[,FDRColumn]),]

  selectedGenesUp <- expData.selected[ expData.selected[,foldChangeColumn]>0,]
  selectedGenesDown <- expData.selected[ expData.selected[,foldChangeColumn]<0,]
  selectedGenes <- c(selectedGenesUp[ 1:n,"Gene"],selectedGenesDown[ 1:n,"Gene"])
  
  #add highlighted genes if significant
  toHighlight <- read.delim("data/toHighlight.txt")
  toHighlight <- toHighlight[ toHighlight$Comparison==title,"Gene"]
  toHighlight <- toHighlight[ toHighlight %in% expData[ expData$test == "yes","Gene"]]
  
  selectedGenes <- c(selectedGenes,toHighlight)
  
  
  g<- ggplot(expData, aes_string(x=foldChangeColumn, y="logFDR")) +
    geom_point(aes(colour=test), size=1, alpha=0.3) +
    scale_colour_manual(values=c('grey', 'red')) +
    geom_vline(xintercept=log2(1.5), colour='blue',linetype=2) +
    geom_vline(xintercept=-log2(1.5), colour='blue',linetype=2) +
    geom_hline(yintercept=-log10(0.05), colour='blue',linetype=2) +
    cowplot::theme_cowplot(font_size = 20) +  theme(legend.position = 'none') +
    xlab("log2 fold change") + ylab("-log10 (adjusted p-value)")  +
    geom_text_repel(data = subset(expData,Gene %in% selectedGenes), size = 5,aes(label=Gene),max.overlaps = Inf,force = 10,min.segment.length = 0,seed = 42) +
    ggtitle(title)
  
  }

#function to plot a MAplot given the differential expression results and the negative probes
MAPlot <- function(diffExp,contrasts=1,FCColumn="Estimate",FDRColumn="FDR",negProbes="NegProbe-WTX",n=8) {
  
  diffExp$Color <- ifelse(diffExp[,FDRColumn]<=0.05,"yes","no")
  
  cutoff <- quantile(diffExp$meanExp, 0.8)
  diffExp$absFC <- abs(diffExp[,FCColumn])
  
  FCColumn <- sym(FCColumn)
  FDRColumn <- sym(FDRColumn)
  
  if(contrasts == 1){
    labelData <- diffExp %>% group_by(Subset) %>% arrange(!!FDRColumn) %>%
      filter(!!FDRColumn<0.01 & absFC >1 & meanExp > cutoff ) %>% slice_head(n=n)
  } else {
    labelData <- diffExp %>% group_by(Subset,Contrast) %>% arrange(!!FDRColumn) %>%
      filter(!!FDRColumn<0.01 & absFC >1 & meanExp > cutoff ) %>% slice_head(n=n)
  }
  
  
  
  g <- ggplot(subset(diffExp, !Gene %in% negProbes),
              aes(x = meanExp, y = !!FCColumn,
                  color = Color, label = Gene)) +
    geom_hline(yintercept = c(0.5, -0.5), lty = "dashed") +
    geom_point(alpha = 0.5) + 
    labs(y = "log2 Fold Change",
         x = "Mean Expression",
         color = "FDR ≤ 0.05") +
    scale_color_manual(values = c(yes = "dodgerblue",
                                  no = "grey")) +
    geom_text_repel(data = labelData,
                    size = 4, point.padding = 0.15, color = "black",
                    min.segment.length = .1, box.padding = .2, lwd = 2,max.overlaps = Inf,seed=42) +
    theme_bw(base_size = 16)
  
  
  if(contrasts == 1){
    g <- g +  facet_wrap(~Subset,ncol = 1)
  } else {
    g <- g + facet_wrap(~Subset + Contrast)
  }
  
  return(g)
}



data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

plotStripChart <- function(GOI,spatialData,elt="q_norm") {
  
  ggplot(pData(spatialData),
         aes(x = Disease, fill = Disease, color=`slide name`,
             y = assayDataElement(spatialData[GOI, ],
                                  elt = elt))) +
    geom_jitter(width = .2) +
    stat_summary(fun.data=data_summary, color="grey") +
    labs(y = paste(GOI,"Expression")) +
    facet_wrap(~segment) +
    theme_cowplot()
}


getReactomePathways <- function(species, outFile){
  
  #the reactome pathways are only provided with entrez IDs
  #so we need to convert them to gene symbols to match the data we have
  
  symbol2eg <- dplyr::case_when(
    species == "human" ~ "org.Hs.egSYMBOL",
    species == "mouse" ~ "org.Mm.egSYMBOL"
  )
  
  symbol2eg <- as.list(get(symbol2eg))
  
  #get eg to reactome pathway
  reactome2eg <- as.list(reactomePATHID2EXTID)
  
  speciesID <- dplyr::case_when(
    species == "human" ~ "R-HSA",
    species == "mouse" ~ "R-MMU"
  )
  
  #filter to obtain pathways for the selected species
  reactome2eg <- reactome2eg[grep(speciesID,names(reactome2eg))]
  
  #function to search through the pathway
  grepREACTOME <- function(id,mapkeys){
    unique(unlist(mapkeys[id],use.names=FALSE))
  }
  
  #convert the entrez ids to gene symbols for each pathway
  reactome <- lapply(reactome2eg,grepREACTOME,symbol2eg)
  
  #get the pathway names rather than the ids
  reactome2name <- as.list(reactomePATHID2NAME)
  reactomeNames <- sapply(names(reactome),grepREACTOME,reactome2name)
  term2Name <- data.frame(term=names(reactome),name=reactomeNames)
  
  reactome <- stack(reactome)[,2:1]
  colnames(reactome) <- c("term","gene")
  
  return(list(TERM2GENE=reactome,TERM2NAME=term2Name))
  
}

# function to get the results table for a limma contrast
getResultsDataFrame <- function(fit2, contrast, numerator, 
                                denominator) {
  data <- topTable(fit2, coef = contrast, number = Inf, 
                   sort.by = "P")
  #data <- independentFiltering(data, filter = data$AveExpr, objectType = "limma")
  #colnames(data) <- paste(paste(numerator, denominator, sep = "vs"), colnames(data), sep = "_")
  return(data)
}

# function to run lme4

lme4BetweenDisease <- function(object,assaySlot="log_q",segments=NULL,cores=7){
  
  if(is.null(segments)){
    segments <- unique(pData(object)$segment)
  }
  
  DEResults <- c()
  for(segment in segments) {
    ind <- pData(object)$segment == segment
    mixedOutmc <-
      mixedModelDE(object[, ind],
                   elt = assaySlot,
                   modelFormula = ~ Region + (1 | PatientID),
                   groupVar = "Region",
                   nCores = cores,
                   multiCore = FALSE)
    
    # format results as data.frame
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    
    # use lapply in case you have multiple levels of your test factor to
    # correctly associate gene name with it's row in the results table
    r_test$Gene <- 
      unlist(lapply(colnames(mixedOutmc),
                    rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Subset <- segment
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                         "Pr(>|t|)", "FDR")]
    DEResults <- rbind(DEResults, r_test)
  }
  return(DEResults)
}


lme4DE <- function(object,assaySlot="log_q",segments=NULL,cores=7){
  
  DEResults <- c()
  mixedOutmc <-
      mixedModelDE(object,
                   elt = assaySlot,
                   modelFormula = ~ Region + (1 | PatientID),
                   groupVar = "Region",
                   nCores = cores,
                   multiCore = FALSE,pAdjust = "BH")
    
    # format results as data.frame
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    
    # use lapply in case you have multiple levels of your test factor to
    # correctly associate gene name with it's row in the results table
    r_test$Gene <- 
      unlist(lapply(colnames(mixedOutmc),
                    rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Contrast", "Estimate", 
                         "Pr(>|t|)", "FDR")]
    DEResults <- rbind(DEResults, r_test)
  
  return(DEResults)
}


makeColours <- function(values,palette="Set1") {
  
  if(palette=="nrc") {
    v <- pal_npg(palette)(length(unique(values)))
  } else{
    v <- brewer.pal(length(unique(values)), palette)
  }
  
  v <- v[1:length(unique(values))]
  names(v) <- sort(unique(values))
  return(v)
}


simplifyGeneOntology <- function(go_analysis,threshold=0.7,orgdb="org.Hs.eg.db",ont="BP",method="Rel") {
  
  simMatrix <- calculateSimMatrix(go_analysis$ID,
                                  orgdb=orgdb,
                                  ont=ont,
                                  method=method)
  scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=threshold,
                                  orgdb=orgdb)
  return(list(reducedTerms,simMatrix))
}

improvedScatterPlot <- function (simMatrix, reducedTerms, size = "score", addLabel = TRUE, 
          labelSize = 3) 
{
  if (!all(sapply(c("ggplot2", "ggrepel","cowplot"), requireNamespace, 
                  quietly = TRUE))) {
    stop("Packages ggplot2, ggrepel, cowplot and/or its dependencies not available. ", 
         "Consider installing them before using this function.", 
         call. = FALSE)
  }
  x <- cmdscale(as.matrix(as.dist(1 - simMatrix)), eig = TRUE, 
                k = 2)
  df <- cbind(as.data.frame(x$points), reducedTerms[match(rownames(x$points), 
                                                          reducedTerms$go), c("term", "parent", "parentTerm", "size")])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = V1, y = V2, color = parentTerm)) + 
    ggplot2::geom_point(ggplot2::aes(size = size), alpha = 0.5) + 
    ggplot2::scale_color_discrete(guide = "none") + ggplot2::scale_size_continuous(guide = "none", 
                                                                                   range = c(0, 25)) + ggplot2::scale_x_continuous(name = "") + 
    ggplot2::scale_y_continuous(name = "") + cowplot::theme_cowplot() + xlab("MDS1") + ylab("MDS2")
  if (addLabel) {
    p + ggrepel::geom_label_repel(aes(label = parentTerm), 
                                  data = subset(df, parent == rownames(df)), box.padding = grid::unit(1, 
                                                                                                      "lines"), size = labelSize,max.overlaps = Inf)
  }
  else {
    p
  }
}

save_pheatmap_png <- function(x, filename, width=12, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height,units = "in",res=300)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#function to plot performance measures for an lme4 based model
checkLME4Model <- function(geneName,object,modelFormula,elt="log_q"){
  
  dat <- assayDataElement(object, elt = elt)
  
  mTerms <- all.vars(modelFormula)
  pDat <- sData(object)[, mTerms]
  for (i in names(pDat)) {
    if (inherits(i, "character")) {
      pDat[, i] <- as.factor(pDat[, i])
    }
  }
  modelFormula <- formula(paste("expr",as.character(modelFormula)[2], sep = " ~ "))
  
  dat <- data.frame(expr = dat[geneName, ], pDat)
  lmOut <- lmerTest::lmer(modelFormula,dat)
  check_model(lmOut,panel = TRUE)
}

getGeneSymbol <- function(ids){
  
  ids <- unlist(strsplit(ids,split = "/"))
  symbols <- bitr(ids,"ENTREZID",  "SYMBOL", OrgDb = org.Hs.eg.db)
  symbols <- paste(sort(symbols$SYMBOL),collapse=" ")
  return(symbols)
}



#function run nichenet analysis given expressed ligands, receptors and differentially expressed genes
runNicheNetAnalysis <- function(expressed_genes_sender, expressed_genes_receiver, geneset_oi, senderName, receiverName, weighted_networks){
  
  
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  #lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds?download=1"))
  
  # If wanted, users can remove ligand-receptor interactions that were predicted based on protein-protein interactions and only keep ligand-receptor interactions that are described in curated databases. To do this: uncomment following line of code:
  # lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  
  ligands = lr_network %>% pull(from) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  
  receptors = lr_network %>% pull(to) %>% unique()
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
  head(lr_network_expressed)
  
  potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
  head(potential_ligands)
  
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities <- ligand_activities %>% arrange(-aupr) 
  
  best_upstream_ligands = ligand_activities %>% top_n(10, aupr) %>% arrange(-aupr) %>% pull(test_ligand)
  
  p_hist_lig_activity = ggplot(ligand_activities, aes(x=aupr)) + 
    geom_histogram(color="black", fill="darkorange")  + 
    # geom_density(alpha=.1, fill="orange") +
    geom_vline(aes(xintercept=min(ligand_activities %>% top_n(10, aupr) %>% pull(aupr))), color="red", linetype="dashed", size=1) + 
    labs(x="ligand activity (PCC)", y = "# ligands") +
    theme_classic()
  p_hist_lig_activity
  
  # active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 100) %>% bind_rows()
  # 
  # active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
  # 
  # order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  # order_targets = active_ligand_target_links_df$target %>% unique()
  # vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  # 
  # p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot(senderName,paste0("Differentially expressed genes in ",receiverName), color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))
  # 
  # p_ligand_target_network
  
  
  # get the ligand-receptor network of the top-ranked ligands
  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  # get the weights of the ligand-receptor interactions as used in the NicheNet model
  # weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
  # weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds?download=1"))
  lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  # convert to a matrix
  lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  # perform hierarchical clustering to order the ligands and receptors
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot(senderName,receiverName, color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
  p_ligand_receptor_network
  
  ligand_pearson_matrix = ligand_activities %>% dplyr::select(aupr) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("aupr")
  p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot(senderName,"Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "aupr \ntarget gene prediction ability)")
  p_ligand_pearson
  
  return(list(ligandActivityPlot=p_hist_lig_activity,ligandActivityTable=ligand_activities,ligandReceptor=p_ligand_receptor_network,ligandPearson=p_ligand_pearson))
  
  
}

# Function to produce circos plot from the outputs of a nichenet analysis
# Code adapted from multinichenetR circos plot by Emily Johnson
# https://github.com/saeyslab/multinichenetr/blob/main/R/plotting.R
# circos_links_oi should be a dataframe containing the columns, ligand, target, weight, id, sender & receiver
# id should be a column that serves as a unique identifier for each interaction 

produceCircosPlot = function(circos_links_oi, colors_sender, colors_receiver, title=""){
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("circlize")
  
  
  # Link each cell type to a color
  grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
  grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
  
  
  # Make the plot for condition of interest - title of the plot
  title = title
  
  
  # Rename 'target' column to receptor
  circos_links = circos_links_oi
  df = circos_links
  
  df$id <- 1:nrow(df)
  
  
  # Each pair of ligand-receptors needs to be unique for circos plot to work
  # Code to make each pair unique by adding spaces after name
  ligand.uni = unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i = df[df$ligand == ligand.uni[i], ]
    sender.uni = unique(df.i$sender)
    for (j in 1:length(sender.uni)) {
      df.i.j = df.i[df.i$sender == sender.uni[j], ]
      df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
      df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
    }
  }
  receptor.uni = unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i = df[df$receptor == receptor.uni[i], ]
    receiver.uni = unique(df.i$receiver)
    for (j in 1:length(receiver.uni)) {
      df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
      df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
      df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
    }
  }
  
  intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
  
  while(length(intersecting_ligands_receptors) > 0){
    df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
    df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
    df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
    df = dplyr::bind_rows(df_unique, df_duplicated)
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
  }
  
  circos_links = df
  
  # Link ligands/Receptors to the colors of senders/receivers
  circos_links = circos_links %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
  links_circle = circos_links %>% dplyr::distinct(ligand,receptor, Sum)
  ligand_color = circos_links %>% dplyr::distinct(ligand,color_ligand_type)
  grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
  receptor_color = circos_links %>% dplyr::distinct(receptor,color_receptor_type)
  grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
  grid_col =c(grid_ligand_color,grid_receptor_color)
  
  # give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
  transparency = circos_links %>% dplyr::mutate(weight =(Sum-min(Sum))/(max(Sum)-min(Sum))) %>% dplyr::mutate(transparency = (1-weight)/1.2) %>% .$transparency
  
  # Define order of the ligands and receptors and the gaps
  ligand_order = circos_links_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    ligands = circos_links %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
  }) %>% unlist()
  
  receptor_order = circos_links_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    receptors = circos_links %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
  }) %>% unlist()
  
  order = c(ligand_order,receptor_order)
  
  width_same_cell_same_ligand_type = 0.275
  width_different_cell = 3
  width_ligand_receptor = 9
  width_same_cell_same_receptor_type = 0.275
  
  sender_gaps = circos_links_oi$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    sector = rep(width_same_cell_same_ligand_type, times = (circos_links %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
    gap = width_different_cell
    return(c(sector,gap))
  }) %>% unlist()
  sender_gaps = sender_gaps[-length(sender_gaps)]
  
  receiver_gaps = circos_links_oi$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    sector = rep(width_same_cell_same_receptor_type, times = (circos_links %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
    gap = width_different_cell
    return(c(sector,gap))
  }) %>% unlist()
  receiver_gaps = receiver_gaps[-length(receiver_gaps)]
  
  gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
  
  if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
    warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
  }
  
  links_circle$weight <- links_circle$Sum
  links_circle$Sum[links_circle$Sum == 0] = 0.01
  circos.clear()
  circos.par(gap.degree = gaps)
  chordDiagram(links_circle,
               directional = 1,
               order=order,
               link.sort = TRUE,
               link.decreasing = TRUE,
               grid.col = grid_col,
               transparency = transparency,
               diffHeight = 0.0075,
               direction.type = c("diffHeight", "arrows"),
               link.visible = links_circle$Sum > 0.01,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.175),
               link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1,
               reduce = 0,
               scale = TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
  }, bg.border = NA) #
  
  title(title)
  p_circos = recordPlot()
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  # grid_col_all = c(colors_receiver, colors_sender)
  legend = ComplexHeatmap::Legend(at = circos_links_oi$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_receiver[circos_links_oi$receiver %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))
  
  legend = ComplexHeatmap::Legend(at = circos_links_oi$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_sender[circos_links_oi$sender %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))
  
  p_legend = grDevices::recordPlot()
  
  circos_plot <- list(circos = p_circos,
                      legend = p_legend)
  
  return(circos_plot)
}
