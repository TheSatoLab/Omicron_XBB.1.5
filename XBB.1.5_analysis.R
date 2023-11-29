library(Biostrings)
library(ape)
library(phangorn)
library(tidyverse)
library(data.table)
library(EnvStats) # Rosner's test 
library(ggplot2)
library(ggtree)
library(ggnewscale)
library(RColorBrewer)

### Change working directory
dir <- "/Users/chainorato/Desktop/XBB.1.5_project/"
setwd(dir)

#################################################################################################### Analyses done in the revised manuscript ####################################################################################################

############################################# Retrieve metadata of XBB subvariants #############################################

pango_note <- read.delim("231024_lineage_notes.txt", header=TRUE, sep="\t")
rownames(pango_note) <- pango_note$Lineage
XBB_start_index <- which(pango_note$Lineage == "XBB")
XBB_end_index <-which(pango_note$Lineage == "XBB.9")
pango_note.XBB_filtered.lineage <- pango_note[XBB_start_index:XBB_end_index,]$Lineage # Checked, all XBB variant name identified by PANGO included

metadata <- fread("metadata_tsv_2023_01_19/metadata.tsv", header=T, sep="\t", quote="", check.names=T)
nextclade <- fread("metadata_tsv_2023_01_19/nextclade.tsv", header=T, sep="\t", quote="", check.names=T)
nextclade <- nextclade %>% arrange(index) %>% mutate(seqName=str_split(seqName, "\\|", simplify=T)[,1])
nextclade <- nextclade %>% mutate(Accession.ID = metadata$Accession.ID)

metadata <- metadata %>% mutate(Nextclade_pango = nextclade$Nextclade_pango)
table((metadata %>% filter(str_detect(Pango.lineage, "XBB")))$Pango.lineage)

# i) only ‘original passage’ sequences
# iii) host labelled as ‘Human’
# iv) sequence length above 28,000 base pairs
# v) proportion of ambiguous bases below 2%.
metadata.filtered <- metadata %>%
  distinct(Accession.ID,.keep_all=T) %>% distinct(Virus.name,.keep_all=T) %>%
  filter(Host == "Human",
         !N.Content > 0.02 | is.na(N.Content),
         str_length(Collection.date) == 10,
         Sequence.length > 28000,
         Passage.details.history == "Original",
         !str_detect(Additional.location.information,"[Qq]uarantine")
  )

metadata.filtered <- metadata.filtered %>%
  mutate(Collection.date = as.Date(Collection.date),
         region = str_split(Location," / ",simplify = T)[,1],
         country = str_split(Location," / ",simplify = T)[,2],
         state = str_split(Location," / ",simplify = T)[,3])

metadata.filtered <- metadata.filtered %>% filter(Pango.lineage %in% pango_note.XBB_filtered.lineage, Nextclade_pango != '') %>% filter(Nextclade_pango %in% pango_note.XBB_filtered.lineage)
nextclade.filtered <- nextclade %>% filter(Accession.ID %in% metadata.filtered$Accession.ID)

# write.table(metadata.filtered, "metadata_tsv_2023_01_19/231101_filtered_metadata.tsv", col.names=T, row.names=F, sep="\t", quote=F)
# write.table(nextclade.filtered, "metadata_tsv_2023_01_19/231101_filtered_nextclade.tsv", col.names=T, row.names=F, sep="\t", quote=F)


############################################################ Sensitivity analysis ############################################################ 

metadata.filtered <- fread("metadata_tsv_2023_01_19/231101_filtered_metadata.tsv", header=T, sep="\t", quote="", check.names=T)
nextclade.filtered <- fread("metadata_tsv_2023_01_19/231101_filtered_nextclade.tsv", header=T, sep="\t", quote="", check.names=T)

set.seed(49091)
 
for (i in 1:20) {
  metadata.filtered.phylogeny <- metadata.filtered %>% group_by(country,Pango.lineage) %>% slice_sample(n = 10)
  # write.table(metadata.filtered.phylogeny, paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_",  i, ".tsv", sep=""), col.names=T, row.names=F, sep="\t", quote=F)
}

########## 1. Remove taxa that may cause long branch attraction using Rosner's test ##########

finished_ten_entry <- c(2,5,7,8,10,12,15,16,17,20) # We picked only the first ten trees completely reconstructed due to the revision submission deadline

for (i in 1:20) {
  if (i %in% finished_ten_entry) {
    metadata.filtered.phylogeny <- fread(paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", i, ".tsv", sep=""), header=T, sep="\t", quote="", check.names=T)
    metadata.filtered.phylogeny <- metadata.filtered.phylogeny %>% mutate(Virus.name = str_replace_all(Virus.name, " ", "_"))
    xbb_tree <- read.tree(paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", i, ".fasta.edited.aln.trimal.treefile", sep=""))
    
    fasta_file <- readDNAStringSet(paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", i, ".fasta.edited.aln.trimal", sep=""), "fasta")
    fasta_seq_name <- names(fasta_file)
    fasta_seq <- paste(fasta_file)
    data_fasta <- data.frame(label=fasta_seq_name, seq=fasta_seq)
    data_fasta$label <- gsub(" ", "_", data_fasta$label)
    data_fasta_filtered <- data_fasta
    
    loop_index <- TRUE
    while (loop_index == TRUE) {
      ### Run Rosner's generalized extreme Studentized deviate test to remove branch length outliers
      xbb_tree.info.df <- ggtree(xbb_tree)$data
      xbb_tree.info.df <- xbb_tree.info.df[,c('label','branch.length')]
      xbb_tree.info.df$branch.length.log <- log(xbb_tree.info.df$branch.length+1,10)
      
      rosner_test <- rosnerTest(xbb_tree.info.df$branch.length.log, k=10, alpha=1e-4)
      rosner_test
      
      if (TRUE %in% rosner_test$all.stats$Outlier) {
        ### Filter outlier taxa out and write output FASTA file for phylogenetic tree reconstruction
        xbb_tree <- drop.tip(xbb_tree, xbb_tree.info.df[rosner_test$all.stats[which(rosner_test$all.stats$Outlier==TRUE),]$Obs.Num,]$label)
        data_fasta_filtered <- data_fasta_filtered[which(!data_fasta_filtered$label %in% xbb_tree.info.df[rosner_test$all.stats[which(rosner_test$all.stats$Outlier==TRUE),]$Obs.Num,]$label),]
        if (all(xbb_tree.info.df[rosner_test$all.stats[which(rosner_test$all.stats$Outlier==TRUE),]$Obs.Num,]$label %in% metadata.filtered.phylogeny$Virus.name)==FALSE) {loop_index <- FALSE} # Solve loop problem
      } else {
        loop_index <- FALSE
    }
  }
  
  data_fasta_filtered$label <- paste('>', data_fasta_filtered$label, sep="")
  # write.table(data_fasta_filtered, paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", i, ".fasta.edited.aln.trimal.outlier_filtered", sep=""), col.names=F, row.names=F, sep="\n", quote=F)
  }
}

########## 1.5 Run IQTree again ##########

########## 2. Create the FASTA file of fully aligned genomic sequences of SARS-CoV-2 for ancestral sequence reconstruction ##########

for (i in 1:20) {
  if (i %in% finished_ten_entry) {
    metadata.filtered.phylogeny <- fread(paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", i, ".tsv", sep=""), header=T, sep="\t", quote="", check.names=T)
    metadata.filtered.phylogeny <- metadata.filtered.phylogeny %>% mutate(Virus.name = str_replace_all(Virus.name, " ", "_"))
    xbb_tree <- read.tree(paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", i, ".fasta.edited.aln.trimal.outlier_filtered.treefile", sep=""))

    fasta_file <- readDNAStringSet(paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", i, ".fasta.edited.aln.full", sep=""), "fasta")
    fasta_seq_name <- names(fasta_file)
    fasta_seq <- paste(fasta_file)
    data_fasta <- data.frame(label=fasta_seq_name, seq=fasta_seq)
    data_fasta$label <- gsub(" ", "_", data_fasta$label)
    data_fasta_filtered <- data_fasta %>% filter(label %in% xbb_tree$tip.label)
    data_fasta_filtered$label <- paste('>', data_fasta_filtered$label, sep="")
    # write.table(data_fasta_filtered, paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", i, ".fasta.edited.aln.full.outlier_filtered", sep=""), col.names=F, row.names=F, sep="\n", quote=F)
  }
}

########## 3. Plot ancestral nucleotide to the tree and determine state changes ##########

final.anc.tree.plot.l <- list()

min.num.descendants <- 20
min.prop.descendants <- 0.5

mut.interest.info.name <- "metadata_tsv_2023_01_19/mut.interest.info.txt"
mut.interest.info <- read.table(mut.interest.info.name, header=T)

root_node.v <- as.numeric(c("X",2358,"X","X",2354,"X",2532,2357,"X",2357,"X",2364,"X","X",2354,2896,2365,"X","X",2397))

for (a in 1:20) {
  if (a %in% finished_ten_entry) {
    metadata.filtered.phylogeny <- fread(paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", a, ".tsv", sep=""), header=T, sep="\t", quote="", check.names=T)
    metadata.filtered.phylogeny <- metadata.filtered.phylogeny %>% mutate(Virus.name = str_replace_all(Virus.name, " ", "_"))
    metadata.filtered.phylogeny <- metadata.filtered %>% mutate(Pango.lineage = ifelse(Pango.lineage=="XBB.1.4.1", "XBB.1.4",
                                                                                ifelse(Pango.lineage=="XBB.3.1", "XBB.3",
                                                                                ifelse(Pango.lineage=="XBB.4.1", "XBB.4", Pango.lineage))))

    ### Read the tree after removing outliers
    xbb_tree <- read.tree(paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", a, ".fasta.edited.aln.trimal.outlier_filtered.treefile", sep=""))

    ### Read tip sequence file
    asr.tip <- readDNAStringSet(paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", a, ".fasta.edited.aln.full.outlier_filtered", sep=""), "fasta")
    asr.tip.fasta_seq_name <- names(asr.tip)
    asr.tip.fasta_seq <- paste(asr.tip)
    asr.tip.df <- data.frame(label=asr.tip.fasta_seq_name, seq=asr.tip.fasta_seq)
    asr.tip.df <- cbind(asr.tip.df, do.call(rbind, strsplit(asr.tip.df$seq, "")))
    
    site.interested <- c('22316','22317','22318','23018','23019','23020','27915','27916','27917')
    asr.tip.df$"Spike:252" <- apply(asr.tip.df[,site.interested[1:3]], 1, paste, collapse="")
    asr.tip.df$"Spike:486" <- apply(asr.tip.df[,site.interested[4:6]], 1, paste, collapse="")
    asr.tip.df$"NS8:8" <- apply(asr.tip.df[,site.interested[7:9]], 1, paste, collapse="")
    asr.tip.df$add.label <- apply(asr.tip.df[,c('Spike:252','Spike:486','NS8:8')], 1, paste, collapse="/")
    asr.tip.df.selected <- asr.tip.df[,c('label','add.label')]
    asr.tip.df.selected <- asr.tip.df.selected %>% arrange(factor(label, levels = xbb_tree$tip.label))
    asr.tip.df.selected <- asr.tip.df.selected %>% mutate(node.label = xbb_tree$tip.label)
    asr.tip.df.selected$new.label <- apply(asr.tip.df.selected[,c("label","add.label")], 1, paste, collapse="|")
    
    ### Read ancestral sequence file
    asr.iqtree <- readDNAStringSet(paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", a, ".fasta.edited.aln.full.outlier_filtered.asr.fasta", sep=""), "fasta")
    asr.iqtree.fasta_seq_name <- names(asr.iqtree)
    asr.iqtree.fasta_seq <- paste(asr.iqtree)
    asr.iqtree.df <- data.frame(label=asr.iqtree.fasta_seq_name, seq=asr.iqtree.fasta_seq)
    asr.iqtree.df <- cbind(asr.iqtree.df, do.call(rbind, strsplit(asr.iqtree.df$seq, "")))
    
    site.interested <- c('22316','22317','22318','23018','23019','23020','27915','27916','27917')
    asr.iqtree.df$"Spike:252" <- apply(asr.iqtree.df[,site.interested[1:3]], 1, paste, collapse="")
    asr.iqtree.df$"Spike:486" <- apply(asr.iqtree.df[,site.interested[4:6]], 1, paste, collapse="")
    asr.iqtree.df$"NS8:8" <- apply(asr.iqtree.df[,site.interested[7:9]], 1, paste, collapse="")
    asr.iqtree.df$add.label <- apply(asr.iqtree.df[,c('Spike:252','Spike:486','NS8:8')], 1, paste, collapse="/")
    asr.iqtree.df.selected <- asr.iqtree.df[,c('label','add.label')]
    asr.iqtree.df.selected <- asr.iqtree.df.selected %>% mutate(Node = str_replace_all(label, "Node", "")) %>% arrange(as.numeric(Node))
    asr.iqtree.df.selected <- asr.iqtree.df.selected %>% mutate(node.label = xbb_tree$node.label)
    asr.iqtree.df.selected$new.label <- apply(asr.iqtree.df.selected[,c("node.label","add.label")], 1, paste, collapse="|")

    ### Root the tree
    xbb_tree$tip.label <- asr.tip.df.selected$new.label
    xbb_tree$node.label <- asr.iqtree.df.selected$new.label
    xbb_tree <- phytools::reroot(xbb_tree, node=root_node.v[a], position=0.5*xbb_tree$edge.length[which(xbb_tree$edge[,2]==root_node.v[a])]) # Root at midpoint of specified and earlier nodes
    
    #detect transition branches
    lineage.interest <- "XBB"
    xbb_tree.info.df <- ggtree(xbb_tree)$data
    
    # Count number of all descendant tips under each node
    count_descendants <- function(node) {
      num_descendants <- lengths(Descendants(xbb_tree, node, type="tips"))
      return(num_descendants)}
    
    max.date <- max(as.Date(metadata.filtered.phylogeny$Collection.date))
    max.x <- max(xbb_tree.info.df$x)
    
    xbb_tree.info.df$num.descendants <- as.numeric(count_descendants(xbb_tree.info.df$node))
    mut.interest.v <- c("Spike:G252V","Spike:S486P","NS8:G8stop")

    #detect transition branches
    node.transition.df <- data.frame()

    for (j in 1:length(mut.interest.v)) {
      #extract transition branch
      xbb_tree.info.df.interest <- xbb_tree.info.df %>% mutate(state = str_split(str_split(label, "\\|", simplify=TRUE)[,2], "/", simplify=TRUE)[,j]) %>% select(parent, node, num.descendants, state) 
      xbb_tree.info.df.interest.parent <- xbb_tree.info.df.interest %>% select(node,state) %>% dplyr::rename(parent = node, state.parent = state)
      xbb_tree.info.df.interest.merged <- merge(xbb_tree.info.df.interest,xbb_tree.info.df.interest.parent,by="parent")
      xbb_tree.info.df.interest.merged.change <- xbb_tree.info.df.interest.merged %>% filter(!str_detect(state, "N|-"), !str_detect(state.parent, "N|-")) %>% filter(state != "", state.parent != "") %>% filter(state != state.parent)
      xbb_tree.info.df.interest.merged.change <- xbb_tree.info.df.interest.merged.change %>% filter(state %in% names(GENETIC_CODE), state.parent %in% names(GENETIC_CODE))
      xbb_tree.info.df.interest.merged.change <- xbb_tree.info.df.interest.merged.change %>% mutate(state = as.character(translate(DNAStringSet(state))), state.parent = as.character(translate(DNAStringSet(state.parent))))
      xbb_tree.info.df.interest.merged.change$state <- ifelse(xbb_tree.info.df.interest.merged.change$state == "*", "stop", xbb_tree.info.df.interest.merged.change$state)
      xbb_tree.info.df.interest.merged.change$state.parent <- ifelse(xbb_tree.info.df.interest.merged.change$state.parent == "*", "stop", xbb_tree.info.df.interest.merged.change$state.parent)
      xbb_tree.info.df.interest.merged.change <- xbb_tree.info.df.interest.merged.change %>% mutate(mut = paste(str_split(mut.interest.v[j], "\\:", simplify=TRUE)[,1], ":", state.parent, gsub("[A-Z]", "", str_split(mut.interest.v[j], "\\:", simplify=TRUE)[,2], ignore.case=TRUE), state, sep=""))
      xbb_tree.info.df.interest.merged.change <- xbb_tree.info.df.interest.merged.change %>% filter(mut %in% mut.interest.v)
      
      transition.node.info.df <- data.frame()
      for (i in 1:nrow(xbb_tree.info.df.interest.merged.change)) {
        node.interest <- xbb_tree.info.df.interest.merged.change$node[i]
        date.interest <- max.date - round(365 * (max.x - xbb_tree.info.df$x[node.interest]))

        tip.v <- Descendants(xbb_tree, node.interest, type = "tips")[[1]]
        tip.df.with_mut <- xbb_tree.info.df.interest %>% filter(state %in% names(GENETIC_CODE)) %>% mutate(state = as.character(translate(DNAStringSet(state))))
        tip.df.with_mut$state <- ifelse(tip.df.with_mut$state == "*", "stop", tip.df.with_mut$state)
        tip.df.with_mut <- tip.df.with_mut %>% filter(state == str_split(str_split(mut.interest.v[j], "\\:", simplify=TRUE)[,2], gsub("[A-Z]", "", str_split(mut.interest.v[j], "\\:", simplify=TRUE)[,2], ignore.case=TRUE), simplify=TRUE)[,2], node %in% tip.v)
        num.tip.with_mut <-  nrow(tip.df.with_mut)

        temp.df <- data.frame(date = date.interest, num.tip.with_mut)
        transition.node.info.df <- rbind(transition.node.info.df, temp.df)
      }

      xbb_tree.info.df.interest.merged.change <- cbind(xbb_tree.info.df.interest.merged.change, transition.node.info.df)
      xbb_tree.info.df.interest.merged.change <- xbb_tree.info.df.interest.merged.change %>% mutate(prop.num.tip.with_mut = num.tip.with_mut / num.descendants)
 
      # filtering transition nodes
      xbb_tree.info.df.interest.merged.change.filtered <- xbb_tree.info.df.interest.merged.change %>% filter(num.tip.with_mut >= min.num.descendants, prop.num.tip.with_mut >= min.prop.descendants)
      node.transition.df <- rbind(node.transition.df, xbb_tree.info.df.interest.merged.change.filtered)
    }
    
    #save transition branch info
    out.name <- paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", a, "_transition_branch.", lineage.interest,".txt", sep="")
    # write.table(node.transition.df, out.name, col.names=T, row.names=F, sep="\t", quote=F)
    
    node.df <- node.transition.df %>% select(node, mut)

    color.mut.v <- brewer.pal(9, "Set1")[c(4,6,8)]
    names(color.mut.v) <- c("Spike:G252V", "NS8:G8stop", "Spike:S486P")

    xbb_tree <- read.tree(paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_", a, ".fasta.edited.aln.trimal.outlier_filtered.treefile", sep=""))
    xbb_tree <- phytools::reroot(xbb_tree, node=root_node.v[a], position=0.5*xbb_tree$edge.length[which(xbb_tree$edge[,2]==root_node.v[a])]) # Root at midpoint of specified and earlier nodes

    g <- ggtree(xbb_tree, size=0.25, aes(color=Pango.lineage)) + labs(shape="PANGO lineage", color="PANGO lineage")
    g <- g + geom_text2(aes(subset=node %in% node.df$node, label=label), size=2, hjust=1.3)
    g <- g %<+% metadata.filtered.phylogeny + geom_tippoint(aes(subset = Pango.lineage %in% c("XBB.1","XBB.1.5"), color = Pango.lineage), alpha=0.5)
    g <- g + scale_color_discrete(na.value="black") + new_scale_color()
    g <- g %<+% node.df + geom_point2(shape=18, size=7, aes(color=mut))
    g <- g + scale_color_manual("Mutation", values=color.mut.v, breaks=names(color.mut.v), na.value = NA)
    g <- g + scale_x_continuous(limits=c(0,0.001))
    g <- g + geom_treescale(width=1e-4)
    
    if (a == 20) {
      h <- ggtree(xbb_tree, size=0.25, aes(color=Pango.lineage)) + labs(shape="PANGO lineage", color="PANGO lineage") #+ geom_tiplab(aes(label = Pango.lineage), size=2, align=FALSE, linesize=.5)
      h <- h + geom_text2(aes(subset=node %in% node.df$node, label=label), size=2, hjust=1.3)
      h <- h %<+% metadata.filtered.phylogeny #+ geom_tippoint(aes(color = Pango.lineage), alpha=0.5)
      h <- h + scale_color_discrete(na.value="black") + new_scale_color()
      h <- h %<+% node.df + geom_point2(shape=18, size=7, aes(color=mut))
      h <- h + scale_color_manual("Mutation", values=color.mut.v, breaks=names(color.mut.v), na.value = NA)
      h <- h + geom_treescale(width=1e-4)
      
      h  <- h  + geom_cladelab(node=3271, label="XBB.1")
      h  <- h  + geom_cladelab(node=4167, label="XBB.1.1")
      h  <- h  + geom_cladelab(node=4363, label="XBB.1.2")
      h  <- h  + geom_cladelab(node=4561, label="XBB.1.3")
      h  <- h  + geom_cladelab(node=3882, label="XBB.1.4")
      h  <- h  + geom_cladelab(node=3401, label="XBB.1.5")
      h  <- h  + geom_cladelab(node=2951, label="XBB.2")
      h  <- h  + geom_cladelab(node=2495, label="XBB.3")
      h  <- h  + geom_cladelab(node=2769, label="XBB.4")
      h  <- h  + geom_cladelab(node=4635, label="XBB.5")
      h  <- collapse(h , 4167, 'max', fill='green') #%>% expand(3142)		# XBB.1.1
      h  <- collapse(h , 4363, 'max', fill='black') #%>% expand(3464) 		# XBB.1.2
      h  <- collapse(h , 4561, 'max', fill='red') #%>% expand(4030) 		# XBB.1.3
      h  <- collapse(h , 3882, 'max', fill='orange') #%>% expand(3471)		# XBB.1.4 + XBB.1.4.1
      h  <- collapse(h , 3401, 'max', fill='steelblue') #%>% expand(2519)	# XBB.1.5
      h  <- collapse(h , 2951, 'max', fill='blue') #%>% expand(4050)		# XBB.2
      h  <- collapse(h , 2495, 'max', fill='yellow') #%>% expand(4641)		# XBB.3 + XBB.3.1
      h  <- collapse(h , 2769, 'max', fill='pink') #%>% expand(4487)		# XBB.4
      h  <- collapse(h , 4635, 'max', fill='brown') #%>% expand(4568)		# XBB.5

      # pdf.name <- paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_final_anc_tree_", a, "_main.pdf", sep="")
      # pdf(pdf.name, width=8, height=15)
      # plot(h)
      # dev.off()
    }
     
    # pdf.name <- paste("metadata_tsv_2023_01_19/sensitivity_analysis/231101_filtered_metadata_phylogeny_final_anc_tree_", a, ".pdf", sep="")
    # pdf(pdf.name, width=8, height=12)
    # plot(g)
    # dev.off()

    final.anc.tree.plot.l[[pdf.name]] <- g
  }
}







#################################################################################################### Analyses done in the original manuscript ####################################################################################################

########## 1. Remove taxa that may cause long branch attraction using Rosner's test ##########

### Read metadata file
xbb_metadata <- read.delim("20230202_XBB_phylogeny/20230202_XBB_selected_metadata.tsv", header=TRUE, sep="\t")
xbb_metadata$Virus.name <- str_replace_all(xbb_metadata$Virus.name, " ", "_")
xbb_metadata$Location <- str_split(xbb_metadata$Location, " / ", simplify=TRUE)[,1]
rownames(xbb_metadata) <- xbb_metadata$Virus.name

### Read the preliminary tree
xbb_tree <- read.tree("20230202_XBB_phylogeny/20230202_XBB_selected_sequences.fasta.edited.aln.trimal.trimmed.treefile")
xbb_tree <- root(xbb_tree, "NC_045512.2", resolve.root=TRUE)
is.rooted(xbb_tree)

### Plot the preliminary tree
xbb_tree_plot <- ggtree(xbb_tree, size=0.25) + theme_tree2() + ggtitle("XBB lineages") + geom_text2(aes(subset=!isTip, label=node), size=2, hjust=1.3) + scale_shape_manual(values=seq(5,18)) #+ geom_tiplab(size=2, align=TRUE, linesize=.5) # use to guide rooting
xbb_tree_plot <- xbb_tree_plot %<+% xbb_metadata + geom_tippoint(aes(shape=Pango.lineage, color=Pango.lineage), alpha=0.75)
xbb_tree_plot

### Gather FASTA sequences of all viruses used in preliminary phylogenetic tree reconstruction
fasta_file <- readDNAStringSet("20230202_XBB_phylogeny/20230202_XBB_selected_sequences.fasta.edited.aln.trimal.trimmed", "fasta")
fasta_seq_name <- names(fasta_file)
fasta_seq <- paste(fasta_file)
data_fasta <- data.frame(label=fasta_seq_name, seq=fasta_seq)
data_fasta$label <- gsub(" ", "_", data_fasta$label)

### Run Rosner's generalized extreme Studentized deviate test to remove branch length outliers
xbb_tree.info.df <- ggtree(xbb_tree)$data
xbb_tree.info.df <- xbb_tree.info.df[,c('label','branch.length')]
xbb_tree.info.df$branch.length.log <- log(xbb_tree.info.df$branch.length+1,10)

log_branch_length_boxplot <- ggplot(xbb_tree.info.df, aes(branch.length.log)) + geom_boxplot()
log_branch_length_boxplot

rosner_test <- rosnerTest(xbb_tree.info.df$branch.length.log, k=10, alpha=1e-4)
rosner_test

### Filter outlier taxa out and write output FASTA file for phylogenetic tree reconstruction
data_fasta_filtered <- data_fasta[which(!data_fasta$label %in% xbb_tree.info.df[rosner_test$all.stats[which(rosner_test$all.stats$Outlier==TRUE),]$Obs.Num,]$label),]
data_fasta_filtered$label <- paste('>', data_fasta_filtered$label, sep="")
#write.table(data_fasta_filtered, "20230202_XBB_selected_sequences.fasta.edited.aln.trimal.trimmed.outlier_filtered", col.names=F, row.names=F, sep="\n", quote=F)

### Initially check the tree whose outliers were dropped
xbb_tree <- drop.tip(xbb_tree, xbb_tree.info.df[rosner_test$all.stats[which(rosner_test$all.stats$Outlier==TRUE),]$Obs.Num,]$label)
xbb_tree_plot <- ggtree(xbb_tree, size=0.25) + theme_tree2() + ggtitle("XBB lineages") + geom_text2(aes(subset=!isTip, label=node), size=2, hjust=1.3) + scale_shape_manual(values=seq(5,18)) #+ geom_tiplab(size=2, align=TRUE, linesize=.5) # use to guide rooting
xbb_tree_plot <- xbb_tree_plot %<+% xbb_metadata + geom_tippoint(aes(shape=Pango.lineage, color=Pango.lineage), alpha=0.75)
xbb_tree_plot

########## 1.5 Run IQTree again ##########

########## 2. Check the final tree and plot it with metadata ##########

### Read the tree after removing outliers
xbb_tree <- read.tree("20230202_XBB_phylogeny/20230202_XBB_selected_sequences.fasta.edited.aln.trimal.trimmed.outlier_filtered.treefile")
xbb_tree <- phytools::reroot(xbb_tree, node=3653, position=0.5*xbb_tree$edge.length[which(xbb_tree$edge[,2]==3653)]) # Root at midpoint of node 3653 and node 3652 # similar root to time-calibrated tree, dividing XBB and other related lineages from XBB.2
is.rooted(xbb_tree)

xbb_metadata <- xbb_metadata[which(rownames(xbb_metadata) %in% xbb_tree$tip.lab),]

### Plot the final tree
xbb_tree_plot <- ggtree(xbb_tree, size=0.25) + theme_tree2() + scale_shape_manual(values=seq(1,14)) + labs(shape="Pango lineages", col="Pango lineages") + geom_text2(aes(subset=!isTip, label=label), size=2, hjust=1.3, alpha=0.75) #+ geom_hilight(node=2522, fill="steelblue", alpha=0.05)  # + geom_tiplab(size=2, align=FALSE, linesize=.5) # use to guide rooting
xbb_tree_plot <- xbb_tree_plot %<+% xbb_metadata + geom_tippoint(aes(shape=Pango.lineage, color=Pango.lineage), size=1.25, alpha=0.75)
xbb_tree_plot

### Read information of Spike and NS8 alignment files
##### Spike:G252 --> Position 282 of the alignment
##### Spike:F486 --> Position 516 of the alignment
spike_prot_info <- readAAStringSet("20230202_XBB_phylogeny/20230202_XBB_selected_Spike_mafft.output", "fasta")
spike_prot_seq_name <- names(spike_prot_info)
spike_prot_seq <- paste(spike_prot_info)
spike_prot_df <- data.frame(accession_id=spike_prot_seq_name, fasta_seq=spike_prot_seq)
spike_prot_df <- cbind(spike_prot_df, do.call(rbind, strsplit(spike_prot_df$fasta_seq, "")))

##### NS8:G8 --> Position 8 of the alignment
NS8_prot_info <- readAAStringSet("20230202_XBB_phylogeny/20230202_XBB_selected_NS8_mafft.output", "fasta")
NS8_prot_seq_name <- names(NS8_prot_info)
NS8_prot_seq <- paste(NS8_prot_info)
NS8_prot_df <- data.frame(accession_id=NS8_prot_seq_name, fasta_seq=NS8_prot_seq)
NS8_prot_df <- cbind(NS8_prot_df, do.call(rbind, strsplit(NS8_prot_df$fasta_seq, "")))

xbb_metadata$"Spike:G252V" <- 0
xbb_metadata$"Spike:F486S" <- 0
xbb_metadata$"Spike:F486P" <- 0
xbb_metadata$"NS8:G8stop" <- 0

for (a in c("Spike:G252V", "Spike:F486S", "Spike:F486P", "NS8:G8stop")) {
  prot_name <- str_split(a, ":", simplify=TRUE)[,1]
  position <- gsub("[A-Z]", "", str_split(a, ":", simplify=TRUE)[,2], ignore.case=TRUE)
  if (prot_name=="Spike") {
    if (position=="252") {position <- "282"}
    if (position=="486") {position <- "516"}
    xbb_metadata[, a] <- spike_prot_df[match(xbb_metadata$Accession.ID, spike_prot_df$accession_id), position]
  }
  if (prot_name=="NS8") {
    if (position=="8") {position <- "8"}
    if (position=="486") {position <- "516"}
    xbb_metadata[, a] <- NS8_prot_df[match(xbb_metadata$Accession.ID, NS8_prot_df$accession_id), position]
  }
}

########## 3. Plot ancestral node to the tree ##########

### Need to create metadata.mut_long.txt at first using the script ... from Ito-san

min.num.descendants <- 5
min.prop.descendants <- 0

mut.interest.info.name <- "20230202_XBB_phylogeny/mut.interest.info.txt"
mut.interest.info <- read.table(mut.interest.info.name, header=T)

metadata <- xbb_metadata

mut.info.name <- "20230202_XBB_phylogeny/metadata.mut_long.txt"
mut.info <- read.table(mut.info.name, header=T,sep="\t", check.names=T)

#preprocessing
mut.info$mut.mod <- mut.info$mut
mut.info.spread <- mut.info %>% select(Id, mut.mod) %>% mutate(value = 1) %>% spread(key=mut.mod, value=value)
mut.info.spread <- mut.info.spread %>% select(Id, Spike_F486S, Spike_G252V, NS8_G8stop, Spike_F486P)
mut.info.spread <- mut.info.spread[which(mut.info.spread$Id %in% xbb_tree$tip.lab),]

#detect transition branches
lineage.interest <- "XBB"
tree <- xbb_tree
tree.info.df <- ggtree(tree)$data

# Count number of all descendant tips under each node
count_descendants <- function(node) {
  num_descendants <- lengths(Descendants(tree, node, type="tips"))
  return(num_descendants)}

max.date <- max(as.Date(metadata$Collection.date))
max.x <- max(tree.info.df$x)

tree.info.df$num.descendants <- as.numeric(count_descendants(tree.info.df$node))
tip.df <- data.frame(tip_Id = 1:length(tree$tip), Id = tree$tip)
tip.df.merged <- tip.df %>% left_join(mut.info.spread, by="Id")
tip.df.merged[is.na(tip.df.merged)] <- 0

mut.interest.v <- mut.interest.info %>% filter(lineage==lineage.interest) %>% pull(mut)

mut.info.mat <- tip.df.merged[,mut.interest.v]

#detect transition branches
node.transition.df <- data.frame()
branch_state.l <- list()

for (mut.name in mut.interest.v) {
  #ancestral state reconstruction
  state.v <- mut.info.mat %>% pull(mut.name)
  names(state.v) <- tip.df.merged$Id
  for (a in 1:length(state.v)) {
    if (!is.na(metadata[names(state.v[a]), str_replace_all(mut.name, "_", ":")]) && metadata[names(state.v[a]), str_replace_all(mut.name, "_", ":")]=="X") {state.v[a] <- NA}} # ace can deal with uncertainty using NA
  
  fit.asr <- ace(state.v+1, tree, model="ER", type="discrete")
  
  state.mat <- fit.asr$lik.anc
  state.node.v <- state.mat[,1]
  state.node.v <- ifelse(state.node.v > 0.5, 0, 1)
  
  state.color <- c(state.v, state.node.v)
  
  branch_state.l[[mut.name]] <- state.color
  
  #extract transition branch
  tree.info.df.interest <- tree.info.df %>% select(parent, node, num.descendants) %>% mutate(state = state.color)
  
  tree.info.df.interest.parent <- tree.info.df.interest %>% select(node,state) %>% dplyr::rename(parent = node, state.parent = state)
  tree.info.df.interest.merged <- merge(tree.info.df.interest,tree.info.df.interest.parent,by="parent")
  tree.info.df.interest.merged.0to1 <- tree.info.df.interest.merged %>% filter(state==1, state.parent==0)
  
  transition.node.info.df <- data.frame()
  for (i in 1:nrow(tree.info.df.interest.merged.0to1)) {
    
    node.interest <- tree.info.df.interest.merged.0to1$node[i]
    date.interest <- max.date - round(365 * (max.x - tree.info.df$x[node.interest]))
    
    tip.v <- Descendants(tree, node.interest, type = "tips")[[1]]
    
    tip.df.with_mut <- tip.df.merged %>% filter(get({{mut.name}}) == 1, tip_Id %in% tip.v)
    
    num.tip.with_mut <-  nrow(tip.df.with_mut)
    
    mut.info.interest <- tip.df.with_mut %>% left_join(mut.info %>% filter(mut.mod == mut.name),by="Id")
    
    mut_type.major <- mut.info.interest %>% group_by(mut) %>% summarize(count = n()) %>% slice_max(count,n=1,with_ties = F) %>% pull(mut)
    
    temp.df <- data.frame(date = date.interest, mut_type = mut_type.major,num.tip.with_mut)
    transition.node.info.df <- rbind(transition.node.info.df,temp.df)
    
  }
  
  tree.info.df.interest.merged.0to1 <- cbind(tree.info.df.interest.merged.0to1,transition.node.info.df)
  tree.info.df.interest.merged.0to1 <- tree.info.df.interest.merged.0to1 %>% mutate(prop.num.tip.with_mut = num.tip.with_mut / num.descendants)
  print(tree.info.df.interest.merged.0to1)
  
  # filtering transition nodes
  tree.info.df.interest.merged.0to1.filtered <- tree.info.df.interest.merged.0to1 %>% filter(num.tip.with_mut >= min.num.descendants, prop.num.tip.with_mut >= min.prop.descendants) %>% mutate(mut = mut.name)
  
  node.transition.df <- rbind(node.transition.df,tree.info.df.interest.merged.0to1.filtered)
  
}

#save transition branch info
out.name <- paste("20230202_XBB_phylogeny/output/transition_branch.",lineage.interest,".txt",sep="")
# write.table(node.transition.df, out.name, col.names=T, row.names=F, sep="\t", quote=F)

branch_state.df <- as.data.frame(branch_state.l)
branch_state.df <- branch_state.df %>% mutate(total = apply(branch_state.df,1,sum), node = 1:nrow(branch_state.df))
node.df <- node.transition.df %>% select(node,mut)
node.df$mut <- str_replace_all(node.df$mut, "_", ":")

color.mut.v <- brewer.pal(9, "Set1")[c(3,4,6,8)]
names(color.mut.v) <- c("Spike:F486S", "Spike:G252V", "NS8:G8stop", "Spike:F486P")

g <- ggtree(tree, size=0.25) + geom_treescale(width=1e-4) + theme_tree2() + scale_shape_manual(values=seq(1,14)) + labs(shape="Pango lineages", col="Pango lineages") #+ geom_text2(aes(subset=!isTip, label=node), size=2, hjust=1.3, alpha=0.75) #+ geom_tiplab(size=2, align=FALSE, linesize=.5)
g <- g %<+% metadata + geom_tippoint(aes(shape=Pango.lineage, color=Pango.lineage), size=1.25, alpha=0.75)
g <- g + geom_cladelab(node=2377, label="XBB.1")
g <- g + geom_cladelab(node=3142, label="XBB.1.1")
g <- g + geom_cladelab(node=3464, label="XBB.1.2")
g <- g + geom_cladelab(node=4030, label="XBB.1.3")
g <- g + geom_cladelab(node=3471, label="XBB.1.4 & XBB.1.4.1")
g <- g + geom_cladelab(node=2519, label="XBB.1.5")
g <- g + geom_cladelab(node=4050, label="XBB.2")
g <- g + geom_cladelab(node=4641, label="XBB.3 & XBB.3.1")
g <- g + geom_cladelab(node=4487, label="XBB.4")
g <- g + geom_cladelab(node=4502, label="XBB.4 & XBB.4.1")
g <- g + geom_cladelab(node=4568, label="XBB.5")
g <- g + new_scale_color()
g <- g %<+% node.df + geom_point2(shape=18, size=7, aes(color=mut))
g <- g + scale_color_manual("Mutation", values=color.mut.v, breaks=names(color.mut.v), na.value = NA)
g <- collapse(g, 3142, 'max', fill='green') #%>% expand(3142)		# XBB.1.1
g <- collapse(g, 3464, 'max', fill='black') #%>% expand(3464) 		# XBB.1.2
g <- collapse(g, 4030, 'max', fill='red') #%>% expand(4030) 		# XBB.1.3
g <- collapse(g, 3471, 'max', fill='orange') #%>% expand(3471)		# XBB.1.4 + XBB.1.4.1
g <- collapse(g, 2519, 'max', fill='steelblue') #%>% expand(2519)	# XBB.1.5
g <- collapse(g, 4050, 'max', fill='blue') #%>% expand(4050)		# XBB.2
g <- collapse(g, 4641, 'max', fill='yellow') #%>% expand(4641)		# XBB.3 + XBB.3.1
g <- collapse(g, 4487, 'max', fill='pink') #%>% expand(4487)		# XBB.4
g <- collapse(g, 4502, 'max', fill='pink') #%>% expand(4502)		# XBB.4 + XBB.4.1
g <- collapse(g, 4568, 'max', fill='brown') #%>% expand(4568)		# XBB.5

# pdf("xbb.1.5_phylogeny.pdf", width=10, height=15)
# g
# dev.off()