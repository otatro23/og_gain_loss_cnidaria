library(dplyr)

if(file.exists("eggnog_annotations.Rdata")) {
  load("eggnog_annotations.Rdata")
} else {
  eggnog <- read.table("eggnog_keep.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE,)
  save(eggnog, file = "eggnog_annotations.Rdata")
}

load("../../orthodollo/loss_gain_pres_abs_H1.Rds")
gl <- lossgain.table
rm(lossgain.table)


gl <- as.data.frame(t(gl))
colnames(gl) <- as.character(1:ncol(gl))
gl$og <- rownames(gl)

load("eggnog_categories.Rdata") 


GO_setup <- function(NODE, DATA = GOs) {
  
  gains <- DATA$g_GOs[[NODE]]
  losses <- DATA$l_GOs[[NODE]]
  
  # take a row of the GOs dataframe and make count tables with each GO term and how many are in gains and losses (may need to do this separately and then cbind)
  if (all(is.na(gains))) {
    g_counts <- data.frame(GO_term = character(), gained = numeric())
  } else {
    g_counts <- as.data.frame(table(gains))
    colnames(g_counts) <- c("GO_term", "gained")
  }

  if (all(is.na(losses))) {
    l_counts <- data.frame(GO_term = character(), lost = numeric())
  } else {
    l_counts <- as.data.frame(table(losses))
    colnames(l_counts) <- c("GO_term", "lost")
  }
  #g_counts <- as.data.frame(table(DATA[NODE,2]))
  #colnames(g_counts) <- c("GO_term", "gained")
  #l_counts <- as.data.frame(table(DATA[NODE,3]))
  #colnames(l_counts) <- c("GO_term", "lost")
  
  GO_counts <- merge(g_counts, l_counts, by = "GO_term", all = TRUE)
  GO_counts$GO_term <- factor(GO_counts$GO_term, levels = GO_counts$GO_term)
  GO_counts$gained <- ifelse(is.na(GO_counts$gained), 0, GO_counts$gained)
  GO_counts$lost <- ifelse(is.na(GO_counts$lost), 0, GO_counts$lost)
  
  total_gains <- sum(GO_counts$gained)
  total_losses <- sum(GO_counts$lost)
  GO_counts$prop_gained <- GO_counts$gained / total_gains
  GO_counts$prop_lost <- GO_counts$lost / total_losses
  GO_counts$diff <- GO_counts$gained - GO_counts$lost
  GO_counts$prop_diff <- GO_counts$prop_gained - GO_counts$prop_lost
  GO_counts$sum <- GO_counts$gained + GO_counts$lost
  
  if(nrow(filter(GO_counts, gained != 0)) != 0) {
    # fishertest: IF EACH GO TERM IS GAINED MORE IN NODE THAN IN OTHER NODES
    # a = number of GO_term gains in anthozoa
    # b = number of other GO terms gained in anthozoa
    # c = number of GO_term gains in other nodes
    # d = number of other GO term gains in other nodes 
    b = as.numeric(length(DATA[[NODE,2]]))
    all_GO_gains <- unlist(DATA$g_GOs, recursive = TRUE)
    d = as.numeric(length(all_GO_gains))
    all_GO_gain_counts <- table(all_GO_gains)
    
    GO_counts <- GO_counts %>% rowwise() %>% 
      mutate(pval_gain = {
        term <- GO_term
        if (term %in% names(all_GO_gain_counts)) {
          a <- gained
          c <- all_GO_gain_counts[[term]]
          mat <- matrix(c(a, b - a, c, d - c), nrow = 2, byrow = TRUE)
          fisher.test(mat, alternative = "greater")$p.value
        } else {
          1
        }
      }) %>% 
      ungroup()
    GO_counts$padj_gain <- p.adjust(GO_counts$pval_gain, method = "bonferroni")
  }  
  if(nrow(filter(GO_counts, lost != 0)) != 0) {
    # fishertest: IF EACH GO TERM IS LOST MORE IN NODE THAN IN OTHER NODES
    # a = number of GO_term losses in anthozoa
    # b = number of other GO terms lost in anthozoa
    # c = number of GO_term losses in other nodes
    # d = number of other GO term losses in other nodes 
    b = as.numeric(length(DATA[[NODE,3]]))
    all_GO_losses <- unlist(DATA$l_GOs, recursive = TRUE)
    d = as.numeric(length(all_GO_losses))
    all_GO_loss_counts <- table(all_GO_losses)
    
    GO_counts <- GO_counts %>% rowwise() %>% 
      mutate(pval_loss = {
        term <- GO_term
        if (term %in% names(all_GO_loss_counts)) {
          a <- lost
          c <- all_GO_loss_counts[[term]]
          mat <- matrix(c(a, b - a, c, d - c), nrow = 2, byrow = TRUE)
          fisher.test(mat, alternative = "greater")$p.value
        } else {
          1
        }
      }) %>% 
      ungroup()
    GO_counts$padj_loss <- p.adjust(GO_counts$pval_loss, method = "bonferroni")
  }  
  
  
  return(GO_counts)
}

# Basal nodes 
print("processing node 86")
anth_GO <- GO_setup(86) # anthozoa 

print("processing node 87")
anth2_GO <- GO_setup(87) # anthozoa2 (hexacorals and ceriantheria)

print("processing node 84")
endomed_GO <- GO_setup(84) # endocnidozoa and medusozoa

print("processing node 83")
med_GO <- GO_setup(83) # medusazoa

print("processing node 129")
endo_GO <- GO_setup(129) # endocnidozoa

# clade gains/losses
print("processing node 130")
myx_GO <- GO_setup(130) # myxozoa

print("processing node 82")
hyd_GO <- GO_setup(82) # hydrozoa

print("processing node 132")
scs_GO <- GO_setup(132) # ancestor of staurozoa, cubozoa, and scyphozoa (scs)

print("processing node 117")
oct_GO <- GO_setup(117) # octocorallia

print("processing node 113")
cer_GO <- GO_setup(113) # ceriantheria

print("processing node 88")
hex_GO <- GO_setup(88) #hexacorallia

print("processing node 91")
scler_GO <- GO_setup(91) #scleractinia

print("processing node 136")
scyph_GO <- GO_setup(136) #scyphozoa

print("processing node 133")
cubscyph_GO <- GO_setup(133) # cubozoa/scyphozoa ancestor


save(anth_GO, anth2_GO, endomed_GO, med_GO, endo_GO, file = "early_nodes.Rdata")
save(myx_GO, hyd_GO, scs_GO, hex_GO, oct_GO, cer_GO, scler_GO, scyph_GO, cubscyph_GO, file = "recent_nodes.Rdata")


