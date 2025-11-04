
library(dplyr)

load("deepfri.Rdata")
print("deepfri loaded")
bp <- bp %>% rename(GO = GO_term.EC_number, name = GO_term.EC_number.name) %>%
	filter(Score > 0.5)
cc <- cc %>% rename(GO = GO_term.EC_number, name = GO_term.EC_number.name) %>%
        filter(Score > 0.5)
mf <- mf %>% rename(GO = GO_term.EC_number, name = GO_term.EC_number.name) %>%
        filter(Score > 0.5)

save(bp, cc, mf, file="deepfri50.Rdata") 

focal_nodes <- data.frame(number = c(132, 133, 103, 104, 129, 130, 105, 116, 163, 159, 134, 137, 120, 108, 102), name = c("anth", "anth2", "endomed", "med", "endo", "myx", "hyd", "acr", "oct", "cer", "hex", "scler", "scyph", "siph", "cnid"))

load("../../orthodollo/tree.Rds")
tree <- tr1
rm(tr1)
print("tree loaded")

load("../../orthodollo/loss_gain_pres_abs_H1.Rds")
gl <- lossgain.table
rm(lossgain.table)
print("lossgain loaded")

gl <- as.data.frame(t(gl))
colnames(gl) <- as.character(1:ncol(gl))
gl$og <- rownames(gl)
print("gl prepped")

bp_OG_GO <- bp %>% mutate(og = sub("_.*", "", Protein)) %>% 
  dplyr::select(og, GO) %>% 
  group_by(og) %>% 
  summarize(GOs = paste(GO, collapse = ","), .groups = "drop")
print("bp collapsed")

cc_OG_GO <- cc %>% mutate(og = sub("_.*", "", Protein)) %>% 
  dplyr::select(og, GO) %>% 
  group_by(og) %>% 
  summarize(GOs = paste(GO, collapse = ","), .groups = "drop")
print("cc collapsed")

mf_OG_GO <- mf %>% mutate(og = sub("_.*", "", Protein)) %>% 
  dplyr::select(og, GO) %>% 
  group_by(og) %>% 
  summarize(GOs = paste(GO, collapse = ","), .groups = "drop")
print("mf collapsed")

bp_GOs <- data.frame(node = character(), g_GOs = vector(), l_GOs = vector())
cc_GOs <- data.frame(node = character(), g_GOs = vector(), l_GOs = vector())
mf_GOs <- data.frame(node = character(), g_GOs = vector(), l_GOs = vector())

for (NODE in 1:193) {
  
  if (NODE %in% focal_nodes$number) {
    print(paste("starting node: ", NODE))

    gained <- gl[gl[[NODE]] == "g", NODE, drop = FALSE] 
    gained$og <- rownames(gained) 
    lost <- gl[gl[[NODE]] == "l", NODE, drop = FALSE] 
    lost$og <- rownames(lost) 
  
    #bp
    gained_GOs <- merge(bp_OG_GO, gained) %>% select(-og) %>% group_by(!!sym(as.character(NODE))) %>% 
      summarize(GOs = paste(GOs, collapse = ","), .groups = "drop")
    lost_GOs <- merge(bp_OG_GO, lost) %>% select(-og) %>% group_by(!!sym(as.character(NODE))) %>% 
      summarize(GOs = paste(GOs, collapse = ","), .groups = "drop")
  
    bp_GOs <- rbind(bp_GOs, data.frame(node = NODE, g_GOs = I(list(strsplit(gained_GOs$GOs[1], ",")[[1]])), l_GOs = I(list(strsplit(lost_GOs$GOs[1], ",")[[1]]))))
  
    #cc
    gained_GOs <- merge(cc_OG_GO, gained) %>% select(-og) %>% group_by(!!sym(as.character(NODE))) %>% 
      summarize(GOs = paste(GOs, collapse = ","), .groups = "drop")
    lost_GOs <- merge(cc_OG_GO, lost) %>% select(-og) %>% group_by(!!sym(as.character(NODE))) %>% 
      summarize(GOs = paste(GOs, collapse = ","), .groups = "drop")
    cc_GOs <- rbind(cc_GOs, data.frame(node = NODE, g_GOs = I(list(strsplit(gained_GOs$GOs[1], ",")[[1]])), l_GOs = I(list(strsplit(lost_GOs$GOs[1], ",")[[1]]))))
  
    #mf
    gained_GOs <- merge(mf_OG_GO, gained) %>% select(-og) %>% group_by(!!sym(as.character(NODE))) %>% 
      summarize(GOs = paste(GOs, collapse = ","), .groups = "drop")
    lost_GOs <- merge(mf_OG_GO, lost) %>% select(-og) %>% group_by(!!sym(as.character(NODE))) %>% 
      summarize(GOs = paste(GOs, collapse = ","), .groups = "drop")
    mf_GOs <- rbind(mf_GOs, data.frame(node = NODE, g_GOs = I(list(strsplit(gained_GOs$GOs[1], ",")[[1]])), l_GOs = I(list(strsplit(lost_GOs$GOs[1], ",")[[1]]))))
  }
  else {
  print(paste("Skipping node:", NODE))
  bp_GOs <- rbind(bp_GOs, data.frame(node = NODE, g_GOs = NA, l_GOs = NA))
  cc_GOs <- rbind(cc_GOs, data.frame(node = NODE, g_GOs = NA, l_GOs = NA))
  mf_GOs <- rbind(mf_GOs, data.frame(node = NODE, g_GOs = NA, l_GOs = NA))
  }
}

print("setup done")

safe_fisher <- function(a, b, c, d) {
  mat <- matrix(c(a, b - a, c, d - c), nrow = 2, byrow = TRUE)
  if (any(mat < 0) || any(is.na(mat)) || any(!is.finite(mat))) return(1)
  fisher.test(mat, alternative = "greater")$p.value
}

GO_setup <- function(DATA, NODE) {
  
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
    
    GO_counts <- GO_counts %>% 
  	filter(!is.na(GO_term) & grepl("^GO:", GO_term))
    
    GO_counts$pval_gain <- apply(GO_counts, 1, function(row) {
      a <- as.numeric(row["gained"])
      go <- row[["GO_term"]]
      c <- as.numeric(all_GO_gain_counts[go])
      #mat <- matrix(c(a,b-a,c,d-c),nrow = 2, byrow = TRUE)
      #fisher.test(mat, alternative = "greater")$p.value
      safe_fisher(a,b,c,d)
    })
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
    
    GO_counts$pval_loss <- apply(GO_counts, 1, function(row) {
      a <- as.numeric(row["lost"])
      go <- row[["GO_term"]]
      c <- as.numeric(all_GO_loss_counts[go])
      #mat <- matrix(c(a,b-a,c,d-c),nrow = 2, byrow = TRUE)
      #fisher.test(mat, alternative = "greater")$p.value
      safe_fisher(a,b,c,d)
    })
    GO_counts$padj_loss <- p.adjust(GO_counts$pval_loss, method = "bonferroni")
  }
  
  
  return(GO_counts)
}
print("GO setup initialized")

process <- function(NODE, DIR) {
  dir.create(DIR)
  GO_df <- GO_setup(bp_GOs, NODE)
  GO_df <- GO_df %>% mutate(padj_loss = if_else(padj_loss == 0, 1e-300, padj_loss),
			    padj_gain = if_else(padj_gain == 0, 1e-300, padj_gain))
  write.table(GO_df, file.path(DIR, "BP_GO_df.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  if(nrow(filter(GO_df, gained != 0)) == 0){
    print(paste("No BP gains in node: ", NODE))
  } else {
    gain <- GO_df %>% 
      dplyr::select(c("GO_term", "padj_gain", "gained"))
    write.table(gain, file.path(DIR, "BP_gain.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
	#write_tsv(gain, file.path(DIR, "BP_gain.tsv"), col_names = FALSE)
  }
  if(nrow(filter(GO_df, lost != 0)) == 0){
    print(paste("No BP losses in node: ", NODE))
  } else {
    loss <- GO_df %>% 
      dplyr::select(c("GO_term", "padj_loss", "lost"))  
    #write_tsv(loss, file.path(DIR, "BP_loss.tsv"), col_names = FALSE)
    write.table(loss, file.path(DIR, "BP_loss.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  GO_df <- GO_setup(cc_GOs, NODE)
  GO_df <- GO_df %>% mutate(padj_loss = if_else(padj_loss == 0, 1e-300, padj_loss),
                            padj_gain = if_else(padj_gain == 0, 1e-300, padj_gain))
  write.table(GO_df, file.path(DIR, "CC_GO_df.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  if(nrow(filter(GO_df, gained != 0)) == 0){
    print(paste("No CC gains in node: ", NODE))
  } else {
    gain <- GO_df %>% 
      dplyr::select(c("GO_term", "padj_gain", "gained"))  
      #filter(padj_gain <= 0.01) %>% 
      #mutate(padj_gain = if_else(padj_gain == 0, 1e-300, padj_gain))
    write.table(gain, file.path(DIR, "CC_gain.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    #write_tsv(gain, file.path(DIR, "CC_gain.tsv"), col_names = FALSE)
  }
  if(nrow(filter(GO_df, lost != 0)) == 0){
    print(paste("No CC losses in node: ", NODE))
  } else {
    loss <- GO_df %>% 
      dplyr::select(c("GO_term", "padj_loss", "lost"))  
      #filter(padj_loss <= 0.01) %>% 
      #mutate(padj_loss = if_else(padj_loss == 0, 1e-300, padj_loss))
    #write_tsv(loss, file.path(DIR, "CC_loss.tsv"), col_names = FALSE)
    write.table(loss, file.path(DIR, "CC_loss.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  GO_df <- GO_setup(mf_GOs, NODE)
  GO_df <- GO_df %>% mutate(padj_loss = if_else(padj_loss == 0, 1e-300, padj_loss),
                            padj_gain = if_else(padj_gain == 0, 1e-300, padj_gain))
  write.table(GO_df, file.path(DIR, "MF_GO_df.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  if(nrow(filter(GO_df, gained != 0)) == 0){
    print(paste("No MF gains in node: ", NODE))
  } else {
    gain <- GO_df %>% 
      dplyr::select(c("GO_term", "padj_gain", "gained"))  
      #filter(padj_gain <= 0.01) %>% 
      #mutate(padj_gain = if_else(padj_gain == 0, 1e-300, padj_gain))
    write.table(gain, file.path(DIR, "MF_gain.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    #write_tsv(gain, file.path(DIR, "MF_gain.tsv"), col_names = FALSE)
  }
  if(nrow(filter(GO_df, lost != 0)) == 0){
    print(paste("No MF losses in node: ", NODE))
  } else {
    loss <- GO_df %>% 
      dplyr::select(c("GO_term", "padj_loss", "lost")) 
      #filter(padj_loss <= 0.01) %>% 
      #mutate(padj_loss = if_else(padj_loss == 0, 1e-300, padj_loss))
    #write_tsv(loss, file.path(DIR, "MF_loss.tsv"), col_names = FALSE)
    write.table(loss, file.path(DIR, "MF_loss.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}
print("process initialized")

dir="go_figure_tsvs50"

dir.create(dir)

for (i in 1:nrow(focal_nodes)) {
  print(paste("processing node:", focal_nodes$number[i]))
  process(focal_nodes$number[i], file.path(dir,focal_nodes$name[i]))
}

print("gofigure setup done")




