library(ggplot2)
print("ggplot loaded")
library(dplyr)
print("dplyr loaded")

ontology_colors <- c("BP" = "palevioletred", "CC" = "steelblue3", "MF" = "seagreen3", "all" = "orange2")

dir = "./deepfri_output"

bp <- read.csv(file.path(dir, "BP.csv"), header = FALSE)
print("bp loaded")
cc <- read.csv(file.path(dir, "CC.csv"), header = FALSE)
print("cc loaded")
mf <- read.csv(file.path(dir, "MF.csv"), header = FALSE)
print("mf loaded")

colnames(bp) <- c("Protein", "GO", "Score", "Description")
colnames(cc) <- c("Protein", "GO", "Score", "Description")
colnames(mf) <- c("Protein", "GO", "Score", "Description")

bp$Ontology = "BP"
cc$Ontology = "CC"
mf$Ontology = "MF"

data <- rbind(bp, cc, mf)
#colnames(data) <- c("Protein", "GO", "Score", "Description", "Ontology")
print("rbind done")

pdf("./score_histogram.pdf", height = 5, width = 8)
ggplot(data = data, aes(x = Score, fill = Ontology)) +
  geom_histogram(bins = 30, col = "black") +
  scale_fill_manual(values = ontology_colors) +
  xlab("DeepFRI score") +
  ylab("Number of annotations") +
  theme_bw()
dev.off()
print("histogram done")


scores <- data.frame(filter = character(), count = numeric(), ontology = character())

for (thisfilter in c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) {
  # total of all the ontologies
  uniq_count <- data %>% filter(Score >= thisfilter) %>% distinct(Protein) %>% nrow()
  scores <- rbind(scores, data.frame(filter = thisfilter, ontology = "all", count = uniq_count))
  # biological process
  uniq_count <- bp %>% filter(Score >= thisfilter) %>% distinct(Protein) %>% nrow()
  scores <- rbind(scores, data.frame(filter = thisfilter, ontology = "BP", count = uniq_count))
  # cellular component
  uniq_count <- cc %>% filter(Score >= thisfilter) %>% distinct(Protein) %>% nrow()
  scores <- rbind(scores, data.frame(filter = thisfilter, ontology = "CC", count = uniq_count))
  # molecular function
  uniq_count <- mf %>% filter(Score >= thisfilter) %>% distinct(Protein) %>% nrow()
  scores <- rbind(scores, data.frame(filter = thisfilter, ontology = "MF", count = uniq_count))
  write.csv(scores, file="scores.csv")
}

pdf("barplot.pdf", height = 5, width = 8)
ggplot(data = scores, aes(x = as.character(filter), y = count/1000000, fill = ontology)) +
  geom_bar(stat = "identity", position = position_dodge(), col = "black") +
  geom_hline(yintercept = 2684787/1000000, linetype = "dashed") +
  scale_fill_manual(values = ontology_colors) +
  xlab("DeepFRI Score Filter") +
  ylab("Number of centroids annotated (millions)") +
  annotate("text", x = 7.5, y = 2.8, label = "Total number of centroids", size = 5) +
  theme_bw()
dev.off()
print("script ended")


