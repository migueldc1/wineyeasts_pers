library(igraph)
library(hillR)
library(ggplot2)
library(reshape2)
library(Biostrings)
library(msa)
library(phangorn)
library(ggtree)

setwd("PATH/TO/YOUR/FILES")

#### FIGURE 1

### PANEL A

## Data loading
asv_df <- t(read.table("asv_df.txt", sep = "\t", header = TRUE, row.names = 1))
meta_df <- read.table("meta_df.txt", sep = "\t", header = TRUE)
tax_df <- read.table("tax_ws.txt", sep = "\t", header = TRUE)
source("veech_prob.R")

## Data filtering
# Remove Unidentified at the genus level taxons
asv_df <- asv_df[,grepl("Unidentified", colnames(asv_df)) == FALSE] 

# Keep genera with more than 0.5% abundance in at least a 5% of alcoholic samples
asv_df.alc <- asv_df[subset(meta_df, Stage == "Alcoholic")$Sample_ID,]
asv_df.alc <- asv_df.alc[, colSums(asv_df.alc > 0.005) > 0.05*nrow(asv_df.alc)] 

# Keep only Alcholic fermentation and Must samples
asv_df.net <- asv_df[subset(meta_df, Stage == "Alcoholic" | Stage == "Must")$Sample_ID, colnames(asv_df.alc)]
asv_df.net[asv_df.net > 0] <- 1

## Network analysis - Probabilistic Co-occurrence (Veech, 2013)

veech_df <- veech_prob(asv_df.net)

# Filter significant co-occurring pairs
edge_df <- subset(veech_df, p.val_Pgt <= 0.05)[,c(1,2,3)]

# Generate an igraph object
cooc.net <- graph_from_data_frame(as.matrix(edge_df)[,2:1], directed = FALSE)

# Create a node data frame with taxonomic and topological information
node_df <- data.frame(Genus = names(V(cooc.net)))
mod_df <- cbind.data.frame(Genus = cluster_walktrap(cooc.net)$names, module = cluster_walktrap(cooc.net)$membership)
node_df <- Reduce(function(...) merge(..., all = FALSE, by = "Genus"), list(node_df, tax_df, mod_df))
row.names(node_df) <- node_df$Genus

node_df <- node_df[names(V(cooc.net)),]

# Introduce the information in the igraph object
V(cooc.net)$Genus <- node_df$Genus
V(cooc.net)$Order <- node_df$Order
V(cooc.net)$community <- node_df$module

# Layout
set.seed(12)
l_ws <- layout_with_fr(cooc.net, grid = "nogrid")
l_ws <- norm_coords(layout_with_fr(cooc.net, grid = "nogrid"), ymin = -1, ymax = 1, xmin = -1, xmax = 1)

# Improve the plot
colrs <- adjustcolor(c("#00509c", "#b54702", "#2fa13a", "#b3b312", "#570e55"), alpha = 0.8)

# Select edges connecting Saccharomyces
idx <- match("Saccharomyces", V(cooc.net)$Genus)
E(cooc.net)$color[E(cooc.net)[!from(idx)]] <- "gray70"
E(cooc.net)$color[E(cooc.net)[from(idx)]] <- "#8a223c"

# Select nodes connected with Saccharomyces
nodes_sacc <- union(head_of(cooc.net, which(E(cooc.net)$color == "#8a223c")),
                    tail_of(cooc.net, which(E(cooc.net)$color == "#8a223c")))

# Plot the network
plot(cooc.net, 
     vertex.size = ifelse(V(cooc.net)$Order == "Saccharomycetales", 8, 6), 
     vertex.label = ifelse(V(cooc.net)$Genus == "Saccharomyces", V(cooc.net)$Genus, NA),
     vertex.color = colrs[V(cooc.net)$community], 
     vertex.frame.color = ifelse(V(cooc.net) %in% nodes_sacc, "#8a223c", "black"),
     rescale = FALSE, layout = l_ws*1.3, edge.curved = 0.25)

### PANEL B
# Filter ASV table
asv_net <- asv_df[,node_df$Genus]

# Calculate Richness in the whole community and in the community forming the network
alpha_net <- cbind.data.frame(Richness = hill_taxa(asv_df, q = 0),
                              Nw.Richness = hill_taxa(asv_net, q = 0))

alpha_net$Sample_ID <- row.names(alpha_net)
alpha_net <- merge(alpha_net, meta_df, by = "Sample_ID")

# Calculate mean values
alpha_net.plot <- cbind(aggregate(alpha_net[,2:3], list(alpha_net$Stage), mean))

colnames(alpha_net.plot) <- c("Stage", "Richness", "Nw.Richness")

alpha_net.plot$Stage <- factor(alpha_net.plot$Stage, levels = c("Must", "Alcoholic", "Post-Fermentative"))

# Plot
ggplot(data = alpha_net.plot) +
  geom_area(aes(x = as.numeric(Stage), y = Richness, fill = "Richness"), stat = "identity", position = "stack", 
            alpha = 0.8) +
  geom_area(aes(x = as.numeric(Stage), y = Nw.Richness, fill = "Nw.Richness"), stat = "identity", position = "stack", 
            alpha = 0.8) +
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("1" = "Must", "2" = "Alcoholic", "3" = "Post-Fermentative")) +
  scale_fill_manual(name = "", values = c(Richness = "#f2c14e", Nw.Richness = "#4ef2d4"), 
                    labels = c("Total", "Network")) +
  theme_bw() + 
  xlab("") + ylab("Average Richness") + 
  theme(axis.text.y = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))

### PANEL C
# Calculate the abundance of the ASVs corresponding to the different Network modules
asv_net.mod <- cbind.data.frame(Genus = colnames(asv_net), t(asv_net))
asv_net.mod <- merge(asv_net.mod, mod_df, by = "Genus")

asv_net.mod.plot <- aggregate(asv_net.mod[,2:273], list(asv_net.mod$module), sum)
asv_net.mod.plot <- melt(t(asv_net.mod.plot[,-1]))

colnames(asv_net.mod.plot) <- c("Sample_ID", "Module", "value")

asv_net.mod.plot <- merge(asv_net.mod.plot, meta_df, by = "Sample_ID")

asv_net.mod.plot <- cbind(aggregate(asv_net.mod.plot$value, list(asv_net.mod.plot$Module, 
                                                                 asv_net.mod.plot$Stage), mean))

colnames(asv_net.mod.plot) <- c("Module", "Stage", "value")


asv_net.mod.plot$Stage <- factor(asv_net.mod.plot$Stage, levels = c("Must", "Alcoholic", "Post-Fermentative"))

# Complete the abundance to 1
asv_net.mod.plot[nrow(asv_net.mod.plot) + 3 , ] <- NA
asv_net.mod.plot[16:18,] <- cbind(Stage = "Out", aggregate(asv_net.mod.plot$value, list(asv_net.mod.plot$Stage),
                                                           function(x) 1-sum(x)))

asv_net.mod.plot$Module <- factor(asv_net.mod.plot$Module, levels = c("Out", 4, 5, 3, 2, 1))

# Plot
ggplot(data = asv_net.mod.plot) +
  geom_area(aes(x = as.numeric(Stage), y = value, fill = Module), stat = "identity", position = "stack", 
            alpha = 0.8) +
  theme_bw() + 
  xlab("") + ylab("Module mean proportion") + 
  labs(fill = "Module") + ylim(0, 1) +
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("1" = "Must", "2" = "Alcoholic", "3" = "Post-Fermentative")) +
  scale_fill_manual(values = c("#00509c", "#b54702", "#2fa13a", "#570e55", "#b3b312", "gray70"),
                    breaks = c(1, 2, 3, 5, 4, "Out")) +
  theme(axis.text.y = element_text(size = 19, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 17, color = "black"),
        legend.title = element_text(size = 19, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"))


#### FIGURE 2

### PHYlOGENETIC TREE

## Align full ITS sequences
seqs_its <- readDNAStringSet("full_ITS.fasta")

mult_its <- msa(seqs_its, method = "ClustalW", type = "dna", order = "input")

## Create phylogenetic tree object
phang.align_its <- as.phyDat(mult_its, type = "DNA", names = seqs_its)
dm_its <- dist.ml(phang.align_its)
treeNJ_its <- NJ(dm_its)

## Fit GTR phylogenetic tree
fit_its <- pml(treeNJ_its, data = phang.align_its)

set.seed(123)
fitGTR_its <- update(fit_its, k = 4, inv = 0.2)
fitGTR_its <- optim.pml(fitGTR_its, model = "GTR", optInv = TRUE, optGamma = TRUE,
                         rearrangement = "stochastic", control = pml.control(trace = 0))

## Set aesthetics
tax_net <- tax_df[tax_df$Genus %in% asv_net.mod$Genus,]

ord.list_its <- vector("list", length(levels(as.factor(tax_net$Order))))
names(ord.list_its) <- levels(as.factor(tax_net$Order))
for(i in 1:length(ord.list_its)) {
  ord.list_its[[i]] <- subset(tax_net, Order == names(ord.list_its)[i])$Genus
}

tree_its <- groupOTU(fitGTR_its$tree, ord.list_its)

phy_tree.color <- cbind.data.frame(Order = c("0", unique(tax_net$Order)),
                                   Color = c("gray50", "#e3ad5d", "#8ab849", "#7678d6",
                                             "#75a8c7", "#186b05", "#80654a", "#577fb3",
                                             "#dbd053", "#b0773a", "#a60857", "#d01e3c",
                                             "#488896", "#65a6c9"))

# Plot the tree
ggtree(tree_its, aes(color = group), size = 1) +
  scale_color_manual(values = phy_tree.color$Color) +
  geom_tiplab(aes(), size = 4)

### OCCURRENCE / ABUNDANCE PLOTs
asv_net.data <- as.data.frame(asv_net)
asv_net.data$Sample_ID <- row.names(asv_net.data)

asv_net.data <- merge(asv_net.data, meta_df, by = "Sample_ID")

asv_net.data_df <- cbind.data.frame(melt(aggregate(asv_net.data[,2:33], list(asv_net.data$Stage), mean)),
                                    melt(aggregate(asv_net.data[,2:33], list(asv_net.data$Stage), 
                                                   function(x) sum(x > 0)/length(x)))[,3])

colnames(asv_net.data_df) <- c("Stage", "Genus", "Mean", "Prevalence")

max.net_df <- aggregate(asv_net.data_df$Mean, list(asv_net.data_df$Genus), function(x) max(x))
colnames(max.net_df) <- c("Genus", "Max")

asv_net.data_df <- merge(asv_net.data_df, max.net_df, by = "Genus")
asv_net.data_df$Mean.n <- asv_net.data_df$Mean/asv_net.data_df$Max

asv_net.data_df$Stage <- factor(asv_net.data_df$Stage, levels = c("Must", "Alcoholic", "Post-Fermentative"))

sacc.gen <- subset(tax_net, Order == "Saccharomycetales")$Genus

for (gen in sacc.gen) {
  print(ggplot(data = subset(asv_net.data_df, Genus == gen)) +
          geom_bar(aes(x = -as.numeric(Stage), y = Prevalence), stat = "identity", position = "stack",
                   fill = "#912437", alpha = 0.8, width = 0.2) +
          geom_area(aes(x = -as.numeric(Stage), y = Mean.n), stat = "identity", position = "stack",  
                    fill = "#c97785", alpha = 0.8) +
          coord_flip() +
          theme_void() +
          theme(plot.margin = margin(15, 10, 0, 10)) +
          labs(title = gen))
}
