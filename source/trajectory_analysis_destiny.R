library(destiny)
library(Matrix)

destinySCEset <- filteredSCEset[ordering_genes, ]
dm <- DiffusionMap(t(as.matrix(logcounts(destinySCEset))))

# dpt <- DPT(dm)
# 
# plot(dpt, pch = 20)
# plot(dpt, col_by = "branch")

tmp <- data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  clusters = Seuratset@active.ident)

class(tmp)
dim(tmp)
head(tmp)

ggplot(tmp, aes(x = DC1, y = DC2, colour = clusters)) +
  geom_point() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()

destinySCEset$pseudotime_diffusionmap <- rank(-eigenvectors(dm)[,1])
destinySCEset$clusters <- Seuratset@active.ident

library(ggbeeswarm)
ggplot(as.data.frame(colData(destinySCEset)), 
       aes(x = pseudotime_diffusionmap, 
           y = clusters, colour = clusters)) +
  geom_quasirandom(groupOnX = FALSE) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(size=1), 
        axis.ticks=element_line(size=1), 
        # legend.title=element_blank(), 
        # legend.key=element_blank(), 
  ) +
  # scale_color_tableau() + 
  # theme_classic() +
  xlab("Diffusion map pseudotime (first diffusion map component)") +
  ylab("clusters") +
  ggtitle("Cells ordered by diffusion map pseudotime")




