### plotting genexp vs pseudotime
plotgene_exp <- function(gene) { #gene is a string eg "Myh6"

#select pseudotime vecetor
pst <- as.matrix(HSMM_wtsub2@phenoData@data$Pseudotime)
#convert to normalized expression to log2 
expv <- log2(as.matrix(HSMM_wtsub2@assayData$exprs_sfnorm[gene,])+1)

#merge pseudotime and expressoin
p_exp <- as.data.frame(cbind(as.numeric(pst), as.numeric(expv)))
colnames(p_exp) <- c("pst","expv") #set column names

#plot using ggplot2
library(ggplot2)
ggplot(data = p_exp, aes(x= pst, y = expv)) +
  geom_point() +
  geom_smooth() +
  ggtitle(gene) +
  ylab("Log2 Expression") +
  xlab("Pseudotime") +
  title(gene) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=4),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")
  )
}