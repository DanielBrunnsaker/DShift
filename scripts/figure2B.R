
library(ggplot2)
library(rstudioapi)

cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

rme_volcano = read.csv('../results/phase/phase.csv')[c('X','P.Value','logFC')]
rme_volcano$change = 'ns'
for (row in 1:nrow(rme_volcano)){
  if ((rme_volcano$logFC[row] <= -0.58) & (rme_volcano$P.Value[row] < 0.1)){
    rme_volcano$change[row] = 'p-value and log2FC'
  }
  
  if (rme_volcano$logFC[row] >= 0.58 & rme_volcano$P.Value[row] < 0.1){
    rme_volcano$change[row] = 'p-value and log2FC'
  }
  
  if (abs(rme_volcano$logFC[row]) <= 0.58 & rme_volcano$P.Value[row] < 0.1){
    rme_volcano$change[row] = 'p-value'
  }
}

cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey")
cols <- c("p-value and log2FC" = "#ffad73", "p-value" = "#26b3ff", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
sizes <- c("p-value and log2FC" = 4, "p-value" = 3, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
alphas <- c("p-value and log2FC" = 1, "p-value" = 0.75, "ns" = 0.5)

ggplot(rme_volcano,aes(x = logFC,
                       y = -log10(P.Value),
                       fill = change,    
                       size = change,
                       alpha = change)) + 
  geom_point(shape = 21, colour = "black") + 
  geom_hline(yintercept = -log10(0.1),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-6, 6, 2)),       
                     limits = c(-6, 6))+theme_bw()+ theme(legend.position="bottom",axis.line = element_line(colour = "black"),text = element_text(size=20))
