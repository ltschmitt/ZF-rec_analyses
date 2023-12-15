library(tidyverse)
library(ggrepel)

d <- read_tsv("reads_vs_vika_aa.exonerate.txt.gz", col_names=F)

gap5 <- d %>% group_by(X1) %>% filter(sum(X2=="frameshift") < 7) %>% ungroup() %>% filter(X2=="gap") %>% group_by(X1, X3) %>% summarize(n=n(), .groups="drop") %>% filter(n==5) %>% group_by(X3) %>% summarize(n=n())

gap5_max <- gap5 %>% filter(n > 100) %>% mutate(Group=dist(X3) %>% hclust(method="complete") %>% cutree(h=30)) %>% group_by(Group) %>% filter(n>= 0.5*max(n)) %>% ungroup()

p <- ggplot(gap5, aes(x=X3, y=n))+geom_spoke(aes(radius=n, angle=-pi/2))+ggrepel::geom_text_repel(aes(label=X3), data=gap5_max, direction="x", vjust=0, color="darkgray")+scale_x_continuous("Residue position", limits=c(1,362), breaks=c(1, seq(30, 362, by=30)), expand=c(0, 30))+scale_y_continuous("Count of sequencing reads", expand=expansion(mult=c(0, 0.05)))+theme_classic(12)+labs(title="Distribution of pentapeptide insetions in front of given residue positions")

ggsave("Distribution of pentapeptide insetions in front of given residue positions.pdf", p, width=8, height=4.5, device=cairo_pdf)
