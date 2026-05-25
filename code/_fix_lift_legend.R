suppressPackageStartupMessages({library(tidyverse); library(ggsci)})
npg10 <- pal_npg("nrc")(10)
rules <- read.csv("D:/Users/Desktop/R_Work/results/association/association_rules.csv")
top_rules <- rules %>% arrange(desc(lift)) %>% head(20)

p <- ggplot(top_rules, aes(x=support, y=confidence, size=lift, color=lift)) +
  geom_point(alpha=0.8) +
  scale_color_gradient(low=npg10[5], high=npg10[3]) +
  scale_size_continuous(range=c(3,12)) +
  labs(title="Association Rules: Gene Expression to Clinical Features",
       subtitle=sprintf("Top 20 of %d rules (lift-ranked)", nrow(rules)),
       x="Support", y="Confidence") +
  theme_classic(base_size=15) +
  theme(plot.title=element_text(face="bold",hjust=0.5)) +
  guides(size = "none")
ggsave("D:/Users/Desktop/R_Work/results/figures_pub/fig15_association_rules.png",
       p, width=10, height=7, dpi=300, bg="white")
cat("Done\n")
