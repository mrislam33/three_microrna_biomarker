#load required package
library(ggpubr)

#load required data
load("./data/Figure3Bdata.RData")

message("data contains cumulative cognitive score for each behaviourial performance day
        Cognitive score from mouse behaviourial experiment was calculated as described in detail in Methods and Material ")
        


names(data) = c("cognition", "age")

my_comparison = list(c("13.5m", "15m"), c("13.5m", "16.5m"), c("15m", "16.5m"))

ggboxplot(data, x = "age", 
          y = "cognition", 
          fill = "age",
          yticks.b = 0.2,
          short.panel.labs = FALSE, 
          outlier.shape = NA,
          ylab = "cognitive score", 
          title = "Figure 3B")+ stat_compare_means(comparisons = my_comparison)

