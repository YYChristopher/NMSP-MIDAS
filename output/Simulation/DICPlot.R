library(data.table)
library(ggplot2)
library(latex2exp)

work_dir = 'F:/24561/Documents/ResearchWrite/NHMM-MIDAS/Submission Requirements/Code in Github'
#Set the Work Space of the data
setwd(work_dir)

# DIC
Rep = 100
DIC_path <- './output/Simulation/DIC_Result'
# NS = 1
file_name <- list.files(DIC_path, pattern = 'DIC1.txt', full.names = TRUE)
dic1.res <- readLines(file_name)
all_dic1 <- c(as.numeric(strsplit(dic1.res[1], split = ",")[[1]]))

# NS = 2
file_name <- list.files(DIC_path, pattern = 'DIC2.txt', full.names = TRUE)
dic2.res <- readLines(file_name)
all_dic2 <- c(as.numeric(strsplit(dic2.res[1], split = ",")[[1]]))

# NS = 3
file_name <- list.files(DIC_path, pattern = 'DIC3.txt', full.names = TRUE)
dic3.res <- readLines(file_name)
all_dic3 <- c(as.numeric(strsplit(dic3.res[1], split = ",")[[1]]))

# NS = 4
file_name <- list.files(DIC_path, pattern = 'DIC4.txt', full.names = TRUE)
dic4.res <- readLines(file_name)
all_dic4 <- c(as.numeric(strsplit(dic4.res[1], split = ",")[[1]]))

dic.df <- data.frame(
  dic = c(all_dic1, all_dic2, all_dic3, all_dic4),
  states = c(rep('1', Rep), rep('2', Rep), rep('3', Rep), rep('4', Rep))
)

DIC_plot <- dic.df %>%  
  ggplot(aes(x = states, y = dic, fill = states)) + 
  stat_boxplot(geom = "errorbar", width = 0.35) +  
  geom_boxplot() + 
  # theme(
  #   panel.background = element_rect(fill = "white", color = NA),  
  #   panel.grid.major = element_line(color = "grey80"),            
  #   panel.grid.minor = element_line(color = "grey90"),            
  #   axis.line = element_line(color = "black")                     
  # ) +
  # theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_boxplot(fill = "grey", color = "black") +
  xlab(expression(paste("Number of states ", italic(S)))) +
  ylab("DIC")
ggsave("./output/Simulation/Plots/DIC.png", plot = DIC_plot, width = 14, height = 8, units = "in", dpi = 600)