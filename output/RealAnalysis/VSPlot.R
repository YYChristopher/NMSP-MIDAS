library("midasr")
library(gtools)
library(data.table)
library(ggplot2)
library(latex2exp)
library(tools)

work_dir = 'F:/24561/Documents/ResearchWrite/NHMM-MIDAS/Submission Requirements/Code in Github'
#Set the Work Space of the data
setwd(work_dir)

H <- 7
NT_test <- 51
NS <- 2
Nbvar <- 122

group1 = c(1, 2) # Income
group2 = c(seq(6, 19)) # Output
group3 = c(seq(20, 47), seq(111, 113)) # Labour
group4 = c(seq(48, 57)) # Housing
group5 = c(seq(3, 5), 114) # Consumption
group6 = c(seq(58, 61)) # Order
group7 = c(62, 63) # Inventories
group8 = c(seq(64, 73), seq(115, 117)) # Money & Credit
group9 = c(seq(78, 85)) # Interest Rates
group10 = c(seq(86, 90)) # Exchange Rates
group11 = c(seq(91, 96)) # PPI
group12 = c(seq(97, 106)) # CPI
group13 = c(seq(107, 110)) # PCE
group14 = c(seq(74, 77)) # Stock Market
group15 = c(seq(118, 121)) # FF3
group16 = c(122) # ADS

Groups = vector("list", length = numGroup)
for (n in 1:numGroup) {
  var_name <- paste0("group", n)
  Groups[[n]] <- get(var_name)
}

pi1Heat.Record = vector("list", length = H + 1)
for (h in 0:H) {
  pi1.Heatmap = array(0, dim = c(NS, NT_test - 1, numGroup))
  for (nt in 1:(NT_test - 1)) {
    
    directory_path <- sprintf("./output/RealAnalysis/Parameter_Estimation_Stream%d", h)
    pattern <- sprintf("Paras_est_Number%d.txt", nt)
    file_name <- list.files(directory_path, pattern = pattern, full.names = TRUE)
    
    est.res_list <- readLines(file_name)
    beta.estimate <- matrix(as.numeric(strsplit(est.res_list[1], split = ",")[[1]]), nrow = NS, ncol = Nbvar, byrow = FALSE)
    # beta.estimate <- beta.estimate[ ,-Nbvar]
    intercept.estimate <- c(as.numeric(strsplit(est.res_list[2], split = ",")[[1]]))
    sig2.estimate <- c(as.numeric(strsplit(est.res_list[3], split = ",")[[1]]))
    phi.estimate <- c(as.numeric(strsplit(est.res_list[4], split = ",")[[1]]))
    zeta.estimate <- matrix(as.numeric(strsplit(est.res_list[5], split = ",")[[1]]), nrow = NS, ncol = NS, byrow = FALSE)
    pi1.estimate <- matrix(as.numeric(strsplit(est.res_list[6], split = ",")[[1]]), nrow = NS, ncol = Nbvar, byrow = FALSE)
    
    for (s in 1:NS) {
      pi1.ns.est <- 1 - pi1.estimate[s, ]
      for (n in 1:numGroup) {
        # pi1Sum <- sum(pi1.ns.est[Groups[[n]]])
        pi1Max <- max(pi1.ns.est[Groups[[n]]])
        pi1.Heatmap[s, nt, n] <- pi1Max
      }
    }
  }
  pi1Heat.Record[[h + 1]] <- pi1.Heatmap
}

################## Plot ##########################
# State 1
data.toplot = as.matrix(pi1Heat.Record[[1]][1, , ])
colnames(data.toplot) <- c("Income", "Output", "Labour", "Housing", "Consumption", "Orders",
                           "Inventories", "Money&Credit", "Interest_Rates", "Exchange_Rates",
                           "PPI_Prices", "CPI_Prices", "PCE_Prices", "Stock_Market", "FF3", "ADS")
rownames(data.toplot) <- c(1:50)
data_long <- melt(data.toplot, varnames = c("Row", "Column"), value.name = "Value")

heatmap.h0s1 <- ggplot(data_long, aes(x = Row, y = Column, fill = Value)) +
  geom_tile(color = "black", width = 1, height = 1, position = "identity") +
  scale_fill_gradientn(
    colors = c("#FFEBEE", "#FFCDD2", "#EF9A9A", "#E57373", "#EF5350", "#B71C1C"),  
    values = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1),  
    limits = c(0, 1) 
  ) +
  # labs(title = "Out-of-sample h=0, S = 1", x = "Out-of-sample Quarter", y = "", fill = "Prob") +
  labs(x = "Out-of-sample Quarter", y = "", fill = "Prob") +
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_discrete(expand = c(0, 0)) +    
  scale_y_discrete(labels = c(
    "Output" = "Output",
    "Income" = "Income",
    "Labour" = "Labour",
    "Housing" = "Housing",
    "Consumption" = "Consumption",
    "Orders" = "Orders",
    "Inventories" = "Inventories",
    "Money&Credit" = "Money",
    "Interest_Rates" = "Interest",
    "Exchange_Rates" = "Exchange",
    "PPI_Prices" = "PPI",
    "CPI_Prices" = "CPI",
    "PCE_Prices" = "PCE",
    "Stock_Market" = "Stock", 
    "FF3" = "FF3",
    "ADS" = "ADS"
  )) +
  theme_minimal() +
  theme( 
    # plot.background = element_rect(fill = "white"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 16),
    text = element_text(size = 14)
    # legend.position = "none"
  ) +
  coord_fixed()
show(heatmap.h0s1)
# ggsave("./output/RealAnalysis/Plots/h0s1.png", plot = heatmap.h0s1, width = 8, height = 7, units = "in", dpi = 600)

# State 2
data.toplot = as.matrix(pi1Heat.Record[[1]][2, , ])
colnames(data.toplot) <- c("Output", "Income", "Labour", "Housing", "Consumption", "Orders",
                           "Inventories", "Money&Credit", "Interest_Rates", "Exchange_Rates",
                           "PPI_Prices", "CPI_Prices", "PCE_Prices", "Stock_Market", "FF3", "ADS")
rownames(data.toplot) <- c(1:50)
data_long <- melt(data.toplot, varnames = c("Row", "Column"), value.name = "Value")

heatmap.h0s2 <- ggplot(data_long, aes(x = Row, y = Column, fill = Value)) +
  geom_tile(color = "black", width = 1, height = 1, position = "identity") +
  scale_fill_gradientn(
    colors = c("#FFEBEE", "#FFCDD2", "#EF9A9A", "#E57373", "#EF5350", "#B71C1C"),
    values = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1),  # 映射关键值点
    limits = c(0, 1)  
  ) +
  # labs(title = "Out-of-sample h=0, S = 1", x = "Out-of-sample Quarter", y = "", fill = "Prob") +
  labs(x = "Out-of-sample h=0, S = 1", y = "", fill = "Prob") +
  scale_x_continuous(expand = c(0, 0)) +   
  scale_y_discrete(expand = c(0, 0)) +     
  scale_y_discrete(labels = c(
    "Output" = "Output",
    "Income" = "Income",
    "Labour" = "Labour",
    "Housing" = "Housing",
    "Consumption" = "Consumption",
    "Orders" = "Orders",
    "Inventories" = "Inventories",
    "Money&Credit" = "Money",
    "Interest_Rates" = "Interest",
    "Exchange_Rates" = "Exchange",
    "PPI_Prices" = "PPI",
    "CPI_Prices" = "CPI",
    "PCE_Prices" = "PCE",
    "Stock_Market" = "Stock", 
    "FF3" = "FF3",
    "ADS" = "ADS"
  )) +
  theme_minimal() +
  theme( 
    # plot.background = element_rect(fill = "white"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 16),
    text = element_text(size = 14)
    # legend.position = "none"
  ) +
  coord_fixed()
show(heatmap.h0s2)
# ggsave("./output/RealAnalysis/Plots/h0s2.png", plot = heatmap.h0s2, width = 8, height = 7, units = "in", dpi = 600)

h_values = c("0", "1/3", "2/3", "1", "4/3", "5/3", "2", "4")
for (h in 0:7) {
  h.Used = h_values[h + 1]
  for (s in 1:NS) {
    data.toplot = as.matrix(pi1Heat.Record[[h + 1]][s, , ])
    colnames(data.toplot) <- c("Output", "Income", "Labour", "Housing", "Consumption", "Orders",
                               "Inventories", "Money&Credit", "Interest_Rates", "Exchange_Rates",
                               "PPI_Prices", "CPI_Prices", "PCE_Prices", "Stock_Market", "FF3", "ADS")
    rownames(data.toplot) <- c(1:50)
    data_long <- melt(data.toplot, varnames = c("Row", "Column"), value.name = "Value")
    
    heatmap.toPlot <- ggplot(data_long, aes(x = Row, y = Column, fill = Value)) +
      geom_tile(color = "black", width = 1, height = 1, position = "identity") +
      scale_fill_gradientn(
        colors = c("#FFEBEE", "#FFCDD2", "#EF9A9A", "#E57373", "#EF5350", "#B71C1C"),  # 渐变颜色
        values = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1), 
        limits = c(0, 1)  
      ) +
      # labs(title = "Out-of-sample h=0, S = 1", x = "Out-of-sample Quarter", y = "", fill = "Prob") +
      labs(x = paste("h =", h.Used, ", s =", s), y = "", fill = "Prob") +
      scale_x_continuous(expand = c(0, 0)) +  
      scale_y_discrete(expand = c(0, 0)) +     
      scale_y_discrete(labels = c(
        "Output" = "Output",
        "Income" = "Income",
        "Labour" = "Labour",
        "Housing" = "Housing",
        "Consumption" = "Consumption",
        "Orders" = "Orders",
        "Inventories" = "Inventories",
        "Money&Credit" = "Money",
        "Interest_Rates" = "Interest",
        "Exchange_Rates" = "Exchange",
        "PPI_Prices" = "PPI",
        "CPI_Prices" = "CPI",
        "PCE_Prices" = "PCE",
        "Stock_Market" = "Stock", 
        "FF3" = "FF3",
        "ADS" = "ADS"
      )) +
      theme_minimal() +
      theme( 
        # plot.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 16),
        text = element_text(size = 14)
        # legend.text = element_text(size = 14),    
        # legend.title = element_text(size = 16)    
        # legend.position = "none"
      ) +
      coord_fixed()
    # show(heatmap.toPlot)
    filename <- file.path("./output/RealAnalysis/Plots", paste0("h", h, "s", s, ".png"))
    ggsave(filename = filename, plot = heatmap.toPlot, width = 12, height = 8, units = "in", dpi = 600)
  }
}

