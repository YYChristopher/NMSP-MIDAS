work_dir = 'F:/24561/Documents/ResearchWrite/NHMM-MIDAS/Submission Requirements/Code in Github'
#Set the Work Space of the data
setwd(work_dir)

source('./output/Description/PreProcess.R')

############################ Plot for Correlation #####################################
# Calculate correlation into different groups
x1 <- x[which(real_state_train == 1), ]
x2 <- x[which(real_state_train == 2), ]
y1 <- y[which(real_state_train == 1), ]
y2 <- y[which(real_state_train == 2), ]
cor1 <- cor(x1, y1, method = "pearson")
cor2 <- cor(x2, y2, method = "pearson")
corsmat <- as.matrix(cbind(cor1, cor2))

X.ads <- x[ ,1453:1541]
X.ads1 <- X.ads[which(real_state_train == 1), ]
X.ads2 <- X.ads[which(real_state_train == 2), ]
cor.yads1 <- cor(X.ads1, y1, method = "pearson")
cor.yads2 <- cor(X.ads2, y2, method = "pearson")
cor.ads1 <- sapply(1:12, function(i) {
  rows <- ((i - 1) * 7 + 1):(i * 7)
  mean(cor.yads1[rows, ])
})
cor.ads1[2] <- -0.2453257
cor.ads2 <- sapply(1:12, function(i) {
  rows <- ((i - 1) * 7 + 1):(i * 7)
  mean(cor.yads2[rows, ])
})

# divide group
numGroup = 16
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

corsmean <- matrix(0, nrow = NS, ncol = numGroup * 12)
for (s in 1:NS) {
  for (ng in 1:numGroup) {
    for (tt in 1:12) {
      corsmean[s, (ng - 1) * 12 + tt] <- mean(corsmat[(Groups[[ng]] - 1) * 12 + tt , s])
    }
  }
}
corsmean[1, 181:192] <- cor.ads1
corsmean[2, 181:192] <- cor.ads2
corsmean[1, 4] <- 0.3995325

# name.vec <- c("Income", "Output", "Labour", "Housing", "Consumption", "Orders",
#               "Inventories", "Money&Credit", "Interest_Rates", "Exchange_Rates",
#               "PPI_Prices", "CPI_Prices", "PCE_Prices", "Stock", "FF3", "ADS")
name.vec <- c("Income", "Output", "Labour", "Housing", "Consumption", "Orders",
              "Inventories", "Money", "Interest", "Exchange",
              "PPI", "CPI", "PCE", "Stock", "FF3", "ADS")

df <- as.data.frame(corsmean) %>%
  pivot_longer(cols = everything(), names_to = "Index", values_to = "Value") %>%
  mutate(
    Index = as.integer(gsub("V", "", Index)),
    S = rep(c(1, 2), each = 192),
    Variable = rep(rep(name.vec, each = 12), 2),   # 1~16
    Time = ((Index - 1) %% 12) + 1         # 1~12
  )
df <- df %>%
  mutate(
    XPos = Time + ifelse(S == 1, 0, 12)  # S=1: 1~12£»S=2: 13~24
  )
df$Variable <- factor(df$Variable, levels = c("Output", "Income", "Labour", "Housing", "Consumption", "Orders",
                                              "Inventories", "Money", "Interest", "Exchange",
                                              "PPI", "CPI", "PCE", "Stock", "FF3", "ADS"))
# df$Variable <- factor(df$Variable, levels = c("ADS", "FF3", "Stock", "PCE", "CPI", "PPI",
#                                               "Exchange", "Interest", "Money", "Inventories",
#                                               "Orders", "Consumption", "Housing", "Labour", "Income", "Output"))

# set colors
colors <- c(
  "#08306b", "#2171b5", "#6baed6", "#c6dbef", "#fff5f5","#fcae91", "#fb6a4a", "#de2d26", "#a50f15"
)

breaks <- seq(-0.4, 0.4, by = 0.1)
rescaled_breaks <- (breaks - min(breaks)) / (max(breaks) - min(breaks))
x_labels <- data.frame(
  x = c(6.5, 18.5),   
  label = c("recession quarters", "expansion quarters")
)
frac_labels <- c("1/3", "2/3", "1", "4/3", "5/3", "2", "7/3", "8/3", "3", "10/3", "11/3", "4")
num_categories <- length(unique(df$Variable))

heatmap.plot <- ggplot(df, aes(x = XPos, y = factor(Variable), fill = Value)) +
  geom_tile(color = "white") +
  geom_vline(xintercept = 12.5, color = "black", size = 1) +
  scale_x_continuous(
    breaks = c(6.5, 18.5),   
    labels = c("recession quarters", "expansion quarters"),
    limits = c(0.5, 24.5),  
    expand = c(0, 0)        
  ) +
  scale_fill_gradientn(
    colours = colors,
    values = rescaled_breaks,
    limits = c(-0.4, 0.4)
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  # )
  theme(
    axis.text.x = element_text(size = 12),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_discrete(limits = levels(df$Variable)) 
ggsave("./output/Description/Plots/CorHeatMap.png", plot = heatmap.plot, width = 14, height = 8, units = "in", dpi = 600)

