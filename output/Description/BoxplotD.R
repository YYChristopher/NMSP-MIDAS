work_dir = 'F:/24561/Documents/ResearchWrite/NHMM-MIDAS/Submission Requirements/Code in Github'
#Set the Work Space of the data
setwd(work_dir)

source('./output/Description/PreProcess.R')

m_mon = 3
NBER = read.csv('./data/USREC.csv', header = TRUE)
NBER_noindex <- NBER['USREC']
dates_nber_all <- seq(as.Date("1854-12-01"), by = "month", length.out = dim(NBER)[1])
NBER_xts <- xts(NBER_noindex, order.by = dates_nber_all)
dates_nber.quarter <- seq(as.Date("1960-01-01"), by = "quarter", length.out = NT_train + 1)
dates_nber.month <- seq(as.Date("1960-01-01"), by = "month", length.out = (NT_train + 1) * m_mon)
NBER_used <- NBER_xts['1960-01-01/2023-12-01']


D.std <- scale(D)
df.D <- data.frame(
  Date = c(rep(dates_nber.quarter, 4)),
  Value = c(D.std[ ,1], D.std[ ,2], D.std[ ,3], D.std[ ,4]),
  State = c(rep(c("SP.500", "PPI_Oil", "PPI_Aco", "ADS"), each = 256))
)

##### Two Boxplot

# d_{t}
df.state1 <- data.frame(
  SP500 = D.std[which(real_state_train == 1), 1],
  PPI_Oil = D.std[which(real_state_train == 1), 2],
  PPI_ACO = D.std[which(real_state_train == 1), 3],
  ADS = D.std[which(real_state_train == 1), 4]
)

df.state2 <- data.frame(
  SP500 = D.std[which(real_state_train == 2), 1],
  PPI_Oil = D.std[which(real_state_train == 2), 2],
  PPI_ACO = D.std[which(real_state_train == 2), 3],
  ADS = D.std[which(real_state_train == 2), 4]
)

df1_long <- df.state1 %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = "recession")

df2_long <- df.state2 %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = "expansion")

combined_data <- bind_rows(df1_long, df2_long)
combined_data$Group <- factor(combined_data$Group, levels = c("recession", "expansion"))
combined_data$Variable <- factor(combined_data$Variable, levels = c("SP500", "PPI_Oil", "PPI_ACO", "ADS"))

combined_data_bar <- combined_data %>%
  group_by(Variable, Group) %>%
  summarise(
    ymin = quantile(Value, 0.25) - 1.5 * IQR(Value),
    ymax = quantile(Value, 0.75) + 1.5 * IQR(Value),
    .groups = "drop"
  )

combined_data <- combined_data %>%
  group_by(Variable, Group) %>%
  mutate(
    ymin = quantile(Value, 0.25) - 1.5 * IQR(Value),
    ymax = quantile(Value, 0.75) + 1.5 * IQR(Value)
  ) %>%
  ungroup()

combined_data <- combined_data %>%
  group_by(Variable, Group) %>%
  mutate(
    ydown = quantile(Value, 0.25),
    yup = quantile(Value, 0.75)
  ) %>%
  ungroup()

dodge <- position_dodge(width = 0.75)
p <- ggplot(combined_data, mapping = aes(x = Variable, y = Value, fill = Group))+ 
  stat_boxplot(mapping = aes(x = Variable, y = Value, fill = Group),
               geom = "errorbar",                            
               width = 0.15, 
               position = position_dodge(0.8)) +    
  geom_boxplot(aes(fill = Group),                             
               position = position_dodge(0.8),                 
               width = 0.6,                                   
               outlier.color = "black") +
  # stat_boxplot(geom = "errorbar", width = 0.5, size = 0.5)+
  scale_fill_manual(values = c("recession" = "grey", "expansion" = "white")) +
  scale_x_discrete(labels = c("SP500" = "S&P500", "PPI_Oil" = "PPI:Oil", "PPI_ACO" = "PPI:All", "ADS" = "ADS")) +
  theme_minimal() +
  labs(title = NULL, x = NULL, y = TeX("$d_{t}")) +
  theme(legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line()) + 
  coord_cartesian(ylim = c(-7, 5))
ggsave("./output/Description/Plots/boxplotDt.png", plot = p, width = 14, height = 8, units = "in", dpi = 600)


# d_{t-1}
df.state1 <- data.frame(
  SP500 = D.std[which(real_state_train == 1) - 1, 1],
  PPI_Oil = D.std[which(real_state_train == 1) - 1, 2],
  PPI_ACO = D.std[which(real_state_train == 1) - 1, 3],
  ADS = D.std[which(real_state_train == 1) - 1, 4]
)

df.state2 <- data.frame(
  SP500 = D.std[which(real_state_train == 1) - 2, 1],
  PPI_Oil = D.std[which(real_state_train == 1) - 2, 2],
  PPI_ACO = D.std[which(real_state_train == 1) - 2, 3],
  ADS = D.std[which(real_state_train == 1) - 2, 4]
)

df1_long <- df.state1 %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = "recession")

df2_long <- df.state2 %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = "expansion")

combined_data <- bind_rows(df1_long, df2_long)
combined_data$Group <- factor(combined_data$Group, levels = c("recession", "expansion"))
combined_data$Variable <- factor(combined_data$Variable, levels = c("SP500", "PPI_Oil", "PPI_ACO", "ADS"))

combined_data_bar <- combined_data %>%
  group_by(Variable, Group) %>%
  summarise(
    ymin = quantile(Value, 0.25) - 1.5 * IQR(Value),
    ymax = quantile(Value, 0.75) + 1.5 * IQR(Value),
    .groups = "drop"
  )

combined_data <- combined_data %>%
  group_by(Variable, Group) %>%
  mutate(
    ymin = quantile(Value, 0.25) - 1.5 * IQR(Value),
    ymax = quantile(Value, 0.75) + 1.5 * IQR(Value)
  ) %>%
  ungroup()

combined_data <- combined_data %>%
  group_by(Variable, Group) %>%
  mutate(
    ydown = quantile(Value, 0.25),
    yup = quantile(Value, 0.75)
  ) %>%
  ungroup()

dodge <- position_dodge(width = 0.75)
p1 <- ggplot(combined_data, mapping = aes(x = Variable, y = Value, fill = Group))+ 
  stat_boxplot(mapping = aes(x = Variable, y = Value, fill = Group),
               geom = "errorbar",                             
               width = 0.15, 
               position = position_dodge(0.8)) +     
  geom_boxplot(aes(fill = Group),                             
               position = position_dodge(0.8),                 
               width = 0.6,                                    
               outlier.color = "black") +
  # stat_boxplot(geom = "errorbar", width = 0.5, size = 0.5)+
  scale_fill_manual(values = c("recession" = "grey", "expansion" = "white")) +
  scale_x_discrete(labels = c("SP500" = "S&P500", "PPI_Oil" = "PPI:Oil", "PPI_ACO" = "PPI:All", "ADS" = "ADS")) +
  theme_minimal() +
  labs(title = NULL, x = NULL, y = TeX("$d_{t-1}")) +
  theme(legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line()) + 
  coord_cartesian(ylim = c(-7, 5))
ggsave("./output/Description/Plots/boxplotDt1.png", plot = p1, width = 14, height = 8, units = "in", dpi = 600)

# d_{t-2}
df.state1 <- data.frame(
  SP500 = D.std[which(real_state_train == 1) - 2, 1],
  PPI_Oil = D.std[which(real_state_train == 1) - 2, 2],
  PPI_ACO = D.std[which(real_state_train == 1) - 2, 3],
  ADS = D.std[which(real_state_train == 1) - 2, 4]
)

df.state2 <- data.frame(
  SP500 = D.std[which(real_state_train == 2)[-1] - 2, 1],
  PPI_Oil = D.std[which(real_state_train == 2)[-1] - 2, 2],
  PPI_ACO = D.std[which(real_state_train == 2)[-1] - 2, 3],
  ADS = D.std[which(real_state_train == 2)[-1] - 2, 4]
)

df1_long <- df.state1 %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = "recession")

df2_long <- df.state2 %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  mutate(Group = "expansion")

combined_data <- bind_rows(df1_long, df2_long)
combined_data$Group <- factor(combined_data$Group, levels = c("recession", "expansion"))
combined_data$Variable <- factor(combined_data$Variable, levels = c("SP500", "PPI_Oil", "PPI_ACO", "ADS"))

combined_data_bar <- combined_data %>%
  group_by(Variable, Group) %>%
  summarise(
    ymin = quantile(Value, 0.25) - 1.5 * IQR(Value),
    ymax = quantile(Value, 0.75) + 1.5 * IQR(Value),
    .groups = "drop"
  )

combined_data <- combined_data %>%
  group_by(Variable, Group) %>%
  mutate(
    ymin = quantile(Value, 0.25) - 1.5 * IQR(Value),
    ymax = quantile(Value, 0.75) + 1.5 * IQR(Value)
  ) %>%
  ungroup()

combined_data <- combined_data %>%
  group_by(Variable, Group) %>%
  mutate(
    ydown = quantile(Value, 0.25),
    yup = quantile(Value, 0.75)
  ) %>%
  ungroup()

dodge <- position_dodge(width = 0.75)
p2 <- ggplot(combined_data, mapping = aes(x = Variable, y = Value, fill = Group))+ 
  stat_boxplot(mapping = aes(x = Variable, y = Value, fill = Group),
               geom = "errorbar",                            
               width = 0.15, 
               position = position_dodge(0.8)) +    
  geom_boxplot(aes(fill = Group),                            
               position = position_dodge(0.8),                 
               width = 0.6,                                   
               outlier.color = "black") +
  # stat_boxplot(geom = "errorbar", width = 0.5, size = 0.5)+
  scale_fill_manual(values = c("recession" = "grey", "expansion" = "white")) +
  scale_x_discrete(labels = c("SP500" = "S&P500", "PPI_Oil" = "PPI:Oil", "PPI_ACO" = "PPI:All", "ADS" = "ADS")) +
  theme_minimal() +
  labs(title = NULL, x = NULL, y = TeX("$d_{t-2}")) +
  theme(legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line()) + 
  coord_cartesian(ylim = c(-7, 5))
ggsave("./output/Description/Plots/boxplotDt2.png", plot = p2, width = 14, height = 8, units = "in", dpi = 600)


