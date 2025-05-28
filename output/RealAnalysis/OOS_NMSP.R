library(data.table)
library(ggplot2)
library(latex2exp)

work_dir = 'F:/24561/Documents/ResearchWrite/NHMM-MIDAS/Submission Requirements/Code in Github'
#Set the Work Space of the data
setwd(work_dir)

NT_test <- 51
NS <- 2

GDP = read.csv('./data/GDP2.csv', header = TRUE)
GDP_noindex <- GDP['GDP']
dates_gdp <- seq(as.Date("1947-03-01"), by = "quarter", length.out = 308)
GDP_xts <- xts(GDP_noindex, order.by = dates_gdp)
gdp <- GDP_xts['1959-03-01/2023-12-01']
gdp.test <- GDP_xts['2011-03-01/2023-12-01']
len.test = as.numeric(length(gdp.test))
gdp_rate_test = rep(0, len.test - 1)
for (i in 1:(len.test - 1)) {
  gdp_rate_test[i] = 100 * log(as.numeric(gdp.test[i + 1]) / as.numeric(gdp.test[i]))
}
growth_dates_test = dates_test <- seq(as.Date("2011-03-01"), by = "quarter", length.out = len.test - 1)
gdp_growth_test = xts(gdp_rate_test, order.by = growth_dates_test)
gdp_growth_test_num = as.numeric(gdp_growth_test)
gdp_growth.test <- gdp_growth_test["2011-01-01/2023-12-01"]
Y.test <- gdp_growth.test

predmat = read.table("./output/RealAnalysis/PredState/predsave_NMSP.txt", header = FALSE, sep = " ")
state.forecast.res <- matrix(0, nrow = 1, ncol = NT_test)
for(j in 1:NT_test){
  state.forecast.res[1, j] = which.max(predmat[j, ])
}

m_mon <- 3
NBER = read.csv('./data/USREC.csv', header = TRUE)
NBER_noindex <- NBER['USREC']
dates_nber_all <- seq(as.Date("1854-12-01"), by = "month", length.out = dim(NBER)[1])
NBER_xts <- xts(NBER_noindex, order.by = dates_nber_all)
dates_nber.quarter <- seq(as.Date("2011-01-01"), by = "quarter", length.out = NT_test)
dates_nber.month <- seq(as.Date("2011-01-01"), by = "month", length.out = NT_test * m_mon)
NBER_used <- NBER_xts['2011-01-01/2023-09-01']

df.nber <- data.frame(
  Date = dates_nber.month, 
  Recession = c(as.numeric(NBER_used))
)

recession_periods <- subset(df.nber, NBER == 1)
df.nber$Date <- as.Date(df.nber$Date)

# Define begin and over
recession_periods <- df.nber %>%
  mutate(Recession_Change = Recession != lag(Recession, default = 0)) %>%
  mutate(Period_ID = cumsum(Recession_Change)) %>%
  group_by(Period_ID) %>%
  filter(Recession == 1) %>%
  summarize(Start = min(Date), End = max(Date))

# define color
df.gdp.growth.color <- data.frame(
  Date = dates_nber.quarter,
  GDP.growth = c(Y.test),
  state = as.character(c(state.forecast.res))
)

plot.gdp.growth.test.color <- ggplot() +
  geom_rect(data = recession_periods, aes(xmin = Start, xmax = End, ymin = -9, ymax = 9),
            fill = "grey", alpha = 0.5) + 
  geom_point(data = df.gdp.growth.color, aes(x = Date, y = GDP.growth, color = state)) +
  scale_color_manual(values = c("1" = "red", "2" = "blue")) +
  labs(y = TeX("$y_t"), x = "Date") + 
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16), 
    axis.text = element_text(size = 12),     
    plot.title = element_text(size = 18, face = "bold"),  
    legend.text = element_text(size = 14),   
    legend.title = element_text(size = 16)  
  ) +
  scale_x_date(date_breaks = "2 year", date_labels = "%Y") +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(axis.text.y = element_text())  
show(plot.gdp.growth.test.color)
ggsave("./output/RealAnalysis/Plots/NMSP-OUT.png", plot = plot.gdp.growth.test.color, width = 12, height = 7, units = "in", dpi = 600)
