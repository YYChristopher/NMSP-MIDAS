library("midasr")
library(gtools)
library(data.table)
library(ggplot2)
library(latex2exp)
library(tools)

work_dir = 'F:/24561/Documents/ResearchWrite/NHMM-MIDAS/Submission Requirements/Code in Github'
#Set the Work Space of the data
setwd(work_dir)

N = 1 
NT = 500 
NT_oos = 50
Rep = 100
K = 30

Res_path <- './output/Simulation/Result_compare'
# LASSO-MIDAS
file_name <- list.files(Res_path, pattern = 'LASSO-500Sample.txt', full.names = TRUE)
lasso.res <- readLines(file_name)
rmse.vec <- c(as.numeric(strsplit(lasso.res[1], split = ",")[[1]]))
mae.vec <- c(as.numeric(strsplit(lasso.res[3], split = ",")[[1]]))
crps.vec <- c(as.numeric(strsplit(lasso.res[4], split = ",")[[1]]))
rmse.res = mean(rmse.vec)
mae.res = mean(mae.vec)
crps.res = mean(crps.vec)

# MCP-MIDAS
file_name <- list.files(Res_path, pattern = 'MCP-500Sample.txt', full.names = TRUE)
mcp.res <- readLines(file_name)
rmse.vec <- c(as.numeric(strsplit(mcp.res[1], split = ",")[[1]]))
mae.vec <- c(as.numeric(strsplit(mcp.res[3], split = ",")[[1]]))
crps.vec <- c(as.numeric(strsplit(mcp.res[4], split = ",")[[1]]))
rmse.res = mean(rmse.vec)
mae.res = mean(mae.vec)
crps.res = mean(crps.vec)

# SCAD-MIDAS
file_name <- list.files(Res_path, pattern = 'SCAD-500Sample.txt', full.names = TRUE)
scad.res <- readLines(file_name)
rmse.vec <- c(as.numeric(strsplit(scad.res[1], split = ",")[[1]]))
mae.vec <- c(as.numeric(strsplit(scad.res[3], split = ",")[[1]]))
crps.vec <- c(as.numeric(strsplit(scad.res[4], split = ",")[[1]]))
rmse.res = mean(rmse.vec)
mae.res = mean(mae.vec)
crps.res = mean(crps.vec)

# EN-MIDAS
file_name <- list.files(Res_path, pattern = 'EN-500Sample.txt', full.names = TRUE)
scad.res <- readLines(file_name)
rmse.vec <- c(as.numeric(strsplit(scad.res[1], split = ",")[[1]]))
mae.vec <- c(as.numeric(strsplit(scad.res[3], split = ",")[[1]]))
crps.vec <- c(as.numeric(strsplit(scad.res[4], split = ",")[[1]]))
rmse.res = mean(rmse.vec)
mae.res = mean(mae.vec)
crps.res = mean(crps.vec)

# factor-MIDAS
file_name <- list.files(Res_path, pattern = 'factorMIDAS-500Sample.txt', full.names = TRUE)
scad.res <- readLines(file_name)
rmse.vec <- c(as.numeric(strsplit(scad.res[1], split = ",")[[1]]))
mae.vec <- c(as.numeric(strsplit(scad.res[3], split = ",")[[1]]))
crps.vec <- c(as.numeric(strsplit(scad.res[4], split = ",")[[1]]))
rmse.res = mean(rmse.vec)
mae.res = mean(mae.vec)
crps.res = mean(crps.vec)

# BMA-MIDAS
file_name <- list.files(Res_path, pattern = 'BMAMIDAS-500Sample.txt', full.names = TRUE)
scad.res <- readLines(file_name)
rmse.vec <- c(as.numeric(strsplit(scad.res[1], split = ",")[[1]]))
mae.vec <- c(as.numeric(strsplit(scad.res[3], split = ",")[[1]]))
crps.vec <- c(as.numeric(strsplit(scad.res[4], split = ",")[[1]]))
rmse.res = mean(rmse.vec)
mae.res = mean(mae.vec)
crps.res = mean(crps.vec)

# BMIDAS-AGL
file_name <- list.files(Res_path, pattern = 'BMIDAS-AGL-500Sample.txt', full.names = TRUE)
scad.res <- readLines(file_name)
rmse.vec <- c(as.numeric(strsplit(scad.res[1], split = ",")[[1]]))
mae.vec <- c(as.numeric(strsplit(scad.res[3], split = ",")[[1]]))
crps.vec <- c(as.numeric(strsplit(scad.res[4], split = ",")[[1]]))
rmse.res = mean(rmse.vec)
mae.res = mean(mae.vec)
crps.res = mean(crps.vec)

# BMIDAS-AGL-SS
file_name <- list.files(Res_path, pattern = 'BMIDAS-AGL-SS-500Sample.txt', full.names = TRUE)
scad.res <- readLines(file_name)
rmse.vec <- c(as.numeric(strsplit(scad.res[1], split = ",")[[1]]))
mae.vec <- c(as.numeric(strsplit(scad.res[3], split = ",")[[1]]))
crps.vec <- c(as.numeric(strsplit(scad.res[4], split = ",")[[1]]))
rmse.res = mean(rmse.vec)
mae.res = mean(mae.vec)
crps.res = mean(crps.vec)

# MIDASml
file_name <- list.files(Res_path, pattern = 'MIDASml-500Sample.txt', full.names = TRUE)
scad.res <- readLines(file_name)
rmse.vec <- c(as.numeric(strsplit(scad.res[1], split = ",")[[1]]))
mae.vec <- c(as.numeric(strsplit(scad.res[3], split = ",")[[1]]))
crps.vec <- c(as.numeric(strsplit(scad.res[4], split = ",")[[1]]))
rmse.res = mean(rmse.vec)
mae.res = mean(mae.vec)
crps.res = mean(crps.vec)

# HMSP-MIDAS
file_name <- list.files(Res_path, pattern = 'HMSP-MIDAS_500Sample.txt', full.names = TRUE)
scad.res <- readLines(file_name[1])
rmse.vec <- c(as.numeric(strsplit(scad.res[1], split = ",")[[1]]))
mae.vec <- c(as.numeric(strsplit(scad.res[3], split = ",")[[1]]))
crps.vec <- c(as.numeric(strsplit(scad.res[4], split = ",")[[1]]))
rmse.res = mean(rmse.vec)
mae.res = mean(mae.vec)
crps.res = mean(crps.vec)

# NMSP-MIDAS
file_name <- list.files(Res_path, pattern = 'NMSP-MIDAS_500Sample.txt', full.names = TRUE)
scad.res <- readLines(file_name[1])
rmse.vec <- c(as.numeric(strsplit(scad.res[1], split = ",")[[1]]))
mae.vec <- c(as.numeric(strsplit(scad.res[3], split = ",")[[1]]))
crps.vec <- c(as.numeric(strsplit(scad.res[4], split = ",")[[1]]))
rmse.res = mean(rmse.vec)
mae.res = mean(mae.vec)
crps.res = mean(crps.vec)

# Boxplot
comparison_path <- './output/Simulation/Result_Compare'
# LASSO-MIDAS
file_name <- list.files(comparison_path, pattern = 'LASSO-500Sample.txt', full.names = TRUE)
lasso.MIDAS.res <- readLines(file_name)
all_y_rmse_lasso <- c(as.numeric(strsplit(lasso.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_lasso <- c(as.numeric(strsplit(lasso.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_lasso <- c(as.numeric(strsplit(lasso.MIDAS.res[4], split = ",")[[1]]))

# SCAD-MIDAS
file_name <- list.files(comparison_path, pattern = 'SCAD-500Sample.txt', full.names = TRUE)
scad.MIDAS.res <- readLines(file_name)
all_y_rmse_scad <- c(as.numeric(strsplit(scad.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_scad <- c(as.numeric(strsplit(scad.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_scad <- c(as.numeric(strsplit(scad.MIDAS.res[4], split = ",")[[1]]))

# MCP-MIDAS
file_name <- list.files(comparison_path, pattern = 'MCP-500Sample.txt', full.names = TRUE)
mcp.MIDAS.res <- readLines(file_name)
all_y_rmse_mcp <- c(as.numeric(strsplit(mcp.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_mcp <- c(as.numeric(strsplit(mcp.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_mcp <- c(as.numeric(strsplit(mcp.MIDAS.res[4], split = ",")[[1]]))

# EN-MIDAS
file_name <- list.files(comparison_path, pattern = 'EN-500Sample.txt', full.names = TRUE)
en.MIDAS.res <- readLines(file_name)
all_y_rmse_en <- c(as.numeric(strsplit(en.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_en <- c(as.numeric(strsplit(en.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_en <- c(as.numeric(strsplit(en.MIDAS.res[4], split = ",")[[1]]))

# factor-MIDAS
file_name <- list.files(comparison_path, pattern = 'factorMIDAS-500Sample.txt', full.names = TRUE)
factor.MIDAS.res <- readLines(file_name)
all_y_rmse_factor <- c(as.numeric(strsplit(factor.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_factor <- c(as.numeric(strsplit(factor.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_factor <- c(as.numeric(strsplit(factor.MIDAS.res[4], split = ",")[[1]]))

# BMA-MIDAS
file_name <- list.files(comparison_path, pattern = 'BMAMIDAS-500Sample.txt', full.names = TRUE)
bma.MIDAS.res <- readLines(file_name)
all_y_rmse_bma <- c(as.numeric(strsplit(bma.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_bma <- c(as.numeric(strsplit(bma.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_bma <- c(as.numeric(strsplit(bma.MIDAS.res[4], split = ",")[[1]]))

# BMIDAS-AGL
file_name <- list.files(comparison_path, pattern = 'BMIDAS-AGL-500Sample.txt', full.names = TRUE)
bagl.MIDAS.res <- readLines(file_name)
all_y_rmse_bagl <- c(as.numeric(strsplit(bagl.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_bagl <- c(as.numeric(strsplit(bagl.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_bagl <- c(as.numeric(strsplit(bagl.MIDAS.res[4], split = ",")[[1]]))

# BMIDAS-AGL-SS
file_name <- list.files(comparison_path, pattern = 'BMIDAS-AGL-SS-500Sample.txt', full.names = TRUE)
baglss.MIDAS.res <- readLines(file_name)
all_y_rmse_baglss <- c(as.numeric(strsplit(baglss.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_baglss <- c(as.numeric(strsplit(baglss.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_baglss <- c(as.numeric(strsplit(baglss.MIDAS.res[4], split = ",")[[1]]))

# HMSP-MIDAS
file_name <- list.files(comparison_path, pattern = 'HMSP-MIDAS_500Sample.txt', full.names = TRUE)
hmm.MIDAS.res <- readLines(file_name[1])
all_y_rmse_hmm <- c(as.numeric(strsplit(hmm.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_hmm <- c(as.numeric(strsplit(hmm.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_hmm <- c(as.numeric(strsplit(hmm.MIDAS.res[4], split = ",")[[1]]))

# MIDASml
file_name <- list.files(comparison_path, pattern = 'MIDASml-500Sample.txt', full.names = TRUE)
ml.MIDAS.res <- readLines(file_name)
all_y_rmse_ml <- c(as.numeric(strsplit(ml.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_ml <- c(as.numeric(strsplit(ml.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_ml <- c(as.numeric(strsplit(ml.MIDAS.res[4], split = ",")[[1]]))

# NMSP-MIDAS
file_name <- list.files(comparison_path, pattern = 'NMSP-MIDAS_500Sample.txt', full.names = TRUE)
nhmm.MIDAS.res <- readLines(file_name[1])
all_y_rmse_nhmm <- c(as.numeric(strsplit(nhmm.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_nhmm <- c(as.numeric(strsplit(nhmm.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_nhmm <- c(as.numeric(strsplit(nhmm.MIDAS.res[4], split = ",")[[1]]))

# Plot
rmse.all.df <- data.frame(
  rmse = c(all_y_rmse_nhmm, all_y_rmse_bagl, all_y_rmse_baglss, all_y_rmse_ml, all_y_rmse_hmm, 
           all_y_rmse_lasso, all_y_rmse_mcp, all_y_rmse_scad, all_y_rmse_en, all_y_rmse_factor, 
           all_y_rmse_bma),
  name = factor(c(rep('NMSP-MIDAS', Rep), rep('BMIDAS-AGL', Rep), rep('BMIDAS-AGL-SS', Rep), rep('SG-Lasso-MIDAS', Rep), 
                  rep('HMSP-MIDAS', Rep), rep('LASSO-MIDAS', Rep), rep('MCP-MIDAS', Rep), rep('SCAD-MIDAS', Rep), 
                  rep('EN-MIDAS', Rep), rep('factor-MIDAS', Rep), rep('BMA-MIDAS', Rep)), 
                levels = c('NMSP-MIDAS', 'BMIDAS-AGL', 'BMIDAS-AGL-SS', 'SG-Lasso-MIDAS', 'HMSP-MIDAS', 'LASSO-MIDAS', 
                           'MCP-MIDAS', 'SCAD-MIDAS', 'EN-MIDAS', 'factor-MIDAS', 'BMA-MIDAS'))
)

RMSFE.plot <- rmse.all.df %>%  
  ggplot(aes(x = name, y = rmse, fill = name)) + 
  stat_boxplot(geom = "errorbar", width = 0.35) + 
  geom_boxplot() + 
  geom_boxplot(fill = "white", color = "black") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  xlab("Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("RMSFE")
ggsave("./output/Simulation/Plots/RMSFE.png", plot = RMSFE.plot, width = 14, height = 7, units = "in", dpi = 600)

mae.all.df <- data.frame(
  mae = c(all_y_mae_nhmm, all_y_mae_bagl, all_y_mae_baglss, all_y_mae_ml, all_y_mae_hmm, 
          all_y_mae_lasso, all_y_mae_mcp, all_y_mae_scad, all_y_mae_en, all_y_mae_factor, 
          all_y_mae_bma),
  name = factor(c(rep('NMSP-MIDAS', Rep), rep('BMIDAS-AGL', Rep), rep('BMIDAS-AGL-SS', Rep), rep('SG-Lasso-MIDAS', Rep), 
                  rep('HMSP-MIDAS', Rep), rep('LASSO-MIDAS', Rep), rep('MCP-MIDAS', Rep), rep('SCAD-MIDAS', Rep), 
                  rep('EN-MIDAS', Rep), rep('factor-MIDAS', Rep), rep('BMA-MIDAS', Rep)), 
                levels = c('NMSP-MIDAS', 'BMIDAS-AGL', 'BMIDAS-AGL-SS', 'SG-Lasso-MIDAS', 'HMSP-MIDAS', 'LASSO-MIDAS', 
                           'MCP-MIDAS', 'SCAD-MIDAS', 'EN-MIDAS', 'factor-MIDAS', 'BMA-MIDAS'))
)

MAFE.plot <- mae.all.df %>%  
  ggplot(aes(x = name, y = mae, fill = name)) + 
  stat_boxplot(geom = "errorbar", width = 0.35) + 
  geom_boxplot() + 
  geom_boxplot(fill = "white", color = "black") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  xlab("Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("MAFE")
ggsave("./output/Simulation/Plots/MAFE.png", plot = MAFE.plot, width = 14, height = 7, units = "in", dpi = 600)

crps.all.df <- data.frame(
  crps = c(all_y_crps_nhmm, all_y_crps_bagl, all_y_crps_baglss, all_y_crps_ml, all_y_crps_hmm, 
           all_y_crps_lasso, all_y_crps_mcp, all_y_crps_scad, all_y_crps_en, all_y_crps_factor, 
           all_y_crps_bma),
  name = factor(c(rep('NMSP-MIDAS', Rep), rep('BMIDAS-AGL', Rep), rep('BMIDAS-AGL-SS', Rep), rep('SG-Lasso-MIDAS', Rep), 
                  rep('HMSP-MIDAS', Rep), rep('LASSO-MIDAS', Rep), rep('MCP-MIDAS', Rep), rep('SCAD-MIDAS', Rep), 
                  rep('EN-MIDAS', Rep), rep('factor-MIDAS', Rep), rep('BMA-MIDAS', Rep)), 
                levels = c('NMSP-MIDAS', 'BMIDAS-AGL', 'BMIDAS-AGL-SS', 'SG-Lasso-MIDAS', 'HMSP-MIDAS', 'LASSO-MIDAS', 
                           'MCP-MIDAS', 'SCAD-MIDAS', 'EN-MIDAS', 'factor-MIDAS', 'BMA-MIDAS'))
)

CRPS.plot <- crps.all.df %>%  
  ggplot(aes(x = name, y = crps, fill = name)) + 
  stat_boxplot(geom = "errorbar", width = 0.35) + 
  geom_boxplot() + 
  geom_boxplot(fill = "white", color = "black") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  xlab("Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("CRPS")
ggsave("./output/Simulation/Plots/CRPS.png", plot = CRPS.plot, width = 14, height = 7, units = "in", dpi = 600)


## T = 200
# Boxplot
comparison_path <- './output/Simulation/Result_Compare2'
# with 0.2 SNR
# LASSO-MIDAS
file_name <- list.files(comparison_path, pattern = 'LASSO.txt', full.names = TRUE)
lasso.MIDAS.res <- readLines(file_name)
all_y_rmse_lasso <- c(as.numeric(strsplit(lasso.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_lasso <- c(as.numeric(strsplit(lasso.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_lasso <- c(as.numeric(strsplit(lasso.MIDAS.res[4], split = ",")[[1]]))

# SCAD-MIDAS
file_name <- list.files(comparison_path, pattern = 'SCAD.txt', full.names = TRUE)
scad.MIDAS.res <- readLines(file_name)
all_y_rmse_scad <- c(as.numeric(strsplit(scad.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_scad <- c(as.numeric(strsplit(scad.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_scad <- c(as.numeric(strsplit(scad.MIDAS.res[4], split = ",")[[1]]))

# MCP-MIDAS
file_name <- list.files(comparison_path, pattern = 'MCP.txt', full.names = TRUE)
mcp.MIDAS.res <- readLines(file_name)
all_y_rmse_mcp <- c(as.numeric(strsplit(mcp.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_mcp <- c(as.numeric(strsplit(mcp.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_mcp <- c(as.numeric(strsplit(mcp.MIDAS.res[4], split = ",")[[1]]))

# EN-MIDAS
file_name <- list.files(comparison_path, pattern = 'EN.txt', full.names = TRUE)
en.MIDAS.res <- readLines(file_name)
all_y_rmse_en <- c(as.numeric(strsplit(en.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_en <- c(as.numeric(strsplit(en.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_en <- c(as.numeric(strsplit(en.MIDAS.res[4], split = ",")[[1]]))

# factor-MIDAS
file_name <- list.files(comparison_path, pattern = 'factorMIDAS.txt', full.names = TRUE)
factor.MIDAS.res <- readLines(file_name)
all_y_rmse_factor <- c(as.numeric(strsplit(factor.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_factor <- c(as.numeric(strsplit(factor.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_factor <- c(as.numeric(strsplit(factor.MIDAS.res[4], split = ",")[[1]]))

# BMA-MIDAS
file_name <- list.files(comparison_path, pattern = 'BMAMIDAS.txt', full.names = TRUE)
bma.MIDAS.res <- readLines(file_name)
all_y_rmse_bma <- c(as.numeric(strsplit(bma.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_bma <- c(as.numeric(strsplit(bma.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_bma <- c(as.numeric(strsplit(bma.MIDAS.res[4], split = ",")[[1]]))

# BMIDAS-AGL
file_name <- list.files(comparison_path, pattern = 'BMIDAS-AGL.txt', full.names = TRUE)
bagl.MIDAS.res <- readLines(file_name)
all_y_rmse_bagl <- c(as.numeric(strsplit(bagl.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_bagl <- c(as.numeric(strsplit(bagl.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_bagl <- c(as.numeric(strsplit(bagl.MIDAS.res[4], split = ",")[[1]]))

# BMIDAS-AGL-SS
file_name <- list.files(comparison_path, pattern = 'BMIDAS-AGL-SS.txt', full.names = TRUE)
baglss.MIDAS.res <- readLines(file_name)
all_y_rmse_baglss <- c(as.numeric(strsplit(baglss.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_baglss <- c(as.numeric(strsplit(baglss.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_baglss <- c(as.numeric(strsplit(baglss.MIDAS.res[4], split = ",")[[1]]))

# HMSP-MIDAS
file_name <- list.files(comparison_path, pattern = 'HMSP-MIDAS.txt', full.names = TRUE)
hmm.MIDAS.res <- readLines(file_name[1])
all_y_rmse_hmm <- c(as.numeric(strsplit(hmm.MIDAS.res[2], split = ",")[[1]]))
all_y_mae_hmm <- c(as.numeric(strsplit(hmm.MIDAS.res[4], split = ",")[[1]]))
all_y_crps_hmm <- c(as.numeric(strsplit(hmm.MIDAS.res[5], split = ",")[[1]]))

# MIDASml
file_name <- list.files(comparison_path, pattern = 'MIDASml.txt', full.names = TRUE)
ml.MIDAS.res <- readLines(file_name)
all_y_rmse_ml <- c(as.numeric(strsplit(ml.MIDAS.res[1], split = ",")[[1]]))
all_y_mae_ml <- c(as.numeric(strsplit(ml.MIDAS.res[3], split = ",")[[1]]))
all_y_crps_ml <- c(as.numeric(strsplit(ml.MIDAS.res[4], split = ",")[[1]]))

# NMSP-MIDAS
file_name <- list.files(comparison_path, pattern = 'NMSP-MIDAS.txt', full.names = TRUE)
nhmm.MIDAS.res <- readLines(file_name[1])
all_y_rmse_nhmm <- c(as.numeric(strsplit(nhmm.MIDAS.res[2], split = ",")[[1]]))
all_y_mae_nhmm <- c(as.numeric(strsplit(nhmm.MIDAS.res[4], split = ",")[[1]]))
all_y_crps_nhmm <- c(as.numeric(strsplit(nhmm.MIDAS.res[5], split = ",")[[1]]))

# Plot
rmse.all.df <- data.frame(
  rmse = c(all_y_rmse_nhmm, all_y_rmse_bagl, all_y_rmse_baglss, all_y_rmse_ml, all_y_rmse_hmm, 
           all_y_rmse_lasso, all_y_rmse_mcp, all_y_rmse_scad, all_y_rmse_en, all_y_rmse_factor, 
           all_y_rmse_bma),
  name = factor(c(rep('NMSP-MIDAS', Rep), rep('BMIDAS-AGL', Rep), rep('BMIDAS-AGL-SS', Rep), rep('SG-Lasso-MIDAS', Rep), 
                  rep('HMSP-MIDAS', Rep), rep('LASSO-MIDAS', Rep), rep('MCP-MIDAS', Rep), rep('SCAD-MIDAS', Rep), 
                  rep('EN-MIDAS', Rep), rep('factor-MIDAS', Rep), rep('BMA-MIDAS', Rep)), 
                levels = c('NMSP-MIDAS', 'BMIDAS-AGL', 'BMIDAS-AGL-SS', 'SG-Lasso-MIDAS', 'HMSP-MIDAS', 'LASSO-MIDAS', 
                           'MCP-MIDAS', 'SCAD-MIDAS', 'EN-MIDAS', 'factor-MIDAS', 'BMA-MIDAS'))
)

RMSFE.plot2 <- rmse.all.df %>%  
  ggplot(aes(x = name, y = rmse, fill = name)) + 
  stat_boxplot(geom = "errorbar", width = 0.35) + 
  geom_boxplot() + 
  geom_boxplot(fill = "white", color = "black") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  xlab("Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("RMSFE")
ggsave("./output/Simulation/Plots/RMSFE2.png", plot = RMSFE.plot2, width = 14, height = 7, units = "in", dpi = 600)

mae.all.df <- data.frame(
  mae = c(all_y_mae_nhmm, all_y_mae_bagl, all_y_mae_baglss, all_y_mae_ml, all_y_mae_hmm, 
          all_y_mae_lasso, all_y_mae_mcp, all_y_mae_scad, all_y_mae_en, all_y_mae_factor, 
          all_y_mae_bma),
  name = factor(c(rep('NMSP-MIDAS', Rep), rep('BMIDAS-AGL', Rep), rep('BMIDAS-AGL-SS', Rep), rep('SG-Lasso-MIDAS', Rep), 
                  rep('HMSP-MIDAS', Rep), rep('LASSO-MIDAS', Rep), rep('MCP-MIDAS', Rep), rep('SCAD-MIDAS', Rep), 
                  rep('EN-MIDAS', Rep), rep('factor-MIDAS', Rep), rep('BMA-MIDAS', Rep)), 
                levels = c('NMSP-MIDAS', 'BMIDAS-AGL', 'BMIDAS-AGL-SS', 'SG-Lasso-MIDAS', 'HMSP-MIDAS', 'LASSO-MIDAS', 
                           'MCP-MIDAS', 'SCAD-MIDAS', 'EN-MIDAS', 'factor-MIDAS', 'BMA-MIDAS'))
)

MAFE.plot2 <- mae.all.df %>%  
  ggplot(aes(x = name, y = mae, fill = name)) + 
  stat_boxplot(geom = "errorbar", width = 0.35) + 
  geom_boxplot() + 
  geom_boxplot(fill = "white", color = "black") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  xlab("Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("MAFE")
ggsave("./output/Simulation/Plots/MAFE2.png", plot = MAFE.plot2, width = 14, height = 7, units = "in", dpi = 600)

crps.all.df <- data.frame(
  crps = c(all_y_crps_nhmm, all_y_crps_bagl, all_y_crps_baglss, all_y_crps_ml, all_y_crps_hmm, 
           all_y_crps_lasso, all_y_crps_mcp, all_y_crps_scad, all_y_crps_en, all_y_crps_factor, 
           all_y_crps_bma),
  name = factor(c(rep('NMSP-MIDAS', Rep), rep('BMIDAS-AGL', Rep), rep('BMIDAS-AGL-SS', Rep), rep('SG-Lasso-MIDAS', Rep), 
                  rep('HMSP-MIDAS', Rep), rep('LASSO-MIDAS', Rep), rep('MCP-MIDAS', Rep), rep('SCAD-MIDAS', Rep), 
                  rep('EN-MIDAS', Rep), rep('factor-MIDAS', Rep), rep('BMA-MIDAS', Rep)), 
                levels = c('NMSP-MIDAS', 'BMIDAS-AGL', 'BMIDAS-AGL-SS', 'SG-Lasso-MIDAS', 'HMSP-MIDAS', 'LASSO-MIDAS', 
                           'MCP-MIDAS', 'SCAD-MIDAS', 'EN-MIDAS', 'factor-MIDAS', 'BMA-MIDAS'))
)

CRPS.plot2 <- crps.all.df %>%  
  ggplot(aes(x = name, y = crps, fill = name)) + 
  stat_boxplot(geom = "errorbar", width = 0.35) + 
  geom_boxplot() + 
  geom_boxplot(fill = "white", color = "black") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  xlab("Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("CRPS")
ggsave("./output/Simulation/Plots/CRPS2.png", plot = CRPS.plot2, width = 14, height = 7, units = "in", dpi = 600)