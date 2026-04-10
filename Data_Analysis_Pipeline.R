# ==============================================================================
# Title: Bayesian Network Analysis and Simulation for SCL-90
# Description: This script contains the full analytical pipeline, including 
#              structure learning, parameter extraction, network visualization, 
#              GGM-based structural attenuation, and Bayesian do-simulations.
# Note: Raw mental health data is not provided due to privacy and ethical 
#       restrictions. Users can substitute "your_dataset.sav" with their own data.
# ==============================================================================

# ------------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(haven)     # For reading data
  library(dplyr)     # For data manipulation
  library(bnlearn)   # For Bayesian Network learning and do-calculus
  library(qgraph)    # For network visualization
  library(corpcor)   # For shrinkage covariance in GGM
  library(ggplot2)   # For plotting simulation results
})

set.seed(20251012) # Set global seed for reproducibility

# ------------------------------------------------------------------------------
# MODULE 1: Data Preparation & Structure Learning
# ------------------------------------------------------------------------------
# Load the pre-processed analytic dataset
# (Raw data is subsetted and anonymized to protect participant privacy)

data2022 <- read_sav("SCL90_Final_Analysis_Data.sav") 
node_names <- c("SOM", "O-C", "INT", "DEP", "ANX", "HOS", "PHOB", "PAR", "PSY")
colnames(data2022) <- node_names

# 1.1 Bootstrap robust structure learning (R = 1000)
boot_res <- boot.strength(
  data = data2022,
  R = 1000,
  algorithm = "hc",
  algorithm.args = list(score = "bic-g")
)

# 1.2 Average network construction
thr <- 0.85  # Inclusion probability threshold
bn_avg <- averaged.network(boot_res, threshold = thr)

# ------------------------------------------------------------------------------
# MODULE 2: Parameter Extraction (Standardized Coefficients for Table S1)
# ------------------------------------------------------------------------------
# 2.1 Standardize data for path coefficient extraction
data_std <- as.data.frame(scale(data2022))

# 2.2 Fit the averaged DAG to standardized data
fit_std <- bn.fit(bn_avg, data_std)
final_arcs <- as.data.frame(arcs(bn_avg))
final_arcs$Beta <- NA

# 2.3 Extract standardized coefficients (Beta) for each directed edge
for (i in 1:nrow(final_arcs)) {
  from_node <- final_arcs$from[i]
  to_node <- final_arcs$to[i]
  # Extract coefficient of the parent node from the child's regression equation
  final_arcs$Beta[i] <- fit_std[[to_node]]$coefficients[from_node]
}

# 2.4 Merge with bootstrap inclusion probabilities and format
table_s1 <- merge(final_arcs, boot_res, by = c("from", "to"), all.x = TRUE)
colnames(table_s1) <- c("From", "To", "Standardized_Beta", "Strength", "Direction")
table_s1 <- table_s1[order(table_s1$From, -table_s1$Strength), ]
rownames(table_s1) <- NULL

write.csv(table_s1, "Table_S1_Path_Coefficients.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# MODULE 3: Network Visualization (qgraph - Figure 1)
# ------------------------------------------------------------------------------
# 3.1 Define hierarchical layout (Waterfall flow)
flow_layout <- matrix(c(
  0.00, -0.70,  # 1. SOM 
  -1.00, -0.35,  # 2. O-C 
  0.70,  0.70,  # 3. INT 
  0.00,  1.00,  # 4. DEP 
  -0.70,  0.40,  # 5. ANX 
  0.50, -0.35,  # 6. HOS 
  -0.80, -1.00,  # 7. PHOB
  0.90,  0.10,  # 8. PAR 
  0.00,  0.10   # 9. PSY 
), ncol = 2, byrow = TRUE)

# 3.2 Generate high-fidelity network plot
pdf("Figure1_Directed_Network.pdf", width = 8, height = 8)
final_network_plot <- qgraph(
  arcs(bn_avg),
  layout = flow_layout,
  directed = TRUE,
  curveAll = TRUE, curve = 0, fade = FALSE, edge.width = 1.2,
  posCol = "#0000EF", negCol = "#C43C39",
  shape = "circle", vsize = 6,
  borders = TRUE, border.color = "grey20", border.width = 1.5,
  labels = node_names, label.cex = 1, label.scale = FALSE
)
dev.off()

# ------------------------------------------------------------------------------
# MODULE 4: Structural Flexibility Analysis (GGM & Virtual Attenuation)
# ------------------------------------------------------------------------------
# 4.1 Fit baseline GGM via shrinkage covariance
X <- as.matrix(data2022)
S <- corpcor::cov.shrink(X, verbose = FALSE)
K_base <- corpcor::pseudoinverse(S)

pcor_from_precision <- function(K) {
  d <- diag(K)
  P <- -K / sqrt(outer(d, d))
  diag(P) <- 1; return(P)
}
global_strength <- function(P) sum(abs(P[upper.tri(P, diag = FALSE)]))

GS_base <- global_strength(pcor_from_precision(K_base))

# 4.2 Function to simulate intervention (attenuation)
intervene_precision <- function(K, j, factor) {
  K2 <- K; K2[j, j] <- K[j, j] * factor; return(K2)
}

# 4.3 Run sensitivity analysis
parent_candidates <- c("DEP", "INT")
intervene_factors <- c(1.5, 2.0, 3.0)
results_list <- list()

for (tg in parent_candidates) {
  j <- match(tg, node_names)
  for (fac in intervene_factors) {
    K_int <- intervene_precision(K_base, j, factor = fac)
    GS_int <- global_strength(pcor_from_precision(K_int))
    results_list[[length(results_list) + 1]] <- data.frame(
      target_node = tg, factor = fac,
      GS_base = GS_base, GS_after = GS_int,
      dGS = GS_int - GS_base
    )
  }
}
res_df <- dplyr::bind_rows(results_list)
res_df$GS_base_ratio <- (res_df$dGS / res_df$GS_base) * 100

# 4.4 Plot Structural Flexibility
ggplot(res_df, aes(x = factor, y = GS_base_ratio, color = target_node, group = target_node)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  labs(x = "Attenuation factor (Kjj multiplier)", y = expression(Delta*" Global Strength (%)"), color = "Target node") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

# ------------------------------------------------------------------------------
# MODULE 5: Functional Impact via Do-calculus Simulations
# ------------------------------------------------------------------------------
# 5.1 Baseline expectation
fit_base <- bn.fit(bn_avg, data = data2022, method = "mle-g")
simulate_expected <- function(fit, nodes, n = 50000L) {
  sim <- rbn(fit, n = n)
  list(per_node_mean = colMeans(sim[, nodes]), total_mean = mean(rowSums(sim[, nodes])))
}
base_res <- simulate_expected(fit_base, node_names, n = 50000L)

# 5.2 Do-calculus function
do_gaussian <- function(dag, data, target, nodes, shift_sd = -1, n = 50000L) {
  dag2 <- dag
  par_t <- parents(dag2, node = target)
  if (length(par_t) > 0) for (p in par_t) dag2 <- drop.arc(dag2, from = p, to = target)
  
  fit_tmp <- bn.fit(dag2, data = data, method = "mle-g")
  new_mu <- mean(data[[target]], na.rm = TRUE) + shift_sd * sd(data[[target]], na.rm = TRUE)
  new_sd <- fit_tmp[[target]]$sd
  
  spec <- vector("list", length(nodes)); names(spec) <- nodes
  for (v in nodes) {
    if (v == target) {
      spec[[v]] <- list(coef = c("(Intercept)" = new_mu), sd = new_sd)
    } else {
      spec[[v]] <- list(coef = fit_tmp[[v]]$coef, sd = fit_tmp[[v]]$sd)
    }
  }
  fit_do <- custom.fit(dag2, dist = spec)
  sim <- rbn(fit_do, n = n)
  
  list(total_mean = mean(rowSums(sim[, nodes])), per_node_mean = colMeans(sim[, nodes]))
}

# 5.3 Run do-simulations
shift_grid <- c(-0.5, -1.0, -1.5)
res_tot <- list()
k <- 0

for (tg in parent_candidates) {
  for (sh in shift_grid) {
    k <- k + 1
    ans <- do_gaussian(dag = bn_avg, data = data2022, target = tg, nodes = node_names, shift_sd = sh, n = 50000L)
    res_tot[[k]] <- data.frame(
      target_node = tg, shift_sd = sh,
      d_total = ans$total_mean - base_res$total_mean,
      rel_change_pct = (ans$total_mean - base_res$total_mean) / base_res$total_mean * 100
    )
  }
}
res_do_total <- bind_rows(res_tot)

# 5.4 Plot Functional Impact
ggplot(res_do_total, aes(x = abs(shift_sd), y = rel_change_pct, color = target_node, group = target_node)) +
  geom_line(linewidth = 1.1) + geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  labs(x = "Intervention intensity (SD shift)", y = expression(Delta*" Total Load (%)"), color = "Target node") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

# ================================ END OF SCRIPT ===============================