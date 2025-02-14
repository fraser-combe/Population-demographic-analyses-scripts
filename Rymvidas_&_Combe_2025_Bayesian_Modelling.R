Supplementary Material: R Scripts for Bayesian Survival Probability Analysis of Dormice

#########################
# Load necessary libraries
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(rjags)
 
#########################
# Load data must have a 'Sex' column
# Filter data for males and females to create separate datasets for analysis
data_males <- cmr_data %>% filter(Sex == "m")
data_females <- cmr_data %>% filter(Sex == "f")
 
#########################
# Function to generate capture history matrix for given data and study years
# This matrix tracks the presence (1) or absence (0) of each individual per year
generate_capture_history <- function(data, study_years) {
  capture_history_wide <- data %>%
    group_by(`Ring No`, Year) %>%
    summarise(Captured = 1, .groups = "drop") %>%
    pivot_wider(names_from = Year, values_from = Captured, values_fill = list(Captured = 0))
  
  # Ensure all study years are represented in the data
  for (year in study_years) {
    if (!as.character(year) %in% names(capture_history_wide)) {
      capture_history_wide[[as.character(year)]] <- 0
    }
  }
  
  # Order columns by year
  capture_history_wide <- capture_history_wide %>%
    select(`Ring No`, all_of(as.character(study_years)))
  
  return(capture_history_wide)
}
 
#########################
# Generate capture histories for both males and females
capture_history_males <- generate_capture_history(data_males, study_years)
capture_history_females <- generate_capture_history(data_females, study_years)
 
#########################
# Convert processed history data into matrices and setup data for JAGS modeling
capture_matrix_males <- as.matrix(capture_history_males[,-1])   # Exclude 'Ring No'
capture_matrix_females <- as.matrix(capture_history_females[,-1]) # Exclude 'Ring No'
 
# Setup data list for JAGS model; these lists feed into the Bayesian model
data_list_males <- list(
  observed = capture_matrix_males,
  N = nrow(capture_matrix_males),
  T = ncol(capture_matrix_males) - 1,
  initial_survival_prob = 0.6
)
data_list_females <- list(
  observed = capture_matrix_females,
  N = nrow(capture_matrix_females),
  T = ncol(capture_matrix_females) - 1,
  initial_survival_prob = 0.6
)
 
#########################
# Define initial values for the JAGS model based on the data list
# These values are starting points for the Bayesian estimation process
init_values <- function(data_list) {
  list(
    mu_Phi = rnorm(1, 0, 0.1),
    tau_Phi = rgamma(1, 1, 0.1),
    p_global = rep(0.5, data_list$T),
    individual_effect = rnorm(data_list$N, 0, 0.1)
  )
}
 
#########################
# Define the JAGS model
model_code <- "
model {
  mu_Phi ~ dnorm(0, 0.01)  # Normal prior for the average survival probability
  tau_Phi ~ dgamma(0.01, 0.01)  # Gamma prior for the precision of survival probability
 
  for (t in 1:T) {
    p_global[t] ~ dbeta(1, 1)  # Uniform prior for global capture probabilities
  }
 
  for (i in 1:N) {
    individual_effect[i] ~ dnorm(0, tau_Phi)  # Normal distribution for individual effects
 
    Phi[i,1] ~ dbern(initial_survival_prob)  # Initial survival probability
 
    for (t in 2:T+1) {
      logit(Phi[i,t]) <- mu_Phi + individual_effect[i]  # Logit link function for survival probability
      observed[i,t] ~ dbern(Phi[i,t] * p_global[t-1])  # Observation model
    }
  }
}
"
 
# Save and load the model, specify initial values and number of chains
writeLines(model_code, con = "capture_model.jags")
jags_model <- jags.model(file = "capture_model.jags", data = data_list_females, inits = init_values(data_list_females), n.chains = 3, n.adapt = 5000)
 
 
# Run the burn-in phase to reach model convergence
update(jags_model, 5000)
 
# Sample from the posterior distribution of the parameters
samples <- coda.samples(jags_model, variable.names = c("mu_Phi", "tau_Phi", "p_global", "Phi"), n.iter = 10000)

