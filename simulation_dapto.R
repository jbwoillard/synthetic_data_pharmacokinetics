# simulations Alex modele simple vanco pour article POPPK avatar

#Loading of the packages


library(mrgsolve)
library(tidyverse)
library(mapbayr)
library(MESS)
library(truncnorm)
library(FNN)
library(ggcorrplot)
library(GGally)

# Exemple modele 1 daptomycine dvorchik


code <-
  "
$PARAM @annotated
TVCL:  0.807 : Clearance
TVV1: 4.80 : Central volume
TVV2  : 3.13 : Peripheral volume of distribution
TVQ   :  3.46 : Intercompartmental clearance
CRCL: 0.00514 : effect of CLCR on CL
TCL: 0.14 : effect of body temperature on CL
WTQ: 0.0593 : linear effect of WT deviation on Q
WTV2: 0.0458 : linear effect of WT deviation on V2

ETA1: 0 : Clearance (L/h)
ETA2: 0 : central volume (L)
ETA3: 0 : intercompartmental Clearance (L/h)
ETA4: 0 : peripheral volume (L)

$PARAM @annotated @covariates
CREATCL : 91.2 : estimated creat clearance ml/min
TEMP: 37.2 : Body temperature
WT: 75.1 : Body weight
SEX : 1 : men (1) women(0)


$OMEGA 0.093636 0.3249 0.425104 0.036481

$SIGMA
0.000001 // proportional
2.068 // additive

$CMT @annotated
CENT  : Central compartment (mg/L)[ADM, OBS]
PERIPH: Peripheral compartment ()

$TABLE
capture DV = (CENT/V1) *(1 + EPS(1)) + EPS(2);

int i = 0;
while(DV<0 && i <100) {
simeps();
DV = (CENT/V1) *(1 + EPS(1)) + EPS(2);
++i;
}

$MAIN
double CL = ((TVCL+ CRCL*(CREATCL - 91.2) + TCL*(TEMP - 37.2))* (0.8 + 0.2*SEX)) * exp(ETA1 + ETA(1)) ;
double V1 = TVV1  * exp(ETA2 + ETA(2)) ;
double Q = TVQ + WTQ * (WT - 75.1) * exp(ETA3 + ETA(3)) ;
double V2 = TVV2 + WTV2 * (WT - 75.1) * exp(ETA4 + ETA(4)) ;

$ODE
dxdt_CENT   =  -(CL+Q)*CENT/V1 + Q*PERIPH/V2 ;
dxdt_PERIPH =  Q*CENT/V1 - Q*PERIPH/V2 ;

$CAPTURE DV CL V1 V2 Q
"
dapto_dvorchik <- mcode("dapto_dvorchik", code)


#Simulation patient 


set.seed(1234)
WT_data = tibble(ID = 1:50) %>% mutate(WT = rtruncnorm(n(),a=48, b=153, mean=75.1, sd=30))
set.seed(1234)
CREATCL_data = tibble(ID = 1:50) %>% mutate(CREATCL = rtruncnorm(n(),a=14, b=150, mean=91.2, sd=30))
set.seed(1234)
temp_data = tibble(ID = 1:50) %>% mutate(TEMP  = rtruncnorm(n(),a=36.1, b=40.1, mean=37.2, sd=1))
set.seed(1234)
SEX_data = tibble(ID = 1:50) %>% mutate(SEX = rbinom(n(),1,0.59)) 

data_ev <- as_tibble(ev(ID = 1:50, amt = 10, ii=24, addl=5, ss=1, tinf=0.5)) %>%
  arrange(ID) %>%
  left_join(WT_data) %>%
  left_join(CREATCL_data) %>%
  left_join(SEX_data) %>% 
  left_join(temp_data) %>% 
  mutate(amt = WT*amt, rate = amt/tinf) 

set.seed(12)
out_ex <- dapto_dvorchik %>% 
  data_set(data_ev) %>%
  mrgsim(end = 48, tgrid = c(0,0.5, 1, 2, 3, 4, 6, 8, 12, 15, 24, 24.5, 25, 26, 27, 28, 30, 32, 36, 39, 48) )

#Simulations graphique
simul_datpo_all <- as_tibble(out_ex) %>% 
  select(-CENT, -PERIPH) %>% 
  left_join(data_ev, join_by(ID)) %>% 
  rename(time = time.x) %>% 
  filter(between(time, 24, 48))

simul_datpo_all %>% ggplot(aes(x = time, y = DV, color = as.factor(ID))) + geom_point() + geom_line() + labs(x = "Time (h)", y = "Daptomycin concentration (mg/L)") + theme_bw() + theme(legend.position = "none")

write_csv(simul_datpo_all, file = "fichier-original-simulation_dapto_param_pk.csv")

# avatar on wide data

# wide format for avatar

original_pk_wide_dapto <- simul_datpo_all %>% 
  select(-CL:-Q, -time.y, -ID) %>% 
  pivot_wider( names_from = "time", values_from = "DV", names_prefix = "conc_") %>% 
  select(-ii:-tinf)

write_csv(original_pk_wide_dapto, file = "fichier-original-simulation_dapto_param_pk_wide.csv")
## Algorithm with knn=5

data_normalized <- scale(original_pk_wide_dapto)
pca <- prcomp(data_normalized, scale. = FALSE)# pour selecitonner le nombre de cp rank. = 3
# Number of neighbors
k <- 5  # Adjust this based on your requirement
pca_transformed_data <- pca$x
knn_result <- get.knn(pca_transformed_data, k)

generate_avatar_weights <- function(knn_result, pca_transformed_data, k) {
  n <- nrow(pca_transformed_data)
  avatar_weights <- matrix(nrow = n, ncol = k)
  
  for (i in 1:n) {
    # Step 1: Inverse of Distances
    distances <- sqrt(rowSums((pca_transformed_data[knn_result$nn.index[i, ], ] - pca_transformed_data[i, ])^2))
    inverse_distances <- 1 / distances
    
    # Step 2: Random Weights
    set.seed(12345)
    random_weights <- rexp(k, rate = 1)
    
    # Step 3: Contribution Factors
    set.seed(12345)
    shuffled_indices <- sample(k)
    contribution_factors <- 1 / (2^shuffled_indices)
    
    # Step 4: Calculate Weights
    weights <- inverse_distances * random_weights * contribution_factors
    
    # Step 5: Normalize Weights
    normalized_weights <- weights / sum(weights)
    
    avatar_weights[i, ] <- normalized_weights
  }
  
  return(avatar_weights)
}



# Generate avatar weights
avatar_weights <- generate_avatar_weights(knn_result, pca_transformed_data, k)

# Assuming pca_result, avatar_weights, knn_result$nn.index, and pca_transformed_data are already defined

# Function to generate avatars in PCA space based on weights
generate_avatars_pca_space <- function(pca_transformed_data, knn_indices, weights) {
  n <- nrow(pca_transformed_data)
  avatars_pca <- matrix(nrow = n, ncol = ncol(pca_transformed_data))
  
  for (i in 1:n) {
    weighted_avatars <- pca_transformed_data[knn_indices[i, ], ] * weights[i, ]
    avatars_pca[i, ] <- colSums(weighted_avatars)
  }
  
  return(avatars_pca)
}
# Generate avatars in PCA space
avatars_pca_space <- generate_avatars_pca_space(pca_transformed_data, knn_result$nn.index, avatar_weights)

# Assuming 'aids_pca' is the PCA object and 'avatars_pca_space' contains the avatars in PCA space
# Inverse PCA transformation
inverse_pca <- function(pca_object, pca_data) {
  return(pca_data %*% t(pca_object$rotation) + matrix(pca_object$center, nrow = nrow(pca_data), ncol = ncol(pca_object$rotation), byrow = TRUE))
}
avatars_original_scale <- inverse_pca(pca, avatars_pca_space)

# Assuming 'aids_data_normalized' contains the scaling attributes of the original data
# Inverse normalization (if the original data was normalized)
avatars_rescaled <- scale(avatars_original_scale, center = FALSE, scale = 1/attr(data_normalized, "scaled:scale"))
avatars_rescaled <- sweep(avatars_rescaled, 2, attr(data_normalized, "scaled:center"), "+")

avatars_tibble_knn5 <- as_tibble(avatars_rescaled) %>% 
  mutate(SEX   = if_else(SEX>0.7 ,1,0))


#Plot of the synthetic and original in the latent space


# Combine original and synthetic data for visualization
combined_data <- rbind(
  original_pk_wide_dapto %>% mutate(DataType = 'Original'),
  avatars_tibble_knn5 %>% mutate(DataType = 'Synthetic')
)

# Perform PCA on combined data
combined_data_normalized <- scale(combined_data[, -which(names(combined_data) %in% c("DataType", "id"))])
combined_pca <- prcomp(combined_data_normalized, scale. = FALSE)

# Extract the first two principal components
combined_pca_data <- data.frame(combined_pca$x[, 1:2])
combined_pca_data$DataType <- combined_data$DataType

# Plot PCA with color differentiation
ggplot(combined_pca_data, aes(x = PC1, y = PC2, color = DataType)) +
  geom_point(alpha = 1) +
  theme_minimal() +
  labs(title = "PCA Plot", x = "Principal Component 1", y = "Principal Component 2", color = "Data Type")

## Compare the data

### summary distribution

summary(avatars_tibble_knn5)
summary(original_pk_wide_dapto)

### plot correlation

##Correlation Analysis
cor_real <- cor(original_pk_wide_dapto, use = "complete.obs")
cor_synthetic <- cor(avatars_tibble_knn5, use = "complete.obs")

# plots
ggcorrplot(cor_real, hc.order = TRUE, type = "lower",
           lab = TRUE,  pch.cex = 5,
           tl.cex = 6, lab_size = 2)
# plots
ggcorrplot(cor_synthetic, hc.order = TRUE, type = "lower",
           lab = TRUE,  pch.cex = 5,
           tl.cex = 6, lab_size = 2)

# ggpair compairons
group_comparison <-
  original_pk_wide_dapto %>% 
  mutate(group="original") %>% 
  bind_rows(avatars_tibble_knn5 %>% 
              mutate(group="synthetic")) 

pm <- group_comparison %>% ggpairs(
  ggplot2::aes(colour = group,alpha = 0.5),
  upper = list(continuous = wrap("cor", size = 1.5)),
  lower=list(combo=wrap("facethist", binwidth=0.5))) + 
  theme(strip.text.x = element_text(size = 5),
        strip.text.y = element_text(size = 5),axis.text = element_text(size = 5))
pm


#Process and plot the data

avatar_knn5 <- avatars_tibble_knn5 %>%
  mutate(ID = 1:nrow(avatars_tibble_knn5)) %>% 
  pivot_longer(cols = starts_with("conc_"), names_to = "time", values_to = "DV", names_prefix = "conc_") %>%
  mutate(time = as.numeric(time)) %>% 
  arrange(ID, time) %>% 
  relocate(ID) %>% 
  mutate(ii = 24, addl = 5, cmt = 1, ss=1, tinf = 0.5, evid = 0) 

avatar_knn5 <- avatar_knn5 %>% 
 bind_rows(avatar_knn5 %>% group_by(ID) %>% slice_head(n=1) %>% mutate(time=0, evid=1))%>% 
  arrange(ID, time) 
write_csv(avatar_knn5, file = "avatar_simul_dapto_knn5.csv")

#Plot the data

simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(avatar_knn5 %>% mutate(group="synthetic")) %>% 
  ggplot(aes(x = time, y = DV, color = group)) + geom_point() +
  geom_point(alpha = 0.2, size = 3) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw()


simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(avatar_knn5 %>% 
              mutate(group="synthetic")) %>% 
  mutate(ID = as.factor(ID)) %>% 
  ggplot(aes(x = time, y = DV, color = ID)) +  geom_line(show.legend = FALSE)+
  geom_point(alpha = 0.2, size = 3, show.legend = FALSE) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw() + facet_wrap(~group)

## Algorithm with knn=15

data_normalized <- scale(original_pk_wide_dapto)
pca <- prcomp(data_normalized, scale. = FALSE)# pour selecitonner le nombre de cp rank. = 3
# Number of neighbors
k <- 15  # Adjust this based on your requirement
pca_transformed_data <- pca$x
knn_result <- get.knn(pca_transformed_data, k)

generate_avatar_weights <- function(knn_result, pca_transformed_data, k) {
  n <- nrow(pca_transformed_data)
  avatar_weights <- matrix(nrow = n, ncol = k)
  
  for (i in 1:n) {
    # Step 1: Inverse of Distances
    distances <- sqrt(rowSums((pca_transformed_data[knn_result$nn.index[i, ], ] - pca_transformed_data[i, ])^2))
    inverse_distances <- 1 / distances
    
    # Step 2: Random Weights
    set.seed(12345)
    random_weights <- rexp(k, rate = 1)
    
    # Step 3: Contribution Factors
    set.seed(12345)
    shuffled_indices <- sample(k)
    contribution_factors <- 1 / (2^shuffled_indices)
    
    # Step 4: Calculate Weights
    weights <- inverse_distances * random_weights * contribution_factors
    
    # Step 5: Normalize Weights
    normalized_weights <- weights / sum(weights)
    
    avatar_weights[i, ] <- normalized_weights
  }
  
  return(avatar_weights)
}



# Generate avatar weights
avatar_weights <- generate_avatar_weights(knn_result, pca_transformed_data, k)

# Assuming pca_result, avatar_weights, knn_result$nn.index, and pca_transformed_data are already defined

# Function to generate avatars in PCA space based on weights
generate_avatars_pca_space <- function(pca_transformed_data, knn_indices, weights) {
  n <- nrow(pca_transformed_data)
  avatars_pca <- matrix(nrow = n, ncol = ncol(pca_transformed_data))
  
  for (i in 1:n) {
    weighted_avatars <- pca_transformed_data[knn_indices[i, ], ] * weights[i, ]
    avatars_pca[i, ] <- colSums(weighted_avatars)
  }
  
  return(avatars_pca)
}
# Generate avatars in PCA space
avatars_pca_space <- generate_avatars_pca_space(pca_transformed_data, knn_result$nn.index, avatar_weights)

# Assuming 'aids_pca' is the PCA object and 'avatars_pca_space' contains the avatars in PCA space
# Inverse PCA transformation
inverse_pca <- function(pca_object, pca_data) {
  return(pca_data %*% t(pca_object$rotation) + matrix(pca_object$center, nrow = nrow(pca_data), ncol = ncol(pca_object$rotation), byrow = TRUE))
}
avatars_original_scale <- inverse_pca(pca, avatars_pca_space)

# Assuming 'aids_data_normalized' contains the scaling attributes of the original data
# Inverse normalization (if the original data was normalized)
avatars_rescaled <- scale(avatars_original_scale, center = FALSE, scale = 1/attr(data_normalized, "scaled:scale"))
avatars_rescaled <- sweep(avatars_rescaled, 2, attr(data_normalized, "scaled:center"), "+")

avatars_tibble_knn15 <- as_tibble(avatars_rescaled) %>% 
  mutate(SEX   = if_else(SEX>0.7 ,1,0))


#Plot of the synthetic and original in the latent space


# Combine original and synthetic data for visualization
combined_data <- rbind(
  original_pk_wide_dapto %>% mutate(DataType = 'Original'),
  avatars_tibble_knn15 %>% mutate(DataType = 'Synthetic')
)

# Perform PCA on combined data
combined_data_normalized <- scale(combined_data[, -which(names(combined_data) %in% c("DataType", "id"))])
combined_pca <- prcomp(combined_data_normalized, scale. = FALSE)

# Extract the first two principal components
combined_pca_data <- data.frame(combined_pca$x[, 1:2])
combined_pca_data$DataType <- combined_data$DataType

# Plot PCA with color differentiation
ggplot(combined_pca_data, aes(x = PC1, y = PC2, color = DataType)) +
  geom_point(alpha = 1) +
  theme_minimal() +
  labs(title = "PCA Plot", x = "Principal Component 1", y = "Principal Component 2", color = "Data Type")

## Compare the data

### summary distribution

summary(avatars_tibble_knn15)
summary(original_pk_wide_dapto)

### plot correlation

##Correlation Analysis
cor_real <- cor(original_pk_wide_dapto, use = "complete.obs")
cor_synthetic <- cor(avatars_tibble_knn15, use = "complete.obs")

# plots
ggcorrplot(cor_real, hc.order = TRUE, type = "lower",
           lab = TRUE,  pch.cex = 5,
           tl.cex = 6, lab_size = 2)
# plots
ggcorrplot(cor_synthetic, hc.order = TRUE, type = "lower",
           lab = TRUE,  pch.cex = 5,
           tl.cex = 6, lab_size = 2)

# ggpair compairons
group_comparison <-
  original_pk_wide_dapto %>% 
  mutate(group="original") %>% 
  bind_rows(avatars_tibble_knn15 %>% 
              mutate(group="synthetic")) 

pm <- group_comparison %>% ggpairs(
  ggplot2::aes(colour = group,alpha = 0.5),
  upper = list(continuous = wrap("cor", size = 1.5)),
  lower=list(combo=wrap("facethist", binwidth=0.5))) + 
  theme(strip.text.x = element_text(size = 5),
        strip.text.y = element_text(size = 5),axis.text = element_text(size = 5))
pm


#Process and plot the data

avatar_knn15 <- avatars_tibble_knn15 %>%
  mutate(ID = 1:nrow(avatars_tibble_knn15)) %>% 
  pivot_longer(cols = starts_with("conc_"), names_to = "time", values_to = "DV", names_prefix = "conc_") %>%
  mutate(time = as.numeric(time)) %>% 
  arrange(ID, time) %>% 
  relocate(ID) %>% 
  mutate(ii = 24, addl = 5, cmt = 1, ss=1, tinf = 0.5, evid = 0)

avatar_knn15 <- avatar_knn15 %>% 
  bind_rows(avatar_knn15 %>% group_by(ID) %>% slice_head(n=1) %>% mutate(time=0, evid=1))%>% 
  arrange(ID, time) 
write_csv(avatar_knn15, file = "avatar_simul_dapto_knn15.csv")

#Plot the data

simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(avatar_knn15 %>% mutate(group="synthetic")) %>% 
  ggplot(aes(x = time, y = DV, color = group)) + geom_point() +
  geom_point(alpha = 0.2, size = 3) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw()


simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(avatar_knn15 %>% 
              mutate(group="synthetic")) %>% 
  mutate(ID = as.factor(ID)) %>% 
  ggplot(aes(x = time, y = DV, color = ID)) +  geom_line(show.legend = FALSE)+
  geom_point(alpha = 0.2, size = 3, show.legend = FALSE) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw() + facet_wrap(~group)

## data augmentaiotn knn = 5

## data augmentation knn5
data_augment_avatar <- function(x) {
  data_normalized <- scale(original_pk_wide_dapto)
  pca <- prcomp(data_normalized, scale. = FALSE)# pour selecitonner le nombre de cp rank. = 3
  # Number of neighbors
  k <- 5  # Adjust this based on your requirement
  pca_transformed_data <- pca$x
  knn_result <- get.knn(pca_transformed_data, k)
  
  generate_avatar_weights <- function(knn_result, pca_transformed_data, k) {
    n <- nrow(pca_transformed_data)
    avatar_weights <- matrix(nrow = n, ncol = k)
    
    for (i in 1:n) {
      # Step 1: Inverse of Distances
      distances <- sqrt(rowSums((pca_transformed_data[knn_result$nn.index[i, ], ] - pca_transformed_data[i, ])^2))
      inverse_distances <- 1 / distances
      
      # Step 2: Random Weights
      set.seed( str_c(1,x))
      random_weights <- rexp(k, rate = 1)
      
      # Step 3: Contribution Factors
      set.seed( str_c(1,x))
      shuffled_indices <- sample(k)
      contribution_factors <- 1 / (2^shuffled_indices)
      
      # Step 4: Calculate Weights
      weights <- inverse_distances * random_weights * contribution_factors
      
      # Step 5: Normalize Weights
      normalized_weights <- weights / sum(weights)
      
      avatar_weights[i, ] <- normalized_weights
    }
    
    return(avatar_weights)
  }
  
  
  
  # Generate avatar weights
  avatar_weights <- generate_avatar_weights(knn_result, pca_transformed_data, k)
  
  # Assuming pca_result, avatar_weights, knn_result$nn.index, and pca_transformed_data are already defined
  
  # Function to generate avatars in PCA space based on weights
  generate_avatars_pca_space <- function(pca_transformed_data, knn_indices, weights) {
    n <- nrow(pca_transformed_data)
    avatars_pca <- matrix(nrow = n, ncol = ncol(pca_transformed_data))
    
    for (i in 1:n) {
      weighted_avatars <- pca_transformed_data[knn_indices[i, ], ] * weights[i, ]
      avatars_pca[i, ] <- colSums(weighted_avatars)
    }
    
    return(avatars_pca)
  }
  # Generate avatars in PCA space
  avatars_pca_space <- generate_avatars_pca_space(pca_transformed_data, knn_result$nn.index, avatar_weights)
  # Assuming 'aids_pca' is the PCA object and 'avatars_pca_space' contains the avatars in PCA space
  # Inverse PCA transformation
  inverse_pca <- function(pca_object, pca_data) {
    return(pca_data %*% t(pca_object$rotation) + matrix(pca_object$center, nrow = nrow(pca_data), ncol = ncol(pca_object$rotation), byrow = TRUE))
  }
  avatars_original_scale <- inverse_pca(pca, avatars_pca_space)
  
  # Assuming 'aids_data_normalized' contains the scaling attributes of the original data
  # Inverse normalization (if the original data was normalized)
  avatars_rescaled <- scale(avatars_original_scale, center = FALSE, scale = 1/attr(data_normalized, "scaled:scale"))
  avatars_rescaled <- sweep(avatars_rescaled, 2, attr(data_normalized, "scaled:center"), "+")
  
  avatars_tibble <- as_tibble(avatars_rescaled) %>% 
    mutate(
      SEX   = if_else(SEX>0.7 ,1,0))
}

iteration <- c(1:4)

augmented_data_5 <- map_dfr(iteration, data_augment_avatar, .id = "iter_")


## Compare the data

### summary dsitribution
summary(augmented_data_5)
summary(original_pk_wide_dapto)


##Correlation Analysis
cor_real <- cor(original_pk_wide_dapto, use = "complete.obs")
cor_synthetic <- cor(augmented_data_5 %>% select(-iter_), use = "complete.obs")

# plots
ggcorrplot(cor_real, hc.order = TRUE, type = "lower",
           lab = TRUE,  pch.cex = 5,
           tl.cex = 6, lab_size = 2)
# plots
ggcorrplot(cor_synthetic, hc.order = TRUE, type = "lower",
           lab = TRUE,  pch.cex = 5,
           tl.cex = 6, lab_size = 2)


#Graphical exploration of distribution


group_comparison <-
  original_pk_wide_dapto %>% 
  mutate(group="original") %>% 
  bind_rows(augmented_data_5 %>% 
              mutate(group="synthetic")) %>% 
  select(-iter_)

pm <- group_comparison %>% ggpairs(
  ggplot2::aes(colour = group,alpha = 0.5),
  upper = list(continuous = wrap("cor", size = 1.5)),
  lower=list(combo=wrap("facethist", binwidth=0.5))) + 
  theme(strip.text.x = element_text(size = 5),
        strip.text.y = element_text(size = 5),axis.text = element_text(size = 5))
pm
#ggsave("comparaison_distribution_PK_augmented_knn5.pdf")


#Process and plot the data

avatar_augmented_knn5 <- augmented_data_5 %>%
  mutate(ID = 1:nrow(augmented_data_5)) %>% 
  pivot_longer(cols = starts_with("conc_"), names_to = "time", values_to = "DV", names_prefix = "conc_") %>%
  mutate(time = as.numeric(time)) %>% 
  arrange(ID, time) %>% 
  relocate(ID) %>% 
  mutate(ii = 24, addl = 5, cmt = 1, ss=1, tinf = 0.5, evid=0)

avatar_augmented_knn5 <- avatar_augmented_knn5 %>% 
  bind_rows(avatar_knn15 %>% group_by(ID) %>% slice_head(n=1) %>% mutate(time=0, evid=1))%>% 
  arrange(ID, time) 
write_csv(avatar_augmented_knn5, file = "avatar_simul_dapto_augmented_knn5.csv")


#Plots

## plots
simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(avatar_augmented_knn5 %>% mutate(group="synthetic")) %>% 
  ggplot(aes(x = time, y = DV, color = group)) + 
  geom_point(alpha = 0.2, size = 3) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw()

simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(avatar_augmented_knn5 %>% 
              mutate(group="synthetic")) %>% 
  mutate(ID = as.factor(ID)) %>% 
  ggplot(aes(x = time, y = DV, color = ID)) +  geom_line(show.legend = FALSE)+
  geom_point(alpha = 0.2, size = 3, show.legend = FALSE) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw() + facet_wrap(~group)


## data augmentaiotn knn = 15

data_augment_avatar <- function(x) {
  data_normalized <- scale(original_pk_wide_dapto)
  pca <- prcomp(data_normalized, scale. = FALSE)# pour selecitonner le nombre de cp rank. = 3
  # Number of neighbors
  k <- 15  # Adjust this based on your requirement
  pca_transformed_data <- pca$x
  knn_result <- get.knn(pca_transformed_data, k)
  
  generate_avatar_weights <- function(knn_result, pca_transformed_data, k) {
    n <- nrow(pca_transformed_data)
    avatar_weights <- matrix(nrow = n, ncol = k)
    
    for (i in 1:n) {
      # Step 1: Inverse of Distances
      distances <- sqrt(rowSums((pca_transformed_data[knn_result$nn.index[i, ], ] - pca_transformed_data[i, ])^2))
      inverse_distances <- 1 / distances
      
      # Step 2: Random Weights
      set.seed( str_c(1,x))
      random_weights <- rexp(k, rate = 1)
      
      # Step 3: Contribution Factors
      set.seed( str_c(1,x))
      shuffled_indices <- sample(k)
      contribution_factors <- 1 / (2^shuffled_indices)
      
      # Step 4: Calculate Weights
      weights <- inverse_distances * random_weights * contribution_factors
      
      # Step 5: Normalize Weights
      normalized_weights <- weights / sum(weights)
      
      avatar_weights[i, ] <- normalized_weights
    }
    
    return(avatar_weights)
  }
  
  
  
  # Generate avatar weights
  avatar_weights <- generate_avatar_weights(knn_result, pca_transformed_data, k)
  
  # Assuming pca_result, avatar_weights, knn_result$nn.index, and pca_transformed_data are already defined
  
  # Function to generate avatars in PCA space based on weights
  generate_avatars_pca_space <- function(pca_transformed_data, knn_indices, weights) {
    n <- nrow(pca_transformed_data)
    avatars_pca <- matrix(nrow = n, ncol = ncol(pca_transformed_data))
    
    for (i in 1:n) {
      weighted_avatars <- pca_transformed_data[knn_indices[i, ], ] * weights[i, ]
      avatars_pca[i, ] <- colSums(weighted_avatars)
    }
    
    return(avatars_pca)
  }
  # Generate avatars in PCA space
  avatars_pca_space <- generate_avatars_pca_space(pca_transformed_data, knn_result$nn.index, avatar_weights)
  # Assuming 'aids_pca' is the PCA object and 'avatars_pca_space' contains the avatars in PCA space
  # Inverse PCA transformation
  inverse_pca <- function(pca_object, pca_data) {
    return(pca_data %*% t(pca_object$rotation) + matrix(pca_object$center, nrow = nrow(pca_data), ncol = ncol(pca_object$rotation), byrow = TRUE))
  }
  avatars_original_scale <- inverse_pca(pca, avatars_pca_space)
  
  # Assuming 'aids_data_normalized' contains the scaling attributes of the original data
  # Inverse normalization (if the original data was normalized)
  avatars_rescaled <- scale(avatars_original_scale, center = FALSE, scale = 1/attr(data_normalized, "scaled:scale"))
  avatars_rescaled <- sweep(avatars_rescaled, 2, attr(data_normalized, "scaled:center"), "+")
  
  avatars_tibble <- as_tibble(avatars_rescaled) %>% 
    mutate(
      SEX   = if_else(SEX>0.7 ,1,0))
}

iteration <- c(1:4)

augmented_data_15 <- map_dfr(iteration, data_augment_avatar, .id = "iter_")



## Compare the data

### summary dsitribution
summary(augmented_data_15)
summary(original_pk_wide_dapto)


##Correlation Analysis
cor_real <- cor(original_pk_wide_dapto, use = "complete.obs")
cor_synthetic <- cor(augmented_data_15 %>% select(-iter_), use = "complete.obs")

# plots
ggcorrplot(cor_real, hc.order = TRUE, type = "lower",
           lab = TRUE,  pch.cex = 5,
           tl.cex = 6, lab_size = 2)
# plots
ggcorrplot(cor_synthetic, hc.order = TRUE, type = "lower",
           lab = TRUE,  pch.cex = 5,
           tl.cex = 6, lab_size = 2)


#Graphical exploration of distribution


group_comparison <-
  original_pk_wide_dapto %>% 
  mutate(group="original") %>% 
  bind_rows(augmented_data_15 %>% 
              mutate(group="synthetic")) %>% 
  select(-iter_)

pm <- group_comparison %>% ggpairs(
  ggplot2::aes(colour = group,alpha = 0.5),
  upper = list(continuous = wrap("cor", size = 1.5)),
  lower=list(combo=wrap("facethist", binwidth=0.5))) + 
  theme(strip.text.x = element_text(size = 5),
        strip.text.y = element_text(size = 5),axis.text = element_text(size = 5))
pm
#ggsave("comparaison_distribution_PK_augmented_knn5.pdf")



#Process and plot the data


avatar_augmented_knn15 <- augmented_data_15 %>%
  mutate(ID = 1:nrow(augmented_data_15)) %>% 
  pivot_longer(cols = starts_with("conc_"), names_to = "time", values_to = "DV", names_prefix = "conc_") %>%
  mutate(time = as.numeric(time)) %>% 
  arrange(ID, time) %>% 
  relocate(ID) %>% 
  mutate(ii = 24, addl = 5, cmt = 1, ss=1, tinf = 0.5, evid = 0)

avatar_augmented_knn15 <- avatar_augmented_knn15 %>% 
  bind_rows(avatar_knn15 %>% group_by(ID) %>% slice_head(n=1) %>% mutate(time=0, evid=1))%>% 
  arrange(ID, time) 
write_csv(avatar_augmented_knn15, file = "avatar_simul_dapto_augmented_knn15.csv")


#Plots

## plots
simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(avatar_augmented_knn15 %>% mutate(group="synthetic")) %>% 
  ggplot(aes(x = time, y = DV, color = group)) + 
  geom_point(alpha = 0.2, size = 3) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw()

simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(avatar_augmented_knn15 %>% 
              mutate(group="synthetic")) %>% 
  mutate(ID = as.factor(ID)) %>% 
  ggplot(aes(x = time, y = DV, color = ID)) +  geom_line(show.legend = FALSE)+
  geom_point(alpha = 0.2, size = 3, show.legend = FALSE) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw() + facet_wrap(~group)


# ouverture et plot simul dapto ctgan clement

#ctgan
library(readr)
ctgan_dapto <- read_csv("~/Downloads/dapto_ctgan_data_3_multi000.dat")
#Process and plot the data


ctgan <- ctgan_dapto %>%
  mutate(ID = 1:nrow(ctgan_dapto)) %>% 
  pivot_longer(cols = starts_with("conc_"), names_to = "time", values_to = "DV", names_prefix = "conc_") %>%
  mutate(time = as.numeric(time)) %>% 
  arrange(ID, time) %>% 
  relocate(ID) %>% 
  mutate(ii = 24, addl = 5, cmt = 1, ss=1, tinf = 0.5, evid = 0) 

ctgan <- ctgan %>% 
  bind_rows(ctgan %>% group_by(ID) %>% slice_head(n=1) %>% mutate(time=0, evid=1))%>% 
  arrange(ID, time) 
write_csv(ctgan, file = "ctgan_simul_dapto.csv")

#Plot the data

simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(ctgan %>% mutate(group="synthetic")) %>% 
  ggplot(aes(x = time, y = DV, color = group)) + geom_point() +
  geom_point(alpha = 0.2, size = 3) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw()


simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(ctgan %>% 
              mutate(group="synthetic")) %>% 
  mutate(ID = as.factor(ID)) %>% 
  ggplot(aes(x = time, y = DV, color = ID)) +  geom_line(show.legend = FALSE)+
  geom_point(alpha = 0.2, size = 3, show.legend = FALSE) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw() + facet_wrap(~group)

#ctgan large
#Process and plot the data
library(readr)
ctgan_dapto_large <- read_csv("~/Downloads/dapto_ctgan_data_large3_multi_000.dat")

ctgan_large <- ctgan_dapto_large %>%
  mutate(ID = 1:nrow(ctgan_dapto_large)) %>% 
  pivot_longer(cols = starts_with("conc_"), names_to = "time", values_to = "DV", names_prefix = "conc_") %>%
  mutate(time = as.numeric(time)) %>% 
  arrange(ID, time) %>% 
  relocate(ID) %>% 
  mutate(ii = 24, addl = 5, cmt = 1, ss=1, tinf = 0.5, evid = 0) 

ctgan_large <- ctgan_large %>% 
  bind_rows(ctgan_large %>% group_by(ID) %>% slice_head(n=1) %>% mutate(time=0, evid=1))%>% 
  arrange(ID, time) 
write_csv(ctgan_large, file = "ctgan_large_simul_dapto.csv")

#Plot the data

simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(ctgan_large %>% mutate(group="synthetic")) %>% 
  ggplot(aes(x = time, y = DV, color = group)) + geom_point() +
  geom_point(alpha = 0.2, size = 3) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw()


simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(ctgan_large %>% 
              mutate(group="synthetic")) %>% 
  mutate(ID = as.factor(ID)) %>% 
  ggplot(aes(x = time, y = DV, color = ID)) +  geom_line(show.legend = FALSE)+
  geom_point(alpha = 0.2, size = 3, show.legend = FALSE) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw() + facet_wrap(~group)


#tvae

library(readr)
tvae_dapto<- read_csv("~/Downloads/dapto_tvae_data3_multi_000.dat")
#Process and plot the data


tvae <- tvae_dapto %>%
  mutate(ID = 1:nrow(tvae_dapto)) %>% 
  pivot_longer(cols = starts_with("conc_"), names_to = "time", values_to = "DV", names_prefix = "conc_") %>%
  mutate(time = as.numeric(time)) %>% 
  arrange(ID, time) %>% 
  relocate(ID) %>% 
  mutate(ii = 24, addl = 5, cmt = 1, ss=1, tinf = 0.5, evid = 0) 

tvae <- tvae %>% 
  bind_rows(tvae %>% group_by(ID) %>% slice_head(n=1) %>% mutate(time=0, evid=1))%>% 
  arrange(ID, time) 
write_csv(tvae, file = "tvae_simul_dapto.csv")

#Plot the data

simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(tvae %>% mutate(group="synthetic")) %>% 
  ggplot(aes(x = time, y = DV, color = group)) + geom_point() +
  geom_point(alpha = 0.2, size = 3) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw()


simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(tvae %>% 
              mutate(group="synthetic")) %>% 
  mutate(ID = as.factor(ID)) %>% 
  ggplot(aes(x = time, y = DV, color = ID)) +  geom_line(show.legend = FALSE)+
  geom_point(alpha = 0.2, size = 3, show.legend = FALSE) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw() + facet_wrap(~group)

#tvae large

library(readr)
tvae_dapto_large <- read_csv("~/Downloads/dapto_tvae_data_large3_multi_000.dat")
#Process and plot the data


tvae_large <- tvae_dapto_large %>%
  mutate(ID = 1:nrow(tvae_dapto_large)) %>% 
  pivot_longer(cols = starts_with("conc_"), names_to = "time", values_to = "DV", names_prefix = "conc_") %>%
  mutate(time = as.numeric(time)) %>% 
  arrange(ID, time) %>% 
  relocate(ID) %>% 
  mutate(ii = 24, addl = 5, cmt = 1, ss=1, tinf = 0.5, evid = 0) 

tvae_large <- tvae_large %>% 
  bind_rows(tvae_large %>% group_by(ID) %>% slice_head(n=1) %>% mutate(time=0, evid=1))%>% 
  arrange(ID, time) 
write_csv(tvae_large, file = "tvae_large_simul_dapto.csv")

#Plot the data

simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(tvae_large %>% mutate(group="synthetic")) %>% 
  ggplot(aes(x = time, y = DV, color = group)) + geom_point() +
  geom_point(alpha = 0.2, size = 3) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw()


simul_datpo_all %>% mutate(group="original") %>% 
  bind_rows(tvae_large %>% 
              mutate(group="synthetic")) %>% 
  mutate(ID = as.factor(ID)) %>% 
  ggplot(aes(x = time, y = DV, color = ID)) +  geom_line(show.legend = FALSE)+
  geom_point(alpha = 0.2, size = 3, show.legend = FALSE) +
  labs(x = "Time (h)", y = "daptomycin concentration(mg/L)")+
  theme_bw() + facet_wrap(~group)
