# Libraries
library(sf)
library(tidyverse)
library(parallel)
library(mapview)
library(spdep)
library(spatialreg)
library(tmap)
library(viridisLite)
library(sandwich)
library(clubSandwich)
library(multiwayvcov)
library(broom)
library(lmtest)

if (Sys.getenv("USER")=="ckoenig") {
  data_raw <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/data_raw/"
  workdata <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/workdata/"
  tabs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/tables/"
  figs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/figures/"
}

## load pre-processed data (analytical sample)
data <- readRDS(paste0(workdata,"data_core_inhabited_02.rds"))%>%
  mutate(rownum=row_number())


## filter on grids for which spatial lags of neighborhood composition as well as environmental variables
## can be generated (contiguity criterion)
listw <- poly2nb(data, queen = TRUE, snap = 10)
list_no_neighbors <- which(card(listw)==0) %>% as.list()

data <- data %>%
  filter(!(rownum %in% list_no_neighbors))

listw <- poly2nb(data, queen = TRUE, snap = 10)
listw <- nb2listw(listw,zero.policy = TRUE, style = "W")

## generate spatial lags of key variables
aux <- data %>%
  mutate(across(.cols=c(share_sgb_total,
                        share_workingage_foreign,
                        share_age_0_14,
                        share_age_15_65,
                        share_age_over_65,
                        airpoll_w,
                        airpoll_w_std,
                        greens18,
                        greens18_std,
                        lack_greens18,
                        lack_greens18_std,
                        noise_exposure,
                        noise_exposure_std,
                        meb,
                        meb_std,
                        grid_pop),
                .fns= ~lag.listw(x=listw,var=.x),
                .names = "{.col}_lag")) %>%
  st_drop_geometry() %>%
  select(id,
         share_sgb_total_lag,
         share_workingage_foreign_lag,
         share_age_0_14_lag,
         share_age_15_65_lag,
         share_age_over_65_lag,
         airpoll_w_lag,
         airpoll_w_std_lag,
         greens18_lag,
         greens18_std_lag,
         lack_greens18_lag,
         lack_greens18_std_lag,
         noise_exposure_lag,
         noise_exposure_std_lag,
         meb_lag,
         meb_std_lag,
         grid_pop_lag)

class(aux) <- "data.frame"

## merge data and spatial lags
data <- data %>%
  inner_join(aux,by="id") %>%
  ## standardization of key variables
  mutate(across(.cols = c(share_sgb_total,share_workingage_foreign,
                          share_sgb_total_lag,share_workingage_foreign_lag,
                          share_age_0_14,share_age_15_65,share_age_over_65,grid_pop,
                          share_age_0_14_lag,share_age_15_65_lag,share_age_over_65_lag,grid_pop_lag),
                .fns = ~ .x/sqrt(Hmisc::wtd.var(.x,weights = data$grid_pop)),
                .names = "{.col}_std")) %>%
  rename(ags8_gen=municipality,
         ags8_key=AGS)


###################################################
################ ANALYSES FUA CORES ###############
###################################################
outcomes <- c(rep("share_sgb_total_std",6),rep("share_workingage_foreign_std",6)) %>% as.list()
treatments <- rep(c("airpoll_w_std","greens18_std","lack_greens18_std","noise_exposure_std","meb","meb_std"),2) %>% as.list()
treatment_lags <- rep(c("airpoll_w_std_lag","greens18_std_lag","lack_greens18_std_lag","noise_exposure_std_lag","meb_lag","meb_std_lag"),2) %>% as.list()

outcomes_treatments_list <- Map(c,outcomes,treatments,treatment_lags)

x <- outcomes_treatments_list[[4]] 

system.time(regression_list <- mclapply(outcomes_treatments_list, function(x){
  outcome <- x[[1]]
  treatment <- x[[2]]
  treatment_lag <- x[[3]]
  
  ## filter to complete cases along main outcomes, main EQ treatments (except noise exposure) and main controls
  df <- data %>% 
    filter(complete.cases(share_age_0_14_std,share_age_15_65_std,share_age_over_65_std,mun_pop,grid_pop,
                          eval(parse(text = outcome)),eval(parse(text = treatment)))) 
  
  ##### SPATIAL REGRESSION MODELS #####
  ## SLX
  formula_m1_slx <- as.formula(paste0(outcome,"~",treatment,"+",treatment_lag))

  ## unweighted
  m1_slx_uw <- lm(formula = formula_m1_slx,
                  data = data)
  
  coef_m1_slx_uw <- tidy(m1_slx_uw, conf.int = TRUE) %>% 
    mutate(model="m1_slx_uw", outcome=outcome, treatment=treatment)%>%
    filter(term==treatment)
  
  ## SAR
  formula_m1_sar <- as.formula(paste0(outcome,"~",treatment))
  m1_sar <- lagsarlm(formula =  formula_m1_sar,
                     data = data,
                     listw = listw,
                     Durbin = FALSE)
  
  coef_m1_sar <- tidy(m1_sar, conf.int = TRUE) %>% 
    mutate(model="m1_sar", outcome=outcome, treatment=treatment)%>%
    filter(term==treatment)
  
  ## SDM
  formula_m1_sdm <- as.formula(paste0(outcome,"~",treatment))
  m1_sdm <- lagsarlm(formula =  formula_m1_sdm,
                     data = data,
                     listw = listw,
                     Durbin = TRUE)
  
  coef_m1_sdm <- tidy(m1_sdm, conf.int = TRUE) %>% 
    mutate(model="m1_sdm", outcome=outcome, treatment=treatment)%>%
    filter(term==treatment)
  
  ## SEM
  formula_m1_sem <- as.formula(paste0(outcome,"~",treatment))
  m1_sem <- errorsarlm(formula =  formula_m1_sem,
                       data = data,
                       listw = listw,
                       Durbin = FALSE)
  
  coef_m1_sem <- tidy(m1_sem, conf.int = TRUE) %>% 
    mutate(model="m1_sem", outcome=outcome, treatment=treatment)%>%
    filter(term==treatment)
  
  ## SDEM
  formula_m1_sdem <- as.formula(paste0(outcome,"~",treatment))
  m1_sdem <- errorsarlm(formula =  formula_m1_sdem,
                        data = data,
                        listw = listw,
                        Durbin = TRUE)
  
  coef_m1_sdem <- tidy(m1_sdem, conf.int = TRUE) %>% 
    mutate(model="m1_sdem", outcome=outcome, treatment=treatment)%>%
    filter(term==treatment)
  
  modelsummary::modelsummary(list(m1_slx_uw,
                                  m1_sar,
                                  m1_sdm,
                                  m1_sem,
                                  m1_sdem),
                             stars = TRUE,
                             output = paste0(tabs,"coefs_spatregs_",outcome,"_",treatment,".md"))
  
  
  coefs_spatialregs <- rbind(coef_m1_slx_uw,
                             coef_m1_sar,
                             coef_m1_sdm,
                             coef_m1_sem,
                             coef_m1_sdem)
  
  
  coefs_spatialregs
  
},mc.cores = 20))

coefs_spatialregs <- do.call(what=rbind, args=regression_list)
#saveRDS(coefs_spatialregs,paste0(tabs,"coefs_spatialregs_all.rds"))

tmp <- coefs_spatialregs %>%
  filter(treatment %in% c("airpoll_w_std","greens18_std","noise_exposure_std","meb_std")) %>%
  mutate(treatment_viz= treatment %>% fct_recode("Air pollution (std.)"="airpoll_w_std",
                                                 "Green space (std.)"="greens18_std",
                                                 "Noise (std.)"="noise_exposure_std",
                                                 "Multiple environmental\nburdens (std.)"="meb_std") %>%
           fct_relevel("Air pollution (std.)","Noise (std.)","Green space (std.)","Multiple environmental\nburdens (std.)"),
         outcome_viz=outcome %>% fct_recode("Poverty rate (std.)"="share_sgb_total_std",
                                            "Foreign minorities (std.)"="share_workingage_foreign_std") %>%
           fct_relevel("Poverty rate (std.)","Foreign minorities (std.)"),
         model=model %>% 
           fct_recode("SAR"="m1_sar",
                      "SEM" = "m1_sem",
                      "SLX"="m1_slx_uw",
                      "SDM"="m1_sdm",
                      "SDEM"="m1_sdem") %>%
           fct_relevel("SDEM","SDM","SLX","SEM","SAR"))



p <- ggplot(data=tmp %>% filter(model %in% c("SAR","SEM","SLX","SDM","SDEM")),
            aes(y=treatment_viz,x=estimate,color=model))+
  geom_pointrange(aes(
    xmin=conf.low,
    xmax=conf.high),
    position=position_dodge2(width = .6)) +
  geom_vline(xintercept = 0)+
  scale_y_discrete(limits=rev)+
  scale_color_viridis_d(limits=rev)+
  labs(title = "",
       y="Environmental quality variable",
       x="Coefficient estimate (w/ 95% CI)",
       color="Spatial regression model")+
  facet_grid(~outcome_viz)+
  coord_cartesian(xlim = c(-.15,.5))+
  theme_bw()+
  theme(text = element_text(face = "bold"))

ggsave(p,filename=paste0(figs,"coefs_spatialregs.svg"),
       width = 10, height = 5)
