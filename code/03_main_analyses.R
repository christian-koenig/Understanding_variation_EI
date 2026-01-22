# Libraries
library(tidyverse)
library(sf)
library(tmap)
library(parallel)
library(estimatr)
library(stats)
library(flextable)
library(officer)
library(lmtest)
library(sandwich)
library(clubSandwich)
library(miceadds)
library(multiwayvcov)
library(fixest)
library(colorspace)
library(writexl)
library(spdep)
library(spatialreg)
library(sp)
library(modelsummary)

## DEFINE DIRECTORIES
## compute/kratos
if (Sys.getenv("USER")=="ckoenig") {
  data_raw <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/data_raw/"
  workdata <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/workdata/"
  tabs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/tables/"
  figs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/figures/"
}


## read pre-processed data on inhabited as well as all grid cells in cities
data <- readRDS(paste0(workdata,"data_core_inhabited_02.rds"))
data_core_full <- readRDS(paste0(workdata,"data_core_full_02.rds"))


## generation of spatial lags of EQ variables;
## based on all urban grid (incl. those not inhabited by any residents)
listw <- poly2nb(data_core_full, queen = TRUE)
listw <- nb2listw(listw,zero.policy = TRUE)

aux1 <- data_core_full %>%
  mutate(across(.cols=c(airpoll_w,
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

class(aux1) <- "data.frame"

## generation of spatial lags of socioeconomic and demgraphic variables;
## excl. those not inhabited by any residents to prevent non-finite values
listw2 <- poly2nb(data,queen = TRUE)
listw2 <- nb2listw(listw2,zero.policy = TRUE)

aux2 <- data %>%
  mutate(across(.cols=c(share_sgb_total,
                        share_workingage_foreign,
                        share_age_0_14,
                        share_age_15_65,
                        share_age_over_65),
                .fns= ~lag.listw(x=listw2,var=.x),
                .names = "{.col}_lag")) %>%
  st_drop_geometry() %>%
  select(id,
         share_sgb_total_lag,
         share_workingage_foreign_lag,
         share_age_0_14_lag,
         share_age_15_65_lag,
         share_age_over_65_lag)

class(aux2) <- "data.frame"

## merge data and spatial lags
data <- data %>%
  inner_join(aux1,by="id") %>%
  inner_join(aux2,by="id") %>%
  ## z-standardization of key variables
  mutate(across(.cols = c(share_sgb_total,share_workingage_foreign,
                          share_sgb_total_lag,share_workingage_foreign_lag,
                          share_age_0_14,share_age_15_65,share_age_over_65,grid_pop,
                          share_age_0_14_lag,share_age_15_65_lag,share_age_over_65_lag,grid_pop_lag),
                .fns = ~ .x/sqrt(Hmisc::wtd.var(.x,weights = data$grid_pop)),
                .names = "{.col}_std")) %>%
  rename(ags8_gen=municipality,
         ags8_key=AGS)



## spatially merge info on FUACODE and FUANAME by largest areal overlap
fua_core <- readRDS(paste0(workdata,"fua_core_shape_adjusted.rds"))
st_crs(fua_core)

system.time(data <- data %>% st_join(fua_core,largest = TRUE))


table(data$ags8_gen) %>% names()

######################################
##### MORAN'S I OF KEY VARIABLES #####
######################################

a <- moran.test(data$share_sgb_total,listw = listw2, zero.policy = TRUE) %>% 
  tidy() %>%
  mutate(variable="Share_SGB_total")

b <- moran.test(data$share_workingage_foreign,listw = listw2, zero.policy = TRUE)%>% 
  tidy() %>%
  mutate(variable="Share_workingage_foreign")

c <- moran.test(data$airpoll_w,listw = listw2, zero.policy = TRUE)%>% 
  tidy() %>%
  mutate(variable="Air_pollution")

d <- moran.test(data$greens18,listw = listw2, zero.policy = TRUE)%>% 
  tidy() %>%
  mutate(variable="Green_space")

e <- moran.test(data$lack_greens18,listw = listw2, zero.policy = TRUE)%>% 
  tidy() %>%
  mutate(variable="Lack_green_space")

f <- moran.test(data$noise_exposure,listw = listw2, zero.policy = TRUE)%>% 
  tidy() %>%
  mutate(variable="Noise")

g <- moran.test(data$meb,listw = listw2, zero.policy = TRUE)%>% 
  tidy() %>%
  mutate(variable="Multiple_environmental_burdens")

moran_estimates <- rbind(a,b,c,d,f,g)%>%
  select(estimate1,statistic,p.value,variable) %>%
  mutate(variable=as.factor(variable) %>% fct_relevel("Share_SGB_total","Share_workingage_foreign",
                                                      "Air_pollution","Green_space","Lack_green_space","Noise","Multiple_environmental_burdens"))
rm(a,b,c,d,e,f)

## MORAN FIGURE
moran_graph <- ggplot(moran_estimates, aes(x=estimate1,y=variable))+
  geom_bar(stat = "identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Moran's I",y="")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,size = 12,face = "bold",hjust = .95,vjust=.2),
        text = element_text(size = 12,face = "bold"))

ggsave(moran_graph,filename = paste0(figs,"moran_main_vars.svg"), width = 4,height = 6)
ggsave(moran_graph,filename = paste0(figs,"moran_main_vars.png"), width = 4,height = 6)

###################################################
############### DESCRIPTIVES ######################
###################################################

# generate functions for minmax and weighted means to be used in datasummary()
minmax <- function(x) paste0("[", round(min(x),2), ", ", round(max(x),2), "]")
w_mean_gridpop <- function(x,weight_vec) Weighted.Desc.Stat::w.mean(x,mu=data$grid_pop)
w_sd_gridpop <- function(x,weight_vec) Weighted.Desc.Stat::w.sd(x, mu=data$grid_pop)

# calculate summary statistics
rows <- tibble::tribble(~term,~n,~mean,~sd,~minmax,~histogram,
                        'a) Grid level','','','','','')

attr(rows,'position') <- c(1)

grid_pop_stats <- datasummary(data = data,
                              add_rows = rows,
                              formula = Heading("Grid population (unweighted)")*grid_pop
                              ~ N + Heading("Mean")*mean + Heading("SD")*sd + Heading("[Min,Max]")*minmax+Heading("Distribution")*Histogram,
                              output = "data.frame")

grid_pop_stats_flex <- datasummary(data = data,
                                   add_rows = rows,
                                   formula = Heading("Grid population (unweighted)")*grid_pop
                                   ~ N + Heading("Mean")*mean + Heading("SD")*sd + Heading("[Min,Max]")*minmax+Heading("Distribution")*Histogram,
                                   output = "flextable")

rows <- tibble::tribble(~term,~n,~mean,~sd,~minmax,~histogram,
                        'b) City level','','','','','')

attr(rows,'position') <- c(8)


grid_stats <-datasummary(data = data,
                         add_rows = rows,
                         formula = Heading("Share of SGBII recipients")*share_sgb_total+
                           Heading("Share of immigrant-origin population")*share_workingage_foreign+
                           Heading("Air pollution, toxicity-weighted (PM25,NO2,SO2)")*airpoll_w+
                           Heading("Green space coverage (in %)")*greens18+
                           Heading("Lack of green space (in %)")*lack_greens18+
                           Heading("Noise exposure (road, rail, aviation)")*noise_exposure+
                           Heading("Multiple_exposure_burden")*meb+
                           Heading("City size (in 100K inhabitants)")*mun_pop+
                           Heading("Segregation, SGBII")*dissim_sgb+
                           Heading("Segregation, immigrant-origin population")*dissim_ausl+
                           Heading("Residential centrality, SGBII")*dichte_arm+
                           Heading("Residential centrality, immigrant-origin population")*dichte_ausl
                         ~ N + Heading("Mean")*w_mean_gridpop + Heading("SD")*w_sd_gridpop + Heading("[Min,Max]")*minmax+Heading("Distribution")*Histogram,
                         output = "data.frame")

grid_stats_flex <-datasummary(data = data,
                              add_rows = rows,
                              formula = Heading("Share of SGBII recipients")*share_sgb_total+
                                Heading("Share of immigrant-origin population")*share_workingage_foreign+
                                Heading("Air pollution, toxicity-weighted (PM25,NO2,SO2)")*airpoll_w+
                                Heading("Green space coverage (in %)")*greens18+
                                Heading("Lack of green space (in %)")*lack_greens18+
                                Heading("Noise exposure (road, rail, aviation)")*noise_exposure+
                                Heading("Multiple_exposure_burden")*meb+
                                Heading("City size (in 100K inhabitants)")*mun_pop+
                                Heading("Segregation, SGBII")*dissim_sgb+
                                Heading("Segregation, immigrant-origin population")*dissim_ausl+
                                Heading("Residential centrality, SGBII")*dichte_arm+
                                Heading("Residential centrality, immigrant-origin population")*dichte_ausl
                              ~ N + Heading("Mean")*w_mean_gridpop + Heading("SD")*w_sd_gridpop + Heading("[Min,Max]")*minmax+Heading("Distribution")*Histogram,
                              output = "flextable")


stats_full <- rbind(grid_pop_stats,grid_stats)

flextable::set_flextable_defaults(font.size = 8,
                                  text.align = "left",
                                  line_spacing = 1,
                                  table.layout = "autofit") 

sect_properties <- prop_section(page_size = page_size(orient = "landscape"))

grid_pop_stats_flex %>%
  flextable::save_as_docx(path = paste0(tabs,"summary_stats_gridpop.docx"),
                          pr_section = sect_properties)

grid_stats_flex %>%
  flextable::save_as_docx(path = paste0(tabs,"summary_stats.docx"),
                          pr_section = sect_properties)

############################################################
################ MAIN SLX ANALYSES FUA CORES ###############
############################################################
outcomes <- c(rep("share_sgb_total_std",5),rep("share_workingage_foreign_std",5)) %>% as.list()
treatments <- rep(c("airpoll_w_std","greens18_std","lack_greens18_std","noise_exposure_std","meb_std"),2) %>% as.list()
treatment_lags <- rep(c("airpoll_w_std_lag","greens18_std_lag","lack_greens18_std_lag","noise_exposure_std_lag","meb_std_lag"),2) %>% as.list()

outcomes_treatments_list <- Map(c,outcomes,treatments,treatment_lags)

x <- outcomes_treatments_list[[2]] 

regression_list <- lapply(outcomes_treatments_list, function(x){
  outcome <- x[[1]]
  treatment <- x[[2]]
  treatment_lag <- x[[3]]
  
  df <- data
  
  ##### SLX #####
  formula_m1_slx <- as.formula(paste0(outcome,"~",treatment,"+",treatment_lag))
  formula_m1_fe_slx <- as.formula(paste0(outcome,"~",treatment,"+",treatment_lag,"+ factor(ags8_key)"))
  
  m1_slx <- lm(data = df,
               formula = formula_m1_slx,
               weights = grid_pop)
  
  coefs_m1_slx <- coeftest(m1_slx, vcov. = cluster.vcov(m1_slx, cluster =~ags8_gen))
  coefs_m1_slx <- tidy(coefs_m1_slx, conf.int = TRUE) %>% 
    mutate(model="m1_slx", outcome=outcome, treatment=treatment)%>%
    filter(term==treatment)
  
  m1_fe_slx <- lm(data = df,
                  formula = formula_m1_fe_slx,
                  weights = grid_pop)
  
  coefs_m1_fe_slx <- coeftest(m1_fe_slx, vcov. = cluster.vcov(m1_fe_slx, cluster =~ags8_gen))
  coefs_m1_fe_slx <- tidy(coefs_m1_fe_slx, conf.int = TRUE) %>% 
    mutate(model="m1_fe_slx", outcome=outcome, treatment=treatment)%>%
    filter(term==treatment)
  
  modelsummary(list(m1_slx,m1_fe_slx),
               stars = TRUE,
               estimate = "{estimate} {stars} ({std.error})",
               statistic = NULL,
               vcov =~ags8_gen,
               gof_omit = "AIC|BIC|F|Lik",
               output = paste0(tabs,"regtable_SLX_",outcome,"_",treatment,".txt"))
  
  coefs_slx <- rbind(coefs_m1_slx,coefs_m1_fe_slx)
  coefs_slx
})

coefs_slx <- do.call(rbind,args = regression_list)

## rename/relevel factors for visualization
coefs_slx <- coefs_slx %>% 
  filter(treatment %in% c("airpoll_w_std","greens18_std","meb_std","noise_exposure_std")) %>%
  mutate(treatment_viz= treatment %>% fct_recode("Air pollution (std.)"="airpoll_w_std",
                                                 "Green space (std.)"="greens18_std",
                                                 "Noise (std.)"="noise_exposure_std",
                                                 "Multiple environmental\nburdens (std.)"="meb_std") %>%
           fct_relevel("Multiple environmental\nburdens (std.)","Green space (std.)","Noise (std.)","Air pollution (std.)"),
         outcome_viz=outcome %>% fct_recode("Poverty rate (std.)"="share_sgb_total_std",
                                            "Foreign minorities (std.)"="share_workingage_foreign_std") %>%
           fct_relevel("Poverty rate (std.)","Foreign minorities (std.)"),
         model=model %>% 
           fct_recode("bivariate SLX"="m1_slx", 
                      "bivariate SLX w/ city FEs"="m1_fe_slx") %>%
           fct_relevel("bivariate SLX w/ city FEs","bivariate SLX")) %>%
  rename(Model=model)


graph <- ggplot(coefs_slx, aes(x = estimate, y = treatment_viz, color = Model, shape = Model)) +
  geom_vline(xintercept = 0, color = "gray40", linetype = "dashed", size = 0.9) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.20, size = 0.8,
                 position = position_dodge(width  = 0.6)) +
  geom_point(size = 3.2, position = position_dodge(width = 0.6)) +
  facet_grid(~ outcome_viz) +
  scale_color_viridis_d(option = "plasma", begin = .7, end=0) +
  scale_shape_manual(values = c(15,19))+
  labs(
    x = "Coefficient estimate (95% CI)",
    y = "Environmental quality treatment",
    color = NULL, shape = NULL) +
  guides(color = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE))+
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),      # remove horizontal grid lines for clarity
    #panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"))



ggsave(graph, filename =paste0(figs,"slx_results_pooled_all.pdf"),
       height = 6, width = 9 )

ggsave(graph, filename =paste0(figs,"slx_results_pooled_all.svg"),
       height = 6, width = 9 )

ggsave(graph, filename =paste0(figs,"slx_results_pooled_all.png"),
       height = 6, width = 9)


###################################
##### CITY SUBSAMPLE ANALYSES #####
###################################

ags8_list <- data$ags8_gen %>% table() %>% names() %>% as.list()

x <- ags8_list[[54]] 

ags8_analyses_list <- lapply(ags8_list, function(x){
  
  df <- data %>% filter(ags8_gen==x)
  
  ## formulas sgb share slx
  formula_m1_slx_sgbtotal_airpoll_w <- as.formula("share_sgb_total_std ~ airpoll_w_std+airpoll_w_std_lag")
  formula_m1_slx_sgbtotal_greens18 <- as.formula("share_sgb_total_std ~ greens18_std+greens18_std_lag")
  formula_m1_slx_sgbtotal_lackgreens18 <- as.formula("share_sgb_total_std ~ lack_greens18_std+lack_greens18_std_lag")
  formula_m1_slx_sgbtotal_noise <- as.formula("share_sgb_total_std ~ noise_exposure_std+noise_exposure_std_lag")
  formula_m1_slx_sgbtotal_meb <- as.formula("share_sgb_total_std ~ meb_std+meb_std_lag")
  
  
  ## formulas foreigner share slx
  formula_m1_slx_foreign_airpoll_w <- as.formula("share_workingage_foreign_std ~ airpoll_w_std+airpoll_w_std_lag")
  formula_m1_slx_foreign_greens18 <- as.formula("share_workingage_foreign_std ~ greens18_std+greens18_std_lag")
  formula_m1_slx_foreign_lackgreens18 <- as.formula("share_workingage_foreign_std ~ lack_greens18_std+lack_greens18_std_lag")
  formula_m1_slx_foreign_noise <- as.formula("share_workingage_foreign_std ~ noise_exposure_std+noise_exposure_std_lag")
  formula_m1_slx_foreign_meb <- as.formula("share_workingage_foreign_std ~ meb_std+meb_std_lag")
  
  
  ## SLX
  m1_slx_sgbtotal_airpoll_w <- lm(data = df,
                                  formula = formula_m1_slx_sgbtotal_airpoll_w,
                                  weights = grid_pop)
  
  coefs_m1_slx_sgbtotal_airpoll_w <- coeftest(m1_slx_sgbtotal_airpoll_w, vcov. = vcovHC(m1_slx_sgbtotal_airpoll_w, type = "HC1"))
  coefs_m1_slx_sgbtotal_airpoll_w <- tidy(coefs_m1_slx_sgbtotal_airpoll_w, conf.int = TRUE) %>% 
    mutate(model="m1",modeltype="SLX", outcome="share_sgb_total_std", treatment="airpoll_w_std")%>%
    filter(str_detect(term,"airpoll_w"))
  
  
  m1_slx_sgbtotal_greens18 <- lm(data = df,
                                 formula = formula_m1_slx_sgbtotal_greens18,
                                 weights = grid_pop)
  
  coefs_m1_slx_sgbtotal_greens18 <- coeftest(m1_slx_sgbtotal_greens18, vcov. = vcovHC(m1_slx_sgbtotal_greens18, type = "HC1"))
  coefs_m1_slx_sgbtotal_greens18 <- tidy(coefs_m1_slx_sgbtotal_greens18, conf.int = TRUE) %>% 
    mutate(model="m1",modeltype="SLX", outcome="share_sgb_total_std", treatment="greens18_std")%>%
    filter(str_detect(term,"greens18"))
  
  m1_slx_sgbtotal_lackgreens18 <- lm(data = df,
                                     formula = formula_m1_slx_sgbtotal_lackgreens18,
                                     weights = grid_pop)
  
  coefs_m1_slx_sgbtotal_lackgreens18 <- coeftest(m1_slx_sgbtotal_lackgreens18, vcov. = vcovHC(m1_slx_sgbtotal_lackgreens18, type = "HC1"))
  coefs_m1_slx_sgbtotal_lackgreens18 <- tidy(coefs_m1_slx_sgbtotal_lackgreens18, conf.int = TRUE) %>% 
    mutate(model="m1",modeltype="SLX", outcome="share_sgb_total_std", treatment="lack_greens18_std")%>%
    filter(str_detect(term,"lack_greens18"))
  
  
  
  
  m1_slx_sgbtotal_noise <- lm(data = df,
                              formula = formula_m1_slx_sgbtotal_noise,
                              weights = grid_pop)
  
  coefs_m1_slx_sgbtotal_noise <- coeftest(m1_slx_sgbtotal_noise, vcov. = vcovHC(m1_slx_sgbtotal_noise, type = "HC1"))
  coefs_m1_slx_sgbtotal_noise <- tidy(coefs_m1_slx_sgbtotal_noise, conf.int = TRUE) %>% 
    mutate(model="m1",modeltype="SLX", outcome="share_sgb_total_std", treatment="noise_exposure_std")%>%
    filter(str_detect(term,"noise"))
  
  
  m1_slx_sgbtotal_meb <- lm(data = df,
                            formula = formula_m1_slx_sgbtotal_meb,
                            weights = grid_pop)
  
  coefs_m1_slx_sgbtotal_meb <- coeftest(m1_slx_sgbtotal_meb, vcov. = vcovHC(m1_slx_sgbtotal_meb, type = "HC1"))
  coefs_m1_slx_sgbtotal_meb <- tidy(coefs_m1_slx_sgbtotal_meb, conf.int = TRUE) %>% 
    mutate(model="m1",modeltype="SLX", outcome="share_sgb_total_std", treatment="meb_std")%>%
    filter(str_detect(term,"meb"))
  
  
  m1_slx_foreign_airpoll_w <- lm(data = df,
                                 formula = formula_m1_slx_foreign_airpoll_w,
                                 weights = grid_pop)
  
  coefs_m1_slx_foreign_airpoll_w <- coeftest(m1_slx_foreign_airpoll_w, vcov. = vcovHC(m1_slx_foreign_airpoll_w, type = "HC1"))
  coefs_m1_slx_foreign_airpoll_w <- tidy(coefs_m1_slx_foreign_airpoll_w, conf.int = TRUE) %>% 
    mutate(model="m1",modeltype="SLX", outcome="share_workingage_foreign_std", treatment="airpoll_w_std")%>%
    filter(str_detect(term,"airpoll_w"))
  
  
  m1_slx_foreign_greens18 <- lm(data = df,
                                formula = formula_m1_slx_foreign_greens18,
                                weights = grid_pop)
  
  coefs_m1_slx_foreign_greens18 <- coeftest(m1_slx_foreign_greens18, vcov. = vcovHC(m1_slx_foreign_greens18, type = "HC1"))
  coefs_m1_slx_foreign_greens18 <- tidy(coefs_m1_slx_foreign_greens18, conf.int = TRUE) %>% 
    mutate(model="m1",modeltype="SLX", outcome="share_workingage_foreign_std", treatment="greens18_std")%>%
    filter(str_detect(term,"greens18"))
  
  m1_slx_foreign_lackgreens18 <- lm(data = df,
                                    formula = formula_m1_slx_foreign_lackgreens18,
                                    weights = grid_pop)
  
  coefs_m1_slx_foreign_lackgreens18 <- coeftest(m1_slx_foreign_lackgreens18, vcov. = vcovHC(m1_slx_foreign_lackgreens18, type = "HC1"))
  coefs_m1_slx_foreign_lackgreens18 <- tidy(coefs_m1_slx_foreign_lackgreens18, conf.int = TRUE) %>% 
    mutate(model="m1",modeltype="SLX", outcome="share_workingage_foreign_std", treatment="lack_greens18_std")%>%
    filter(str_detect(term,"lack_greens18"))
  
  m1_slx_foreign_noise <- lm(data = df,
                             formula = formula_m1_slx_foreign_noise,
                             weights = grid_pop)
  
  coefs_m1_slx_foreign_noise <- coeftest(m1_slx_foreign_noise, vcov. = vcovHC(m1_slx_foreign_noise, type = "HC1"))
  coefs_m1_slx_foreign_noise <- tidy(coefs_m1_slx_foreign_noise, conf.int = TRUE) %>% 
    mutate(model="m1",modeltype="SLX", outcome="share_workingage_foreign_std", treatment="noise_exposure_std")%>%
    filter(str_detect(term,"noise"))
  
  
  m1_slx_foreign_meb <- lm(data = df,
                           formula = formula_m1_slx_foreign_meb,
                           weights = grid_pop)
  
  coefs_m1_slx_foreign_meb <- coeftest(m1_slx_foreign_meb, vcov. = vcovHC(m1_slx_foreign_meb, type = "HC1"))
  coefs_m1_slx_foreign_meb <- tidy(coefs_m1_slx_foreign_meb, conf.int = TRUE) %>% 
    mutate(model="m1",modeltype="SLX", outcome="share_workingage_foreign_std", treatment="meb_std")%>%
    filter(str_detect(term,"meb"))
  
  
  coefs_slx <- rbind(coefs_m1_slx_sgbtotal_airpoll_w,coefs_m1_slx_sgbtotal_noise,coefs_m1_slx_sgbtotal_greens18,coefs_m1_slx_sgbtotal_lackgreens18,coefs_m1_slx_sgbtotal_meb,
                     coefs_m1_slx_foreign_airpoll_w,coefs_m1_slx_foreign_noise,coefs_m1_slx_foreign_greens18,coefs_m1_slx_foreign_lackgreens18,coefs_m1_slx_foreign_meb) %>%
    dplyr:: mutate(ags8_gen=x)
  
  coefs_slx
})

ags8_sample <- readRDS(paste0(workdata,"ags8_Städte_noisedata.rds")) %>%
  select(AGS,GEN,geometry) %>%
  group_by(AGS,GEN) %>%
  summarise(geometry=st_union(geometry))

## combine city-level estimates of EI
coefs_slx <- do.call(rbind, args = ags8_analyses_list)%>%
  left_join(ags8_sample,by=c("ags8_gen"="GEN")) %>% 
  st_as_sf() %>% 
  st_set_crs(32632) %>%
  ungroup() %>%
  group_by(ags8_gen) %>%
  mutate(city_number=cur_group_id()) %>%
  ungroup() %>%
  mutate(ags8_gen_num=paste(ags8_gen,city_number,sep = "  ")) 

lvls <- table(coefs_slx$ags8_gen_num) %>%
  names() %>% sort(decreasing = F)

coefs_slx <- coefs_slx %>%
  filter(treatment %in% c("airpoll_w_std","greens18_std","lack_greens18_std","meb_std","noise_exposure_std")) %>%
  mutate(treatment_viz= treatment %>% fct_recode("Air pollution (std.)"="airpoll_w_std",
                                                 "Green space (std.)"="greens18_std",
                                                 "Lack green space (std.)"="lack_greens18_std",
                                                 "Noise (std.)"="noise_exposure_std",
                                                 "Multiple environmental burdens (std.)"="meb_std") %>%
           fct_relevel("Air pollution (std.)","Noise (std.)","Green space (std.)","Lack green space (std.)","Multiple environmental burdens (std.)"),
         outcome_viz=outcome %>% fct_recode("Share of SGBII recipients (std.)"="share_sgb_total_std",
                                            "Share of foreign minorities (std.)"="share_workingage_foreign_std") %>%
           fct_relevel("Share of SGBII recipients (std.)","Share of foreign minorities (std.)"),
         ags8_gen_num=factor(ags8_gen_num,levels=lvls))

## read NUTS1 shapes for maps
nuts_shape <- st_read(paste0(data_raw,"nuts250/250_NUTS1.shp")) %>%
  st_transform(32632)


## Map serving as legend; assigning city names to polygons on map
graph <- ggplot()+
  geom_sf(data = nuts_shape, fill="gray90",color="NA")+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_sgb_total_std" & term=="airpoll_w_std"),
          color="white",fill="gray70")+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_sgb_total_std" & term=="airpoll_w_std"),
          aes(color=ags8_gen_num),fill=NA)+
  geom_sf_text(data = coefs_slx%>% 
                 filter(outcome=="share_sgb_total_std" & term=="airpoll_w_std"),
               aes(label=city_number),size=2.5,color="black")+
  labs(color="City")+
  scale_color_manual(values = rep("NA",70))+
  theme_void()+
  guides(color=guide_legend(ncol = 2))

ggsave(graph, filename=paste0(figs,"map_city_legend.pdf"),width = 12, height = 11)
ggsave(graph, filename=paste0(figs,"map_city_legend.svg"),width = 12, height = 11)

###########################
##### SLX MAPS CITIES #####
###########################

a <- coefs_slx %>%
  mutate(vartype=case_when(
    str_detect(term,"lag") ~ "Spatial lag coefficient",
    TRUE ~ "SLX main coefficient")%>% as.factor() %>% fct_relevel("SLX main coefficient","Spatial lag coefficient")) %>%
  filter(outcome=="share_sgb_total_std" & term=="airpoll_w_std"| outcome=="share_sgb_total_std" & term=="airpoll_w_std_lag") %>%
  filter(vartype=="SLX main coefficient") %>%
  left_join(coefs_slx %>% st_drop_geometry() %>% filter(outcome=="share_sgb_total_std" & term=="airpoll_w_std")%>%
              mutate(rank=dense_rank(estimate)) %>% select(ags8_gen,rank),by="ags8_gen")%>%
  ggplot(aes(y=reorder(ags8_gen,rank),x=estimate,fill=estimate))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high),shape=21)+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid = "white",midpoint = 0)+
  labs(x="Coefficient estimate and 95% CI",y="Municipality (AGS8)")+
  guides(fill="none")+
  facet_wrap(~vartype,scales="free_x")+
  theme_bw()

b <- ggplot()+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_sgb_total_std" & term=="airpoll_w_std"),aes(fill=estimate))+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid="white",midpoint = 0,guide = guide_colorbar(title = "Air pollution\ncoefficient estimate"))+
  geom_sf(data = nuts_shape, fill=NA,color="black")+
  theme_void()+
  theme(legend.box.spacing = unit(-10, "pt"))



ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=FALSE,
                  legend="right",
                  widths = c(1,2)) %>% 
  ggsave(filename=paste0(figs,"sgbtotal_airpoll_slx.pdf"), height = 8,width = 14)

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=FALSE,
                  legend="right",
                  widths = c(1,1.5)) %>% 
  ggsave(filename=paste0(figs,"sgbtotal_airpoll_slx.svg"), height = 8,width = 13)


a <- coefs_slx %>%
  mutate(vartype=case_when(
    str_detect(term,"lag") ~ "Spatial lag coefficient",
    TRUE ~ "SLX main coefficient")%>% as.factor() %>% fct_relevel("SLX main coefficient","Spatial lag coefficient")) %>%
  filter(outcome=="share_sgb_total_std" & term=="noise_exposure_std"| outcome=="share_sgb_total_std" & term=="noise_exposure_std_lag") %>%
  filter(vartype=="SLX main coefficient") %>%
  left_join(coefs_slx %>% st_drop_geometry() %>% filter(outcome=="share_sgb_total_std" & term=="noise_exposure_std")%>%
              mutate(rank=dense_rank(estimate)) %>% select(ags8_gen,rank),by="ags8_gen")%>%
  ggplot(aes(y=reorder(ags8_gen,rank),x=estimate,fill=estimate))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high),shape=21)+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid = "white",midpoint = 0)+
  labs(x="Coefficient estimate and 95% CI",y="Municipality (AGS8)")+
  guides(fill="none")+
  facet_wrap(~vartype,scales="free_x")+
  theme_bw()

b <- ggplot()+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_sgb_total_std" & term=="noise_exposure_std"),aes(fill=estimate))+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid="white",midpoint = 0,guide = guide_colorbar(title = "Noise exposure\ncoefficient estimate"))+
  geom_sf(data = nuts_shape, fill=NA,color="black")+
  theme_void()+
  theme(legend.box.spacing = unit(-10, "pt"))

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=TRUE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"sgbtotal_noise_slx.pdf"), height = 8,width = 14)

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=FALSE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"sgbtotal_noise_slx.svg"), height = 8,width = 13)

a <- coefs_slx %>%
  mutate(vartype=case_when(
    str_detect(term,"lag") ~ "Spatial lag coefficient",
    TRUE ~ "SLX main coefficient")%>% as.factor() %>% fct_relevel("SLX main coefficient","Spatial lag coefficient")) %>%
  filter(outcome=="share_sgb_total_std" & term=="greens18_std"| outcome=="share_sgb_total_std" & term=="greens18_std_lag") %>%
  filter(vartype=="SLX main coefficient") %>%
  left_join(coefs_slx %>% st_drop_geometry() %>% filter(outcome=="share_sgb_total_std" & term=="greens18_std")%>%
              mutate(rank=dense_rank(estimate)) %>% select(ags8_gen,rank),by="ags8_gen")%>%
  ggplot(aes(y=reorder(ags8_gen,desc(rank)),x=estimate,fill=estimate))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high),shape=21)+
  scale_fill_gradient2(low = "red",high = "darkgreen",mid = "white",midpoint = 0)+
  labs(x="Coefficient estimate and 95% CI",y="Municipality (AGS8)")+
  guides(fill="none")+
  facet_wrap(~vartype,scales="free_x")+
  theme_bw()

b <- ggplot()+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_sgb_total_std" & term=="greens18_std"),aes(fill=estimate))+
  scale_fill_gradient2(low = "red",high = "darkgreen",mid="white",midpoint = 0,guide = guide_colorbar(title = "Green space\ncoefficient estimate"))+
  geom_sf(data = nuts_shape, fill=NA,color="black")+
  theme_void()+
  theme(legend.box.spacing = unit(-10, "pt"))

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=TRUE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"sgbtotal_greens_slx.pdf"), height = 8,width = 14)

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=FALSE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"sgbtotal_greens_slx.svg"), height = 8,width = 13)


a <- coefs_slx %>%
  mutate(vartype=case_when(
    str_detect(term,"lag") ~ "Spatial lag coefficient",
    TRUE ~ "SLX main coefficient")%>% as.factor() %>% fct_relevel("SLX main coefficient","Spatial lag coefficient")) %>%
  filter(outcome=="share_sgb_total_std" & term=="lack_greens18_std"| outcome=="share_sgb_total_std" & term=="lack_greens18_std_lag") %>%
  filter(vartype=="SLX main coefficient") %>%
  left_join(coefs_slx %>% st_drop_geometry() %>% filter(outcome=="share_sgb_total_std" & term=="lack_greens18_std")%>%
              mutate(rank=dense_rank(estimate)) %>% select(ags8_gen,rank),by="ags8_gen")%>%
  ggplot(aes(y=reorder(ags8_gen,rank),x=estimate,fill=estimate))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high),shape=21)+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid = "white",midpoint = 0)+
  labs(x="Coefficient estimate and 95% CI",y="Municipality (AGS8)")+
  guides(fill="none")+
  facet_wrap(~vartype,scales="free_x")+
  theme_bw()

b <- ggplot()+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_sgb_total_std" & term=="lack_greens18_std"),aes(fill=estimate))+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid="white",midpoint = 0,guide = guide_colorbar(title = "Lack of green space\ncoefficient estimate"))+
  geom_sf(data = nuts_shape, fill=NA,color="black")+
  theme_void()+
  theme(legend.box.spacing = unit(-10, "pt"))

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=TRUE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"sgbtotal_lack_greens_slx.pdf"), height = 8,width = 14)

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=FALSE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"sgbtotal_lack_greens_slx.svg"), height = 8,width = 13)

a <- coefs_slx %>%
  mutate(vartype=case_when(
    str_detect(term,"lag") ~ "Spatial lag coefficient",
    TRUE ~ "SLX main coefficient")%>% as.factor() %>% fct_relevel("SLX main coefficient","Spatial lag coefficient")) %>%
  filter(outcome=="share_sgb_total_std" & term=="meb_std"| outcome=="share_sgb_total_std" & term=="meb_std_lag") %>%
  filter(vartype=="SLX main coefficient") %>%
  left_join(coefs_slx %>% st_drop_geometry() %>% filter(outcome=="share_sgb_total_std" & term=="meb_std")%>%
              mutate(rank=dense_rank(estimate)) %>% select(ags8_gen,rank),by="ags8_gen")%>%
  ggplot(aes(y=reorder(ags8_gen,rank),x=estimate,fill=estimate))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high),shape=21)+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid = "white",midpoint = 0)+
  labs(x="Coefficient estimate and 95% CI",y="Municipality (AGS8)")+
  guides(fill="none")+
  facet_wrap(~vartype,scales = "free_x")+
  theme_bw()

b <- ggplot()+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_sgb_total_std" & term=="meb_std"),aes(fill=estimate))+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid="white",midpoint = 0,guide = guide_colorbar(title = "Multiple burdens\ncoefficient estimate"))+
  geom_sf(data = nuts_shape, fill=NA,color="black")+
  theme_void()+
  theme(legend.box.spacing = unit(-10, "pt"))

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=TRUE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"sgbtotal_meb_slx.pdf"), height = 8,width = 14)

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=FALSE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"sgbtotal_meb_slx.svg"), height = 8,width = 13)


a <- coefs_slx %>%
  mutate(vartype=case_when(
    str_detect(term,"lag") ~ "Spatial lag coefficient",
    TRUE ~ "SLX main coefficient")%>% as.factor() %>% fct_relevel("SLX main coefficient","Spatial lag coefficient")) %>%
  filter(outcome=="share_workingage_foreign_std" & term=="airpoll_w_std"| outcome=="share_workingage_foreign_std" & term=="airpoll_w_std_lag") %>%
  filter(vartype=="SLX main coefficient") %>%
  left_join(coefs_slx %>% st_drop_geometry() %>% filter(outcome=="share_workingage_foreign_std" & term=="airpoll_w_std")%>%
              mutate(rank=dense_rank(estimate)) %>% select(ags8_gen,rank),by="ags8_gen")%>%
  ggplot(aes(y=reorder(ags8_gen,rank),x=estimate,fill=estimate))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high),shape=21)+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid = "white",midpoint = 0)+
  labs(x="Coefficient estimate and 95% CI",y="Municipality (AGS8)")+
  guides(fill="none")+
  facet_wrap(~vartype,scales="free_x")+
  theme_bw()

b <- ggplot()+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_workingage_foreign_std" & term=="airpoll_w_std"),aes(fill=estimate))+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid="white",midpoint = 0,guide = guide_colorbar(title = "Air pollution\ncoefficient estimate"))+
  geom_sf(data = nuts_shape, fill=NA,color="black")+
  theme_void()+
  theme(legend.box.spacing = unit(-10, "pt"))

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=TRUE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"foreign_airpoll_slx.pdf"), height = 8,width = 14)

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=FALSE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"foreign_airpoll_slx.svg"), height = 8,width = 13)

a <- coefs_slx %>%
  mutate(vartype=case_when(
    str_detect(term,"lag") ~ "Spatial lag coefficient",
    TRUE ~ "SLX main coefficient")%>% as.factor() %>% fct_relevel("SLX main coefficient","Spatial lag coefficient")) %>%
  filter(outcome=="share_workingage_foreign_std" & term=="noise_exposure_std"| outcome=="share_workingage_foreign_std" & term=="noise_exposure_std_lag") %>%
  filter(vartype=="SLX main coefficient") %>%
  left_join(coefs_slx %>% st_drop_geometry() %>% filter(outcome=="share_workingage_foreign_std" & term=="noise_exposure_std")%>%
              mutate(rank=dense_rank(estimate)) %>% select(ags8_gen,rank),by="ags8_gen")%>%
  ggplot(aes(y=reorder(ags8_gen,rank),x=estimate,fill=estimate))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high),shape=21)+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid = "white",midpoint = 0)+
  labs(x="Coefficient estimate and 95% CI",y="Municipality (AGS8)")+
  guides(fill="none")+
  facet_wrap(~vartype,scales="free_x")+
  theme_bw()

b <- ggplot()+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_workingage_foreign_std" & term=="noise_exposure_std"),aes(fill=estimate))+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid="white",midpoint = 0,guide = guide_colorbar(title = "Noise exposure\ncoefficient estimate"))+
  geom_sf(data = nuts_shape, fill=NA,color="black")+
  theme_void()+
  theme(legend.box.spacing = unit(-10, "pt"))

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=TRUE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"foreign_noise_slx.pdf"), height = 8,width = 14)

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=FALSE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"foreign_noise_slx.svg"), height = 8,width = 13)

a <- coefs_slx %>%
  mutate(vartype=case_when(
    str_detect(term,"lag") ~ "Spatial lag coefficient",
    TRUE ~ "SLX main coefficient")%>% as.factor() %>% fct_relevel("SLX main coefficient","Spatial lag coefficient")) %>%
  filter(outcome=="share_workingage_foreign_std" & term=="greens18_std"| outcome=="share_workingage_foreign_std" & term=="greens18_std_lag") %>%
  filter(vartype=="SLX main coefficient") %>%
  left_join(coefs_slx %>% st_drop_geometry() %>% filter(outcome=="share_workingage_foreign_std" & term=="greens18_std")%>%
              mutate(rank=dense_rank(estimate)) %>% select(ags8_gen,rank),by="ags8_gen")%>%
  ggplot(aes(y=reorder(ags8_gen,desc(rank)),x=estimate,fill=estimate))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high),shape=21)+
  scale_fill_gradient2(low = "red",high = "darkgreen",mid = "white",midpoint = 0)+
  labs(x="Coefficient estimate and 95% CI",y="Municipality (AGS8)")+
  guides(fill="none")+
  facet_wrap(~vartype,scales="free_x")+
  theme_bw()

b <- ggplot()+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_workingage_foreign_std" & term=="greens18_std"),aes(fill=estimate))+
  scale_fill_gradient2(low = "red",high = "darkgreen",mid="white",midpoint = 0,guide = guide_colorbar(title = "Green space\ncoefficient estimate"))+
  geom_sf(data = nuts_shape, fill=NA,color="black")+
  theme_void()+
  theme(legend.box.spacing = unit(-10, "pt"))

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=TRUE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"foreign_greens_slx.pdf"), height = 8,width = 14)

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=FALSE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"foreign_greens_slx.svg"), height = 8,width = 13)

a <- coefs_slx %>%
  mutate(vartype=case_when(
    str_detect(term,"lag") ~ "Spatial lag coefficient",
    TRUE ~ "SLX main coefficient")%>% as.factor() %>% fct_relevel("SLX main coefficient","Spatial lag coefficient")) %>%
  filter(outcome=="share_workingage_foreign_std" & term=="lack_greens18_std"| outcome=="share_workingage_foreign_std" & term=="lack_greens18_std_lag") %>%
  filter(vartype=="SLX main coefficient") %>%
  left_join(coefs_slx %>% st_drop_geometry() %>% filter(outcome=="share_workingage_foreign_std" & term=="lack_greens18_std")%>%
              mutate(rank=dense_rank(estimate)) %>% select(ags8_gen,rank),by="ags8_gen")%>%
  ggplot(aes(y=reorder(ags8_gen,rank),x=estimate,fill=estimate))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high),shape=21)+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid = "white",midpoint = 0)+
  labs(x="Coefficient estimate and 95% CI",y="Municipality (AGS8)")+
  guides(fill="none")+
  facet_wrap(~vartype,scales="free_x")+
  theme_bw()

b <- ggplot()+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_workingage_foreign_std" & term=="lack_greens18_std"),aes(fill=estimate))+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid="white",midpoint = 0,guide = guide_colorbar(title = "Lack of green space\ncoefficient estimate"))+
  geom_sf(data = nuts_shape, fill=NA,color="black")+
  theme_void()+
  theme(legend.box.spacing = unit(-10, "pt"))

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=TRUE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"foreign_lack_greens_slx.pdf"), height = 8,width = 14)

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=FALSE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"foreign_lack_greens_slx.svg"), height = 8,width = 13)

a <- coefs_slx %>%
  mutate(vartype=case_when(
    str_detect(term,"lag") ~ "Spatial lag coefficient",
    TRUE ~ "SLX main coefficient")%>% as.factor() %>% fct_relevel("SLX main coefficient","Spatial lag coefficient")) %>%
  filter(outcome=="share_workingage_foreign_std" & term=="meb_std"| outcome=="share_workingage_foreign_std" & term=="meb_std_lag") %>%
  filter(vartype=="SLX main coefficient") %>%
  left_join(coefs_slx %>% st_drop_geometry() %>% filter(outcome=="share_workingage_foreign_std" & term=="meb_std")%>%
              mutate(rank=dense_rank(estimate)) %>% select(ags8_gen,rank),by="ags8_gen")%>%
  ggplot(aes(y=reorder(ags8_gen,rank),x=estimate,fill=estimate))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high),shape=21)+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid = "white",midpoint = 0)+
  labs(x="Coefficient estimate and 95% CI",y="Municipality (AGS8)")+
  guides(fill="none")+
  facet_wrap(~vartype,scales = "free_x")+
  theme_bw()

b <- ggplot()+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_workingage_foreign_std" & term=="meb_std"),aes(fill=estimate))+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid="white",midpoint = 0,guide = guide_colorbar(title = "Multiple burdens\ncoefficient estimate"))+
  geom_sf(data = nuts_shape, fill=NA,color="black")+
  theme_void()+
  theme(legend.box.spacing = unit(-10, "pt"))

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=TRUE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"foreign_meb_slx.pdf"), height = 8,width = 14)

ggpubr::ggarrange(a,b,nrow = 1, 
                  common.legend=FALSE,
                  legend="right",
                  widths = c(1,1.5)) %>% ggsave(filename=paste0(figs,"foreign_meb_slx.svg"), height = 8,width = 13)


### Combined figure for poverty rate and foreign minorities (MEB variable); Figure 3

a <- coefs_slx %>%
  mutate(vartype=case_when(
    str_detect(term,"lag") ~ "Spatial lag coefficient",
    TRUE ~ "SLX main coefficient")%>% as.factor() %>% fct_relevel("SLX main coefficient","Spatial lag coefficient")) %>%
  filter(outcome=="share_sgb_total_std" & term=="meb_std"| outcome=="share_sgb_total_std" & term=="meb_std_lag") %>%
  filter(vartype=="SLX main coefficient") %>%
  left_join(coefs_slx %>% st_drop_geometry() %>% filter(outcome=="share_sgb_total_std" & term=="meb_std")%>%
              mutate(rank=dense_rank(estimate)) %>% select(ags8_gen,rank),by="ags8_gen")%>%
  ggplot(aes(y=reorder(ags8_gen,rank),x=estimate,fill=estimate))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high),shape=21)+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid = "white",midpoint = 0)+
  labs(x="Coefficient estimate and 95% CI",y="Municipality (AGS8)")+
  guides(fill="none")+
  facet_wrap(~vartype,scales = "free_x")+
  theme_bw()+
  theme(strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12))

b <- ggplot()+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_sgb_total_std" & term=="meb_std"),aes(fill=estimate))+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid="white",midpoint = 0,guide = guide_colorbar(title = "Multiple burdens\ncoefficient estimate"))+
  geom_sf(data = nuts_shape, fill=NA,color="black")+
  theme_void()+
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box.spacing = unit(-20, "pt"))

fig3a <- ggpubr::ggarrange(a,b,nrow = 1, 
                           common.legend=TRUE,
                           legend="right",
                           widths = c(1,1.4),
                           labels = c("(a) Link between multiple environmental burdens and poverty rate"),
                           label.x = 0.3,
                           label.y = 0)



a <- coefs_slx %>%
  mutate(vartype=case_when(
    str_detect(term,"lag") ~ "Spatial lag coefficient",
    TRUE ~ "SLX main coefficient")%>% as.factor() %>% fct_relevel("SLX main coefficient","Spatial lag coefficient")) %>%
  filter(outcome=="share_workingage_foreign_std" & term=="meb_std"| outcome=="share_workingage_foreign_std" & term=="meb_std_lag") %>%
  filter(vartype=="SLX main coefficient") %>%
  left_join(coefs_slx %>% st_drop_geometry() %>% filter(outcome=="share_workingage_foreign_std" & term=="meb_std")%>%
              mutate(rank=dense_rank(estimate)) %>% select(ags8_gen,rank),by="ags8_gen")%>%
  ggplot(aes(y=reorder(ags8_gen,rank),x=estimate,fill=estimate))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high),shape=21)+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid = "white",midpoint = 0)+
  labs(x="Coefficient estimate and 95% CI",y="Municipality (AGS8)")+
  guides(fill="none")+
  facet_wrap(~vartype,scales = "free_x")+
  theme_bw()+
  theme(strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12))

b <- ggplot()+
  geom_sf(data = coefs_slx %>% 
            filter(outcome=="share_workingage_foreign_std" & term=="meb_std"),aes(fill=estimate))+
  scale_fill_gradient2(low = "darkgreen",high = "red",mid="white",midpoint = 0,guide = guide_colorbar(title = "Multiple burdens\ncoefficient estimate"))+
  geom_sf(data = nuts_shape, fill=NA,color="black")+
  theme_void()+
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box.spacing = unit(-20, "pt"))

fig3b <- ggpubr::ggarrange(a,b,nrow = 1, 
                           common.legend=TRUE,
                           legend="right",
                           widths = c(1,1.4),
                           labels = c("(b) Link between multiple environmental burdens and foreign minority rate"),
                           label.x = .1,
                           label.y = 0)


empty <- ggplot() + theme_void()

fig3 <- ggpubr::ggarrange(fig3a,empty,fig3b,empty,
                          nrow = 4,
                          heights = c(1,0.05,1,0.05),
                          common.legend = FALSE)


ggsave(filename = paste0(figs,"fig3_meb_slx.svg"),
       height = 16,
       width = 13)


#################################
##### CITY CONTEXT ANALYSES #####
#################################

## build data set of city-level SLX estimates of EI (main coefs) and city level contexts
EI_context <- coefs_slx %>%
  st_drop_geometry() %>%
  filter(!str_detect(term,"lag")) %>%
  mutate(se_squared=std.error^2,
         weights_uncertainty=1/se_squared)


city_variables <- haven::read_dta(paste0(data_raw,"city-level_indices.dta")) %>%
  mutate(ags=as.character(ags),
         ags=case_when(
           ags_name=="Hamburg" ~ "02000000",
           ags_name=="Berlin" ~ "11000000",
           TRUE ~ ags),
         ags=str_pad(ags,width=8,side="left",pad="0"))%>%
  select(ags_name,ags,bev_log,dissim_sgb,dissim_ausl,sgb_quote_stadt,ausl_quote_stadt,
         dichte_arm,dichte_ausl,akad_quote_stadt,Bev_ges)

municipality_eq <- readRDS(paste0(workdata,"municipality_eq_02.rds"))

city_level_centrality_robustness <- readRDS(paste0(workdata,"city_level_centrality_robustness.rds"))

aux <- city_variables %>% 
  inner_join(municipality_eq,by=c("ags_name"="GEN")) %>%
  inner_join(city_level_centrality_robustness, by=c("ags_name"="municipality", "ags"="AGS")) %>%
  mutate(across(
    .cols = c(bev_log,dissim_sgb,dissim_ausl,sgb_quote_stadt,ausl_quote_stadt,dichte_arm,dichte_ausl,
              airpoll_mun,greens_mun,lack_greens_mun,noise_mun,meb_mun,meb_mun_w),
    .fns = ~ scale(.x,center = T,scale = T),
    .names = "{.col}_std"))


EI_context <- EI_context %>% inner_join(aux,by=c("ags8_gen"="ags_name"))

aux <- data_core_full %>%
  select(AGS,Bev_km2) %>%
  group_by(AGS) %>%
  slice_head() %>%
  st_drop_geometry()

EI_context <- EI_context %>% inner_join(aux,by="AGS")


which(is.na(EI_context),arr.ind = TRUE)

## generate varlist specific to sub-analyses
outcomes <- c(rep("share_sgb_total_std",4),rep("share_workingage_foreign_std",4)) %>% as.list()
treatments <- rep(c("airpoll_w_std","greens18_std","noise_exposure_std","meb_std"),2) %>% as.list()
outcomes_mun <- c(rep("sgb_quote_stadt_std",4),rep("ausl_quote_stadt_std",4)) %>% as.list()
treatments_mun <- rep(c("airpoll_mun_std","greens_mun_std","noise_mun_std","meb_mun_std"),2) %>% as.list()
outcomes_density <- c(rep("dichte_arm_std",4),rep("dichte_ausl_std",4)) %>% as.list()
outcomes_dissim <- c(rep("dissim_sgb_std",4),rep("dissim_ausl_std",4)) %>% as.list()

outcomes_treatments_list <- Map(c,outcomes,treatments,outcomes_mun,treatments_mun,outcomes_density,outcomes_dissim)

x <- outcomes_treatments_list[[8]] 

context_list <- lapply(outcomes_treatments_list, function(x){
  out <- x[[1]] 
  treat <- x[[2]]
  out_mun <- x[[3]]
  treat_mun <- x[[4]] 
  out_density <- x[[5]]
  out_dissim <- x[[6]] 
  
  ## filter
  df <- EI_context %>% filter(outcome==out & treatment==treat)
  
  ## formula (all city level contexts)
  formula_m1 <- as.formula(paste0("estimate ~",treat_mun,"+",out_dissim,"+",out_density,"+","bev_log_std"))
  formula_m1b <- as.formula(paste0("estimate ~",treat_mun,"+",out_dissim,"+",out_density,"+",out_mun,"+","Bev_km2 + bev_log_std"))

  
  ## lm
  context_m1_uncertainty <- lm(data = df,
                               formula = formula_m1,
                               weights = weights_uncertainty) %>% 
    broom::tidy(conf.int=TRUE) %>%
    filter(term!="(Intercept)") %>%
    mutate(group=out,
           environmental_exposure=treat,
           model="multivariate_uncertainty")
  
  context_m1b_uncertainty <- lm(data = df,
                               formula = formula_m1b,
                               weights = weights_uncertainty) %>% 
    broom::tidy(conf.int=TRUE) %>%
    filter(term!="(Intercept)") %>%
    mutate(group=out,
           environmental_exposure=treat,
           model="multivariate_b_uncertainty")
  
  results <- rbind(context_m1_uncertainty,context_m1b_uncertainty) %>%
    mutate(term2=case_when(
      term==out_mun ~ "Group variable (city-level)",
      term==treat_mun ~ "Environmental exposure (city-level)",
      term==out_dissim ~ "Group segregation (city-level DI)",
      term==out_density ~ "Group-density-link (city-level)",
      term=="bev_log_std" ~ "City size",
      TRUE~ as.character(NA)) %>% as.factor() %>% fct_relevel("Group-density-link (city-level)",
                                                              "Group segregation (city-level DI)",
                                                              "City size",
                                                              "Environmental exposure (city-level)",
                                                              "Group variable (city-level)"))
  
  
  
  
  results
}) 

coefs_contexts <- do.call(what = rbind,args = context_list)

coefs_contexts <- coefs_contexts %>%
  filter(environmental_exposure %in% c("airpoll_w_std","greens18_std","meb_std","noise_exposure_std")) %>%
  mutate(environmental_exposure= environmental_exposure %>% fct_recode("Air pollution (std.)"="airpoll_w_std",
                                                                       "Green space (std.)"="greens18_std",
                                                                       "Noise (std.)"="noise_exposure_std",
                                                                       "Multiple env. burdens (std.)"="meb_std") %>%
           fct_relevel("Air pollution (std.)","Noise (std.)","Green space (std.)","Multiple env. burdens (std.)"),
         group=group %>% fct_recode("Poor vs. non-poor"="share_sgb_total_std",
                                    "Foreign minority vs. majority"="share_workingage_foreign_std") %>%
           fct_relevel("Poor vs. non-poor","Foreign minority vs. majority"),
         term2= case_when(
           group=="Poor vs. non-poor" & term2=="Group segregation (city-level DI)" ~ "Residential segregation\nof the poor",
           group=="Poor vs. non-poor" & term2=="Group-density-link (city-level)" ~ "Correlation between poverty share   \nand grid pop. density",
           group=="Poor vs. non-poor" & term2=="Environmental exposure (city-level)" ~ "Environmental exposure\n(aggregated)",
           group=="Foreign minority vs. majority" & term2=="Group segregation (city-level DI)" ~ "Residential segregation\nof foreign minorities",
           group=="Foreign minority vs. majority" & term2=="Group-density-link (city-level)" ~ "Correlation between share foreign   \nand grid pop. density",
           group=="Foreign minority vs. majority" & term2=="Environmental exposure (city-level)" ~ "Environmental exposure\n(aggregated)",
           TRUE ~ as.character(NA)) %>%
           fct_relevel("Correlation between poverty share   \nand grid pop. density","Correlation between share foreign   \nand grid pop. density",
                       "Environmental exposure\n(aggregated)",
                       "Residential segregation\nof the poor","Residential segregation\nof foreign minorities"))


# graph
library(ggh4x)

# uncertainty-weighted graph, additionally adjusting for outcome aggregated to municipality level and city population density
contexts_graph1 <- ggplot(coefs_contexts %>% 
                            filter(model=="multivariate_uncertainty" & 
                                     term2 %in% c("Residential segregation\nof the poor","Residential segregation\nof foreign minorities",
                                                  "Correlation between poverty share   \nand grid pop. density","Correlation between share foreign   \nand grid pop. density",
                                                  "Environmental exposure\n(aggregated)") &
                                     group=="Poor vs. non-poor"), 
                          aes(x=estimate,y=term2))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high,color=term2,fill=term2,shape=term2,group=term2),position = position_dodge2(width = .5),
                  size=1,
                  lwd=1.5)+
  scale_color_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_fill_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_shape_manual(limits=rev,values=c(25,16,25),guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  #scale_shape_discrete(guide=guide_legend(title = "Model"))+
  facet_grid2(group~environmental_exposure,scales = "free_x",switch = "y")+
  scale_x_continuous(limits=c(-.15,.15))+
  ggh4x::scale_x_facet(PANEL==1, limits=c(-.6,.6))+
  labs(x="")+
  theme_bw()+
  theme(text = element_text(size = 14,face="bold"),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=5,r=10,b=5,l=5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        legend.spacing.y = unit(1, "cm"))

contexts_graph2 <- ggplot(coefs_contexts %>% 
                            filter(model=="multivariate_uncertainty" & 
                                     term2 %in% c("Residential segregation\nof the poor","Residential segregation\nof foreign minorities",
                                                  "Correlation between poverty share   \nand grid pop. density","Correlation between share foreign   \nand grid pop. density",
                                                  "Environmental exposure\n(aggregated)") &
                                     group=="Foreign minority vs. majority"), 
                          aes(x=estimate,y=term2))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high,color=term2,fill=term2,shape=term2,group=term2),position = position_dodge2(width = .5),
                  size=1,
                  lwd=1.5)+
  scale_color_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_fill_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_shape_manual(limits=rev,values=c(22,16,22),guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  #scale_shape_discrete(guide=guide_legend(title = "Model"))+
  facet_grid2(group~environmental_exposure,scales = "free_x",switch = "y")+
  scale_x_continuous(limits=c(-.15,.15))+
  ggh4x::scale_x_facet(PANEL==1, limits=c(-.6,.6))+
  labs(x="Change in environmental inequality coefficient (SLX, city FEs) associated with city-level context variables")+
  theme_bw()+
  theme(text = element_text(size = 14,face="bold"),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=5,r=10,b=5,l=5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.x = element_blank(),
        legend.position = "right",
        legend.spacing.y = unit(1, "cm"))

ggpubr::ggarrange(contexts_graph1,contexts_graph2,
                  nrow = 2, align = "v") %>% ggsave(filename=paste0(figs,"city-level_predictors_uncertainty-weighted_EI_SLX_a.svg"),
                                                    width = 14,
                                                    height = 7)

ggpubr::ggarrange(contexts_graph1,contexts_graph2,
                  nrow = 2, align = "v") %>% ggsave(filename=paste0(figs,"city-level_predictors_uncertainty-weighted_EI_SLX_a.png"),
                                                    width = 14,
                                                    height = 7)




contexts_graph1 <- ggplot(coefs_contexts %>% 
                            filter(model=="multivariate_b_uncertainty" & 
                                     term2 %in% c("Residential segregation\nof the poor","Residential segregation\nof foreign minorities",
                                                  "Correlation between poverty share   \nand grid pop. density","Correlation between share foreign   \nand grid pop. density",
                                                  "Environmental exposure\n(aggregated)") &
                                     group=="Poor vs. non-poor"), 
                          aes(x=estimate,y=term2))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high,color=term2,fill=term2,shape=term2,group=term2),position = position_dodge2(width = .5),
                  size=1,
                  lwd=1.5)+
  scale_color_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_fill_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_shape_manual(limits=rev,values=c(25,16,25),guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  #scale_shape_discrete(guide=guide_legend(title = "Model"))+
  facet_grid2(group~environmental_exposure,scales = "free_x",switch = "y")+
  scale_x_continuous(limits=c(-.13,.18))+
  ggh4x::scale_x_facet(PANEL==1, limits=c(-.6,.6))+
  labs(x="")+
  theme_bw()+
  theme(text = element_text(size = 14,face="bold"),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=5,r=10,b=5,l=5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        legend.spacing.y = unit(1, "cm"))



contexts_graph2 <- ggplot(coefs_contexts %>% 
                            filter(model=="multivariate_b_uncertainty" & 
                                     term2 %in% c("Residential segregation\nof the poor","Residential segregation\nof foreign minorities",
                                                  "Correlation between poverty share   \nand grid pop. density","Correlation between share foreign   \nand grid pop. density",
                                                  "Environmental exposure\n(aggregated)") &
                                     group=="Foreign minority vs. majority"), 
                          aes(x=estimate,y=term2))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high,color=term2,fill=term2,shape=term2,group=term2),position = position_dodge2(width = .5),
                  size=1,
                  lwd=1.5)+
  scale_color_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_fill_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_shape_manual(limits=rev,values=c(22,16,22),guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  #scale_shape_discrete(guide=guide_legend(title = "Model"))+
  facet_grid2(group~environmental_exposure,scales = "free_x",switch = "y")+
  scale_x_continuous(limits=c(-.13,.18))+
  ggh4x::scale_x_facet(PANEL==1, limits=c(-.6,.6))+
  labs(x="Change in environmental inequality coefficient (SLX, city FEs) associated with city-level context variables")+
  theme_bw()+
  theme(text = element_text(size = 14,face="bold"),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=5,r=10,b=5,l=5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.x = element_blank(),
        legend.position = "right",
        legend.spacing.y = unit(1, "cm"))

ggpubr::ggarrange(contexts_graph1,contexts_graph2,
                  nrow = 2, align = "v") %>% ggsave(filename=paste0(figs,"city-level_predictors_uncertainty-weighted_EI_SLX_b.svg"),
                                                    width = 14,
                                                    height = 7)

ggpubr::ggarrange(contexts_graph1,contexts_graph2,
                  nrow = 2, align = "v") %>% ggsave(filename=paste0(figs,"city-level_predictors_uncertainty-weighted_EI_SLX_b.png"),
                                                    width = 14,
                                                    height = 7)

## Scatterplott, clusters of cities
gentrification_cities <- readRDS(file=paste0(workdata,"gentrification_indices_cities.rds"))

EI_context2 <- EI_context %>%
  mutate(city_size_akad_cat=case_when(
    Bev_ges>=500000 | akad_quote_stadt>=25 ~ ">500K inhabitants or >25% academics",
    TRUE ~ "All other cities")) %>% 
  inner_join(gentrification_cities, by="AGS")

graph <- ggplot(EI_context2 %>% filter(term=="meb_std" & outcome=="share_workingage_foreign_std"),aes(x=corr_acad_density,y=estimate)) +
  geom_hline(yintercept = 0,linetype="dashed",color="grey60")+
  geom_point(color="black",size=2.5,alpha=1)+
  geom_line(stat="smooth",method="lm",aes(color=after_stat(x)),lwd=1.5,alpha=.75)+
  geom_smooth(method="lm",alpha=.2,color=NA,se=TRUE)+
  scale_color_viridis_c(option = "inferno",begin = .9,end = .2)+
  guides(color="none")+
  labs(x="Correlation between grid-level share of academics and population density within city",y="Multiple environmental burdens (foreign vs German)")+
  theme_bw()+
  theme(legend.position = "bottom",
        text = element_text(face = "bold"))


ggsave(graph,file=paste0(figs,"scatterplot_EI_foreign_gentrification.pdf"),
       height = 5.5,width = 8)
ggsave(graph,file=paste0(figs,"scatterplot_EI_foreign_gentrification.png"),
       height = 5.5,width = 8)
ggsave(graph,file=paste0(figs,"scatterplot_EI_foreign_gentrification.svg"),
       height = 5.5,width = 8)




graph <- ggplot(EI_context2 %>% filter(term=="meb_std" & outcome=="share_sgb_total_std"),aes(x=corr_acad_density,y=estimate)) +
  geom_hline(yintercept = 0,linetype="dashed",color="grey60")+
  geom_point(color="black",size=2.5,alpha=1)+
  geom_line(stat="smooth",method="lm",aes(color=after_stat(x)),lwd=1.5,alpha=.75)+
  geom_smooth(method="lm",alpha=.2,color=NA,se=TRUE)+
  scale_color_viridis_c(option = "inferno",begin = .9,end = .2)+
  guides(color="none")+
  labs(x="Correlation between grid-level share of academics and population density within city",y="Multiple environmental burdens (poor vs non-poor)")+
  theme_bw()+
  theme(legend.position = "bottom",
        text = element_text(face = "bold"))


ggsave(graph,file=paste0(figs,"scatterplot_EI_sgb_gentrification.pdf"),
       height = 5.5,width = 8)

ggsave(graph,file=paste0(figs,"scatterplot_EI_sgb_gentrification.png"),
       height = 5.5,width = 8)

ggsave(graph,file=paste0(figs,"scatterplot_EI_sgb_gentrification.svg"),
       height = 5.5,width = 8)



#############################################################
##### ROBUSTNESS CHECK: ALTERNATIVE CENTRALITY MEASURES #####
#############################################################

## generate varlist specific to sub-analyses
outcomes <- c(rep("share_sgb_total_std",4),rep("share_workingage_foreign_std",4)) %>% as.list()
treatments <- rep(c("airpoll_w_std","greens18_std","noise_exposure_std","meb_std"),2) %>% as.list()
outcomes_mun <- c(rep("sgb_quote_stadt_std",4),rep("ausl_quote_stadt_std",4)) %>% as.list()
treatments_mun <- rep(c("airpoll_mun_std","greens_mun_std","noise_mun_std","meb_mun_std"),2) %>% as.list()
outcomes_density <- c(rep("corr_sgb2_distcityhall",4),rep("corr_foreign_distcityhall",4)) %>% as.list()
outcomes_dissim <- c(rep("dissim_sgb_std",4),rep("dissim_ausl_std",4)) %>% as.list()

outcomes_treatments_list <- Map(c,outcomes,treatments,outcomes_mun,treatments_mun,outcomes_density,outcomes_dissim)

x <- outcomes_treatments_list[[8]] 

context_list_robustness_cityhalls <- lapply(outcomes_treatments_list, function(x){
  out <- x[[1]] 
  treat <- x[[2]]
  out_mun <- x[[3]]
  treat_mun <- x[[4]] 
  out_density <- x[[5]]
  out_dissim <- x[[6]] 
  
  ## filter
  df <- EI_context %>% filter(outcome==out & treatment==treat)
  
  ## formula (all city level contexts)
  formula_m1 <- as.formula(paste0("estimate ~",treat_mun,"+",out_dissim,"+",out_density,"+","bev_log_std"))
  formula_m1b <- as.formula(paste0("estimate ~",treat_mun,"+",out_dissim,"+",out_density,"+",out_mun,"+","Bev_km2 + bev_log_std"))
  
  
  ## lm
  context_m1_uncertainty <- lm(data = df,
                               formula = formula_m1,
                               weights = weights_uncertainty) %>% 
    broom::tidy(conf.int=TRUE) %>%
    filter(term!="(Intercept)") %>%
    mutate(group=out,
           environmental_exposure=treat,
           model="multivariate_uncertainty")
  
  context_m1b_uncertainty <- lm(data = df,
                                formula = formula_m1b,
                                weights = weights_uncertainty) %>% 
    broom::tidy(conf.int=TRUE) %>%
    filter(term!="(Intercept)") %>%
    mutate(group=out,
           environmental_exposure=treat,
           model="multivariate_b_uncertainty")
  
  results <- rbind(context_m1_uncertainty,context_m1b_uncertainty) %>%
    mutate(term2=case_when(
      term==out_mun ~ "Group variable (city-level)",
      term==treat_mun ~ "Environmental exposure (city-level)",
      term==out_dissim ~ "Group segregation (city-level DI)",
      term==out_density ~ "Group-centrality-link, city hall (city-level)",
      term=="bev_log_std" ~ "City size",
      TRUE~ as.character(NA)) %>% as.factor() %>% fct_relevel("Group-centrality-link, city hall (city-level)",
                                                              "Group segregation (city-level DI)",
                                                              "City size",
                                                              "Environmental exposure (city-level)",
                                                              "Group variable (city-level)"))
  
  
  
  
  results
}) 

coefs_contexts <- do.call(what = rbind,args = context_list_robustness_cityhalls)

coefs_contexts <- coefs_contexts %>%
  filter(environmental_exposure %in% c("airpoll_w_std","greens18_std","meb_std","noise_exposure_std")) %>%
  mutate(environmental_exposure= environmental_exposure %>% fct_recode("Air pollution (std.)"="airpoll_w_std",
                                                                       "Green space (std.)"="greens18_std",
                                                                       "Noise (std.)"="noise_exposure_std",
                                                                       "Multiple env. burdens (std.)"="meb_std") %>%
           fct_relevel("Air pollution (std.)","Noise (std.)","Green space (std.)","Multiple env. burdens (std.)"),
         group=group %>% fct_recode("Poor vs. non-poor"="share_sgb_total_std",
                                    "Foreign minority vs. majority"="share_workingage_foreign_std") %>%
           fct_relevel("Poor vs. non-poor","Foreign minority vs. majority"),
         term2= case_when(
           group=="Poor vs. non-poor" & term2=="Group segregation (city-level DI)" ~ "Residential segregation\nof the poor",
           group=="Poor vs. non-poor" & term2=="Group-centrality-link, city hall (city-level)" ~ "Correlation between poverty share   \nand dist. to city hall",
           group=="Poor vs. non-poor" & term2=="Environmental exposure (city-level)" ~ "Environmental exposure\n(aggregated)",
           group=="Foreign minority vs. majority" & term2=="Group segregation (city-level DI)" ~ "Residential segregation\nof foreign minorities",
           group=="Foreign minority vs. majority" & term2=="Group-centrality-link, city hall (city-level)" ~ "Correlation between share foreign   \nand dist. to city hall",
           group=="Foreign minority vs. majority" & term2=="Environmental exposure (city-level)" ~ "Environmental exposure\n(aggregated)",
           TRUE ~ as.character(NA)) %>%
           fct_relevel("Correlation between poverty share   \nand dist. to city hall","Correlation between share foreign   \nand dist. to city hall",
                       "Environmental exposure\n(aggregated)",
                       "Residential segregation\nof the poor","Residential segregation\nof foreign minorities"))


# graph
library(ggh4x)

# uncertainty-weighted graph, additionally adjusting for outcome aggregated to municipality level and city population density
contexts_graph1 <- ggplot(coefs_contexts %>% 
                            filter(model=="multivariate_uncertainty" & 
                                     term2 %in% c("Residential segregation\nof the poor","Residential segregation\nof foreign minorities",
                                                  "Correlation between poverty share   \nand dist. to city hall","Correlation between share foreign   \nand dist. to city hall",
                                                  "Environmental exposure\n(aggregated)") &
                                     group=="Poor vs. non-poor"), 
                          aes(x=estimate,y=term2))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high,color=term2,fill=term2,shape=term2,group=term2),position = position_dodge2(width = .5),
                  size=1,
                  lwd=1.5)+
  scale_color_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_fill_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_shape_manual(limits=rev,values=c(25,16,25),guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  #scale_shape_discrete(guide=guide_legend(title = "Model"))+
  facet_grid2(group~environmental_exposure,scales = "free_x",switch = "y")+
  #scale_x_continuous(limits=c(-.15,.15))+
  #ggh4x::scale_x_facet(PANEL==1, limits=c(-.6,.6))+
  labs(x="")+
  theme_bw()+
  theme(text = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=5,r=10,b=5,l=5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "plain"),
        legend.position = "right",
        legend.spacing.y = unit(1, "cm"))

contexts_graph2 <- ggplot(coefs_contexts %>% 
                            filter(model=="multivariate_uncertainty" & 
                                     term2 %in% c("Residential segregation\nof the poor","Residential segregation\nof foreign minorities",
                                                  "Correlation between poverty share   \nand dist. to city hall","Correlation between share foreign   \nand dist. to city hall",
                                                  "Environmental exposure\n(aggregated)") &
                                     group=="Foreign minority vs. majority"), 
                          aes(x=estimate,y=term2))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high,color=term2,fill=term2,shape=term2,group=term2),position = position_dodge2(width = .5),
                  size=1,
                  lwd=1.5)+
  scale_color_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_fill_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_shape_manual(limits=rev,values=c(22,16,22),guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  #scale_shape_discrete(guide=guide_legend(title = "Model"))+
  facet_grid2(group~environmental_exposure,scales = "free_x",switch = "y")+
  #scale_x_continuous(limits=c(-.15,.15))+
  #ggh4x::scale_x_facet(PANEL==1, limits=c(-.6,.6))+
  labs(x="Change in environmental inequality coefficient (SLX, city FEs) associated with city-level context variables")+
  theme_bw()+
  theme(text = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=5,r=10,b=5,l=5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "plain"),
        strip.text.x = element_blank(),
        legend.position = "right",
        legend.spacing.y = unit(1, "cm"))

ggpubr::ggarrange(contexts_graph1,contexts_graph2,
                  nrow = 2, align = "v") %>% ggsave(filename=paste0(figs,"city-level_predictors_uncertainty-weighted_EI_SLX_a_robustness_cityhalls.svg"),
                                                    width = 16,
                                                    height = 8)

ggpubr::ggarrange(contexts_graph1,contexts_graph2,
                  nrow = 2, align = "v") %>% ggsave(filename=paste0(figs,"city-level_predictors_uncertainty-weighted_EI_SLX_a_robustness_cityhalls.png"),
                                                    width = 16,
                                                    height = 8)




contexts_graph1 <- ggplot(coefs_contexts %>% 
                            filter(model=="multivariate_b_uncertainty" & 
                                     term2 %in% c("Residential segregation\nof the poor","Residential segregation\nof foreign minorities",
                                                  "Correlation between poverty share   \nand dist. to city hall","Correlation between share foreign   \nand dist. to city hall",
                                                  "Environmental exposure\n(aggregated)") &
                                     group=="Poor vs. non-poor"), 
                          aes(x=estimate,y=term2))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high,color=term2,fill=term2,shape=term2,group=term2),position = position_dodge2(width = .5),
                  size=1,
                  lwd=1.5)+
  scale_color_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_fill_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_shape_manual(limits=rev,values=c(25,16,25),guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  #scale_shape_discrete(guide=guide_legend(title = "Model"))+
  facet_grid2(group~environmental_exposure,scales = "free_x",switch = "y")+
  #scale_x_continuous(limits=c(-.13,.18))+
  #ggh4x::scale_x_facet(PANEL==1, limits=c(-.6,.6))+
  labs(x="")+
  theme_bw()+
  theme(text = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=5,r=10,b=5,l=5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "plain"),
        legend.position = "right",
        legend.spacing.y = unit(1, "cm"))

contexts_graph2 <- ggplot(coefs_contexts %>% 
                            filter(model=="multivariate_b_uncertainty" & 
                                     term2 %in% c("Residential segregation\nof the poor","Residential segregation\nof foreign minorities",
                                                  "Correlation between poverty share   \nand dist. to city hall","Correlation between share foreign   \nand dist. to city hall",
                                                  "Environmental exposure\n(aggregated)") &
                                     group=="Foreign minority vs. majority"), 
                          aes(x=estimate,y=term2))+
  geom_vline(xintercept = 0)+
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high,color=term2,fill=term2,shape=term2,group=term2),position = position_dodge2(width = .5),
                  size=1,
                  lwd=1.5)+
  scale_color_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_fill_viridis_d(limits=rev,option = "plasma", begin = .8, end = 0, guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  scale_shape_manual(limits=rev,values=c(22,16,22),guide=guide_legend(title = "City-level context\nvariables (std.)"))+
  #scale_shape_discrete(guide=guide_legend(title = "Model"))+
  facet_grid2(group~environmental_exposure,scales = "free_x",switch = "y")+
  #scale_x_continuous(limits=c(-.13,.18))+
  #ggh4x::scale_x_facet(PANEL==1, limits=c(-.6,.6))+
  labs(x="Change in environmental inequality coefficient (SLX, city FEs) associated with city-level context variables")+
  theme_bw()+
  theme(text = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=5,r=10,b=5,l=5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "plain"),
        strip.text.x = element_blank(),
        legend.position = "right",
        legend.spacing.y = unit(1, "cm"))

ggpubr::ggarrange(contexts_graph1,contexts_graph2,
                  nrow = 2, align = "v") %>% ggsave(filename=paste0(figs,"city-level_predictors_uncertainty-weighted_EI_SLX_b_robustness_cityhalls.svg"),
                                                    width = 16,
                                                    height = 8)

ggpubr::ggarrange(contexts_graph1,contexts_graph2,
                  nrow = 2, align = "v") %>% ggsave(filename=paste0(figs,"city-level_predictors_uncertainty-weighted_EI_SLX_b_robustness_cityhalls.png"),
                                                    width = 16,
                                                    height = 8)



