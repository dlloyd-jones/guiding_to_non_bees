# 

# Dangerous mutualists: on the guiding of humans by greater honeyguides to animals other than bees
# R code for running all analyses included in the manuscript submitted to The American Naturalist

# Data are stored in stored in four files: 
# (1) guided_trips_with_gps.xlsx
# (2) prior_harvests_data.csv
# (3) all_guides_2018.csv"
# (4) hg_selections.csv
# and two folders:
# (1) gpx_guiding_tracks (containing gpx.tracks)
# (2) honeyguide_audio_clips (containing wav.files)

library(readxl)
library(tidyverse)# load packages
library(ggplot2)
library(ggsci)
library(lubridate)
library(hms)
library(forcats)
library(knitr)
library(plotrix)
library(lme4)
library(trajr)
library(ggmap)
library(rstudioapi)
library(here)
library(RCurl)
library(XML)
library(geosphere)
library(ggplot2)
library(ggspatial)
library(rstatix)
library(dplyr)
library(ggfortify)
library(warbleR)
library(seewave)
library(corrplot)

options(digits.secs = 3)
options(scipen = 999)

nonbees <- read_excel("~/Desktop/Am Nat Submission/guided_trips_with_gps.xlsx", 
                            col_types = c("numeric", "text", "text", "text", 
                                          "text", "text", "text", "text", "text", 
                                          "date", "date", "date","numeric", 
                                          "numeric", "date", "date", "date", 
                                          "numeric", "numeric", "numeric", 
                                          "date", "date", "numeric"))


nonbees$date <- as.Date(nonbees$date, origin="1899-12-30")
nonbees$audio_start_indication <- as_hms(ymd_hms(nonbees$audio_start_indication)) # covert time columns to postix:
nonbees$audio_finish_woo <- as_hms(ymd_hms(nonbees$audio_finish_woo))
nonbees$indication_to_woo <- as_hms(ymd_hms(nonbees$indication_to_woo))
nonbees$guiding_starttime <- as_hms(ymd_hms(nonbees$guiding_starttime))
nonbees$finalphasecalls_time <- as_hms(ymd_hms(nonbees$finalphasecalls_time))
nonbees$foundtree_time <- as_hms(ymd_hms(nonbees$foundtree_time))
nonbees$GPSguiding_starttime <- as_hms(ymd_hms(nonbees$GPSguiding_starttime))
nonbees$GPSguiding_endtime <- as_hms(ymd_hms(nonbees$GPSguiding_endtime))

# ---- Figure 2 A and B ---- 

# this step is set up to be manually iterated over the twenty-four GPX tracks in the 'tracks' folder

setwd("~/gpx_guiding_tracks")
# register API
register_google(key = "INSERT API KEY")

# parse GPX file
path_tmp <- paste("~/gpx_guiding_tracks/C20.gpx")
parsed <- htmlTreeParse(file = path_tmp, useInternalNodes = TRUE)
#parsed

# get values via the respective path
coords <- xpathSApply(parsed, path = "//trkpt", xmlAttrs)
elev   <- xpathSApply(parsed, path = "//trkpt/ele", xmlValue)
ts_chr <- xpathSApply(parsed, path = "//trkpt/time", xmlValue)

# combine into dataframe 
dat_df <- data.frame(
  ts_POSIXct = ymd_hms(ts_chr, tz = "Africa/Maputo"),
  lat = as.numeric(coords["lat",]), 
  lon = as.numeric(coords["lon",]), 
  elev = as.numeric(elev)
)
head(dat_df)

dat_df <- 
  dat_df %>%
  mutate(lat_lead = lead(lat)) %>%
  mutate(lon_lead = lead(lon)) %>%
  rowwise() %>%
  mutate(dist_to_lead_m = distm(c(lon, lat), c(lon_lead, lat_lead), fun = distHaversine)[1,1]) %>%
  ungroup()
dat_df

# get the map background (valid GoogleMaps API required)
dat_df_map <- get_googlemap(center = c(mean(range(dat_df$lon)), mean(range(dat_df$lat))), zoom = 17)

# data frame to add distance marks
dat_df_dist_marks <- 
  dat_df %>% 
  mutate(dist_m_cumsum = cumsum(dist_to_lead_m)) %>%
  mutate(dist_m_cumsum_km_floor = floor(dist_m_cumsum / 1000)) %>%
  group_by(dist_m_cumsum_km_floor) %>%
  dplyr::filter(row_number() == 1, dist_m_cumsum_km_floor > 0) 

# generate plot
plt_path_fancy <- 
  ggmap(dat_df_map) + 
  theme_classic()+
  geom_path(data = dat_df, aes(lon,lat,), linewidth = 0.5,colour = "#ff2c1c")+
  #colour = "#0a68cd"
  #scale_colour_gradientn(colours = topo.colors (4))+
  scale_colour_gradientn(colours = terrain.colors (4))+
  #scale_colour_gradientn(colours = heat.colors (4))+
  #geom_path(data = dat_df, aes(lon, lat),
  #          size = 0.3) +
  geom_label(data = dat_df_dist_marks, aes(lon, lat, label = dist_m_cumsum_km_floor),
             size = 3) +
  labs(x = "Longitude", 
       y = "Latitude", 
       color = "Speed km[hr]",
       title = "C20 Track") # change for each track
plt_path_fancy
# SAVE IMAGE AS: 6 x 6 inches

# ---- Figure 3A  ---- 

fcd <- nonbees
fcd <- fcd %>% dplyr::filter(!treeID == "NA")
fcd$duration = fcd$GPSguiding_endtime-fcd$GPSguiding_starttime
fcd$indi = fcd$GPSguiding_endtime-fcd$audio_start_indication
fcd$difference =ifelse(is.na(fcd$indi),fcd$duration,ifelse(is.na(fcd$duration),fcd$indi,fcd$duration-fcd$indi))
fcd$start <- sample(2, size = nrow(fcd), replace = TRUE)
fcd["start"][!is.na(fcd["start"])] <- 0


fcd <- fcd %>% select(follow_GPStrackID,duration,difference,indi,start)
fcd <- fcd %>% mutate(follow_GPStrackID = fct_relevel(follow_GPStrackID,"C34","C43_PY","C20","C16","C45","C39","C13","C22","C21","C33",
                                                          "C6","C19","C47","C14","C36","C49","C27","C3","C7","C5","C10","C46"))

color1 <- c("#0366cc")
color2 <- c("#fba200")
ggplot(fcd, aes(x=start, xend=difference, y=follow_GPStrackID, yend=follow_GPStrackID)) + 
  geom_segment(aes(colour = color1), linewidth =3, alpha = 1, show.legend = FALSE)+
  geom_segment(aes(x=difference, xend=duration, y=follow_GPStrackID, yend=follow_GPStrackID, color = color2),
               linewidth = 3, alpha = 1, show.legend = FALSE) +
  xlab("time (minutes)") +
  theme_classic()
# export as 4 x 6 inches


# ----- Calculation of overall guiding rates to non-bees in 2018 in section 'Rarity of guiding to animals other than bees' ----

all_guides_2018 <- read.csv("~/all_guides_2018.csv")

# calculations of visit rates:
4/108 # 4 guiding events to non-bee animals out of 108 guiding events 


# ---- Figure 4 ----

prior_harvests_data <- read.csv("~/prior_harvests_data.csv") # load data 
prop <- prior_harvests_data
prop$proportionNonHarvest <- replace(prop$proportionNonHarvest, which(prop$proportionNonHarvest == 0), 0.01)

ggplot(prop) +
  aes(x = outcome, y = proportionNonHarvest,fill = outcome) +
  geom_boxplot(width = .5)+
  coord_cartesian(ylim = c(0, 1))+
  theme_classic()+
  geom_jitter(width = 0.15, height = 0)+
  coord_flip()

# Shapiro-Wilk Test for Normality
shapiro_test <- shapiro.test(prop$proportionNonHarvest)

# Histogram
hist(prop$proportionNonHarvest, main = "Histogram of Outcome", col = "lightblue", border = "black")

# Q-Q Plot
qqnorm(prop$proportionNonHarvest)
qqline(prop$proportionNonHarvest, col = 2)

# Display results of Shapiro-Wilk test
cat("Shapiro-Wilk Test for Normality:\n")
cat("Test Statistic =", shapiro_test$statistic, "\n")
cat("p-value =", shapiro_test$p.value, "\n")

# Interpretation: The data is not normally distributed based on the Shapiro-Wilk test.
# Run a permutation test to see whether this is any statistical difference between the proportions of rewarding
# (proportionNonHarvest) prior to being guided to bees vs animals other than bees ("NonBees").

Bees <- c(0.50, 1.00, 0.5, 0.01, 0.01, 1.00, 1.00)  
NonBees <- c(0.50, 0.01, 0.75, 0.01)  

# Combine the data
all_proportions <- c(Bees, NonBees)
group_indicator <- rep(c("Bees", "NonBees"))

# Calculate observed difference in proportions
obs_diff <- mean(Bees) - mean(NonBees)

# Number of permutations
n_permutations <- 1000

# Function to calculate the difference in proportions for each permutation
perm_test <- function(data, indicator) {
  group1_perm <- sample(data[indicator == "Bees"], sum(indicator == "Bees"))
  group2_perm <- sample(data[indicator == "NonBees"], sum(indicator == "NonBees"))
  return(mean(group1_perm) - mean(group2_perm))
}

# Permutation test
set.seed(123)  # Set seed for reproducibility
perm_diffs <- replicate(n_permutations, perm_test(all_proportions, group_indicator))

# Calculate p-value
p_value <- mean(abs(perm_diffs) >= abs(obs_diff))


# ---- Statistics for section: 'Honeyguide behavior when guiding to animals is similar to when guiding to bees' ----

# Calculation of sinuosity
nonbees <- nonbees[complete.cases(nonbees$treeID), ] 
nonbees$direct <- nonbees$follow_trackdist/nonbees$followdist_crowfly # make a new column
ranges_ratio <- nonbees %>%
  group_by(group) %>%
  summarize(min_measure = min(direct),
            max_measure = max(direct),
            range_measure = max_measure - min_measure)
print(ranges_ratio)

# Mann-Whitney U test of a difference between track sinuosity to bees vs other animals

sinu_bees <- nonbees$direct[nonbees$bee == "apis"]
sinu_nonbees <- nonbees$direct[nonbees$bee == "NB"|nonbees$bee == "C14_BB"|nonbees$bee == "MAMBA"|nonbees$bee == "PYTHON"]

mwu_result <- wilcox.test(sinu_bees, sinu_nonbees)
print(mwu_result)

mean(sinu_bees) # calculate mean
mean(sinu_nonbees)
sd(sinu_bees) / sqrt(length(sinu_bees)) # calculate SE
sd(sinu_nonbees) / sqrt(length(sinu_nonbees)) # calculate SE


# ---- Calculations of guiding distances ----

gdist <- nonbees
gdist <- gdist %>% dplyr::filter(group == "B") # guided to bees group 
mean(gdist$follow_trackdist, na.rm = TRUE)
std.error(gdist$follow_trackdist,na.rm = TRUE)
range(gdist$follow_trackdist,na.rm = TRUE)

gdist <- nonbees
gdist <- gdist %>% dplyr::filter(group == "A") # guided to animals other than bees
mean(gdist$follow_trackdist, na.rm = TRUE)
std.error(gdist$follow_trackdist,na.rm = TRUE)
range(gdist$follow_trackdist,na.rm = TRUE)


# ---- Stats for Result section "Non-rewarding behavior by honey-hunters is not more likely to precede being guided to non-bee animals"

prior_h <- prop

prior_h_nh <- prior_h %>% dplyr::filter(outcome == "NonBees") # guided to bees group 
mean(prior_h_nh$proportionNonHarvest, na.rm = TRUE)
std.error(prior_h_nh$proportionNonHarvest,na.rm = TRUE)
range(prior_h_nh$proportionNonHarvest,na.rm = TRUE)

prior_h_h <- prior_h %>% dplyr::filter(outcome == "Bees") # guided to bees group 
mean(prior_h_h$proportionNonHarvest, na.rm = TRUE)
std.error(prior_h_h$proportionNonHarvest,na.rm = TRUE)
range(prior_h_h$proportionNonHarvest,na.rm = TRUE)


# ---- Supplementary material: acoustic analyses ----

setwd("~/honeyguide_audio_clips")
# ls("package:warbleR")
list.files()

# code for making spectrograms of audio files (one at a time)
freq_range_detec(readWave("C7bees_ind.wav"), flim = c(0, 10), bp = c(0, 3), threshold = 5, plot = TRUE) # manually change the name of the wav. audio clips to visualize
w2 <- readWave("34mamba_ind.wav")
spectro(w2, wl = 300, grid = FALSE, contlevels=seq(1,20), scale = FALSE, ovlp = 90, flim=c(0,10),) # visualize with spectrograph only

# load in time parameters for all signals in all audio clips:
hg.ad <- read.csv("~/hg_selections.csv")
hg.ad <- selection_table(X = hg.ad)

# uncomment the following to make jpegs with which to inspect spectrograms to see if selections are correctly applied
# full_spectrograms(hg.ad, wl = 300, flim = c(1, 10), ovlp = 10, sxrow = 3, rows = 3, it = "jpeg")

# measure all freguency parameters for all files
params <- spectro_analysis(hg.ad, bp = c(2, 10), threshold = 15) 

params_cor <- params[, 3:28] # filter columns needed
cor_matrix <- cor(params_cor) # inspect the correlation matrix 
print(cor_matrix)
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust") # plot correlation matrix using library(corplot)

#remove variables for pca
col_rm <- c("sp.ent", "startdom", "enddom", "time.Q25", "time.Q75", "freq.Q25", "freq.Q75") # remove columns which are highly correlated (> 0.75) with other columns
params <- params[, !(names(params) %in% col_rm)]
params <- params[, grep("fun|peakf", colnames(params), invert = TRUE)]
# write.csv(params, "acoustic_parameters.csv", row.names = FALSE)
pca <- prcomp(x = params[, sapply(params, is.numeric)], scale. = TRUE)
summary(pca)

pcascor <- as.data.frame(pca[[5]])

# Plot the 2 first PCs for supplementary material
plot(pcascor[, 1], pcascor[, 2], col = as.numeric(as.factor(params$sound.files)), pch = 20, 
     cex = 1, xlab = "PC1", ylab = "PC2")

params$nonbees <- hg.ad$nonbees
params$chatter <- hg.ad$chatter
params$both <- hg.ad$both

#plotting pca
colors <- c("#0366cc","#fba200", "#0366cc","#fba200")  # Specify your desired colors
shapes <- c(17,17,19,19)  # Specify your desired shapes
linetypes <- c("solid", "solid", "dashed", "dashed")

pca_plot <- autoplot(pca, data = params, colour = "both", shape = "both") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  theme_classic()
pca_plot
# add ellipses

params$chatter <- as.factor(params$chatter)
pca_plot + stat_ellipse(aes(color = both, linetype = both), level = 0.95, geom = "polygon", alpha = 0)+
  scale_linetype_manual(values = linetypes)

# save as 4 x 6 inches










