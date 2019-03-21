####### TEST zone #####


#### Example 1 - Simple analysis ####
# Example can be found in the supplementary material 1 - Bolduc et al 2017
set.seed(91)
library(R2MCDS)
### Import and filter data
data(alcidae)
df1 <-mcds.filter(alcidae[alcidae$WatchID == -1946788232,], transect.id = "WatchID",
                  distance.field = "Distance",
                  distance.labels = c("A", "B", "C", "D"),
                  distance.midpoints = c(25, 75, 150, 250),
                  effort.field = "WatchLenKm",
                  lat.field = "LatStart",
                  long.field = "LongStart",
                  sp.field = "Alpha",
                  date.field = "Date")
### Run analysis with the MCDS engine.
### Here, the 5-minute observation period (WatchID) is used as the sampling unit.
mod1 <- mcds.wrap(df1,
                  SMP_EFFORT="WatchLenKm",
                  DISTANCE="Distance",
                  SIZE="Count",
                  units=list(Type="Line",
                             Distance="Perp",
                             Length_units="Kilometers",
                             Distance_units="Meters",
                             Area_units="Square kilometers"),
                  breaks=c(0,50,100,200,300),
                  SMP_LABEL="WatchID",
                  STR_LABEL="STR_LABEL",
                  STR_AREA="STR_AREA",
                  estimator=list(c("HN","CO")),
                  multiplier = c(1, 0, 0),
                  path="C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/R_examples",
                  pathMCDS="C:/Program Files (x86)/Distance 7",
                  verbose=FALSE)

# File with errors = "log_xxx.tmp"
mod1
summary(mod1)
predicted_hist_new(mod1)

### Run all six detection keys
listmod <- mcds.wrap.point(df1,
                     SMP_EFFORT="WatchLenKm",
                     DISTANCE="Distance",
                     SIZE="Count",
                     units=list(Type="Point",
                                Distance="Radial",
                                Length_units="Kilometers",
                                Distance_units="Meters",
                                Area_units="Square kilometers"),
                     breaks=c(0,50,100,200,300),
                     STR_LABEL="STR_LABEL",
                     STR_AREA="STR_AREA",
                     estimator = NULL, # function fits all six detection models and returns a distanceList object that contains 6 distanceFit objects
                     multiplier = c(1, 0, 0),
                     SMP_LABEL="WatchID",
                     path="C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/R_examples",
                     pathMCDS="C:/Program Files (x86)/Distance 7",
                     verbose=FALSE)
listmod
summary(listmod[[1]])
lapply(listmod, summary)

predicted_hist_new(listmod[[6]])
lapply(listmod, predicted_hist)

##### Keep the 'best' model in the list
best.mod <- keep.best.model(listmod)
best.mod
## Model selected randomly among the model with the lowest AICc values
summary(best.mod)
predicted_hist(best.mod)


#### Example 2 - Multispecies analysis ####
# Example can be found in the supplementary material 2 - Bolduc et al 2017
# When the user has collected data for several species, it is possible to run the same detection models on each of them with the lsub and split arguments

set.seed(91) 
library(R2MCDS) 
### Import and filter data 
data(laridae) 
df2 <- mcds.filter(laridae,
                   transect.id = "WatchID",
                   distance.field = "Distance",
                   distance.labels = c("A", "B", "C", "D"),
                   distance.midpoints = c(25, 75, 150, 250),
                   effort.field = "WatchLenKm",
                   lat.field = "LatStart",
                   long.field = "LongStart",
                   sp.field = "Alpha",
                   date.field = "Date")
### Run analysis with the MCDS engine.
### Here, the 5-minute observation period (WatchID) is used as the sampling unit. 
mod1 <- mcds.wrap.point(df2,
                  SMP_EFFORT="WatchLenKm",
                  DISTANCE="Distance",
                  SIZE="Count",
                  units=list(Type="Point",
                             Distance="Radial",
                             Length_units="Kilometers",
                             Distance_units="Meters",
                             Area_units="Square kilometers"),
                  breaks=c(0,50,100,200,300),
                  #estimator=list(c("HN","CO")),
                  estimator = NULL,
                  multiplier = c(1, 0, 0),
                  factor = c("Observer"), # Covariate
                  monotone = "none", # Important for analysis with covariate
                  lsub=list(Alpha=c("BLKI","GBBG","HERG")), #list of all the species of interest to include in the analysis
                  split=TRUE, #indicates whether the species should be analyzed together split = FALSE
                  STR_LABEL="STR_LABEL",
                  STR_AREA="STR_AREA",
                  SMP_LABEL="WatchID",
                  path="C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/R_examples",
                  pathMCDS="C:/Program Files (x86)/Distance 7",
                  verbose= TRUE)

mod1

#look at the output for the Black-legged kittiwake
summary(mod1[[1]])

# plot the prediction for the Black-legged kittiwake
predicted_hist(mod1[[1]])

#### Example 3 - Rare species analysis ####
# Example can be found in the supplementary material 2 - Bolduc et al 2017
# possible to estimate the density and abundance of a rare species by incorporating the detection probability of a similar species into a multiplier

set.seed(91) 
library(R2MCDS) 
### Import and filter data 
data(laridae) 
df2 <- mcds.filter(laridae,
                   transect.id = "WatchID",
                   distance.field = "Distance",
                   distance.labels = c("A", "B", "C", "D"),
                   distance.midpoints = c(25, 75, 150, 250),
                   effort.field = "WatchLenKm",
                   lat.field = "LatStart",
                   long.field = "LongStart",
                   sp.field = "Alpha",
                   date.field = "Date")
### Run analysis with the MCDS engine. Here, the WatchID is used as the sample.
mod2 <- mcds.wrap.point(df2,
                  SMP_EFFORT="WatchLenKm",
                  DISTANCE="Distance",
                  SIZE="Count",
                  units=list(Type="Point",
                             Distance="Radial",
                             Length_units="Kilometers",
                             Distance_units="Meters",
                             Area_units="Square kilometers"),
                  breaks=c(0,50,100,200,300),
                  estimator=list(c("HN","CO")),
                  #multiplier = c(3, 0, 0),
                  lsub=list(Alpha=c("HERG")),
                  rare= list(Alpha=c("RBGU")),
                  split=TRUE,
                  STR_LABEL="STR_LABEL",
                  STR_AREA="STR_AREA",
                  SMP_LABEL="WatchID",
                  path="C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/R_examples",
                  pathMCDS="C:/Program Files (x86)/Distance 7",
                  verbose=FALSE)
mod2
summary(mod2)
predicted_hist(mod2)

#### Package examples ####
#### Example 1 ####
### Simple models without stratification
### Import and filter data
set.seed(91)
library(R2MCDS)
data(alcidae)
alcids <- mcds.filter(alcidae,
                      transect.id = "WatchID",
                      distance.field = "Distance",
                      distance.labels = c("A", "B", "C", "D"), 
                      distance.midpoints = c(25, 75, 150, 250),
                      effort.field = "WatchLenKm",
                      lat.field = "LatStart", 
                      long.field = "LongStart",
                      sp.field = "Alpha",
                      date.field = "Date") 

### Run analysis with the MCDS engine. Here, the WatchID is used as the sample.
dist.out1 <- mcds.wrap(alcids,
                       SMP_EFFORT="WatchLenKm",
                       DISTANCE="Distance",
                       SIZE="Count",
                       units=list(Type="Line",
                                  Distance="Perp",
                                  Length_units="Kilometers",
                                  Distance_units="Meters",
                                  Area_units="Square kilometers"),
                       breaks=c(0,50,100,200,300),
                       estimator=list(c("HN","CO")),
                       STR_LABEL="STR_LABEL",
                       STR_AREA="STR_AREA",
                       SMP_LABEL="WatchID", 
                       path="C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/R_examples",
                       pathMCDS="C:/Program Files (x86)/Distance 7",
                       verbose=FALSE)

summary(dist.out1)
#### Example 2 ####
### Run separate analysis for years 2007-2008
# To fixe the beug example
utils::View(alcids)
alcids$Year <- substr(alcids$Date, start = 1, stop = 4)
alcids$Year <- as.numeric(alcids$Year)
summary(alcids)


dist.out2 <- mcds.wrap(alcids,
                       SMP_EFFORT="WatchLenKm",
                       DISTANCE="Distance",
                       SIZE="Count",
                       units=list(Type="Line",
                                  Distance="Perp",
                                  Length_units="Kilometers",
                                  Distance_units="Meters",
                                  Area_units="Square kilometers"),
                       breaks=c(0,50,100,200,300),
                       estimator=list(c("HN","CO")),
                       lsub=list(Year=c(2007,2008)),
                       split=TRUE,
                       empty="Year",
                       STR_AREA="STR_AREA",
                       SMP_LABEL="WatchID", 
                       path="C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/R_examples",
                       pathMCDS="C:/Program Files (x86)/Distance 7",
                       verbose=FALSE)

### Get the names of the different models produced
names(dist.out2)
#####summary for the Year 2008 model
summary(dist.out2[["2008"]])
summary(dist.out2[["2007"]])
#### ***** type = "point"
# voir pour argument multiplier = 1 pour point d'écoute, multiplier = c(2, 0, 0) par défaur 
# voir pour Length/MEASURE (ligne juste éteinte pour le moment)
# Uniformiser les exemples entre packages et appendices du papier

# browser() permet de s'arreter dans le code de la fonction# Permet de rentrer dans l'environnement de la fonction


################## example with simulated data ########################
########################### Phase 1 ##################################
library(AHMbook)
ls("package:AHMbook")
?sim.pdata
###################
# Data simulation #
###################
set.seed(12345)
simu.data <- sim.pdata(N = 1000,
                       sigma = 1,
                       B = 3,
                       keep.all = TRUE,
                       show.plot = TRUE)

tmp <- sim.pdata(N = 1000,
                 sigma = 1,
                 keep.all = FALSE,
                 show.plot = FALSE)
summary(tmp)
str(tmp)

delta <- 0.5 # Width of distance bins
B <- 3 # Max count distance
dist.breaks <- seq(0, B, delta) # Make the interval cut points
dclass <- tmp$d %/% delta + 1 # %/% division entiere
nD <- length(dist.breaks) - 1 # How many intervals do you have
# y.obs <- table(dclass)
# y.padded <- rep(0, nD)
# names(y.padded) <- 1:nD
# y.padded[names(y.obs)] <- y.obs
# y.obs <- y.padded

###################
### DF building ###
###################

# Try to build a dataframe to MCDS analysis
transect.id <- "pioupiou"
distance.field <- dclass

effort.field <- 0# ? cause it's the length of transect or watch
lat.field <- 47.0
long.field <- -45.24
sp.field <- "Piou"
date.field <- "2019-03-15"

piou.data <- data.frame(transect.id = rep(transect.id, 163),
                        distance.field = distance.field,
                        effort.field = rep(1, 163),
                        lat.field = rep(lat.field, 163),
                        long.field = rep(long.field, 163),
                        sp.field = rep(sp.field, 163),
                        date.field = rep(date.field, 163),
                        Count = rep(1, 163))

piou.data$distance.field[piou.data$distance.field == 1] <- "A"
piou.data$distance.field[piou.data$distance.field == 2] <- "B"
piou.data$distance.field[piou.data$distance.field == 3] <- "C"
piou.data$distance.field[piou.data$distance.field == 4] <- "D"
piou.data$distance.field[piou.data$distance.field == 5] <- "E"
piou.data$distance.field[piou.data$distance.field == 6] <- "F"

piou.data$distance.field <- as.factor(piou.data$distance.field)

summary(piou.data)

###################
###### Model ######
###################

piou <-mcds.filter(piou.data,
                  transect.id = "transect.id",
                  distance.field = "distance.field",
                  distance.labels <- c("A", "B", "C", "D", "E", "F"),
                  distance.midpoints <- c(0.25, 0.75, 1.25, 1.75, 2.25, 2.75),
                  effort.field = "effort.field",
                  lat.field = "lat.field",
                  long.field = "long.field",
                  sp.field = "sp.field",
                  date.field = "date.field")
### Run analysis with the MCDS engine.
### Here, the 5-minute observation period (WatchID) is used as the sampling unit.
mod1 <- mcds.wrap.point(piou,
                        SMP_EFFORT="WatchLenKm",
                        DISTANCE="Distance",
                        SIZE="Count",
                        units=list(Type="Point",
                                   Distance="Radial",
                                   Length_units="Meters",
                                   Distance_units="Meters",
                                   Area_units="Square meters"),
                        breaks=c(0, 0.5, 1, 1.50, 2, 2.50, 3),
                        SMP_LABEL="WatchID",
                        STR_LABEL="STR_LABEL",
                        STR_AREA="STR_AREA",
                        estimator=list(c("HN","CO")),
                        multiplier = c(1, 0, 0),
                        path="C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/R_examples",
                        pathMCDS="C:/Program Files (x86)/Distance 7",
                        verbose=FALSE)

# File with errors = "log_xxx.tmp"
mod1
summary(mod1)
plot.distanceFit(mod1)

################## example with simulated data ########################
########################### Phase 2 ##################################
library(AHMbook)
ls("package:AHMbook")
?sim.pdata
#########################
# More data simulation  #
########################
ll <- list()
j <- 1:100
#set.seed(64)
al <- sample(300:1000, max(j))

for (i in j){
  # set.seed(i)
  simu.data <- sim.pdata(N = al[i],
                         sigma = 1,
                         B = 3,
                         keep.all = TRUE,
                         show.plot = TRUE)
  print(simu.data$N.real)
  tmp <- sim.pdata(N = al[i],
                   sigma = 1,
                   keep.all = FALSE,
                   show.plot = FALSE)
  ll[[i]] <- tmp
}


delta <- 0.5 # Width of distance bins
B <- 3 # Max count distance
dist.breaks <- seq(0, B, delta) # Make the interval cut points
dclass <- lapply(ll, function(i){
  i$d %/% delta + 1
})
nD <- length(dist.breaks) - 1 # How many intervals do you have

###################
### DF building ###
###################

# Try to build a dataframe to MCDS analysis
piou.data <- data.frame()
for(i in 1:length(ll)){
  transect.id <- rep(paste("piou", i, sep = ""), length(dclass[[i]]))
  d.field <- dclass[[i]]
  
  effort.field <- rep(1, length(dclass[[i]])) # ? cause it's the length of transect or watch
  lat.field <- rep(47.0, length(dclass[[i]]))
  long.field <- rep((-45.24 + i), length(dclass[[i]]))
  sp.field <- rep("Piou", length(dclass[[i]]))
  date.field <- rep("2019-03-15", length(dclass[[i]])) 
  Count <- rep(1, length(dclass[[i]]))
  real.abun <- rep(ll[[i]]$N.real, length(dclass[[i]]))
  
  dt <- data.frame(transect.id ,
                          d.field,
                          effort.field,
                          lat.field,
                          long.field,
                          sp.field,
                          date.field,
                          Count,
                   real.abun)
  piou.data <- rbind(piou.data, dt)
  
}

# Convert numerical distance value to categorical with letters
piou.data$distance.field <- LETTERS[piou.data$d.field]
piou.data$distance.field <- as.factor(piou.data$distance.field)

summary(piou.data)
utils::View(piou.data)

# Real total abundance

real.N <- sum(tapply(piou.data$real.abun, piou.data$transect.id, unique))

real.per.meter <- real.N / (3.141593 * B^2 * length(ll))

###################
###### Model ######
###################
library(R2MCDS)
piou <-mcds.filter(piou.data,
                   transect.id = "transect.id",
                   distance.field = "distance.field",
                   distance.labels <- c("A", "B", "C", "D", "E", "F"),
                   distance.midpoints <- c(0.25, 0.75, 1.25, 1.75, 2.25, 2.75),
                   effort.field = "effort.field",
                   lat.field = "lat.field",
                   long.field = "long.field",
                   sp.field = "sp.field",
                   date.field = "date.field")
### Run analysis with the MCDS engine.
### Here, the 5-minute observation period (WatchID) is used as the sampling unit.
mod.piou <- mcds.wrap.point(piou,
                        SMP_EFFORT="WatchLenKm",
                        DISTANCE="Distance",
                        SIZE="Count",
                        units=list(Type="Point",
                                   Distance="Radial",
                                   Length_units="Meters",
                                   Distance_units="Meters",
                                   Area_units="Square meters"),
                        breaks=c(0, 0.5, 1, 1.50, 2, 2.50, 3),
                        SMP_LABEL="WatchID",
                        STR_LABEL="STR_LABEL",
                        STR_AREA="STR_AREA",
                        estimator=list(c("HN","CO")),
                        #estimator = NULL,
                        multiplier = c(1, 0, 0),
                        path="C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/R_examples",
                        pathMCDS="C:/Program Files (x86)/Distance 7",
                        verbose=FALSE)

# File with errors = "log_xxx.tmp"
mod.piou
summary(mod.piou)
real.per.meter
plot.distanceFit(mod.piou)

################################### example with Distance data ##############################
################### Project 1 - Amakihi analysis -- > Simple test model #####################
###################### *** Correct example **** ###########################
library(R2MCDS)
setwd("C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/DISTANCE 7.2")
list.files()

amak <- read.table("Amakihi.txt", h = T, sep = "\t")
head(amak)
summary(amak)
amak <- na.omit(amak)

amak$sp <- "Honeycreeper"
amak$lat <- 43 
amak$long <- 90

# Conversion of the metric distance to letter categories
amak$Distance_L[amak$Radial_distance < 50] <- "A" 
amak$Distance_L[amak$Radial_distance >= 50 & amak$Radial_distance < 100] <- "B" 
amak$Distance_L[amak$Radial_distance >= 100 & amak$Radial_distance < 150] <- "C"
amak$Distance_L[amak$Radial_distance >= 150 & amak$Radial_distance < 200] <- "D"
amak$Distance_L[amak$Radial_distance >= 200 & amak$Radial_distance < 250] <- "E"
amak$Distance_L[amak$Radial_distance >= 250 & amak$Radial_distance < 300] <- "F"

amak$Distance_L <- as.factor(amak$Distance_L)

summary(amak)
utils::View(amak)

amak.1 <- mcds.filter(amak,
                    transect.id = "point_ID",
                    distance.field = "Distance_L", #Use the distance variable with letters here
                    distance.labels <- c("A", "B", "C", "D", "E", "F"),
                    distance.midpoints <- c(25, 75, 125, 175, 225, 275),
                    effort.field = "Survey_effort",
                    lat.field = "lat",
                    long.field = "long",
                    sp.field = "sp",
                    date.field = "date"
                    ) 

mod1 <- mcds.wrap(amak,
                  path="C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/R_examples",
                  pathMCDS="C:/Program Files (x86)/Distance 7",
                  STR_LABEL="date",
                  STR_AREA="region_ID",
                  SMP_LABEL="point_ID",
                  units=list(Type="Point",
                             Distance="Radial",
                             Length_units="Meters",
                             Distance_units="Meters",
                             Area_units="Hectares"),
                  SMP_EFFORT="Survey_effort",
                  DISTANCE="Radial_distance",
                  SIZE="Count",
                  breaks=seq(0, 250, 50),
                  #covariates = c("obs", "MAS", "HAS"),
                  #factor = "a factor giving the name of factors to be used in analysis",
                  #lsub = list(unique(amak$date)),
                  #stratum = TRUE,
                  #split = TRUE,
                  #rare = "used when a species has few observations to estimate a detection function",
                  #period = "vector of characters of length 2 containing the extreme dates for which the analysis should be restricted",
                  #detection = "see help file",
                  #monotone = "Strict",
                  estimator=list(c("HN","CO")),
                  multiplier = c(1, 0, 0),
                  #empty = "determine how empty transects are to be selected in the analysis",
                  verbose=FALSE)

# File with errors = "log_xxx.tmp"
summary(mod1)
x11()
plot.distanceFit(mod1)

################################### example with Distance data ##############################
################### Project 2 - Transect Line example analysis #####################
######################## *** Correct example *** ##################################
library(R2MCDS)
setwd("C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/DISTANCE 7.2")
list.files()

line <- read.table("Line_Transect.txt", h = T, sep = "\t", dec = ",")
head(line)
summary(line)
line <- na.omit(line)

line$sp <- "sp"
line$lat <- 43 
line$long <- 90

# Conversion of the metric distance to letter categories
line$Distance_L[line$Perp_distance < 5] <- "A" 
line$Distance_L[line$Perp_distance >= 5 & line$Perp_distance < 10] <- "B" 
line$Distance_L[line$Perp_distance >= 10 & line$Perp_distance < 15] <- "C"
line$Distance_L[line$Perp_distance >= 15 & line$Perp_distance < 20] <- "D"
line$Distance_L[line$Perp_distance >= 20 & line$Perp_distance < 25] <- "E"
line$Distance_L[line$Perp_distance >= 25 & line$Perp_distance < 30] <- "F"
line$Distance_L[line$Perp_distance >= 30 & line$Perp_distance < 35] <- "G"
line$Distance_L[line$Perp_distance >= 35 & line$Perp_distance < 40] <- "H"

line$Distance_L <- as.factor(line$Distance_L)

line$date <- "2005-06-07"
line$Count <- 1

summary(line)
utils::View(line)

line.1 <- mcds.filter(line,
                      transect.id = "ID_transect",
                      distance.field = "Distance_L", #Use the distance variable with letters here
                      distance.labels <- c("A", "B", "C", "D", "E", "F", "G", "H"),
                      distance.midpoints <- c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5),
                      effort.field = "Line_length",
                      lat.field = "lat",
                      long.field = "long",
                      sp.field = "sp",
                      date.field = "date"
) 

mod.line <- mcds.wrap(line,
                  path="C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/R_examples",
                  pathMCDS="C:/Program Files (x86)/Distance 7",
                  STR_LABEL= "Label",
                  STR_AREA="Area",
                  SMP_LABEL="ID_transect",
                  units=list(Type="Line",
                             Distance="Perp",
                             Length_units="Kilometers",
                             Distance_units="Meters",
                             Area_units="Square kilometers"),
                  SMP_EFFORT="Line_length",
                  DISTANCE="Perp_distance", # Need to use numerical value if original database
                  SIZE="Cluster_size",
                  breaks=seq(0, 40, 5),
                  #covariates = NA,
                  #factor = "a factor giving the name of factors to be used in analysis",
                  #lsub = NA,
                  #stratum = TRUE,
                  #split = TRUE,
                  #rare = "used when a species has few observations to estimate a detection function",
                  #period = "vector of characters of length 2 containing the extreme dates for which the analysis should be restricted",
                  #detection = "see help file",
                  #monotone = "Strict",
                  estimator=list(c("HN","CO")),
                  multiplier = c(1, 0, 0),
                  #empty = "determine how empty transects are to be selected in the analysis",
                  verbose=FALSE)
summary(mod.line)
#x11()
plot.distanceFit(mod.line)


mod.line.wrap <- mcds.wrap(line.1,
                      path="C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/R_examples",
                      pathMCDS="C:/Program Files (x86)/Distance 7",
                      STR_LABEL= "Label",
                      STR_AREA="Area",
                      SMP_LABEL="WatchID",
                      units=list(Type="Line",
                                 Distance="Perp",
                                 Length_units="Kilometers",
                                 Distance_units="Meters",
                                 Area_units="Square kilometers"),
                      SMP_EFFORT="WatchLenKm",
                      DISTANCE="Distance",
                      SIZE="Cluster_size",
                      breaks=seq(0, 40, 5),
                      #covariates = NA,
                      #factor = "a factor giving the name of factors to be used in analysis",
                      #lsub = NA,
                      #stratum = TRUE,
                      #split = TRUE,
                      #rare = "used when a species has few observations to estimate a detection function",
                      #period = "vector of characters of length 2 containing the extreme dates for which the analysis should be restricted",
                      #detection = "see help file",
                      #monotone = "Strict",
                      estimator=list(c("HN","CO")),
                      multiplier = c(1, 0, 0),
                      #empty = "determine how empty transects are to be selected in the analysis",
                      verbose=FALSE)

summary(mod.line.wrap)
#x11()
plot.distanceFit(mod.line.wrap)
