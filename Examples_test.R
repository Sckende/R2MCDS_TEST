####### TEST zone #####


#### Example 1 - Simple analysis ####
# Example can be found in the supplementary material 1 - Bolduc et al 2017
set.seed(91)
library(R2MCDS)
### Import and filter data
data(alcidae)
df1 <-mcds.filter(alcidae, transect.id = "WatchID",
                  distance.field = "Distance",
                  distance.labels = c("A", "B", "C", "D"),
                  distance.midpoints = c(25, 75, 150, 250),
                  effort.field = "WatchLenKm",
                  lat.field = "LatStart", long.field = "LongStart",
                  sp.field = "Alpha", date.field = "Date")
### Run analysis with the MCDS engine.
### Here, the 5-minute observation period (WatchID) is used as the sampling unit.
mod1 <- mcds.wrap.point(df1,
                  SMP_EFFORT="WatchLenKm",
                  DISTANCE="Distance",
                  SIZE="Count",
                  units=list(Type="Point",
                             Distance="Radial",
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
predicted_hist(mod1)

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
                     SMP_LABEL="WatchID",
                     path="C:/Users/HP_9470m/OneDrive - Université de Moncton/GC job - R2MCDS/R_examples",
                     pathMCDS="C:/Program Files (x86)/Distance 7",
                     verbose=FALSE)
listmod
summary(listmod[[1]])
lapply(listmod, summary)

predicted_hist(listmod[[1]])
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
                  estimator=list(c("HN","CO")),
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



#### ***** type = "point"
# voir pour argument multiplier = 1 pour point d'écoute, multiplier = c(2, 0, 0) par défaur 
# voir pour Length/MEASURE (ligne juste éteinte pour le moment)
# Uniformiser les exemples entre packages et appendices du papier

# browser() permet de s'arreter dans le code de la fonction# Permet de rentrer dans l'environnement de la fonction