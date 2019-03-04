####### TEST zone #####

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


mod1


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
                   date.field = "Date") ### Run analysis with the MCDS engine.

### Here, the 5-minute observation period (WatchID) is used as the sampling unit. 
mod1 <- mcds.wrap.point(df2,
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
                  factor = c("Observer"),
                  monotone = "none",
                  lsub=list(Alpha=c("BLKI","GBBG","HERG")),
                  split=TRUE,
                  STR_LABEL="STR_LABEL",
                  STR_AREA="STR_AREA",
                  SMP_LABEL="WatchID",
                  path="c:/temp/distance",
                  pathMCDS="C:/Program Files (x86)/Distance 7",
                  verbose= TRUE)

#### ***** type = "point"
# voir pour argument multiplier = 1 pour point d'écoute, multiplier = c(2, 0, 0) par défaur 
# voir pour Length/MEASURE (ligne juste éteinte pour le moment)
# Uniformiser les exemples entre packages et appendices du papier

# browser() permet de s'arreter dans le code de la fonction# Permet de rentrer dans l'environnement de la fonction