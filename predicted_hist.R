#' Title Plot a detection histogram of apparent density and the predicted detction curve for the model.
#'
#' @param model A distanceFit object
#'
#'@details
#' Make a histogram of the observations and the predicted detection curve in function of the distance class 
#'@section Author:Christian Roy
#'@examples
#'data(alcidae)
#'alcids <- distance.filter(alcidae, transect.id = "WatchID", distance.field = "Distance", distance.labels = c("A", "B", "C", "D"), 
#'                          distance.midpoints = c(25, 75, 150, 250), effort.field = "WatchLenKm", lat.field = "LatStart", 
#'                          long.field = "LongStart", sp.field = "Alpha", date.field = "Date") 
#'
#'### Run analysis with the MCDS engine. Here, the WatchID is used as the sample.
#'dist.out1 <-distance.wrap(alcids, SMP_EFFORT="WatchLenKm",DISTANCE="Distance",SIZE="Count",
#'                          units=list(Type="Line",Distance="Perp",Length_units="Kilometers",
#'                                     Distance_units="Meters",Area_units="Square kilometers"),
#'                          breaks=c(0,50,100,200,300), estimator=list(c("HN","CO")),
#'                          STR_LABEL="STR_LABEL", STR_AREA="STR_AREA",SMP_LABEL="WatchID", 
#'                          path="c:/temp/distance",
#'                          pathMCDS="C:/Distance 6",verbose=FALSE)
#'
#'summary(dist.out1)
#'predicted_hist(x)
#'
#' #END
predicted_hist <- function(model){
  p <- observation_hist(model[['input_data']][['observations']], count='SIZE', dist.class='DISTANCE',
                        keep.class=as.character(unique(sort(as.numeric(model[['input_data']][['observations']]$DISTANCE)))),
                        breaks=model[['input_data']][['breaks']], color='powderblue', ungroup=F,
                        rescale=model$detection[['Global']][,'predicted'][1]) +
    labs(x = 'Distance', y = 'Detection probability')
  
  pred.df <- data.frame(x=model$detection[['Global']][,'distance'],y=model$detection[['Global']][,'predicted'])
  p + geom_line(data=pred.df, aes(x=x,y=y), linetype=1, size=1.25)
  
}

#' Title Plot a detection histogram of apparent density and the predicted detction curve for the model.
#'
#' @param model A distanceFit object
#'
#'@details
#' Make a histogram of the observations and the predicted detection curve in function of the distance class 
#'@section Author:Christian Roy
#'@examples
#'data(alcidae)
#'alcids <- mcds.filter(alcidae, transect.id = "WatchID", distance.field = "Distance", distance.labels = c("A", "B", "C", "D"), 
#'                          distance.midpoints = c(25, 75, 150, 250), effort.field = "WatchLenKm", lat.field = "LatStart", 
#'                          long.field = "LongStart", sp.field = "Alpha", date.field = "Date") 
#'
#'### Run analysis with the MCDS engine. Here, the WatchID is used as the sample.
#'dist.out1 <- mcds.wrap(alcids, SMP_EFFORT="WatchLenKm",DISTANCE="Distance",SIZE="Count",
#'                          units=list(Type="Line",Distance="Perp",Length_units="Kilometers",
#'                                     Distance_units="Meters",Area_units="Square kilometers"),
#'                          breaks=c(0,50,100,200,300), estimator=list(c("HN","CO")),
#'                          STR_LABEL="STR_LABEL", STR_AREA="STR_AREA",SMP_LABEL="WatchID", 
#'                          path="c:/temp/distance",
#'                          pathMCDS="C:/Distance 6",verbose=FALSE)
#'
#'summary(dist.out1)
#'plot(dist.out1)
#'
#' #END
plot.distanceFit <- function(model){
  p <- observation_hist(model[['input_data']][['observations']], count='SIZE', dist.class='DISTANCE',
                        keep.class=as.character(unique(sort(as.numeric(model[['input_data']][['observations']]$DISTANCE)))),
                        breaks=model[['input_data']][['breaks']], color='powderblue', ungroup=F,
                        rescale=model$detection[['Global']][,'predicted'][1]) +
    labs(x = 'Distance', y = 'Detection probability')
  
  pred.df <- data.frame(x=model$detection[['Global']][,'distance'],y=model$detection[['Global']][,'predicted'])
  p + geom_line(data=pred.df, aes(x=x,y=y), linetype=1, size=1.25)
  
}