##############################################################################
##############################################################################
############## See MCDS Command Language ####################################
################Distance manual page 331 #####################################
##############################################################################

#' @export
#'@title Fit the distance sampling model provided by the Distance software.
#'
#'
#'@description This function fits the distance sampling model provided by the Distance program using its MCDS engine.

#'@param dataset A \code{\link{data.frame}} with observations.

#'@param path The path where to store the input and output files of the MCDS engine.

#'@param pathMCDS The path where the MCDS engine is located.

#'@param STR_LABEL Name of the column to use as a stratum label.

#'@param STR_AREA Name of the column to use for the stratum area.

#'@param SMP_LABEL Name of the column to use for the transect/watch label.

#'@param units list of th units used for the analysis. Contains the Type of analysis (Line or Point), the Distance engine to use (Perp or Radial), 
#'the Length units, the Distance units, and the Area_units.For the possible units of distance and area see Distance 6.2 documentation.     

#'@param SMP_EFFORT Length in of the transect or the transect/watch unit.

#'@param DISTANCE Distance of the observation.

#'@param SIZE Number of individuals in the observation.

#'@param breaks A vector giving the distance intervals in meters to be used in the analysis.

#'@param covariates A vector giving the name of covariates.

#'@param factor A vector giving the name of factors to be used in the analysis.

#'@param lsub A named list giving the subsets to be used in the analysis. The names of the list are the names of columns used to subset the dataset.
#'Each element of the list is a vector giving the values to keep in the analysis for a given column. When a vector is \code{NULL},
#'all values are kept for this column. Default value \code{lsub = NULL}. See examples for further details.

#'@param stratum When \code{stratum = "STR_LABEL"}, the model will be stratified with the label.
#'When \code{stratum} is the name of a column in the data, the model will be
#'post-stratified according to this column (see Thomas et al. 2010). Default value
#'is \code{NULL}.

#'@param split When \code{split = TRUE}, separate models are ran according to the subsets determined by \code{lsub}. Default value \code{FALSE}.

#'@param rare This argument is used when a species has few observations to estimate a detection function.
#'The probability of detection is estimated from a group of similar species and a multiplier. A named list
#'has to be given, with the column name where species are stored as the name of the list and the name of the species as
#'an element (ex: \code{list(Species = "HERG")}). Default value \code{NULL}. See details.

#'@param period A vector of characters of length 2 containing the extreme dates for which the analysis
#'should be restricted. Dates have to be in the "yyyy-mm-dd" format.

#'@param detection Currently, set to \code{detection = "All"}. Can also be \code{detection = "Stratum"}.

#'@param monotone Currently, set to \code{monotone = "Strict"}.

#'@param estimator When set to \code{NULL}, the following key functions and expansion terms will be used: UN-CO, UN-PO, HN-CO, HN-HE, HA-CO and HA-PO.
#'If the user wants to choose the key functions and the expansion terms used, a list has to be given with each element a vector of
#'length 2 with the first element the key function and the second element the expansion term (ex: \code{list(c("HN","CO"),c("HA","PO")}).
#'When \code{rare != NULL}, only one set is used (UN-CO) for the final specific model. 

#'@param multiplier Value by which the estimates of density or abundance are multiplied. The first value is the multiplier, the second the SE and the third the degree of freedom associated with the multiplier (useful when using the probability of detection as a multiplier in a two-step analyses with the second step). Default value \code{multiplier = c(2,0,0)} meaning
#'only one-half of the transect is surveyed and the value is known with certainty with an infinite degree of freedom. When \code{rare != NULL}, the multiplier will be modified to account
#'for the probability of detection.

#'@param empty Determine how empty transects are to be selected in the analysis. When \code{empty = NULL}, all
#'empty transects are included for every element of \code{lsub}. For example, when models are splitted according
#'to species, empty transects and transects where a species was not detected need to be considered in the analysis for that species.
#'When \code{lsub} contains geographic or temporal elements, empty transects need to be restricted to the subsets given. In this case,
#'a vector of character has to be given with the names in \code{lsub} for which the empty transects are to be restricted. When \code{split = TRUE}
#'and \code{empty != NULL}, empty transects will be splitted according to the names in \code{empty}. In any case, it is assumed that
#'the \code{dataset} contains at least a line for every transect executed, either with or without an observation. See examples for
#'further details.

#'@param verbose When set to \code{TRUE}, prints the input file given to the MCDS engine. Default value \code{FALSE}.

#'@details
#'Uses the MCDS engine from program Distance. The function produces an input file and submits it to the MCDS engine
#'through the \code{\link{system}} function. The results are then extracted from the output files and returned as a
#'list object.

#'@return
#'An object of class \code{"distanceFit"}, when \code{split = FALSE}. When \code{split = TRUE} and \code{lsub != NULL},
#'a named list with the different models of class \code{"distanceFit"}.

#'Each object of class \code{"distanceFit"} is a named list with components \code{model_fitting}, \code{parameter_estimates},
#'\code{chi_square_test}, \code{density_estimate}, \code{detection} and \code{path}. Elements
#'of the list are accessible through the \code{$} operator. For each component
#'except \code{detection} and \code{path}, a \code{\link{list}} of length 2 is given with
#'component "\code{Global}" and "\code{Stratum}". Depending on the analysis chosen,
#'one of these component will be empty (\code{NULL}).

#'@references
#'Thomas, L., S.T. Buckland, E.A. Rexstad, J. L. Laake, S. Strindberg, S. L. Hedley, J. R.B. Bishop,
#'T. A. Marques, and K. P. Burnham. 2010. \emph{Distance software: design and analysis of distance sampling
#'surveys for estimating population size.} Journal of Applied Ecology 47: 5-14.
#'DOI: 10.1111/j.1365-2664.2009.01737.x

#'@section Author:Francois Rousseu, Christian Roy, Francois Bolduc

#'@examples
#'########################################
#'### Simple models without stratification
#'### Import and filter data
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
#'### Run separate analysis for years 2008-2009
#'dist.out2 <- mcds.wrap(alcids, SMP_EFFORT="WatchLenKm",DISTANCE="Distance",SIZE="Count",
#'                           units=list(Type="Line",Distance="Perp",Length_units="Kilometers",
#'                                      Distance_units="Meters",Area_units="Square kilometers"),
#'                           breaks=c(0,50,100,200,300), estimator=list(c("HN","CO")),
#'                           lsub=list(Year=c(2008,2009)), split=TRUE, empty="Year",
#'                           STR_AREA="STR_AREA",SMP_LABEL="WatchID", 
#'                           path="c:/temp/distance",
#'                          pathMCDS="C:/Distance 6",verbose=FALSE)
#'
#'### Get the names of the different models produced
#'names(dist.out2)
#'#####summary for the Year 2008 model
#'summary(dist.out2[["2008"]])
#'#END


mcds.wrap.point <-
  function(dataset,
           path,
           pathMCDS,
           SMP_LABEL="SMP_LABEL",
           SMP_EFFORT="SMP_EFFORT",
           DISTANCE="DISTANCE",
           SIZE="SIZE",                												
           STR_LABEL="STR_LABEL",
           STR_AREA="STR_AREA",												
           units=list(Type="Line",Distance="Perp",Length_units="Kilometer",
                      Distance_units="Meter",Area_units="Square kilometer"),
           breaks=c(0,50,100,200,300),
           covariates=NULL,
           factor=NULL,
           lsub=NULL,
           stratum=NULL, 
           split=TRUE,
           rare=NULL, 
           period=NULL,
           detection="All", 
           monotone="Strict",
           estimator=NULL,
           #multiplier=c(2,0,0), 
           multiplier=NULL,
           empty=NULL,
           verbose=FALSE
  ){
    
    ##Make sure input are correct
    #check for detection options
    if(toupper(detection)%in%c("ALL","STRATUM")==FALSE)
      stop("Detection options are 'All' or 'Stratum'")
    
    if(toupper(monotone)%in%c("NONE","WEAK", "STRICT")==FALSE)
      stop("monotone options are 'None', 'Weak' or 'Strict'")
    
    #Check for monotone options
    if(!is.null(covariates) & monotone!="none")
      stop("monotone must be set to 'none' when  factors or covariates are included")
    
    if(!is.null(factor) & monotone!="none")
      stop("monotone must be set to 'none' when  factors or covariates are included")
    
    #Check for units
    if(length(units)!= 5)
      stop("units list must includes values for Type, Distance, Length_units, Distance_units and Area_units")
    
    if(toupper(units$Type)%in%c("LINE","POINT", "CUE")==FALSE)
      stop("Type of sampling available are: Line, Point or Cue")
    
    if(toupper(units$Distance)%in%c("PERP","RADIAL")==FALSE)
      stop("Type of distance available are: Perp or Radial")
    
    # if(toupper(units$Type)%in%c("POINT", "CUE") & toupper(units$Distance)!=c("RADIAL"))
    #   stop("For Points and Cue sampling scheme Distance must be set to radial")
    
    # Automatic Distance value when units$Type = Point or Cue
    if(toupper(units$Type)%in%c("POINT", "CUE") & toupper(units$Distance)!=c("RADIAL")) {
      units$Distance = c("Radial")
      warning("Distance value has been replaced from Perp to Radial")
      }
    
    if(toupper(units$Length_units)%in%c("CENTIMETERS", "METERS", "KILOMETERS", "MILES", 
                                        "INCHES", "FEET", "YARDS", "NAUTICAL MILES")==FALSE)
      stop("Distance units must be one of thse: 'Centimeters', 'Meters', 'Kilometers', 'Miles', 'Inches', 
           'Feet', 'Yards', Or 'Nautical Miles'")  
    
    if(toupper(units$Distance_units)%in%c("CENTIMETERS", "METERS", "KILOMETERS", "MILES", 
                                          "INCHES", "FEET", "YARDS", "NAUTICAL MILES")==FALSE)
      stop("Distance units must be one of thse: 'Centimeters', 'Meters', 'Kilometers', 'Miles', 'Inches', 
           'Feet', 'Yards', Or 'Nautical Miles'")  
    
    if(toupper(units$Area_units)%in%c("SQUARE CENTIMETERS", "SQUARE METERS", "SQUARE KILOMETERS", "SQUARE MILES", 
                                      "SQUARE INCHES", "SQUARE FEET", "SQUARE YARDS", "HECTARES")==FALSE)
      stop("Distance units must be one of thse: 'Centimeters', 'Meters', 'Kilometers', 'Miles', 'Inches', 
           'Feet', 'Yards', Or 'Nautical Miles'")  
    
    #get the list of arguments
    arguments<-as.list(environment())
    if(!is.null(period)){
      dataset<-dataset[dataset$Date>=min(period) & dataset$Date<=max(period),]
    }
    
    #Make sure the directory exist to save output
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    
    
    #Make the factor characters for the function
    if(!is.null(factor)){
      for(i in 1:length(factor)){
        dataset[,factor[i]] <- as.character(dataset[,factor[i]])
      }   
    }
    
    # Automatic values for multiplier argument depending on the type argument
    if(is.null(multiplier)){
      if(units$Type == "Line"){
        multiplier <- c(2, 0, 0)
      }else{
          multiplier <- c(1, 0, 0)
        }
      
    }
    
    # this recursively fits a model to get an estimate of p for a rarer species
    if(!is.null(rare)){		
      if(!names(rare)%in%names(dataset)){stop(paste("The column",names(rare)[1],"not in dataset"))}
      if(length(rare)>1 || length(rare[[1]])>1){stop("Only a list and an element of length 1 is accepted for argument rare")}
      
      #arguments[["dataset"]]<-dataset[dataset[,names(rare)[1]]!=rare[[1]],]
      arguments[["dataset"]]<-dataset
      arguments[["lsub"]]<-lsub # here I could use the lsub argument instead of subsetting in the previous line
      arguments[["split"]]<-FALSE
      arguments[["stratum"]]<-NULL
      arguments[["detection"]]<-"All"
      arguments[["rare"]]<-NULL
      arguments[["verbose"]]<-FALSE
     
       premod<-do.call("mcds.wrap.point",arguments)
      if(class(premod)=="distanceFit"){ # distanceFit = output for a single model. If more than one = distanceList
        best.premod <- premod
      }else{
        best.premod <- keep.best.model(premod)
      }      
      detect_table<-best.premod [["parameter_estimates"]][["Global"]]
      pdetect<-detect_table[detect_table[,"Parameters"]=="p","Estimates"]
      sedetect<- detect_table[detect_table[,"Parameters"]=="p","SE"]
      dfdetect<- detect_table[detect_table[,"Parameters"]=="p","df"]
      warning(paste("Rare species detecion modified by a factor of ", round(1/pdetect,3), " with a ", round(sedetect/pdetect^2,3), " standard error and ",dfdetect," degrees of freedom based on estimation of the global model", sep=""), call. = FALSE)
    }
    
    if(!any("STR_AREA"==names(dataset)) && STR_AREA=="STR_AREA"){
      dataset<-cbind(STR_AREA=rep(1,nrow(dataset)),dataset,stringsAsFactors=F)
    }	
    if(!any("STR_LABEL"==names(dataset)) && STR_LABEL=="STR_LABEL"){
      dataset<-cbind(STR_LABEL=rep(1,nrow(dataset)),dataset,stringsAsFactors=F)
    }	
    if(!is.null(stratum)){
      if(STR_LABEL!="STR_LABEL"){
        dataset[,"STR_LABEL"]<-dataset[,STR_LABEL]	
      }
      dataset[,"STR_LABEL"]<-dataset[,stratum] # stratum as predominance over STR_LABEL
      STR_LABEL<-stratum	
    }
    if(!is.null(stratum) && stratum=="STR_LABEL"){
      if(is.null(STR_AREA)){
        stop("No area given (STR_AREA = NULL) and stratum = \"STR_LABEL\"")
      }else{	
        dataset[,"STR_AREA"]<-dataset[,STR_AREA]
      }
    }
    
    if(!is.null(empty) && !is.null(lsub)){
      if(!all(empty%in%names(lsub))){stop("Names in empty do not correspond to names in lsub")}
    }
    transects<-dataset
    if(!is.null(rare)){
      dataset<-dataset[dataset[,names(rare)[1]]==rare[[1]],]
    }else{
      if(!is.null(lsub)){  
        for(i in 1:length(lsub)){
          if(is.null(lsub[[i]])){next}
          dataset<-dataset[dataset[,names(lsub)[i]]%in%lsub[[i]],]			 
          if(names(lsub)[i]%in%empty){
            transects<-transects[transects[,names(lsub)[i]]%in%lsub[[i]],]
          }			
        }
      }
    }
    transects[,DISTANCE]<-""
    
    ######################## add-on
    if(!is.null(lsub) && split){
      if(!is.null(empty)){			
        transects<-dlply(transects,empty)
        dataset<-dlply(dataset,empty)			
        transects<-transects[names(dataset)]
        l<-1:length(dataset)
        names(l)<-names(transects)
        dataset<-llply(l,function(i){
          ans<-dlply(dataset[[i]],setdiff(names(lsub),empty))
          empty.transects<-transects[[i]]
          n<-names(ans)
          ans<-lapply(ans,function(j){rbind(j,empty.transects)})
          names(ans)<-n
          ans
        })
        if(all(names(lsub)%in%empty)){
          name<-names(dataset)
          dataset<-unlist(dataset,recursive=FALSE,use.names=FALSE)
          names(dataset)<-name
        }else{
          dataset<-unlist(dataset,recursive=FALSE,use.names=TRUE)
        }
      }else{
        dataset<-dlply(dataset,names(lsub))
        dataset<-llply(dataset,function(i){rbind(i,transects)})
      }
    }else{		
      dataset<-list(rbind(dataset,transects))
    }
    #########################
    
    dataset<-lapply(dataset,function(i){
      
      res<-i[!is.na(i[,DISTANCE]),] #add1
      good.label<-unique(res[res[,DISTANCE]!="",SMP_LABEL])#add2
      res<-res[(res[,DISTANCE]!="") | (res[,DISTANCE]=="" & !res[,SMP_LABEL]%in%good.label & !duplicated(res[,SMP_LABEL])),]#add3
      res[,DISTANCE]<-ifelse(is.na(res[,DISTANCE]),"",res[,DISTANCE])
      res[,SIZE]<-ifelse(res[,DISTANCE]=="","",res[,SIZE])
      if(all(res[,DISTANCE]=="")){
        res<-NULL
      }
      res
    })
    dataset<-dataset[!sapply(dataset,is.null)] #removes emptys transects when they are marked with "" for the species name
    
    if(!is.null(names(dataset))){
      names(dataset)<-gsub("\\.","\\_",names(dataset))
      names(dataset)<-gsub("'","",names(dataset))	# this was put for gyre d'anticosti #########
    }
    
    dat.file<-paste(paste("data",names(dataset),sep="_"),".tmp",sep="")
    
    #browser()
    lans<-vector(mode="list",length=length(dataset))
    names(lans)<-names(dataset)
    
    notrun<-NULL
    #browser()
    for(i in seq_along(dataset)){
      dat<-dataset[[i]]
      dat<-dat[!(is.na(dat[,DISTANCE]) | is.na(dat[,SIZE])),]
      #dat<-dat[!(duplicated(dat[,SMP_LABEL]) & dat[,DISTANCE]==""),]
      #######################################################
      ### input
      
      opts2<-list()
      opts2["Options;"]<-""
      opts2["Type="] <- paste(units$Type,";",sep="")
      
      if(units$Type=="Line"){
      opts2["Length /Measure="] <- paste("'",units$Length_units,"';",sep="")
      }
      
      if(units$Distance=="Perp"){
        opts2["Distance=Perp /Measure="] <- paste("'",units$Distance_units,"';",sep="")    
      }else{
        opts2["Distance=Radial /Measure="] <- paste("'",units$Distance_units,"';",sep="")    
      }
      opts2["Area /Units="] <- paste("'",units$Area_units,"';",sep="")
      opts2["Object="]<-"Cluster;"
      opts2["SF="]<-"1;"
      opts2["Selection="]<-"Specify;"
      if(is.null(rare)){
        opts2["Selection="]<-"Sequential;"
        opts2["Lookahead="]<-"1;"
        opts2["Maxterms="]<-"5;"
      }
      opts2["Confidence="]<-"95;"
      opts2["Print="]<-"Selection;"
      opts2["End1;"]<-""
      opts2["Data /Structure="]<-"Flat;"
      dat<-dat[,unique(c(STR_LABEL,STR_AREA,SMP_LABEL,SMP_EFFORT,DISTANCE,SIZE,factor,covariates))]
      names(dat)[1:6]<-c("STR_LABEL","STR_AREA","SMP_LABEL","SMP_EFFORT","DISTANCE","SIZE")
      
      opts2["Fields="]<-paste(paste(names(dat),collapse=", "),";",sep="")
      if(!is.null(factor)){
        for(j in 1:length(factor)){
          labels<-sort(unique(dat[,factor[j]]))
          labels<-labels[!is.na(labels)]
          opts2[[length(opts2)+1]] <- paste(paste(paste(c(" /Name="," /Levels="," /Labels="),c(factor[j],length(labels),paste(labels,collapse=",")),sep=""),collapse=""),";",sep="") 
        }
        names(opts2)[(which(names(opts2)=="Fields=")+1):length(opts2)] <- "FACTOR="
      }else{
        NULL
      }
      opts2["Infile="]<-paste(paste(gsub("/","\\\\\\\\",path),"\\\\",dat.file[i],sep="")," /NoEcho;",sep="")
      opts2["End2;"]<-""
      opts2["Estimate;"]<-""
      opts2["Distance /Intervals="]<-paste(paste(breaks,collapse=",")," /Width=",max(breaks)," /Left=",min(breaks),";",sep="")
      
      if(is.null(stratum)){
        opts2["Density="]<-"All;"
        opts2["Encounter="]<-"All;"
        opts2["Detection="]<-"All;"
        opts2["Size="]<-"All;"
      }else{
        if(stratum=="STR_LABEL"){
          
          opts2["Density="]<-"Stratum /DESIGN=strata /WEIGHT=area;" # in post stratify and stratify only, both lines are the same
          opts2["Encounter="]<-"Stratum;"
          opts2["Detection="]<-paste(detection,";",sep="")
          opts2["Size="]<-"All;"
        }else{
          opts2["Density="]<-"Stratum /DESIGN=strata /WEIGHT=area;"
          opts2["Encounter="]<-"Stratum;"
          opts2["Detection="]<-paste(detection,";",sep="")
          opts2["Size="]<-"All;"
          opts2["Fields="]<-paste(paste(names(dat),collapse=", "),";",sep="")
          
        }
      }	
      va<-if(!is.null(factor) | !is.null(covariates)){TRUE}else{FALSE}
      if(!is.null(covariates) | !is.null(factor)){ 
        factor_list<-paste(factor,collapse=", ")
        covariates_list<-paste(covariates,collapse=", ")
        covariates_list<-paste(factor_list,covariates_list,sep=", ")
      }
      if(is.null(rare)){  #fits a uniform model when rare is not NULL
        #browser()
        if(is.null(estimator)){
          if(is.null(covariates) & is.null(factor)){
            opts2[["Estimator1"]]<-gsub(", ;",";",gsub("=, ","=",paste(paste(paste(c(" /Key="," /Adjust="," /Criterion=",if(va){" /Covariates="}else{NULL}),c("UN","CO","AIC",if(va){covariates_list}else{NULL}),sep=""),collapse=""),";",sep="")))
            opts2[["Estimator2"]]<-gsub(", ;",";",gsub("=, ","=",paste(paste(paste(c(" /Key="," /Adjust="," /Criterion=",if(va){" /Covariates="}else{NULL}),c("UN","PO","AIC",if(va){covariates_list}else{NULL}),sep=""),collapse=""),";",sep="")))
          }
          opts2[["Estimator3"]]<-gsub(", ;",";",gsub("=, ","=",paste(paste(paste(c(" /Key="," /Adjust="," /Criterion=",if(va){" /Covariates="}else{NULL}),c("HN","CO","AIC",if(va){covariates_list}else{NULL}),sep=""),collapse=""),";",sep="")))
          opts2[["Estimator4"]]<-gsub(", ;",";",gsub("=, ","=",paste(paste(paste(c(" /Key="," /Adjust="," /Criterion=",if(va){" /Covariates="}else{NULL}),c("HN","HE","AIC",if(va){covariates_list}else{NULL}),sep=""),collapse=""),";",sep="")))
          opts2[["Estimator5"]]<-gsub(", ;",";",gsub("=, ","=",paste(paste(paste(c(" /Key="," /Adjust="," /Criterion=",if(va){" /Covariates="}else{NULL}),c("HA","CO","AIC",if(va){covariates_list}else{NULL}),sep=""),collapse=""),";",sep="")))
          opts2[["Estimator6"]]<-gsub(", ;",";",gsub("=, ","=",paste(paste(paste(c(" /Key="," /Adjust="," /Criterion=",if(va){" /Covariates="}else{NULL}),c("HA","PO","AIC",if(va){covariates_list}else{NULL}),sep=""),collapse=""),";",sep="")))
        }else{
          if(!is.list(estimator) | !all(sapply(estimator,length)==2)){stop("Wrong list format for argument estimator")}
          for(j in 1:length(estimator)){
            opts2[[paste("Estimator",j,sep="")]]<-gsub(", ;",";",gsub("=, ","=",paste(paste(paste(c(" /Key="," /Adjust="," /Criterion=",if(va){" /Covariates="}else{NULL}),c(estimator[[j]][1],estimator[[j]][2],"AIC",if(va){covariates_list}else{NULL}),sep=""),collapse=""),";",sep="")))
          }
        }	
      }else{
        opts2[["Estimator1"]]<-gsub(", ;",";",gsub("=, ","=",paste(paste(paste(c(" /Key="," /Adjust="," /Criterion="," /NAP=",if(va){" /Covariates="}else{NULL}),c("UN","CO","AIC","0",if(va){covariates_list}else{NULL}),sep=""),collapse=""),";",sep="")))
      }
      opts2["Monotone="]<-paste(monotone,";",sep="")
      opts2["Pick="]<-"AICC;"
      opts2["GOF;"]<-""
      opts2["Cluster /Bias="]<-"GXLOG;"
      opts2["VarN="]<-"Empirical;"
      #opts2["Multiplier1="]<-paste(multiplier," /Label='Sampling fraction';",sep="")
      opts2["Multiplier1="]<-paste(multiplier[1]," /Label='Sampling fraction'"," /SE=",multiplier[2]," /DF=",multiplier[3],";",sep="") 
      if(!is.null(rare)){
        opts2["Multiplier2="]<-paste(1/pdetect," /Label='Rare'"," /SE=",sedetect/pdetect^2," /DF=",dfdetect,";",sep="")  
      }
      opts2["End3;"]<-""
      
      names(opts2)[grep("Estimator",names(opts2))]<-"Estimator" #to get the nb out which are used to prevent overwriting previously
      names(opts2)[grep("End",names(opts2))]<-"End;"
      names(opts2)[grep("Multiplier",names(opts2))]<-"Multiplier="
      opts2<-lapply(opts2,paste,collapse="")
      
      
      n.model <- sum(names(opts2)=="Estimator")
      inp.file<- sapply(1:n.model, function(j){paste(paste("input",names(dataset)[i],
                                                           substr(opts2[names(opts2)=="Estimator"][j][[1]],7,8),
                                                           substr(opts2[names(opts2)=="Estimator"][j][[1]],18,19),
                                                           sep="_"),".tmp",sep="")})
      res.file<- sapply(1:n.model, function(j){paste(paste("results",names(dataset)[i],
                                                           substr(opts2[names(opts2)=="Estimator"][j][[1]],7,8),
                                                           substr(opts2[names(opts2)=="Estimator"][j][[1]],18,19),
                                                           sep="_"),".tmp",sep="")})
      log.file<- sapply(1:n.model, function(j){paste(paste("log",names(dataset)[i],
                                                           substr(opts2[names(opts2)=="Estimator"][j][[1]],7,8),
                                                           substr(opts2[names(opts2)=="Estimator"][j][[1]],18,19),
                                                           sep="_"),".tmp",sep="")})
      det.file<- sapply(1:n.model, function(j){paste(paste("detections",names(dataset)[i],
                                                           substr(opts2[names(opts2)=="Estimator"][j][[1]],7,8),
                                                           substr(opts2[names(opts2)=="Estimator"][j][[1]],18,19),
                                                           sep="_"),".tmp",sep="")})
      
      
      input<-lapply(1:n.model, function(j){
        opts1<-list()
        opts1[""]<-paste(gsub("/","\\\\\\\\",path),"\\\\",res.file[j],sep="")
        opts1[""]<-paste(gsub("/","\\\\\\\\",path),"\\\\",log.file[j],sep="")
        opts1[""]<-"None"
        opts1[""]<-paste(gsub("/","\\\\\\\\",path),"\\\\",det.file[j],sep="")
        opts1[""]<-"None"
        opts1[""]<-"None"
        opts <- c(opts1,opts2)
        model.lines <- which(names(opts)=="Estimator")
        if(length(model.lines)==1){
          out <- paste(names(opts),unlist(opts),sep="") 
        }else{
          out <- paste(names(opts)[-model.lines[-j]],unlist(opts)[-model.lines[-j]],sep="")  
        }
        return(out)
      })
      
      if(verbose){ #prints the input file
        if(is.list(input) == FALSE){
          cat(input,fill=1)
        }else{
          #browser()
          lapply(input, function(x){cat(x, "\n\n", fill = 1)})
        }
      }# Modification pour beug avec print de liste
      for(j in 1:n.model){
        write.table(input[[j]],file.path(path,inp.file[j]),row.names=FALSE,quote=FALSE,col.names=FALSE)   
      }
      
      dat$SMP_LABEL<-paste(as.numeric(factor(dat$SMP_LABEL)),dat$SMP_LABEL,sep=".")
      dat<-dat[order(dat[,"STR_LABEL"],dat[,"SMP_LABEL"]),] #always order according to the first column/ str_label, otherwise MCDS creates too many stratum or samples
      w<-which(is.na(dat[,"DISTANCE"]) | dat[,"DISTANCE"]=="")
      if(any(w)){
        dat[w,"SIZE"]<-"" #previously NA was given, but gives status 3 with distance
        dat[w,"DISTANCE"]<-"" #previously not here but don't know why
        if(!is.null(covariates)){for(j in 1:length(covariates)){dat[w,covariates[j]]<-""}}
        if(!is.null(factor)){for(j in 1:length(factor)){dat[w,factor[j]]<-""}}
      }
      dat<-dat[dat[,"STR_LABEL"]!="",]#temporary to get rid of stratum when it is species and it is an empty transect, has to be clarified what to do
      write.table(dat,file.path(path,dat.file[i]),na="",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
      
      ####################################################
      ### running distance
      run.cmd <- vector(mode="list",length=n.model)
      for(j in 1:n.model){
        if(file.exists(file.path(pathMCDS,"MCDS.exe"))){
        cmd<-paste(shQuote(file.path(pathMCDS,"MCDS"))," 0, ",shQuote(file.path(path,inp.file[j])),sep="")
        run.cmd[j] <- try(system(cmd,wait=TRUE,ignore.stdout=FALSE,ignore.stderr=FALSE,invisible=TRUE), silent = TRUE)
        }else{
          stop("pathMCDS doesn't exist")
        }
      }
      
      ####################################################
      ### output  
      ### subset original observations
      input.data <- list()
      input.data[[1]] <- dat
      input.data[[2]] <- breaks
      names(input.data)<-c("observations","breaks")
      mans <-vector(mode="list",length=n.model)
      names(mans) <-  sapply(1:n.model, function(j){paste(names(dataset)[i],
                                                          substr(opts2[names(opts2)=="Estimator"][j][[1]],7,8),
                                                          substr(opts2[names(opts2)=="Estimator"][j][[1]],18,19), sep="_")})
      
      #list of model
      for(j in 1:n.model){
        if(run.cmd[[j]]%in%c(1,2,3)){
          log<-readLines(file.path(path,log.file[j]))
          res<-readLines(file.path(path,res.file[j]))
          if(any(grep("No fit possible",log)) || length(res)<1L){
            notrun<-c(notrun,names(lans)[j])
            next
          }
          x<-readLines(file.path(path,res.file[j]))
          y<-readLines(file.path(path,det.file[j]))
          ans<-vector(mode="list",length=10)
          ans[[1]]<- input.data
          ans[[2]]<- model_fittingMCDS(x=x, units=units)
          ans[[3]]<- parameter_estimatesMCDS(x)
          ans[[4]]<- ifelse(!is.null(factor) | !is.null(covariates)==T,param_namesMCDS(x),"No covariates in the model") 
          ans[[5]]<- chi_square_testMCDS(x)
          ans[[6]]<- density_estimateMCDS(x)
          ans[[7]]<- ifelse(is.null(rare)==T,cluster_sizeMCDS(x),"No cluster estimation for the rare species model")
          ans[[8]]<- detection_probabilityMCDS(y, covariates=covariates)
          ans[[9]]<- AIC_MCDS(x)
          ans[[10]]<- path
          names(ans)<-c("input_data","model_fitting","parameter_estimates","covar_key","chi_square_test","density_estimate","cluster_size","detection","AIC","path")
          class(ans)<-"distanceFit"
          #list of model
          mans[[j]] <- ans
        }else{
          mans[[j]] <- "Model failed to converge" 
        }
      } 
      if(length(mans)==1){
        mans <- mans[[1]]
      }else{
        class(mans) <- "distanceList" 
      }
      #list of species
      lans[[i]]<-mans
    }
    if(!is.null(notrun)){warning(paste("No models for the following combinations: ",paste(notrun,collapse=" "),". See log files." ,sep=""),immediate.=TRUE)}
    lans<-lans[!sapply(lans,is.null)]
    class(lans) <- "speciesList"
    if(length(lans)==1){lans[[1]]}else{lans}
  }


