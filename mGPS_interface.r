library(geosphere)
library(sp)
library(rworldmap)
library(caret)
library(maps)
library(geosphere)
library(caret)
library(plyr)
library(rgeos)
library(mapplots)
library(shiny)
library(dplyr)
library(xgboost)
library(e1071)
library(rgdal)
library(maptools)
library(sp)


setwd(rprojroot::find_rstudio_root_file())

# Part 1 Function ------------------------------------------------------------------

## Data transformation 
data_normalise <- function(df) {
  return(df/rowSums(df))
}


## Feature selection algorithm
species_select <-
  function(x,
           y,
           remove_correlated = T,
           subsets = NULL,
           cores = 1) {
    doParallel::registerDoParallel(cores)
    
    y <- factor(y)
    
    print("feature elimination --------------------")
    
    if (remove_correlated == T) {
      correlated <- findCorrelation(
        cor(x, use = "complete.obs"),
        cutoff = 0.98,
        verbose = FALSE,
        names = FALSE,
        exact = FALSE
      )
      
      x <- x[,-c(correlated)]
      print(paste0("correlated features removed:",length(correlated)))
      
    }
    
    len <- ncol(x)
    if(is.null(subsets)){
      subsets <-
        c(floor(len/2), floor(len / 4), floor(len / 8), floor(len / 16), floor(len / 32),floor(len / 64))
    }
    
    rfe_ctrl <- rfeControl(
      functions = rfFuncs,
      method = "cv",
      number =  5,
      verbose = FALSE,
      allowParallel = TRUE
    )
    set.seed(123)
    featureElimination <- rfe(
      x = x,
      y = y,
      sizes = subsets,
      rfeControl = rfe_ctrl,
      tuneLength = 2
    )
    doParallel::registerDoParallel(1) 
    
    print("feature elimination has been done ...")
    
    return(featureElimination)
  }


## Main mGPS algorithm 
mGPS <-
  function(training = NULL,
           testing = NULL,
           classTarget,
           hierarchy = c('continent','city','latitude','longitude'),
           variables,
           nthread = 1,
           coast = NULL) {
    
    #Check for training set
    if(is.null(training)){
      return(message("No training set given"))
    } else{
      
      training <- droplevels(training)
      training[,classTarget] <- factor(training[,classTarget])
      training[,hierarchy[1]] <- factor(training[,hierarchy[1]])
      
      #Train mGPS with 5-fold cross validation of training set for hyperparameter tuning. 
      message("Training mGPS...")
      
      set.seed(1234)
      
      folds <- createFolds(training[,classTarget], k = 5, returnTrain = T)
      
      trControlClass <-  trainControl(
        method = "cv",
        number = 5,  
        verboseIter = FALSE,
        returnData = FALSE,
        search = "grid",
        savePredictions = "final",
        classProbs = T,
        allowParallel = T,
        index = folds )
      
      trControl <-  trainControl(
        method = "cv",
        number = 5,  
        verboseIter = FALSE,
        returnData = FALSE,
        search = "grid",
        savePredictions = "final",
        allowParallel = T,
        index = folds)
      
      tune_grid <- expand.grid(
        nrounds = c(300,600),
        eta = c( 0.05, 0.1),
        max_depth = c(3,6,9),
        gamma = 0,
        colsample_bytree = c(0.6,0.8),
        min_child_weight = c(1),
        subsample = (0.7)
      )
      
      if(length(hierarchy) == 4){
        
        print("--------- Continent training ---------")

        Xgb_region <- train(x = training[,variables],y = training[,hierarchy[1]],
                            method = "xgbTree",
                            trControl = trControlClass,
                            tuneGrid = tune_grid,
                            nthread = nthread)
        
        print("Continent training has been done ...")

        
        l1_train <- data.frame(training[,c(variables)],Xgb_region[["pred"]][order(Xgb_region$pred$rowIndex),levels(training[,hierarchy[1]]) ])

      }else{
        
        l1_train <- training[,c(variables)]
        
      }
      

      print("--------- City training --------- ")
      
      Xgb_class <- train(x = l1_train,y = training[,classTarget],
                         method = "xgbTree",
                         trControl = trControlClass,
                         tuneGrid = tune_grid,
                         nthread = nthread)
      
      print("City training has been done ...")
      
      l2_train <- data.frame(l1_train,Xgb_class[["pred"]][order(Xgb_class$pred$rowIndex),levels(training[,classTarget]) ])
      
      print("--------- latitude training ---------")
      
      if (length(hierarchy)==3){
        lat <- hierarchy[2]
      }else if(length(hierarchy)==4){
        lat <- hierarchy[3]
      }
      
      Xgb_latitude <- train(x = l2_train ,y = training[,lat],
                            method = "xgbTree",
                            trControl = trControl,
                            tuneGrid = tune_grid,
                            nthread = nthread)
      print("Latitude training has been done ...")
      
      l3_train <- data.frame(l2_train, "latPred" = Xgb_latitude[["pred"]][order(Xgb_latitude$pred$rowIndex),"pred" ])
      
      print("--------- longitude training ---------")
      
      if (length(hierarchy)==3){
        long <- hierarchy[3]
      }else if(length(hierarchy)==4){
        long <- hierarchy[4]
      }
      
      Xgb_longitude <- train(x = l3_train ,y = training[,long],
                             method = "xgbTree",
                             trControl = trControl,
                             tuneGrid = tune_grid,
                             nthread = nthread)
      
      print("Longitude training has been done ...")
      
    }
    
    #check for test set, return trained model if no test set. 
      model <- function(test,variables){
        
        if(length(hierarchy) == 4){
          regProbs <- predict(Xgb_region, newdata = test[,variables],type ="prob")
          
          l1_test <- data.frame(test[,variables], regProbs)
        }else{
          l1_test <- test[,variables]
        }
        classPred <- predict(Xgb_class, newdata = l1_test)
        classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
        
        l2_test <-  data.frame(l1_test, classProbs) 
        latPred <- predict(Xgb_latitude, newdata = l2_test)
        
        l3_test <- data.frame(l2_test, latPred)
        longPred <- predict(Xgb_longitude, newdata = l3_test)
        return(list(classPred, latPred, longPred))
        
      }

      message("Generating predictions...")
      #generate mGPS predictions for test set

      features <- variables
      newdata_in <- as.data.frame(setNames(replicate(length(features),numeric(0), simplify = F), features))
      for (f in features){
        if (f %in% colnames(testing)){
          newdata_in[c(1:dim(testing)[1]),f] <- testing[,f]
        }else{
          newdata_in[c(1:dim(testing)[1]),f] <- data.frame(0)
        }
      }
      
      if(length(hierarchy) == 4){
        
        regProbs <- predict(Xgb_region, newdata = newdata_in,type ="prob")
        
        l1_test <- data.frame(newdata_in, regProbs)
      }else{
        l1_test <- newdata_in
      }
      
      classPred <- predict(Xgb_class, newdata = l1_test)
      classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
      
      l2_test <-  data.frame(l1_test, classProbs) 
    
      latPred <- predict(Xgb_latitude, newdata = l2_test)
      
      l3_test <- data.frame(l2_test, latPred)
      
      longPred <- predict(Xgb_longitude, newdata = l3_test)
      
      print("Prediction has been done ...")
      
      #adjust out of bounds predictions
      longPred[longPred > 180] <- 180
      longPred[longPred < -180] <- -180
      latPred[latPred > 90] <- 90
      latPred[latPred < -90] <- -90
      
      if(length(hierarchy) == 3){
        model_store <- list(Xgb_class,Xgb_latitude,Xgb_longitude,"model" = model)
      }else  if(length(hierarchy) == 4){
        model_store <- list(Xgb_region,Xgb_class,Xgb_latitude,Xgb_longitude,"model" = model)
      }
      prediction_store <- list(classPred, latPred, longPred)
      
      return(list(model_store, prediction_store))
  }


### Function: Pull points to marine and land

# Pull to nearest coastline

pull_land <- function(land_preds){

  coastlines <- cbind("x"  = maps::SpatialLines2map(rworldmap::coastsCoarse)$x ,"y" =maps::SpatialLines2map(rworldmap::coastsCoarse)$y)
  coastlines <- coastlines[complete.cases(coastlines),]
  coastlines <- coastlines[coastlines[,1] < 180 ,]
  
  find_coast <- function(long, lat) {
    
    distances_from_coastline <-
      sp::spDistsN1(coastlines, c(long, lat), longlat = TRUE)
    
    closest_point <-  which.min(distances_from_coastline)
    new_coords <- coastlines[closest_point,]
    
    return(new_coords)
  }

  toAdjust <-
    which(is.na(maps::map.where(database = "world", land_preds$longPred, land_preds$latPred)))
  
  adjusted <-
    mapply(find_coast, long = land_preds$longPred[toAdjust], lat = land_preds$latPred[toAdjust])
  
  land_preds$longPred[toAdjust] <- adjusted[1,]
  land_preds$latPred[toAdjust] <- adjusted[2,]
  
  return(land_preds)
}


# Pull to nearest water body

pull_marine <- function(marine_preds){


  seas <- rgdal::readOGR(dsn = "Data/Geo/ne_10m_geography_marine_polys", layer = "ne_10m_geography_marine_polys")
  coastlines <- cbind("x"  =maps::SpatialPolygons2map(seas)$x ,"y" =maps::SpatialPolygons2map(seas)$y)
  coastlines <- coastlines[complete.cases(coastlines),]
  coastlines <- coastlines[coastlines[,1] < 180 ,]
  
  find_coast2 <- function(long,lat){
    distances_from_coastline <-  sp::spDistsN1(coastlines , c(long,lat), longlat = TRUE)
    closest_point <-  which.min(distances_from_coastline)
    new_coords <- coastlines[closest_point,]
    return(new_coords)
  }
  
  data(wrld_simpl)
  
  ## Create a SpatialPoints object
  set.seed(0)
  points <- data.frame(marine_preds$longPred, marine_preds$latPred) 
  pts <- SpatialPoints(points, proj4string=CRS(proj4string(wrld_simpl)))
  ## Change CRS if program has error for different CRS in over() step
  proj4string(wrld_simpl) <- pts@proj4string
  ## Find which points fall over land
  ii <- !is.na(over(pts, wrld_simpl)$FIPS)
  toAdjust <- marine_preds[which(ii == TRUE),]
  adjusted <- mapply(find_coast2, long = toAdjust$longPred, lat = toAdjust$latPred )
  
  marine_preds[which(ii == TRUE), "latPred"] <- adjusted[2,]
  marine_preds[which(ii == TRUE), "longPred"] <- adjusted[1,]
  
  return(marine_preds)
}


data_preprocess_f <- function(train_f,target_in,hierarchy,remove_small){

  #Remove control samples
  metasub_data <- droplevels(train_f)
  
  if(length(hierarchy) == 4){
    metasub_data[,hierarchy[1]] <- factor(metasub_data[,hierarchy[1]])
    metasub_data[,hierarchy[2]] <- factor(metasub_data[,hierarchy[2]])
  }else{
    metasub_data[,hierarchy[1]] <- factor(metasub_data[,hierarchy[1]])
  }
  #remove sparse samples locations and dubiously labelled samples. 
  
  remove_small <- as.numeric(remove_small)
  
  if (remove_small > 0){
  
    small_cities <-  names(which(summary(metasub_data[,target_in]) < remove_small))
      
    remove_samples <- which(metasub_data[,target_in] %in%  c("antarctica", small_cities))
  
    if (length(remove_samples) != 0 ){
      metasub_data <- droplevels(metasub_data[-c(remove_samples), ])
    }
  }

  #Correction of identified misslabelling of data 
  
  for (i in 1:length(hierarchy)){
    empty <- which(metasub_data[,hierarchy[i]] == "" | is.na(metasub_data[,hierarchy[i]]))

    if (length(empty) != 0  ){
      metasub_data <- metasub_data[-c(empty),]
    }
  }

  metasub_data[,hierarchy[1]] <- droplevels(metasub_data[,hierarchy[1]])
  
  if(length(hierarchy) == 4){
    metasub_data[,hierarchy[2]] <- droplevels(metasub_data[,hierarchy[2]])
  }

  print("metasub_data has been pre-processed ...")
  return(metasub_data)
}


featureElimination_f <- function(metasub_data,classTarget_in,range_1,range_2,subsets_in){

  ### Find GITs ####
  range_1 <- as.numeric(range_1)
  range_2 <- as.numeric(range_2)
 
  featureElim <- species_select(x = metasub_data[, c(range_1:range_2)],y = metasub_data[,classTarget_in],remove_correlated = F, subsets = subsets_in ,cores = 8)

  return(featureElim)
}


v_f <- function(featureElim){
  optVars <- featureElim$optVariables
  v <- varImp(featureElim$fit, type = 1, scale = F)
  v[,"taxa"] <- row.names(v)
  v <- v[order(v$Overall,decreasing = T),]
  write.csv(v, file = "Outputs/Optimal_features.csv")
  return(v)
}
  

git_subset_f <- function(featureElim){
  git_subset <- data.frame("n_vars" = featureElim$results$Variables, "accuracy" = featureElim$results$Accuracy)
  write.csv(git_subset,file = "Outputs/Features_subset_accuracy.csv")
  return(git_subset)
}
  
 
# Accuracy of model detected on original dataset
model_accuracy_f <- function(metasub_data,optVars,classTarget_in,hierarchy_in){
  set.seed(18)
 
  #generate 5 stratified folds for test predictions.
  trainFolds <-  createFolds(metasub_data[,classTarget_in], k = 5, returnTrain = T)
  
  GeoPreds <- list()
  
  #iteratively train the model on each of the 5 training folds and generate predictions using the coresponding test fold.
  for (i in 1:5){
    
    print(i)
    print("-------------------------------------------")
    
    train <- metasub_data[trainFolds[[i]],]
    test <- metasub_data[-trainFolds[[i]],]
    
    testPreds <-mGPS(training = train, testing = test, classTarget = classTarget_in,variables = optVars,nthread = 8,hierarchy = hierarchy_in, coast=NULL)

    GeoPreds[[i]] <- testPreds[[2]]

    print(i)
    print("has been done .....")
    
  }
  
  model_last_one <- testPreds[[1]]
  
  #Combine these test predictions into one data set 
  add_preds <- list()
  for (i in 1:5){
    
    add_preds[[i]] <- cbind(metasub_data[-trainFolds[[i]],] , 
                            "cityPred"= GeoPreds[[i]][[1]], 
                            "latPred" = GeoPreds[[i]][[2]], 
                            "longPred" = GeoPreds[[i]][[3]] )
  }
  
  MetasubDataPreds <- rbind.fill(add_preds)
  hierarchy = hierarchy_in
  for (i in 1:nrow(MetasubDataPreds)){
    
    if(length(hierarchy) == 3){
    
      MetasubDataPreds[i,"Distance_from_origin"] <- geosphere::distm(c(MetasubDataPreds[i,"longPred"],MetasubDataPreds[i,"latPred"]), c(MetasubDataPreds[i,hierarchy[3]],MetasubDataPreds[i,hierarchy[2]]), fun = geosphere::distHaversine)/1000
  
    }else  if(length(hierarchy) == 4){
      MetasubDataPreds[i,"Distance_from_origin"] <- geosphere::distm(c(MetasubDataPreds[i,"longPred"],MetasubDataPreds[i,"latPred"]), c(MetasubDataPreds[i,hierarchy[4]],MetasubDataPreds[i,hierarchy[3]]), fun = geosphere::distHaversine)/1000
    }
  }
  
  write.csv(MetasubDataPreds,"Outputs/Prediction_results.csv")
  
  save(model_last_one,file="Outputs/Prediction_model.Rda")
  
  return(list(MetasubDataPreds,model_last_one)) 
}


# Function for new model construction and sample origin prediction
New_model_prediction <- function(train,test,classTarget_in,optVars,hierarchy_in){
  
  testPreds <-mGPS(training = train, testing = test, classTarget = classTarget_in,variables = optVars,nthread = 8,hierarchy = hierarchy_in, coast=NULL)
  
  model <- testPreds[[1]]
  GeoPreds <- testPreds[[2]]
  
  MetasubDataPreds <- cbind(test, 
                            "cityPred"= GeoPreds[[1]], 
                            "latPred" = GeoPreds[[2]], 
                            "longPred" = GeoPreds[[3]] )
    
  write.csv(MetasubDataPreds,"Outputs/Prediction_result.csv")
  
  save(model,file="Outputs/Prediction_model.Rda")
  
  return(list(MetasubDataPreds,model))
}


# Extract features that used in built prediction model
feature_filter <- function(test_dataset,model){
  features <- (model$finalModel$feature_names)
  newdata <- as.data.frame(setNames(replicate(length(features),numeric(0), simplify = F), features))
  for (f in features){
    if (f %in% colnames(test_dataset)){
      newdata[c(1:dim(test_dataset)[1]),f] <- test_dataset[,f]
    }else{
      newdata[c(1:dim(test_dataset)[1]),f] <- data.frame(0)
    }
  }
  return(newdata)
}


# Use MetaSub model to predict new samples
MetaSub_prediction <- function(model_store,test_dataset){
  Xgb_region <- model_store[[1]]
  Xgb_class <- model_store[[2]]
  Xgb_latitude <- model_store[[3]]
  Xgb_longitude <- model_store[[4]]
  # prediction_model <- model_store[[5]]
  
  print("---------- Prediction ----------")
  
  model <- function(test_dataset){
    
    newdata <- feature_filter(test_dataset, Xgb_region) 
    
    regProbs <- predict(Xgb_region, newdata, type ="prob")
    
    l1_test <- data.frame(newdata, regProbs)
    
    classPred <- predict(Xgb_class, newdata = l1_test)
    classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
    
    l2_test <-  data.frame(l1_test, classProbs) 
    latPred <- predict(Xgb_latitude, newdata = l2_test)
    
    l3_test <- data.frame(l2_test, latPred)
    longPred <- predict(Xgb_longitude, newdata = l3_test)
    
    longPred[longPred > 180] <- 180
    longPred[longPred < -180] <- -180
    latPred[latPred > 90] <- 90
    latPred[latPred < -90] <- -90
    
    return(list(classPred, latPred, longPred))
  }
  
  testPreds <- model(test_dataset)
  
  DataPreds <- cbind(test_dataset , 
                     "cityPred"= testPreds[[1]], 
                     "latPred" = testPreds[[2]], 
                     "longPred" = testPreds[[3]] )
  
  write.csv(DataPreds,"Outputs/Predicted_origin.csv")
  
  return(DataPreds)
}


# Plot predicted location of samples in training dataset. It can show how many samples are located into their continent.
plot_map <- function(MetasubDataPreds,hierarchy_in,classTarget_in,x_ran,y_ran){
  
  par(mai=c(2,1,0.5,0.5), mar=par()$mar+c(3,0,0,0))
  
  
  
  map <- rworldmap::getMap(resolution = "coarse")
  palette <-c( "darkorchid4","gold2","dodgerblue3","brown","orangered2","mediumspringgreen","deeppink2")
  plot(map, xlim = x_ran,ylim = y_ran,col = "grey",border = "darkgrey", xlab = "", ylab = '', bg = "lightskyblue1")
  title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
  
  for (i in 1:(length(hierarchy_in)-2)){
    MetasubDataPreds[,hierarchy_in[i]] <- factor(MetasubDataPreds[,hierarchy_in[i]])
  }
  
  MetasubDataPreds$cityPred <- factor(MetasubDataPreds$cityPred)
  
  # get the coordination of continent
  continent_list <- c("east_asia","europe","middle_east","north_america","oceania","south_america", "sub_saharan_africa",'Australia','North_America','Europe','Africa','South_america','Asia')
  continent_lats <- c(55,69,8,40,-40,-10,-5,-40,40,69,-5,-10,55)
  continent_longs <- c(125,0,60,-130,140,-80,5,140,-130,0,5,-80,125)
  continent_position <- data.frame(cbind(continent_list,continent_lats,continent_longs))
  
  label_continent <- c()
  flag <- 0
  for ( i in 1:length(levels(MetasubDataPreds[,hierarchy_in[1]]))){
    this_continent <- levels(MetasubDataPreds[,hierarchy_in[1]])[i]
    label_continent <- c(label_continent,this_continent)
    find_lats <- MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,][,"latPred"]
    find_longs <- MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,][,"longPred"]
    
    #plot predicted co-ordinates
    points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.2,xlim = c(-165,168))
    
    #plot city prediction accuravy by continent as pies
    
    correctly_pred <-  mean(as.numeric(as.character(MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,"cityPred"])== 
                                                      as.character(MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,classTarget_in]))) 
    
    incorrectly_pred <- (1 - correctly_pred)
    
    if (this_continent %in% continent_position$continent_list){
      add.pie(z = c(correctly_pred, incorrectly_pred), x = as.numeric(continent_position[continent_position$continent_list==this_continent,3]), y = as.numeric(continent_position[continent_list==this_continent,2])
              ,edges=200,
              radius=10,
              col=c(palette[i],"black") , labels = ""
      )
    }else{
      
      add.pie(z = c(correctly_pred, incorrectly_pred), x = as.numeric(-150+flag), y = 105
              ,edges=200,
              radius=10,
              col=c(palette[i],"black") , labels = ""
      )
      flag <- flag + 25
    }
  }
  map.axes(cex.axis = 1.1)
  
  legend(xpd = T,'bottom',inset=c(0,-0.35), label_continent, pch = "+", col = palette[1:length(label_continent)], cex = 0.8,n=4)
  
  box( col = 'black')
  par(mar=c(5, 4, 4, 2) + 0.1)
  
}


 
# Plot the prediction origin of sample on worldmap
plot_prediction <- function(MetasubDataPreds,x_ran,y_ran){
  
  par(mai=c(1,1,0.5,0.5))

  map <- rworldmap::getMap(resolution = "coarse")

  plot(map, xlim = x_ran,ylim = y_ran, col = "grey",border = "darkgrey", xlab = "", ylab = '', bg = "lightskyblue1")
  title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)

  find_lats <- MetasubDataPreds$latPred
  find_longs <- MetasubDataPreds$longPred
  
  points(find_longs,find_lats, type = "p",col = "purple", pch = "+", cex = 1.3)
  map.axes(cex.axis = 1.1)
  
}


# Plot of accuracy

plot_accuracy <- function(MetasubDataPreds,classTarget_in){
  MetasubDataPreds[,classTarget_in] <- factor(MetasubDataPreds[,classTarget_in])
  
  bar_df1 <- data.frame(row.names = c("Overall",levels(MetasubDataPreds[,classTarget_in])))
  
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 100)
    bar_df1[i+1,"0 - 100km"] <- prop
    
    overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] < 100)
    bar_df1[ 1,"0 - 100km"] <- overall_prop
  }
  
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 100 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 500)
    bar_df1[i+1,"100 - 500km"] <- prop
    
    overall_prop <-mean(MetasubDataPreds[,"Distance_from_origin"] > 100 & MetasubDataPreds[,"Distance_from_origin"] < 500)
    bar_df1[ 1,"100 - 500km"] <- overall_prop
  }
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 500 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 1000)
    bar_df1[i+1,"500 - 1000km"] <- prop
    
    overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 500 & MetasubDataPreds[,"Distance_from_origin"] < 1000)
    bar_df1[ 1,"500 - 1000km"] <- overall_prop
  }
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 1000 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 2000)
    bar_df1[i+1,"1000 - 2000km"] <- prop
    
    overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 1000 & MetasubDataPreds[,"Distance_from_origin"] < 2000)
    bar_df1[ 1,"1000 - 2000km"] <- overall_prop
  }
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 2000 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 3000)
    bar_df1[i+1,"2000 - 3000km"] <- prop
    
    overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 2000 & MetasubDataPreds[,"Distance_from_origin"] < 3000)
    bar_df1[1,"2000 - 3000km"] <- overall_prop
  }
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 3000 )
    bar_df1[i+1,"> 3000km"] <- prop
    
    overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 3000)
    bar_df1[ 1,"> 3000km"] <- overall_prop
  }
  
  size <- c()
  for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
    
    this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
    size[i] <- length(which(MetasubDataPreds[,classTarget_in] == this_city))
  }
  
  par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)
  bp <- barplot(t(bar_df1*100), space = 0,col=c("lightyellow","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), 
                names.arg=c("Overall",paste0(levels(MetasubDataPreds[,classTarget_in]),"  (",size,")"), axes = FALSE) , 
                las =2, cex.names=.6, ylab = "", axisnames = F, axes = F)
  axis(side =2, pos = 0)
  mtext(text = c("Overall",paste0(levels(MetasubDataPreds[,classTarget_in])," (",size,")")), side = 1, at = bp, line = 0, padj = 1, cex = 0.7)
  title(ylab="Proportion of sample predictions %", mgp=c(2,1,0),cex.lab=1)
  legend("topright",inset = c(-0.26,0.4), rev(c(colnames(bar_df1))), 
         fill = rev(c("lightyellow","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) ,
         bty = 1, cex = 0.8)
  par(mar=c(5, 4, 4, 2) + 0.1)
}


# Part 2 shiny ------------------------------------------------------------------
  
# ui

library(shiny)

ui <- fluidPage(
  titlePanel('Geographical origin prediction of microbiome'),
  
  tags$head(
    tags$style(HTML("
    .shiny-output-error-validation {
    color: red;
    }
    "))
  ),
  
  sidebarLayout(
    sidebarPanel(
      
      # Data loading tips
      tags$head(tags$style(type="text/css", "
                #loadmessage {
                top: 0px; left: 0px;
                width: 100%; padding: 5px 0px 5px 0px;
                text-align: center; font-weight: bold;
                font-size: 100%; color: #000000;
                background-color: #FFC1C1; z-index: 105;}"),
                ), 
      conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                       tags$div("Data loading....",id="loadmessage")),
      
      # Program selection
      selectInput("program","Prediction program",
                  c("Build a new prediction model using mGPS" = "cc",
                    "Build a new prediction model using mGPS and predict new samples" = "bb",
                    "Use an existing model to predict new samples" = "aa"
                    )
                  ),
      
      # function : Microbial origin prediction based on MetaSub
      conditionalPanel(
        condition = "input.program == 'aa'",
        fileInput(inputId = "f_new_test_1",
                  label = "Upload sample(s) abundance file",
                  accept = c("text/csv",
                             ".csv")),
        fileInput(inputId = "model",
                  label = "Upload the prediction model (In .Rda format)",
                  accept = ".Rda"),
        actionButton(inputId = "acbutton_1", 
                     label ="Start")
      ),
      
      # function: Build a new prediction model using mGPS and predict new samples
      conditionalPanel(
        condition = "input.program == 'bb'",
        fileInput(inputId = "f_new_test_2",
                  label = "Upload new sample abundance file(s)",
                  accept = c("text/csv",
                             ".csv")),
        radioButtons(inputId = "file_num_2",
                     label = "Upload reference file(s)",
                     c("Merged metadata and abundance file" = "one_file",
                       "Separate metadata and abundance file" = "two_file")
        ),
        conditionalPanel(
          condition = "input.file_num_2 == 'one_file'",
          fileInput(inputId = "f_new_train_2",
                    label = "Upload the reference merged dataset file",
                    accept = c("text/csv",
                               ".csv")),
        ),
        conditionalPanel(
          condition = "input.file_num_2 == 'two_file'",
          fileInput(inputId = "f_metadata_2",
                    label = "Upload the metadata file",
                    accept = c("text/csv",
                               ".csv")),
          fileInput(inputId = "f_abundance_2",
                    label = "Upload the abundance file",
                    accept = c("text/csv",
                               ".csv")),
          textInput(inputId = "by_x_2",
                    label = "Merge column name in metadata file"),
          textInput(inputId = "by_y_2",
                    label = "Merge column name in abundance file")
        ),
        textInput(inputId= "target_2",
                  label ="Enter the main locality level"),
        textAreaInput(inputId = "hierarchy_2",
                      label = "Enter the locality hierarchy",
                      rows = 1),
        textAreaInput(inputId = "abundance_range_2",
                      label = "Column range of abundance data",
                      rows = 1),
        
        
        
        checkboxInput('remove_small_2',"Locality sample size cut off (Optional)", value = FALSE),
        conditionalPanel( 
          condition = "input.remove_small_2 == true",
          textAreaInput(inputId = "remove_small_value_2",
                        label = "Cut off of sample number",
                        rows = 1),
        ),
        
        
        checkboxInput("subsets_2", "Subsets in feature elimination (Optional)", value = FALSE),
        conditionalPanel( 
          condition = "input.subsets_2 == true",
          textAreaInput(inputId = "subsets_value_2",
                        label = "Subsets size",
                        rows = 1),
        ),
        
        actionButton(inputId = "acbutton_2", 
                     label ="Start")
      ),
      
      # function: Build a new prediction model using mGPS
      conditionalPanel(
        condition = "input.program == 'cc'",
        radioButtons(inputId = "file_num_3",
                     label = "Input file",
                     c("Merged metadata and abundance file" = "one_file",
                       "Separate metadata and abundance file" = "two_file")
                      ),
        conditionalPanel(
          condition = "input.file_num_3 == 'one_file'",
          fileInput(inputId = "f_new_train_3",
                    label = "Upload the reference merged dataset file",
                    accept = c("text/csv",
                               ".csv")),
        ),
        conditionalPanel(
          condition = "input.file_num_3 == 'two_file'",
          fileInput(inputId = "f_metadata_3",
                    label = "Upload the metadata file",
                    accept = c("text/csv",
                               ".csv")),
          fileInput(inputId = "f_abundance_3",
                    label = "Upload the abundance file",
                    accept = c("text/csv",
                               ".csv")),
          textInput(inputId = "by_x_3",
                    label = "Merge column name in metadata file"),
          textInput(inputId = "by_y_3",
                    label = "Merge column name in abundance file")
        ),
        textInput(inputId= "target_3",
                  label ="Enter the main locality level"),
        textAreaInput(inputId = "hierarchy_3",
                      label = "Enter the locality hierarchy",
                      rows = 1),
        textAreaInput(inputId = "abundance_range_3",
                      label = "Column range of abundance data",
                      rows = 1),
        
        checkboxInput('remove_small_3',"Locality sample size cut off (Optional)", value = FALSE),
        conditionalPanel( 
          condition = "input.remove_small_3 == true",
          textAreaInput(inputId = "remove_small_value_3",
                        label = "Cut off of sample number",
                        rows = 1),
        ),
        
        
        checkboxInput("subsets_3", "Subsets in feature elimination (Optional)", value = FALSE),
        conditionalPanel( 
          condition = "input.subsets_3 == true",
          textAreaInput(inputId = "subsets_value_3",
                        label = "Subsets size",
                        rows = 1),
        ),
        
        actionButton(inputId = "acbutton_3", 
                     label ="Start")
        
    )
  ),
    
  # main panel part
    mainPanel(
      
      conditionalPanel(
        condition = "input.program == 'aa'",
        tabsetPanel(
          type = "tabs",
          tabPanel("HELP",
                   
                   
                   h3("Use an existing model to predict new samples"),
                   h4("_______________________________________________"),
                   
                   br(),
                   
                   h4("Function description"),
                   
                   br(),
                   p("In this mode, user can predict new sample origin based on an exsiting prediction model."),
                   h4("_______________________________________________"),
                   
                   br(),
                   
                   h4("Usage"),
                   br(),
                   p("In left side bar:",style = "font-size:16px"),
                   p("1. Select ",strong("Prediction program")," as ",em("Use an existing model to predict new samples", style = "color:purple")), 
                   p("2. ",strong("Upload sample(s) abundance file:")," Upload data file(s) (in .csv format) containing new microbial sample abundance data."),  
                   p("3. ",strong("Upload the prediction model:")," Upload a constructed origin prediction model in .Rda format. Model can be downloaded in ", strong("Output")," tab of function:",em("Build a new prediction model using mGPS", style = "color:purple"), "or ",em("Build a new prediction model using mGPS and predict new samples", style = "color:purple")),
                   br(),
                   p("In ",strong("Result Plot")," tab:",style = "font-size:16px"),
                   p("4. ",strong("Change longitude/latitude range in output map"),em(" (Optional) ")),
                   p("5. ",strong("Whether pull points to land/marine:"),em(" (Optional) "),"If checked, predicted origin location will be pull to the nearest land/marine if predicted coordinates are out of the expected boarder."),
                   p(strong("Start the program:")," Click the",strong("Start"),"bar and then click the ",strong("Result Plot")," tab",style = "font-size:16px"),
                   p('Data processing, please wait while output files are being generated. When the prompt bar disappears you can see the results and download files.',style = "font-size:16px"),
                   
                   h4("_______________________________________________"),
                   
                   br(),
                   
                   h4("Result plot and output"),
                   br(),
                   p(strong("Result Plot")," tab:"),
                   p("The exsiting model will be used to predict the origin of new sample."),
                  
                   p("* World map: new samples' prediction origins are plotted on the world map"),
                   br(),
                   
                   p(strong("Output")," tab:"),
                   p("The results of using the exsiting prediction model to predict the source coordinates of new sample(s)."),
                   p("For more introduction of output results, view the tutorial file in ",
                     a("mGPS interface Gitbub",href = "https://github.com/YaliZhang98/mGPS_interface"))
                   
                   ,h4("_______________________________________________"),
                   
                   ),
          tabPanel("Result Plot",
                   fluidRow( column(6, sliderInput("xlim_1", "longitude range:",
                                                   min = -165, max = 168,
                                                   value = c(-165,168))), 
                             column(6, sliderInput("ylim_1", "latitude range:",
                                                   min = -90, max = 90,
                                                   value = c(-90,90)))) ,
                   radioButtons("pull_1", "Whether pull points to land/marine",
                                choices = c("Pull to land" = "land_1",
                                            "Pull to waterbody" = "marine_1",
                                            "Default" = "none"),
                                selected = "none"), 
                   plotOutput(outputId = "predicted_map_1"),
                   downloadButton("downloadMap_1",
                                  "DownloadMap")
                   ),
          
          tabPanel("Output",
                   helpText("Here you can download the predicted original coordinates of samples in your uploaded file. 
                            (Reference microbiome locations are based on MetaSub project microbiome database)"),
                   
                   downloadButton("downloadData_1",
                                  "DownloadData")
                   ),
        )
      ),
    
      
      conditionalPanel(
        condition = "input.program == 'bb'",
        tabsetPanel(
          type = "tabs",
          tabPanel("HELP",
                   
                   h3("Build a new prediction model using mGPS and predict new samples"),
                   h4("_______________________________________________"),
                   
                   br(),
                   h4("Function description"),
                   br(),
                   p("In this mode, user can train the microbial origin prediction model based on 
                   the reference data set uploaded by the user. The constructed prediction 
                   model will be used to predict the new sample to be tested provided by the 
                   user and report the prediction result of the sample source. (If user want 
                   to visualize the accuracy of the model, please use function:",em("Build a new prediction model using mGPS", style = "color:purple"),")"),
                   h4("_______________________________________________"),
                   
                   br(),
                   
                   h4("Usage"),
                   br(),
                   p("In left side bar:",style = "font-size:16px"),
                   p("1. Select ",strong("Prediction program")," as ",em("Build a new prediction model using mGPS and predict new samples", style = "color:purple")), 
                   p("2. ",strong("Upload new sample(s) abundance file:")," Upload file (in .csv format) containing abundance data of new sample(s)."),  
                   
                   p("3. ",strong("Upload reference file(s):")," Upload data file(s) (in .csv format) containing microbial abundance data and metadata."),  
                   p("   In metadata, at least one locality (eg. continent, city) and coordinates (necessary) data columns should be included. The metadata and abundance data of the sample can be merged into one file (",em("Merged metadata and abundance data", style = "color:purple"),"), or uploaded as two files (",em("Separate metadata and abundance data", style = "color:purple"),")"),
                   p("When ",em("Separate metadata and abundance file", style = "color:purple")," is selected. ",strong("Merge column name in metadata/abundance file: "),"Input the header name of column which is the merged column in two files."),
                   p("4. ",strong("Enter the main locality level:")," Input the main locality target. It should same as that column header. (eg. city)"),
                   p("5. ",strong("Enter the locality hierarchy:")," The locality chain used in mGPS to construct the prediction model (same column headers). It should contain one or two locality information, latitude and longitude. Use ',' as the separator. (eg. continent,city,latitude,longitude)"),
                   p("6. ",strong("Columns range of abundance data:")," Input the columns range number of abundance data in the abundance/merged file. Use ':' as separator (eg 44:1000)"),
                   p("7. ",strong("Locality sample size cut off:"),em(" (Optional) "),"Remove locality whose sample size is less than a certain value. If checked, input the cut off number (eg. 8)"),
                   p("8. ",strong("Subsets in feature elimination:"),em(" (Optional) ")," Limit the number of features to a certain value. If unchecked, mGPS will find the optimal subset size of microbiome features. If checked, there are three types of input format: a.input the subsets size with separator as ',' (eg. 50,100,200,300); b. input the subsets size range with separator as '-' (eg. 50-300); c. input a single value."),
                   br(),
                   p("In ",strong("Result Plot")," tab:",style = "font-size:16px"),
                   p("9. ",strong("Change longitude/latitude range in output map"),em(" (Optional) ")),
                   p("10. ",strong("Whether pull points to land/marine:"),em(" (Optional) "),"If checked, predicted origin location will be pull to the nearest land/marine if predicted coordinates are out of the expected boarder."),
                   p(strong("Start the program:")," Click the",strong("Start"),"bar and then click the ",strong("Result Plot")," tab",style = "font-size:16px"),
                   p('Data processing, please wait while output files are being generated. When the prompt bar disappears you can see the results and download files.',style = "font-size:16px"),
                   
                   h4("_______________________________________________"),
                   
                   br(),
                   
                   h4("Result plot and output"),
                   br(),
                   p(strong("Result Plot")," tab:"),
                   p("The reference datasets will be used to construct a origin prediction model by mGPS.
                     Then this model will be used to predict origin of new samples. The prediction origin map and
                     prediction accuracy plot can be downloaded by bottons under the plots."),
                   
                   p("* World map: samples' prediction origins are plotted on the world map"),
                   br(),
                   
                   p(strong("Output")," tab:"),
                   p("The results of using the constructed prediction model to predict the source coordinates of the new samples. In addition, the constructed prediction model can be downloaded (In Rda format)."),
                   p("For more introduction of output results, view the tutorial file in ",
                     a("mGPS interface Gitbub",href = "https://github.com/YaliZhang98/mGPS_interface"))
          ,h4("_______________________________________________"),

                   ),
                                                                                           
                     
          tabPanel("Result Plot",
                   fluidRow( column(6, sliderInput("xlim_2", "longitude range:",
                                                   min = -165, max = 168,
                                                   value = c(-165,168))), 
                             column(6, sliderInput("ylim_2", "latitude range:",
                                                   min = -90, max = 90,
                                                   value = c(-90,90)))) ,
                   radioButtons("pull_2", "Whether pull points to land/marine",
                                choices = c("Pull to land" = "land_2",
                                            "Pull to waterbody" = "marine_2",
                                            "Default" = "none_2"),
                                selected = "none_2"),   
                   plotOutput(outputId = "predicted_map_2"),
                   downloadButton("downloadMap_2",
                                  "DownloadMap")),
          
          tabPanel("Output",
                   helpText("Predicted original coordinates of samples can be downloaded. Prediction model in RDA format can be downloaded"),
                   downloadButton("downloadData_2",
                                  "Download prediction data"),
                  
                   downloadButton("download_featuresub_2",
                                  "Download feature subsets accuracy in feature elimination"),
                   downloadButton("download_optfeatures_2",
                                  "Download optimal features in prediction model"),
                   downloadButton("downloadModel_2",
                                  "DownloadModel")
                   )
        )
      ),
      
      conditionalPanel(
        condition = "input.program == 'cc'",
        tabsetPanel(
          type = "tabs",
          
          tabPanel("Welcome",
                   h2("Welcome to mGPS application"),
                   
                   br(),
                   p("This is a web program based on the mGPS application created by 
      Shiny. It can build a microbial origin prediction model and predict 
        the origin of microbes.", style = "font-size:16px"),
                   p("To learn more about mGPS, please visit:",
                     a(em("mGPS"),
                       href="https://github.com/eelhaik/mGPS"), style = "font-size:16px"),
                   br(),
                   h3("Function"),
                   br(),
                   p("This program contains three functions. These function can be performed by selected ",strong("Prediction program")," in the left side bar. The detail usage 
        will be introduced in ",strong("HELP")," tab for each function.", style = "font-size:16px"),
                   br(),
                   p("1. ",strong("Build a new prediction model using mGPS")),
                   p("In this mode, user can use the mGPS tool to build a microbial source prediction model based
                            on the microbial abundance data and coordinates data uploaded by the user.","To learn more about mGPS, please visit:",
                     a(em("mGPS"),
                       href="https://github.com/eelhaik/mGPS"),
                     p("2. ",strong("Build a new prediction model using mGPS and predict new samples")),
                     p("In this mode, user can train the microbial origin
        prediction model based on the reference data set uploaded by the user. The
        constructed prediction model will be used to predict the new sample to be tested provided
        by the user and report the prediction result of the sample source. (If user want
                                                                            to visualize the accuracy of the model, please use function:",em("Build a new prediction model using mGPS"),")"),
                     p("3. ",strong("Use an existing model to predict new samples")),
                     p("In this mode, user can predict new sample origin based on an existing prediction model. Model can be downloaded in ", strong("Output")," tab of function:",em("Build a new prediction model using mGPS", style = "color:purple"), "or ",em("Build a new prediction model using mGPS and predict new samples", style = "color:purple")),
                     br(),
                     p("For more detail introduction and examples, visit the ",
                       a("mGPS interface on Gitbub", 
                         href = "https://github.com/YaliZhang98/mGPS_interface"), style = "font-size:16px"),
                   )),
          
          tabPanel("HELP",
                   
                   h3("Build a new prediction model using mGPS"),
                   h4("_______________________________________________"),
                   
                   br(),
                   h4("Function description"),
                   br(),
                   p("This model can use the mGPS tool to build a microbial source prediction model based
                            on the microbial abundance data and coordinates data uploaded by the user."),
                   h4("_______________________________________________"),
                   
                   br(),
                   
                   h4("Usage"),
                   br(),
                   p("In left side bar:",style = "font-size:16px"),
                   p("1. Select ",strong("Prediction program")," as ",em("Build a new prediction model using mGPS", style = "color:purple")), 
                   p("2. ",strong("Input file(s):")," Upload data file(s) (in .csv format) containing microbial abundance data and metadata."),  
                   p("   In metadata, at least one locality (eg. continent, city) and coordinates (necessary) data columns should be included. The metadata and abundance data of the sample can be merged into one file (",em("Merged metadata and abundance data"),"), or uploaded as two files (",em("Separate metadata and abundance data"),")"),
                   p("When ",em("Separate metadata and abundance file")," is selected. ",strong("Merge column name in metadata/abundance file: "),"Input the header name of column which is the merged column in two files."),
                   p("3. ",strong("Enter the main locality level:")," Input the main locality target. It should same as that column header. (eg. city)"),
                  p("4. ",strong("Enter the locality hierarchy:")," The locality chain used in mGPS to construct the prediction model (same column headers). It should contain one or two locality information, latitude and longitude. Use ',' as the separator. (eg. continent,city,latitude,longitude)"),
                  p("5. ",strong("Columns range of abundance data:")," Input the columns range number of abundance data in the abundance/merged file. Use ':' as separator (eg 44:1000)"),
                  p("6. ",strong("Locality sample size cut off:"),em(" (Optional) "),"Remove locality whose sample size is less than a certain value. If checked, input the cut off number (eg. 8)"),
                  p("7. ",strong("Subsets in feature elimination:"),em(" (Optional) ")," Limit the number of features to a certain value. If unchecked, mGPS will find the optimal subset size of microbiome features. If checked, there are three types of input format: a.input the subsets size with separator as ',' (eg. 50,100,200,300); b. input the subsets size range with separator as '-' (eg. 50-300); c. input a single value."),
                  br(),
                  p("In ",strong("Result Plot")," tab:",style = "font-size:16px"),
                  p("8. ",strong("Change longitude/latitude range in output map"),em(" (Optional) ")),
                  p("9. ",strong("Whether pull points to land/marine:"),em(" (Optional) "),"If checked, predicted origin location will be pull to the nearest land/marine if predicted coordinates are out of the expected boarder."),
                  p(strong("Start the program:")," Click the",strong("Start"),"bar and then click the ",strong("Result Plot")," tab",style = "font-size:16px"),
                  p('Data processing, please wait while output files are being generated. When the prompt bar disappears you can see the results and download files.',style = "font-size:16px"),
                  h4("_______________________________________________"),
                  
                  br(),
                  
                  h4("Result plot and output"),
                  br(),
                  p(strong("Result Plot")," tab:"),
                  p("Show the accuracy
                            of the prediction model trained by the mGPS tool and based on the reference microbial database
                            uploaded by the user."),
                  p("The original database will be divided into 5 folds, and mGPS will use 4 of
                            these folds to train the model, and the resulting model will be used to predict the microbial source
                            of the remaining fold. Iteratively obtain the prediction result of the original database and compare
                            it with the actual location of the microorganism."
                  ),
                  
                  p("* World map: samples' prediction origins are plotted on the world map"),
                  p("* Accuracy bar plot: model built by mGPS accuracy per site for the original reference dataset. mGPS accuracy is shown per-site as the distances between the predicted and true sampling site for the reference samples. The average prediction accuracy across all samples with each population given equal weight is shown on the left."),
                  br(),
                  
                  p(strong("Output")," tab:"),
                  p("The results of using the constructed prediction model to predict the source coordinates of the original dataset samples. In addition, the constructed prediction model can be downloaded (In Rda format)."),
                  p("For more introduction of output results, view the tutorial file in ",
                    a("mGPS interface Gitbub",href = "https://github.com/YaliZhang98/mGPS_interface"))
                  ,h4("_______________________________________________"),
),
          
          tabPanel("Result Plot",
                   
                   fluidRow( column(6, sliderInput("xlim_3", "longitude range:",
                                                   min = -165, max = 168,
                                                   value = c(-165,168))), 
                             column(6, sliderInput("ylim_3", "latitude range:",
                                                   min = -90, max = 120,
                                                   value = c(-90,120)))),
                   
                   radioButtons("pull_3", "Whether pull points to land/marine",
                                choices = c("Pull to land" = "land_3",
                                            "Pull to waterbody" = "marine_3",
                                            "Default" = "none"),
                                selected = "none"),   
                   plotOutput(outputId = "predicted_map_3"),
                   downloadButton("downloadMap_3",
                                  "DownloadMap"),
                   plotOutput(outputId = "predicted_accuracy_3"),
                   downloadButton("downloadAccuracy_3",
                                  "DownloadPlot")
                   ),
          tabPanel("Output",
                   helpText("Predicted original coordinates of samples can be downloaded. Prediction model in RDA format can be downloaded"),
                   downloadButton("downloadData_3",
                                  "Download prediction data"),
                   downloadButton("download_featuresub_3",
                                  "Download feature subsets accuracy in feature elimination"),
                   downloadButton("download_optfeatures_3",
                                  "Download optimal features in prediction model"),
                   downloadButton("downloadModel",
                                  "DownloadModel"))
      
      )
    )
  )
))


# -------------------------------------------------------------------

# Server

options(shiny.maxRequestSize=200*1024^2)

server <- function(input,output){
    
      model_store <- reactive({
        req(input$model)
        path <- input$model$datapath
        validate(
          need('Rda' %in%  strsplit(input$model$datapath,".",fixed = T)[[1]], "The model file you upload should be in .Rda format.")
        )
        load(path)
          return(model_store)
      })
      
      test_f <- reactive({
        if(input$program == "aa" ){
          req(input$f_new_test_1)
          test_f <- read.csv(input$f_new_test_1$datapath, 
                             header = TRUE)
          print(input$f_new_test_1$datapath)
          validate(
            need('csv' %in%  strsplit(input$f_new_test_1$datapath,".",fixed = T)[[1]], "The sample abundance file you upload should be in.csv format.")
          )
        }
        if(input$program == "bb" ){
          req(input$f_new_test_2)
          test_f <- read.csv(input$f_new_test_2$datapath, 
                             header = TRUE)
          validate(
            need('csv' %in%  strsplit(input$f_new_test_2$datapath,".",fixed = T)[[1]], "The sample abundance file you upload should be in.csv format.")
          )}
          return(test_f)
      })
      
      MetasubDataPreds <- reactive({
          MetasubDataPreds <- MetaSub_prediction(model_store(),test_f())
          
        return(MetasubDataPreds)    
      })

    pull_MetasubDataPreds <- reactive({
      if (input$pull_1 == "land_1"){
        pull_MetasubDataPreds <- pull_land(MetasubDataPreds())
      }else if (input$pull_1 == "marine_1"){
        pull_MetasubDataPreds <- pull_marine(MetasubDataPreds())  
      }else{
        pull_MetasubDataPreds <- MetasubDataPreds()
      }
      return(pull_MetasubDataPreds)
    })
      
    observeEvent(input$acbutton_1,{
    output$predicted_map_1 <- renderPlot({
      x_ran <- input$xlim_1
      y_ran <- input$ylim_1
      plot_prediction(pull_MetasubDataPreds(),x_ran,y_ran)
    })
    
    output$downloadMap_1 <- downloadHandler(
      
      filename = function(){
        return("Predicted_map.png")
      },
      content = function(file){
        
        MetasubDataPreds <- pull_MetasubDataPreds()
        
        x_ran <- input$xlim_1
        y_ran <- input$ylim_1
        
        png(file, width = 13,height = 8, units = 'in', res = 600)
        
        par(mai=c(1,1,0.5,0.5))
        
        map <- rworldmap::getMap(resolution = "coarse")
        
        plot(map, xlim = x_ran,ylim = y_ran, col = "grey",border = "darkgrey", xlab = "", ylab = '', bg = "lightskyblue1")
        title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
        
        find_lats <- MetasubDataPreds$latPred
        find_longs <- MetasubDataPreds$longPred
        
        points(find_longs,find_lats, type = "p",col = "purple", pch = "+", cex = 1.3)

        map.axes()
        dev.off()
      })
    
    
    
    output$downloadData_1 <- downloadHandler(
      filename = function(){
        return("Prediction_results.csv")
      },
      content = function(file){
        write.csv(pull_MetasubDataPreds(),file)
      }
    )
  })
    

    observeEvent(input$acbutton_2,{
      output$predicted_map_2 <- renderPlot({
        df <- prediction_output()[[1]]
        x_ran <- input$xlim_2
        y_ran <- input$ylim_2
        plot_prediction(df,x_ran,y_ran)
      })
      
      output$downloadMap_2 <- downloadHandler(
        
        filename = function(){
          return("Predicted_map.png")
        },
        content = function(file){

          MetasubDataPreds <- prediction_output()[[1]]
         
          x_ran <- input$xlim_2
          y_ran <- input$ylim_2
          
          png(file, width = 13,height = 8, units = 'in', res = 600)
          
          
          par(mai=c(1,1,0.5,0.5))
          
          map <- rworldmap::getMap(resolution = "coarse")
          
          plot(map, xlim = x_ran,ylim = y_ran, col = "grey",border = "darkgrey", xlab = "", ylab = '', bg = "lightskyblue1")
          title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
          
          find_lats <- MetasubDataPreds$latPred
          find_longs <- MetasubDataPreds$longPred
          
          points(find_longs,find_lats, type = "p",col = "purple", pch = "+", cex = 1.3)

          map.axes()
          dev.off()
        })
      
      output$downloadData_2 <- downloadHandler(
        
        filename = function(){
          return("Predicted_origin.csv")
        },
        content = function(file){
          write.csv(prediction_output()[[1]],file)
        })
      
      output$download_optfeatures_2 <- downloadHandler(
        
        filename = function(){
          return("Optimal_features.csv")
        },
        content = function(file){
          write.csv(v(),file)
        })
      
      output$download_featuresub_2 <- downloadHandler(
        
        filename = function(){
          return("Features_subsets_accuracy.csv")
        },
        content = function(file){
          write.csv(git_subset(),file)
        })
      
      output$downloadModel_2 <- downloadHandler(
        filename = function(){
          return("Prediction_model.Rda")
        },
        content = function(file){
          output_model <- prediction_output()[[2]]
          save(output_model,file=file)
        })
    })

# function b and c
    
    by_x_in <- reactive({
      if(input$program == "bb" ){
        req(input$by_x_2)
        by_x_in <- input$by_x_2
        
      }
      if(input$program == "cc" ){
        req(input$by_x_3)
        by_x_in <- input$by_x_3
      }
      return(by_x_in)
    })
    
    by_y_in <- reactive({
      if(input$program == "bb" ){
        req(input$by_y_2)
        by_y_in <- input$by_y_2
      }
      if(input$program == "cc" ){
        req(input$by_y_3)
        by_y_in <- input$by_y_3}
      return(by_y_in)
    })
    
    classTarget_in <- reactive({
      if(input$program == "bb" ){
        req(input$target_2)
        classTarget_in <- input$target_2}
      if(input$program == "cc" ){
        req(input$target_3)
        classTarget_in <- input$target_3}
      return(classTarget_in)
    })
    
    hierarchy_in <- reactive({
      if(input$program == "bb" ){
        req(input$hierarchy_2)
        text <- input$hierarchy_2
        hierarchy_in <- strsplit(text,",")[[1]]}
      if(input$program == "cc" ){
        req(input$hierarchy_3)
        text <- input$hierarchy_3
        hierarchy_in <- strsplit(text,",")[[1]]}
      return(hierarchy_in)
    })
    
    
    subsets_in <- reactive({
      if(input$program == "bb" ){
        if (input$subsets_2 == T){
          req(input$subsets_value_2)
          text <- input$subsets_value_2
          
          if ('-' %in% strsplit(text,"")[[1]]){
            subsets_range <- as.numeric(strsplit(text,"-")[[1]])
            
            subsets_in <- c(subsets_range[1],
                            round((subsets_range[2]-subsets_range[1])/5+subsets_range[1]),
                            round((subsets_range[2]-subsets_range[1])/5*2+subsets_range[1]),
                            round((subsets_range[2]-subsets_range[1])/5*3+subsets_range[1]),
                            round((subsets_range[2]-subsets_range[1])/5*4+subsets_range[1]),
                            subsets_range[2])
            
          }else{
            subsets_in <- as.numeric(strsplit(text,",")[[1]])
          }
          
          
         }else{
          subsets_in <- NULL
          }}
      if(input$program == "cc" ){
        if (input$subsets_3 == T){
          req(input$subsets_value_3)
          text <- input$subsets_value_3
          if ('-' %in% strsplit(text,"")[[1]]){
            subsets_range <- as.numeric(strsplit(text,"-")[[1]])
            
            subsets_in <- c(subsets_range[1],
                            round((subsets_range[2]-subsets_range[1])/5+subsets_range[1]),
                            round((subsets_range[2]-subsets_range[1])/5*2+subsets_range[1]),
                            round((subsets_range[2]-subsets_range[1])/5*3+subsets_range[1]),
                            round((subsets_range[2]-subsets_range[1])/5*4+subsets_range[1]),
                            subsets_range[2])
            
          }else{
            subsets_in <- as.numeric(strsplit(text,",")[[1]])
          }
        }else{
          subsets_in <- NULL
        }}
      return(subsets_in)
    })
    
    remove_small <- reactive({
      if(input$remove_small_2 == T ){
        req(input$remove_small_value_2)
        remove_small <- input$remove_small_value_2
        remove_small <- as.numeric(remove_small)
      }else{
        remove_small <- 0
      }
      
      if(input$remove_small_3 == T ){
        req(input$remove_small_value_3)
        remove_small <- input$remove_small_value_3
        remove_small <- as.numeric(remove_small)
      }else{
        remove_small <- 0
      }
      return(remove_small)
    })
    
    abundance_r <- reactive({
      if(input$program == "bb" ){
        req(input$abundance_range_2)
        text <- input$abundance_range_2
        range_1 <- as.numeric(strsplit(text,":")[[1]][1])
        range_2 <- as.numeric(strsplit(text,":")[[1]][2])
        abundance_r <- list(range_1,range_2)}
      if(input$program == "cc" ){
        req(input$abundance_range_3)
        text <- input$abundance_range_3
        range_1 <- as.numeric(strsplit(text,":")[[1]][1])
        range_2 <- as.numeric(strsplit(text,":")[[1]][2])
        abundance_r <- list(range_1,range_2)}
      return(abundance_r)
    })
    
    
    
      train_f <- reactive({
        
        if(input$program == "bb" ){
          
          if(input$file_num_2 == "one_file" ){
            
            req(input$f_new_train_2)
            train_f <- read.csv(input$f_new_train_2$datapath, 
                                header = TRUE)
            validate(
              need('csv' %in%  strsplit(input$f_new_train_2$datapath,".",fixed = T)[[1]], "The training file you upload should be in.csv format.")
            )
          
            return(train_f)
          }
          
          if (input$file_num_2 == "two_file"){
            
            req(input$f_metadata_2)
            req(input$f_abundance_2)
            data_metadata <- read.csv(input$f_metadata_2$datapath, 
                                      header = TRUE)
            validate(
              need('csv' %in%  strsplit(input$f_metadata_2$datapath,".",fixed = T)[[1]], "The metadata file you upload should be in.csv format.")
            )
            data_abundance <- read.csv(input$f_abundance_2$datapath,
                                       header = T)
            validate(
              need('csv' %in%  strsplit(input$f_abundance_2$datapath,".",fixed = T)[[1]], "The abundance file you upload should be in.csv format.")
            )
            data_abundance <- unique(data_abundance)
            by_x <- by_x_in()
            by_y <- by_y_in()
            
            metasub_data <- merge(data_metadata,data_abundance,by.x=by_x,by.y=by_y) # #merge bacterial and meta data

            
            train_f <- metasub_data
            return(train_f)
          }
        }
        
        if(input$program == "cc" ){
        
        if(input$file_num_3 == "one_file" ){
          
          req(input$f_new_train_3)
          train_f <- read.csv(input$f_new_train_3$datapath, 
                              header = TRUE)
          validate(
            need('csv' %in%  strsplit(input$f_new_train_3$datapath,".",fixed = T)[[1]], "The training file you upload should be in.csv format.")
          )
         
          return(train_f)
        }
        
        if (input$file_num_3 == "two_file"){
          
          req(input$f_metadata_3)
          req(input$f_abundance_3)
          data_metadata <- read.csv(input$f_metadata_3$datapath, 
                                    header = TRUE)
          validate(
            need('csv' %in%  strsplit(input$f_metadata_3$datapath,".",fixed = T)[[1]], "The metadata file you upload should be in.csv format.")
          )
          data_abundance <- read.csv(input$f_abundance_3$datapath,
                                     header = T)
          validate(
            need('csv' %in%  strsplit(input$f_abundance_3$datapath,".",fixed = T)[[1]], "The abundance file you upload should be in.csv format.")
          )
          data_abundance <- unique(data_abundance)
          metasub_data <- merge(data_metadata,data_abundance,by.x=by_x_in(),by.y=by_y_in()) # #merge bacterial and meta data
          
          
          train_f <- metasub_data
          return(train_f)
        }
        }
      })
     
      metasub_data <- reactive({
          metasub_data <- data_preprocess_f(train_f(),classTarget_in(),hierarchy_in(),remove_small())
          return(metasub_data)
      })
      
      featureElim <- reactive({
        range_1 <- abundance_r()[1]
        range_2 <- abundance_r()[2]
        
        if(input$program == "bb" ){
          
          if(input$file_num_2 == "two_file" ){
            req(input$f_metadata_2)
            data_metadata <- read.csv(input$f_metadata_2$datapath, 
                                      header = TRUE)
            validate(
              need('csv' %in%  strsplit(input$f_metadata_2$datapath,".",fixed = T)[[1]], "The metadata file you upload should be in.csv format.")
            )
            range_1 <- dim(data_metadata)[2] + as.numeric(range_1) - 1
            range_2 <- dim(data_metadata)[2] + as.numeric(range_2) - 1
            
          }}
        if(input$program == "cc" ){
          
          if(input$file_num_3 == "two_file" ){
            req(input$f_metadata_3)
            data_metadata <- read.csv(input$f_metadata_3$datapath, 
                                      header = TRUE)
            validate(
              need('csv' %in%  strsplit(input$f_metadata_3$datapath,".",fixed = T)[[1]], "The metadata file you upload should be in.csv format.")
            )
            range_1 <- dim(data_metadata)[2] + as.numeric(range_1) - 1
            range_2 <- dim(data_metadata)[2] + as.numeric(range_2) - 1
            
          }}
        
        featureElim <- featureElimination_f(metasub_data(),classTarget_in(),range_1,range_2, subsets_in()) # 44:3712 for MetaSub
        return(featureElim)
      })
      
      prediction_output_before <- reactive({
        if(input$program == "bb" ){
          test_file <- test_f()
          prediction_output_before <- New_model_prediction(train_f(),test_f(),classTarget_in(),featureElim()$optVariables,hierarchy_in())
        }
        if(input$program == "cc" ){
          prediction_output_before <- model_accuracy_f(metasub_data(),featureElim()$optVariables,classTarget_in(),hierarchy_in())
        }
        return(prediction_output_before)
      })
      
      prediction_output <- reactive({
        prediction_output <- prediction_output_before()
        if(input$program == "bb" ){
          if (input$pull_2 == "marine_2"){
            df <- prediction_output[[1]]
            df_new <- pull_marine(df)
            prediction_output[[1]] <- df_new
          }else if(input$pull_2 == "land_2"){
            df <- prediction_output[[1]]
            df_new <- pull_land(df)
            prediction_output[[1]] <- df_new
          }
        }
        if(input$program == "cc" ){
          if (input$pull_3 == "marine_3"){
            df <- prediction_output[[1]]
            df_new <- pull_marine(df)
            prediction_output[[1]] <- df_new
          }else if(input$pull_3 == "land_3"){
            df <- prediction_output[[1]]
            df_new <- pull_land(df)
            prediction_output[[1]] <- df_new
          }
        }
        return(prediction_output)
      })
      
      
      v <- reactive({
        text <- featureElim()
        v <- v_f(text)
        return(v)
      })
      
      
      git_subset <- reactive({
        text <- featureElim()
        git_subset <- git_subset_f(text)
        return(git_subset)
      })
      
      
      observeEvent(input$acbutton_3,{
        output$predicted_map_3 <- renderPlot({
          x_ran <- input$xlim_3
          y_ran <- input$ylim_3
          df <- prediction_output()[[1]]
          hierarchy_in <- hierarchy_in()
          classTarget_in <- classTarget_in()
          
          plot_map(df,hierarchy_in,classTarget_in,x_ran,y_ran)
        })
        
        output$predicted_accuracy_3 <- renderPlot({
          classTarget_in <- classTarget_in()
          plot_accuracy(prediction_output()[[1]],classTarget_in)
        })
        
        output$downloadData_3 <- downloadHandler(
          filename = function(){
            return("Predicted_origin.csv")
          },
          content = function(file){
            write.csv(prediction_output()[[1]],file)
          })
        
        
        output$download_optfeatures_3 <- downloadHandler(
          filename = function(){
            return("Optimal_features.csv")
          },
          content = function(file){
            write.csv(v(),file)
          })
        
        output$download_featuresub_3 <- downloadHandler(
          filename = function(){
            return("Features_subsets_accuracy.csv")
          },
          content = function(file){
            write.csv(git_subset(),file)
          })
      
        output$downloadMap_3 <- downloadHandler(
          filename = function(){
            return("Prediction_map.png")
          },
          content = function(file){

            MetasubDataPreds <- prediction_output()[[1]]
            hierarchy_in <- hierarchy_in()
            classTarget_in <- classTarget_in()
            x_ran <- input$xlim_3
            y_ran <- input$ylim_3
            
            png(file, width = 13,height = 8, units = 'in', res = 600)

            par(mai=c(2,1,0.5,0.5), mar=par()$mar+c(3,0,0,0))
            
            map <- rworldmap::getMap(resolution = "coarse")
            palette <-c( "darkorchid4","gold2","dodgerblue3","brown","orangered2","mediumspringgreen","deeppink2")
            
            plot(map, xlim = x_ran,ylim = y_ran,col = "grey",border = "darkgrey", xlab = "", ylab = '', bg = "lightskyblue1")
            title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
            
            for (i in 1:(length(hierarchy_in)-2)){
              MetasubDataPreds[,hierarchy_in[i]] <- factor(MetasubDataPreds[,hierarchy_in[i]])
            }
            
            MetasubDataPreds$cityPred <- factor(MetasubDataPreds$cityPred)
            
            # get the coordination of continent
            continent_list <- c("east_asia","europe","middle_east","north_america","oceania","south_america", "sub_saharan_africa",'Australia','North_America','Europe','Africa','South_america','Asia')
            continent_lats <- c(55,69,8,40,-40,-10,-5,-40,40,69,-5,-10,55)
            continent_longs <- c(125,0,60,-130,140,-80,5,140,-130,0,5,-80,125)
            continent_position <- data.frame(cbind(continent_list,continent_lats,continent_longs))
            
            
            label_continent <- c()
            flag <- 0
            for ( i in 1:length(levels(MetasubDataPreds[,hierarchy_in[1]]))){
              this_continent <- levels(MetasubDataPreds[,hierarchy_in[1]])[i]
              label_continent <- c(label_continent,this_continent)
              find_lats <- MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,][,"latPred"]
              find_longs <- MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,][,"longPred"]
              
              #plot predicted co-ordinates
              points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.2,xlim = c(-165,168))
              
              #plot city prediction accuravy by continent as pies
              
              correctly_pred <-  mean(as.numeric(as.character(MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,"cityPred"])== 
                                                                as.character(MetasubDataPreds[MetasubDataPreds[,hierarchy_in[1]] == this_continent,classTarget_in]))) 
              
              incorrectly_pred <- (1 - correctly_pred)
              
              if (this_continent %in% continent_position$continent_list){
                add.pie(z = c(correctly_pred, incorrectly_pred), x = as.numeric(continent_position[continent_position$continent_list==this_continent,3]), y = as.numeric(continent_position[continent_list==this_continent,2])
                        ,edges=200,
                        radius=10,
                        col=c(palette[i],"black") , labels = ""
                )
              }else{
                
                add.pie(z = c(correctly_pred, incorrectly_pred), x = as.numeric(-150+flag), y = 105
                        ,edges=200,
                        radius=10,
                        col=c(palette[i],"black") , labels = ""
                )
                flag <- flag + 25
              }
            }

            legend(xpd=T,'bottom',inset=c(0,-0.25), label_continent, pch = "+", col = palette[1:length(label_continent)], cex = 1,n=5)
            
            box( col = 'black')
            
            map.axes()
            dev.off()
          })
        
        
        
        
        output$downloadAccuracy_3 <- downloadHandler(
          filename = function(){
            return("Prediction_accuracy.png")
          },
          content = function(file){

            MetasubDataPreds <- prediction_output()[[1]]

            classTarget_in <- classTarget_in()
            
            MetasubDataPreds[,classTarget_in] <- factor(MetasubDataPreds[,classTarget_in])
            
            bar_df1 <- data.frame(row.names = c("Overall",levels(MetasubDataPreds[,classTarget_in])))
            
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 100)
              bar_df1[i+1,"0 - 100km"] <- prop
              
              overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] < 100)
              bar_df1[ 1,"0 - 100km"] <- overall_prop
            }
            
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 100 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 500)
              bar_df1[i+1,"100 - 500km"] <- prop
              
              overall_prop <-mean(MetasubDataPreds[,"Distance_from_origin"] > 100 & MetasubDataPreds[,"Distance_from_origin"] < 500)
              bar_df1[ 1,"100 - 500km"] <- overall_prop
            }
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 500 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 1000)
              bar_df1[i+1,"500 - 1000km"] <- prop
              
              overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 500 & MetasubDataPreds[,"Distance_from_origin"] < 1000)
              bar_df1[ 1,"500 - 1000km"] <- overall_prop
            }
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 1000 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 2000)
              bar_df1[i+1,"1000 - 2000km"] <- prop
              
              overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 1000 & MetasubDataPreds[,"Distance_from_origin"] < 2000)
              bar_df1[ 1,"1000 - 2000km"] <- overall_prop
            }
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 2000 & MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] < 3000)
              bar_df1[i+1,"2000 - 3000km"] <- prop
              
              overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 2000 & MetasubDataPreds[,"Distance_from_origin"] < 3000)
              bar_df1[1,"2000 - 3000km"] <- overall_prop
            }
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              prop <- mean(MetasubDataPreds[MetasubDataPreds[,classTarget_in] == this_city,][,"Distance_from_origin"] > 3000 )
              bar_df1[i+1,"> 3000km"] <- prop
              
              overall_prop <- mean(MetasubDataPreds[,"Distance_from_origin"] > 3000)
              bar_df1[ 1,"> 3000km"] <- overall_prop
            }
            
            size <- c()
            for (i in 1: length(levels(MetasubDataPreds[,classTarget_in]))){
              
              this_city <- levels(MetasubDataPreds[,classTarget_in])[i]
              size[i] <- length(which(MetasubDataPreds[,classTarget_in] == this_city))
            }
            

            png(file, width = 13,height = 8, units = 'in', res = 600)
            par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)

            bp <- barplot(t(bar_df1*100), space = 0,col=c("lightyellow","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), 
                          names.arg=c("Overall",paste0(levels(MetasubDataPreds[,classTarget_in]),"  (",size,")"), axes = FALSE) , 
                          las =2, cex.names=.6, ylab = "", axisnames = F, axes = F)
            axis(side =2, pos = 0)
            mtext(text = c("Overall",paste0(levels(MetasubDataPreds[,classTarget_in])," (",size,")")), side = 1, at = bp, line = 0, padj = 1, cex = 0.7)
            title(ylab="Proportion of sample predictions %", mgp=c(2,1,0),cex.lab=1)
            legend("topright",inset = c(-0.12,0.4), rev(c(colnames(bar_df1))), 
                   fill = rev(c("lightyellow","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) ,
                   bty = 1, cex = 0.8)
           
            par(mar=c(5, 4, 4, 2) + 0.1)
           
            dev.off()
          })
        
        output$downloadModel <- downloadHandler(
          filename = function(){
            return("Prediction_model.Rda")
          },
          content = function(file){
            output_model <- prediction_output()[[2]]
            save(output_model,file=file)
          })
      })
}


shinyApp(ui=ui,server=server)






