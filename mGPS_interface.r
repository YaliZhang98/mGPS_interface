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

##Data transformation 
data_normalise <- function(df) {
  
  return(df/rowSums(df))
  
}




##Feature selection algorithm
species_select <-
  function(x,
           y,
           remove_correlated = T,
           subsets = NULL,
           cores = 1) {
    doParallel::registerDoParallel(cores)
    
    #!!!!! transfer value to factor
    y <- factor(y)
    
    message("feature elimination ....................")
    
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
    
    message("feature elimination has been done ....................")
    
    return(featureElimination)
    
    
  }

##Main mGPS algorithm 
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
        
        message("Continent training .........................")
        


        Xgb_region <- train(x = training[,variables],y = training[,hierarchy[1]],
                            method = "xgbTree",
                            trControl = trControlClass,
                            tuneGrid = tune_grid,
                            nthread = nthread)
        
        message("Continent training has been done .........................")
        
        l1_train <- data.frame(training[,c(variables)],Xgb_region[["pred"]][order(Xgb_region$pred$rowIndex),levels(training[,hierarchy[1]]) ])

      }else{
        
        l1_train <- training[,c(variables)]
        
      }
      

      message("City training .........................")
      
      Xgb_class <- train(x = l1_train,y = training[,classTarget],
                         method = "xgbTree",
                         trControl = trControlClass,
                         tuneGrid = tune_grid,
                         nthread = nthread)
      
      message("City training has been done .........................")
      
      l2_train <- data.frame(l1_train,Xgb_class[["pred"]][order(Xgb_class$pred$rowIndex),levels(training[,classTarget]) ])
      
      message("latitude training .........................")
      
      Xgb_latitude <- train(x = l2_train ,y = training[,"latitude"],
                            method = "xgbTree",
                            trControl = trControl,
                            tuneGrid = tune_grid,
                            nthread = nthread)
      message("Latitude training has been done .........................")
      
      l3_train <- data.frame(l2_train, "latPred" = Xgb_latitude[["pred"]][order(Xgb_latitude$pred$rowIndex),"pred" ])
      
      message("longitude training .........................")
      
      Xgb_longitude <- train(x = l3_train ,y = training[,"longitude"],
                             method = "xgbTree",
                             trControl = trControl,
                             tuneGrid = tune_grid,
                             nthread = nthread)
      
      message("Longitude training has been done .........................")
      
    }
    #check for test set, return trained model if no test set. 
    
      
      model <- function(test,variables){
        regProbs <- predict(Xgb_region, newdata = test[,variables],type ="prob")
        
        l1_test <- data.frame(test[,variables], regProbs)
        
        classPred <- predict(Xgb_class, newdata = l1_test)
        classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
        
        l2_test <-  data.frame(l1_test, classProbs) 
        latPred <- predict(Xgb_latitude, newdata = l2_test)
        
        l3_test <- data.frame(l2_test, latPred)
        longPred <- predict(Xgb_longitude, newdata = l3_test)
        return(list(classPred, latPred, longPred))
        
      }

     
      message("Generating predictions")
      #generate mGPS predictions for test set
      
      
      
      #########!!!!!!!!!!!!!!!!!!ZYL
      features <- variables
      newdata_in <- as.data.frame(setNames(replicate(length(features),numeric(0), simplify = F), features))
      for (f in features){
        if (f %in% colnames(testing)){
          newdata_in[c(1:dim(testing)[1]),f] <- testing[,f]
        }else{
          newdata_in[c(1:dim(testing)[1]),f] <- data.frame(0)
        }
      }
      
      write.csv(testing,file = "Outputs/tesing_in.csv")
      write.csv(newdata_in,file = "Outputs/newdata_in_tesing.csv")
      
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
      
      message("Prediction has been done ..........................")
      
      #adjust out of bounds predictions
      longPred[longPred > 180] <- 180
      longPred[longPred < -180] <- -180
      latPred[latPred > 90] <- 90
      latPred[latPred < -90] <- -90
      #Pull to nearest coastline if provided
      
      find_coast <- function(long, lat) {
        distances_from_coastline <-
          sp::spDistsN1(coast, c(long, lat), longlat = TRUE)
        
        closest_point <-  which.min(distances_from_coastline)
        new_coords <- coast[closest_point,]
        
        return(new_coords)
        
      }
      
      if (!is.null(coast)) {
        toAdjust <-
          which(is.na(maps::map.where(database = "world", longPred, latPred)))
        
        adjusted <-
          mapply(find_coast, long = longPred[toAdjust], lat = latPred[toAdjust])
        
        
        longPred[toAdjust] <- adjusted[1,]
        latPred[toAdjust] <- adjusted[2,]
        
        
      }
      
      model_store <- list(Xgb_region,Xgb_class,Xgb_latitude,Xgb_longitude,"model" = model)
      
      prediction_store <- list(classPred, latPred, longPred)
      
      return(list(model_store, prediction_store))
      
    # }
    
  }





#Import data sets 
#setwd( rprojroot::find_rstudio_root_file())

remove_control <- function(metasub_data,text){
  for (i in range(length(strsplit(text,";")[[1]]))){
    v <- strsplit(text,";")[[1]][i]
    col <- strsplit(v,":")[[1]][1]
    v_remove <-  strsplit(strsplit(v,":")[[1]][2],",")[[1]]
    
    if (col %in% colnames(metasub_data)){
      control_samples <- c( which(metasub_data[,col] %in% v_remove)) 
      if (length(control_samples) != 0){
        metasub_data <- droplevels(metasub_data[-c(control_samples), ])
      }
    }
  }
  write.csv(metasub_data,"Outputs/after_remove_control.csv")
  return(metasub_data)
}

data_preprocess_f <- function(train_f,target_in,hierarchy){

  #Remove control samples
  metasub_data <- droplevels(train_f)
  
  if(length(hierarchy) == 4){
    metasub_data[,hierarchy[1]] <- factor(metasub_data[,hierarchy[1]])
    metasub_data[,hierarchy[2]] <- factor(metasub_data[,hierarchy[2]])
  }else{
    metasub_data[,hierarchy[1]] <- factor(metasub_data[,hierarchy[1]])
  }
  #remove sparse samples locations and dubiously labelled samples. 
  
  
  small_cities <-  names(which(summary(metasub_data[,target_in]) < 8))
    
  remove_samples <- which(metasub_data[,target_in] %in%  c("antarctica", small_cities))
  

  
  if (length(remove_samples) != 0 ){
    metasub_data <- droplevels(metasub_data[-c(remove_samples), ])
  }

  #Correction of identified misslabelling of data 
  
  for (i in range(length(hierarchy))){
    empty <- which(metasub_data[,hierarchy[i]] == "" | is.na(metasub_data[,hierarchy[i]]))

    if (length(empty) != 0  ){
      metasub_data <- metasub_data[-c(empty),]
    }
  
  }

 

  metasub_data[,hierarchy[1]] <- droplevels(metasub_data[,hierarchy[1]])
  
  if(length(hierarchy) == 4){
    metasub_data[,hierarchy[2]] <- droplevels(metasub_data[,hierarchy[2]])
  }

  message("Data has been pre-processed")
  write.csv(metasub_data,file = "Outputs/after_preprocess_dataset.csv")
  return(metasub_data)
}



featureElimination_f <- function(metasub_data,classTarget_in,range_1,range_2,subsets_in){

  ### Find GITs ####
  range_1 <- as.numeric(range_1)
  range_2 <- as.numeric(range_2)
 
  featureElim <- species_select(x = metasub_data[, c(range_1:range_2)],y = metasub_data[,classTarget_in],remove_correlated = F, subsets_in ,cores = 8)

  return(featureElim)
}

v_f <- function(featureElim){
  optVars <- featureElim$optVariables
  v <- varImp(featureElim$fit, type = 1, scale = F)
  v[,"taxa"] <- row.names(v)
  v <- v[order(v$Overall,decreasing = T),]
  write.csv(v, file = "Outputs/optimal_features.csv")
  return(v)
}
  
git_subset_f <- function(featureElim){
  git_subset <- data.frame("n_vars" = featureElim$results$Variables, "accuracy" = featureElim$results$Accuracy)
  write.csv(git_subset,file = "Outputs/features_subset_accuracy.csv")
  return(git_subset)
}
  
  
 
coastlines_f <- function(){
  ### Make predictions using mGPS ###
  coastlines <- cbind("x"  = maps::SpatialLines2map(rworldmap::coastsCoarse)$x ,"y" =maps::SpatialLines2map(rworldmap::coastsCoarse)$y)
  coastlines <- coastlines[complete.cases(coastlines),]
  coastlines <- coastlines[coastlines[,1] < 180 ,]
  
  return(coastlines)
}

find_coast <- function(long, lat,coastlines) {
  
  distances_from_coastline <-
    sp::spDistsN1(coastlines, c(long, lat), longlat = TRUE)
  
  closest_point <-  which.min(distances_from_coastline)
  new_coords <- coastlines[closest_point,]
  
  return(new_coords)
}

  #generate 5 stratified folds for test predictions.

# Accuracy of model detected on original dataset
model_accuracy_f <- function(metasub_data,optVars,classTarget_in,hierarchy_in,coastlines){
  set.seed(18)
 
  trainFolds <-  createFolds(metasub_data[,classTarget_in], k = 5, returnTrain = T)
  
 
  
  GeoPreds <- list()
  
  #iteratively train the model on each of the 5 training folds and generate predictions using the coresponding test fold.
  for (i in 1:5){
    
    print(i)
    print("-------------------------------------------")
    
    train <- metasub_data[trainFolds[[i]],]
    test <- metasub_data[-trainFolds[[i]],]
    
    # testPreds <-mGPS(training = train, testing = test, classTarget = "city",variables = optVars,nthread = 8,hierarchy = c('continent','city','latitude','longitude'), coast=coastlines)
    testPreds <-mGPS(training = train, testing = test, classTarget = classTarget_in,variables = optVars,nthread = 8,hierarchy = hierarchy_in, coast=coastlines)
    
    GeoPreds[[i]] <- testPreds[[2]]

    print(i)
    print("has been done -------------------------------------------")
    
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
 
  for (i in 1:nrow(MetasubDataPreds)){
    MetasubDataPreds[i,"Distance_from_origin"] <- geosphere::distm(c(MetasubDataPreds[i,"longPred"],MetasubDataPreds[i,"latPred"]), c(MetasubDataPreds[i,"longitude"],MetasubDataPreds[i,"latitude"]), fun = geosphere::distHaversine)/1000
  }
  
  message("Distance_from_origin has been done")
  
  write.csv(MetasubDataPreds,"Outputs/Prediction_results.csv")
  
  save(model_last_one,file="Outputs/Prediction_model.Rda")
  
  return(list(MetasubDataPreds,model_last_one)) 
}






# writing..................
New_model_prediction <- function(train,test,classTarget_in,optVars,hierarchy_in,coastlines){
  
  write.csv(test,"Outputs/New_model_prediction_test.csv")
  
  testPreds <-mGPS(training = train, testing = test, classTarget = classTarget_in,variables = optVars,nthread = 8,hierarchy = hierarchy_in, coast=coastlines)
  
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





feature_filter <- function(test_dataset,model){
  features <- (model$finalModel$feature_names)[1:200]
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
MetaSub_prediction <- function(model_store,test_dataset,coastlines){
  Xgb_region <- model_store[[1]]
  Xgb_class <- model_store[[2]]
  Xgb_latitude <- model_store[[3]]
  Xgb_longitude <- model_store[[4]]
  # prediction_model <- model_store[[5]]
  
  print("Prediction ..................................")
  
  model <- function(test_dataset,coastlines){
    
    newdata <- feature_filter(test_dataset, Xgb_class) # !!!!!!!!!!
    
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
    
    find_coast <- function(long, lat) {
      distances_from_coastline <-
        sp::spDistsN1(coastlines, c(long, lat), longlat = TRUE)
      
      closest_point <-  which.min(distances_from_coastline)
      new_coords <- coastlines[closest_point,]
      
      return(new_coords)
      
    }
    
    if (!is.null(coastlines)) {
      
      toAdjust <-
        which(is.na(maps::map.where(database = "world", longPred, latPred)))
      
      
      if (length(toAdjust) != 0){
        
        adjusted <-
          mapply(find_coast, long = longPred[toAdjust], lat = latPred[toAdjust])
        
        longPred[toAdjust] <- adjusted[1,]
        latPred[toAdjust] <- adjusted[2,]
      }
    }
    
    return(list(classPred, latPred, longPred))
  }
  
  testPreds <- model(test_dataset,coastlines)
  
  DataPreds <- cbind(test_dataset , 
                     "cityPred"= testPreds[[1]], 
                     "latPred" = testPreds[[2]], 
                     "longPred" = testPreds[[3]] )
  
  
  write.csv(DataPreds,"Outputs/Predicted_result.csv")
  
  message("Prediction result file has been built")
  
  return(DataPreds)
}


# Plot predicted location of samples in training dataset. It can show how many samples are located into their continent.
plot_map <- function(MetasubDataPreds){
  
  par(mai=c(1,1,0.5,0.5))
  
  map <- rworldmap::getMap(resolution = "coarse")
  palette <-c( "darkorchid4","gold2","dodgerblue3","brown","orangered2","mediumspringgreen","deeppink2")
  plot(map, xlim = c(-165,168), col = "grey",border = "darkgrey", xlab = "", ylab = '', bg = "lightskyblue1")
  title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)
  #find coord preds by region
  
 
  MetasubDataPreds$continent <- factor(MetasubDataPreds$continent)
  MetasubDataPreds$city <- factor(MetasubDataPreds$city)
  MetasubDataPreds$cityPred <- factor(MetasubDataPreds$cityPred)
  
  # get the coordination of continent
  continent_list <- c("east_asia","europe","middle_east","north_america","oceania","south_america", "sub_saharan_africa")
  continent_lats <- c(55,69,8,40,-40,-10,-5)
  continent_longs <- c(125,0,60,-130,140,-80,5)
  continent_position <- data.frame(cbind(continent_list,continent_lats,continent_longs))
  
  
  
  for ( i in 1:length(levels(MetasubDataPreds$continent))){
    this_continent <- levels(MetasubDataPreds$continent)[i]

    find_lats <- MetasubDataPreds[MetasubDataPreds[,"continent"] == this_continent,][,"latPred"]
    find_longs <- MetasubDataPreds[MetasubDataPreds[,"continent"] == this_continent,][,"longPred"]
    
    #plot predicted co-ordinates
    points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.2,xlim = c(-165,168))
    
    #plot city prediction accuravy by continent as pies
    correctly_pred <-  mean(MetasubDataPreds[MetasubDataPreds$continent == this_continent,"cityPred"]== 
                              MetasubDataPreds[MetasubDataPreds$continent == this_continent,"city"]) 
    
    incorrectly_pred <- (1 - correctly_pred) 
    
    
    
    add.pie(z = c(correctly_pred, incorrectly_pred), x = as.numeric(continent_position[continent_position$continent_list==this_continent,3]), y = as.numeric(continent_position[continent_list==this_continent,2])
            ,edges=200,
            radius=10,
            col=c(palette[i],"black") , labels = ""
    )
  }
  
  #Plot city sampling locations
  map.axes(cex.axis = 1.1)
  par(fig = c(0,0.4,0.0,0.5), new = T) 
  plot(map,xlim = c(-165,168), col = "grey", border = "darkgrey", bg ="lightskyblue1")
  label_continent <- c()
  for ( i in 1:length(levels(MetasubDataPreds$continent))){
    this_continent <- levels(MetasubDataPreds$continent)[i]
    label_continent <- c(label_continent,this_continent)
    find_lats <- MetasubDataPreds[MetasubDataPreds$continent == this_continent,]$city_latitude
    find_longs <- MetasubDataPreds[MetasubDataPreds$continent == this_continent,]$city_longitude
    
    points(find_longs, find_lats, col = palette[i], pch = 17, cex = 1,xlim = c(-165,168))
  }
  legend(-165,-15, label_continent, pch = 17, col = palette, cex = 0.5, bg ="lightskyblue1")
  box( col = 'black')
  
}

# Plot the prediction origin of sample on worldmap
plot_prediction <- function(MetasubDataPreds){
  par(mai=c(1,1,0.5,0.5))

  map <- rworldmap::getMap(resolution = "coarse")

  plot(map, xlim = c(-165,168), col = "grey",border = "darkgrey", xlab = "", ylab = '', bg = "lightskyblue1")
  title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)

  find_lats <- MetasubDataPreds$latPred
  find_longs <- MetasubDataPreds$longPred
  
  points(find_longs,find_lats, type = "p",col = "purple", pch = "+", cex = 1.2)
}

# --------------------------------------------------------------------

# Plot of accuracy

plot_accuracy <- function(MetasubDataPreds,classTarget_in){
  
  MetasubDataPreds[,classTarget_in] <- factor(MetasubDataPreds[,classTarget_in])
  
  for (i in 1:nrow(MetasubDataPreds)){
    MetasubDataPreds[i,"Distance_from_origin"] <- geosphere::distm(c(MetasubDataPreds[i,"longPred"],MetasubDataPreds[i,"latPred"]), c(MetasubDataPreds[i,"longitude"],MetasubDataPreds[i,"latitude"]), fun = geosphere::distHaversine)/1000
  }
  
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
  title(ylab="Proportion of sample predictions %", mgp=c(0,0,0),cex.lab=1)
  legend("topright",inset = c(-0.1,0.4), rev(c(colnames(bar_df1))), 
         fill = rev(c("lightyellow","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) ,
         bty = 1, cex = 0.8)
  par(mar=c(5, 4, 4, 2) + 0.1)
}




# ------------------------------------------------------------------
  
# ui

library(shiny)



ui <- fluidPage(
  titlePanel('Geographical origin prediction of microbiome'),
  sidebarLayout(
    sidebarPanel(
      
    
      selectInput("program","Prediction program",
                  c("Microbial origin prediction based on MetaSub" = "aa",
                    "Build a new prediction model using mGPS" = "cc",
                    "Build a new prediction model using mGPS and predict new samples" = "bb"
                    )
                  ),
      
      conditionalPanel(
        condition = "input.program == 'aa'",
        fileInput(inputId = "f_new_test_1",
                  label = "upload the sample abundance file",
                  accept = c("text/csv",
                             ".csv")),
        actionButton(inputId = "acbutton_1", 
                     label ="Start")
      ),
      
      
      conditionalPanel(
        condition = "input.program == 'bb'",
        fileInput(inputId = "f_new_test_2",
                  label = "upload the testing file (sample abundance file)",
                  accept = c("text/csv",
                             ".csv")),
        radioButtons(inputId = "file_num_2",
                     label = "upload the training file (reference file)",
                     c("Merged metadata and abundance file" = "one_file",
                       "Separate metadata and abundance file" = "two_file")
        ),
        conditionalPanel(
          condition = "input.file_num_2 == 'one_file'",
          fileInput(inputId = "f_new_train_2",
                    label = "upload the reference abundance dataset file",
                    accept = c("text/csv",
                               ".csv")),
        ),
        conditionalPanel(
          condition = "input.file_num_2 == 'two_file'",
          fileInput(inputId = "f_metadata_2",
                    label = "upload the metadata file",
                    accept = c("text/csv",
                               ".csv")),
          fileInput(inputId = "f_abundance_2",
                    label = "upload the abundance file",
                    accept = c("text/csv",
                               ".csv")),
          textInput(inputId = "by_x_2",
                    label = "merge column name in metadata file"),
          textInput(inputId = "by_y_2",
                    label = "merge column name in abundance file")
        ),
        textInput(inputId= "target_2",
                  label ="Please enter the target (same as column name) "),
        textAreaInput(inputId = "hierarchy_2",
                      label = "hierarchy in prediction model (separator:',' )",
                      rows = 1),
        textAreaInput(inputId = "abundance_range_2",
                      label = "column range of abundance data with separator as ':' (eg. 47:100)",
                      rows = 1),
        checkboxInput("control_remove_2", "Remove values (eg. control values)", value = FALSE),
        conditionalPanel( 
          condition = "input.control_remove_2 == true",
          textAreaInput(inputId = "control_value_2",
                        label = "Values need to be removed with format as 'remove_column1:value1,value2;remove_column2:value1,value2'",
                        rows = 1),
        ),
        
        checkboxInput("subsets_2", "Subsets in feature elimination", value = FALSE),
        conditionalPanel( 
          condition = "input.subsets_2 == true",
          textAreaInput(inputId = "subsets_value_2",
                        label = "Subsets with separator as ',' (eg. 50,100,200,300,500,1500)",
                        rows = 1),
        ),
        
        actionButton(inputId = "acbutton_2", 
                     label ="Start")
      ),
      conditionalPanel(
        condition = "input.program == 'cc'",
        radioButtons(inputId = "file_num_3",
                     label = "Input file",
                     c("merged Metadata and abundance file" = "one_file",
                       "separate Metadata and abundance file" = "two_file")
                      ),
        conditionalPanel(
          condition = "input.file_num_3 == 'one_file'",
          fileInput(inputId = "f_new_train_3",
                    label = "upload the reference abundance dataset file",
                    accept = c("text/csv",
                               ".csv")),
        ),
        conditionalPanel(
          condition = "input.file_num_3 == 'two_file'",
          fileInput(inputId = "f_metadata_3",
                    label = "upload the metadata file",
                    accept = c("text/csv",
                               ".csv")),
          fileInput(inputId = "f_abundance_3",
                    label = "upload the abundance file",
                    accept = c("text/csv",
                               ".csv")),
          textInput(inputId = "by_x_3",
                    label = "merge column name in metadata file"),
          textInput(inputId = "by_y_3",
                    label = "merge column name in abundance file")
        ),
        textInput(inputId= "target_3",
                  label ="Please enter the target (same as column name) "),
        textAreaInput(inputId = "hierarchy_3",
                      label = "hierarchy in prediction model (separator:',' )",
                      rows = 1),
        textAreaInput(inputId = "abundance_range_3",
                      label = "column range of abundance data with separator as ':' (eg. 47:100)",
                      rows = 1),
        checkboxInput("control_remove_3", "Remove values (eg. control values)", value = FALSE),
        conditionalPanel( 
          condition = "input.control_remove_3 == true",
          textAreaInput(inputId = "control_value_3",
                        label = "Values need to be removed with format as 'remove_column1:value1,value2;remove_column2:value1,value2'",
                        rows = 1),
        ),
        
        checkboxInput("subsets_3", "Subsets in feature elimination", value = FALSE),
        conditionalPanel( 
          condition = "input.subsets_3 == true",
          textAreaInput(inputId = "subsets_value_3",
                        label = "Subsets with separator as ',' (eg. 50,100,200,300,500,1500)",
                        rows = 1),
        ),
        
        actionButton(inputId = "acbutton_3", 
                     label ="Start")
    )
  ),
    
    mainPanel(
      
      conditionalPanel(
        condition = "input.program == 'aa'",
        tabsetPanel(
          type = "tabs",
          tabPanel("HELP",
                   helpText("This is a web program based on the mGPS application created by Shiny. It can 
                            build a microbial origin prediction model or predict the origin of microbes."),
                   helpText(strong("Microbial origin prediction based on MetaSub")),

                   helpText(strong("Function description: "),"This model can predict the source of microorganisms
                            according to the microorganism database from MetaSub. ","MetaSUB is a mapping project
                            which carried out a large worldwide metagenomic investigation in urban areas and provided
                            a huge microbial sample dataset with metadata.","To learn more about MetaSub, please visit:",
                            a(em("MetaSUB International Consortium inaugural meeting report"),
                              href="https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0168-z"),
                   ),
                   helpText(strong("Usage: "),"The user uploads the microbial abundance data of the sample. (The file can contain multiple samples)"
                   ),
                   helpText(strong("Function and output:"),"The original geographic location predicted by the sample will be displayed on the world
                            map in the Map tab of the program. In the Data tab of the program, users can directly download the predicted source coordinate
                            data of the sample in the uploaded file. (Prelatitude column: predicted source latitude; Prelongitude column: predicted source longitude)"
                   )
          ),
          tabPanel("Map",
                   plotOutput(outputId = "predicted_map_1"),
                   width = "100%"
                   ),
          
          tabPanel("Data",
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
                   helpText("This is a web program based on the mGPS application created by Shiny. It can 
                            build a microbial origin prediction model or predict the origin of microbes."
                   ),
                   helpText(strong("Build a new prediction model using mGPS and predict new samples")),

                   helpText(strong("Function description: "),"This mode can train the microbial source
                            prediction model based on the reference data set uploaded by the user. The
                            prediction model will be used to predict the new sample to be tested provided
                            by the user and report the prediction result of the sample source. (If you need
                            to visualize the accuracy of the model, please use",strong("Build a new prediction model using mGPS"),")"),
                   
                   helpText(strong("Usage: ")),
                   helpText("User need to upload a data file containing microbial abundance. The metadata and abundance data of the sample can be merged into one file (Merged metadata and abundance data), or uploaded as two files (Separate metadata and abundance data)"),
                   helpText("target: the target (marker) of the model selected by the user. It should same as the name of column header. (eg. city)"
                   ),
                   helpText("hierarchy: The marker chain used in mGPS to construct the prediction model (corresponds to the header of the column). In addition to latitude and longitude, there can be at most two. Use ',' as separator  (eg. continent,city,latitude,longitude)"),
                   helpText("column range of abundance data: the number of columns in the file for the column containing abundance data. Use ':' as separator (eg 44:1000)"),
                   helpText("Remove values: user can select to remove some special values in column such as control samples. If checked, user need to input the value and correponding column name. (eg. city:control,neg_control,pos_control;control_type:ctrl cities,negative_control,positive_control)"),
                   helpText(strong("Function and output:")),
                   helpText("The",strong( "Map bar"),"The original geographic location predicted by the sample will be displayed on the world
                            map"
                   ),
                    helpText("The",strong( "Data bar"),"of the program will provide users with the results of using the constructed prediction model to predict the source coordinates of the original data set samples. In addition, users can download the constructed prediction model (In Rda format).")
          ),
          tabPanel("Map",
                   plotOutput(outputId = "predicted_map_2")),
          
          tabPanel("Data",
                   helpText("Predicted original coordinates of samples can be downloaded. Prediction model in RDA format can be downloaded"),
                   downloadButton("downloadData_2",
                                  "Download prediction data"),
                   downloadButton("download_optfeatures_2",
                                  "Download optimal features in prediction model"),
                   downloadButton("download_featuresub_2",
                                  "Download feature subsets accuracy in feature elimination"),
                   downloadButton("downloadModel_2",
                                  "DownloadModel")
                   )
        )
      ),
      
      conditionalPanel(
        condition = "input.program == 'cc'",
        tabsetPanel(
          type = "tabs",
          
          tabPanel("HELP",
                   helpText("This is a web program based on the mGPS application created by Shiny. It can 
                            build a microbial origin prediction model or predict the origin of microbes."
                           ),
                   helpText(strong("Build a new prediction model using mGPS")),

                   helpText(strong("Function description: "),"This model can use the mGPS tool to build a microbial source prediction model based
                            on the microbial abundance data uploaded by the user.","To learn more about mGPS, please visit:",
                            a(em("mGPS"),
                              href="https://github.com/eelhaik/mGPS"),
                   ),
                   helpText(strong("Usage: ")),
                   helpText("User need to upload a data file containing microbial abundance. The metadata and abundance data of the sample can be merged into one file (Merged metadata and abundance data), or uploaded as two files (Separate metadata and abundance data)"),
                  helpText("target: the target (marker) of the model selected by the user. It should same as the name of column header. (eg. city)"
                   ),
                  helpText("hierarchy: The marker chain used in mGPS to construct the prediction model (corresponds to the header of the column). In addition to latitude and longitude, there can be at most two. Use ',' as separator  (eg. continent,city,latitude,longitude)"),
                  helpText("column range of abundance data: the number of columns in the file for the column containing abundance data. Use ':' as separator (eg 44:1000)"),
                  helpText("Remove values: user can select to remove some special values in column such as control samples. If checked, user need to input the value and correponding column name. (eg. city:control,neg_control,pos_control;control_type:ctrl cities,negative_control,positive_control)"),
                  helpText(strong("Function and output:")),
                  helpText("The",strong( "Accuracy bar"),"of the program will show the user the accuracy
                            of the prediction model trained by the mGPS tool and based on the reference microbial database
                            uploaded by the user. The original database will be divided into 5 folds, and mGPS will use 4 of
                            these folds to train the model, and the resulting model will be used to predict the microbial source
                            of the remaining fold. Iteratively obtain the prediction result of the original database and compare
                            it with the actual location of the microorganism."
                  ),
                  helpText("world map: samples' prediction origins are plotted on the world map"),
                  helpText("accuracy bar plot: model built by mGPS accuracy per site for the original reference dataset. mGPS accuracy is shown per-site as the distances between the predicted and true sampling site for the reference samples. The average prediction accuracy across all samples with each population given equal weight is shown on the left. The panel shows the results for QCed site"),
                  helpText("The",strong( "Data bar"),"of the program will provide users with the results of using the constructed prediction model to predict the source coordinates of the original data set samples. In addition, users can download the constructed prediction model (In Rda format).")
                  ),
          
          tabPanel("Accuracy",
                   plotOutput(outputId = "predicted_map_3"),
                   plotOutput(outputId = "predicted_accuracy_3")
                   ),
          tabPanel("Data",
                   helpText("Predicted original coordinates of samples can be downloaded. Prediction model in RDA format can be downloaded"),
                   downloadButton("downloadData_3",
                                  "Download prediction data"),
                   downloadButton("download_optfeatures_3",
                                  "Download optimal features in prediction model"),
                   downloadButton("download_featuresub_3",
                                  "Download feature subsets accuracy in feature elimination"),
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
  
  
  
  data_input <- reactive({
    req(input$program)
    data_input <- input$program
    
  })

  coastlines <- reactive({
    coastlines <- coastlines_f()
    return(coastlines)
  })
  

    
      model_store <- reactive({
          load("model_Metasub/20211011/MetaSub_model_2.Rda")
          return(model_store)
      })
      
      test_f <- reactive({
        if(input$program == "aa" ){
          req(input$f_new_test_1)
          test_f <- read.csv(input$f_new_test_1$datapath, 
                             header = TRUE)
        }
        if(input$program == "bb" ){
          req(input$f_new_test_2)
          test_f <- read.csv(input$f_new_test_2$datapath, 
                             header = TRUE)}
          return(test_f)
      })
      
      MetasubDataPreds <- reactive({
          MetasubDataPreds <- MetaSub_prediction(model_store(),test_f(),coastlines())
        return(MetasubDataPreds)    
      })

    observeEvent(input$acbutton_1,{
    
    output$predicted_map_1 <- renderPlot({
      
      plot_prediction(MetasubDataPreds())
    })
    
    output$downloadData_1 <- downloadHandler(
      filename = function(){
        return("Predicted_origin.csv")
      },
      content = function(file){
        write.csv(MetasubDataPreds(),file)
      }
    )
    
    })
    
    #-------------------------------------------------------

    observeEvent(input$acbutton_2,{
      output$predicted_map_2 <- renderPlot({
        df <- prediction_output()[[1]]
        plot_prediction(df)
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
          return("optimal_features.csv")
        },
        content = function(file){
          write.csv(v(),file)
        })
      
      output$download_featuresub_2 <- downloadHandler(
        
        filename = function(){
          return("features_subsets_accuracy.csv")
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

  # -------------------------------------------------------------------------
  # 
  # 
    
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
    
    control_text <- reactive({
      if(input$program == "bb" ){
        req(input$control_value_2)
        control_text <- input$control_value_2
      }
      if(input$program == "cc" ){
        req(input$control_value_3)
        control_text <- input$control_value_3
      }
      return(control_text)
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
          subsets_in <- as.numeric(strsplit(text,",")[[1]])
         }else{
          subsets_in <- NULL
          }}
      if(input$program == "cc" ){
        if (input$subsets_2 == T){
          req(input$subsets_value_3)
          text <- input$subsets_value_3
          subsets_in <- as.numeric(strsplit(text,",")[[1]])
        }else{
          subsets_in <- NULL
        }}
      
      return(subsets_in)
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
            if (input$control_remove_2 == T){
              train_f <- remove_control(train_f,control_text()) 
            }
            return(train_f)
          }
          
          if (input$file_num_2 == "two_file"){
            
            req(input$f_metadata_2)
            req(input$f_abundance_2)
            data_metadata <- read.csv(input$f_metadata_2$datapath, 
                                      header = TRUE)
            data_abundance <- read.csv(input$f_abundance_2$datapath,
                                       header = T)
            data_abundance <- unique(data_abundance)
            by_x <- by_x_in()
            by_y <- by_y_in()
            print(by_x)
            print(by_y)
            metasub_data <- merge(data_metadata,data_abundance,by.x=by_x,by.y=by_y) # #merge bacterial and meta data

            if (input$control_remove_2 == T){
              metasub_data <- remove_control(metasub_data,control_text()) 
            }
            
            train_f <- metasub_data
            return(train_f)
          }
        }
        
        if(input$program == "cc" ){
        
        if(input$file_num_3 == "one_file" ){
          
          req(input$f_new_train_3)
          train_f <- read.csv(input$f_new_train_3$datapath, 
                              header = TRUE)
          if (input$control_remove_3 == T){
            train_f <- remove_control(train_f,control_text()) 
          }
          return(train_f)
        }
        
        if (input$file_num_3 == "two_file"){
          
          req(input$f_metadata_3)
          req(input$f_abundance_3)
          data_metadata <- read.csv(input$f_metadata_3$datapath, 
                                    header = TRUE)
          data_abundance <- read.csv(input$f_abundance_3$datapath,
                                     header = T)
          data_abundance <- unique(data_abundance)
          metasub_data <- merge(data_metadata,data_abundance,by.x=by_x_in(),by.y=by_y_in()) # #merge bacterial and meta data
          
          if (input$control_remove_3 == T){
            metasub_data <- remove_control(metasub_data,control_text()) 
          }
          
          train_f <- metasub_data
          return(train_f)
        }
        }
      })
      
      
     
      metasub_data <- reactive({
          metasub_data <- data_preprocess_f(train_f(),classTarget_in(),hierarchy_in())
          return(metasub_data)

      })
      
      featureElim <- reactive({
        range_1 <- abundance_r()[1]
        range_2 <- abundance_r()[2]
        featureElim <- featureElimination_f(metasub_data(),classTarget_in(),range_1,range_2, subsets_in()) # 44:3712 for MetaSub
        return(featureElim)
      })
      
      prediction_output <- reactive({
        if(input$program == "bb" ){
          test_file <- test_f()
          write.csv(test_file,"Outputs/prediction_output_test_f.csv")
          prediction_output <- New_model_prediction(train_f(),test_f(),classTarget_in(),featureElim()$optVariables,hierarchy_in(),coastlines())
        }
        if(input$program == "cc" ){
          prediction_output <- model_accuracy_f(metasub_data(),featureElim()$optVariables,classTarget_in(),hierarchy_in(),coastlines())
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
          df <- prediction_output()[[1]]
          plot_map(df)
        })
        
        output$predicted_accuracy_3 <- renderPlot({
          plot_accuracy(prediction_output()[[1]],"city")
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
            return("optimal_features.csv")
          },
          content = function(file){
            write.csv(v(),file)
          })
        
        output$download_featuresub_3 <- downloadHandler(
          
          filename = function(){
            return("features_subsets_accuracy.csv")
          },
          content = function(file){
            write.csv(git_subset(),file)
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






