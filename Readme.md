# mGPS Interface

## Introduction

This is a web program based on the mGPS application created by Shiny. It can build a microbial origin prediction model or predict the origin of microbes. A detailed explanation of the output result is reported in the tutorial file: mGPS Interface Tutorial.pdf 
To learn more about mGPS, please visit: [mGPS](https://github.com/eelhaik/mGPS)  

## Usage

Open the R script **mGPS_interface.r** in R and click label ```Run App``` in R . Then the user can select the function mode and operate on the interface. Do not jump to other function modes when using one function mode (that is, don't click on other “Prediction programs” during the running process). It needs to reload the app after using one function mode and then can use another function mode.
The output files that the user can download from the interface will be automatically stored in the **Outputs** folder. Detailed usage is introduced in each function part in **Prediction programs** part.

## Dependencies

Required packages can be found in ```packages.r```

## Prediction programs

### 1 Build a new prediction model using mGPS

- Function description  
  In this mode, users can use the mGPS tool to build a microbial source prediction model based on the microbial abundance data uploaded by the user. To learn more about mGPS, please visit: [mGPS](https://github.com/eelhaik/mGPS)  

- Usage    
  
  - A. **In left sidebar**:
  
  - Select ```Prediction program``` as `Build a new prediction model using mGPS`  
  
  - ```Input file(s)```: Upload data file(s) (in .csv format) containing microbial abundance data and metadata.   
    In metadata, at least one locality (eg. continent, city) and coordinates (necessary) data columns should be included. The metadata and abundance data of the sample can be merged into one file ( `Merged metadata and abundance data` ) or uploaded as two files ( `Separate metadata and abundance data` )   
    When the `Separate metadata and abundance file` is selected. `Merge column name in metadata/abundance file`: Input the header name of the column which is the merged column in two files.    
  
  - ```Enter the main locality level``` -  Input the main locality target. It should same as that column header. (eg. city)  
  
  - ```Enter the locality hierarchy``` -  The locality chain used in mGPS to construct the prediction model (same column headers). It should contain one or two locality information, latitude and longitude. Use ',' as the separator. (eg. continent,city,latitude,longitude)
  
  - ```Column range of abundance data``` -  the number of columns in the file for the column containing abundance data. Use ':' as separator (eg 44:1000)    
  
  - ```Locality sample size cut off``` -  *(Optional)* Remove locality whose sample size is less than a certain value. If checked, input the cut off number (eg. 8)  
  
  - ```Subsets in feature elimination```: *(Optional)* Limit the number of features to a certain value. If unchecked, mGPS will find the optimal subset size of microbiome features. If checked, there are three types of input format: a. Input sizes of the subset with separator as ',' (eg. 50,100,200,300); b. Input the size range of subset with separator as '-' (eg. 50-300); c. Input a single value.  
       Hint: In addition to checking the given subset sizes, the algorithm also tries to use all features(taxa) to predict the origin. If the accuracy of the subset sizes uploaded by the user is lower than the accuracy of using all features, the algorithm will still choose to use all features for prediction to improve the accuracy. If the user still wants to use a specific number of features, the feature columns in the input file can be manually filtered according to the order of importance of the features in the output file "*Optimal_features.csv*" (Download optimal features in prediction model).
  
  - **B. In `Result Plot` tab**:    
  
  - ```Change longitude/latitude range in output map``` *(Optional)*  
  
  - ```Whether pull points to land/waterbody```: *(Optional)* If checked, predicted origin location will be pulled to the nearest land/waterbody if predicted coordinates are out of the expected border.   
  
  - **C. Start the program:** Click the ```Start``` bar and then click the ```Result Plot``` tab   
  
  - **D. Data processing:** Please wait while output files are being generated. When the prompt bar disappears you can see the results and download files.  


- Result plot and output  
  For more explanation of output results, view the tutorial file "*mGPS Interface Tutorial.pdf*"  
  
  - ```Result Plot``` tab:  
    Show the accuracy of the prediction model trained by the mGPS tool and based on the reference microbial database uploaded by the user.  
    The original database will be divided into 5 folds, and mGPS will use 4 of these folds to train the model, and the resulting model will be used to predict the microbial source of the remaining fold. Iteratively obtain the prediction result of the original database and compare it with the actual location of the microorganism.    
    
    - *World map*: Samples' prediction origins are plotted on the world map  
    - *Accuracy bar plot*: model built by mGPS accuracy per site for the original reference dataset. mGPS accuracy is shown per-site as the distances between the predicted and true sampling sites for the reference samples. The average prediction accuracy across all samples with each population given equal weight is shown on the left.   
  
  - ```Output``` tab:  
    The results of using the constructed prediction model to predict the source coordinates of the original dataset samples. In addition, the constructed prediction model can be downloaded (In .Rda format).  
    - ```Download prediction data``` -  The prediction results of using the constructed prediction model to predict the source coordinates of the original dataset samples  
    - ```Download optimal features in prediction model``` -  The optimal features used in prediction model construction.  
    - ```Download feature subsets accuracy in feature elimination``` -  The prediction accuracy for different sizes of feature subsets.  
    - ```DownloadModel``` -  The constructed prediction model (In Rda format). Users can load this model into R to check the detailed model information through  

### 2 Build a new prediction model using mGPS and predict new samples

- Function description  
  In this mode, users can train the microbial source prediction model based on the reference data set uploaded by the user. The prediction model will be used to predict the new sample to be tested provided by the user and report the prediction result of the sample source. (If you need to visualize the accuracy of the model, please use *Build a new prediction model using mGPS* model function )  

- Usage  
  
  - A. **In left sidebar**:  
  
  - Select ```Prediction program``` as `Build a new prediction model using mGPS and predict new samples`   
  
  - ```Upload new sample(s) abundance file```: Upload file (in .csv format) containing abundance data of new sample(s).   
  
  - ```Upload reference file(s)```: Upload data file(s) (in .csv format) containing microbial abundance data and metadata.   
    In metadata, at least one locality (eg. continent, city) and coordinates (necessary) data columns should be included. The metadata and abundance data of the sample can be merged into one file ( `Merged metadata and abundance data` ) or uploaded as two files ( `Separate metadata and abundance data` )   
    When the `Separate metadata and abundance file` is selected. ```Merge column name in metadata/abundance file```: Input the header name of the column which is the merged column in two files.   
  
  - ``` Enter the main locality level``` -  Input the main locality target. It should same as that column header. (eg. city)     
    ```Enter the locality hierarchy``` -  The locality chain used in mGPS to construct the prediction model (same column headers). It should contain one or two locality information, latitude and longitude. Use ',' as the separator. (eg. continent,city,latitude,longitude)   
  
  - ```Column range of abundance data``` -  the number of columns in the file for the column containing abundance data. Use ':' as separator (eg 44:1000)  
  
  - ```Locality sample size cut off``` -  *(Optional)* Remove locality whose sample size is less than a certain value. If checked, input the cut off number (eg. 8)   
  
 - ```Subsets in feature elimination```: *(Optional)* Limit the number of features to a certain value. If unchecked, mGPS will find the optimal subset size of microbiome features. If checked, there are three types of input format: a. Input sizes of the subset with separator as ',' (eg. 50,100,200,300); b. Input the size range of subset with separator as '-' (eg. 50-300); c. input a single value.  
       Hint: In addition to checking the given subset sizes, the algorithm also tries to use all features(taxa) to predict the origin. If the accuracy of the subset sizes uploaded by the user is lower than the accuracy of using all features, the algorithm will still choose to use all features for prediction to improve the accuracy. If the user still wants to use a specific number of features, the feature columns in the input file can be manually filtered according to the order of importance of the features in the output file "*Optimal_features.csv*" (Download optimal features in prediction model).
  
  - **B. In `Result Plot` tab**:    
  
  - ```Change longitude/latitude range in output map``` *(Optional)*  
  
  - ```Whether pull points to land/waterbody```: *(Optional)* If checked, predicted origin location will be pulled to the nearest land/waterbody if predicted coordinates are out of the expected border.   
  
  - **C. Start the program:** Click the ```Start``` bar and then click the ```Result Plot``` tab   
  - **D. Data processing:** Please wait while output files are being generated. When the prompt bar disappears you can see the results and download files.  

- Result plot and output 
  For more explanation of output results, view the tutorial file "*mGPS Interface Tutorial.pdf*"  
  
  - ```Result Plot``` tab:  
    The reference datasets will be used to construct an origin prediction model by mGPS. Then this model will be used to predict the origin of new samples.  
    - *World map*:  samples' prediction origins are plotted on the world map  
    
  - ```Output``` tab:   
    The results of using the constructed prediction model to predict the source coordinates of the new samples. In addition, the constructed prediction model can be downloaded (In .Rda format).  
    - ```Download prediction data``` -  The origin prediction result of sample in csv table.  
    - ```Download optimal features in prediction model``` -  The optimal features used in prediction model construction.  
    - ```Download feature subsets accuracy in feature elimination``` -  The prediction accuracy for different sizes of feature subsets.  
    - ```DownloadModel``` -  The constructed prediction model (In Rda format). Users can load this model into R to check the detailed model information through `load("Prediction_model.Rda")`

### 3 Use an existing model to predict new samples

- Function description  
  In this mode, users can predict new sample origin based on an existing prediction model.   

- Usage  
  
  - A. **In left sidebar**:  
  
  - Select ```Prediction program``` as `Use an existing model to predict new samples`   
  
  - ```Upload new sample(s) abundance file```: Upload data file(s) (in .csv format) containing new microbial sample abundance data.   
  
  - ```Upload the prediction model```: Upload a constructed origin prediction model in .Rda format. The model can be downloaded in the `Output` tab of function: *Build a new prediction model using mGPS* or *Build a new prediction model using mGPS and predict new samples*   
  
  - **B. In `Result Plot` tab**:    
  
  - ```Change longitude/latitude range in output map``` *(Optional)*  
  
  - ```Whether pull points to land/waterbody```: *(Optional)* If checked, predicted origin location will be pulled to the nearest land/marine if predicted coordinates are out of the expected border.   
  
  - **C. Start the program:** Click the ```Start``` bar and then click the ```Result Plot``` tab   
  - **D. Data processing:** Please wait while output files are being generated. When the prompt bar disappears you can see the results and download files.  

- Result plot and output  
  For more explanation of output results, view the tutorial file "*mGPS Interface Tutorial.pdf*"     

  - ```Result Plot``` tab:  
    The existing model will be used to predict the origin of new sample(s).  
    - *World map*:  new samples' prediction origins are plotted on the world map     
    
  - ```Output``` tab:   
    The results of using the existing prediction model to predict the source coordinates of new sample(s). Predicted source coordinate data of the sample in the uploaded file. (*Prelatitude column*: predicted original latitude; *Prelongitude column*: predicted original longitude)      

## Example

Users can find the test files (example files) in the folder **Example_file**.  

### 1 Build a new prediction model using mGPS  

For short runtimes, the sample database contains only a small amount of sample data, so models built from this reference dataset are less accurate. It runs for a different time on different computers, between about 10-40min.  

- *Merged_training_dataset.csv*: Merged metadata and abundance data file.  

- *Separate_training_abundance.csv*: Abundance data of samples.  

- *Separate_training_metadata.csv*: Metadata of samples.  

1. ```Input file```
   - **Merged metadata and abundance file**: upload file *Merged_training_dataset.csv*
   - **Separate metadata and abundance file**: upload file *Separate_training_metadata.csv* in ```Upload the metadata file``` and upload file *Separate_training_abundance.csv* in ```Upload the abundance file```   
     ```merge column name in metadata file``` - uuid  
     ```merge column name in abundance file``` - uuid  

2. ```Enter the main locality level``` - city  
   ```Enter the locality hierarchy``` - continent,city,latitude,longitude  
   ```Column range of abundance data``` - 43:70  

3. If select ```Locality sample size cut off (Optional)```:  
   ```Cut off of sample number ```  - 3     

4. If select ```Subsets in feature elimination (Optional)``` :
   ```Subsets size``` - 5,15,20,25

5. Start the program: Click the ```Start``` bar and then click the ```Result Plot``` tab   

6. In the ```Result Plot``` tab, the latitude and longitude range can be adjusted. The prediction points can be pulled to land and waterbody.  

7. For more explanation of output results, view the tutorial file "*mGPS Interface Tutorial.pdf*" 

### 2 Build a new prediction model using mGPS and predict new samples

For short runtimes, the sample database contains only a small amount of sample data, so models built from this reference dataset are less accurate. It runs for a different time on different computers, between about 10-40min.  

- *Sample_prediction_MetaSub.csv*: New samples with taxa abundance data.  

1. ```Upload new sample abundance file```-  upload file *Sample_prediction_MetaSub.csv*  

2. ```Upload reference file(s)```  
   - **Merged metadata and abundance file**: upload file *Merged_training_dataset.csv*
   - **Separate metadata and abundance file**: upload file *Separate_training_metadata.csv* in ```Upload the metadata file``` and upload file *Separate_training_abundance.csv* in ```Upload the abundance file```   
     ```merge column name in metadata file``` - uuid  
     ```merge column name in abundance file``` - uuid  

3. Similar enter information and operations as function 2 **Build a new prediction model using mGPS**.  
4. For more explanation of output results, view the tutorial file "*mGPS Interface Tutorial.pdf*"   

### 3 Use an existing model to predict new samples

- *Sample_prediction_MetaSub.csv*: samples that can be used to do the prediction  
- *Data/model_Metasub/MetaSub_model.Rda*: A model built from the MetaSub dataset. It can be used to predict the sample origin in urban. To learn more about MetaSub, please visit: [MetaSUB International Consortium inaugural meeting report](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0168-z)   
1. ```Upload sample(s) abundance file```: *Sample_prediction_MetaSub.csv*  
2. ```Upload the prediction model (In .Rda format)```: *MetaSub_model.Rda* 
3. Start the program: Click the ```Start``` bar and then click the ```Result Plot``` tab   
4. In the ```Result Plot``` tab, the latitude and longitude range can be adjusted. The prediction points can be pulled to land and waterbody.  
5. For more explanation of output results, view the tutorial file "*mGPS Interface Tutorial.pdf*"   
