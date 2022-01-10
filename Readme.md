# mGPS Interface

## Introduction

This is a web program based on the mGPS application created by Shiny. It can build a microbial origin prediction model or predict the origin of microbes. The detail explanation of output result is reported in the tutorial file: mGPS Interface Tutorial.pdf 
To learn more about mGPS, please visit: [mGPS](https://github.com/eelhaik/mGPS)  

### Usage

Open the R script **mGPS_interface.r** in R and click label ```Run App``` in R . Then user can select the function model and operate on interface. The output files that user can download from interface will be automaticlly stored in **Outputs** folder.

## Model Function

### 1 Build a new prediction model using mGPS

- Function description  
  This model can use the mGPS tool to build a microbial source prediction model based on the microbial abundance data uploaded by the user. To learn more about mGPS, please visit: [mGPS](https://github.com/eelhaik/mGPS)  

- Usage    
  
  - A. **In left side bar**:
  
  - Select ```Prediction program``` as `Build a new prediction model using mGPS`  
  
  - ```Input file(s)```: Upload data file(s) (in .csv format) containing microbial abundance data and metadata.   
    In metadata, at least one locality (eg. continent, city) and coordinates (necessary) data columns should be included. The metadata and abundance data of the sample can be merged into one file ( `Merged metadata and abundance data` ), or uploaded as two files ( `Separate metadata and abundance data` )   
    When `Separate metadata and abundance file` is selected. `Merge column name in metadata/abundance file`: Input the header name of column which is the merged column in two files.    
  
  - ```Enter the main locality level``` -  Input the main locality target. It should same as that column header. (eg. city)  
  
  - ```Enter the locality hierarchy``` -  The locality chain used in mGPS to construct the prediction model (same column headers). It should contain one or two locality information, latitude and longitude. Use ',' as separator. (eg. continent,city,latitude,longitude)
  
  - ```Column range of abundance data``` -  the number of columns in the file for the column containing abundance data. Use ':' as separator (eg 44:1000)    
  
  - ```Locality sample size cut off``` -  *(Optional)* Remove locality whose sample size is less than a certain value. If checked, input the cut off number (eg. 8)  
  
  - ``` Remove values``` -  *(Optional)* Remove some special values in column such as control samples. If checked, input column name(s) and corresponding value(s). (format eg. city:control,neg_control,pos_control;control_type:ctrl cities,negative_control,positive_control)   
  
  - ```Subsets in feature elimination```: *(Optional)* Set the size of subsets used in feature elimination (find the optimal subset size of microbiome features). If checked, input the subsets size with separator as ',' (eg. 50,100,200,300)  
  
  - **B. In `Result Plot` tab**:    
  
  - ```Change longitude/latitude range in output map``` *(Optional)*  
  
  - ```Whether pull points to land/marine```: *(Optional)* If checked, predicted origin location will be pull to the nearest land/marine if predicted coordinates are out of the expected boarder.   
  
  - **C. Start the program:** Click the ```Start``` bar and then click the ```Result Plot``` tab   

- Result plot and output  
  
  - ```Result Plot``` tab:  
    Show the accuracy of the prediction model trained by the mGPS tool and based on the reference microbial database uploaded by the user.  
    The original database will be divided into 5 folds, and mGPS will use 4 of these folds to train the model, and the resulting model will be used to predict the microbial source of the remaining fold. Iteratively obtain the prediction result of the original database and compare it with the actual location of the microorganism.    
    
    - *World map*: Samples' prediction origins are plotted on the world map  
    - *Accuracy bar plot*: model built by mGPS accuracy per site for the original reference dataset. mGPS accuracy is shown per-site as the distances between the predicted and true sampling site for the reference samples. The average prediction accuracy across all samples with each population given equal weight is shown on the left.   
  
  - ```Output``` tab:  
    The results of using the constructed prediction model to predict the source coordinates of the original dataset samples. In addition, the constructed prediction model can be downloaded (In Rda format).  
    
    - ```Download prediction data``` -  The prediction results of using the constructed prediction model to predict the source coordinates of the original dataset samples  
    - ```Download optimal features in prediction model``` -  The optimal features used in prediction model construction.  
    - ```Download feature subsets accuracy in feature elimination``` -  The prediction accuracy for different size of feature subsets.  
    - ```DownloadModel``` -  The constructed prediction model (In Rda format). User can load this model into R to check the detailed model information through  

### 2 Build a new prediction model using mGPS and predict new samples

- Function description  
  This mode can train the microbial source prediction model based on the reference data set uploaded by the user. The prediction model will be used to predict the new sample to be tested provided by the user and report the prediction result of the sample source. (If you need to visualize the accuracy of the model, please use *Build a new prediction model using mGPS* model function )  

- Usage  
  
  - A. **In left side bar**:  
  
  - Select ```Prediction program``` as `Build a new prediction model using mGPS and predict new samples`   
  
  - ```Upload new sample(s) abundance file```: Upload file (in .csv format) containing abundance data of new sample(s).   
  
  - ```Upload reference file(s)```: Upload data file(s) (in .csv format) containing microbial abundance data and metadata.   
    In metadata, at least one locality (eg. continent, city) and coordinates (necessary) data columns should be included. The metadata and abundance data of the sample can be merged into one file ( `Merged metadata and abundance data` ), or uploaded as two files ( `Separate metadata and abundance data` )   
    When `Separate metadata and abundance file` is selected. ```Merge column name in metadata/abundance file```: Input the header name of column which is the merged column in two files.   
  
  - ``` Enter the main locality level``` -  Input the main locality target. It should same as that column header. (eg. city)     
    ```Enter the locality hierarchy``` -  The locality chain used in mGPS to construct the prediction model (same column headers). It should contain one or two locality information, latitude and longitude. Use ',' as separator. (eg. continent,city,latitude,longitude)   
  
  - ```Column range of abundance data``` -  the number of columns in the file for the column containing abundance data. Use ':' as separator (eg 44:1000)  
  
  - ```Locality sample size cut off``` -  *(Optional)* Remove locality whose sample size is less than a certain value. If checked, input the cut off number (eg. 8)   
  
  - ```Remove values``` -  *(Optional)* Remove some special values in column such as control samples. If checked, input column name(s) and corresponding value(s). (format eg. city:control,neg_control,pos_control;control_type:ctrl cities,negative_control,positive_control)    
  
  - ```Subsets in feature elimination```: *(Optional)* Set the size of subsets used in feature elimination (find the optimal subset size of microbiome features). If checked, input the subsets size with separator as ',' (eg. 50,100,200,300)   
  
  - **B. In `Result Plot` tab**:    
  
  - ```Change longitude/latitude range in output map``` *(Optional)*  
  
  - ```Whether pull points to land/marine```: *(Optional)* If checked, predicted origin location will be pull to the nearest land/marine if predicted coordinates are out of the expected boarder.   
  
  - **C. Start the program:** Click the ```Start``` bar and then click the ```Result Plot``` tab   

- Result plot and output  
  
  - ```Result Plot``` tab:  
    The reference datasets will be used to construct a origin prediction model by mGPS. Then this model will be used to predict origin of new samples.  
    
    - *World map*:  samples' prediction origins are plotted on the world map  
  
  - ```Output``` tab:   
    The results of using the constructed prediction model to predict the source coordinates of the new samples. In addition, the constructed prediction model can be downloaded (In Rda format).  
    
    - ```Download prediction data``` -  The origin prediction result of sample in csv table.  
    - ```Download optimal features in prediction model``` -  The optimal features used in prediction model construction.  
    - ```Download feature subsets accuracy in feature elimination``` -  The prediction accuracy for different size of feature subsets.  
    - ```DownloadModel``` -  The constructed prediction model (In Rda format). User can load this model into R to check the detailed model information through `load("Prediction_model.Rda")`

### 3 Use existing model to predict new samples

- Function description  
  This mode can predict new sample origin based on an exsiting prediction model.   

- Usage  
  
  - A. **In left side bar**:  
  
  - Select ```Prediction program``` as `Use existing model to predict new samples`   
  
  - ```Upload new sample(s) abundance file```: Upload data file(s) (in .csv format) containing new microbial sample abundance data.   
  
  - ```Upload the prediction model```: Upload a constructed origin prediction model in .Rda format. Model can be downloaded in `Output` tab of function: *Build a new prediction model using mGPS* or *Build a new prediction model using mGPS and predict new samples*   
  
  - **B. In `Result Plot` tab**:    
  
  - ```Change longitude/latitude range in output map``` *(Optional)*  
  
  - ```Whether pull points to land/marine```: *(Optional)* If checked, predicted origin location will be pull to the nearest land/marine if predicted coordinates are out of the expected boarder.   
  
  - **C. Start the program:** Click the ```Start``` bar and then click the ```Result Plot``` tab   

- Result plot and output  
  
  - ```Result Plot``` tab:  
    The exsiting model will be used to predict the origin of new sample.  
    
    - *World map*:  new samples' prediction origins are plotted on the world map     
  
  - ```Output``` tab:   
    The results of using the exsiting prediction model to predict the source coordinates of new sample(s). Predicted source coordinate data of the sample in the uploaded file. (*Prelatitude column*: predicted original latitude; *Prelongitude column*: predicted original longitude)      

## Example

User can find the test files (example files) in folder **Example_file**.  

### 1 Build a new prediction model using mGPS

- *Merged_training_dataset.csv*: Merged metadata and abundance data file.  

- *Separate_training_abundance.csv*: Abundance data of samples.  

- *Separate_training_metadata.csv*: Metadata of samples.  
1. ```Input file```
   
   - **Merged metadata and abundance file**: upload file *Merged_training_dataset.csv*
   - **Separate metadata and abundance file**:upload file *Separate_training_metadata.csv* in ```Upload the metadata file``` and upload file *Separate_training_abundance.csv* in ```Upload the abundance file```   
     ```merge column name in metadata file``` - uuid  
     ```merge column name in abundance file``` - uuid  

2. ```Enter the main locality level``` - city  
   ```Enter the locality hierarchy``` - continent,city,latitude,longitude  
   ```Column range of abundance data``` - 43:70  

3. If select ```Locality sample size cut off (Optional)```:  
   ```Cut off of sample number ```  - 3  

4. If select ```Cut off of sample number```:  
   ```Values need to be removed ``` - city:control,other_control,neg_control,other,pos_control;control_type:ctrl cities,negative_control,positive_control    

5. If select ```Subsets in feature elimination (Optional)``` :
   ```Subsets size``` - 5,15,20,25

6. Start the program: Click the ```Start``` bar and then click the ```Result Plot``` tab   

7. In ```Result Plot``` tab, the latitude and logitude range can be adjusted. The prediction points can be pull to land and marine.

### 2 Build a new prediction model using mGPS and predict new samples

Similar enter information and oprations as function 2 **Build a new prediction model using mGPS**.  

### 3 Use existing model to predict new samples

*Sample_prediction_MetaSub.csv* : Select ```Use existing model to predict new samples function``` in prediction program (left side bar) and upload this file on interface. Then click ```start``` and click the ```Result Plot``` tab. In ```Result Plot``` tab, the latitude and logitude range can be adjusted. The prediction points can be pull to land and marine.
