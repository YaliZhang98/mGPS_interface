# mGPS Interface



## Introduction

This is a web program based on the mGPS application created by Shiny. It can build a microbial origin prediction model or predict the origin of microbes.  
To learn more about mGPS, please visit: [mGPS](https://github.com/eelhaik/mGPS)  

### Usage
Open the R script **mGPS_interface.r** in R and click label ```Run App``` in R . Then user can select the function model and operate on interface. The output files that user can download from interface will be automaticlly stored in **Outputs** folder.



## Model Function



### 1 Microbial origin prediction based on MetaSub  
- Function description  
  This model can predict the source of microorganisms according to the microorganism database from MetaSub. MetaSUB is a mapping project which carried out a large worldwide metagenomic investigation in urban areas and provided a huge microbial sample dataset with metadata. To learn more about MetaSub, please visit: [MetaSUB International Consortium inaugural meeting report](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0168-z)  
- Usage  
  The user uploads the microbial abundance data of the sample (The file can contain multiple samples). Then click start and click the **Map** tab.  
- Output  
  **Map** bar: The original geographic location predicted by the sample.  
  **Data** bar: Predicted source coordinate data of the sample in the uploaded file. (*Prelatitude column*: predicted original latitude; *Prelongitude column*: predicted original longitude)  



### 2 Build a new prediction model using mGPS

- Function description  
  This model can use the mGPS tool to build a microbial source prediction model based on the microbial abundance data uploaded by the user. To learn more about mGPS, please visit: [mGPS](https://github.com/eelhaik/mGPS)  

- Usage  
  User need to upload a reference training file with microbial abundance and corresponding coordinates. The metadata and abundance data of samples can be merged into one file (Merged metadata and abundance data), or uploaded as two files (Separate metadata and abundance data).  

  - If user upload separate metadata and abundance data, the merge column name should be entered in `merge column name in metadata file` and `merge column name in abundance file` textarea.  

  - ```target``` -  the target (marker) of the model selected by the user. It should same as the name of column header. (eg. city)  
    ```hierarchy``` -  The marker chain used in mGPS to construct the prediction model (corresponds to the header of the column). In addition to latitude and longitude, there can be at most two. Use ',' as separator (eg. continent,city,latitude,longitude). Hierarchy should same as the name of column header.  
  - ```column range of abundance data``` -  the number of columns in the file for the column containing abundance data. Use ':' as separator (eg 44:1000)  
  - ```Remove values``` -  user can select to remove some special values in column such as control samples. If checked, user need to input the value and correponding column name with format as "colomn_name_1:value_1,value_2;colomn_name_2:value_1,value_2" (eg. city:control,neg_control,pos_control;control_type:ctrl cities,negative_control,positive_control)  
  - ```Subsets in feature elimination``` -  user can select set the size of subsets used in feature elimination. If checked, user need to input the subsets size with separator as ',' (eg. 50,100,200,300,500,1500)  
  - After enter all the values, please click start and click the **Map** tab.  

- Output  

  - **Accuracy** bar:  
    The accuracy of the prediction model trained by the mGPS tool based on the reference microbial database uploaded by the user.   
    The original database will be divided into 5 folds, and mGPS will use 4 of these folds to train the model. The built model will be used to predict the microbial source of the remaining fold. Iteratively obtain the prediction result of the original database and compare it with the actual location of the microorganism.  
    - *World map*: Samples' prediction origins are plotted on the world map  
    - *Accuracy bar plot*: model built by mGPS accuracy per site for the original reference dataset. mGPS accuracy is shown per-site as the distances between the predicted and true sampling site for the reference samples. The average prediction accuracy across all samples with each population given equal weight is shown on the left.   
  - **Data** bar:  
    - ```Download prediction data``` -  The prediction results of using the constructed prediction model to predict the source coordinates of the original dataset samples  
    - ```Download optimal features in prediction model``` -  The optimal features used in prediction model construction.  
    - ```Download feature subsets accuracy in feature elimination``` -  The prediction accuracy for different size of feature subsets.  
    - ```DownloadModel``` -  The constructed prediction model (In Rda format). User can load this model into R to check the detailed model information through  



### 3 Build a new prediction model using mGPS and predict new samples  

- Function description  
  This mode can train the microbial source prediction model based on the reference data set uploaded by the user. The prediction model will be used to predict the new sample to be tested provided by the user and report the prediction result of the sample source. (If you need to visualize the accuracy of the model, please use **Build a new prediction model using mGPS** model function )  

- Usage  
  User need to upload a testing file with samples and a reference training file with microbial abundance and corresponding coordinates. For training file, the metadata and abundance data of the sample can be merged into one file (Merged metadata and abundance data), or uploaded as two files (Separate metadata and abundance data).  
  - If user upload separate metadata and abundance data, the merge column name should be entered in ```merge column name in metadata file``` and ```merge column name in abundance file``` textarea.  
  - ```target``` -  the target (marker) of the model selected by the user. It should same as the name of column header. (eg. city)  
    ```hierarchy``` -  The marker chain used in mGPS to construct the prediction model (corresponds to the header of the column). In addition to latitude and longitude, there can be at most two. Use ',' as separator (eg. continent,city,latitude,longitude). Hierarchy should same as the name of column header.    
  - ```column range of abundance data``` -  the number of columns in the file for the column containing abundance data. Use ':' as separator (eg 44:1000)  
  - ```Remove values``` -  user can select to remove some special values in column such as control samples. If checked, user need to input the value and correponding column name with format as "colomn_name_1:value_1,value_2;colomn_name_2:value_1,value_2" (eg. city:control,neg_control,pos_control;control_type:ctrl cities,negative_control,positive_control)  
  - ```Subsets in feature elimination``` -  user can select set the size of subsets used in feature elimination. If checked, user need to input the subsets size with separator as ',' (eg. 50,100,200,300,500,1500)  
  - After enter all the values, please click start and click the **Map** tab.  

- Output  
  - **Map** bar: The original geographic location predicted by the sample will be displayed on the world map  
  - **Data** bar:   
    - ```Download prediction data``` -  The origin prediction result of sample in csv table.  
    - ```Download optimal features in prediction model``` -  The optimal features used in prediction model construction.  
    - ```Download feature subsets accuracy in feature elimination``` -  The prediction accuracy for different size of feature subsets.  
    - ```DownloadModel``` -  The constructed prediction model (In Rda format). User can load this model into R to check the detailed model information through `load("Prediction_model.Rda")`




## Example  

User can find the test files (example files) in folder **Example_file**.  



### 1 Microbial origin prediction based on MetaSub  

*Sample_prediction_MetaSub.csv* : Select **Microbial origin prediction based on MetaSub** function in prediction program and upload this file on interface. Then click start and click the **Map** tab.  



### 2 Build a new prediction model using mGPS  

- *Merged_training_dataset.csv*: Merged metadata and abundance data file.  
- *Separate_training_abundance.csv*: Abundance data of samples.  
- *Separate_training_metadata.csv*: Metadata of samples.  

1. If select **Separate metadata and abundance file**:  
   ```merge column name in metadata file``` - uuid  
   ```merge column name in abundance file``` - uuid  
2. ```Please enter the target (same as column name)``` - city  
   ```hierarchy in prediction model (separator:',' )``` - continent,city,latitude,longitude  
   ```column range of abundance data with separator as ':' (eg. 47:100)``` - 43:70  
3. If select **Remove values (eg. control values)**:  
   ```Values need to be removed with format as ```  - city:control,other_control,neg_control,other,pos_control;control_type:ctrl cities,negative_control,positive_control  
4. If select **Subsets in feature elimination**:  
   ```Values need to be removed with format as ``` - 5,15,20,25  

5. After enter all the values, please click start and click the **Accuracy** tab.  



### 3 Build a new prediction model using mGPS and predict new samples  

Similar enter information as function 2 **Build a new prediction model using mGPS**. After enter all the values, please click start and click the **Map** tab.  
