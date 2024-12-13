# OLINK-Target-96-SPOC-Correction-procedure

There is a need for a  robust method that allows all OLINK Target 96 studies to compare proteomics data accurately and cost-efficiently. To meet these goals we developed the Synthetic Plasma Pool Cohort Correction, the ‘SPOC correction’ approach, based on the use of an OLINK-composed synthetic plasma sample. The method can easily be implemented in a federated data-sharing context.

## Running the correction locally

1. Open the terminal/command prompt on your pc.
2. change the directory to the folder you want to store the application in.
```bash
cd Documents\OLINK-Target-96-cohort-bridging
```
3. Clone this GitHub repository into your folder.
```bash
git clone https://github.com/driesheylen123/OLINK-Target-96-cohort-bridging.git
```
4. In order to have your own cohort data bridged need to insert (each batch) of the cohort data in raw npx format at the "./data/" level. Be aware of the fact that if two seperate cohorts are placed in the "./data/" folder, the file that is alphabetically ranked first is considered as the reference dataset and the other file(s) will be corrected towards this reference dataset. You can rename your files in order to correct the required cohort. The reference_dataset contains a small set of synthetic plasma samples that is  in "./testdata/" can be used as a general reference pool if needed.
   
5. Now running the bridging_run.R file (with the briding_fun.R file in the same directory) in any R environment will apply the cohort bridging procedure on your data to make it comparable. In the ".data/output" directory the complete bridged dataset will be saved together with a file called 'intervals.csv' that computes an upper and lower interval bound for each protein based on your SPOC corrected cohort data, indicating the variation present in your cohort for each specific protein.

6. OPTONAL: Federated_control_substraction.R contains the Code designed to isolate external control sample data from your cohort or that of your collaborators, enabling the application of the SPOC correction method in compliance with federated data principles. The path for any complete raw NPX data file that you wish to substract the external sample controls from should be placed in line 289 of the script. 


