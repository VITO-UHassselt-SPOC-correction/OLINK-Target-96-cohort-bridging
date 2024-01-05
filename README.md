# OLINK-Target-96-cohort-bridging procedure


## Running the app online
shiny app to be deployed based on IRI shiny app from murih
## Running the app locally

1. Open the terminal/command promt on your pc.
2. change the directory to the folder you want to store the application in.
```bash
cd Documents\OLINK-Target-96-cohort-bridging
```
3. Clone this GitHub repository into your folder.
```bash
git clone https://github.com/driesheylen123/OLINK-Target-96-cohort-bridging.git
```
4. The repository includes the required reference files but in orde to have your own cohort data bridged towards this reference pool you need to insert (each batch) of the cohort data in raw npx format at the "./data/" level.
   

5. Once your data is in the right place the last thing you need to do before you can run the script is insert the correct filename in the bridging_run.R file where "./data/yourdata_NPX.xlsx"




6. Now running the bridging_run.R file in any R environemnt will apply the cohort bridging procedure on you data to make it comparable. In the ".data/output" directory the complete bridged dataset will be saved as well as a file called 'intervals.csv' that computes an upper and lower interval bound for each protein based on your data, indicating the variation present in your cohort for that specific protein.

