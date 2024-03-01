# QSAR-MODELLING
DATA &amp; MODELLING PROCEDURE.

 This is a repository to facilitate ChEMBL data curation and classification QSAR modelling procedure. Please read this paper : doi... before using the programs.
 If you use the programs, make sure you have the requirements listed below, and rename your files in the same format of the given examples in the repo.
 For the math treatment program, first make sure you have a csv file with molecular descriptors, or use our example : "Descriptors.csv" (INSR). Second, make sure you save the two columns csv file from the cleaned matrix (containing only the chembl_id and the activity labels (in our example (INSR) there is also the IC50 column). 


## Requirements
2. Clone the repository
2. install python3
3. pip install requirements.txt

## Data treatment
1. Download the csv raw data from https://ebi.ac.uk/chembl in the same directory, or use the IGF1R data example in the repo 
2. In the linux terminal, type : python3 ChEMBL_data_preprocessing.py, you will be invited to type the name of the raw data file and your new cleaned file will be generated in the same directory

## Math treatment
1. Generate your descriptors file in the same directory, or use the INSR example
2. In the linux terminal, type : python3 MPT.py, you will be invited to type the name of the descriptors file and your MPT will be done. You will find the correlation matrix and the train/test files in the a new directory called "preparation_output".

## Machine learning
1. In the linux terminal, type : python3 Machine_learning.py, you will be invited to type the abbreviation of the kinase and the adequate number of decision trees (see Mousser et al. paper). A report on the QSAR model will appear. 
