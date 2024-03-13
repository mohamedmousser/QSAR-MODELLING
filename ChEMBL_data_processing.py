'''
Program for pre-processing the raw data from ChEMBL, intending to build a Qualitative SAR
Please, read the 'README.md' file before using this program
'''

import pandas as pd
import os as os
import numpy as np

os.chdir(os.getcwd())

def rename_columns(df, column_to_rename, new_name):
    """
    If yu want to change the name of a column, note that for multiple column names changing just use a comma ',' and give
    the other column name and new name as for the first one
    @param df: dataframe to search in
    @param column_to_rename: 'Name of the column to rename'
    @param new_name: 'New desired name of the column'
    @return:
    """
    return df.rename(columns={column_to_rename: new_name}, inplace=True)


def drop_col(df, col_to_drop):
    """

    @param df: df to search in
    @param col_to_drop: 'Name of the column to be dropped'
    @return:
    """
    return (df.drop(columns=[col_to_drop], inplace=True))


def search_index(df, column_to_check, tracked_value):
    """
    If you want to find line numbers containing a certain value in a certain column :D
    @param df: dataframe to search in
    @param column_to_check: 'column name'
    @param tracked_value: the cell to search
    @return:
    """
    return df.loc[df[column_to_check] == tracked_value]

# Requesting the abbrev of the kinase target
TARGET = str(input('Enter the abbrev of the kinase target please '))

# load our csv input file
df = pd.read_csv(TARGET + '.csv', sep=';')

# rename the first column as CHEMBL_ID
rename_columns(df, 'Molecule ChEMBL ID', 'CHEMBL_ID')

# set variables for replacing simply by calling them
old = float('nan')
no_value = pd.NaT

# set cells to eliminate as containing missing values to eliminate them (with dropna)
to_eliminate = ["'<'", "'<='", "'>='", "'~'", old]
for el in to_eliminate:
    df.replace(el, no_value, inplace=True)

# executing the drop for our undesired lines
df.dropna(subset=['Standard Relation', 'Standard Value', 'Standard Units', 'Document Journal', 'Document Year'],
          inplace=True)

# Convert all columns that contains numbers to numeric series
ls_num_cols = ['Molecular Weight', 'Standard Value', 'Document Year']
df[ls_num_cols] = df[ls_num_cols].apply(pd.to_numeric)

# convert the uM IC50 to nM
df.loc[df['Standard Units'] == 'uM', 'Standard Value'] = df['Standard Value'] * 1000

# convert the ug.mL-1 IC50 to nM (plz check the convertion operation)
df.loc[df['Standard Units'] == 'ug.mL-1', 'Standard Value'] = (df['Standard Value'] / df['Molecular Weight']) * 1000000

# select the cells of standard units for which the standard value was converted and replace them by nM
selection_cond_to_replace_units = (df['Standard Units'] == 'ug.mL-1') | (df['Standard Units'] == 'uM')
df.loc[selection_cond_to_replace_units, 'Standard Units'] = 'nM'

# select '>>' and '>' relation cells with IC50 values under 10 000 and replace them by NaT to remove them by dropna
selection_cond = ((df['Standard Relation'] == "'>'") | (df['Standard Relation'] == "'>>'")) & (
        df['Standard Value'] < 10000)
df.loc[selection_cond, 'Standard Value'] = no_value
df.dropna(subset=['Standard Value'], inplace=True)

# Drop some unwanted columns
drop_col(df, 'Compound Key')
drop_col(df, 'Molecule Max Phase')
drop_col(df, '#RO5 Violations')
drop_col(df, 'Assay Organism')
drop_col(df, 'BAO Label')
drop_col(df, 'Assay Tissue ChEMBL ID')
drop_col(df, 'Assay Tissue Name')
drop_col(df, 'Assay Cell Type')
drop_col(df, 'Assay Subcellular Fraction')
drop_col(df, 'Assay Parameters')
drop_col(df, 'Assay Variant Accession')
drop_col(df, 'Assay Variant Mutation')
drop_col(df, 'Target ChEMBL ID')
drop_col(df, 'Target Type')
drop_col(df, 'Document ChEMBL ID')
drop_col(df, 'Source ID')
drop_col(df, 'Cell ChEMBL ID')
drop_col(df, 'Properties')
drop_col(df, 'pChEMBL Value')
drop_col(df, 'Data Validity Comment')
drop_col(df, 'Comment')
drop_col(df, 'Ligand Efficiency BEI')
drop_col(df, 'Ligand Efficiency LE')
drop_col(df, 'Ligand Efficiency LLE')
drop_col(df, 'Ligand Efficiency SEI')
drop_col(df, 'Potential Duplicate')
drop_col(df, 'AlogP')
drop_col(df, 'Assay Description')
drop_col(df, 'Uo Units')
drop_col(df, 'Assay ChEMBL ID')
drop_col(df, 'Assay Type')
drop_col(df, 'BAO Format ID')
drop_col(df, 'Target Name')
drop_col(df, 'Target Organism')
drop_col(df, 'Source Description')
drop_col(df, 'Action Type')
drop_col(df, 'Standard Text Value')



# Insert the activity column
l = []
(i, j) = df.shape
for n in range(i):
    l.append(pd.NaT)
df.insert(loc=1, column='Activity', value=np.array(l))

# Group CHEMBL_ID's and take for each one the minimal value of IC50
grouped = df.groupby('CHEMBL_ID')['Standard Value'].min()

# Replace the IC50 values by the minimal ones for each CHEMBL_ID in df
for cbl in grouped.keys():
    df.loc[df['CHEMBL_ID'] == cbl, 'Standard Value'] = grouped[cbl]

# Tell if a ligand is active, Intermediate or inactive
df.loc[(df['Standard Value'] <= 200), 'Activity'] = 'Active'
df.loc[((df['Standard Value'] > 200) & (df['Standard Value'] <= 1000)), 'Activity'] = 'Intermediate'
df.loc[(df['Standard Value'] > 1000), 'Activity'] = 'Inactive'

# Drop the repeated CHEMBL_ID rows (now it doesn't matter which one we take because all duplicates have the same ic50)
# the minimal one
df.drop_duplicates(subset=['CHEMBL_ID'], inplace=True)

# Considering the intermidates as inactives
df['Activity'] = df['Activity'].replace('Intermediate', 'Inactive')

# Get the number of ligands
ligand_number = df.shape[0]

# Write the cleaned data into a new csv output file
df.to_excel('treated_' + TARGET + '.xlsx', index=False)

print(100*'*', '\nYour raw data have been cleaned. A file named "' + 'treated_' + TARGET + '.xlsx" was successfully generated !\n', 'Please, if these programs help your scientific research, cite us : doi ...')
print(100 * '*')
