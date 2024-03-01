'''
Code for math pre-treatment. Please: take a look at the 'README.md' file before using any of these programs.
Before you generate the descriptors, make sure you save the two columns csv file from the cleaned matrix
 (containing only the chembl_id and the activity labels (in our example (IGF1R) there is also the IC50 column)
'''

import pandas as pd
import numpy as np
import os

os.chdir(os.getcwd())

Activity_file = str(input('Enter the name of activity file please ("act.csv" in our example) '))
Descriptors_filePATH = str(input('Enter the name of the descriptors file please ("Descriptors.csv" in our example) '))

IC50_filePATH = Activity_file
Descriptors_filePATH = Descriptors_filePATH


""" MPT : Mathematical preprocessing!!!"""

"""Read data (descriptors and IC50)"""

print(100 * '*', '\nReading the data.')

# Read data
X_filename = 'Descriptors.csv'
y_filename = 'act.csv'
X_df = pd.read_csv(X_filename, sep=';')  # , dtype={'':float})
y_df = pd.read_csv(y_filename, sep=';')
print('Ok ..')
# y_df = y.drop(['Activity'], axis=1)
print(X_df.shape, y_df.shape)

df = pd.merge(y_df, X_df, how='inner', on='CHEMBL_ID')

"""Sort and split to train and test
*   (0,1,2) for train
*   (3) for test
"""

if 'IC50' in df.columns:
    df = df.sort_values(by=['IC50'], ascending=False)

# Drop the ic50 column
df = df.drop(['IC50'], axis=1)

print('Splitting data into train/test.')

y_train = df['Activity']

idx = list(range(0, df.shape[0]))
idx_train = [indx for indx in idx if indx % 4 != 3]
idx_test = [indx for indx in idx if indx % 4 == 3]
train_df = df.iloc[idx_train]
test_df = df.iloc[idx_test]
print('Ok ..')
print(train_df.shape, test_df.shape)

"""Remove columns that contain 'NA'"""

print('Removing columns that contain NA values.')

train_df = train_df.dropna(axis=1)
print('Ok ..')
"""Remove Constant Descriptors"""

print('Removing constant Descriptors.')

train_df = train_df.loc[:, (train_df != train_df.iloc[0]).any()]
print('Ok ..')

"""Get the column "CHEMBL_ID" values and drop it to calculate the correlation"""

print('Calculating correlation.')

train_rows_names = train_df["CHEMBL_ID"]
train_df = train_df.drop(["CHEMBL_ID"], axis=1)

X_train = train_df.drop(["Activity"], axis=1)
corr_matrix = pd.DataFrame(np.corrcoef(X_train.values, rowvar=False), columns=X_train.columns)
#print(corr_matrix.shape)
#print(corr_matrix.head)
print('Ok ..')
"""Remove columns based on correlation"""

def drop_columns_corr(df, corr_matrix, threshould_corr):
    corr_matrix = corr_matrix.abs()
    # Select upper triangle of correlation matrix
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    # Find index of feature columns with correlation greater than threshould
    to_drop = [column for column in upper.columns if any(upper[column] > threshould_corr)]
    print('Ok ..')
    print("Number of dropped columns %s\n" % str(len(to_drop)))

    # Drop features
    new_df = df.drop(to_drop, axis=1)
    return new_df, to_drop

print('Applying the correlation filter.')

threshould_corr = 0.9
train_df, _ = drop_columns_corr(train_df, corr_matrix, threshould_corr)

"""Return the molecules names to the train matrix"""

train_df.insert(0, 'CHEMBL_ID', train_rows_names)

"""Filter test dataframe"""

test_df = test_df[train_df.columns]

"""Save the train/test *dataframes* in new folder"""

folder_name = "preparation_output"

print('Creating the preparation_output directory.')

path = os.path.join(folder_name)

# Check whether the specified path exists or not
if not os.path.exists(path):
    # Create a new directory because it does not exist
    os.makedirs(path)
    print("New directory is created !")

corr_matrix.to_csv(path + "/correlation.csv", index=False)
train_df.to_csv(path + '/train_IGF1R.csv', index=False)
test_df.to_csv(path + '/test_IGF1R.csv', index=False)

print('Correlation matrix, train and test files were successfully created !')
print('Please, if this program helps your research, cite this article : doi ... ')
print(100 * '*')
