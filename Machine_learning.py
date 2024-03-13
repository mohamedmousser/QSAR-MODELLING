'''
Code for Xgboost machine learning over the --- data example
Please, read the 'README.md' file before using these programs.
'''

import pandas as pd
import numpy as np
import math
import xgboost as xgb
import os
from sklearn import preprocessing
# from sklearn.metrics import classification_report, confusion_matrix
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
os.chdir(os.getcwd())
# Define the csv file separator & the scaling function
separator = ','
def preapre_matrix(X, min=None, max=None):
    if (min is None):
        min = X.min(axis=0)
    if (max is None):
        max = X.max(axis=0)
    res = (X.T - min[:, None]).T
    diff = (max[:, None] - min[:, None])
    res = np.divide(res.T, diff).T
    return X, min, max

# Chose one abreviation among the given kinases training data & take its corresponding number of decision trees
# Requesting the abbrev of the kinase target
TARGET = str(input('Enter the abbrev of the kinase target please '))
nbr_dec_trees = int(input('Enter the adequate number of decision trees please (if you do not know, read this paper please : ...) '))


# Reading the train & test sets
train_filename = 'training_data/train_' + TARGET + '.csv'
test_filename = 'training_data/test_' + TARGET + '.csv'

train_df = pd.read_csv(train_filename, sep=separator)
test_df = pd.read_csv(test_filename, sep=separator)

# Print shapes
print(100*'*','\nFor our ' + TARGET + ' QSAR model we have the following train/test shapes : ', train_df.shape, test_df.shape)

# Keeping a copy of the molecules IDs
train_molecules_names = train_df['CHEMBL_ID']
test_molecules_names = test_df['CHEMBL_ID']

# Defining the train and test dfs & the descriptors names
train_df = train_df.drop(columns=['CHEMBL_ID'])
test_df = test_df.drop(columns=['CHEMBL_ID'])
descriptors_names = list(train_df.columns)
#print(train_df.columns[0])

# Now defining the ML target and features vectors
y_train = train_df['Activity'].values
X_train = train_df.drop(columns='Activity')
#train_descriptors = list(X_train.columns)

# Scaling
X_train = X_train.values
X_train, min, max = preapre_matrix(X_train)

y_test = test_df['Activity'].values
X_test = test_df.drop(columns='Activity').values
X_test, min, max = preapre_matrix(X_test, min, max)

#print('Training dataset shape:', X_train.shape, y_train.shape)
#print('Testing dataset shape:', X_test.shape, y_test.shape)

# Configure Xgboost classifier
xg_cl = xgb.XGBClassifier(colsample_bytree=0.4, learning_rate=0.15, gamma=0.3, max_depth=12, min_child_weight=1,
                          n_estimators=nbr_dec_trees)

# Encode the labels for binary classification and transform train and test targets
le = preprocessing.LabelEncoder()
le.fit(y_train)
le.classes_

y_train_transformed = le.transform(y_train)
y_test_transformed = le.transform(y_test)

# Fit the classifier to the training set
xg_cl.fit(X_train, y_train_transformed)


# Compute the evaluation metrics using the test set

# Predict the labels of the test set: preds
preds = xg_cl.predict(X_test)

# confusion matrix
tn, fp, fn, tp = confusion_matrix(y_test_transformed, preds).ravel()
cm = confusion_matrix(y_test_transformed, preds)

# Set up k-fold cross-validation
k_folds = 10  # You can adjust the number of folds
kf = KFold(n_splits=k_folds, shuffle=True, random_state=42)

# Perform k-fold cross-validation and compute scores
scores = cross_val_score(xg_cl, X_train, y_train_transformed, cv=kf,
                         scoring='accuracy')  # Replace 'accuracy' with your desired metric

# Print the cross-validation scores
print('And the following evaluation metrics : ')
print("10-fold CV :", scores.mean().round(2))

# accuracy
accuracy = (tp + tn) / (tp + tn + fp + fn)

# sensitivity (recall)
sensitivity = tp / (tp + fn)

# specificity
specificity = tn / (tn + fp)

# AUC
auc = roc_auc_score(y_test_transformed, preds)

print("Accuracy :", round(accuracy, 2))
print("Sensitivity :", round(sensitivity, 2))
print("Specificity :", round(specificity, 2))
print("AUC :", round(auc, 2))

print(100 * '*')

print('Please, if you use this program to build your model for scientific research, cite us : doi ...')
print(100 * '*')
