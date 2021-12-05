# 19F_BIS335_KAIST - Biomedical Statistics & Statistical Learning project

## Goal of this project

"The breast cancer patient's data" is given. It contains clinical labels, somatic mutation status, and gene expression of the patients. The goal of this project is to identify the significant features and use them for the classification and clustering of the patients.

## HW1

+ Calculate 5-year survival probabilities of all subtypes of the breast cancer patients in the clinical dataset.

+ For each subtype of breast cancer patients, suggest highly possible driving mutation of the subtypes with estimated probability using conditional probability concept.

## HW2

+ Test whether there is a significant difference in the patient's survival time on the stage and subtypes of breast cancer. For significant changes, plot box plot to visualize the difference.

+ Test whether there is a significant difference in the patients' survival with or without mutation using student t-test and multiple testing correction

+ Using expression data, derive the differentially expressed genes in the group with different survival time in the Problem 1. Estimate the differential expression threshold of the genes, and test whether the survival of the patient is significantly different among the groups divided by the obtained thresholds.

## HW3

+ Test whether patient's survival time is related with gene expression using Pearson correlation coefficient. Suggest candidate genes which are highly correlated with patient survival

+ Perform linear regression model and genes which can predict survival time of patients well.

## HW4

+ Find the best class label for prediction of survival time among 3 class labels

	+ Perform Kaplan-Meier test with given class labels
	+ What's the best class label which minimize p-value of KM-test
	+ Perform KM-test using samples exist in all three datasets. Are there any differences from the result in b)?

## HW5

+ Obtain estimated accuracy using 5-fold cross-validation of k-nearest neighbor model of breast cancer patient classification. And find the optimal k value.

+ Repeat 5-fold cross-validation after removing class imbalance calibration step. Find the tendency of accuracy according to the k. Discuss the reason of the difference with the characteristics of K-NN algorithm.

## HW6

+ Find optimal number of predictors which minimize adjusted R-square using forward stepwise feature selection with 10-fold cross-validation.

+ Compare Ridge regression and Lasso regression model. Which one was better the predict response variable? Discuss the reason in terms of the regularization intensity.

## HW7

+ Build SVM models with changing kernel functions to classify patient stage. Compare their performances.

+ In the case of the model with radial basis kernel. Find its optimal cost and gamma parameter. Write its mathematical description and relation with model's bias-variance tradeoff.

## HW8

+ Build hierarchical clustering model using whole genes. According to silhouette index, what is optimal cluster number? Assess model by calculating Jaccard index between cluster membership and tumor stage.

+ Repeat upper question using PCA transformed data and discuss the results.