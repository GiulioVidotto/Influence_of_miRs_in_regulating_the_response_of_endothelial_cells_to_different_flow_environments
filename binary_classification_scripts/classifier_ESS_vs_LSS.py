# Import modules
import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt
import seaborn as sns
import os
import shap
from sklearn.model_selection import train_test_split, GridSearchCV, learning_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import RobustScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_curve, roc_auc_score
shap.initjs()

# 1) Import the data
# Before running the script, ensure that the working directory is the cloned github repository
file_path = os.getcwd() + "/binary_classification_scripts/ESS_vs_LSS_classifier_data_final.csv"
classifier_data = pd.read_csv(file_path)[["y", "y_1", "y_2", "y_3", "net_score", "Kd", "miRNA", "mRNA"]]

# 2) Visualize the data with an heat map
classifier_data_for_matrix = classifier_data[["net_score", "Kd", "miRNA", "mRNA"]]
corr_matrix = classifier_data_for_matrix.corr()
corr_matrix.columns = ["net effect score", "Kd", "miRNA concentration", "mRNA concentration"]
corr_matrix.index = ["net effect score", "Kd", "miRNA concentration", "mRNA concentration"]

# Plot heatmap
sns.set_theme(font_scale=1.5)  # Adjust font scale for better readability
heatmap = sns.heatmap(
    corr_matrix, 
    cmap="coolwarm",
    center=0,
    annot=True, 
    annot_kws={"size": 30},  # Annotation (number) size in heatmap
    cbar=True
)

# Customize x and y ticks
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

# Access and customize the colorbar
colorbar = heatmap.collections[0].colorbar  # Access the colorbar from the heatmap
colorbar.ax.tick_params(labelsize=30)       # Adjust the label size
plt.show()

# Divide the data in X and y
X = classifier_data.drop(["y", "y_1", "y_2", "y_3"], axis = 1)
y = classifier_data[["y", "y_1", "y_2", "y_3"]]

# Create a train/test split for the miRNAs (because there are few miRNAs and they should be unique in both train and test)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, train_size = 0.8, random_state=42)

# Starting from the y variable, obtain the information of the specific groups to test the best classifier later
filter_test_group_1 = np.where((y_test["y_2"] == 1) | (y_test["y_3"] == 1), False, True)
filter_test_group_2 = np.where((y_test["y_1"] == 1) | (y_test["y_3"] == 1), False, True)
filter_test_group_3 = np.where((y_test["y_1"] == 1) | (y_test["y_2"] == 1), False, True)

# Select only the general y to train and test the models
y_train = np.array(y_train[["y"]]).ravel()
y_test = np.array(y_test[["y"]]).ravel()

# Initialize the scaler
scaler = RobustScaler()

# Fit and transform the training data
X_train_scaled = scaler.fit_transform(X_train)
X_train_scaled = pd.DataFrame(X_train_scaled, columns=X_train.columns)  # Convert back to DataFrame
X_train = X_train_scaled

# Apply the same transformation to the test data
X_test_scaled = scaler.transform(X_test)
X_test_scaled = pd.DataFrame(X_test_scaled, columns=X_test.columns) 
X_test = X_test_scaled

# Dimenions of the X train set
print("Dimension of the X train set: {}".format(X_train.shape))

# Dimensions of the X test set
print("Dimension of the X test set: {}".format(X_test.shape))

# Dimenions of the Y train set
print("Dimension of the y train set: {}".format(y_train.shape))

# Dimensions of the Y test set
print("Dimension of the y test set: {}".format(y_test.shape))

#__________________________________________________________________________________________________________________

# LOGISTIC REGRESSION

# Train the model (logistic regression classifier)
# Define a set of hyperparameters that wil be test to find the best combination for the model
C = [0.01, 0.1, 1, 10]

# Write hyperparameters as a dictionary for elastic logistic regression
elastic_logistic_grid = {
    "C": C,
    "penalty": ["elasticnet"],
    "solver": ["saga"],
    "class_weight": ["balanced"],
    "l1_ratio" : np.arange(0.0, 1.0, 0.1),
    "max_iter" : [1000000]
}

# Search for the best combinations of hyperparameters
logistic_model = LogisticRegression()

# Elastic
elastic_logistic_grid_search = GridSearchCV(
    logistic_model,
    elastic_logistic_grid,
    cv=10,
    refit=True,
    scoring="roc_auc",
    n_jobs=-1
)

# Start the timer
start_time = time.time()

# Fit the grid search model to the training data
elastic_logistic_grid_search.fit(X_train, y_train)

# End the timer
end_time = time.time()

# Print the time taken to run the grid search
print(f"Time taken to run the grid search on the hyperparameters: {(end_time - start_time) / 60 :4f} minutes")

# Print the best hyperparameters
print(f"Best hyperparameters: \n {elastic_logistic_grid_search.best_params_}")

# Print the best model (estimator)
print(f"Best estimator: \n {elastic_logistic_grid_search.best_estimator_}")

# Get the best model from each search
elastic_logistic_best_model = elastic_logistic_grid_search.best_estimator_

# LOGISTIC REGRESSION ASSUMPTIONS

# Multicollinearity
from statsmodels.stats.outliers_influence import variance_inflation_factor

# Assuming X is your feature dataframe
vif = pd.DataFrame()
vif["Feature"] = X_train.columns
vif["VIF"] = [variance_inflation_factor(X_train.values, i) for i in range(X_train.shape[1])]
print(vif) # ALL MUST BE LOWER THAN 5 OR 10

# Linearity of Independent Variables
import numpy as np
import matplotlib.pyplot as plt

for col in X_train.columns:
    plt.figure()
    transformed_X = np.array((X_train[col])).reshape(-1, 1)

    # Calculate the logit after applying the transformation
    log_odds = elastic_logistic_grid_search.predict_proba(np.array(X_train))[:, 1]
    logit = np.log10(log_odds / (1 - log_odds))

    # Plot the transformed feature vs. logit
    plt.scatter(x = transformed_X, y = logit)
    plt.ylabel('Logit of the response variable', fontsize = 20)
    plt.xlabel(f"{col}", fontsize = 20)
    plt.show()

# Compute the predicted probabilities and residuals
pred_probs = elastic_logistic_best_model.predict_proba(X_train)[:, 1]  # Probability for the positive class
residuals = y_train - pred_probs

# Calculate leverage (diagonal elements of the hat matrix)
# Leverage can be approximated as the dot product of X and model coefficients
X_with_intercept = np.hstack([np.ones((X_train.shape[0], 1)), X_train])  # Add intercept term
hat_matrix_diag = np.sum((X_with_intercept @ np.linalg.pinv(X_with_intercept.T @ X_with_intercept)) * X_with_intercept, axis=1)

# Cook's Distance calculation
cooks_d = (residuals**2 / (2 * np.mean(residuals**2))) * (hat_matrix_diag / (1 - hat_matrix_diag)**2)

# Set a threshold for identifying outliers (influence threshold, e.g., 4/n)
n = len(X_train)
influence_threshold = 4 / n

# Plotting Cook's Distance
plt.figure(figsize=(10, 6))
plt.stem(np.arange(n), cooks_d, markerfmt=",", basefmt=" ", use_line_collection=True)
plt.axhline(y=influence_threshold, color='red', linestyle='--', label=f'Influence threshold (4/n = {influence_threshold:.4f})')
plt.xlabel('Observation index', fontsize = 20)
plt.ylabel("Cook's Distance", fontsize = 20)
plt.legend()
plt.show()

# Display potential outliers based on Cook's Distance
outliers = np.where(cooks_d > influence_threshold)[0]
print("Potential outliers (high influence):", outliers)
len(outliers)

#__________________________________________________________________________________________________________________

# GRADIENT BOOSTING

# Define the hyperparameters you want to tune for Gradient Boosting
n_estimators = [100, 200, 300, 400]  # Number of trees
learning_rate = [0.005]  # Learning rate
max_depth = [2,3]  # Maximum depth of individual trees
min_samples_split = [16,17]  # Minimum number of samples required to split a node
subsample = [0.4, 0.5]  # Fraction of samples used for fitting trees

# Write all the hyperparameters as a dictionary
gb_param_grid = [{
    "n_estimators": n_estimators,
    "learning_rate": learning_rate,
    "max_depth": max_depth,
    "min_samples_split": min_samples_split,
    "subsample": subsample
}]

# Initialize the model
gb_model = GradientBoostingClassifier(random_state=42)

# Initialize GridSearchCV with the hyperparameter grid and cross-validation
gb_grid_search = GridSearchCV(
    gb_model, 
    gb_param_grid, 
    cv=10,  # 5-fold cross-validation
    scoring="roc_auc",  # Use ROC AUC as the scoring metric
    n_jobs=-1  # Use all cores for parallel processing
)

# Start the timer
start_time = time.time()

# Fit the grid search model to the training data
gb_grid_search.fit(X_train, y_train)

# End the timer
end_time = time.time()

# Print the time taken to run the grid search
print(f"Time taken to run the grid search on the hyperparameters: {(end_time - start_time) / 60 :4f} minutes")

# Print the best hyperparameters
print(f"Best hyperparameters: \n {gb_grid_search.best_params_}")

# Print the best model (estimator)
print(f"Best estimator: \n {gb_grid_search.best_estimator_}")

# Get the best model from the grid search
gb_best_model = gb_grid_search.best_estimator_

#__________________________________________________________________________________________________________________

# RANDOM FOREST

# Train the model (random forest)
# Define a set of hyperparameters that wil be test to find the best combination for the model
n_estimators = [int(x) for x in np.linspace(start = 200, stop = 500, num= 5)]
max_depth = range(2,10)
max_features = ["sqrt"]
min_samples_leaf = range(8,10)
class_weight = ["balanced"]

# Write all the hyperparameters as a dictionary
rf_param_grid = [{
    "n_estimators" : n_estimators,
    "max_depth" : max_depth,
    "max_features" : max_features,
    "min_samples_leaf" : min_samples_leaf,
    #"min_samples_split" :min_samples_split,
    "class_weight" : class_weight
    }]

# Search for the best combinations of hyperparameters
# Initialize the model and grid search
rf_model = RandomForestClassifier(random_state=123)
rf_grid_search = GridSearchCV(
    rf_model, 
    rf_param_grid, 
    cv=10, 
    scoring="roc_auc",
    n_jobs=-1
)
start_time = time.time()
rf_grid_search.fit(X_train, y_train)
end_time = time.time()

# Print the time taken to run this search
print(f"Time taken to run the grid search on the hyperparameters: {(end_time - start_time) / 60 :4f} minutes")

# Search for the best parameters
print(f"Best hyperparameters: \n {rf_grid_search.best_params_}")

# Search for the best estimator
print(f"Best estimator: \n {rf_grid_search.best_estimator_}")
rf_best_model = rf_grid_search.best_estimator_

#__________________________________________________________________________________________________________________

# LEARNING CURVES

# Generate learning curves (training set sizes vary)
train_sizes, train_scores, val_scores = learning_curve(
    rf_best_model, X_train, y_train, cv=10, scoring='neg_log_loss', n_jobs=-1,
    train_sizes=np.linspace(0.1, 1.0, 10) # Change the model based on which one you are analysing
)

# Convert negative log loss to positive values for easier interpretation
train_scores_mean = -np.mean(train_scores, axis=1)  # Convert negative loss to positive
val_scores_mean = -np.mean(val_scores, axis=1)  # Convert negative loss to positive
train_scores_std = np.std(train_scores, axis=1)
val_scores_std = np.std(val_scores, axis=1)

# Plotting the learning curves
plt.figure(figsize=(10, 6))

# Plot training and validation losses
plt.plot(train_sizes, train_scores_mean, label='Training Loss', color='blue', marker='o')
plt.plot(train_sizes, val_scores_mean, label='Validation Loss', color='red', marker='x')

# Add shaded area for the standard deviation
plt.fill_between(train_sizes, train_scores_mean - train_scores_std, train_scores_mean + train_scores_std, color='blue', alpha=0.2)
plt.fill_between(train_sizes, val_scores_mean - val_scores_std, val_scores_mean + val_scores_std, color='red', alpha=0.2)

# Add labels and title
plt.xlabel('Training Set Size')
plt.ylabel('Log Loss')
plt.legend()
plt.grid(True)
plt.show()

#__________________________________________________________________________________________________________________

# COMPARE THE MODELS AND EVALUATE WHICH IS THE BEST ONE

# ROC CURVES

# Define the function to plot the ROC curves
def plot_roc(model_name, predict_fn, X, y):
    # Calculate ROC curve
    fprs, tprs, thresholds = roc_curve(y, predict_fn(X)[:, 1])
    # Calculate ROC AUC score
    roc_auc = roc_auc_score(y, predict_fn(X)[:, 1])
    # Plot the ROC
    plt.plot(fprs, tprs, label = f"{model_name} AUC = {round(roc_auc, ndigits = 2)}")
    plt.plot(fprs, fprs, linestyle = "--", color = "lightgrey") # Add the random classifier line
    plt.xlabel("FPR = 1 - specificity", fontsize = 30)
    plt.ylabel("TPR = sensitivity", fontsize = 30)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.01])
    plt.legend(loc = 'lower center', bbox_to_anchor = (2.5, 0.8), fontsize = 30)

# On the test data
plot_roc("Gradient Boosting classifier", gb_best_model.predict_proba, X_test, y_test)
plot_roc("Random Forest classifier", rf_best_model.predict_proba, X_test, y_test)
plot_roc("Logistic Regression classifier (Elastic-net)", elastic_logistic_best_model.predict_proba, X_test, y_test)
plt.tick_params(axis='both', which='major', labelsize=30) 
plt.show()

# PRECISION AND RECALL CURVES
# Define the function for the precision and recall curves
def plot_precision_recall_curve(model_name, predict_fn, X, y):
    # Calculate precision and recall values
    precision, recall, thresholds = precision_recall_curve(y, predict_fn(X)[:, 1])
    # Plot the ROC
    plt.plot(recall, precision, label=f"{model_name}")
    plt.plot([0, 1], [0, 0], linestyle="--", color="lightgrey")  # Baseline for no skill classifier
    plt.xlabel("Recall", fontsize = 20)
    plt.ylabel("Precision", fontsize = 20)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.title("Precision-Recall Curve")
    plt.legend(loc='lower center', bbox_to_anchor = (2.5, 0.8), fontsize = 20)

# On the test data
plot_precision_recall_curve("Gradient Boosting classifier", gb_best_model.predict_proba, X_test, y_test)
plot_precision_recall_curve("Random Forest classifier", rf_best_model.predict_proba, X_test, y_test)
plot_precision_recall_curve("Logistic Regression classifier (Elastic-net)", elastic_logistic_best_model.predict_proba, X_test, y_test)
plt.tick_params(axis='both', which='major', labelsize=20) 
plt.show()

#__________________________________________________________________________________________________________________

# Analyse the best model based on the curved and then decide which one is the best model
# The best model is...
best_model = None

# Calculate the different model matrices (specifically specificity and sensitivity. In this project we were not intersted in looking at precision, recall and f1 score)
# predict the values
y_predicted_train = best_model.predict(X_train)
y_predicted_test = best_model.predict(X_test)

# compute the probabilities
y_pred_prob = best_model.predict_proba(X_test)[:, 1]
fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)

# Calculate specificity (TNR) = 1 - FPR
specificity = 1 - fpr

# Sensitivity is the same as TPR (True Positive Rate)
sensitivity = tpr

# Find the index where the sum of sensitivity and specificity is highest
best_threshold_idx = np.argmax(sensitivity + specificity)
best_threshold = thresholds[best_threshold_idx]

# Get the best sensitivity and specificity at this threshold
best_sensitivity = sensitivity[best_threshold_idx]
best_specificity = specificity[best_threshold_idx]

# Get the best sensitivity and specificity at this threshold
best_sensitivity = sensitivity[best_threshold_idx]
best_specificity = specificity[best_threshold_idx]

# Output the results
print(f"Best Threshold: {round(best_threshold, 2)}")
print(f"Best Sensitivity (Recall): {round(best_sensitivity, 2)}")
print(f"Best Specificity: {round(best_specificity, 2)}")

#__________________________________________________________________________________________________________________

# Roc curves with only random forest for the different specific test groups
plot_roc("rf - general model", best_model.predict_proba, X_test, y_test)
plot_roc("rf - test group 1", best_model.predict_proba, X_test[filter_test_group_1], y_test[filter_test_group_1])
plot_roc("rf - test group 2", best_model.predict_proba, X_test[filter_test_group_2], y_test[filter_test_group_2])
plot_roc("rf - test group 3", rf_best_model.predict_proba, X_test[filter_test_group_3], y_test[filter_test_group_3])
plt.show()

#__________________________________________________________________________________________________________________

# Shap values
explainer = shap.Explainer(best_model)
shap_values = explainer(pd.DataFrame(X_train))

np.shape(shap_values) # Number of samples: 2634, # Number of features: 4, Number of classes: 2

# Mean shap plot
fig = shap.plots.beeswarm(shap_values[:, :, 1])

# Access the current axes from the figure
ax = fig.gcf()

# Customize the tick label size for x and y axes
ax.tick_params(axis='x', labelsize=30)  # Set x-axis label size
ax.tick_params(axis='y', labelsize=30)  # Set y-axis label size

# Show the plot
plt.show()

# Create the SHAP beeswarm plot and get the axis object
shap.plots.beeswarm(shap_values[:, :, 1])

# Now modify the font size of the axis labels
plt.gcf().tick_params(axis='x', labelsize=14)  # Adjust font size for x-axis labels
plt.gcf().tick_params(axis='y', labelsize=14)  # Adjust font size for y-axis labels

# Optionally adjust the font size of the title or other parts
plt.xlabel('Feature', fontsize=20)
plt.ylabel('SHAP Value', fontsize=20)

plt.show()

# - X axis
# A positive SHAP value means the feature is pushing the prediction towards a positive outcome (class 1)
# A negative SHAP value pushes the prediction towards the negative outcome (class 0)
# The further a point is from 0, the stronger its influence on the prediction.

# - Y axis
# Each row corresponds to one of the features
# Ordered by importance

# - Color
# The color of each dot represents the actual value of that feature for a given sample (red high values and blue low ones)

# We cannot use shap values to get to conclusions about the pathology.
# SHAP is not a measure of how important a given feature is in the real world, it is simply how important a feature is to the model.
