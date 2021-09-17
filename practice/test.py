# standard lib imports
import sys
import warnings
from collections import Counter
from collections import defaultdict

# 3rd party imports
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# sklearn impots
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
                             average_precision_score,
                             roc_auc_score)
from scipy.stats import sem

# external packages in local directory
sys.path.append("../simulate-groups")
from simulate_groups import simulate_ll

# removing warning message
warnings.filterwarnings("ignore")

sample_sizes = [10, 50, 100, 500, 1000, 5000, 7500, 10000]
l1_ratios = [0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1]
random_state_vals = [0, 2, 4, 6, 8, 10] # 12, 14, 16, 18, 20]# 20 different values


for sample_size in sample_sizes:
	X, y, info_dict = simulate_ll(n=sample_size, p=20, uncorr_frac=0.1, num_groups=5)
	for idx, random_state_val in enumerate(random_state_vals):

		# splitting and trianing data
		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=random_state_val)
		model_cv = LogisticRegressionCV(cv=3, penalty="elasticnet", solver="saga", l1_ratios=l1_ratios).fit(X_train, y_train)

		# testing model
		y_pred = model_cv.predict(X_test)
		aupr_score = round(average_precision_score(y_test, y_pred), 3)
		try:
			auroc_score = round(roc_auc_score(y_test, y_pred), 3)
		except ValueError:
			# aupr_score = 0
			auroc_score = 0

		print(f"Sample size={sample_size} Random_seed={random_state_val}, Best_l1_ratio={model_cv.l1_ratio_[0]}, AUPR={aupr_score} AUORC={auroc_score}")