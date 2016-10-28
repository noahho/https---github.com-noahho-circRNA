import sys
sys.path.append("feature")
from intervaltree import IntervalTree
from sklearn import svm
from parse import *
import random
import time
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sklearn.model_selection import cross_val_score
from Bio import SeqIO
import features

print("Fetching the training data...")
# Fetch the data
crnas = []
not_crnas = []

POSITIVE_EXAMPLES_FILE = 'data/positive.fa';
NEGATIVE_EXAMPLES_FILE = 'data/negative.fa';

fasta = SeqIO.parse(POSITIVE_EXAMPLES_FILE, "fasta")
for record in fasta:
    crnas.append([str(record.seq), record.name])
	
fasta = SeqIO.parse(NEGATIVE_EXAMPLES_FILE, "fasta")
for record in fasta:
    not_crnas.append([str(record.seq), record.name])
	
del not_crnas[len(crnas):]
	
# Concatenate data
data = crnas + not_crnas
labels = [1] * len(crnas) + [0] * len(not_crnas)

true_inputs = len(crnas)
false_inputs = len(not_crnas)
print("## Loaded "+str(len(crnas)) + " positive examples and "+str(len(not_crnas)) + " negatives")

print("Building features...")
# Build features

features = features.get_features(data)
	
# SVM
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
import numpy as np
from sklearn import metrics
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import f1_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.model_selection import ShuffleSplit

X = np.array(features)
y = np.array(labels)
ss = ShuffleSplit(n_splits=10, test_size=0.1, random_state=0)
i = 0;
scores = []
f1scores = []
auc = []
prs = []
for train, test in ss.split(X):
	X_train, X_test, y_train, y_test = X[train], X[test], y[train], y[test]
	clf = GradientBoostingClassifier(n_estimators=100, learning_rate=0.5, random_state=0).fit(X_train, y_train)  
	#clf = SVC(probability = True)
	#clf = clf.fit(X_train, y_train)  
	test_predict = clf.predict(X_test);
	test_score = clf.score(X_test, y_test)
	test_predict_proba = clf.predict_proba(X_test)[:,1]
	
	print "Set "+str(i)+" with score: "+str(test_score)
	
	fpr, tpr, thresholds = metrics.roc_curve(y_test, test_predict, pos_label=1)
	roc = metrics.auc(fpr, tpr)
	
	precision, recall, thresholds = precision_recall_curve(y_test, test_predict_proba)
	pr = metrics.auc(recall, precision)
	
	f1 = f1_score(y_test, test_predict, average='weighted')
	scores.append(test_score)
	f1scores.append(f1);
	auc.append(roc);
	prs.append(pr)
	i += 1;
aclu = np.average(np.array(prs))
aauc = np.average(np.array(auc))
af1 = np.average(np.array(f1scores))
print "ACLU: " + str(aclu) + ", Average AUC (ROC): " + str(aauc) + " Average F1 Score: "+str(af1);