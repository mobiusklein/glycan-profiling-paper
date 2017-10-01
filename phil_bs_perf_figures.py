import matplotlib
try:
    matplotlib.use("agg")
except Exception:
    pass
from glycan_profiling import (
    serialize, chromatogram_tree, plotting, trace, profiler,
    database, scoring, models, composition_distribution_model)

import numpy as np
from sklearn.metrics import auc as sk_auc, roc_curve, roc_auc_score, precision_recall_curve

from matplotlib import pyplot as plt, rcParams
rcParams['figure.figsize'] = 8, 6

Formate = chromatogram_tree.Formate

db = database.GlycanCompositionDiskBackedStructureDatabase("./analysis/hypothesis/native.db", 1)

ads = serialize.AnalysisDeserializer("./analysis/results/phil-bs-native-unregularized.db")
ads2 = serialize.AnalysisDeserializer("./analysis/results/phil-bs-native-fit-prior.db")
ads3 = serialize.AnalysisDeserializer("./analysis/results/phil-bs-native-grid.db")

print("loading Unregularized")
gcs = ads.load_glycan_composition_chromatograms()
print("loading Partially Regularized")
gcs2 = ads2.load_glycan_composition_chromatograms()
print("loading Fully Regularized")
gcs3 = ads3.load_glycan_composition_chromatograms()

true_positives = [x.key for x in gcs if x.glycan_composition['Neu5Ac'] == 0]
false_positives = [x.key for x in gcs if x.glycan_composition['Neu5Ac'] > 0]

fpr, tpr, thresholds = roc_curve([x.key in true_positives for x in gcs], [x.score for x in gcs])
auc = roc_auc_score([x.key in true_positives for x in gcs], [x.score for x in gcs])
fpr2, tpr2, thresholds2 = roc_curve([x.key in true_positives for x in gcs2], [x.score for x in gcs2])
auc2 = roc_auc_score([x.key in true_positives for x in gcs2], [x.score for x in gcs2])
fpr3, tpr3, thresholds3 = roc_curve([x.key in true_positives for x in gcs3], [x.score for x in gcs3])
auc3 = roc_auc_score([x.key in true_positives for x in gcs3], [x.score for x in gcs3])

fig = plt.figure()
plt.plot(fpr, tpr, alpha=0.8, label='Unregularized (%0.3f AUC)' % auc)
plt.plot(fpr2, tpr2, alpha=0.8, label='Partially Regularized (%0.3f AUC)' % auc2)
plt.plot(fpr3, tpr3, alpha=0.8, label='Fully Regularized (%0.3f AUC)' % auc3)
plt.plot([0, 1], [0, 1], color='black', linestyle='--', lw=0.8, alpha=0.8)
plt.xlabel("FPR", fontsize=14)
plt.ylabel("TPR", fontsize=14)
plt.legend()
plt.xlim(-0.01, 1)
plt.ylim(0, 1.01)
plt.title("ROC Curve Comparing Regularization Performance\n"
          "for 20141101-04-Phil-BS", fontsize=16)
ax = plt.gca()
ax.axes.spines['right'].set_visible(False)
ax.axes.spines['top'].set_visible(False)
fig.set_figwidth(8)
fig.set_figheight(6)
plt.savefig("figure/phil_bs_native_roc.pdf", bbox_inches='tight')


fig = plt.figure()

prec, rec, _ = precision_recall_curve([x.key in true_positives for x in gcs], [x.score for x in gcs])
prec2, rec2, _ = precision_recall_curve([x.key in true_positives for x in gcs2], [x.score for x in gcs2])
prec3, rec3, _ = precision_recall_curve([x.key in true_positives for x in gcs3], [x.score for x in gcs3])

plt.step(rec, prec, label="Unregularized (%0.3f AUC)" % (sk_auc(rec, prec)))
plt.step(rec2, prec2, label="Partially Regularized (%0.3f AUC)" % (sk_auc(rec2, prec2)))
plt.step(rec3, prec3, label="Fully Regularized (%0.3f AUC)" % (sk_auc(rec3, prec3)))
plt.legend()
plt.xlabel("Recall", fontsize=14)
plt.ylabel("Precision", fontsize=14)
plt.title("Precision-Recall Curve Comparing Regularization Performance\n"
          "for 20141101-04-Phil-BS", fontsize=16)
ax = plt.gca()
ax.axes.spines['right'].set_visible(False)
ax.axes.spines['top'].set_visible(False)
fig.set_figwidth(8)
fig.set_figheight(6)
plt.savefig("figure/phil_bs_native_prec_rec.pdf", bbox_inches='tight')
