import csv

import matplotlib
try:
    matplotlib.use("agg")
except Exception:
    pass

from glycan_profiling import (
    serialize, chromatogram_tree, plotting, trace, profiler,
    database, scoring, models, composition_distribution_model)
from glycan_profiling.plotting import summaries
import glypy
import ms_deisotope
import ms_peak_picker
import brainpy

import numpy as np
from sklearn.metrics import auc as sk_auc, roc_curve, roc_auc_score, precision_recall_curve

from matplotlib import pyplot as plt, rcParams
rcParams['figure.figsize'] = 8, 6
Formate = chromatogram_tree.Formate
Ammonium = chromatogram_tree.Ammonium
Unmodified = chromatogram_tree.Unmodified


db = database.GlycanCompositionDiskBackedStructureDatabase("./analysis/hypothesis/reduced-permethylated.db", 1)


ads = serialize.AnalysisDeserializer("./analysis/test/serum-rp-unregularized.db")
ads2 = serialize.AnalysisDeserializer("./analysis/test/serum-rp-prior.db")
ads3 = serialize.AnalysisDeserializer("./analysis/test/serum-rp-grid.db")


with open("./analysis/multiglycan-match.csv", "rb") as f:
    obs = list(csv.reader(f))
true_positives = []
false_positives = []

for row in obs:
    gc = glypy.glycan_composition.HashableGlycanComposition.parse(row[0])
    status = row[1]
    try:
        remark = row[2]
    except IndexError:
        pass
    if status == "hit":
        true_positives.append(gc)
    elif status.startswith("ambiguous with unassign"):
        false_positives.append(gc)
    elif status.startswith("ambiguous with"):
        true_positives.append(gc)
    elif status.startswith("hit"):
        true_positives.append(gc)

with open("./analysis/perm_bs_other_sane_matches.csv") as f:
    for line in f:
        if not line.strip():
            continue
        gc = glypy.glycan_composition.HashableGlycanComposition.parse(line.strip())
        true_positives.append(gc)


print("loading Unregularized")
gcs = ads.load_glycan_composition_chromatograms()
print("loading Partially Regularized")
gcs2 = ads2.load_glycan_composition_chromatograms()
print("loading Fully Regularized")
gcs3 = ads3.load_glycan_composition_chromatograms()


def is_true_match(x):
    gc = x.glycan_composition
    if gc in true_positives:
        return True
    return False


gcs_true = map(is_true_match, gcs)

fpr, tpr, thresholds = roc_curve(gcs_true, [x.score for x in gcs])
auc = roc_auc_score(gcs_true, [x.score for x in gcs])


gcs2_true = map(is_true_match, gcs2)

fpr2, tpr2, thresholds2 = roc_curve(gcs2_true, [x.score for x in gcs2])
auc2 = roc_auc_score(gcs2_true, [x.score for x in gcs2])


gcs3_true = map(is_true_match, gcs3)

fpr3, tpr3, thresholds3 = roc_curve(gcs3_true, [x.score for x in gcs3])
auc3 = roc_auc_score(gcs3_true, [x.score for x in gcs3])


plt.figure()
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
          "for Perm-BS-070111-04-Human-Serum", fontsize=16)
ax = plt.gca()
ax.axes.spines['right'].set_visible(False)
ax.axes.spines['top'].set_visible(False)
plt.savefig("figure/multiglycan_serum_roc.pdf", bbox_inches='tight')


prec, rec, _ = precision_recall_curve(gcs_true, [x.score for x in gcs])
prec2, rec2, _ = precision_recall_curve(gcs2_true, [x.score for x in gcs2])
prec3, rec3, _ = precision_recall_curve(gcs3_true, [x.score for x in gcs3])


plt.figure()
plt.step(rec, prec, label="Unregularized (%0.3f)" % (sk_auc(rec, prec)))
plt.step(rec2, prec2, label="Partially Regularized (%0.3f)" % (sk_auc(rec2, prec2)))
plt.step(rec3, prec3, label="Fully Regularized (%0.3f)" % (sk_auc(rec3, prec3)))
plt.legend()
plt.xlabel("Recall", fontsize=14)
plt.ylabel("Precision", fontsize=14)
plt.title("Precision-Recall Curve Comparing Regularization Performance\n"
          "for Perm-BS-070111-04-Human-Serum", fontsize=16)
ax = plt.gca()
ax.axes.spines['right'].set_visible(False)
ax.axes.spines['top'].set_visible(False)
plt.savefig("figure/multiglycan_serum_prec_rec.pdf", bbox_inches='tight')


print "Generating EIC plot"
und = trace.ChromatogramFilter(
    ads2.query(serialize.UnidentifiedChromatogram).filter(
        serialize.UnidentifiedChromatogram.analysis_id == ads2.analysis_id).all())


summary_plot = summaries.GlycanChromatographySummaryGraphBuilder(
    filter(lambda x: x.score > 5, gcs3 + und))
lcms_plot, composition_abundance_plot = summary_plot.draw(min_score=5)
lcms_plot.ax.legend_.set_visible(False)
# lcms_plot.ax.set_title("Glycan Composition\nLC-MS Aggregated EICs", fontsize=24)

lcms_plot.ax.set_xlabel("Retention Time (Min)", fontsize=16)
lcms_plot.ax.set_ylabel("Relative Abundance", fontsize=16)

fig = lcms_plot.ax.figure
# fig.set_figwidth(fig.get_figwidth() * 2.)
# fig.set_figheight(fig.get_figheight() * 2.)

fig.savefig("figure/multiglycan_rp_serum_chromatograms.pdf", bbox_inches='tight')