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
from glypy.composition.composition_transform import strip_derivatization
import numpy as np
import pandas as pd
from sklearn.metrics import auc as sk_auc, roc_curve, roc_auc_score, precision_recall_curve

from matplotlib import pyplot as plt, rcParams
rcParams['figure.figsize'] = 8, 6
Formate = chromatogram_tree.Formate
Ammonium = chromatogram_tree.Ammonium
Unmodified = chromatogram_tree.Unmodified

files = {
    "combinatorial_unregularized": "./analysis/results/igg-native-unregularized.db",
    "combinatorial_partial": "./analysis/results/igg-native-fit-prior.db",
    "combinatorial_grid": "./analysis/results/igg-native-grid.db",
}

deserializers = {
    key: serialize.AnalysisDeserializer(value) for key, value in files.items()
}

print("Loading Datasets")
gcs_map = {
    key: trace.ChromatogramFilter(value.query(serialize.GlycanCompositionChromatogram))
    for key, value in deserializers.items()
}

all_keys = set()
for group, members in gcs_map.items():
    for member in members:
        all_keys.add(member.key)
true_positives = set([
    gc.glycan_composition for gc in gcs_map["combinatorial_partial"]])

datasets = {}
for group, members in gcs_map.items():
    universal_list = dict()
    for case in sorted(all_keys, key=lambda x: x.mass()):
        match = members.find_all_instances(case)
        if not match:
            universal_list[(case)] = 0
        else:
            assert len(match) == 1
            universal_list[(case)] = match[0].score
    datasets[group] = universal_list
datasets = pd.DataFrame(datasets)
true_positive_mask = datasets.index.isin(true_positives)
detected = datasets.sum(axis=1) != 0
datasets = datasets[detected]
true_positive_mask = true_positive_mask[detected]


# fig = plt.figure()
# print("Drawing ROC")
# for group in datasets:
#     fpr, tpr, _ = roc_curve(true_positive_mask, datasets[group])
#     plt.plot(fpr, tpr, alpha=0.8, label="%s (%0.3f AUC)" % (' '.join(group.split("_")).title(),
#              round(roc_auc_score(true_positive_mask, datasets[group]), 3)))
# plt.plot([0, 1], [0, 1], color='black', linestyle='--', lw=0.8, alpha=0.8)

# plt.xlabel("FPR", fontsize=16)
# plt.ylabel("TPR", fontsize=16)
# plt.legend()
# plt.xlim(-0.01, 1)
# plt.ylim(0, 1.01)
# plt.title("Receiver-Operator Characteristic Curve", fontsize=18)
# ax = plt.gca()
# ax.axes.spines['right'].set_visible(False)
# ax.axes.spines['top'].set_visible(False)
# fig.set_figwidth(8)
# fig.set_figheight(6)
# plt.savefig("figure/igg_native_roc.pdf", bbox_inches='tight')

print("Drawing PR")
fig = plt.figure()
for group in datasets:
    prec, rec, _ = precision_recall_curve(true_positive_mask, datasets[group])
    plt.step(rec, prec, alpha=0.8, label="%s (%0.3f AUC)" % (' '.join(group.split("_")).title(), sk_auc(
        rec, prec)))

plt.xlim(-0.01, 1)
plt.ylim(0, 1.01)
plt.legend()
plt.xlabel("Recall", fontsize=16)
plt.ylabel("Precision", fontsize=16)
plt.title("Precision-Recall Curve", fontsize=18)
ax = plt.gca()
ax.axes.spines['right'].set_visible(False)
ax.axes.spines['top'].set_visible(False)
fig.set_figwidth(8)
fig.set_figheight(6)
plt.savefig("figure/igg_native_prec_rec.pdf", bbox_inches='tight')


native_path = "analysis/results/igg-native-fit-prior.db"
ads = serialize.AnalysisDeserializer(native_path)

und = trace.ChromatogramFilter(
    ads.query(serialize.UnidentifiedChromatogram).all())
gcs = trace.ChromatogramFilter(
    ads.query(serialize.GlycanCompositionChromatogram).all())


summary_plot = summaries.GlycanChromatographySummaryGraphBuilder(
    filter(lambda x: x.score > 5 and not x.used_as_adduct, gcs + und))
lcms_plot, composition_abundance_plot = summary_plot.draw(min_score=5)
lcms_plot.ax.legend_.set_visible(False)

lcms_plot.ax.set_xlabel("Retention Time (Min)", fontsize=14)
lcms_plot.ax.set_ylabel("Relative Abundance", fontsize=14)
lcms_plot.ax.figure.set_figheight(3)
composition_abundance_plot.ax.figure.set_figheight(3)
composition_abundance_plot.ax.set_title("")

tls = lcms_plot.ax.xaxis.get_ticklabels()
for tl in tls:
    tl.set(fontsize=14)
lcms_plot.ax.xaxis.set_ticklabels(tls)
lcms_plot.ax.tick_params(axis='y', labelsize=14)
lcms_plot.ax.figure.savefig("figure/native_igg_chromatograms.pdf", bbox_inches='tight')
composition_abundance_plot.ax.figure.savefig(
    "figure/native_igg_abundances.pdf", bbox_inches='tight')
