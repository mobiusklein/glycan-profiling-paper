import matplotlib
try:
    matplotlib.use("agg")
except Exception:
    pass
from glycan_profiling import (
    serialize, chromatogram_tree, plotting, trace, profiler,
    database, scoring, models, composition_distribution_model)
from glycan_profiling.plotting import summaries
import numpy as np
import pandas as pd
import glypy
from glypy.composition.composition_transform import strip_derivatization
from sklearn.metrics import auc as sk_auc, roc_curve, roc_auc_score, precision_recall_curve
import variable_writer
from matplotlib import pyplot as plt, rcParams
rcParams['figure.figsize'] = 8, 6

Formate = chromatogram_tree.Formate

db = database.GlycanCompositionDiskBackedStructureDatabase("./analysis/hypothesis/native.db", 1)

files = {
    "combinatorial_unregularized": "./analysis/results/sulfated/phil-bs-native-unregularized.db",
    "combinatorial_partial": "./analysis/results/sulfated/phil-bs-native-fit-prior.db",
    "combinatorial_grid": "./analysis/results/sulfated/phil-bs-native-grid.db",
    "glyspace_unregularized": "./analysis/test/phil-bs-humanish-unregularized.db",
    "glyspace_partial": "./analysis/test/phil-bs-humanish-prior.db",
    "glyspace_grid": "./analysis/test/phil-bs-humanish-grid.db",
    "krambeck_unregularized": "./analysis/test/phil-bs-krambeck-unregularized.db",
    "krambeck_partial": "./analysis/test/phil-bs-krambeck-prior.db",
    "krambeck_grid": "./analysis/test/phil-bs-krambeck-grid.db"
}

deserializers = {
    key: serialize.AnalysisDeserializer(value) for key, value in files.items()
}

print("Loading Datasets")
gcs_map = {
    key: chromatogram_tree.ChromatogramFilter(value.query(serialize.GlycanCompositionChromatogram))
    for key, value in deserializers.items()
}

all_keys = set()
for group, members in gcs_map.items():
    for member in members:
        all_keys.add(member.key)
true_positives = set([
    gc.glycan_composition for gc in gcs_map["combinatorial_partial"] + gcs_map[
        "krambeck_grid"
    ]
    if gc.glycan_composition['Neu5Ac'] == 0])
true_positives.add("{Hex:5; HexNAc:4; Neu5Ac:1}")

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

glyspace_identified = set(datasets[true_positive_mask & (datasets.glyspace_unregularized > 0)].index)
rest_identified = set(datasets[true_positive_mask & (
    datasets.combinatorial_partial > 0)].index) - glyspace_identified

print "Writing missing glycan list"

with open("philbs_glycans_not_in_glyspace.txt", 'w') as fh:
    for gc in (map(strip_derivatization, map(glypy.GlycanComposition.parse, rest_identified))):
        fh.write("%s\n" % gc)


fig = plt.figure()
print("Drawing ROC")
for group in datasets:
    fpr, tpr, _ = roc_curve(true_positive_mask, datasets[group])
    plt.plot(fpr, tpr, alpha=0.8, label="%s (%0.3f AUC)" % (' '.join(group.split("_")).title(),
             round(roc_auc_score(true_positive_mask, datasets[group]), 3)))
plt.plot([0, 1], [0, 1], color='black', linestyle='--', lw=0.8, alpha=0.8)

plt.xlabel("FPR", fontsize=16)
plt.ylabel("TPR", fontsize=16)
plt.legend()
plt.xlim(-0.01, 1)
plt.ylim(0, 1.01)
plt.title("Receiver-Operator Characteristic Curve", fontsize=18)
ax = plt.gca()
ax.axes.spines['right'].set_visible(False)
ax.axes.spines['top'].set_visible(False)
fig.set_figwidth(8)
fig.set_figheight(6)
plt.savefig("figure/sulfated_phil_bs_native_roc.pdf", bbox_inches='tight')

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
plt.savefig("figure/sulfated_phil_bs_native_prec_rec.pdf", bbox_inches='tight')


print("Drawing EIC")
und = chromatogram_tree.ChromatogramFilter(
    deserializers['combinatorial_partial'].query(serialize.UnidentifiedChromatogram).filter(
        serialize.UnidentifiedChromatogram.analysis_id == deserializers['combinatorial_partial'].analysis_id).all())


summary_plot = summaries.GlycanChromatographySummaryGraphBuilder(
    filter(lambda x: x.score > 5, gcs_map['combinatorial_partial'] + und))
lcms_plot, composition_abundance_plot = summary_plot.draw(min_score=5)
lcms_plot.ax.legend_.set_visible(False)
lcms_plot.ax.set_title("Aggregated Extracted Ion Chromatograms", fontsize=18)
lcms_plot.ax.set_xlabel("Retention Time (Min)", fontsize=16)
lcms_plot.ax.set_ylabel("Relative Abundance", fontsize=16)

fig = lcms_plot.ax.figure

fig.savefig("figure/sulfated_phil_bs_native_chromatograms.pdf", bbox_inches='tight')
composition_abundance_plot.ax.set_title("Total Abundance", fontsize=18)
composition_abundance_plot.ax.set_xlabel(
    composition_abundance_plot.ax.get_xlabel().replace("@sulfate", 'SO3'), fontsize=16)
composition_abundance_plot.ax.set_ylabel(
    "Relative Abundance", fontsize=16)
xtick_labels = composition_abundance_plot.ax.get_xticklabels()
for label in xtick_labels:
    label.set_fontsize(7.5)
composition_abundance_plot.ax.set_xticklabels(xtick_labels)

fig = composition_abundance_plot.ax.figure
fig.savefig("figure/sulfated_phil_bs_native_abundances.pdf", bbox_inches='tight')

variables = variable_writer.VariableCollection('PhilBSStats')
for group in datasets:
    score = roc_auc_score(true_positive_mask, datasets[group])
    label = group.title().replace("_", '')
    variables[label + 'ROCAUC'] = round(score, 3)
    variables[label + 'Total'] = (datasets[group][true_positive_mask] > 5.0).sum()
variables.write()
