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

import variable_writer

db = database.GlycanCompositionDiskBackedStructureDatabase("./analysis/hypothesis/reduced-permethylated.db", 1)


ads = serialize.AnalysisDeserializer("./analysis/results/serum-rp-unregularized.db")
ads2 = serialize.AnalysisDeserializer("./analysis/results/serum-rp-prior.db")
ads3 = serialize.AnalysisDeserializer("./analysis/results/serum-rp-grid.db")


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

ads_map = {
    "combinatorial_unregularized": "./analysis/results/serum-rp-unregularized.db",
    "combinatorial_partial": "./analysis/results/serum-rp-prior.db",
    "combinatorial_grid": "./analysis/results/serum-rp-grid.db",
    "krambeck_unregularized": "./analysis/test/serum-rp-unregularized.db",
    "krambeck_partial": "./analysis/test/serum-rp-prior.db",
    "krambeck_grid": "./analysis/test/serum-rp-grid.db",
    # "glyspace_unregularized": "./analysis/test/serum-rp-unregularized-glyspace.db",
    # "glyspace_partial": "./analysis/test/serum-rp-prior-glyspace.db",
    # "glyspace_grid": "./analysis/test/serum-rp-grid-glyspace.db"
    "glyspace_unregularized": "./analysis/test/serum-rp-unregularized-humanish.db",
    "glyspace_partial": "./analysis/test/serum-rp-prior-humanish.db",
    "glyspace_grid": "./analysis/test/serum-rp-grid-humanish.db"
}

ads_map = {k: serialize.AnalysisDeserializer(v) for k, v in ads_map.items()}


print "Loading data"

gcs_map = {
    # k: v.load_glycan_composition_chromatograms()
    k: trace.ChromatogramFilter(
        v.query(serialize.GlycanCompositionChromatogram).all())
    for k, v in ads_map.items()
}

datasets = {}
used_as_adduct = set()
for group, members in gcs_map.items():
    universal_list = dict()
    for case in db:
        match = members.find_all_instances(case)
        if not match:
            universal_list[(case)] = 0
        else:
            assert len(match) == 1
            if match[0].used_as_adduct:
                used_as_adduct.add(match[0].glycan_composition)
            universal_list[(case)] = match[0].score
    datasets[group] = universal_list
datasets = pd.DataFrame(datasets)
true_positive_mask = datasets.index.isin(true_positives)
detected = datasets.sum(axis=1) != 0
datasets = datasets[detected]
true_positive_mask = true_positive_mask[detected]

not_used_as_adduct_mask = ~datasets.index.isin(used_as_adduct) & true_positive_mask

glyspace_identified = set(datasets[not_used_as_adduct_mask & (datasets.glyspace_unregularized > 0)].index)
rest_identified = set(datasets[not_used_as_adduct_mask & (
    datasets.combinatorial_partial > 0)].index) - glyspace_identified

with open("serum_glycans_not_in_glyspace.txt", 'w') as fh:
    for gc in (map(strip_derivatization, map(glypy.GlycanComposition.parse, rest_identified))):
        gc.reducing_end = None
        fh.write("%s\n" % gc)

print("Drawing ROC")
for group in datasets:
    fpr, tpr, _ = roc_curve(true_positive_mask, datasets[group])
    plt.plot(fpr, tpr, alpha=0.8, label="%s (%0.3f AUC)" % (' '.join(group.split("_")).title(), roc_auc_score(
        true_positive_mask, datasets[group])))
plt.plot([0, 1], [0, 1], color='black', linestyle='--', lw=0.8, alpha=0.8)
plt.xlabel("FPR", fontsize=16)
plt.ylabel("TPR", fontsize=16)
plt.xlim(-0.01, 1)
plt.ylim(0, 1.01)
plt.title("ROC Curve", fontsize=18)
ax = plt.gca()
plt.legend(loc=4)
ax.axes.spines['right'].set_visible(False)
ax.axes.spines['top'].set_visible(False)
plt.savefig("figure/serum_roc.pdf", bbox_inches='tight')


print("Drawing PR")
plt.figure()
for group in datasets:
    prec, rec, _ = precision_recall_curve(true_positive_mask, datasets[group])
    plt.step(rec, prec, alpha=0.8, label="%s (%0.3f AUC)" % (' '.join(group.split("_")).title(), sk_auc(
             rec, prec)))
plt.legend()
plt.xlabel("Recall", fontsize=16)
plt.ylabel("Precision", fontsize=16)
plt.title("Precision-Recall Curve", fontsize=18)
ax = plt.gca()
ax.axes.spines['right'].set_visible(False)
ax.axes.spines['top'].set_visible(False)
plt.savefig("figure/serum_prec_rec.pdf", bbox_inches='tight')


print "Generating Tetra-antennary Plot"
labeler = plotting.NGlycanLabelProducer(
    glypy.GlycanComposition.parse("{Fuc^Me:1; Hex^Me:7; HexNAc^Me:6; Neu5NAc^Me:1}$C1H4"))
color_cycle = iter(("red", "green", "blue", "orange")).next

art = plotting.SmoothingChromatogramArtist(
    [gcs_map['combinatorial_partial'].find_key("{Fuc^Me:1; Hex^Me:7; HexNAc^Me:6; Neu5NAc^Me:1}$C1H4"),
     gcs_map['combinatorial_partial'].find_key("{Fuc^Me:1; Hex^Me:7; HexNAc^Me:6; Neu5NAc^Me:2}$C1H4"),
     gcs_map['combinatorial_partial'].find_key("{Fuc^Me:1; Hex^Me:7; HexNAc^Me:6; Neu5NAc^Me:3}$C1H4"),
     gcs_map['combinatorial_partial'].find_key("{Fuc^Me:1; Hex^Me:7; HexNAc^Me:6; Neu5NAc^Me:4}$C1H4")],
    colorizer=lambda *a, **k: color_cycle()).draw(
        label_function=labeler)

art.ax.set_xlim(*map(sum, zip(art.ax.get_xlim(), (-3, 2))))
art.ax.annotate(
    "Low Abundance Form\nSaved By Smoothing",
    xy=(art[0].apex_time, max(art[0].as_arrays()[1])), xytext=(art[0].apex_time, 2e7),
    ha='center', arrowprops=(dict(facecolor='black', shrink=0.05, width=2, headwidth=10)),
    fontsize=12)
# art.ax.set_title("Smoothing Recovers\nLow Abundance Features", fontsize=20)
art.ax.set_xlabel("Retention Time (Min)", fontsize=16)
art.ax.set_ylabel("Relative Abundance", fontsize=16)
axbox = art.ax.get_position()
art.ax.legend(
    # loc=(axbox.x0 - 0.1, axbox.y0 + 0.7),
    loc=2,
    ncol=2)
art.ax.legend_.set_title(labeler.label_key.replace("^Me", ""))
art.ax.figure.set_figheight(3.5)
art.ax.figure.savefig("figure/fucosylated_tetra_antennary_structures.pdf", bbox_inches='tight')


print "Generating Ammonium Adduct Plot"
art = plotting.SmoothingChromatogramArtist([
    gcs_map['combinatorial_partial'].find_key("{Hex^Me:6; HexNAc^Me:5; Neu5NAc^Me:3}$C1H4"),
    gcs_map['combinatorial_partial'].find_key("{Fuc^Me:1; Hex^Me:7; HexNAc^Me:5; Neu5NAc^Me:2}$C1H4"),
    gcs_map['combinatorial_partial'].find_key("{Fuc^Me:2; Hex^Me:8; HexNAc^Me:5; Neu5NAc^Me:1}$C1H4"),
]).draw(label_function=labeler)
case = art[0]
x = case.apex_time + 0.1
y1 = np.max(case.as_arrays()[1])
y2 = np.max(art[1].as_arrays()[1])
y3 = np.max(art[2].as_arrays()[1])
art.ax.annotate("Ammonium\nAdduct", xy=(x, y2 * 1.0), xytext=(x, (y1 + y2) / 2), ha='center',
                arrowprops=(dict(facecolor='black', shrink=0.05, width=2, headwidth=10)), fontsize=12)
# art.ax.annotate("", xy=(x, ((y1 + y2) / 2) * 1.15), xytext=(x, y1 * .95),
#                 arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=10))
art.ax.annotate("", xy=(x, y3 * 1.05), xytext=(x, y2 * .95),
                arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=10))
art.ax.legend_.remove()
art.ax.legend(loc=1, ncol=2, fontsize=10)
art.ax.legend_.set_title(labeler.label_key.replace("^Me", ""))

art.ax.set_xlabel("Retention Time (Min)", fontsize=16)
art.ax.set_ylabel("Relative Abundance", fontsize=16)
art.ax.figure.set_figheight(3.5)
art.ax.figure.savefig("figure/ammonium_adduct_ambiguity.pdf", bbox_inches='tight')

print "Generating EIC plot"
und = trace.ChromatogramFilter(
    ads2.query(serialize.UnidentifiedChromatogram).filter(
        serialize.UnidentifiedChromatogram.analysis_id == ads2.analysis_id).all())


summary_plot = summaries.GlycanChromatographySummaryGraphBuilder(
    filter(lambda x: x.score > 5, gcs_map['krambeck_partial'] + und))
lcms_plot, composition_abundance_plot = summary_plot.draw(min_score=5)
lcms_plot.ax.legend_.set_visible(False)

# lcms_plot.ax.set_title("Aggregated Extracted Ion Chromatograms", fontsize=18)
lcms_plot.ax.set_xlabel("Retention Time (Min)", fontsize=14)
lcms_plot.ax.set_ylabel("Relative Abundance", fontsize=14)

fig = lcms_plot.ax.figure
# fig.set_figwidth(fig.get_figwidth() * 2.)
# fig.set_figheight(2ig.get_figheight() * 2.)
lcms_plot.ax.figure.set_figheight(2.5)
lcms_plot.ax.figure.set_figwidth(12)

tls = lcms_plot.ax.xaxis.get_ticklabels()
for tl in tls:
    tl.set(fontsize=14)
lcms_plot.ax.xaxis.set_ticklabels(tls)
lcms_plot.ax.tick_params(axis='y', labelsize=14)

fig.savefig("figure/rp_serum_chromatograms.pdf", bbox_inches='tight')


variables = variable_writer.VariableCollection('Serum')
for group in datasets:
    score = roc_auc_score(true_positive_mask, datasets[group])
    label = group.title().replace("_", '')
    variables[label + 'ROCAUC'] = round(score, 3)
    variables[label + 'Total'] = (datasets[group][true_positive_mask] > 5.0).sum()
    variables[label + 'TotalSimplified'] = (datasets[group][not_used_as_adduct_mask] > 5.0).sum()
variables.write()
