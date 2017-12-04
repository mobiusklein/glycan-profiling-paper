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


print "Native Plots"
native_path = "analysis/results/sulfated/phil-82-native-fit-prior.db"
ads = serialize.AnalysisDeserializer(native_path)

und = trace.ChromatogramFilter(
    ads.query(serialize.UnidentifiedChromatogram).all())
gcs = trace.ChromatogramFilter(
    ads.query(serialize.GlycanCompositionChromatogram).all())

print len(filter(lambda x: x.score > 5 and not x.used_as_adduct, gcs + und))

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
lcms_plot.ax.figure.savefig("figure/native_phil82_chromatograms.pdf", bbox_inches='tight')
composition_abundance_plot.ax.figure.savefig(
    "figure/native_phil82_abundances.pdf", bbox_inches='tight')


print "Deuteroreduced and Permethylated Plots"
dp_path = "analysis/results/phil-82-dr-fit-prior.db"
ads = serialize.AnalysisDeserializer(dp_path)

und = trace.ChromatogramFilter(
    ads.query(serialize.UnidentifiedChromatogram).all())
gcs = trace.ChromatogramFilter(
    ads.query(serialize.GlycanCompositionChromatogram).all())

print len(filter(lambda x: x.score > 5 and not x.used_as_adduct, gcs + und))

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
lcms_plot.ax.figure.savefig("figure/dp_phil82_chromatograms.pdf", bbox_inches='tight')
composition_abundance_plot.ax.figure.savefig(
    "figure/dp_phil82_abundances.pdf", bbox_inches='tight')
