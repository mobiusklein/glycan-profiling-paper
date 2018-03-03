import glypy
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import matplotlib_venn
from glypy.structure.glycan_composition import HashableGlycanComposition

with open("combinatorial-glycans.txt") as combfile:
    combinatorial = []
    for line in combfile:
        combinatorial.append(HashableGlycanComposition.parse(line.split("\t")[0]))


with open("krambeck_glycan_compositions.txt") as krambeck_file:
    krambeck = []
    for line in krambeck_file:
        krambeck.append(HashableGlycanComposition.parse(line.split("\t")[0]))


with open("glyspace_glycan_compositions.txt") as glyspace_file:
    glyspace = []
    for line in glyspace_file:
        glyspace.append(HashableGlycanComposition.parse(line.split("\t")[0]))


combinatorial = set(combinatorial)
krambeck = set(krambeck)
glyspace = set(glyspace)

c = len(combinatorial)
k = len(krambeck)
g = len(glyspace)

ck = len(combinatorial & krambeck)
cg = len(combinatorial & glyspace)
kg = len(krambeck & glyspace)
ckg = len(combinatorial & krambeck & glyspace)

print("combinatorial", c)
print("krambeck", k)
print("glyspace", g)

print("combinatorial, krambeck", ck)
print("combinatorial, glyspace", cg)
print("krambeck, glyspace", kg)
print("all", ckg)


fig = plt.figure()
matplotlib_venn.venn3((combinatorial, krambeck, glyspace), set_labels=("Combinatorial", "Krambeck", "glySpace"))
ax = plt.gca()
ax.set_title("$N$-Glycan Database Overlaps")
fig.savefig("venndiagram.pdf")
