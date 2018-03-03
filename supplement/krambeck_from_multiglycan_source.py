import csv
from glypy.structure.glycan_composition import HashableGlycanComposition


with open("./Default_Combination_V2.csv") as fh, open("krambeck_glycan_compositions.txt", 'wb') as out:
    reader = csv.reader(fh)
    # skip header
    next(reader)
    for i, row in enumerate(reader):
        if i % 100 == 0:
            print("%d glycan compositions processed" % (i,))
        gc = HashableGlycanComposition({
            k: v for k, v in dict(HexNAc=row[0], Hex=row[1], Fuc=row[2], NeuAc=row[3]).items() if v > 0})
        out.write("%s\tN-Linked\n" % (gc,))
print("%d glycan compositions parsed" % (i,))
