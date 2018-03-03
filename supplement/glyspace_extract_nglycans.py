from collections import Counter
import glypy
from glypy.io import glyspace, glycoct
from glypy.structure.glycan_composition import GlycanComposition, HashableGlycanComposition


# Execute a SPARQL query against the SPARQL endpoint hosted by GlyTouCan to
# download all N-glycans, which are identified by having the motif "glycoinfo:G00026MO"
# and which were associated with a taxonomic source of any kind
result = glyspace.query("""
    SELECT DISTINCT ?saccharide ?glycoct ?motif WHERE {
        ?saccharide a glycan:saccharide .
        ?saccharide glycan:has_glycosequence ?sequence .
        ?saccharide skos:exactMatch ?gdb .
        ?gdb glycan:has_reference ?ref .
        ?ref glycan:is_from_source ?source .
        ?source glycan:has_taxon ?taxon
        FILTER CONTAINS(str(?sequence), "glycoct") .
        ?sequence glycan:has_sequence ?glycoct .
        ?saccharide glycan:has_motif ?motif .
        FILTER(?motif in (glycoinfo:G00026MO))
    }
""")

structures = []
for i, bind in enumerate(result.bindings):
    # result.vars[1] contains the RDF key for the GlycoCT Condensed encoding of the
    # glycan structure
    text = bind[result.vars[1]]
    if i % 100 == 0:
        print("Parsed %d glycan structures" % (i,))
    try:
        structure = glycoct.loads(text)
        structures.append(structure)
    except Exception as ex:
        print(i, bind[result.vars[0]], ex)
        continue


def detatch_monosaccharide_substituents(composition, substituents=None):
    if substituents is None:
        substituents = []
    if not substituents:
        return composition

    gc = GlycanComposition()

    for key, value in composition.items():
        node_substituents = [s for p, s in key.substituents()]
        counts = Counter()
        for node_sub in node_substituents:
            for sub in substituents:
                if node_sub.name == sub.name:
                    key.drop_substituent(-1, sub.name)
                    counts[sub.name] += 1
        gc[key] = value
        for sub_key, count in counts.items():
            gc['@%s' % sub_key] += count * value
    return gc


substituents_to_detatch = [
    glypy.Substituent("sulfate"),
    glypy.Substituent("phosphate")
]

compositions = set()
for i, structure in enumerate(structures):
    if i % 100 == 0:
        print("Converted %d glycan structures. %d glycan compositions." % (i, len(compositions)))
    glycan_comp = GlycanComposition.from_glycan(structure)
    glycan_comp.drop_configurations()
    glycan_comp.drop_stems()
    glycan_comp.drop_positions()
    glycan_comp = detatch_monosaccharide_substituents(glycan_comp, substituents_to_detatch)
    try:
        compositions.add(HashableGlycanComposition(glycan_comp))
    except ValueError as ex:
        print(ex)
        continue


valid_components = {"Hex", "HexNAc", "Neu5Ac", "Fuc", "@sulfate"}

filtered_compositions = set()
for i, composition in enumerate(compositions):
    if i % 100 == 0:
        print("%d glycan compositions filtered. %d glycan compositions accepted." % (i, len(filtered_compositions)))
    components = {str(k) for k in composition.keys()}
    if len(components - valid_components) > 0:
        continue
    filtered_compositions.add(composition)

print("%d glycan compositions after filtering" % (len(filtered_compositions),))

with open("glyspace_glycan_compositions.txt", 'wb') as fh:
    for gc in filtered_compositions:
        fh.write("%s\tN-Linked\n" % gc)
