import glypy
from glycan_profiling.database import composition_network

with open('./combinatorial-glycans.txt') as fh:
    glycan_compositions = []
    for line in fh:
        gc = line.split("\t")[0]
        glypy.GlycanComposition.parse(gc)
        glycan_compositions.append(gc)


neighborhoods = composition_network.make_n_glycan_neighborhoods()

graph = composition_network.CompositionGraph(glycan_compositions, neighborhoods=neighborhoods)

graph.create_edges()

walker = composition_network.NeighborhoodWalker(graph, neighborhoods, assign=True)


for neighborhood, members in walker.neighborhood_maps.items():
    print(neighborhood, len(members))
