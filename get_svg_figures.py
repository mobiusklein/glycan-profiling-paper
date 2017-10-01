import os
import sys
import glob
import json
from urllib import unquote
from lxml import etree


def parse_html(path):
    parser = etree.HTMLParser(encoding='utf-16')
    tree = etree.parse(path, parser)
    return tree


def get_chromatogram_figure(tree):
    uri = unquote(tree.find(
        ".//section[@id='lc-ms-overview']/img").attrib['src'])
    return uri.split(",", 1)[1]


def get_abundance_plot(tree):
    uri = unquote(tree.find(
        ".//section[@id='composition-abundance-overview']/img").attrib['src'])
    return uri.split(",", 1)[1]


loc = os.path.dirname(__file__)
search_dir = os.path.join(loc, "analysis", "results")
figure_dir = os.path.join(loc, "figure")


def main(sample_name):
    query = os.path.join(search_dir, "%s*-prior-report.html" % sample_name)
    files = glob.glob(query)
    if not files:
        raise ValueError("Could not find files for %r" % query)
    else:
        tree = parse_html(files[0])
        name_map = json.load(open(os.path.join(loc, "prefix_map.json")))
        prefix = name_map[sample_name]
        chrom_svg = get_chromatogram_figure(tree)
        with open(os.path.join(figure_dir, "%s_chromatograms.svg" % prefix), 'w') as f:
            f.write(chrom_svg)
        abund_svg = get_abundance_plot(tree)
        with open(os.path.join(figure_dir, "%s_abundances.svg" % prefix), 'w') as f:
            f.write(abund_svg)


if __name__ == '__main__':
    try:
        name = sys.argv[1]
        main(name)
    except IndexError:
        name_map = json.load(open(os.path.join(loc, "prefix_map.json")))
        for key in name_map:
            main(key)
