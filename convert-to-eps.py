import os
import re
import glob

exepath = r"C:\Program Files\Inkscape\inkscape.com"

has_figure = re.compile(r"\includegraphic")
figure_file = re.compile(r"\{(figure/.*\.pdf)\}")

with open('build/_supplement.tex') as fh:
    for line in fh:
        line = line.strip()
        if has_figure.search(line):
            figure_pdf = figure_file.search(line)
            if figure_pdf:
                figure_pdf = figure_pdf.group(1)
                print("%s %s --export-eps=%s" % (exepath, figure_pdf, figure_pdf[:-3] + 'eps'))
                os.system("\"%s\" %s --export-eps=%s" % (exepath, figure_pdf, figure_pdf[:-3] + 'eps'))
