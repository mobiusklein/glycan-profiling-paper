import os
import re
import shutil

os.system("flatex formatted.tex build/_formatted.tex")
os.system("flatex supplement.tex build/_supplement.tex")
shutil.copy("bibliography.bib", "build/bibliography.bib")
shutil.copy("natbib.bst", "build/natbib.bst")

has_figure = re.compile(r"\includegraphic")
figure_file = re.compile(r"\{(figure/.*)\}")

with open("build/_formatted.tex") as fh, open("build/formatted.tex", 'w') as outfh:
    for line in fh:
        if has_figure.search(line):
            figure = figure_file.search(line)
            figure_path = figure.group(1)
            figure_name = os.path.basename(figure_path)
            line = figure_file.sub("{%s}" % figure_name, line)
            shutil.copy(figure_path, os.path.join("build", figure_name))
        outfh.write(line)


with open("build/_supplement.tex") as fh, open("build/supplement.tex", 'w') as outfh:
    for line in fh:
        if has_figure.search(line):
            figure = figure_file.search(line)
            figure_path = figure.group(1)
            figure_name = os.path.basename(figure_path)
            line = figure_file.sub("{%s}" % figure_name, line)
            shutil.copy(figure_path, os.path.join("build", figure_name))
        outfh.write(line)


os.chdir("build/")
print os.getcwd()

script = '''
pdflatex supplement.tex
bibtex supplement
pdflatex supplement.tex
pdflatex supplement.tex
pdflatex formatted.tex
bibtex formatted
pdflatex formatted.tex
pdflatex formatted.tex
'''

for line in script.splitlines():
    print line
    os.system(line)

try:
    os.remove("application-of-network-smoothing-supplement.pdf")
except OSError:
    pass

try:
    os.remove("application-of-network-smoothing.pdf")
except OSError:
    pass

os.rename("supplement.pdf", "application-of-network-smoothing-supplement.pdf")
os.rename("formatted.pdf", "application-of-network-smoothing.pdf")
