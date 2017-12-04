import sys
import os

from glycan_profiling.composition_distribution_model import GridPointSolution


def main(point_solution_files, outfile=None):
    point_solutions = []
    for solution_file in point_solution_files:
        with open(solution_file) as f:
            point = GridPointSolution.load(f)
            point_solutions.append(point)

    # n_points = len(point_solutions)
    # n_measures = len(point_solutions[0].neighborhood_names) + 2

    def label_cleaner(file_name):
        base = os.path.basename(file_name)
        label = base.replace("-grid-regularization-parameters.txt", '')
        return "{%s}" % label

    table = [
        r'\begin{table}',
        r'    \centering',
        r'    \small',
        r'    \begin{tabular}{%s}' % ("l " + 'S' * (len(point_solutions))),
        r'        \toprule',
        r"        $\tau_i$ & " + " & ".join(map(label_cleaner,
                                            point_solution_files)) + r"\\",
        r'        \midrule'
    ]
    for i, neighborhood_name in enumerate(point_solutions[0].neighborhood_names):
        row = [(" " * 8) + neighborhood_name]
        for j, point in enumerate(point_solutions):
            row.append("%0.3f" % (point.tau[i],))
        line = " & ".join(row) + r"\\"
        table.append(line)

    table.append(r"        \midrule")
    row = [r'        ${\hat \lambda}$']
    for j, point in enumerate(point_solutions):
        row.append("%0.2f" % (point.lmbda,))
    line = " & ".join(row) + r"\\"
    table.append(line)
    row = [r'        ${\hat \gamma}$']
    for j, point in enumerate(point_solutions):
        row.append("%0.2f" % (point.threshold,))
    line = " & ".join(row) + r"\\"
    table.append(line)
    table.append(r"        \bottomrule")
    table.append(r"    \end{tabular}")
    table.append(r"\end{table}")
    if outfile is None:
        outfile = "tau-table.tex"
    with open(outfile, 'w') as f:
        f.write('\n'.join(table))


if __name__ == '__main__':
    main(sys.argv[1:-1], sys.argv[-1])
