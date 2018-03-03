import matplotlib
try:
    matplotlib.use("agg")
except Exception:
    pass
from matplotlib import pyplot as plt, rcParams
import numpy as np
from scipy import stats
from glycan_profiling.scoring.chromatogram_solution import logitsum

rcParams['figure.figsize'] = 12, 8


rng = np.random.RandomState(500)

simulation_data = (rng.uniform(size=(100000, 5)).T)

simulation_scores = logitsum(simulation_data)
dist = stats.norm(loc=np.mean(simulation_scores), scale=np.std(simulation_scores))

plt.hist(simulation_scores, bins='fd', density=1)
plt.vlines(0, 0, .11, color='black', linestyle='--', label='Never Considered (p = %0.3f)' % (1 - dist.cdf(0)))
plt.vlines(5, 0, .11, color='red', linestyle='--', label='Low Confidence (p = %.3f)' % (1 - dist.cdf(5)))
plt.vlines(8, 0, .11, color='gold', linestyle='--', label='Medium Confidence (p = %.3f)' % (1 - dist.cdf(8)))
plt.vlines(15, 0, .11, color='green', linestyle='--', label='High Confidence (p = %.3e)' % (1 - dist.cdf(15)))
plt.xticks(np.arange(-20, 22, 2))
plt.xlim(-20, 20)
plt.ylim(0, 0.11)
plt.title("Simulation of Scoring Distribution", fontsize=18)
plt.ylabel("Frequency", fontsize=14)
plt.xlabel("Summarization Score", fontsize=14)
plt.text(-18, 0.09, "Normal Distribution\n $\mu: %0.3f$\n$\sigma: %0.3f$" % (
    np.mean(simulation_scores), np.std(simulation_scores)))
legend = plt.legend()
legend.set_title("Thresholds and 1 tail p-value")
plt.savefig("figure/simulation_of_scoring_distribution.pdf", bbox_inches='tight')


thresh = [0, 5, 8, 15]
p_value = 1 - dist.cdf(thresh)

print(zip(thresh, p_value))
