import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
fig = plt.figure()
from matplotlib_venn import venn2
venn2(subsets=(38390, 18341, 49194), set_labels=('CTCF+TCF-1-', 'CTCF-TCF-1+'))
plt.savefig("Figure1C.venn.pdf")
