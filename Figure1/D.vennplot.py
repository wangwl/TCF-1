import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
fig = plt.figure()
from matplotlib_venn import venn2
venn2(subsets=(39183, 45796, 11021), set_labels=('CTCF+TCF1-', 'CTCF-TCF1+'))
plt.savefig("Figure1D.venn.pdf")
