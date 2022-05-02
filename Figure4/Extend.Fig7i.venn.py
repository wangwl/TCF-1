import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
fig = plt.figure()
from matplotlib_venn import venn2
venn2(subsets=(367, 324, 3083), set_labels=('WT', 'KO'))
plt.savefig("FigureS7B.pdf")
