import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
fig = plt.figure()
from matplotlib_venn import venn2
venn2(subsets=(322, 504, 3198), set_labels=('Norm 3T3', 'TCF-1 3T3'))
plt.savefig("FigureS5D.pdf")
