import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D
import matplotlib.cm as cm
import matplotlib.colors as col
import matplotlib.gridspec as gridspec

#Set figure parameters

matplotlib.rcParams['xtick.major.size'] = 2
matplotlib.rcParams['xtick.major.width'] = .3
matplotlib.rcParams['ytick.major.size'] = 2
matplotlib.rcParams['ytick.major.width'] = .3

matplotlib.rcParams['xtick.labelsize'] = 7
matplotlib.rcParams['ytick.labelsize'] = 7
matplotlib.rcParams['legend.numpoints'] = 1
matplotlib.rcParams['legend.frameon'] = 'False'
matplotlib.rcParams['legend.fontsize'] = 10
matplotlib.rcParams['legend.markerscale'] = 1.0
matplotlib.rcParams['legend.loc'] = 'lower right'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['patch.edgecolor'] = 'k'
# matplotlib.rcParams['patch.linewidth'] = 0.1
matplotlib.rcParams['axes.linewidth'] = 1.
matplotlib.rcParams['lines.markeredgewidth'] = 0.0    # the line width around the marker symbol
matplotlib.rcParams['lines.markersize'] = 5           # markersize, in points
matplotlib.rcParams['figure.figsize'] = (10,5)
matplotlib.rcParams['font.size'] = 12

matplotlib.rcParams['pdf.fonttype'] = 42
