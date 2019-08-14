import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col
import numpy
from matplotlib.ticker import FormatStrFormatter
import sys

class ColorCycler:
	def __init__(self, num_colors = 8, **kwargs):
		self.number_of_colors_in_wheel = num_colors
		
		if 'cmap' in kwargs.keys():
			self.color_set = plt.get_cmap(kwargs['cmap'])
		else:
			self.color_set = plt.get_cmap('Set1')

		if 'start_index' in kwargs.keys():
			self.start_index = kwargs['start_index']
		else:
			self.start_index = 0

		if 'increment' in kwargs.keys():
			self.increment = kwargs['increment']
		else:
			self.increment = 1

		cNorm  = col.Normalize(vmin=0, vmax= self.number_of_colors_in_wheel*(1))
		self.scalarMap = cm.ScalarMappable(norm=cNorm, cmap=self.color_set)

		self.color_cycle_position = self.start_index

	def get_new_color(self, alpha = 1., increment_cycle = True):
		R, G, B, A = self.scalarMap.to_rgba(self.color_cycle_position%self.number_of_colors_in_wheel)
		if increment_cycle:
			self.color_cycle_position += self.increment
		return (R,G,B, alpha)


def logit(x):
	return numpy.log10(x/(1.0-x))

def logistic(x):
	return 10**x/(1 + 10**x)


def logitlabelentry(x):
	if x > 1:
		return r'$1 - 10^{-%d}$' % x
	elif abs(x) <= 1:
			return r'$%.1g$' % ((10**x/(1.0 + 10**x)))
	else:
		return r'$10^{%d}$' %  x


class BarChart:
	def __init__(self,figure,axes,xvals):
		self.fig = figure
		self.axes = axes
		self.xvals = xvals
		self.bar_width = (xvals[1] - xvals[0])
		self.tolerance = self.bar_width
		self.bar_lower = numpy.zeros(len(self.xvals))
		self.temp_bar = numpy.zeros(len(self.xvals))
		self.axes.set_xlim([min(xvals),max(xvals)])
		start, end = self.axes.get_xlim()
		self.axes.xaxis.set_ticks([start, 0, end])
		self.axes.xaxis.set_ticklabels([start,0,end])
		self.axes.xaxis.set_ticks(numpy.arange(start,end,1),minor = True)
		# Hide the right and top spines
		self.axes.spines['right'].set_visible(False)
		self.axes.spines['top'].set_visible(False)
		self.axes.xaxis.set_major_formatter(FormatStrFormatter('%g'))
		self.axes.yaxis.set_major_formatter(FormatStrFormatter('%g'))
		# Only show ticks on the left and bottom spines
		self.axes.yaxis.set_ticks_position('left')
		self.axes.xaxis.set_ticks_position('bottom')

	def add_bar(self,xval,x_lower,x_upper,height,color,no_error_bar = False,outline = False,update = True):
		# adds bar of appropriate height and color at xval, 
		# and plots an error bar with the 95% confidence interval for xval
		self.temp_bar = numpy.zeros(len(self.xvals))
		spot = numpy.arange(len(self.xvals))[abs(self.xvals - xval) == min(abs(self.xvals - xval))][0]
		if outline:
			lw = 1
		else:
			lw = 0

		self.axes.bar(self.xvals[spot],height,width = self.bar_width,bottom = self.bar_lower[spot],color = color ,align = 'center',lw = lw, edgecolor = (0.0,0.0,0.0,1.0))

		if not no_error_bar:
			errorbar_y = (self.bar_lower[spot]*2 + height)/2.
			self.add_error_bar(xval,x_lower,x_upper,errorbar_y,color,outline)
		if update:
			self.bar_lower[spot] += height

		if max(self.bar_lower) > 1.0:
			sys.stderr.write('bar exceeded 1.0\n')

	def add_lineage(self, lineage,  use_epoch, environment = "evolution", color = None,no_error_bar = False,outline = False,scale = None,update = True):
		if color == None:
			color = lineage.color
		xval = lineage.relative_fitness[environment][use_epoch]
		x_lower = lineage.relative_fitness_CI[environment][use_epoch][0]
		x_upper = lineage.relative_fitness_CI[environment][use_epoch][1]
		height = lineage.freqs[11*(use_epoch) +5]
		if scale is not None:
			height = height/scale[11*(use_epoch) + 5]
		self.add_bar(xval,x_lower,x_upper,height,color,no_error_bar,outline, update )

	def add_error_bar(self,xval,x_lower,x_upper,yval,color,outline):
		# plots an error bar with the 95% confidence interval for xval
		self.axes.errorbar([xval], yval, xerr=[[xval - x_lower]*1,[x_upper - xval+10**-5]*1], fmt='--',color = color)
	def clear(self):
		# clears lower bar info
		self.bar_lower = numpy.zeros(len(self.xvals))

	def save(self,filename,**kwargs):
		if 'dpi' in kwargs.keys():
			pass
		else:
			kwargs.update({'dpi':400})
		self.fig.savefig(filename,bbox_inches = 'tight',**kwargs)

class MullerDiagram:
	def __init__(self, figure, axes,xvals):
		self.fig = figure
		self.axes = axes
		self.xvals = xvals
		self.axes.set_ylim([0.0,1.0])
		self.axes.set_xlim([0, max(self.xvals)])

	def place(self,lineage,environment = "evolution",ignore_epochs = 0,color = None, single_epoch = False, outline = False,zorder = 0, scale = None,**kwargs):
		if color is None:
			color = lineage.color

		if environment == "evolution":
			begin = ignore_epochs*11
			single_epoch_length = 10
		else:
			begin  = ignore_epochs*11 + 10
			single_epoch_length = 1


		xvals = self.xvals[begin:]
		lower = lineage.muller_lower[begin:]

		if scale is None:	
			upper = lineage.muller_lower[begin:] + lineage.freqs[begin:]
		else:
			upper = lineage.muller_lower[begin:] + lineage.freqs[begin:]/scale[begin:]
		if single_epoch:
			xvals = xvals[:single_epoch_length+1]
			lower = lower[:single_epoch_length+1]
			upper = upper[:single_epoch_length+1]
		if outline:
			lw = 1.
		else:
			lw = 0.0
		self.axes.fill_between(xvals,lower, upper, lw = lw, facecolor = color,zorder=zorder,**kwargs)

	def save(self,filename,**kwargs):
		sys.stderr.write('saving figure to '+ filename + '\n')
		if 'dpi' in kwargs.keys():
			pass
		else:
			kwargs.update({'dpi':300})
		self.fig.savefig(filename,bbox_inches = 'tight', **kwargs)
