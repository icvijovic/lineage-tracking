import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col
import numpy
from matplotlib.ticker import FormatStrFormatter
import sys, itertools

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

		if 'swap_colors' in kwargs.keys():
			swap_colors = kwargs['swap_colors']
		else:
			swap_colors = []
		cNorm  = col.Normalize(vmin=0, vmax= self.number_of_colors_in_wheel*(1))
		self.scalarMap = cm.ScalarMappable(norm=cNorm, cmap=self.color_set)

		self.color_cycle_position = self.start_index

		self.index_map = {}
		for idx_pair in swap_colors:
			idx1, idx2 = idx_pair
			self.index_map.update({idx1:idx2,idx2:idx1})
		for idx in range(0,self.number_of_colors_in_wheel+1):
			if idx in self.index_map.keys():
				pass
			else:
				self.index_map.update({idx:idx})


	def map_index_to_color(self,index):
		mapped_index = self.index_map[index]
		return self.scalarMap.to_rgba(mapped_index)

	def get_new_color(self, alpha = 1., increment_cycle = True):
		R, G, B, A = self.map_index_to_color(self.color_cycle_position%self.number_of_colors_in_wheel)#self.scalarMap.to_rgba(self.color_cycle_position%self.number_of_colors_in_wheel)
		if increment_cycle:
			self.color_cycle_position += self.increment
		return (R,G,B, alpha)
	
	def get_color_by_index(self,color_index,alpha = 1.):
		R, G, B, A = self.scalarMap.to_rgba(color_index%self.number_of_colors_in_wheel)
		return (R,G,B,alpha)


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

def get_new_color(color_cycle,alpha = 1.):
	R, G, B, A = scalarMap.to_rgba(color_cycle%number_of_colors_in_wheel)
	return (R,G,B, alpha)

class TreeBranch:
	def __init__(self,lineage = None,yval= 0,parent = '',parent_y = 0):
		if lineage == None:
			self.x_begin = 0
			self.x_end = 0
			self.y = yval
			self.color = None
			self.parent = parent
			self.ID = ID
			self.parent_y = parent_y
		else:

			if lineage.ID == "":
				self.x_begin = 0
				self.parent_y = yval
			else:
				self.parent_y = parent_y
				barcode_epoch = len(lineage.ID.split("_"))
				# self.x_begin = barcode_epoch*110

				first_epoch_with_positive_fitness  = 1 + next(i for i in range(0,len(lineage.fitness)) if numpy.sum(abs(lineage.fitness[:i+1]))>5*10**-3)
				if first_epoch_with_positive_fitness == barcode_epoch:
					pos = 11*(barcode_epoch-1) 
				else:
					pos = 11*(barcode_epoch-1) + 5
				self.x_begin = 10*pos - 5 * (first_epoch_with_positive_fitness == barcode_epoch)
				
			last_time_point_with_nonzero_frequency = 1 + next(i for i in reversed(range(0,len(lineage.freqs))) if lineage.freqs[i] > 10**-4)
			self.x_end = 10 * last_time_point_with_nonzero_frequency

			self.y = yval

			self.color = lineage.color
			self.parent = parent
			self.ID = lineage.ID
			


class Tree:
	def __init__(self, figure,axes,xvals):
		self.fig = figure
		self.axes = axes
		self.xvals = xvals
		self.axes.set_xlim([min(xvals),max(xvals)])
		# Hide the right and top spines
		self.axes.spines['left'].set_visible(False)
		self.axes.spines['right'].set_visible(False)
		self.axes.spines['top'].set_visible(False)
		self.axes.xaxis.set_major_formatter(FormatStrFormatter('%g'))
		# Only show ticks on the left and bottom spines
		self.axes.yaxis.set_ticks_position('none')
		self.axes.yaxis.set_ticks([])
		self.axes.xaxis.set_ticks_position('bottom')
		self.last_y = 0

		self.branch_list = []

	def count_leaves(self,lineage_tree):
		total_leaves= len(lineage_tree.keys())
		for key in lineage_tree.keys():
			total_leaves += self.count_leaves(lineage_tree[key])
		return total_leaves

	def mut_time(self,lineage):
		barcode_epoch = len(lineage.ID.split("_"))

		first_epoch_with_positive_fitness  = 1 + next(i for i in range(0,len(lineage.fitness)) if numpy.sum(abs(lineage.fitness[:i+1]))>5*10**-3)
		if first_epoch_with_positive_fitness == barcode_epoch:
			pos = 11*(barcode_epoch-1) 
		else:
			pos = 11*(barcode_epoch-1) + 5
		
		return 10*pos - 5 * (first_epoch_with_positive_fitness == barcode_epoch)

	def death_time(self,lineage):
		return 1 + next(i for i in reversed(range(0,len(lineage.freqs))) if lineage.freqs[i] > 10**-4)

	def list_ancestors(self,child_ID, tree):
		ancestor = next(item for item in tree.keys() if child_ID.startswith(item))
		if ancestor == child_ID:
			return []
		else:
			return itertools.chain([ancestor],self.list_ancestors(child_ID,tree[ancestor]))

	def add_subtree(self,lineage_tree,dictionary,all_lineages,current_ID = "",parent_y = None):

		root_branches = lineage_tree.keys()
		if len(root_branches) > 0:

			self.last_y += 1

			mut_time = [self.mut_time(dictionary[brch]) for brch in root_branches]
			leaf_counts = [self.count_leaves(lineage_tree[brch]) for brch in root_branches]
			death_times = [self.death_time(dictionary[brch]) for brch in root_branches]
			sorted_root_branches = [x for y,z,w,x in reversed(sorted( zip(mut_time,leaf_counts,death_times,root_branches) ) )]
			
			self.branch_list.append(TreeBranch(dictionary[current_ID],yval= self.last_y + 1,parent_y = parent_y))
			self.last_y += 1	
			
			parent_y = self.last_y 
			for brch in sorted_root_branches:
				if self.count_leaves(lineage_tree[brch]) > 0 or all_lineages or current_ID != "":
					self.add_subtree(lineage_tree[brch],dictionary,all_lineages,current_ID = brch,parent_y = parent_y)

			self.last_y += 1
		else:
			self.branch_list.append(TreeBranch(dictionary[current_ID],parent_y = parent_y, yval= self.last_y + 1))
			self.last_y += 1

	def plot_tree(self,lineage_tree,dictionary,all_lineages = True,labels = False,resize = False):
		self.add_subtree(lineage_tree,dictionary,all_lineages)
		for brch in self.branch_list:
			if brch.x_end >= self.axes.get_xlim()[1]:
				brch.x_end = self.axes.get_xlim()[1] - 40
			self.axes.plot([brch.x_begin,brch.x_end],[brch.y,brch.y],color = 'k',lw = 1)
			self.axes.plot([brch.x_begin,brch.x_begin],[brch.parent_y,brch.y],color = 'k',lw = 1)
			if brch.ID != "":
				list_of_ancestors = self.list_ancestors(brch.ID,lineage_tree)
				for anc in list_of_ancestors:
					self.axes.scatter([brch.x_end],[brch.y],marker = 's', color = dictionary[anc].color, lw = 1.,zorder = 20,edgecolor= 'k')
			self.axes.scatter([brch.x_end],[brch.y],marker = 's', color = dictionary[brch.ID].color, lw = 1.,zorder = 20,edgecolor= 'k')
			if labels:
				if brch.ID == "":
					self.axes.text(brch.x_end + 20, brch.y, 'ancestor',va = 'center', ha = 'left', fontsize = 6)
				else:
					self.axes.text(brch.x_end + 20, brch.y, brch.ID,va = 'center', ha = 'left', fontsize = 6)
		if resize:
			self.fig.set_size_inches(9, 8/44.*len(self.branch_list))
			# self.axes.text(brch.x_end + 5, brch.y, brch.ID,va = 'center', ha = 'left', fontsize = 6)






class BarChart:
	def __init__(self,figure,axes,xvals):
		self.fig = figure
		self.axes = axes
		self.xvals = xvals
		self.bar_width = (xvals[1] - xvals[0])*0.7
		self.tolerance = self.bar_width
		self.bar_lower = numpy.zeros(len(self.xvals))
		self.temp_bar = numpy.zeros(len(self.xvals))
		self.axes.set_xlim([min(xvals),max(xvals)])
		start, end = self.axes.get_xlim()
		self.axes.xaxis.set_ticks([start, 0, end])
		self.axes.xaxis.set_ticklabels([start,0,end])
		self.axes.xaxis.set_ticks(numpy.arange(start,end,0.01),minor = True)
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
			lw = 0.3
		else:
			lw = 0

		axes_zorder = self.axes.patch.get_zorder()
		self.axes.bar(self.xvals[spot],height,width = self.bar_width,bottom = self.bar_lower[spot],color = color ,align = 'center',lw = lw, edgecolor = (0.0,0.0,0.0,1.0),zorder = axes_zorder)

		if not no_error_bar:
			errorbar_y = (self.bar_lower[spot]*2 + height)/2.
			self.add_error_bar(xval,x_lower,x_upper,errorbar_y,color,outline)
		if update:
			self.bar_lower[spot] += height

		if max(self.bar_lower) > 1.0:
			sys.stderr.write('bar exceeded 1.0\n')

	def add_lineage(self, lineage, use_epoch,color = None,no_error_bar = False,outline = False,scale = None,update = True, t = 0):
		if color == None:
			color = lineage.color

		xval = lineage.relative_fitness[0]
		x_lower = lineage.relative_fitness_CI[0][0]
		x_upper = lineage.relative_fitness_CI[0][1]

		height = lineage.freqs[11*(use_epoch) + t]
		if scale is not None:
			height = height/scale[11*(use_epoch + 1) - 1]
			print "adjusting scale"
		self.add_bar(xval,x_lower,x_upper,height,color,no_error_bar,outline, update)

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

	def place(self,lineage,ignore_epochs = 0,color = None, single_epoch = False, outline = False,zorder = 0, scale = None):
		if color is None:
			color = lineage.color

		ignore_epochs = (len(lineage.ID.split("_")) - 1)
		start_timepoint = ignore_epochs*11 
		xvals = self.xvals[start_timepoint:]
		lower = lineage.muller_lower[start_timepoint:]

		if scale is None:	
			upper = lineage.muller_lower[start_timepoint:] + lineage.freqs[start_timepoint:]
		else:
			upper = lineage.muller_lower[start_timepoint:] + lineage.freqs[start_timepoint:]/scale[start_timepoint:]
		if single_epoch:
			xvals = xvals[:12]
			lower = lower[:12]
			upper = upper[:12]
		if outline:
			lw = .2
		else:
			lw = 0.0
		R,G,B,A = color
		self.axes.fill_between(xvals,lower, upper, lw = lw, facecolor = color,zorder=zorder,edgecolor = (R,G,B,1.))


	def save(self,filename,**kwargs):
		sys.stderr.write('saving figure to '+ filename + '\n')
		if 'dpi' in kwargs.keys():
			pass
		else:
			kwargs.update({'dpi':400})
		self.fig.savefig(filename,bbox_inches = 'tight', **kwargs)
