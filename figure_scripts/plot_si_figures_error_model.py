import sys, os, glob, csv, re
import string, math, numpy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col
import matplotlib.patches as patches
import argparse
from scipy.special import erf
import matplotlib.gridspec as gridspec
import scipy.stats

# add custom modules to path
sys.path.insert(0,'../modules/') 

# import custom modules
import config
import local_matplotlibrc
import lineage.inference_params as inference_params
import lineage.file_parser as file_parser
from lineage.fitness_estimator import *


def plot_hist(best_kappa,counts_before, counts_after, total_reads_ratio, epoch = 0, t = 0,plot = True):
	
	bin_min = 1
	bin_max = 200
	bin_num = 100
	dx = 1./(numpy.linspace(bin_min,bin_max,bin_num)[1] - numpy.linspace(bin_min,bin_max,bin_num)[0])
	x_vals = numpy.linspace(bin_min,bin_max,bin_max - bin_min)

	before_filter = numpy.logical_and(counts_before <= inference_params.drifting_lineage_size_range[1],
									  counts_before >= inference_params.drifting_lineage_size_range[0])
	counts_after = counts_after[before_filter]
	counts_before = counts_before[before_filter]

	hist_after , bin_edges= numpy.histogram(counts_after, bins = numpy.linspace(bin_min,bin_max,bin_num))
	bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.

	unique_before = numpy.unique(counts_before)
	weights = numpy.asarray([float(len(counts_before[abs(counts_before - value) < 0.5])) for value in unique_before])

	exponent =numpy.sqrt(numpy.outer(unique_before*total_reads_ratio,numpy.ones(len(x_vals)))) - numpy.sqrt(x_vals)
	dist = numpy.exp(-(exponent**2)/best_kappa)/(x_vals)**0.75 
	dist = dist.T * dx
	dist *= numpy.sqrt(unique_before*total_reads_ratio/4./math.pi/best_kappa)
	dist *= weights
	dist = numpy.sum(dist,axis = 1)
	best_histogram = scipy.stats.binned_statistic(x_vals, dist, statistic='sum', bins=bin_edges, range=None)[0]

	ax_hists[epoch,t].scatter(bin_centres,hist_after, marker = 'o',color = config.pop_colors[population],lw = 0, s= 8,alpha = 0.8)
	ax_hists[epoch,t].plot(bin_centres,best_histogram,lw = 1,color = config.pop_colors[population],alpha = 0.8)

	ax_hists[epoch,t].text(0.2,0.75,r'$t = %d.%d$'%(epoch + 1, t * 10), transform = ax_hists[epoch,t].transAxes,fontsize = 6)

	ax_hists[epoch,t].set_xticks([0,100,200])
	ax_hists[epoch,t].set_xlim([bin_edges[0],bin_edges[-1]])

	if epoch == 3 and t == 0:
		ax_inset.scatter(bin_centres,hist_after, marker = 'o',color = config.pop_colors[population],lw = 0, s= 8)
		ax_inset.plot(bin_centres,best_histogram,lw = 1, color= config.pop_colors[population])


fig,ax = plt.subplots(figsize = (7,2.7))
outer_grid = gridspec.GridSpec(2,2, width_ratios = [5,2],height_ratios = [0.7,2],hspace = 0.1,wspace = 0.1)
ax = plt.subplot(outer_grid[3])
ax_kappas = plt.subplot(outer_grid[2])
ax_no_samples = plt.subplot(outer_grid[0])

fig_inset, ax_inset = plt.subplots(1,1,figsize = (1,1))

fig_hists,ax_hists = plt.subplots(10,11,figsize = (7.7,9))
for population in config.populations:
	max_barcode = config.max_barcode[population]


	timepoints, data, counts = file_parser.get_data(population, config.barcode_data_root_directory, as_matrix = True)
	count_data = [data[i]*counts[11*i:] for i in range(0,len(data))]

	kappa_filename = config.error_model_directory+'%s-kappas.tsv' % (population)

	kappa_file = open(kappa_filename,'r')
	kappa_file.readline()	
	kappas = numpy.asarray([float(item) for item in kappa_file.readline().strip().split('\t')], dtype = float)
	kappa_CI_lower = numpy.asarray([float(item) for item in kappa_file.readline().strip().split('\t')], dtype = float)
	kappa_CI_upper = numpy.asarray([float(item) for item in kappa_file.readline().strip().split('\t')], dtype = float)
	ns = numpy.asarray([float(item) for item in kappa_file.readline().strip().split('\t')], dtype = float)
	kappa_file.close()

	# sort count data on min size in first epoch
	arrinds= [numpy.argsort(numpy.min(count_data[i].T[:11],axis = 0)) for i in range(0,len(data))]
	sorted_mins = [numpy.min(count_data[i].T[:11],axis = 0)[arrinds[i][::-1]] for i in range(0,len(data))]
	sorted_count_data = [count_data[i][arrinds[i][::-1]] for i in range(0,len(data))]
	sorted_count_data = numpy.asarray(sorted_count_data)

	centers_before= numpy.asarray([40])
	for timepoint in xrange(0,inference_params.INTERVALS_PER_EPOCH*max_barcode-1):
		epoch = timepoint / inference_params.INTERVALS_PER_EPOCH
		t = timepoint % inference_params.INTERVALS_PER_EPOCH
		counts_before, counts_after = sorted_count_data[epoch].T[t],sorted_count_data[epoch].T[t+1]
		plot_hist(kappas[timepoint],counts_before,counts_after, total_reads_ratio = 1.*counts[11*epoch+t+1]/counts[11*epoch+t], epoch = epoch, t = t)

	counts = counts[:11*max_barcode]
	kappas = kappas[:11*max_barcode-1]
	kappa_CI_lower = kappa_CI_lower[:11*max_barcode-1]
	kappa_CI_upper = kappa_CI_upper[:11*max_barcode-1]
	timepoints = timepoints[0][:11*max_barcode]

	ratios = counts[1:]*1./counts[:-1]+1.
	next_counts = counts[1:]

	no_ev_intervals = inference_params.EVOLUTION_INTERVALS_PER_EPOCH*max_barcode
	no_bc_intervals = max_barcode - 1
	
	ev_kappas, ev_next_counts, ev_ratios = numpy.zeros(no_ev_intervals), numpy.zeros(no_ev_intervals), numpy.zeros(no_ev_intervals)
	ev_kappas_CI_upper, ev_kappas_CI_lower = numpy.zeros(no_ev_intervals), numpy.zeros(no_ev_intervals)
	ev_timepoints = numpy.zeros(no_ev_intervals)
	ev_ns = numpy.zeros(no_ev_intervals)
	bc_kappas, bc_next_counts, bc_ratios = numpy.zeros(no_bc_intervals), numpy.zeros(no_bc_intervals), numpy.zeros(no_bc_intervals)
	bc_kappas_CI_upper, bc_kappas_CI_lower = numpy.zeros(no_bc_intervals), numpy.zeros(no_bc_intervals)
	bc_timepoints = numpy.zeros(no_bc_intervals)
	bc_ns= numpy.zeros(no_bc_intervals)

	error_bar_kwargs = dict(markersize = 4, capthick=0.5,color = config.pop_colors[population],alpha = 0.8)


	for ep in range(0,max_barcode):
		begin_copy, end_copy = ep*inference_params.EVOLUTION_INTERVALS_PER_EPOCH, (1+ep)*inference_params.EVOLUTION_INTERVALS_PER_EPOCH
		begin, end = ep*inference_params.INTERVALS_PER_EPOCH, ep*inference_params.INTERVALS_PER_EPOCH + inference_params.EVOLUTION_INTERVALS_PER_EPOCH
		
		ev_kappas[begin_copy:end_copy] = kappas[begin:end]
		ev_kappas_CI_upper[begin_copy:end_copy] = kappa_CI_upper[begin:end]
		ev_kappas_CI_lower[begin_copy:end_copy] = kappa_CI_lower[begin:end]
		ev_next_counts[begin_copy:end_copy] = next_counts[begin:end]
		ev_ratios[begin_copy:end_copy] = ratios[begin:end]
		ev_timepoints[begin_copy:end_copy] = timepoints[begin:end]

		ev_ns[begin_copy:end_copy] = ns[begin:end]

		if ep < max_barcode - 1:
			bc_kappas[ep] = kappas[end]
			bc_kappas_CI_upper[ep] = kappa_CI_upper[end]
			bc_kappas_CI_lower[ep] = kappa_CI_lower[end]
			bc_next_counts[ep] = next_counts[end]
			bc_ratios[ep] = ratios[end]
			bc_timepoints[ep] = timepoints[end]
			bc_ns[ep] = ns[end]

	ax_kappas.errorbar(ev_timepoints, ev_kappas, 
						yerr=[ev_kappas - ev_kappas_CI_lower+10**-6, ev_kappas_CI_upper - ev_kappas],
						fmt = '.',**error_bar_kwargs)
	ax_kappas.errorbar(bc_timepoints, bc_kappas, 
						yerr=[bc_kappas - bc_kappas_CI_lower+10**-6, bc_kappas_CI_upper - bc_kappas],
						fmt = 'o',**error_bar_kwargs)

	ax_no_samples.scatter(ev_timepoints,ev_ns,marker = '.', color = config.pop_colors[population],s = 16,lw = 0, alpha = 0.8)
	ax_no_samples.scatter(bc_timepoints,bc_ns,marker = 'o', color = config.pop_colors[population],s = 16,lw = 0, alpha = 0.8)

	xs_ev = numpy.asarray([ev_next_counts,ev_ratios])
	xs_ev = xs_ev.T
	weights_ev = (4./(ev_kappas_CI_upper - ev_kappas_CI_lower))**2.

	xsw_ev = xs_ev * numpy.sqrt(weights_ev[:,numpy.newaxis])
	ysw_ev = ev_kappas * numpy.sqrt(weights_ev)

	xs_bc = numpy.asarray([bc_next_counts,bc_ratios])
	xs_bc = xs_bc.T
	weights_bc = (4./(bc_kappas_CI_upper - bc_kappas_CI_lower))**2.

	xsw_bc = xs_bc * numpy.sqrt(weights_bc[:,numpy.newaxis])
	ysw_bc = bc_kappas * numpy.sqrt(weights_bc)

	# bover_Nb_ev, beta_ev = numpy.linalg.lstsq(xsw_ev,ysw_ev,rcond=None)[0]
	bover_Nb_ev, beta_ev = np.linalg.solve(xsw_ev.T.dot(xsw_ev), xsw_ev.T.dot(ysw_ev))
	bover_Nb_ev_error, beta_ev_error = np.sqrt(np.linalg.inv(xsw_ev.T.dot(xsw_ev)).diagonal())

	print population, "evolution: [%.3g (%.3g,%.3g)," % (2*bover_Nb_ev, 
														(2*bover_Nb_ev - 4*bover_Nb_ev_error),
														(2*bover_Nb_ev + 4*bover_Nb_ev_error))
	print "%.3g (%.3g,%.3g)]" %  (2*beta_ev, 
								 (2*beta_ev - 4*beta_ev_error), 
								 (2*beta_ev + 4*beta_ev_error))

	# print "  ", "    error: (%.3g,%.3g)" % (2*bover_Nb_ev_error, 2*beta_ev_error)

	# bover_Nb_bc, beta_bc = numpy.linalg.lstsq(xsw_bc,ysw_bc,rcond = None)[0]
	bover_Nb_bc, beta_bc = np.linalg.solve(xsw_bc.T.dot(xsw_bc), xsw_bc.T.dot(ysw_bc))
	bover_Nb_bc_error, beta_bc_error = np.sqrt(np.linalg.inv(xsw_bc.T.dot(xsw_bc)).diagonal())

	# print population, "barcoding: (%.3g,%.3g)" % (2*bover_Nb_bc, 2*beta_bc)
	# print "  ", "    error: (%.3g,%.3g)" % (2*bover_Nb_bc_error, 2*beta_bc_error)
	
	print population, "barcoding: [%.3g (%.3g,%.3g)," % (2*bover_Nb_bc, 
														(2*bover_Nb_bc - 4*bover_Nb_bc_error),
														(2*bover_Nb_bc + 4*bover_Nb_bc_error))
	print "%.3g (%.3g,%.3g)]" %  (2*beta_bc, 
								 (2*beta_bc - 4*beta_bc_error),
								 (2*beta_bc + 4*beta_bc_error))

	kappa_pred_ev = bover_Nb_ev*ev_next_counts + beta_ev*ev_ratios
	kappa_pred_bc = bover_Nb_bc*bc_next_counts + beta_bc*bc_ratios

	ax.errorbar(kappa_pred_ev, (ev_kappas), 
				yerr=[ev_kappas - ev_kappas_CI_lower+10**-6, ev_kappas_CI_upper - ev_kappas],
				fmt = '.',**error_bar_kwargs)
	ax.errorbar(kappa_pred_bc, (bc_kappas), 
				yerr=[bc_kappas - bc_kappas_CI_lower+10**-6, bc_kappas_CI_upper - bc_kappas],
				fmt = 'o',**error_bar_kwargs)

	ax.plot(numpy.linspace(1,100,1000),numpy.linspace(1,100,1000),lw = 1,color = '0.5')


for population in config.populations:
	ax.plot([],[],lw = 5, color = config.pop_colors[population],label = config.pop_labels[population])
ax.scatter([],[], marker = '.', s = 16,color = 'k', alpha = 0.8, label = 'evolution',lw = 0)
ax.scatter([],[], marker = 'o', s = 16, color = 'k',alpha = 0.8, label = 'barcoding',lw = 0)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'$\kappa_{t,\mathrm{model}} = \hat{\beta}(1 + \frac{R_{t + \Delta t}}{R_{t}}) +  \frac{\hat{b}}{N_b}  R_{t + \Delta t}$',fontsize = 8)
# ax.set_ylabel(r'Best fit $\kappa_t$',fontsize = 8)

ax.legend(bbox_to_anchor = (0.75,1.00),scatterpoints = 1,handlelength = 0.5,fontsize = 7)
ax.set_ylim(1,200)
ax.set_ylim(1,200)


ax_kappas.set_ylim(1,200)
ax_kappas.set_yscale('log')
ax_no_samples.set_yscale('log')
xticks = []
xticklabels = []
xtick_dt = 110
for i in range(0,max_barcode+1):
	if i == 0:
		xticks.extend([xtick_dt*i,xtick_dt*i+100])
	else:
		xticks.extend([xtick_dt*i,xtick_dt*i+100])
	xticklabels.extend(['%d.0'%(i+1),'%d.100'%(i+1)])
ax_kappas.set_xticks(xticks)
ax_kappas.set_xticklabels(xticklabels)
ax_kappas.set_xlim([0,xtick_dt*max_barcode-10])

i = 0
for tick in ax_kappas.get_xticklabels():
	tick.set_rotation(90)
	if i % 2 == 0:
		tick.set_horizontalalignment('left')
	else:
		tick.set_horizontalalignment('right')
	i+=1
ax_no_samples.set_xticks(xticks)
ax_no_samples.set_xticklabels([])
ax_no_samples.set_xlim([0,xtick_dt*max_barcode-10])

ax_kappas.set_xlabel(r'Generations, $t$',fontsize = 8)
ax_kappas.set_ylabel(r'Best fit $\kappa_t$',fontsize =8)
ax_no_samples.set_ylabel('No. lineages\nused in fit', fontsize = 8)

fig.savefig(config.figure_directory + 'si/kappa_errors.pdf',bbox_inches = 'tight')

for i in range(0,len(ax_hists)):
	for j in range(0,len(ax_hists[i])):
		ax_hists[i,j].patch.set_alpha(0.) 
		ax_hists[i,j].set_ylim([0,ax_hists[i,j].get_ylim()[1]])
		ax_hists[i,j].tick_params(axis='both', which='major', labelsize=5,length = 1,pad = 2)
		for axis_position in ['top','right']:
			ax_hists[i,j].axes.spines[axis_position].set_visible(False)
		ax_hists[i,j].axes.xaxis.set_ticks_position('bottom')
		ax_hists[i,j].axes.yaxis.set_ticks_position('left')

		if i == len(ax_hists) - 1 and j == len(ax_hists[i]) - 1:
			for axis_position in ['bottom','left']:
				ax_hists[i,j].axes.spines[axis_position].set_visible(False)
			ax_hists[i,j].axes.xaxis.set_ticks_position('none')
			ax_hists[i,j].axes.yaxis.set_ticks_position('none')
			ax_hists[i,j].axes.xaxis.set_ticks([])
			ax_hists[i,j].axes.yaxis.set_ticks([])

		if (i == len(ax_hists) - 1 and j != len(ax_hists[i]) - 1) or i == len(ax_hists) - 2 and j == len(ax_hists[i]) - 1:
			ax_hists[i,j].set_xlabel(r'Reads, $r_{i,T + \Delta T}$',fontsize = 5,labelpad = 0)
		if j == 0:
			ax_hists[i,j].set_ylabel('Number',fontsize = 5,labelpad =0)

ax_hists[0,10].add_patch(patches.Rectangle(
		    	(-0.3, 1.1), 1.6, -11.25,
		    	transform = ax_hists[0,10].transAxes,
		    	facecolor= 'k',
		    	alpha = 0.1,
		    	edgecolor = "none",
		    	clip_on = False))
ax_hists[9,10].text(0.5,0.8,'Barcoding', transform = ax_hists[9,10].transAxes,fontsize = 7,ha = 'center')


OUTFILE = config.figure_directory + 'si/kappa_fit.pdf'
fig_hists.savefig(OUTFILE,bbox_inches = 'tight')

	
		

