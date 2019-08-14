echo "Generating left subpanel of Fig 1d, and equivalent for YPA population..."
python plot_muller_raw.py --stretched # generates subpanel of Figure 1d and equivalent for YPA population
echo "Generating right subpanel of Fig 1d, and equivalent for YPA population..."
python plot_muller_fitness.py --stretched # generates subpanel of Figure 1d and equivalent for YPA population
echo "Generating Figure 2..."
python plot_figure2.py --stretched
echo "Generating Fig 3 using average fitness..."
python plot_figure3.py
echo "Generating Fig 4, Extended Data Figs 3 and 4, and Fig S14..."
python plot_figure4_modified_layout.py # generates Figure 4, Extended Data Figures 3 and 4
echo "Generating Extended Data Fig 1..."
python plot_ed_figure_wgs_trajectories.py
echo "Generating Extended Data Fig 2..."
python plot_ed_figure_mean_fitness.py
echo "Generating Fig S4..."
python plot_si_figure_count_distribution_after_barcoding.py
echo "Generating Figs S5 and S6..."
python plot_si_figures_error_model.py
echo "Generating Fig S7..."
python plot_si_figure_time_to_detect.py
echo "Generating Fig S8..."
python plot_si_figure_jdfe.py
echo "Generating Fig S9..."
python plot_si_figure_log_frequency_change.py
echo "Generating Fig S10..."
python plot_si_figure_fitness_in_two_environments.py
echo "Generating Fig S11..."
python plot_figure3.py --evolution
echo "Generating Fig S12..."
python plot_figure3.py --barcoding
echo "Generating Fig S13..."
python plot_mullers_replicates.py
echo "Generating Figs S15 and S16..."
python plot_si_figure_segment_pvalues.py
echo "Generating Fig S17..."
python plot_si_figure_clone_pvalues.py
echo "Generating Fig S18..."
python plot_simulated_wgs_results.py

echo "Generating Tables S3 and S4..."
python print_clone_table.py
