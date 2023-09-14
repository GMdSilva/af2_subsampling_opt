from angra.workflow import *
from angra.utilities.utilities import load_from_pickle


prefix = config.SYSTEM_NAME
build_msa(prefix)
optimizer_manager = OptimizerManager(config.SYSTEM_NAME)
variations = optimizer_manager.get_variation_regions()
optimizer_manager.get_parameter_set_variations()
optimizer_manager.plot_variation_regions()
best = optimizer_manager.get_best_parameter_set()
optimizer_manager.plot_parameter_set()
tester = MutantTester(config.SYSTEM_NAME)
tester.test_mutants()
tester.plot_results()
#all_trials = ['2_4', '4_8', '16_32', '32_64', '64_128']
all_trials = ["16_32", "32_64", '64_128', '128_256', "256_512", "2048_4096", "1024_2048"]
get_representative_structures(config.SYSTEM_NAME, all_trials)
mut_manager = MutationClusterManager(config.SYSTEM_NAME)
mut_manager.build_wt_model()
mut_manager.measure_effects()
mut_manager.plot_mutant_state_finder_results()
mut_manager.plot_cluster_results()
mut_manager.plot_clustering_diffs()