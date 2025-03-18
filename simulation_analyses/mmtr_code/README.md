<pre>
After running the 'construct_mmtr_sims.R' script, you should have a 
directory with multiple sets of simulated data, structured as follows:

  ./simulated_mmtr_data
    |
    |_sim_set_1
    | |
    | |_simulation_caseID_1_group_sizeID_1_obs_sizeID_1.RData
    | |_simulation_caseID_1_group_sizeID_1_obs_sizeID_2.RData
    | ...
    | |_simulation_caseID_2_group_sizeID_4_obs_sizeID_2.RData
    |
    |_sim_set_2
    | |
    | |_simulation_caseID_1_group_sizeID_1_obs_sizeID_1.RData
    | ...
    | |_simulation_caseID_2_group_sizeID_4_obs_sizeID_2.RData
    |
    ...
    |_sim_set_100
      |
      |_simulation_caseID_1_group_sizeID_1_obs_sizeID_1.RData
      ...
      |_simulation_caseID_2_group_sizeID_4_obs_sizeID_2.RData

In order to use the mmtr library to estimate model parameters on the
simulated data in sim_set_1 and save the results to some dirctory
"./results_save_dir/", call 'run_mmtr_for_mmtr_sim_analysis.R' from
the command line as follows:

  > Rscript run_mmtr_for_mmtr_sim_analysis.R "./simulated_mmtr_data" 1 "./results_save_dir/"

For sim_set_2:

  > Rscript run_mmtr_for_mmtr_sim_analysis.R "./simulated_mmtr_data" 2 "./results_save_dir/"

so on and so forth.

Similarly, after running the 'construct_equicorr_sims.R' script, you should have a 
directory with multiple sets of simulated data, structured as follows:

  ./simulated_equicorr_data
    |
    |_sim_set_1
    | |
    | |_simulation_caseID_1_group_sizeID_1_obs_sizeID_1.RData
    | |_simulation_caseID_1_group_sizeID_2_obs_sizeID_1.RData
    | ...
    | |_simulation_caseID_2_group_sizeID_4_obs_sizeID_1.RData
    |
    |_sim_set_2
    | |
    | |_simulation_caseID_1_group_sizeID_1_obs_sizeID_1.RData
    | ...
    | |_simulation_caseID_2_group_sizeID_4_obs_sizeID_1.RData
    |
    ...
    |_sim_set_100
      |
      |_simulation_caseID_1_group_sizeID_1_obs_sizeID_1.RData
      ...
      |_simulation_caseID_2_group_sizeID_4_obs_sizeID_1.RData

In order to use the mmtr library to estimate model parameters on the
simulated data in sim_set_1 and save the results to some dirctory
"./results_save_dir/", call 'run_mmtr_for_equicorr_sim_analysis.R' from
the command line as follows:

  > Rscript run_mmtr_for_equicorr_sim_analysis.R "./simulated_equicorr_data" 1 "./results_save_dir/"

For sim_set_2:

  > Rscript run_mmtr_for_equicorr_sim_analysis.R "./simulated_equicorr_data" 2 "./results_save_dir/"

so on and so forth.
</pre>
