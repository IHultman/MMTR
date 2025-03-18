<pre>
After running the 'construct_mmtr_sims.R' and 'mmmtr_sim_rdata_to_mat.R' 
scripts, you should have a directory with multiple sets of simulated data, 
structured as follows:

  ./simulated_mmtr_matlab_data
    |
    |_sim_set_1
    | |
    | |_simulation_caseID_1_group_sizeID_1_obs_sizeID_1.mat
    | |_simulation_caseID_1_group_sizeID_1_obs_sizeID_2.mat
    | ...
    | |_simulation_caseID_2_group_sizeID_4_obs_sizeID_2.mat
    |
    |_sim_set_2
    | |
    | |_simulation_caseID_1_group_sizeID_1_obs_sizeID_1.mat
    | ...
    | |_simulation_caseID_2_group_sizeID_4_obs_sizeID_2.mat
    |
    ...
    |_sim_set_100
      |
      |_simulation_caseID_1_group_sizeID_1_obs_sizeID_1.mat
      ...
      |_simulation_caseID_2_group_sizeID_4_obs_sizeID_2.mat

In order to use the Kruskal library to estimate model parameters with max 
rank 10 on the simulated data in sim_set_1 and save the results to some 
dirctory "./results_save_dir/", call 'run_kruskal_for_mmtr_sim_analysis.m' 
from the command line as follows:

  > matlab -nodisplay -nodesktop -nosplash -r "run_kruskal_for_mmtr_sim_analysis ./simulated_mmtr_matlab_data ./results_save_dir/ 1 10"

For sim_set_2:

  > matlab -nodisplay -nodesktop -nosplash -r "run_kruskal_for_mmtr_sim_analysis ./simulated_mmtr_matlab_data ./results_save_dir/ 2 10"

so on and so forth.

Similarly, after running the 'construct_equicorr_sims.R' and 
'equicorr_sim_rdata_to_mat.R' scripts, you should have a directory with 
multiple sets of simulated data, structured as follows:

  ./simulated_equicorr_matlab_data
    |
    |_sim_set_1
    | |
    | |_simulation_caseID_1_group_sizeID_1_obs_sizeID_1.mat
    | |_simulation_caseID_1_group_sizeID_2_obs_sizeID_1.mat
    | ...
    | |_simulation_caseID_2_group_sizeID_4_obs_sizeID_1.mat
    |
    |_sim_set_2
    | |
    | |_simulation_caseID_1_group_sizeID_1_obs_sizeID_1.mat
    | ...
    | |_simulation_caseID_2_group_sizeID_4_obs_sizeID_1.mat
    |
    ...
    |_sim_set_100
      |
      |_simulation_caseID_1_group_sizeID_1_obs_sizeID_1.mat
      ...
      |_simulation_caseID_2_group_sizeID_4_obs_sizeID_1.mat

In order to use the Kruskal library to estimate model parameters with max 
rank 10 on the simulated data in sim_set_1 and save the results to some 
dirctory "./results_save_dir/", call 'run_kruskal_for_equicorr_sim_analysis.m' 
from the command line as follows:

  > matlab -nodisplay -nodesktop -nosplash -r "run_kruskal_for_equicorr_sim_analysis ./simulated_equicorr_matlab_data ./results_save_dir/ 1 10"

For sim_set_2:

  > matlab -nodisplay -nodesktop -nosplash -r "run_kruskal_for_equicorr_sim_analysis ./simulated_equicorr_matlab_data ./results_save_dir/ 2 10"

so on and so forth.
</pre>
