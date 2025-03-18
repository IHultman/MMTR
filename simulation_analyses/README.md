<pre>
Use the 'construct_mmtr_sims.R' and 'construct_equicorr_sims.R' scripts to create datasets 
simulated according to the MMTR model and the equicorrelation model respectively.

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
</pre>
