#pragma once
#include <iostream>

#include "functions_library.cuh"
#include "parameter_load.h"
#include "node_within_host.cuh"
#include <cstdlib>

#include <curand.h>
#include <curand_kernel.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>

#include <sstream>

#include <algorithm>
#include <random>
#include <chrono>
#include <iomanip>
#include <string>
#include <map>

#include <string>
#include <vector>
#include <queue>

#include <thread>
#include <mutex>
#include <shared_mutex>

using namespace std;

class simulator_Master
{
private:
    shared_mutex g_mutex;
    mutex gi_mutex;

    string output_Folder_location;
    string intermediate_Folder_location;

    string intermediary_Sequence_location;
    string intermediary_Index_location;

    int CUDA_device_number;

    int *CUDA_device_IDs;
    int num_Cuda_devices;
    int gpu_Limit;
    int *tot_Blocks;
    int *tot_ThreadsperBlock;

    int CPU_cores;
    int max_Cells_at_a_time = 0;
    string output_Folder;
    string output_Folder_Sequences;
    string intermediate_Folders;
    string intermediate_sequence_Store;

    string enable_Folder_management = "NO";
    string enable_Compression = "NO";

    int max_sequences_per_File = 0;

    string node_Master_location;
    string sequence_Master_location;

    string multi_Read;

    string start_Date;

    string network_Model = "NA";
    int Total_number_of_Nodes = 0;

    // vector<vector<pair<int, int>>> each_Nodes_Connection;
    // TODO: change
    vector<vector<int>> each_Nodes_Connection_INT;
    vector<pair<int, int>> all_node_IDs;

    // vector<int> search_Indexes;
    // int overall_Found = 0;

    string connection_Model = "FIXED";

    int BA_FIXED = 1;

    int BA_NB_sucesses = 0;
    float BA_NB_probability = 0;

    float BA_Poisson_mean = 0;

    int SCM_number_of_caves = 0;
    int SCM_number_of_nodes_per_cave = 0;
    // int Total_num_Nodes_SCM = 0;

    int DCM_number_of_caves = 0;
    int *per_cave_Stride;

    int DC_ND_sucesses = 0;
    float DC_ND_probability = 0;

    float DC_Poisson_mean = 0;

    float DC_percent_Neighbouring = 0;
    float DC_percent_Global_freedom = 0;

    float ER_link_probability = 0;

    float generation_Time = 0;

    float shape_days_in_Host = 0;
    float scale_days_in_Host = 0;

    int trials_Sampling = -1;
    string sampled_sequences_Folder;
    // 1 equals active
    int resampling = -1;
    set<int> sampled_Nodes;
    float sampling_trials = 0;
    float sampling_probability = 0;

    // col 0 = 0= FIXED, 1 = Binomial
    // co 1 = fixed/ trials
    // col 2 = prob
    float *per_Node_sampling;

    int limit_Sampled = -1;
    int count_Sampling_instances = 0;

    string progeny_distribution_Model = "NA";

    float *progeny_distribution_parameters_Array;

    // int progeny_NB_sucesses = 0;
    // float progeny_NB_probability = 0;

    // float progeny_Poisson_mean = 0;

    // float progeny_Gamma_shape = 0;
    // float progeny_Gamma_scale = 0;

    int number_of_node_Profiles = 0;
    string node_Profile_folder_Location = "";
    float *node_profile_Distributions;
    vector<string> profile_names;

    float **each_Node_Profile_Configuration;

    // 0 = No Change
    // 1 = Removed
    //-1 = Less infectious col 1 = alpha, col 2 = beta
    float **node_sampling_effect;

    int num_tissues_per_Node = 0;
    vector<string> tissue_Names;

    int entry_tissues = 0;
    int *entry_array;

    int infectious_tissues = 0;
    int *infectious_array;

    int terminal_tissues = 0;
    int *terminal_array;

    int exit_tissues = 0;
    int *exit_array;

    int sampling_tissues = 0;
    int *sampling_array;

    string first_Infection = "Random";

    // rows = profiles
    // columns
    // 0 = distribution_type; 0=Binomial, -1 = Fixed
    // 1 = trials/ fixed value
    // 2 = prob
    float **infectious_load_Profiles_param;
    float **terminal_load_Profiles_param;

    string viral_Migration = "No";
    float **viral_Migration_Values;
    int *migration_start_Generation;

    // column 0 = Yes = 1, NO =0;
    // column 1 = trials
    // column 2 = prob
    float **profile_tissue_Limits;

    // rows = profiles
    //  columns = tissues;
    // 0 = Binomial
    // 1 = Gamma
    int **cell_Distribution_Type;
    vector<vector<pair<float, float>>> ALL_profiles_Tissue_cell_disribution;

    int *num_replication_phases;
    int *tissue_param_profile_Stride;
    float **tissue_replication_data;

    string parent_Sequence_Folder;
    // row 0 = mutations
    // row 1 = recombination
    // row 2 = proof reading
    // 0 = inactive 1 = activated
    int *mutation_recombination_proof_Reading_availability;

    // 0 = no reinfection
    // 1 = reinfection YES
    int reinfection_Availability = 0;

    // 0 = Fixed (0) or Binomial (1)
    // 1 = Fixed value / Trials bim,
    // 2 = Probability binom
    float *transmission_parameters;

    // row = profile
    //  0 = Fixed (0) or Binomial (1)
    //  1 = Fixed value / Trials bim,
    //  2 = Probability binom
    float **infection_parameters;

    float *Reference_fitness_survivability_proof_reading;

    // columns
    // 0 = position-1;
    // 1 = A
    // 2 = T
    // 3 = G
    // 4 = C
    float **sequence_Fitness_changes;
    float **sequence_Survivability_changes;
    float **sequence_Proof_reading_changes;

    // 0 = fitness
    // 1 = survivability
    // 2 = proof reading
    int *num_effect_Segregating_sites;

    int mutation_Hotspots = 0;
    int recombination_Hotspots = 0;

    float **A_0_mutation;
    float **T_1_mutation;
    float **G_2_mutation;
    float **C_3_mutation;

    float **mutation_hotspot_parameters;
    float **recombination_hotspot_parameters;

    int *tot_prob_selectivity;

    int *recombination_prob_Stride;
    int *recombination_select_Stride;

    float **recombination_Prob_matrix;
    float **recombination_Select_matrix;

    string output_Network_location = "";
    string network_File_location = "";
    string Host_source_target_network_location = "";
    string output_Node_location = "";

    // vector<string> raw_parent_Sequences;
    int genome_Length = 0;

    string infected_to_Recovered = "NO";

    string stop_after_generations = "NO";
    // 0 = Generations 1 = Date
    int stop_gen_Mode = 0;
    int stop_generations_Count = 0;
    float stop_Date = 0;

public:
    simulator_Master(string parameter_Master_Location);

    void configure_Network_Profile(string network_Profile_File, parameter_load &Parameters);

    void ingress();

    void network_Manager(functions_library &functions);

    void BA_Model_Engine(functions_library &functions);
    void SCM_Model_Engine(functions_library &functions);
    void DCM_Model_Engine(functions_library &functions);
    void RANDOM_Model_Engine(functions_library &functions);
    void ER_RANDOM_Model_Engine(functions_library &functions);

    void node_Master_Manager(functions_library &functions);

    void sequence_Master_Manager(functions_library &functions);

    vector<node_within_host> node_Profile_assignment_Manager(functions_library &functions);
    void node_Profile_assignment_thread(int start_Node, int stop_Node);

    void apollo(functions_library &functions, vector<node_within_host> &Hosts);

    int get_first_Infected(vector<int> &susceptible_Population,
                           vector<int> &infected_Population, functions_library &functions);

    vector<string> read_Reference_Sequences(vector<int> &tissue_Sequence_Count);
    vector<string> read_Reference_Sequence_Files(vector<string> &reference_Files);

    // void Node_search(vector<pair<int, int>> &host_Connections);
    // void thread_Node_search(int start_Node, int stop_Node, vector<pair<int, int>> host_Connections);

    vector<int> get_new_Hosts_Indexes(int &node_Profile, mt19937 &gen, vector<int> &possible_Infections);

    void thread_compress_Folders(int new_dead_Host, vector<node_within_host> &Hosts);
};