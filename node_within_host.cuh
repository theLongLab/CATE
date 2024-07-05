#pragma once
#include <iostream>

#include "functions_library.cuh"
#include "parameter_load.h"
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

class node_within_host
{
private:
    // shared_mutex g_mutex;
    //  mutex gi_mutex;

    int host_Index, cave_ID, host_ID, profile_ID;
    int num_Generation;

    int infectious_Load;
    int terminal_Load;

    float infection_probability = 1;
    float sampling_Effect;
    // 0=unsampled
    // 1= sampled
    // int sampled_status = 0;

    int num_Tissues;
    int *cell_Limit;

    string status = "Susceptible";

    int current_Generation = -1;
    int *current_Viral_load_per_Tissue;

    int *dead_Particle_count;

    int *parents_Prev_generation;

    // Remember to clear after getting the indexes in th current tissue;
    vector<pair<string, string>> to_write_Sequence_Store;
    vector<string> converted_Sequences;

public:
    node_within_host();

    void setHost(int host_Index, int cave_ID, int host_ID, int profile_ID, int num_Tissues);
    void setNum_Generation(int num_Generation);
    void setInfectious_Load(int infectious_Load);
    void setTerminal_Load(int terminal_Load);
    void setSampling_Effect(float sampling_Effect);
    void setCell_Limit(vector<int> cell_Limit_vec);

    int get_host_Index();
    string get_Name();
    string get_Status();
    int get_Profile();
    int get_Generation();
    int *get_current_Viral_load_per_Tissue();
    float get_infection_probability();

    void set_Infected();
    void set_Infectious();
    void set_Removed();
    void set_Dead();
    void set_Susceptible();

    int get_Load(int &num_tissues_Calc, int *tissue_array);

    int infectious_status(int &num_tissues_Calc, int *tissue_array);
    int terminal_status(int &num_tissues_Calc, int *tissue_array, string source_sequence_Data_folder,
                        string enable_Folder_management, string enable_Compression);

    int terminal_status(int &num_tissues_Calc, int *tissue_array);

    void print_All();

    void begin_Infection(functions_library &functions, string &intermediary_Sequence_location,
                         int entry_tissues, int *entry_array, int &max_sequences_per_File,
                         string &output_Node_location,
                         vector<string> &tissue_Names, string first_Infection);

    vector<set<int>> removed_by_Transfer_Indexes;

    string transfer_Infection(functions_library &functions, string &intermediary_Sequence_location, string &source_Target_file_Location,
                              int &source_Index, int &source_Generation, string &source_Name, int *source_current_Viral_load_per_Tissue,
                              int num_viruses_to_transfer,
                              int &entry_tissues, int *entry_array, int exit_Load, int &exit_tissues, int *exit_array,
                              vector<set<int>> &source_removed_by_Transfer_Indexes,
                              int &max_sequences_per_File,
                              vector<vector<pair<int, int>>> &indexed_Source_Folders,
                              float &decimal_Date,
                              string &Host_source_target_network_location,
                              string &output_Node_location,
                              vector<string> &tissue_Names,
                              mt19937 &gen);

    void intialize_Tissues(string &host_Folder, vector<vector<string>> &tissue_Sequences, functions_library &functions);

    void run_Generation(functions_library &functions, string &multi_Read, int &max_Cells_at_a_time, int &gpu_Limit, int *CUDA_device_IDs, int &num_Cuda_devices, int &genome_Length,
                        int &CPU_cores, int &max_sequences_per_File,
                        string source_sequence_Data_folder, string &output_Node_location,
                        vector<string> &tissue_Names,
                        int *num_replication_phases, float **tissue_replication_data, int *tissue_param_profile_Stride,
                        int terminal_tissues, int *terminal_array,
                        int **cell_Distribution_Type, vector<pair<float, float>> &viral_distribution_per_Tissue_param,
                        float *Reference_fitness_survivability_proof_reading,
                        int *mutation_recombination_proof_Reading_availability,
                        int *num_effect_Segregating_sites,
                        float **sequence_Fitness_changes,
                        float **sequence_Survivability_changes,
                        float **sequence_Proof_reading_changes,
                        int &mutation_Hotspots,
                        float **mutation_hotspot_parameters,
                        float **A_0_mutation,
                        float **T_1_mutation,
                        float **G_2_mutation,
                        float **C_3_mutation,
                        int &recombination_Hotspots,
                        float **recombination_hotspot_parameters,
                        int *tot_prob_selectivity,
                        int *recombination_prob_Stride,
                        int *recombination_select_Stride,
                        float **recombination_Prob_matrix,
                        float **recombination_Select_matrix,
                        float *progeny_distribution_parameters_Array,
                        string &viral_Migration,
                        float **viral_Migration_Values,
                        int *migration_start_Generation,
                        int &overall_Generations,
                        string &infected_to_Recovered,
                        string enable_Folder_management,
                        string enable_Compression,
                        mt19937 &gen);

    int get_generation_Phase(int generation, int *num_replication_phases, float **tissue_replication_data, int *tissue_param_profile_Stride, int &tissue,
                             float &variable_1, float &variable_2);

    vector<int> assign_Cells(functions_library &functions, int **parents_in_Tissue, int num_Viral_particles, int &tissue,
                             int distribution_Type, float &parameter_1, float &parameter_2,
                             set<int> &check_to_Remove,
                             int &gen_Phase, float &variable_1, float &variable_2,
                             mt19937 &gen);

    void simulate_Cell_replication(functions_library &functions, string &multi_Read, int &gpu_Limit, int *CUDA_device_IDs, int &num_Cuda_devices, string &source_sequence_Data_folder, vector<pair<int, int>> &indexed_Tissue_Folder,
                                   int &CPU_cores, int &max_sequences_per_File,
                                   int &genome_Length,
                                   int &tissue, int **parents_in_Tissue, string &tissue_Name,
                                   vector<int> &start_Stop_cells, int &start_Cell, int &stop_Cell, int &num_Cells,
                                   int &Total_seqeunces_to_Process,
                                   int &sequence_Count,
                                   float *Reference_fitness_survivability_proof_reading,
                                   int *mutation_recombination_proof_Reading_availability,
                                   int *num_effect_Segregating_sites,
                                   float **sequence_Fitness_changes,
                                   float **sequence_Survivability_changes,
                                   float **sequence_Proof_reading_changes,
                                   int &mutation_Hotspots,
                                   float **mutation_hotspot_parameters,
                                   float **A_0_mutation,
                                   float **T_1_mutation,
                                   float **G_2_mutation,
                                   float **C_3_mutation,
                                   int &recombination_Hotspots,
                                   float **recombination_hotspot_parameters,
                                   int *tot_prob_selectivity,
                                   int *recombination_prob_Stride,
                                   int *recombination_select_Stride,
                                   float **recombination_Prob_matrix,
                                   float **recombination_Select_matrix,
                                   float *progeny_distribution_parameters_Array,
                                   string &cells_of_parents_location,
                                   string &dead_List, string &sequence_Profiles, string &sequence_parent_Progeny_relationships, string &cells_of_progeny_location,
                                   int &index_Last_Written,
                                   mt19937 &gen);

    void process_Sequences_get_Configuration(functions_library &functions,
                                             vector<string> &collected_Sequences, int *CUDA_device_IDs, int &num_Cuda_devices, int &genome_Length,
                                             float *Reference_fitness_survivability_proof_reading,
                                             int *mutation_recombination_proof_Reading_availability,
                                             int *num_effect_Segregating_sites,
                                             float **sequence_Fitness_changes,
                                             float **sequence_Proof_reading_changes,
                                             int recombination_Hotspots,
                                             float **recombination_hotspot_parameters,
                                             int *tot_prob_selectivity,
                                             int *recombination_prob_Stride,
                                             int *recombination_select_Stride,
                                             float **recombination_Prob_matrix,
                                             float **recombination_Select_matrix,
                                             float *progeny_distribution_parameters_Array,
                                             int **parent_Sequences,
                                             float **sequence_Configuration_standard,
                                             int &start_Index);

    void progeny_Configurator(functions_library &functions,
                              float **cuda_sequence_Configuration_standard, int recombination_Hotspots,
                              int start_Index, int num_Parents_to_Process,
                              int **progeny_Configuration, int *cuda_progeny_Stride, int progeny_Total, int remove_Back);

    void progeny_Completion(functions_library &functions,
                            int *CUDA_device_IDs, int &num_Cuda_devices,
                            int &CPU_cores, int &max_sequences_per_File,
                            int genome_Length, float *Reference_fitness_survivability_proof_reading, int *mutation_recombination_proof_Reading_availability,
                            int *num_effect_Segregating_sites,
                            float **sequence_Survivability_changes,
                            int recombination_Hotspots,
                            float **recombination_hotspot_parameters,
                            int *tot_prob_selectivity,
                            int mutation_Hotspots,
                            float **A_0_mutation,
                            float **T_1_mutation,
                            float **G_2_mutation,
                            float **C_3_mutation,
                            float **mutation_hotspot_parameters,
                            int **parent_Sequences, int num_Parent_sequence, float **sequence_Configuration_standard, int **parent_IDs, int num_Cells, int *cell_Index, int &start_Cell,
                            int **progeny_Configuration, int num_Progeny_being_Processed,
                            float **totals_Progeny_Selectivity,
                            int start_Progeny, int stop_Progeny, int &sequence_Count, string source_sequence_Data_folder,
                            int &tissue, string &tissue_Name,
                            string &dead_List, string &sequence_Profiles, string &sequence_parent_Progeny_relationships, string &cells_of_progeny_location,
                            int &index_Last_Written);

    void write_Full_Sequences_Progeny(functions_library &functions,
                                      int &CPU_cores,
                                      string &tissue_Name,
                                      int &num_Progeny_being_Processed,
                                      int &genome_Length, int &recombination_Hotspots,
                                      int &sequence_Count,
                                      int **parent_IDs,
                                      int **progeny_Sequences, int *Dead_or_Alive, int **progeny_Configuration_Filled,
                                      string intermediary_Tissue_folder,
                                      string &dead_List, string &sequence_Profiles, string &sequence_parent_Progeny_relationships,
                                      string &cells_of_progeny_location, int &start_Cell,
                                      int &max_sequences_per_File, int &last_index_Seq_Written,
                                      int &tissue);

    void thread_Sequence_to_String(int start, int stop, int **progeny_Sequences, int genome_Length, shared_mutex &g_mutex);

    void write_Partial_Sequence_Progeny(functions_library &functions,
                                        string &intermediary_Tissue_folder,
                                        string &dead_List,
                                        int &last_index_Seq_Written,
                                        int &tissue);

    void particle_Migration_between_Tissues(functions_library &functions,
                                            float **viral_Migration_Values,
                                            int *migration_start_Generation,
                                            string &source_sequence_Data_folder,
                                            vector<string> &tissue_Names,
                                            string &sequence_parent_Progeny_relationships, string &sequence_Profiles,
                                            mt19937 &gen);

    int sample_Host(functions_library &functions, float &decimal_Date,
                    vector<string> &tissue_Names,
                    string source_sequence_Data_folder, int &tissue, int &num_Samples,
                    string &sampled_sequences_Folder,
                    mt19937 &gen);

    void compress_Folder(string path, string &enable_Compression);
    void compress_Folder(string path, string enable_Compression, int thread);

    void clear_Arrays_end();

    void set_Infection_prob_Zero(string source_sequence_Data_folder,
                                 string enable_Folder_management,
                                 string enable_Compression);
};