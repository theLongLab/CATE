#include "node_within_host.cuh"

node_within_host::node_within_host()
{
    cout << "Intializing host: ";
}

void node_within_host::setHost(int host_Index, int cave_ID, int host_ID, int profile_ID, int num_Tissues)
{
    this->host_Index = host_Index;
    this->cave_ID = cave_ID;
    this->host_ID = host_ID;
    this->profile_ID = profile_ID;
    this->num_Tissues = num_Tissues;
    cout << this->cave_ID << "_" << host_ID << endl;
}

void node_within_host::setNum_Generation(int num_Generation)
{
    this->num_Generation = num_Generation;
}

void node_within_host::setInfectious_Load(int infectious_Load)
{
    this->infectious_Load = infectious_Load;
}

void node_within_host::setTerminal_Load(int terminal_Load)
{
    this->terminal_Load = terminal_Load;
}

void node_within_host::setSampling_Effect(float sampling_Effect)
{
    this->sampling_Effect = sampling_Effect;
}

void node_within_host::setCell_Limit(vector<int> cell_Limit_vec)
{
    num_Tissues = cell_Limit_vec.size();
    this->cell_Limit = (int *)malloc(sizeof(int) * num_Tissues);

    for (size_t i = 0; i < num_Tissues; i++)
    {
        cell_Limit[i] = cell_Limit_vec[i];
    }
}

void node_within_host::print_All()
{
    cout << host_Index << "\t"
         << cave_ID << "_" << host_ID << "\t"
         << profile_ID << "\t"
         << num_Generation << "\t"
         << infectious_Load << "\t"
         << terminal_Load << "\t"
         << sampling_Effect;

    for (size_t i = 0; i < num_Tissues; i++)
    {
        cout << "\t" << cell_Limit[i];
    }
    cout << endl;
}

void node_within_host::begin_Infection(functions_library &functions, string &intermediary_Sequence_location,
                                       int entry_tissues, int *entry_array, int &max_sequences_per_File,
                                       string &output_Node_location,
                                       vector<string> &tissue_Names, string first_Infection)
{
    // FIRST NODE OF INFECTION IN THE HOST

    string host_Folder = intermediary_Sequence_location + "/" + to_string(host_Index);
    functions.config_Folder(host_Folder, to_string(cave_ID) + "_" + to_string(host_ID));

    vector<vector<string>> tissue_Sequences;
    intialize_Tissues(host_Folder, tissue_Sequences, functions);

    string reference_Sequences = intermediary_Sequence_location + "/reference_Sequences";

    cout << "\nFirst infection mode: " << first_Infection << endl;

    if (first_Infection == "RANDOM")
    {

        vector<string> files;

        for (const auto &entry : filesystem::directory_iterator(reference_Sequences))
        {
            if (entry.is_regular_file() && entry.path().extension() == ".nfasta")
            {
                files.push_back(entry.path().string());
            }
        }

        vector<string> Sequences;
        vector<string> Sequence_IDs;

        cout << endl;
        vector<char> seq_Status;
        for (int file = 0; file < files.size(); file++)
        {
            cout << "Reading file: " << files[file] << endl;
            fstream nfasta;
            nfasta.open(files[file], ios::in);

            if (nfasta.is_open())
            {
                string line;
                string sequence = "";

                while (getline(nfasta, line))
                {
                    if (line.at(0) != '>')
                    {
                        sequence.append(line);
                    }
                    else
                    {
                        Sequence_IDs.push_back(line);
                        if (sequence != "")
                        {
                            Sequences.push_back(sequence);
                            sequence = "";
                        }
                    }
                }

                if (sequence != "")
                {
                    Sequences.push_back(sequence);
                    sequence = "";
                }
                nfasta.close();
            }
            else
            {
                cout << "ERROR UNABLE TO OPEN NFATSA FILE: " << files[file] << endl;
                exit(-1);
            }
        }

        random_device rd; // Will be used to obtain a seed for the random number engine
        mt19937 gen(rd());
        uniform_int_distribution<int> entry_Tissue_select(0, entry_tissues - 1);

        cout << endl;

        for (int sequence = 0; sequence < Sequences.size(); sequence++)
        {
            int tissue_Index = entry_array[entry_Tissue_select(gen)];
            cout << "Sequence " << sequence + 1 << " infects tissue: " << tissue_Index << endl;
            tissue_Sequences[tissue_Index].push_back(Sequences[sequence]);
        }

        for (int tissue = 0; tissue < entry_tissues; tissue++)
        {
            if (tissue_Sequences[entry_array[tissue]].size() > 0)
            {
                current_Viral_load_per_Tissue[entry_array[tissue]] = tissue_Sequences[entry_array[tissue]].size();
                functions.config_Folder(host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), to_string(cave_ID) + "_" + to_string(host_ID) + " Tissue " + to_string(entry_array[tissue]) + " Generation 0");

                if (!filesystem::exists(output_Node_location + "/" + get_Name()))
                {
                    functions.config_Folder(output_Node_location + "/" + get_Name(), get_Name() + " node");
                    functions.create_File(output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", "Sequence_ID\tHost\tTissue");
                    functions.create_File(output_Node_location + "/" + get_Name() + "/sequence_parent_Progeny_relationships.csv", "Source\tTarget\tType");
                }

                vector<string> sequence_Write_Store_All;
                int last_seq_Num = 0;
                functions.sequence_Write_Configurator(sequence_Write_Store_All, tissue_Sequences[entry_array[tissue]],
                                                      max_sequences_per_File, host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), last_seq_Num, seq_Status,
                                                      output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", get_Name(), tissue_Names[entry_array[tissue]], current_Generation);
                functions.partial_Write_Check(sequence_Write_Store_All,
                                              host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), last_seq_Num, seq_Status,
                                              output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", get_Name(), tissue_Names[entry_array[tissue]], current_Generation);
            }
        }
    }
    else
    {
        cout << "\nInfecting tissues\n";
        vector<string> line_Data;

        //// Write to sequence profiles file

        for (int tissue = 0; tissue < tissue_Names.size(); tissue++)
        {
            string reference_Sequences_tissue = reference_Sequences + "/" + tissue_Names[tissue];
            int last_seq_Num = 0;

            if (filesystem::exists(reference_Sequences_tissue) && filesystem::is_directory(reference_Sequences_tissue))
            {
                cout << "\nTissue: " << tissue_Names[tissue] << endl;
                for (const auto &entry : filesystem::directory_iterator(reference_Sequences_tissue))
                {
                    if (entry.is_regular_file() && entry.path().extension() == ".nfasta")
                    {
                        string file_Name = entry.path().stem();
                        functions.split(line_Data, file_Name, '_');

                        int num_Particles_Tissue = stoi(line_Data[1]) - stoi(line_Data[0]) + 1;

                        cout << "Sequences migrating: " << num_Particles_Tissue << endl;

                        current_Viral_load_per_Tissue[tissue] = current_Viral_load_per_Tissue[tissue] + num_Particles_Tissue;

                        functions.config_Folder(host_Folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation), to_string(cave_ID) + "_" + to_string(host_ID) + " Tissue " + to_string(tissue) + " Generation 0");

                        if (!filesystem::exists(output_Node_location + "/" + get_Name()))
                        {
                            functions.config_Folder(output_Node_location + "/" + get_Name(), get_Name() + " node");
                            functions.create_File(output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", "Sequence_ID\tHost\tTissue");
                            functions.create_File(output_Node_location + "/" + get_Name() + "/sequence_parent_Progeny_relationships.csv", "Source\tTarget\tType");
                        }

                        filesystem::copy_file(entry.path().string(), host_Folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation) + "/" + file_Name + ".nfasta");

                        fstream sequence_Profile;
                        sequence_Profile.open(output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", ios::app);

                        for (int sequence_Num = 0; sequence_Num < num_Particles_Tissue; sequence_Num++)
                        {
                            sequence_Profile << get_Name() << "_" << tissue_Names[tissue] << "_" << current_Generation << "_" << last_seq_Num << "\t" << get_Name() << "\t" << tissue_Names[tissue] << endl;
                            last_seq_Num++;
                        }

                        sequence_Profile.close();
                    }
                }
            }
        }
    }

    set_Infected();

    // exit(-1);
    //  for (int tissue = 0; tissue < num_Tissues; tissue++)
    //  {
    //      cout << current_Viral_load_per_Tissue[tissue] << endl;
    //  }
    //  exit(-1);
}

string node_within_host::transfer_Infection(functions_library &functions, string &intermediary_Sequence_location, string &source_Target_file_Location,
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
                                            mt19937 &gen)
{
    /**
     * First we determine of the exit tissues have a viral population to be transmitted
     **/
    if (exit_Load > 0)
    {
        cout << "\nNode " << this->cave_ID << "_" << this->host_ID << " is being infected by " << source_Name << endl;

        if (num_viruses_to_transfer > exit_Load)
        {
            num_viruses_to_transfer = exit_Load;
        }

        if (num_viruses_to_transfer > 0)
        {
            string host_Folder = intermediary_Sequence_location + "/" + to_string(host_Index);

            int year, month, day;
            functions.decimal_to_Date(decimal_Date, year, month, day);
            string date_String = to_string(year) + "-" + to_string(month) + "-" + to_string(day);

            if (current_Generation == -1)
            {
                functions.config_Folder(host_Folder, to_string(cave_ID) + "_" + to_string(host_ID));
                vector<vector<string>> tissue_Sequences;
                intialize_Tissues(host_Folder, tissue_Sequences, functions);
            }

            cout << "Attempting to transfer " << num_viruses_to_transfer << " viral particle(s)\n";

            uniform_int_distribution<> distribution_exit_Tissue(0, exit_tissues - 1);

            vector<set<int>> unique_indexes_to_Remove_Tissues;

            for (int tissue = 0; tissue < exit_tissues; tissue++)
            {
                set<int> init_Set;
                unique_indexes_to_Remove_Tissues.push_back(init_Set);
            }

            // cout << "Attempting to transfer " << num_viruses_to_transfer << " viral particle(s)\n";

            for (int particle = 0; particle < num_viruses_to_transfer; particle++)
            {
                int exit_tissue_Index = distribution_exit_Tissue(gen);
                if (source_current_Viral_load_per_Tissue[exit_array[exit_tissue_Index]] > 0)
                {
                    uniform_int_distribution<> distribution_particle(0, source_current_Viral_load_per_Tissue[exit_array[exit_tissue_Index]] - 1);
                    unique_indexes_to_Remove_Tissues[exit_tissue_Index].insert(distribution_particle(gen));
                }
            }

            cout << "Viral particle(s) and their exit tissue(s) have been indentifed\n";

            // vector<vector<int>> indexes_to_Remove;

            vector<vector<string>> seq_to_Write;
            vector<vector<string>> source_Seq_Data;
            for (int init = 0; init < entry_tissues; init++)
            {
                vector<string> initialize;
                seq_to_Write.push_back(initialize);
                source_Seq_Data.push_back(initialize);
            }

            uniform_int_distribution<> entry_Select(0, entry_tissues - 1);

            for (int tissue = 0; tissue < exit_tissues; tissue++)
            {
                vector<int> init_Tissue(unique_indexes_to_Remove_Tissues[tissue].begin(), unique_indexes_to_Remove_Tissues[tissue].end());

                if (init_Tissue.size() > 0)
                {
                    cout << "Exit tissue: " << exit_array[tissue] + 1 << endl;

                    vector<int> indexes_of_Seq_write;

                    for (int transfer_Cell = 0; transfer_Cell < init_Tissue.size(); transfer_Cell++)
                    {
                        auto it = source_removed_by_Transfer_Indexes[exit_array[tissue]].find(init_Tissue[transfer_Cell]);

                        if (it == source_removed_by_Transfer_Indexes[exit_array[tissue]].end())
                        {
                            // not present
                            indexes_of_Seq_write.push_back(init_Tissue[transfer_Cell]);
                            source_removed_by_Transfer_Indexes[exit_array[tissue]].insert(init_Tissue[transfer_Cell]);
                        }
                    }
                    if (indexes_of_Seq_write.size() > 0)
                    {
                        int valid_Sequences = 0;
                        // cout << "Collecting " << indexes_of_Seq_write.size() << " sequence(s)\n";
                        vector<string> collected_Sequences = functions.find_Sequences_Master(source_Target_file_Location, indexes_of_Seq_write, exit_array[tissue], indexed_Source_Folders[exit_array[tissue]], source_Generation, valid_Sequences);
                        cout << "Assinging sequence(s) to entry tissue(s)\n";

                        if (valid_Sequences != 0)
                        {
                            if (!filesystem::exists(output_Node_location + "/" + get_Name()))
                            {
                                functions.config_Folder(output_Node_location + "/" + get_Name(), get_Name() + " node");
                                functions.create_File(output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", "Sequence_ID\tHost\tTissue");
                                functions.create_File(output_Node_location + "/" + get_Name() + "/sequence_parent_Progeny_relationships.csv", "Source\tTarget\tType");
                            }
                        }

                        for (int check_Seq = 0; check_Seq < collected_Sequences.size(); check_Seq++)
                        {
                            if (collected_Sequences[check_Seq] != "")
                            {
                                int entry_Tissue_index = entry_Select(gen);
                                seq_to_Write[entry_Tissue_index].push_back(collected_Sequences[check_Seq]);
                                // sequence_Profile << host << "_" << tissue << "_" << last_seq_Num << "\t" << host << tissue<<endl;
                                source_Seq_Data[entry_Tissue_index].push_back(source_Name + "_" + tissue_Names[exit_array[tissue]] + "_" + to_string(source_Generation) + "_" + to_string(indexes_of_Seq_write[check_Seq]));
                            }
                        }
                    }
                }
            }
            vector<char> seq_Status;
            cout << "Writing sequence(s) to entry tissue(s)\n";
            int infected_Check = 0;
            for (int tissue = 0; tissue < entry_tissues; tissue++)
            {
                if (seq_to_Write[tissue].size() > 0)
                {
                    if (!filesystem::exists(host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation)))
                    {
                        functions.config_Folder(host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), to_string(cave_ID) + "_" + to_string(host_ID) + " Tissue " + to_string(entry_array[tissue]) + " Generation 0");
                    }

                    vector<string> sequence_Write_Store_All;
                    vector<int> indexes_Written;
                    functions.sequence_Write_Configurator_transfer(sequence_Write_Store_All, seq_to_Write[tissue],
                                                                   max_sequences_per_File, host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), current_Viral_load_per_Tissue[entry_array[tissue]], seq_Status,
                                                                   output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", get_Name(), tissue_Names[entry_array[tissue]], current_Generation,
                                                                   indexes_Written);
                    functions.partial_Write_Check_transfer(sequence_Write_Store_All,
                                                           host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), current_Viral_load_per_Tissue[entry_array[tissue]], seq_Status,
                                                           output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", get_Name(), tissue_Names[entry_array[tissue]], current_Generation, indexes_Written);
                    infected_Check = 1;

                    //(output_Node_location + "/" + get_Name() + "/sequence_parent_Progeny_relationships.csv", "Source\tTarget\tType");
                    fstream source;
                    fstream target;
                    source.open(output_Node_location + "/" + source_Name + "/sequence_parent_Progeny_relationships.csv", ios::app);
                    target.open(output_Node_location + "/" + get_Name() + "/sequence_parent_Progeny_relationships.csv", ios::app);
                    // functions.create_File(output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", "Sequence_ID\tHost\tTissue");
                    fstream target_Profiles;
                    target_Profiles.open(output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", ios::app);
                    fstream source_Profiles;
                    source_Profiles.open(output_Node_location + "/" + source_Name + "/sequence_Profiles.csv", ios::app);
                    for (int transfers = 0; transfers < indexes_Written.size(); transfers++)
                    {
                        source << source_Seq_Data[tissue][transfers] << "\t" << get_Name() << "_" << tissue_Names[entry_array[tissue]] << "_" << current_Generation << "_" << to_string(indexes_Written[transfers]) << "\tTransmission" << endl;
                        target << source_Seq_Data[tissue][transfers] << "\t" << get_Name() << "_" << tissue_Names[entry_array[tissue]] << "_" << current_Generation << "_" << to_string(indexes_Written[transfers]) << "\tTransmission" << endl;
                        target_Profiles << source_Seq_Data[tissue][transfers] << "\t" << source_Name << "\t" << tissue_Names[entry_array[tissue]] << endl;
                        source_Profiles << get_Name() << "_" << tissue_Names[entry_array[tissue]] << "_" << current_Generation << "_" << to_string(indexes_Written[transfers]) << "\t" << get_Name() << "\t" << tissue_Names[entry_array[tissue]] << endl;
                    }
                    source.close();
                    target.close();
                    target_Profiles.close();
                }
            }
            if (infected_Check == 1)
            {
                if (status == "Susceptible")
                {
                    set_Infected();
                }

                fstream write_source_Target;
                write_source_Target.open(Host_source_target_network_location, ios::app);

                if (write_source_Target.is_open())
                {
                    cout << "Writing host's source target relationship\n";
                    write_source_Target << source_Name << "\t" << this->get_Name() << "\t" << to_string(decimal_Date) << "\t" << date_String << endl;
                    write_source_Target.close();
                }
                else
                {
                    cout << "ERROR: UNABLE TO OPEN SOURCE TARGET FILE: " << Host_source_target_network_location << "\n";
                    exit(-1);
                }
            }
        }
    }
    else
    {
        cout << source_Name << " has no viral particles in the exit tissues\n";
    }

    return status;
}

int node_within_host::get_generation_Phase(int generation, int *num_replication_phases, float **tissue_replication_data, int *tissue_param_profile_Stride, int &tissue,
                                           float &variable_1, float &variable_2)
{
    /**
     * *Extract the generation of the current tissue phase to determine the parent viral population.
     **/
    cout << "Getting generation phase\n";
    int gen_Phase = -1;

    int num_Phases = num_replication_phases[(profile_ID * num_Tissues) + tissue];

    int num_phases_per_tissue = 0;
    int tissue_Check = 0;

    for (int param_Index = tissue_param_profile_Stride[profile_ID]; param_Index < tissue_param_profile_Stride[profile_ID + 1]; param_Index++)
    {
        if (tissue_Check == tissue)
        {
            float time_Check = 0;
            float current_Generation_Ratio = (float)generation / (float)num_Generation;
            for (int phases = param_Index; phases < (param_Index + num_Phases); phases++)
            {
                time_Check = time_Check + tissue_replication_data[phases][0];
                if (current_Generation_Ratio < time_Check)
                {
                    variable_1 = tissue_replication_data[phases][2];
                    variable_2 = tissue_replication_data[phases][3];
                    gen_Phase = tissue_replication_data[phases][1];
                    return gen_Phase;
                    // break;
                }
            }

            //   break;
        }

        num_phases_per_tissue++;

        if (num_phases_per_tissue == num_replication_phases[(profile_ID * num_Tissues) + tissue_Check])
        {
            num_phases_per_tissue = 0;
            tissue_Check++;
        }
    }

    return gen_Phase;
}

void node_within_host::run_Generation(functions_library &functions, string &multi_Read, int &max_Cells_at_a_time, int &gpu_Limit, int *CUDA_device_IDs, int &num_Cuda_devices, int &genome_Length,
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
                                      mt19937 &gen)
{
    /**
     * ! Main function handling the processing of the within host viral infection dynamics.
     **/

    cout << "\nSimulating generation " << current_Generation << " of " << num_Generation << " for " << get_Name() << endl
         << endl;

    if (current_Generation < num_Generation)
    {
        /**
         * First we determine if there is a within host viral population in the host and the counts present within each tissue.
         **/
        cout << "Calculating actual particles in each tissue: \n";
        ////clear array
        int *real_Particle_count_per_Tissue = (int *)malloc(sizeof(int) * num_Tissues);
        int sum_Check = 0;
        for (int tissue = 0; tissue < num_Tissues; tissue++)
        {
            real_Particle_count_per_Tissue[tissue] = current_Viral_load_per_Tissue[tissue] - removed_by_Transfer_Indexes[tissue].size() - dead_Particle_count[tissue];
            cout << tissue_Names[tissue] << " tissue: " << real_Particle_count_per_Tissue[tissue] << endl;
            sum_Check = sum_Check + real_Particle_count_per_Tissue[tissue];
        }

        if (sum_Check > 0)
        {
            if (terminal_status(terminal_tissues, terminal_array, source_sequence_Data_folder, enable_Folder_management, enable_Compression) != 1)
            {
                cout << "\nInitiating simulation\n";
                // cout << profile_ID << endl;

                vector<vector<pair<int, int>>> indexed_Source_Folders = functions.index_sequence_Folders(source_sequence_Data_folder, num_Tissues, current_Generation, multi_Read);

                string sequence_Profiles = output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv";
                string sequence_parent_Progeny_relationships = output_Node_location + "/" + get_Name() + "/sequence_parent_Progeny_relationships.csv";
                string generational_Summary = output_Node_location + "/" + get_Name() + "/node_generational_Summary.csv";

                if (!filesystem::exists(output_Node_location + "/" + get_Name()))
                {
                    functions.config_Folder(output_Node_location + "/" + get_Name(), get_Name() + " node");
                }
                if (!filesystem::exists(sequence_Profiles))
                {
                    functions.create_File(sequence_Profiles, "Sequence_ID\tHost\tTissue");
                }
                if (!filesystem::exists(sequence_parent_Progeny_relationships))
                {
                    functions.create_File(sequence_parent_Progeny_relationships, "Source\tTarget\tType");
                }
                if (!filesystem::exists(generational_Summary))
                {
                    functions.create_File(generational_Summary, "overall_Generation\tnode_Generation\tTissue\tnum_Parents\tnum_Progeny\tdead_Progeny");
                }

                for (int tissue = 0; tissue < num_Tissues; tissue++)
                {
                    if (real_Particle_count_per_Tissue[tissue] > 0)
                    {
                        // real_Particle_count_per_Tissue[tissue] = 123;

                        cout << "\nSimulating " << real_Particle_count_per_Tissue[tissue] << " particle(s) for " << tissue_Names[tissue] << " tissue\n"
                             << endl;

                        // cout << profile_ID << endl;

                        // for (int generation = current_Generation; generation < num_Generation; generation++)
                        // {
                        //     float variable_1, variable_2;
                        //     int gen_Phase = get_generation_Phase(generation, num_replication_phases, tissue_replication_data, tissue_param_profile_Stride, tissue,
                        //                                          variable_1, variable_2);
                        //     cout << "Gen phase " << generation << ": " << gen_Phase << endl;
                        //     cout << variable_1 << "\t" << variable_2 << endl;
                        // }
                        // exit(-1);

                        cout << "Identifying indexes to remove\n";
                        set<int> check_to_Remove;
                        //// Account for dead file
                        if (dead_Particle_count[tissue] > 0)
                        {
                            cout << "\nIdentifying dead viral index(es)\n";
                            // indexes_of_Dead = (int *)malloc(sizeof(int) * dead_Particle_count[tissue]);

                            fstream dead_File;
                            dead_File.open(source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation) + "/dead_List.txt");
                            if (dead_File.is_open())
                            {
                                string line;
                                // int index = 0;
                                while (getline(dead_File, line))
                                {
                                    check_to_Remove.insert(stoi(line));
                                    // index++;
                                }
                                dead_File.close();
                            }
                            else
                            {
                                cout << "ERROR: UNABLE TO OPEN DEAD LIST FILE: " << source_sequence_Data_folder << "/" << tissue << "/generation_" << current_Generation << "/dead_List.txt" << endl;
                                exit(-1);
                            }
                        }

                        if (removed_by_Transfer_Indexes[tissue].size() > 0)
                        {
                            cout << "Identifying transferred viral index(es)\n";
                            for (auto it = removed_by_Transfer_Indexes[tissue].begin(); it != removed_by_Transfer_Indexes[tissue].end(); ++it)
                            {
                                int value = *it; // Dereference the iterator to get the value
                                check_to_Remove.insert(value);
                            }
                        }

                        //// clear 2d
                        // int *parents_in_Tissue = (int *)malloc(sizeof(int) * real_Particle_count_per_Tissue[tissue]);
                        //= functions.create_INT_2D_arrays(2, real_Particle_count_per_Tissue[tissue]);

                        /**
                         * @param parents_in_Tissue : 2D array storing the parent viral ID and its cell
                         * ROW 0: PARTICLE ID
                         * ROW 1: CELL ID of the cell the viral particle occupies.
                         **/
                        int **parents_in_Tissue;
                        parents_in_Tissue = (int **)malloc(2 * sizeof(int *));
                        for (int row = 0; row < 2; row++)
                        {
                            parents_in_Tissue[row] = (int *)malloc(real_Particle_count_per_Tissue[tissue] * sizeof(int));
                        }

                        // test
                        // check_to_Remove.insert(0);
                        // check_to_Remove.insert(1);
                        // check_to_Remove.insert(5);
                        // check_to_Remove.insert(99);
                        // cout << profile_ID << endl;
                        // for (int generation = current_Generation; generation < num_Generation; generation++)
                        //{

                        /**
                         * Extract the generation of the current tissue phase to determine the parent viral population.
                         **/
                        float variable_1, variable_2;
                        int gen_Phase = get_generation_Phase(current_Generation, num_replication_phases, tissue_replication_data, tissue_param_profile_Stride, tissue,
                                                             variable_1, variable_2);

                        // cout << "Gen phase " << generation << ": " << gen_Phase << endl;
                        // cout << variable_1 << "\t" << variable_2 << endl;

                        // cout << "Gen phase " << generation << ": " << gen_Phase << endl;
                        // cout << variable_1 << "\t" << variable_2 << endl;

                        /**
                         * Viral particles infect the cells.
                         **/
                        vector<int> start_Stop_cells = assign_Cells(functions, parents_in_Tissue, real_Particle_count_per_Tissue[tissue], tissue,
                                                                    cell_Distribution_Type[profile_ID][tissue], viral_distribution_per_Tissue_param[tissue].first, viral_distribution_per_Tissue_param[tissue].second,
                                                                    check_to_Remove,
                                                                    gen_Phase, variable_1, variable_2,
                                                                    gen);

                        check_to_Remove.clear();
                        removed_by_Transfer_Indexes[tissue].clear();
                        dead_Particle_count[tissue] = 0;
                        current_Viral_load_per_Tissue[tissue] = 0;

                        cout << "Total number of cell(s) infected: " << start_Stop_cells.size() - 1 << endl;

                        // if (start_Stop_cells.size() - 1 > 0)
                        // {
                        //     for (int i = 0; i < start_Stop_cells.size() - 1; i++)
                        //     {
                        //         // cout << start_Stop_cells[i] << " : \t" << start_Stop_cells[i + 1] << endl;
                        //         for (int particle = start_Stop_cells[i]; particle < start_Stop_cells[i + 1]; particle++)
                        //         {
                        //             // cout << parents_in_Tissue[0][particle] << " :\t" << parents_in_Tissue[1][particle] << endl;
                        //             cout << parents_in_Tissue[1][particle] << "_" << parents_in_Tissue[0][particle] << ", ";
                        //         }
                        //         cout << endl;
                        //     }
                        // }
                        // }
                        // exit(-1);
                        // vector<int> start_Stop_cells;
                        if ((start_Stop_cells.size() - 1) > 0)
                        {
                            // cout << "check\n";
                            //  for (int i = 0; i < start_Stop_cells.size(); i++)
                            //  {
                            //      cout << start_Stop_cells[i] << endl;
                            //  }

                            // for (int i = 0; i < start_Stop_cells.size() - 1; i++)
                            // {
                            //     // cout << start_Stop_cells[i] << " : \t" << start_Stop_cells[i + 1] << endl;
                            //     for (int particle = start_Stop_cells[i]; particle < start_Stop_cells[i + 1]; particle++)
                            //     {
                            //         // cout << parents_in_Tissue[0][particle] << " :\t" << parents_in_Tissue[1][particle] << endl;
                            //         cout << parents_in_Tissue[1][particle] << "_" << parents_in_Tissue[0][particle] << ", ";
                            //     }
                            //     cout << endl;
                            // }

                            // exit(-1);

                            cout << "\nOrganising cell processing schedule\n"
                                 << endl;

                            vector<pair<int, int>> cells_Rounds_start_stop;

                            int full_Rounds = (start_Stop_cells.size() - 1) / max_Cells_at_a_time;
                            int partial_Rounds = (start_Stop_cells.size() - 1) % max_Cells_at_a_time;

                            for (int full = 0; full < full_Rounds; full++)
                            {
                                int start = full * max_Cells_at_a_time;
                                int stop = start + max_Cells_at_a_time;
                                cells_Rounds_start_stop.push_back(make_pair(start, stop));
                            }

                            if (partial_Rounds != 0)
                            {
                                int start = (start_Stop_cells.size() - 1) - partial_Rounds;
                                cells_Rounds_start_stop.push_back(make_pair(start, (start_Stop_cells.size() - 1)));
                            }

                            int sequence_Count = 0;
                            int index_Last_Written = 0;

                            // source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation + 1),

                            string intermediary_Tissue_folder = source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation + 1);
                            string dead_List = intermediary_Tissue_folder + "/dead_List.txt";

                            if (!filesystem::exists(intermediary_Tissue_folder))
                            {
                                functions.config_Folder(intermediary_Tissue_folder, to_string(current_Generation + 1) + " generation Tissue " + tissue_Names[tissue] + " sequences");
                            }

                            if (!filesystem::exists(dead_List))
                            {
                                functions.create_File(dead_List);
                            }

                            string cells_of_parents_location = output_Node_location + "/" + get_Name() + "/cells_of_Parents.csv";
                            string cells_of_progeny_location = output_Node_location + "/" + get_Name() + "/cells_of_Progeny.csv";

                            if (!filesystem::exists(cells_of_parents_location))
                            {
                                functions.create_File(cells_of_parents_location, "Sequence_ID\tParent_Cell_ID");
                            }
                            if (!filesystem::exists(cells_of_progeny_location))
                            {
                                functions.create_File(cells_of_progeny_location, "Sequence_ID\tProgeny_Cell_ID");
                            }

                            // exit(-1);
                            //  size_t arraySize = sizeof(parents_in_Tissue) / sizeof(parents_in_Tissue[0]);

                            // exit(-1);

                            for (int cell_Round = 0; cell_Round < cells_Rounds_start_stop.size(); cell_Round++)
                            {
                                int num_of_Cells = cells_Rounds_start_stop[cell_Round].second - cells_Rounds_start_stop[cell_Round].first;
                                cout << "\nProcessing round " << cell_Round + 1 << " of " << cells_Rounds_start_stop.size() << ": " << num_of_Cells << " cell(s)" << endl;

                                int seqeunces_to_Process = start_Stop_cells[cells_Rounds_start_stop[cell_Round].second] - start_Stop_cells[cells_Rounds_start_stop[cell_Round].first];
                                cout << "Processing " << seqeunces_to_Process << " sequence(s) in total\n";

                                simulate_Cell_replication(functions, multi_Read, gpu_Limit, CUDA_device_IDs, num_Cuda_devices, source_sequence_Data_folder, indexed_Source_Folders[tissue],
                                                          CPU_cores, max_sequences_per_File,
                                                          genome_Length,
                                                          tissue, parents_in_Tissue, tissue_Names[tissue],
                                                          start_Stop_cells, cells_Rounds_start_stop[cell_Round].first, cells_Rounds_start_stop[cell_Round].second, num_of_Cells,
                                                          seqeunces_to_Process,
                                                          sequence_Count,
                                                          Reference_fitness_survivability_proof_reading,
                                                          mutation_recombination_proof_Reading_availability,
                                                          num_effect_Segregating_sites,
                                                          sequence_Fitness_changes,
                                                          sequence_Survivability_changes,
                                                          sequence_Proof_reading_changes,
                                                          mutation_Hotspots,
                                                          mutation_hotspot_parameters,
                                                          A_0_mutation,
                                                          T_1_mutation,
                                                          G_2_mutation,
                                                          C_3_mutation,
                                                          recombination_Hotspots,
                                                          recombination_hotspot_parameters,
                                                          tot_prob_selectivity,
                                                          recombination_prob_Stride,
                                                          recombination_select_Stride,
                                                          recombination_Prob_matrix,
                                                          recombination_Select_matrix,
                                                          progeny_distribution_parameters_Array,
                                                          cells_of_parents_location,
                                                          dead_List, sequence_Profiles, sequence_parent_Progeny_relationships, cells_of_progeny_location,
                                                          index_Last_Written,
                                                          gen);
                            }
                            // free(parents_in_Tissue);

                            write_Partial_Sequence_Progeny(functions,
                                                           intermediary_Tissue_folder,
                                                           dead_List,
                                                           index_Last_Written,
                                                           tissue);

                            // cout << "\nCount\n";
                            // cout << dead_Particle_count[tissue] << endl
                            //      << current_Viral_load_per_Tissue[tissue] << endl;
                            // // TODO: COMPRESS THE PREVIOUS GENERAIONS (Current generations) SEQUENCES per tissue FOLDER
                            if (enable_Folder_management == "YES")
                            {
                                compress_Folder(source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation), enable_Compression);
                            }
                        }
                        cout << "Clearing parent cell array: ";
                        functions.clear_Array_int_CPU(parents_in_Tissue, 2);
                        cout << "cleared" << endl;
                    }
                    else
                    {
                        removed_by_Transfer_Indexes[tissue].clear();
                        dead_Particle_count[tissue] = 0;
                        current_Viral_load_per_Tissue[tissue] = 0;
                        if (enable_Folder_management == "YES")
                        {
                            if (filesystem::exists(source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation)))
                            {
                                compress_Folder(source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation), enable_Compression);
                            }
                        }
                    }
                    // cout << "Cell Limit: " << cell_Limit[tissue] << endl;

                    // cout << "Distribution type: " << cell_Distribution_Type[profile_ID][tissue] << endl;
                    // cout << viral_distribution_per_Tissue_param[tissue].first << "\t" << viral_distribution_per_Tissue_param[tissue].second << endl;

                    // // TODO: Write per node generational summary
                    fstream generational_Summary_File;
                    generational_Summary_File.open(generational_Summary, ios::app);
                    if (generational_Summary_File.is_open())
                    {
                        // overall_Generation\tnode_Generation\tTissue\tnum_Parents\tnum_Progeny\tdead_Progeny
                        generational_Summary_File << to_string(overall_Generations)
                                                  << "\t" << to_string(current_Generation)
                                                  << "\t" << tissue_Names[tissue]
                                                  << "\t" << to_string(parents_Prev_generation[tissue])
                                                  << "\t" << to_string(current_Viral_load_per_Tissue[tissue])
                                                  << "\t" << to_string(dead_Particle_count[tissue]) << endl;
                        generational_Summary_File.close();
                    }
                    else
                    {
                        cout << "ERROR: UNABLE TO OPEN NODE GENERATIONAL SUMMARY FILE: " << generational_Summary << endl;
                        exit(-1);
                    }
                }

                // particle migration between tissues
                current_Generation++;

                // for (int tissue = 0; tissue < num_Tissues; tissue++)
                // {
                //     cout << "Tissue: " << tissue << endl;
                //     cout << "\nCount\n";
                //     cout << dead_Particle_count[tissue] << endl
                //          << current_Viral_load_per_Tissue[tissue] << endl;
                // }

                // exit(-1);

                if (viral_Migration == "YES")
                {
                    // TEST SOME more
                    particle_Migration_between_Tissues(functions,
                                                       viral_Migration_Values,
                                                       migration_start_Generation,
                                                       source_sequence_Data_folder,
                                                       tissue_Names,
                                                       sequence_parent_Progeny_relationships, sequence_Profiles,
                                                       gen);
                }

                // exit(-1);
            }
        }
        else
        {
            if (infected_to_Recovered == "NO")
            {
                set_Removed();
                if (enable_Folder_management == "YES")
                {
                    compress_Folder(source_sequence_Data_folder, enable_Compression);
                }

                clear_Arrays_end();
            }
            else
            {
                set_Susceptible();
            }
        }

        free(real_Particle_count_per_Tissue);
    }
    else
    {
        set_Removed();
        if (enable_Folder_management == "YES")
        {
            compress_Folder(source_sequence_Data_folder, enable_Compression);
        }
    }

    // get each tissues generational phase
}

void node_within_host::compress_Folder(string path, string &enable_Compression)
{
    if (filesystem::exists(path))
    {
        cout << "\nCompressing folder: " << path << endl;

        string tar_Folder;
        string command_Tar;

        if (enable_Compression == "YES")
        {
            tar_Folder = path + ".tar.gz";
            command_Tar = "tar -czf " + tar_Folder + " " + path + " && rm -R " + path;
        }
        else
        {
            tar_Folder = path + ".tar";
            command_Tar = "tar -cf " + tar_Folder + " " + path + " && rm -R " + path;
        }

        int result = system(command_Tar.c_str());

        if (result == 0)
        {
            cout << "Compression successful" << endl;
        }
        else
        {
            cout << "Failed to compress the folder: " << path << endl;
            exit(-1);
        }
    }
    else
    {
        cout << "COMPRESSION ERROR: UNABLE TO FIND THE FOLDER: " << path << endl;
        exit(-1);
    }
}

void node_within_host::compress_Folder(string path, string enable_Compression, int thread)
{
    if (filesystem::exists(path))
    {
        cout << "\nCompressing folder: " << path << endl;

        string tar_Folder;
        string command_Tar;

        if (enable_Compression == "YES")
        {
            tar_Folder = path + ".tar.gz";
            command_Tar = "tar -czf " + tar_Folder + " " + path + " && rm -R " + path;
        }
        else
        {
            tar_Folder = path + ".tar";
            command_Tar = "tar -cf " + tar_Folder + " " + path + " && rm -R " + path;
        }

        int result = system(command_Tar.c_str());

        if (result == 0)
        {
            cout << "Compression successful: " << path << endl;
        }
        else
        {
            cout << "Failed to compress the folder: " << path << endl;
            exit(-1);
        }
    }
    else
    {
        cout << "COMPRESSION ERROR: UNABLE TO FIND THE FOLDER: " << path << endl;
        exit(-1);
    }
}

void node_within_host::set_Infection_prob_Zero(string source_sequence_Data_folder,
                                               string enable_Folder_management,
                                               string enable_Compression)
{
    set_Removed();
    clear_Arrays_end();
    if (enable_Folder_management == "YES")
    {
        compress_Folder(source_sequence_Data_folder, enable_Compression);
    }
}

int node_within_host::sample_Host(functions_library &functions, float &decimal_Date,
                                  vector<string> &tissue_Names,
                                  string source_sequence_Data_folder, int &tissue, int &num_Samples,
                                  string &sampled_sequences_Folder,
                                  mt19937 &gen)
{
    // set_Removed if effect of sampling is to remove

    int success_Sampling = 0;

    int tissue_Particle_Check = current_Viral_load_per_Tissue[tissue] - removed_by_Transfer_Indexes[tissue].size();

    if (tissue_Particle_Check > 0)
    {
        if (num_Samples > tissue_Particle_Check)
        {
            num_Samples = tissue_Particle_Check;
        }

        success_Sampling = 1;

        infection_probability = 1 * sampling_Effect;
        // if (infection_probability <= 0)
        // {
        //     set_Removed();
        //     clear_Arrays_end();
        //     compress_Folder(source_sequence_Data_folder);
        // }

        uniform_int_distribution<> sample_Indexes_draw(0, current_Viral_load_per_Tissue[tissue] - 1);
        set<int> sequences_to_Sample;

        cout << "Identifying " << num_Samples << " sequences\n";

        do
        {
            int potential_Sequence = sample_Indexes_draw(gen);

            auto it = removed_by_Transfer_Indexes[tissue].find(potential_Sequence);

            if (it == removed_by_Transfer_Indexes[tissue].end())
            {
                sequences_to_Sample.insert(potential_Sequence);
            }

        } while (sequences_to_Sample.size() < num_Samples);

        cout << "Collecting sequences\n";
        vector<pair<int, int>> indexed_Source_Folder = functions.index_Source_folder(source_sequence_Data_folder, tissue, current_Generation);

        vector<int> indexes_of_Seq_write(sequences_to_Sample.begin(), sequences_to_Sample.end());
        sequences_to_Sample.clear();

        vector<string> collected_Sequences = functions.find_Sequences_Master(source_sequence_Data_folder, indexes_of_Seq_write, tissue, indexed_Source_Folder, current_Generation);

        cout << "Writing sequences to folder: " << sampled_sequences_Folder << "\n";

        fstream nFASTA_file;
        nFASTA_file.open(sampled_sequences_Folder + "/sampled_Sequences_FASTA.nfasta", ios::app);
        fstream sequence_summary_File;
        sequence_summary_File.open(sampled_sequences_Folder + "/sampled_Sequences_summary.csv", ios::app);
        //"Host_ID\tTissue\tSequence_ID\tSampling_time\t_Sampling_date"

        if (nFASTA_file.is_open() && sequence_summary_File.is_open())
        {

            int year, month, day;
            functions.decimal_to_Date(decimal_Date, year, month, day);

            for (int sequence = 0; sequence < collected_Sequences.size(); sequence++)
            {
                nFASTA_file << ">" << get_Name() << "_" << tissue_Names[tissue] << "_" << to_string(current_Generation) << "_" << to_string(indexes_of_Seq_write[sequence])
                            << "_collection_date_" << to_string(decimal_Date) << "_" << to_string(year) << "-" << to_string(month) << "-" << to_string(day) << endl;
                nFASTA_file << collected_Sequences[sequence] << endl;

                sequence_summary_File << get_Name()
                                      << "\t" << tissue_Names[tissue]
                                      << "\t" << get_Name() << "_" << tissue_Names[tissue] << "_" << to_string(current_Generation) << "_" << to_string(indexes_of_Seq_write[sequence])
                                      << "\t" << to_string(decimal_Date)
                                      << "\t" << to_string(year) << "-" << to_string(month) << "-" << to_string(day) << endl;
            }

            nFASTA_file.close();
            sequence_summary_File.close();
        }
        else
        {
            cout << "ERROR: UNABLE TO OPEN ONE OF THE FOLLOWING FILES\nNFASTA SAMPLE FILE: " << sampled_sequences_Folder << "/sampled_Sequences_FASTA.nfasta\n";
            cout << "SEQUENCE SUMARY FILE: " << sampled_sequences_Folder << "/sampled_Sequences_summary.csv\n";
            exit(-1);
        }

        // Write collected sequences to a Folder

        // get_Name() + "_" + tissue_Names[destination] + "_" + to_string(current_Generation) + "_";
    }
    else
    {
        cout << "No sequences in tissue to sample\n";
    }

    return success_Sampling;
}

void node_within_host::particle_Migration_between_Tissues(functions_library &functions,
                                                          float **viral_Migration_Values,
                                                          int *migration_start_Generation,
                                                          string &source_sequence_Data_folder,
                                                          vector<string> &tissue_Names,
                                                          string &sequence_parent_Progeny_relationships, string &sequence_Profiles,
                                                          mt19937 &gen)
{
    cout << "\nIntializing viral migration\n";

    for (int migration_Check = 0; migration_Check < (num_Tissues * (num_Tissues - 1)); migration_Check++)
    {
        if (viral_Migration_Values[migration_Check][0] != -1)
        {
            if (current_Generation > migration_start_Generation[migration_Check])
            {
                int source = migration_Check / (num_Tissues - 1);
                int destination = migration_Check % (num_Tissues - 1);

                if (destination >= source)
                {
                    destination = destination + 1;
                }

                int tissue_Particle_Check = current_Viral_load_per_Tissue[source] - dead_Particle_count[source] - removed_by_Transfer_Indexes[source].size();

                if (tissue_Particle_Check > 0)
                {
                    cout << "\nViral particle(s) migrating from " << tissue_Names[source] << " tissue to " << tissue_Names[destination] << " tissue" << endl;

                    binomial_distribution<int> num_Particles(viral_Migration_Values[migration_Check][0], viral_Migration_Values[migration_Check][1]);

                    int num_viruses_to_transfer = num_Particles(gen);

                    if (num_viruses_to_transfer >= tissue_Particle_Check)
                    {
                        num_viruses_to_transfer = tissue_Particle_Check;
                    }

                    cout << "Attempting to transfer " << num_viruses_to_transfer << " viral particle(s)\n";

                    set<int> potenital_Candidates_for_transfer;
                    uniform_int_distribution<> distribution_exit_Tissue(0, current_Viral_load_per_Tissue[source] - 1);

                    cout << "Identifying sequence(s) to transfer\n";
                    for (int trial = 0; trial < num_viruses_to_transfer; trial++)
                    {
                        potenital_Candidates_for_transfer.insert(distribution_exit_Tissue(gen));
                    }

                    vector<int> Candidates_for_transfer(potenital_Candidates_for_transfer.begin(), potenital_Candidates_for_transfer.end());
                    potenital_Candidates_for_transfer.clear();

                    vector<int> indexes_of_Seq_write;

                    for (int particle = 0; particle < Candidates_for_transfer.size(); particle++)
                    {
                        auto it = removed_by_Transfer_Indexes[source].find(Candidates_for_transfer[particle]);

                        if (it == removed_by_Transfer_Indexes[source].end())
                        {
                            indexes_of_Seq_write.push_back(Candidates_for_transfer[particle]);
                            removed_by_Transfer_Indexes[source].insert(Candidates_for_transfer[particle]);
                        }
                    }

                    if (indexes_of_Seq_write.size() > 0)
                    {
                        cout << "Collecting sequence(s) to transfer\n";
                        vector<pair<int, int>> indexed_Source_Folder = functions.index_Source_folder(source_sequence_Data_folder, source, current_Generation);

                        int valid_Sequences = 0;
                        vector<string> collected_Sequences = functions.find_Sequences_Master(source_sequence_Data_folder, indexes_of_Seq_write, source, indexed_Source_Folder, current_Generation, valid_Sequences);

                        if (valid_Sequences > 0)
                        {
                            cout << "Transferring " << valid_Sequences << " sequence(s)\n";
                            string destination_Path = source_sequence_Data_folder + "/" + to_string(destination) + "/generation_" + to_string(current_Generation);
                            if (!filesystem::exists(destination_Path))
                            {
                                functions.config_Folder(destination_Path, "Destination tissue");
                            }

                            fstream sequence_parent_Progeny_relationships_File;
                            sequence_parent_Progeny_relationships_File.open(sequence_parent_Progeny_relationships, ios::app);
                            fstream sequence_Profiles_File;
                            sequence_Profiles_File.open(sequence_Profiles, ios::app);

                            if (sequence_parent_Progeny_relationships_File.is_open() && sequence_Profiles_File.is_open())
                            {
                                fstream transfer_nFasta_File;
                                transfer_nFasta_File.open(destination_Path + "/" + to_string(current_Viral_load_per_Tissue[destination]) + "_" + to_string(current_Viral_load_per_Tissue[destination] + valid_Sequences - 1) + ".nfasta", ios::out);

                                if (transfer_nFasta_File.is_open())
                                {
                                    string viral_prefix_Progeny = get_Name() + "_" + tissue_Names[destination] + "_" + to_string(current_Generation) + "_";
                                    string viral_prefix_Parent = get_Name() + "_" + tissue_Names[source] + "_" + to_string(current_Generation) + "_";

                                    for (int sequence = 0; sequence < collected_Sequences.size(); sequence++)
                                    {
                                        if (collected_Sequences[sequence] != "")
                                        {
                                            transfer_nFasta_File << ">" << current_Viral_load_per_Tissue[destination] << "_A\n";
                                            transfer_nFasta_File << collected_Sequences[sequence] << endl;

                                            sequence_Profiles_File << viral_prefix_Progeny << current_Viral_load_per_Tissue[destination]
                                                                   << "\t" << get_Name()
                                                                   << "\t" << tissue_Names[destination] << endl;

                                            sequence_parent_Progeny_relationships_File << viral_prefix_Parent << to_string(indexes_of_Seq_write[sequence])
                                                                                       << "\t" << viral_prefix_Progeny << current_Viral_load_per_Tissue[destination]
                                                                                       << "\tTransmission\n";

                                            current_Viral_load_per_Tissue[destination]++;
                                        }
                                    }
                                }
                                else
                                {
                                    cout << "ERROR: UNABLE TO OPEN NFASTA FILE: " << destination_Path << "/" << to_string(current_Viral_load_per_Tissue[destination]) << "_" + to_string(current_Viral_load_per_Tissue[destination] + valid_Sequences - 1) << ".nfasta" << endl;
                                    exit(-1);
                                }
                                sequence_parent_Progeny_relationships_File.close();
                                sequence_Profiles_File.close();
                                transfer_nFasta_File.close();
                            }
                            else
                            {
                                cout << "ERROR: UNABLE TO OPEN ONE OF THE FOLLOWING FILE:\n"
                                     << sequence_parent_Progeny_relationships << endl
                                     << sequence_Profiles << endl;
                                exit(-1);
                            }
                        }
                    }
                }
            }
        }
    }
}

void node_within_host::write_Full_Sequences_Progeny(functions_library &functions,
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
                                                    int &tissue)
{
    cout << "\nWriting parent progeny_Relationships\n";

    // exit(-1);

    fstream sequence_Profile_File;
    sequence_Profile_File.open(sequence_Profiles, ios::app);
    fstream parent_Progeny_Relationships_File;
    parent_Progeny_Relationships_File.open(sequence_parent_Progeny_relationships, ios::app);
    fstream cells_of_progeny_File;
    cells_of_progeny_File.open(cells_of_progeny_location, ios::app);

    string viral_prefix_Progeny = get_Name() + "_" + tissue_Name + "_" + to_string(current_Generation + 1) + "_";
    string viral_prefix_Parent = get_Name() + "_" + tissue_Name + "_" + to_string(current_Generation) + "_";

    //// Cells of parent and progeny

    if (sequence_Profile_File.is_open() && parent_Progeny_Relationships_File.is_open() && cells_of_progeny_File.is_open())
    {
        for (int progeny = 0; progeny < num_Progeny_being_Processed; progeny++)
        {
            converted_Sequences.push_back("");
            sequence_Profile_File << viral_prefix_Progeny << to_string(sequence_Count)
                                  << "\t" << get_Name()
                                  << "\t" << tissue_Name << endl;

            parent_Progeny_Relationships_File << viral_prefix_Parent << to_string(parent_IDs[0][progeny_Configuration_Filled[progeny][0]])
                                              << "\t" << viral_prefix_Progeny << to_string(sequence_Count)
                                              << "\tprimary_parent" << endl;

            cells_of_progeny_File << viral_prefix_Progeny << to_string(sequence_Count)
                                  << "\t" << get_Name() << "_" << tissue_Name << "_" << to_string(current_Generation) << "_" << to_string(parent_IDs[1][progeny_Configuration_Filled[progeny][0]] + start_Cell) << endl;

            for (int hotspot_parent = 1; hotspot_parent < (1 + recombination_Hotspots); hotspot_parent++)
            {
                if (progeny_Configuration_Filled[progeny][hotspot_parent] != -1)
                {
                    parent_Progeny_Relationships_File << viral_prefix_Parent << to_string(parent_IDs[0][progeny_Configuration_Filled[progeny][hotspot_parent]])
                                                      << "\t" << viral_prefix_Progeny << to_string(sequence_Count)
                                                      << "\trecombination_parent_" << to_string(hotspot_parent) << endl;
                }
            }

            sequence_Count++;
            current_Viral_load_per_Tissue[tissue] = current_Viral_load_per_Tissue[tissue] + 1;
        }

        sequence_Profile_File.close();
        parent_Progeny_Relationships_File.close();
        cells_of_progeny_File.close();
    }
    else
    {
        cout << "ERROR IN OPENING EITHER ONE OF THE FOLLOWING PROGENY SEQUENCE INORMATION FILES.\n";
        cout << "SEQUENCE PROFILE FILE:  " << sequence_Profiles << endl;
        cout << "SEQUENCE PROGENY PARENT RELATIONSHIPS FILE:  " << sequence_parent_Progeny_relationships << endl;
        cout << "SEQUENCE PROGENY CELL TRACK FILE: " << cells_of_progeny_location << endl;
        exit(-1);
    }

    cout << "Converting sequences to nFASTA\n";

    int num_per_Core = num_Progeny_being_Processed / CPU_cores;
    int remainder = num_Progeny_being_Processed % CPU_cores;

    vector<thread> threads_vec;

    shared_mutex g_mutex;
    // mutex gi_mutex;

    for (int core_ID = 0; core_ID < CPU_cores; core_ID++)
    {
        int start_Node = core_ID * num_per_Core;
        int stop_Node = start_Node + num_per_Core;

        threads_vec.push_back(thread{&node_within_host::thread_Sequence_to_String, this, start_Node, stop_Node, progeny_Sequences, genome_Length, ref(g_mutex)});
    }

    if (remainder != 0)
    {
        int start_Node = num_Progeny_being_Processed - remainder;
        int stop_Node = num_Progeny_being_Processed;

        threads_vec.push_back(thread{&node_within_host::thread_Sequence_to_String, this, start_Node, stop_Node, progeny_Sequences, genome_Length, ref(g_mutex)});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    for (int progeny = 0; progeny < num_Progeny_being_Processed; progeny++)
    {
        to_write_Sequence_Store.push_back(make_pair(to_string(Dead_or_Alive[progeny]), converted_Sequences[progeny]));
    }

    converted_Sequences.clear();
    //  int **progeny_Sequences, int *Dead_or_Alive, int **progeny_Configuration_Filled,
    functions.clear_Array_int_CPU(progeny_Sequences, num_Progeny_being_Processed);
    free(Dead_or_Alive);
    functions.clear_Array_int_CPU(progeny_Configuration_Filled, num_Progeny_being_Processed);

    if (to_write_Sequence_Store.size() >= max_sequences_per_File)
    {
        cout << "Writing nFASTA sequences: " << intermediary_Tissue_folder << endl;
        int full_Write_Count = to_write_Sequence_Store.size() / max_sequences_per_File;

        for (int full = 0; full < full_Write_Count; full++)
        {
            string fasta_file_Location = intermediary_Tissue_folder + "/" + to_string(last_index_Seq_Written) + "_" + to_string(last_index_Seq_Written + max_sequences_per_File - 1) + ".nfasta";
            fstream fasta_File;
            fasta_File.open(fasta_file_Location, ios::out);

            fstream dead_List_File;
            dead_List_File.open(dead_List, ios::app);

            if (fasta_File.is_open() && dead_List_File.is_open())
            {
                for (int write_Seq = (full * max_sequences_per_File); write_Seq < ((full * max_sequences_per_File) + max_sequences_per_File); write_Seq++)
                {
                    fasta_File << ">" << last_index_Seq_Written << "_";
                    if (to_write_Sequence_Store[write_Seq].first == "0")
                    {
                        fasta_File << "D";
                        dead_List_File << last_index_Seq_Written << endl;
                        dead_Particle_count[tissue] = dead_Particle_count[tissue] + 1;
                    }
                    else
                    {
                        fasta_File << "A";
                    }
                    fasta_File << endl;
                    fasta_File << to_write_Sequence_Store[write_Seq].second << endl;

                    last_index_Seq_Written++;
                }
                fasta_File.close();
                dead_List_File.close();
            }
            else
            {
                cout << "ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
                cout << "OR DEAD LIST FILE COULD NOT BE OPENED: " << dead_List << endl;
                exit(-1);
            }
        }

        vector<pair<string, string>> to_write_Sequence_Store_temp;
        for (int fill = full_Write_Count * max_sequences_per_File; fill < to_write_Sequence_Store.size(); fill++)
        {
            to_write_Sequence_Store_temp.push_back(make_pair(to_write_Sequence_Store[fill].first, to_write_Sequence_Store[fill].second));
        }

        to_write_Sequence_Store.clear();
        to_write_Sequence_Store = to_write_Sequence_Store_temp;
    }

    // // CREATE A FUNCTION FOR PARTIAL WRITE remainders
}

void node_within_host::write_Partial_Sequence_Progeny(functions_library &functions,
                                                      string &intermediary_Tissue_folder,
                                                      string &dead_List,
                                                      int &last_index_Seq_Written,
                                                      int &tissue)
{
    if (to_write_Sequence_Store.size() > 0)
    {
        cout << "Partial writing nFASTA sequences: " << intermediary_Tissue_folder << endl;
        string fasta_file_Location = intermediary_Tissue_folder + "/" + to_string(last_index_Seq_Written) + "_" + to_string(last_index_Seq_Written + to_write_Sequence_Store.size() - 1) + ".nfasta";
        fstream fasta_File;
        fasta_File.open(fasta_file_Location, ios::out);

        fstream dead_List_File;
        dead_List_File.open(dead_List, ios::app);

        if (fasta_File.is_open() && dead_List_File.is_open())
        {
            for (int write_Seq = 0; write_Seq < to_write_Sequence_Store.size(); write_Seq++)
            {
                fasta_File << ">" << last_index_Seq_Written << "_";
                if (to_write_Sequence_Store[write_Seq].first == "0")
                {
                    fasta_File << "D";
                    dead_List_File << last_index_Seq_Written << endl;
                    dead_Particle_count[tissue] = dead_Particle_count[tissue] + 1;
                }
                else
                {
                    fasta_File << "A";
                }
                fasta_File << endl;
                fasta_File << to_write_Sequence_Store[write_Seq].second << endl;

                last_index_Seq_Written++;
            }

            fasta_File.close();
            dead_List_File.close();
        }
        else
        {
            cout << "PARTIAL WRITE ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
            cout << "OR DEAD LIST FILE COULD NOT BE OPENED: " << dead_List << endl;
            exit(-1);
        }

        to_write_Sequence_Store.clear();
    }
}

void node_within_host::thread_Sequence_to_String(int start, int stop, int **progeny_Sequences, int genome_Length, shared_mutex &g_mutex)
{
    vector<string> converted_Sequences_Store;

    for (int progeny = start; progeny < stop; progeny++)
    {
        string sequence = "";
        for (int base = 0; base < genome_Length; base++)
        {
            sequence.append(to_string(progeny_Sequences[progeny][base]));
        }
        converted_Sequences_Store.push_back(sequence);
    }

    unique_lock<shared_mutex> ul(g_mutex);
    int index = 0;
    for (int progeny = start; progeny < stop; progeny++)
    {
        converted_Sequences[progeny] = converted_Sequences_Store[index];
        index++;
    }
}

void node_within_host::simulate_Cell_replication(functions_library &functions, string &multi_Read, int &gpu_Limit, int *CUDA_device_IDs, int &num_Cuda_devices, string &source_sequence_Data_folder, vector<pair<int, int>> &indexed_Tissue_Folder,
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
                                                 mt19937 &gen)
{
    // gpu_Limit = 5;

    int full_Rounds = Total_seqeunces_to_Process / gpu_Limit;
    int partial_Rounds = Total_seqeunces_to_Process % gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * gpu_Limit;
        int stop = start + gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = Total_seqeunces_to_Process - partial_Rounds;
        start_stops.push_back(make_pair(start, Total_seqeunces_to_Process));
    }

    cout << "Retrieving parent sequences and configuring their profiles\n";

    // // clear 2d array
    int **parent_Sequences;
    //= functions.create_INT_2D_arrays(Total_seqeunces_to_Process, genome_Length);
    parent_Sequences = (int **)malloc(Total_seqeunces_to_Process * sizeof(int *));
    for (int row = 0; row < Total_seqeunces_to_Process; row++)
    {
        parent_Sequences[row] = (int *)malloc(genome_Length * sizeof(int));
    }
    // // clear 2d array
    float **sequence_Configuration_standard;
    // int columns = 3 + (2 * recombination_Hotspots);
    // sequence_Configuration_standard = functions.create_FLOAT_2D_arrays(Total_seqeunces_to_Process, 2 + (2 * recombination_Hotspots));

    sequence_Configuration_standard = (float **)malloc(Total_seqeunces_to_Process * sizeof(float *));
    for (int row = 0; row < Total_seqeunces_to_Process; row++)
    {
        sequence_Configuration_standard[row] = (float *)malloc((2 + (2 * recombination_Hotspots)) * sizeof(float));
    }

    // float **sequence_Configuration_recombination;
    // if (mutation_recombination_proof_Reading_availability[1] == 1)
    // {
    //     sequence_Configuration_recombination = functions.create_FLOAT_2D_arrays(Total_seqeunces_to_Process, 2 * recombination_Hotspots);
    // }

    // // clear 2d array
    // int *parent_IDs = (int *)malloc(sizeof(int) * Total_seqeunces_to_Process);
    int **parent_IDs;
    // functions.create_INT_2D_arrays(2, Total_seqeunces_to_Process);
    parent_IDs = (int **)malloc(2 * sizeof(int *));
    for (int row = 0; row < 2; row++)
    {
        parent_IDs[row] = (int *)malloc(Total_seqeunces_to_Process * sizeof(int));
    }

    // cout << progeny_distribution_parameters_Array[0] << endl;
    // cout << progeny_distribution_parameters_Array[1] << endl;
    // cout << progeny_distribution_parameters_Array[2] << endl;
    // exit(-1);

    int cell_ID = 0;
    int *cell_Index = (int *)malloc(sizeof(int) * (num_Cells + 1));
    cell_Index[0] = 0;

    int augment_Value = parents_in_Tissue[1][start_Stop_cells[start_Cell]];

    int check_Cell = -1;

    // if (start_stops.size() > 1)
    // {
    //     cout << "Check: " << start_stops.size() << endl;
    //     exit(-1);
    // }

    fstream cells_of_parents_File;
    cells_of_parents_File.open(cells_of_parents_location, ios::app);
    if (cells_of_parents_File.is_open())
    {
        for (int round = 0; round < start_stops.size(); round++)
        {
            cout << "\nParent sequence processing round " << round + 1 << " of " << start_stops.size() << endl;

            vector<int> sequence_List;
            for (int parent = start_stops[round].first; parent < start_stops[round].second; parent++)
            {
                sequence_List.push_back(parents_in_Tissue[0][start_Stop_cells[start_Cell] + parent]);
                parent_IDs[0][parent] = parents_in_Tissue[0][start_Stop_cells[start_Cell] + parent];
                int current_Cell = parents_in_Tissue[1][start_Stop_cells[start_Cell] + parent];
                // parent_IDs[1][parent] = parents_in_Tissue[1][start_Stop_cells[start_Cell] + parent];

                cells_of_parents_File << get_Name() + "_" + tissue_Name + "_" + to_string(current_Generation) + "_" + to_string(parent_IDs[0][parent])
                                      << "\t" << get_Name() << "_" << tissue_Name << "_" << to_string(current_Generation) << "_" << to_string(current_Cell) << endl;

                if (parent != 0)
                {
                    if (check_Cell != current_Cell)
                    {
                        // cell_ID++;
                        // cell_Index[cell_ID] = parent;
                        cell_Index[current_Cell - augment_Value] = parent;
                        check_Cell = current_Cell;
                    }
                    // parent_IDs[1][parent] = cell_ID;
                    parent_IDs[1][parent] = current_Cell - augment_Value;
                }
                else
                {
                    // parent_IDs[1][parent] = cell_ID;
                    parent_IDs[1][parent] = current_Cell - augment_Value;
                    check_Cell = current_Cell;
                }

                if ((parent + 1) == start_stops[round].second)
                {
                    // cell_ID++;
                    // cell_Index[cell_ID] = parent + 1;
                    cell_Index[current_Cell - augment_Value + 1] = parent + 1;
                    cell_ID = current_Cell - augment_Value + 1;
                }
            }

            vector<string> collected_Sequences = functions.find_Sequences_Master(source_sequence_Data_folder, sequence_List, tissue, indexed_Tissue_Folder, current_Generation);

            if (collected_Sequences.size() == sequence_List.size())
            {
                sequence_List.clear();
                // for (int test = 0; test < collected_Sequences.size(); test++)
                // {
                //     cout << collected_Sequences[test] << endl;
                // }
                // cout << endl;

                process_Sequences_get_Configuration(functions,
                                                    collected_Sequences, CUDA_device_IDs, num_Cuda_devices, genome_Length,
                                                    Reference_fitness_survivability_proof_reading,
                                                    mutation_recombination_proof_Reading_availability,
                                                    num_effect_Segregating_sites,
                                                    sequence_Fitness_changes,
                                                    sequence_Proof_reading_changes,
                                                    recombination_Hotspots,
                                                    recombination_hotspot_parameters,
                                                    tot_prob_selectivity,
                                                    recombination_prob_Stride,
                                                    recombination_select_Stride,
                                                    recombination_Prob_matrix,
                                                    recombination_Select_matrix,
                                                    progeny_distribution_parameters_Array,
                                                    parent_Sequences,
                                                    sequence_Configuration_standard,
                                                    start_stops[round].first);
            }
            else
            {
                cout << "ERROR: WAS UNABLE TO FIND ALL REQUIRED SEQUENCES.\n";
                exit(-1);
            }

            cout << endl;
        }

        cells_of_parents_File.close();
    }
    else
    {
        cout << "ERROR: UNABLE TO OPEN PARENT CELLS TRACK FILE: " << cells_of_parents_location << endl;
        exit(-1);
    }

    // for (int i = 0; i < num_Cells; i++)
    // {
    //     cout << "Cell: " << cell_Index[i] << " to " << cell_Index[i + 1] << endl;
    //     for (int parent = cell_Index[i]; parent < cell_Index[i + 1]; parent++)
    //     {
    //         cout << parent_IDs[1][parent] << "_" << parent_IDs[0][parent] << ",";
    //     }
    //     cout << endl;
    // }

    cout << "Cells check: " << cell_ID << " : " << num_Cells << endl
         << endl;

    // exit(-1);

    // cout << endl;

    // for (int row = 0; row < Total_seqeunces_to_Process; row++)
    // {
    //     for (int col = 0; col < genome_Length; col++)
    //     {
    //         cout << parent_Sequences[row][col];
    //     }
    //     cout << endl;
    // }

    // cout << endl;

    if (cell_ID != num_Cells)
    {
        cout << "Check: " << start_stops.size() << endl;
        cout << "Cells check: " << cell_ID << " : " << num_Cells << endl
             << endl;
        exit(-1);
    }

    // exit(-1);

    // for (int row = 0; row < Total_seqeunces_to_Process; row++)
    // {
    //     for (int col = 0; col < (2 + (2 * recombination_Hotspots)); col++)
    //     {
    //         cout << sequence_Configuration_standard[row][col] << "\t";
    //     }
    //     cout << endl;
    // }

    // exit(-1);

    cout << "All parent sequences configured\n";
    int total_Progeny = 0;
    // // clear 2d array
    // float *totals_Progeny_Selectivity = (float *)malloc(sizeof(float) * (1 + recombination_Hotspots));
    float **totals_Progeny_Selectivity = functions.create_Fill_2D_array_FLOAT(num_Cells, recombination_Hotspots, 0);
    ////  clear 1d array
    int *progeny_Stride = (int *)malloc(sizeof(int) * (Total_seqeunces_to_Process + 1));
    progeny_Stride[0] = 0;

    // for (int fill = 0; fill < (1 + recombination_Hotspots); fill++)
    // {
    //     totals_Progeny_Selectivity[fill] = 0;
    // }

    cout << "\nDetermining total progeny and configuring recombination hotspots\n";

    for (int cell = 0; cell < num_Cells; cell++)
    {
        for (int parent = cell_Index[cell]; parent < cell_Index[cell + 1]; parent++)
        {
            // cout << parent_IDs[1][parent] << "_" << parent_IDs[0][parent] << ",";
            progeny_Stride[parent + 1] = progeny_Stride[parent] + sequence_Configuration_standard[parent][0];

            for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
            {
                totals_Progeny_Selectivity[cell][hotspot] = totals_Progeny_Selectivity[cell][hotspot] + sequence_Configuration_standard[parent][(hotspot * 2) + 3];
            }
        }
        // cout << endl;
    }

    // for (int row = 0; row < Total_seqeunces_to_Process; row++)
    // {
    //     // totals_Progeny_Selectivity[0] = totals_Progeny_Selectivity[0] + sequence_Configuration_standard[row][0];
    //     progeny_Stride[row + 1] = progeny_Stride[row] + sequence_Configuration_standard[row][0];

    //     for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
    //     {
    //         totals_Progeny_Selectivity[hotspot + 1] = totals_Progeny_Selectivity[hotspot + 1] + sequence_Configuration_standard[row][(hotspot * 2) + 3];
    //     }
    // }

    total_Progeny = progeny_Stride[Total_seqeunces_to_Process];
    cout << "Total progeny to be simulated: " << total_Progeny << endl;

    // cout << endl;
    // for (int i = 0; i < num_Cells; i++)
    // {
    //     for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
    //     {
    //         cout << totals_Progeny_Selectivity[i][hotspot] << "\t";
    //     }
    //     cout << endl;
    // }

    cout << "\nIntiating Progeny configurations\n";
    cout << "Intializing GPU memory structures\n";

    cudaSetDevice(CUDA_device_IDs[0]);

    // //clear array
    // int **cuda_progeny_Configuration = functions.create_CUDA_2D_int(total_Progeny, 1 + recombination_Hotspots);
    int **progeny_Configuration;
    // = functions.create_INT_2D_arrays(total_Progeny, 1 + recombination_Hotspots);
    progeny_Configuration = (int **)malloc(total_Progeny * sizeof(int *));
    for (int row = 0; row < total_Progeny; row++)
    {
        progeny_Configuration[row] = (int *)malloc((1 + recombination_Hotspots) * sizeof(int));
    }

    //// clear array
    int *cuda_progeny_Stride;
    cudaMallocManaged(&cuda_progeny_Stride, (Total_seqeunces_to_Process + 1) * sizeof(int));
    cudaMemcpy(cuda_progeny_Stride, progeny_Stride, (Total_seqeunces_to_Process + 1) * sizeof(int), cudaMemcpyHostToDevice);
    free(progeny_Stride);

    // // ! clear
    float **cuda_sequence_Configuration_standard;
    // functions.float_2D_Array_load_to_CUDA(sequence_Configuration_standard, Total_seqeunces_to_Process, 2 + (2 * recombination_Hotspots));

    cudaMallocManaged(&cuda_sequence_Configuration_standard, Total_seqeunces_to_Process * sizeof(float *));
    for (int row = 0; row < Total_seqeunces_to_Process; row++)
    {
        cudaMalloc((void **)&(cuda_sequence_Configuration_standard[row]), (2 + (2 * recombination_Hotspots)) * sizeof(float));
        cudaMemcpy(cuda_sequence_Configuration_standard[row], sequence_Configuration_standard[row], (2 + (2 * recombination_Hotspots)) * sizeof(float), cudaMemcpyHostToDevice);
    }

    // for (int row = 0; row < Total_seqeunces_to_Process; row++)
    // {
    //     cudaMemcpy(cuda_sequence_Configuration_standard[row], sequence_Configuration_standard[row], (2 + (2 * recombination_Hotspots)) * sizeof(float), cudaMemcpyHostToDevice);
    // }

    for (int round = 0; round < start_stops.size(); round++)
    {
        cout << "\nParent sequence processing round " << round + 1 << " of " << start_stops.size() << endl;

        progeny_Configurator(functions,
                             cuda_sequence_Configuration_standard, recombination_Hotspots,
                             start_stops[round].first, start_stops[round].second - start_stops[round].first,
                             progeny_Configuration, cuda_progeny_Stride, cuda_progeny_Stride[start_stops[round].second] - cuda_progeny_Stride[start_stops[round].first], cuda_progeny_Stride[start_stops[round].first]);
    }

    // cout << "\nCopying test\n";
    //  int **progeny_Configuration = functions.load_to_Host(cuda_progeny_Configuration, total_Progeny, 1 + recombination_Hotspots);

    // for (int row = 0; row < total_Progeny; row++)
    // {
    //     for (int col = 0; col < (1 + recombination_Hotspots); col++)
    //     {
    //         cout << progeny_Configuration[row][col] << "\t";
    //     }
    //     cout << endl;
    // }

    // functions.clear_Array_INT(cuda_progeny_Configuration, total_Progeny);
    cudaFree(cuda_progeny_Stride);

    // Free each row
    for (int row = 0; row < Total_seqeunces_to_Process; row++)
    {
        cudaFree(cuda_sequence_Configuration_standard[row]);
    }
    // see
    cudaFree(cuda_sequence_Configuration_standard);

    // Free the array of pointers

    // create and save sequence
    start_stops.clear();

    // exit(-1);

    cout << "\nScheduling progeny processing schedule\n";

    full_Rounds = total_Progeny / gpu_Limit;
    partial_Rounds = total_Progeny % gpu_Limit;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * gpu_Limit;
        int stop = start + gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = total_Progeny - partial_Rounds;
        start_stops.push_back(make_pair(start, total_Progeny));
    }

    for (int round = 0; round < start_stops.size(); round++)
    {
        cout << "\nProgeny completition round " << round + 1 << " of " << start_stops.size() << endl;
        progeny_Completion(functions,
                           CUDA_device_IDs, num_Cuda_devices,
                           CPU_cores, max_sequences_per_File,
                           genome_Length, Reference_fitness_survivability_proof_reading, mutation_recombination_proof_Reading_availability,
                           num_effect_Segregating_sites,
                           sequence_Survivability_changes,
                           recombination_Hotspots,
                           recombination_hotspot_parameters,
                           tot_prob_selectivity,
                           mutation_Hotspots,
                           A_0_mutation,
                           T_1_mutation,
                           G_2_mutation,
                           C_3_mutation,
                           mutation_hotspot_parameters,
                           parent_Sequences, Total_seqeunces_to_Process, sequence_Configuration_standard, parent_IDs, num_Cells, cell_Index, start_Cell,
                           progeny_Configuration, start_stops[round].second - start_stops[round].first,
                           totals_Progeny_Selectivity,
                           start_stops[round].first, start_stops[round].second, sequence_Count, source_sequence_Data_folder,
                           tissue, tissue_Name,
                           dead_List, sequence_Profiles, sequence_parent_Progeny_relationships, cells_of_progeny_location,
                           index_Last_Written);
    }

    cout << "Unloading memory data points\n";
    functions.clear_Array_int_CPU(parent_Sequences, Total_seqeunces_to_Process);
    functions.clear_Array_float_CPU(sequence_Configuration_standard, Total_seqeunces_to_Process);
    functions.clear_Array_int_CPU(parent_IDs, 2);

    free(cell_Index);

    // functions.clear_Array_float_CPU(totals_Progeny_Selectivity, num_Cells);
    for (int row = 0; row < num_Cells; row++)
    {
        free(totals_Progeny_Selectivity[row]);
    }
    free(totals_Progeny_Selectivity);

    functions.clear_Array_int_CPU(progeny_Configuration, total_Progeny);
}

__global__ void cuda_Progeny_Complete_Configuration(int genome_Length,
                                                    float *cuda_Reference_fitness_survivability_proof_reading,
                                                    int *cuda_num_effect_Segregating_sites,
                                                    float **cuda_sequence_Survivability_changes,
                                                    int recombination_Hotspots,
                                                    float **cuda_recombination_hotspot_parameters,
                                                    int *cuda_tot_prob_selectivity,
                                                    int mutation_Hotspots,
                                                    float **cuda_A_0_mutation,
                                                    float **cuda_T_1_mutation,
                                                    float **cuda_G_2_mutation,
                                                    float **cuda_C_3_mutation,
                                                    float **cuda_mutation_hotspot_parameters,
                                                    int **cuda_parent_Sequences, int **cuda_parent_IDs,
                                                    float **cuda_sequence_Configuration_standard,
                                                    int *cuda_cell_Index, int num_Cells,
                                                    float **cuda_totals_Progeny_Selectivity,
                                                    int **cuda_progeny_Configuration,
                                                    int **cuda_progeny_Sequences,
                                                    int *cuda_Dead_or_Alive,
                                                    int per_gpu_Progeny)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < per_gpu_Progeny)
    {
        for (int base = 0; base < genome_Length; base++)
        {
            cuda_progeny_Sequences[tid][base] = cuda_parent_Sequences[cuda_progeny_Configuration[tid][0]][base];
        }

        curandState localState;
        curand_init(clock64(), tid, 0, &localState);

        if (recombination_Hotspots > 0)
        {
            int get_Cell = cuda_parent_IDs[1][cuda_progeny_Configuration[tid][0]];

            for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
            {
                if (cuda_progeny_Configuration[tid][hotspot + 1] != -1)
                {
                    // TEST BLOCK 1
                    float rand_num = curand_uniform(&localState);
                    float cumulative_prob = 0.0f;

                    //! CHECK LATER
                    int recomb_parent = cuda_progeny_Configuration[tid][0];

                    for (int check = cuda_cell_Index[get_Cell]; check < cuda_cell_Index[get_Cell + 1]; check++)
                    {
                        cumulative_prob += (cuda_sequence_Configuration_standard[check][(hotspot * 2) + 3] / cuda_totals_Progeny_Selectivity[get_Cell][hotspot]);
                        if (rand_num < cumulative_prob)
                        {
                            recomb_parent = check;
                            break;
                        }
                    }

                    cuda_progeny_Configuration[tid][hotspot + 1] = recomb_parent;

                    // TEST BLOCK 2
                    if (recomb_parent != cuda_progeny_Configuration[tid][0])
                    {
                        for (int base = ((int)cuda_recombination_hotspot_parameters[hotspot][0] - 1); base < (int)cuda_recombination_hotspot_parameters[hotspot][1]; base++)
                        {
                            cuda_progeny_Sequences[tid][base] = cuda_parent_Sequences[recomb_parent][base];
                        }
                    }
                }
            }
        }

        // TEST BLOCK 3
        if (mutation_Hotspots > 0)
        {
            for (int hotspot = 0; hotspot < mutation_Hotspots; hotspot++)
            {
                int num_Mutations = -1;

                if (cuda_mutation_hotspot_parameters[hotspot][2] == 0)
                {
                    // Poisson
                    num_Mutations = curand_poisson(&localState, cuda_mutation_hotspot_parameters[hotspot][3]);
                }
                else if (cuda_mutation_hotspot_parameters[hotspot][2] == 1)
                {
                    // neg binomial
                    int failures = 0;
                    int successes = 0;

                    while (successes < cuda_mutation_hotspot_parameters[hotspot][3])
                    {
                        float rand_num = curand_uniform(&localState);
                        if (rand_num < cuda_mutation_hotspot_parameters[hotspot][4])
                        {
                            successes++;
                        }
                        else
                        {
                            failures++;
                        }
                    }

                    num_Mutations = failures;
                }
                else
                {
                    // fixed or binomial distribution
                    int count = 0;

                    int bases_in_Region = cuda_mutation_hotspot_parameters[hotspot][1] - (cuda_mutation_hotspot_parameters[hotspot][0] - 1);

                    for (int trial = 0; trial < bases_in_Region; trial++)
                    {
                        if (curand_uniform(&localState) < cuda_mutation_hotspot_parameters[hotspot][3])
                        {
                            count++;
                        }
                    }

                    num_Mutations = count;
                }

                if (num_Mutations > 0)
                {
                    if (cuda_sequence_Configuration_standard[cuda_progeny_Configuration[tid][0]][1] != -1)
                    {
                        int count = 0;

                        // int bases_in_Region = cuda_mutation_hotspot_parameters[hotspot][1] - (cuda_mutation_hotspot_parameters[hotspot][0] - 1);

                        for (int trial = 0; trial < num_Mutations; trial++)
                        {
                            if (curand_uniform(&localState) < cuda_sequence_Configuration_standard[cuda_progeny_Configuration[tid][0]][1])
                            {
                                count++;
                            }
                        }
                        num_Mutations = num_Mutations - count;
                    }

                    if (num_Mutations > 0)
                    {
                        for (int mutation = 0; mutation < num_Mutations; mutation++)
                        {
                            int position = (int)(curand_uniform(&localState) * (((int)cuda_mutation_hotspot_parameters[hotspot][1] - 1) - ((int)cuda_mutation_hotspot_parameters[hotspot][0] - 1) + 1)) + ((int)cuda_mutation_hotspot_parameters[hotspot][0] - 1);

                            float rand_num = curand_uniform(&localState);
                            float cumulative_prob = 0.0f;

                            int original_BASE = cuda_progeny_Sequences[tid][position];
                            int new_Base = 0;

                            for (int base = 0; base < 4; base++)
                            {
                                cumulative_prob += (original_BASE == 0)   ? cuda_A_0_mutation[hotspot][base]
                                                   : (original_BASE == 1) ? cuda_T_1_mutation[hotspot][base]
                                                   : (original_BASE == 2) ? cuda_G_2_mutation[hotspot][base]
                                                   : (original_BASE == 3) ? cuda_C_3_mutation[hotspot][base]
                                                                          : 0.0f;

                                if (rand_num < cumulative_prob)
                                {
                                    new_Base = base;
                                    break;
                                }
                            }

                            cuda_progeny_Sequences[tid][position] = new_Base;
                        }
                    }
                }
            }
        }

        // TEST BLOCK 4
        // Determine survivability
        float survivability = cuda_Reference_fitness_survivability_proof_reading[1];

        for (int pos = 0; pos < cuda_num_effect_Segregating_sites[1]; pos++)
        {
            if (cuda_progeny_Sequences[tid][(int)cuda_sequence_Survivability_changes[pos][0] - 1] == 0)
            {
                survivability = survivability + cuda_sequence_Survivability_changes[pos][1];
            }
            else if (cuda_progeny_Sequences[tid][(int)cuda_sequence_Survivability_changes[pos][0] - 1] == 1)
            {
                survivability = survivability + cuda_sequence_Survivability_changes[pos][2];
            }
            else if (cuda_progeny_Sequences[tid][(int)cuda_sequence_Survivability_changes[pos][0] - 1] == 2)
            {
                survivability = survivability + cuda_sequence_Survivability_changes[pos][3];
            }
            else if (cuda_progeny_Sequences[tid][(int)cuda_sequence_Survivability_changes[pos][0] - 1] == 3)
            {
                survivability = survivability + cuda_sequence_Survivability_changes[pos][4];
            }
        }

        if (survivability >= 1)
        {
            cuda_Dead_or_Alive[tid] = 1;
        }
        else if (survivability <= 0)
        {
            cuda_Dead_or_Alive[tid] = 0;
        }
        else
        {
            float survivability_Check = curand_uniform(&localState);
            cuda_Dead_or_Alive[tid] = (survivability_Check < survivability) ? 1 : 0;
        }

        tid += blockDim.x * gridDim.x;
    }
}

void node_within_host::progeny_Completion(functions_library &functions,
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
                                          int &index_Last_Written)
{
    int num_of_Sequences_current = num_Progeny_being_Processed;
    cout << "\nConfiguring multi gpu distribution of " << num_of_Sequences_current << " sequence(s)\n";

    int standard_num_per_GPU = num_of_Sequences_current / num_Cuda_devices;
    int remainder = num_of_Sequences_current % num_Cuda_devices;

    vector<pair<int, int>> start_stop_Per_GPU;

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        int start = gpu * standard_num_per_GPU;
        int stop = start + standard_num_per_GPU;

        start_stop_Per_GPU.push_back(make_pair(start, stop));
    }

    start_stop_Per_GPU[num_Cuda_devices - 1].second = start_stop_Per_GPU[num_Cuda_devices - 1].second + remainder;

    cudaStream_t streams[num_Cuda_devices];
    cudaDeviceProp deviceProp;

    float *cuda_Reference_fitness_survivability_proof_reading[num_Cuda_devices];

    int *cuda_num_effect_Segregating_sites[num_Cuda_devices];
    float **cuda_sequence_Survivability_changes[num_Cuda_devices];

    float **cuda_recombination_hotspot_parameters[num_Cuda_devices];
    int *cuda_tot_prob_selectivity[num_Cuda_devices];

    float **cuda_A_0_mutation[num_Cuda_devices];
    float **cuda_T_1_mutation[num_Cuda_devices];
    float **cuda_G_2_mutation[num_Cuda_devices];
    float **cuda_C_3_mutation[num_Cuda_devices];

    float **cuda_mutation_hotspot_parameters[num_Cuda_devices];

    int **cuda_parent_Sequences[num_Cuda_devices];
    float **cuda_sequence_Configuration_standard[num_Cuda_devices];
    int **cuda_parent_IDs[num_Cuda_devices];

    int *cuda_cell_Index[num_Cuda_devices];

    float **cuda_totals_Progeny_Selectivity[num_Cuda_devices];

    int **cuda_progeny_Configuration[num_Cuda_devices];

    int **cuda_progeny_Sequences[num_Cuda_devices];
    int *cuda_Dead_or_Alive[num_Cuda_devices];

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cudaGetDeviceProperties(&deviceProp, gpu);
        cout << "Intializing GPU " << CUDA_device_IDs[gpu] << "'s stream: " << deviceProp.name << endl;

        // cout << "Check 1" << endl;

        cudaMallocManaged(&cuda_progeny_Sequences[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(int *));
        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaMalloc((void **)&cuda_progeny_Sequences[gpu][row], genome_Length * sizeof(int));
        }

        // cout << "Check 2" << endl;

        cudaMallocManaged(&cuda_Dead_or_Alive[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(int));

        // cout << "Check 3" << endl;

        cudaMallocManaged(&cuda_Reference_fitness_survivability_proof_reading[gpu], 3 * sizeof(float));
        cudaMemcpy(cuda_Reference_fitness_survivability_proof_reading[gpu], Reference_fitness_survivability_proof_reading, 3 * sizeof(float), cudaMemcpyHostToDevice);

        // cout << "Check 4" << endl;

        cudaMallocManaged(&cuda_num_effect_Segregating_sites[gpu], 3 * sizeof(int));
        cudaMemcpy(cuda_num_effect_Segregating_sites[gpu], num_effect_Segregating_sites, 3 * sizeof(int), cudaMemcpyHostToDevice);

        // cout << "Check 5" << endl;

        cudaMallocManaged(&cuda_sequence_Survivability_changes[gpu], num_effect_Segregating_sites[1] * sizeof(float *));
        for (int row = 0; row < num_effect_Segregating_sites[1]; row++)
        {
            cudaMalloc((void **)&cuda_sequence_Survivability_changes[gpu][row], 5 * sizeof(float));
            cudaMemcpy(cuda_sequence_Survivability_changes[gpu][row], sequence_Survivability_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
        }

        // cout << "Check 6" << endl;

        cudaMallocManaged(&cuda_tot_prob_selectivity[gpu], 2 * sizeof(int));

        // cout << "Check 7" << endl;

        cudaMallocManaged(&cuda_recombination_hotspot_parameters[gpu], recombination_Hotspots * sizeof(float *));
        if (recombination_Hotspots > 0)
        {
            cudaMemcpy(cuda_tot_prob_selectivity[gpu], tot_prob_selectivity, 2 * sizeof(int), cudaMemcpyHostToDevice);

            for (int row = 0; row < recombination_Hotspots; row++)
            {
                cudaMalloc((void **)&cuda_recombination_hotspot_parameters[gpu][row], 4 * sizeof(float));
                cudaMemcpy(cuda_recombination_hotspot_parameters[gpu][row], recombination_hotspot_parameters[row], 4 * sizeof(float), cudaMemcpyHostToDevice);
            }
        }

        // cout << "Check 8" << endl;

        cudaMallocManaged(&cuda_A_0_mutation[gpu], mutation_Hotspots * sizeof(float *));
        cudaMallocManaged(&cuda_T_1_mutation[gpu], mutation_Hotspots * sizeof(float *));
        cudaMallocManaged(&cuda_G_2_mutation[gpu], mutation_Hotspots * sizeof(float *));
        cudaMallocManaged(&cuda_C_3_mutation[gpu], mutation_Hotspots * sizeof(float *));

        // cout << "Check 9" << endl;

        cudaMallocManaged(&cuda_mutation_hotspot_parameters[gpu], mutation_Hotspots * sizeof(float *));

        if (mutation_Hotspots > 0)
        {
            for (int row = 0; row < mutation_Hotspots; row++)
            {
                cudaMalloc((void **)&cuda_A_0_mutation[gpu][row], 4 * sizeof(float));
                cudaMemcpy(cuda_A_0_mutation[gpu][row], A_0_mutation[row], 4 * sizeof(float), cudaMemcpyHostToDevice);

                cudaMalloc((void **)&cuda_T_1_mutation[gpu][row], 4 * sizeof(float));
                cudaMemcpy(cuda_T_1_mutation[gpu][row], T_1_mutation[row], 4 * sizeof(float), cudaMemcpyHostToDevice);

                cudaMalloc((void **)&cuda_G_2_mutation[gpu][row], 4 * sizeof(float));
                cudaMemcpy(cuda_G_2_mutation[gpu][row], G_2_mutation[row], 4 * sizeof(float), cudaMemcpyHostToDevice);

                cudaMalloc((void **)&cuda_C_3_mutation[gpu][row], 4 * sizeof(float));
                cudaMemcpy(cuda_C_3_mutation[gpu][row], C_3_mutation[row], 4 * sizeof(float), cudaMemcpyHostToDevice);

                cudaMalloc((void **)&cuda_mutation_hotspot_parameters[gpu][row], 5 * sizeof(float));
                cudaMemcpy(cuda_mutation_hotspot_parameters[gpu][row], mutation_hotspot_parameters[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
            }
        }
        // cout << "Check 10" << endl;

        cudaMallocManaged(&cuda_parent_Sequences[gpu], num_Parent_sequence * sizeof(int *));
        // cout << "Check 11" << endl;
        cudaMallocManaged(&cuda_sequence_Configuration_standard[gpu], num_Parent_sequence * sizeof(float *));
        for (int row = 0; row < num_Parent_sequence; row++)
        {
            cudaMalloc((void **)&cuda_parent_Sequences[gpu][row], genome_Length * sizeof(int));
            cudaMemcpy(cuda_parent_Sequences[gpu][row], parent_Sequences[row], genome_Length * sizeof(int), cudaMemcpyHostToDevice);

            cudaMalloc((void **)&cuda_sequence_Configuration_standard[gpu][row], (2 + (2 * recombination_Hotspots)) * sizeof(float));
            cudaMemcpy(cuda_sequence_Configuration_standard[gpu][row], sequence_Configuration_standard[row], (2 + (2 * recombination_Hotspots)) * sizeof(float), cudaMemcpyHostToDevice);
        }
        // cout << "Check 12" << endl;
        cudaMallocManaged(&cuda_parent_IDs[gpu], 2 * sizeof(int *));
        for (int row = 0; row < 2; row++)
        {
            cudaMalloc((void **)&cuda_parent_IDs[gpu][row], num_Parent_sequence * sizeof(int));
            cudaMemcpy(cuda_parent_IDs[gpu][row], parent_IDs[row], num_Parent_sequence * sizeof(int), cudaMemcpyHostToDevice);
        }

        // cout << "Check 13" << endl;

        // cudaMalloc(&cuda_cell_Index[gpu], (num_Cells + 1) * sizeof(int));
        cudaMallocManaged(&cuda_cell_Index[gpu], (num_Cells + 1) * sizeof(int));
        cudaMemcpy(cuda_cell_Index[gpu], cell_Index, (num_Cells + 1) * sizeof(int), cudaMemcpyHostToDevice);

        // cout << "Check 14" << endl;

        // num_Cells, recombination_Hotspots,

        cudaMallocManaged(&cuda_totals_Progeny_Selectivity[gpu], num_Cells * sizeof(float *));
        for (int row = 0; row < num_Cells; row++)
        {
            cudaMalloc((void **)&cuda_totals_Progeny_Selectivity[gpu][row], recombination_Hotspots * sizeof(float));
            cudaMemcpy(cuda_totals_Progeny_Selectivity[gpu][row], totals_Progeny_Selectivity[row], recombination_Hotspots * sizeof(float), cudaMemcpyHostToDevice);
        }
        // cout << "Check 15" << endl;

        cudaMallocManaged(&cuda_progeny_Configuration[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(int *));
        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaMalloc((void **)&cuda_progeny_Configuration[gpu][row], (1 + recombination_Hotspots) * sizeof(int));
            cudaMemcpy(cuda_progeny_Configuration[gpu][row], progeny_Configuration[row + start_stop_Per_GPU[gpu].first + start_Progeny], (1 + recombination_Hotspots) * sizeof(int), cudaMemcpyHostToDevice);
        }

        // cout << "Check 16" << endl;

        cudaStreamCreate(&streams[gpu]);
    }

    cout << "Loaded " << num_Progeny_being_Processed << " sequence(s) and all pre-requisites to the GPU(s)\nInitiating GPU(s) execution\n";

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cuda_Progeny_Complete_Configuration<<<functions.tot_Blocks_array[gpu], functions.tot_ThreadsperBlock_array[gpu], 0, streams[gpu]>>>(genome_Length,
                                                                                                                                            cuda_Reference_fitness_survivability_proof_reading[gpu],
                                                                                                                                            cuda_num_effect_Segregating_sites[gpu],
                                                                                                                                            cuda_sequence_Survivability_changes[gpu],
                                                                                                                                            recombination_Hotspots,
                                                                                                                                            cuda_recombination_hotspot_parameters[gpu],
                                                                                                                                            cuda_tot_prob_selectivity[gpu],
                                                                                                                                            mutation_Hotspots,
                                                                                                                                            cuda_A_0_mutation[gpu],
                                                                                                                                            cuda_T_1_mutation[gpu],
                                                                                                                                            cuda_G_2_mutation[gpu],
                                                                                                                                            cuda_C_3_mutation[gpu],
                                                                                                                                            cuda_mutation_hotspot_parameters[gpu],
                                                                                                                                            cuda_parent_Sequences[gpu], cuda_parent_IDs[gpu],
                                                                                                                                            cuda_sequence_Configuration_standard[gpu],
                                                                                                                                            cuda_cell_Index[gpu], num_Cells,
                                                                                                                                            cuda_totals_Progeny_Selectivity[gpu],
                                                                                                                                            cuda_progeny_Configuration[gpu],
                                                                                                                                            cuda_progeny_Sequences[gpu],
                                                                                                                                            cuda_Dead_or_Alive[gpu],
                                                                                                                                            start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first);
    }

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cudaStreamSynchronize(streams[gpu]);
    }

    cout << "GPU(s) streams completed and synchronized\nCopying data from GPU to Host memory\n";

    //// cleared all 3
    int **progeny_Configuration_Filled;
    // = functions.create_INT_2D_arrays_for_GPU(num_Progeny_being_Processed, (1 + recombination_Hotspots));
    progeny_Configuration_Filled = (int **)malloc(num_Progeny_being_Processed * sizeof(int *));
    for (int row = 0; row < num_Progeny_being_Processed; row++)
    {
        progeny_Configuration_Filled[row] = (int *)malloc((1 + recombination_Hotspots) * sizeof(int));
    }

    int **progeny_Sequences;
    // functions.create_INT_2D_arrays_for_GPU(num_Progeny_being_Processed, genome_Length);
    progeny_Sequences = (int **)malloc(num_Progeny_being_Processed * sizeof(int *));
    for (int row = 0; row < num_Progeny_being_Processed; row++)
    {
        progeny_Sequences[row] = (int *)malloc(genome_Length * sizeof(int));
    }

    int *Dead_or_Alive = (int *)malloc(sizeof(int) * num_Progeny_being_Processed);

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);

        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaMemcpy(progeny_Configuration_Filled[start_stop_Per_GPU[gpu].first + row], cuda_progeny_Configuration[gpu][row], (1 + recombination_Hotspots) * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(progeny_Sequences[start_stop_Per_GPU[gpu].first + row], cuda_progeny_Sequences[gpu][row], genome_Length * sizeof(int), cudaMemcpyDeviceToHost);
        }

        cudaMemcpy(Dead_or_Alive + start_stop_Per_GPU[gpu].first, cuda_Dead_or_Alive[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(int), cudaMemcpyDeviceToHost);
    }

    cout << "Data received by host\n";

    cout << "Terminating GPU streams: ";
    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);

        cudaFree(cuda_Reference_fitness_survivability_proof_reading[gpu]);

        for (int row = 0; row < num_effect_Segregating_sites[1]; row++)
        {
            cudaFree(cuda_sequence_Survivability_changes[gpu][row]);
        }
        cudaFree(cuda_sequence_Survivability_changes[gpu]);

        for (int row = 0; row < recombination_Hotspots; row++)
        {
            cudaFree(cuda_recombination_hotspot_parameters[gpu][row]);
        }
        cudaFree(cuda_recombination_hotspot_parameters[gpu]);

        cudaFree(cuda_tot_prob_selectivity[gpu]);

        for (int row = 0; row < mutation_Hotspots; row++)
        {
            cudaFree(cuda_A_0_mutation[gpu][row]);
            cudaFree(cuda_T_1_mutation[gpu][row]);
            cudaFree(cuda_G_2_mutation[gpu][row]);
            cudaFree(cuda_C_3_mutation[gpu][row]);

            cudaFree(cuda_mutation_hotspot_parameters[gpu][row]);
        }
        cudaFree(cuda_A_0_mutation[gpu]);
        cudaFree(cuda_T_1_mutation[gpu]);
        cudaFree(cuda_G_2_mutation[gpu]);
        cudaFree(cuda_C_3_mutation[gpu]);

        cudaFree(cuda_mutation_hotspot_parameters[gpu]);

        for (int row = 0; row < num_Parent_sequence; row++)
        {
            cudaFree(cuda_parent_Sequences[gpu][row]);
            cudaFree(cuda_sequence_Configuration_standard[gpu][row]);
        }
        cudaFree(cuda_parent_Sequences[gpu]);
        cudaFree(cuda_sequence_Configuration_standard[gpu]);

        for (int row = 0; row < 2; row++)
        {
            cudaFree(cuda_parent_IDs[gpu][row]);
        }
        cudaFree(cuda_parent_IDs[gpu]);

        cudaFree(cuda_cell_Index[gpu]);

        for (int row = 0; row < num_Cells; row++)
        {
            cudaFree(cuda_totals_Progeny_Selectivity[gpu][row]);
        }
        cudaFree(cuda_totals_Progeny_Selectivity[gpu]);

        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaFree(cuda_progeny_Configuration[gpu][row]);
            cudaFree(cuda_progeny_Sequences[gpu][row]);
        }
        cudaFree(cuda_progeny_Configuration[gpu]);
        cudaFree(cuda_progeny_Sequences[gpu]);

        cudaFree(cuda_Dead_or_Alive[gpu]);

        cudaStreamDestroy(streams[gpu]);
    }
    cout << "Completed\n";

    // for (int test = 0; test < num_Progeny_being_Processed; test++)
    // {
    //     for (size_t i = 0; i < recombination_Hotspots + 1; i++)
    //     {
    //         cout << progeny_Configuration_Filled[test][i] << "\t";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

    // for (int test = 0; test < 1; test++)
    // {
    //     for (size_t i = 0; i < genome_Length; i++)
    //     {
    //         cout << progeny_Sequences[test][i];
    //     }
    //     cout << endl;
    // }
    // cout << endl;
    // for (int test = 0; test < num_Progeny_being_Processed; test++)
    // {
    //     cout << Dead_or_Alive[test] << endl;
    // }

    // test when more than one sequence in cell
    // get_Name() << "_" << tissue_Names[entry_array[tissue]] << "_" << current_Generation << "_" << to_string(indexes_Written[transfers])
    write_Full_Sequences_Progeny(functions,
                                 CPU_cores,
                                 tissue_Name,
                                 num_Progeny_being_Processed,
                                 genome_Length, recombination_Hotspots,
                                 sequence_Count,
                                 parent_IDs,
                                 progeny_Sequences, Dead_or_Alive, progeny_Configuration_Filled,
                                 source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation + 1),
                                 dead_List, sequence_Profiles, sequence_parent_Progeny_relationships,
                                 cells_of_progeny_location, start_Cell,
                                 max_sequences_per_File, index_Last_Written,
                                 tissue);
    // exit(-1);
}

__global__ void cuda_Progeny_Configurator(int num_Parents_to_Process, int start_Index,
                                          float **cuda_sequence_Configuration_standard, int recombination_Hotspots,
                                          int **cuda_progeny_Configuration, int *cuda_progeny_Stride, int remove_Back)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Parents_to_Process)
    {
        int parent_Index = tid + start_Index;
        int progeny_Fill_start = cuda_progeny_Stride[parent_Index] - remove_Back;
        int progeny_Fill_end = cuda_progeny_Stride[parent_Index + 1] - remove_Back;

        for (int progeny = 0; progeny < cuda_sequence_Configuration_standard[parent_Index][0]; progeny++)
        {
            cuda_progeny_Configuration[progeny_Fill_start + progeny][0] = parent_Index;

            for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
            {
                if (progeny < cuda_sequence_Configuration_standard[parent_Index][(hotspot * 2) + 2])
                {
                    cuda_progeny_Configuration[progeny_Fill_start + progeny][hotspot + 1] = parent_Index;
                }
                else
                {
                    cuda_progeny_Configuration[progeny_Fill_start + progeny][hotspot + 1] = -1;
                }
            }
        }

        for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
        {
            if (cuda_sequence_Configuration_standard[parent_Index][(hotspot * 2) + 2] > 0)
            {
                if (cuda_sequence_Configuration_standard[parent_Index][(hotspot * 2) + 2] != cuda_sequence_Configuration_standard[parent_Index][0])
                {
                    curandState state;
                    curand_init(clock64(), tid, 0, &state);
                    for (int i = progeny_Fill_start; i < progeny_Fill_end - 1; i++)
                    {

                        int j = curand(&state) % (progeny_Fill_end - i) + i;

                        int temp = cuda_progeny_Configuration[i][hotspot + 1];
                        cuda_progeny_Configuration[i][hotspot + 1] = cuda_progeny_Configuration[j][hotspot + 1];
                        cuda_progeny_Configuration[j][hotspot + 1] = temp;
                    }
                }
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

void node_within_host::progeny_Configurator(functions_library &functions,
                                            float **cuda_sequence_Configuration_standard, int recombination_Hotspots,
                                            int start_Index, int num_Parents_to_Process,
                                            int **progeny_Configuration, int *cuda_progeny_Stride, int progeny_Total, int remove_Back)
{
    cout << "Configuring " << num_Parents_to_Process << " parents' " << progeny_Total << " progeny\n";

    int **cuda_progeny_Configuration;
    //= functions.create_CUDA_2D_int(progeny_Total, 1 + recombination_Hotspots);
    cudaMallocManaged(&cuda_progeny_Configuration, progeny_Total * sizeof(int *));
    for (int row = 0; row < progeny_Total; row++)
    {
        cudaMalloc((void **)&(cuda_progeny_Configuration[row]), (1 + recombination_Hotspots) * sizeof(int));
    }

    cuda_Progeny_Configurator<<<functions.tot_Blocks_array[0], functions.tot_ThreadsperBlock_array[0]>>>(num_Parents_to_Process, start_Index,
                                                                                                         cuda_sequence_Configuration_standard, recombination_Hotspots,
                                                                                                         cuda_progeny_Configuration, cuda_progeny_Stride, remove_Back);
    cudaDeviceSynchronize();

    for (int row = 0; row < progeny_Total; row++)
    {
        cudaMemcpy(progeny_Configuration[row + remove_Back], cuda_progeny_Configuration[row], (recombination_Hotspots + 1) * sizeof(cuda_progeny_Configuration[0][0]), cudaMemcpyDeviceToHost);
    }

    functions.clear_Array_INT(cuda_progeny_Configuration, progeny_Total);
}

__device__ float generateExponential(curandState *state, float lambda)
{
    float u = curand_uniform(state);
    return -logf(u) / lambda;
}

__global__ void cuda_Parent_configuration(int num_Sequences, int **sequence_INT, int genome_Length, char *sites, float **cuda_sequence_Configuration_standard,
                                          float *cuda_Reference_fitness_survivability_proof_reading, int *cuda_num_effect_Segregating_sites,
                                          float **cuda_sequence_Fitness_changes, float **cuda_sequence_Proof_reading_changes,
                                          int recombination_Hotspots, float **cuda_recombination_hotspot_parameters,
                                          int *cuda_recombination_prob_Stride, float **cuda_recombination_Prob_matrix,
                                          int *cuda_recombination_select_Stride, float **cuda_recombination_Select_matrix,
                                          float *cuda_progeny_distribution_parameters_Array)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Sequences)
    {
        int site_Start = tid * genome_Length;
        int site_End = site_Start + genome_Length;

        int bp_Pos = 0;

        for (int site = site_Start; site < site_End; site++)
        {
            if (sites[site] == 'A' || sites[site] == 'a' || sites[site] == '0')
            {
                sequence_INT[tid][bp_Pos] = 0;
            }
            else if (sites[site] == 'T' || sites[site] == 't' || sites[site] == '1')
            {
                sequence_INT[tid][bp_Pos] = 1;
            }
            else if (sites[site] == 'G' || sites[site] == 'g' || sites[site] == '2')
            {
                sequence_INT[tid][bp_Pos] = 2;
            }
            else if (sites[site] == 'C' || sites[site] == 'c' || sites[site] == '3')
            {
                sequence_INT[tid][bp_Pos] = 3;
            }

            bp_Pos++;
        }

        // Fitness
        float fitness = cuda_Reference_fitness_survivability_proof_reading[0];
        for (int pos = 0; pos < cuda_num_effect_Segregating_sites[0]; pos++)
        {
            if (sequence_INT[tid][(int)cuda_sequence_Fitness_changes[pos][0] - 1] == 0)
            {
                fitness = fitness * cuda_sequence_Fitness_changes[pos][1];
            }
            else if (sequence_INT[tid][(int)cuda_sequence_Fitness_changes[pos][0] - 1] == 1)
            {
                fitness = fitness * cuda_sequence_Fitness_changes[pos][2];
            }
            else if (sequence_INT[tid][(int)cuda_sequence_Fitness_changes[pos][0] - 1] == 2)
            {
                fitness = fitness * cuda_sequence_Fitness_changes[pos][3];
            }
            else if (sequence_INT[tid][(int)cuda_sequence_Fitness_changes[pos][0] - 1] == 3)
            {
                fitness = fitness * cuda_sequence_Fitness_changes[pos][4];
            }
        }

        curandState localState;
        curand_init(clock64(), tid, 0, &localState);

        int progeny = 1;

        if (cuda_progeny_distribution_parameters_Array[0] == 0)
        {
            int failures = 0;
            int successes = 0;

            while (successes < cuda_progeny_distribution_parameters_Array[1])
            {
                float rand_num = curand_uniform(&localState);
                if (rand_num < cuda_progeny_distribution_parameters_Array[2])
                {
                    successes++;
                }
                else
                {
                    failures++;
                }
            }

            progeny = failures;
        }
        else if (cuda_progeny_distribution_parameters_Array[0] == 1)
        {
            // progeny = (int)rand_gamma_node(&localState, cuda_progeny_distribution_parameters_Array[1], cuda_progeny_distribution_parameters_Array[2]);

            float sum = 0.0f;
            for (int j = 0; j < cuda_progeny_distribution_parameters_Array[1]; ++j)
            {
                sum += generateExponential(&localState, 1.0f / cuda_progeny_distribution_parameters_Array[2]);
            }
            progeny = (int)sum;
        }
        else if (cuda_progeny_distribution_parameters_Array[0] == 2)
        {
            progeny = curand_poisson(&localState, cuda_progeny_distribution_parameters_Array[1]);
        }

        cuda_sequence_Configuration_standard[tid][0] = (int)(progeny * fitness);

        // proof reading
        if (cuda_Reference_fitness_survivability_proof_reading[2] != -1)
        {
            float proof_Reading = cuda_Reference_fitness_survivability_proof_reading[2];

            for (int pos = 0; pos < cuda_num_effect_Segregating_sites[2]; pos++)
            {
                if (sequence_INT[tid][(int)cuda_sequence_Proof_reading_changes[pos][0] - 1] == 0)
                {
                    proof_Reading = proof_Reading + cuda_sequence_Proof_reading_changes[pos][1];
                }
                else if (sequence_INT[tid][(int)cuda_sequence_Proof_reading_changes[pos][0] - 1] == 1)
                {
                    proof_Reading = proof_Reading + cuda_sequence_Proof_reading_changes[pos][2];
                }
                else if (sequence_INT[tid][(int)cuda_sequence_Proof_reading_changes[pos][0] - 1] == 2)
                {
                    proof_Reading = proof_Reading + cuda_sequence_Proof_reading_changes[pos][3];
                }
                else if (sequence_INT[tid][(int)cuda_sequence_Proof_reading_changes[pos][0] - 1] == 3)
                {
                    proof_Reading = proof_Reading + cuda_sequence_Proof_reading_changes[pos][4];
                }
            }
            if (proof_Reading > 1)
            {
                proof_Reading = 1;
            }
            else if (proof_Reading < 0)
            {
                proof_Reading = 0;
            }
            cuda_sequence_Configuration_standard[tid][1] = proof_Reading;
        }
        else
        {
            cuda_sequence_Configuration_standard[tid][1] = -1;
        }

        if (recombination_Hotspots > 0)
        {
            for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
            {
                int index_Progeny = (hotspot * 2) + 2;
                int index_Selectivity = index_Progeny + 1;

                float probability = cuda_recombination_hotspot_parameters[hotspot][2];
                float selectivity = cuda_recombination_hotspot_parameters[hotspot][3];

                for (int stride = cuda_recombination_prob_Stride[hotspot]; stride < cuda_recombination_prob_Stride[hotspot + 1]; stride++)
                {
                    if (sequence_INT[tid][(int)cuda_recombination_Prob_matrix[stride][0] - 1] == 0)
                    {
                        probability = probability + cuda_recombination_Prob_matrix[stride][1];
                    }
                    else if (sequence_INT[tid][(int)cuda_recombination_Prob_matrix[stride][0] - 1] == 1)
                    {
                        probability = probability + cuda_recombination_Prob_matrix[stride][2];
                    }
                    else if (sequence_INT[tid][(int)cuda_recombination_Prob_matrix[stride][0] - 1] == 2)
                    {
                        probability = probability + cuda_recombination_Prob_matrix[stride][3];
                    }
                    else if (sequence_INT[tid][(int)cuda_recombination_Prob_matrix[stride][0] - 1] == 3)
                    {
                        probability = probability + cuda_recombination_Prob_matrix[stride][4];
                    }
                }

                for (int stride = cuda_recombination_select_Stride[hotspot]; stride < cuda_recombination_select_Stride[hotspot + 1]; stride++)
                {
                    if (sequence_INT[tid][(int)cuda_recombination_Select_matrix[stride][0] - 1] == 0)
                    {
                        selectivity = selectivity * cuda_recombination_Select_matrix[stride][1];
                    }
                    else if (sequence_INT[tid][(int)cuda_recombination_Select_matrix[stride][0] - 1] == 1)
                    {
                        selectivity = selectivity * cuda_recombination_Select_matrix[stride][2];
                    }
                    else if (sequence_INT[tid][(int)cuda_recombination_Select_matrix[stride][0] - 1] == 2)
                    {
                        selectivity = selectivity * cuda_recombination_Select_matrix[stride][3];
                    }
                    else if (sequence_INT[tid][(int)cuda_recombination_Select_matrix[stride][0] - 1] == 3)
                    {
                        selectivity = selectivity * cuda_recombination_Select_matrix[stride][4];
                    }
                }

                if (probability > 1)
                {
                    probability = 1;
                }
                else if (probability < 0)
                {
                    probability = 0;
                }

                int hotspot_Progeny = 0;

                for (int trial = 0; trial < cuda_sequence_Configuration_standard[tid][0]; trial++)
                {
                    if (curand_uniform(&localState) < probability)
                    {
                        hotspot_Progeny++;
                    }
                }

                cuda_sequence_Configuration_standard[tid][index_Progeny] = hotspot_Progeny;

                cuda_sequence_Configuration_standard[tid][index_Selectivity] = selectivity;
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

void node_within_host::process_Sequences_get_Configuration(functions_library &functions,
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
                                                           int &start_Index)
{
    int num_of_Sequences_current = collected_Sequences.size();
    cout << "\nConfiguring multi gpu distribution of " << num_of_Sequences_current << " sequence(s)\n";

    int standard_num_per_GPU = num_of_Sequences_current / num_Cuda_devices;
    int remainder = num_of_Sequences_current % num_Cuda_devices;

    vector<pair<int, int>> start_stop_Per_GPU;

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        int start = gpu * standard_num_per_GPU;
        int stop = start + standard_num_per_GPU;

        start_stop_Per_GPU.push_back(make_pair(start, stop));
    }

    // cout<<"Test 1\n";

    start_stop_Per_GPU[num_Cuda_devices - 1].second = start_stop_Per_GPU[num_Cuda_devices - 1].second + remainder;

    string all_Sequences = "";

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        for (int sequence = start_stop_Per_GPU[gpu].first; sequence < start_stop_Per_GPU[gpu].second; sequence++)
        {
            all_Sequences.append(collected_Sequences[sequence]);
        }
    }

    // cout<<"Test 2\n";

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        for (int sequence = start_stop_Per_GPU[gpu].first; sequence < start_stop_Per_GPU[gpu].second; sequence++)
        {
            all_Sequences.append(collected_Sequences[sequence]);
        }
    }

    // cout<<"Test 3\n";

    char *full_Char;
    full_Char = (char *)malloc((all_Sequences.size() + 1) * sizeof(char));
    strcpy(full_Char, all_Sequences.c_str());

    // cout<<"Test 4\n";

    cudaStream_t streams[num_Cuda_devices];
    cudaDeviceProp deviceProp;

    char *cuda_full_Char[num_Cuda_devices];
    int **cuda_Sequence[num_Cuda_devices];
    float **cuda_sequence_Configuration_standard[num_Cuda_devices];

    float *cuda_Reference_fitness_survivability_proof_reading[num_Cuda_devices];
    int *cuda_mutation_recombination_proof_Reading_availability[num_Cuda_devices];
    int *cuda_num_effect_Segregating_sites[num_Cuda_devices];

    float **cuda_sequence_Fitness_changes[num_Cuda_devices];
    float **cuda_sequence_Proof_reading_changes[num_Cuda_devices];

    float **cuda_recombination_hotspot_parameters[num_Cuda_devices];
    int *cuda_tot_prob_selectivity[num_Cuda_devices];
    int *cuda_recombination_prob_Stride[num_Cuda_devices];
    int *cuda_recombination_select_Stride[num_Cuda_devices];
    float **cuda_recombination_Prob_matrix[num_Cuda_devices];
    float **cuda_recombination_Select_matrix[num_Cuda_devices];

    float *cuda_progeny_distribution_parameters_Array[num_Cuda_devices];

    // cout<<"Test 5\n";

    // cout << "Sites: " << num_effect_Segregating_sites[0] << endl;
    // for (int row = 0; row < num_effect_Segregating_sites[0]; row++)
    // {
    //     for (size_t i = 0; i < 5; i++)
    //     {
    //         cout << sequence_Fitness_changes[row][i] << "\t";
    //     }
    //     cout << endl;
    // }

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cudaGetDeviceProperties(&deviceProp, gpu);
        cout << "Intializing GPU " << CUDA_device_IDs[gpu] << "'s stream: " << deviceProp.name << endl;
        // cout << "Test 1\n";
        cudaMalloc(&cuda_full_Char[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * genome_Length * sizeof(char));
        cudaMemcpy(cuda_full_Char[gpu], full_Char + (start_stop_Per_GPU[gpu].first * genome_Length), (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * genome_Length * sizeof(char), cudaMemcpyHostToDevice);
        // cout << "Test 2\n";
        cudaMallocManaged(&cuda_Sequence[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(int *));
        cudaMallocManaged(&cuda_sequence_Configuration_standard[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(float *));
        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaMalloc((void **)&cuda_Sequence[gpu][row], genome_Length * sizeof(int));
            cudaMalloc((void **)&cuda_sequence_Configuration_standard[gpu][row], (3 + (2 * recombination_Hotspots)) * sizeof(float));
        }
        // cout << "Test 3\n";
        cudaMallocManaged(&cuda_Reference_fitness_survivability_proof_reading[gpu], 3 * sizeof(float));
        cudaMemcpy(cuda_Reference_fitness_survivability_proof_reading[gpu], Reference_fitness_survivability_proof_reading, 3 * sizeof(float), cudaMemcpyHostToDevice);
        // cout << "Test 4\n";
        cudaMallocManaged(&cuda_mutation_recombination_proof_Reading_availability[gpu], 3 * sizeof(int));
        cudaMemcpy(cuda_mutation_recombination_proof_Reading_availability[gpu], mutation_recombination_proof_Reading_availability, 3 * sizeof(int), cudaMemcpyHostToDevice);
        // cout << "Test 5\n";
        cudaMallocManaged(&cuda_num_effect_Segregating_sites[gpu], 3 * sizeof(int));
        cudaMemcpy(cuda_num_effect_Segregating_sites[gpu], num_effect_Segregating_sites, 3 * sizeof(int), cudaMemcpyHostToDevice);
        // cout << "Test 6\n";
        cudaMallocManaged(&cuda_progeny_distribution_parameters_Array[gpu], 3 * sizeof(float));
        cudaMemcpy(cuda_progeny_distribution_parameters_Array[gpu], progeny_distribution_parameters_Array, 3 * sizeof(float), cudaMemcpyHostToDevice);
        // cout << "Test 7\n";
        cudaMallocManaged(&cuda_sequence_Fitness_changes[gpu], num_effect_Segregating_sites[0] * sizeof(float *));
        // cout << "Test 7.1\n";
        for (int row = 0; row < num_effect_Segregating_sites[0]; row++)
        {
            cudaMalloc((void **)&cuda_sequence_Fitness_changes[gpu][row], 5 * sizeof(float));
            cudaMemcpy(cuda_sequence_Fitness_changes[gpu][row], sequence_Fitness_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
        }
        // cout << "Test 8\n";
        //  cudaMallocManaged(&cuda_sequence_Survivability_changes[gpu], num_effect_Segregating_sites[1] * sizeof(float *));
        //  for (int row = 0; row < num_effect_Segregating_sites[1]; row++)
        //  {
        //      cudaMalloc((void **)&cuda_sequence_Survivability_changes[gpu][row], 5 * sizeof(float));
        //      cudaMemcpy(cuda_sequence_Survivability_changes[gpu][row], sequence_Survivability_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
        //  }
        // cout << "Test 9\n";
        cudaMallocManaged(&cuda_sequence_Proof_reading_changes[gpu], num_effect_Segregating_sites[2] * sizeof(float *));
        for (int row = 0; row < num_effect_Segregating_sites[2]; row++)
        {
            cudaMalloc((void **)&cuda_sequence_Proof_reading_changes[gpu][row], 5 * sizeof(float));
            cudaMemcpy(cuda_sequence_Proof_reading_changes[gpu][row], sequence_Proof_reading_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
        }
        // cout << "Test 10\n";
        cudaMallocManaged(&cuda_tot_prob_selectivity[gpu], 2 * sizeof(int));
        cudaMallocManaged(&cuda_recombination_prob_Stride[gpu], (recombination_Hotspots + 1) * sizeof(int));
        cudaMallocManaged(&cuda_recombination_select_Stride[gpu], (recombination_Hotspots + 1) * sizeof(int));
        // cout << "Test 11\n";
        cudaMallocManaged(&cuda_recombination_hotspot_parameters[gpu], recombination_Hotspots * sizeof(float *));
        // cout << "Test 12\n";
        if (recombination_Hotspots > 0)
        {
            cudaMallocManaged(&cuda_recombination_Prob_matrix[gpu], tot_prob_selectivity[0] * sizeof(float *));
            cudaMallocManaged(&cuda_recombination_Select_matrix[gpu], tot_prob_selectivity[1] * sizeof(float *));

            cudaMemcpy(cuda_tot_prob_selectivity[gpu], tot_prob_selectivity, 2 * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(cuda_recombination_prob_Stride[gpu], recombination_prob_Stride, (recombination_Hotspots + 1) * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(cuda_recombination_select_Stride[gpu], recombination_select_Stride, (recombination_Hotspots + 1) * sizeof(int), cudaMemcpyHostToDevice);

            for (int row = 0; row < recombination_Hotspots; row++)
            {
                cudaMalloc((void **)&cuda_recombination_hotspot_parameters[gpu][row], 4 * sizeof(float));
                cudaMemcpy(cuda_recombination_hotspot_parameters[gpu][row], recombination_hotspot_parameters[row], 4 * sizeof(float), cudaMemcpyHostToDevice);
            }

            for (int row = 0; row < tot_prob_selectivity[0]; row++)
            {
                cudaMalloc((void **)&cuda_recombination_Prob_matrix[gpu][row], 5 * sizeof(float));
                cudaMemcpy(cuda_recombination_Prob_matrix[gpu][row], recombination_Prob_matrix[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
            }

            for (int row = 0; row < tot_prob_selectivity[1]; row++)
            {
                cudaMalloc((void **)&cuda_recombination_Select_matrix[gpu][row], 5 * sizeof(float));
                cudaMemcpy(cuda_recombination_Select_matrix[gpu][row], recombination_Select_matrix[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
            }
        }
        else
        {
            cudaMallocManaged(&cuda_recombination_Prob_matrix[gpu], 0 * sizeof(float *));
            cudaMallocManaged(&cuda_recombination_Select_matrix[gpu], 0 * sizeof(float *));
        }
        // cout << "Test 13\n";
        cudaStreamCreate(&streams[gpu]);
    }

    free(full_Char);

    cout << "Loaded " << num_of_Sequences_current << " sequence(s) and all pre-requisites to the GPU(s)\nInitiating GPU(s) execution\n";

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        // (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first, cuda_Sequence[gpu], genome_Length, cuda_full_Char[gpu]);

        cuda_Parent_configuration<<<functions.tot_Blocks_array[gpu], functions.tot_ThreadsperBlock_array[gpu], 0, streams[gpu]>>>(start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first, cuda_Sequence[gpu], genome_Length, cuda_full_Char[gpu], cuda_sequence_Configuration_standard[gpu],
                                                                                                                                  cuda_Reference_fitness_survivability_proof_reading[gpu], cuda_num_effect_Segregating_sites[gpu],
                                                                                                                                  cuda_sequence_Fitness_changes[gpu], cuda_sequence_Proof_reading_changes[gpu],
                                                                                                                                  recombination_Hotspots, cuda_recombination_hotspot_parameters[gpu],
                                                                                                                                  cuda_recombination_prob_Stride[gpu], cuda_recombination_Prob_matrix[gpu],
                                                                                                                                  cuda_recombination_select_Stride[gpu], cuda_recombination_Select_matrix[gpu],
                                                                                                                                  cuda_progeny_distribution_parameters_Array[gpu]);
    }

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cudaStreamSynchronize(streams[gpu]);
    }

    cout << "GPU(s) streams completed and synchronized\nCopying data from GPU to Host memory\n";

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);

        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaMemcpy(parent_Sequences[start_stop_Per_GPU[gpu].first + row + start_Index], cuda_Sequence[gpu][row], genome_Length * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(sequence_Configuration_standard[start_stop_Per_GPU[gpu].first + row + start_Index], cuda_sequence_Configuration_standard[gpu][row], (2 + (2 * recombination_Hotspots)) * sizeof(float), cudaMemcpyDeviceToHost);
        }
    }
    cout << "Data received by host\n";

    cout << "Terminating GPU streams: ";
    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);

        cudaFree(cuda_full_Char[gpu]);

        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaFree(cuda_Sequence[gpu][row]);
            cudaFree(cuda_sequence_Configuration_standard[gpu][row]);
        }
        cudaFree(cuda_Sequence[gpu]);
        cudaFree(cuda_sequence_Configuration_standard[gpu]);

        cudaFree(cuda_Reference_fitness_survivability_proof_reading[gpu]);
        cudaFree(cuda_mutation_recombination_proof_Reading_availability[gpu]);
        cudaFree(cuda_num_effect_Segregating_sites[gpu]);

        for (int row = 0; row < num_effect_Segregating_sites[0]; row++)
        {
            cudaFree(cuda_sequence_Fitness_changes[gpu][row]);
        }
        cudaFree(cuda_sequence_Fitness_changes[gpu]);

        for (int row = 0; row < num_effect_Segregating_sites[2]; row++)
        {
            cudaFree(cuda_sequence_Proof_reading_changes[gpu][row]);
        }
        cudaFree(cuda_sequence_Proof_reading_changes[gpu]);

        cudaFree(cuda_tot_prob_selectivity[gpu]);
        cudaFree(cuda_recombination_prob_Stride[gpu]);
        cudaFree(cuda_recombination_select_Stride[gpu]);
        if (recombination_Hotspots > 0)
        {
            for (int row = 0; row < recombination_Hotspots; row++)
            {
                cudaFree(cuda_recombination_hotspot_parameters[gpu][row]);
            }

            for (int row = 0; row < tot_prob_selectivity[0]; row++)
            {
                cudaFree(cuda_recombination_Prob_matrix[gpu][row]);
            }

            for (int row = 0; row < tot_prob_selectivity[1]; row++)
            {
                cudaFree(cuda_recombination_Select_matrix[gpu][row]);
            }
        }
        cudaFree(cuda_recombination_hotspot_parameters[gpu]);
        cudaFree(cuda_recombination_Prob_matrix[gpu]);
        cudaFree(cuda_recombination_Select_matrix[gpu]);

        cudaFree(cuda_progeny_distribution_parameters_Array[gpu]);

        cudaStreamDestroy(streams[gpu]);
    }
    cout << "Completed\n";
}

vector<int> node_within_host::assign_Cells(functions_library &functions, int **parents_in_Tissue, int num_Viral_particles, int &tissue,
                                           int distribution_Type, float &parameter_1, float &parameter_2,
                                           set<int> &check_to_Remove,
                                           int &gen_Phase, float &variable_1, float &variable_2,
                                           mt19937 &gen)
{
    /**
     * *Function handling the infection of of cells by the viral particles.
     **/

    cout << "\nAssigning cell(s) their virulant particle(s)\n";

    // if (parents_Prev_generation != 0)
    // {
    //     num_Viral_particles = parents_Prev_generation;
    // }

    vector<int> start_Stop_cells;

    int cells_Assigned = 0;
    int particles_Assigned = 0;

    start_Stop_cells.push_back(0);

    // int cells_Full = 0;

    int index_Track_removed = 0;
    int particle_ID = 0;

    // test
    // cell_Limit[tissue] = 2;

    vector<int> removals(check_to_Remove.begin(), check_to_Remove.end());
    /**
     * First all viral particles are assigned to a cell.
     **/
    do
    {
        // if (cell_Limit[tissue] != -1 && cells_Assigned >= cell_Limit[tissue])
        // {
        //     cells_Full = 1;
        //     break;
        // }

        int num_Particles_in_Cell;

        if (distribution_Type == 0)
        {
            binomial_distribution<int> num_Particles(parameter_1, parameter_2);
            num_Particles_in_Cell = num_Particles(gen);
        }
        else
        {
            gamma_distribution<float> num_Particles(parameter_1, parameter_2);
            num_Particles_in_Cell = (int)num_Particles(gen);
        }

        if (num_Particles_in_Cell > 0)
        {
            for (int cell = 0; cell < num_Particles_in_Cell; cell++)
            {
                parents_in_Tissue[1][particles_Assigned] = cells_Assigned;
                //// Account for dead
                if (index_Track_removed < removals.size())
                {
                    while (particle_ID == removals[index_Track_removed])
                    {
                        particle_ID++;
                        index_Track_removed++;
                    }
                }

                parents_in_Tissue[0][particles_Assigned] = particle_ID;

                particle_ID++;
                particles_Assigned++;

                if (particles_Assigned >= num_Viral_particles)
                {
                    break;
                }
            }

            start_Stop_cells.push_back(particles_Assigned);
            cells_Assigned++;
        }

    } while (particles_Assigned < num_Viral_particles);

    // cout << cells_Assigned - 1 << endl;
    srand(time(0));
    random_shuffle(parents_in_Tissue[0], parents_in_Tissue[0] + num_Viral_particles);

    // for (int i = 0; i < num_Viral_particles; i++)
    // {
    //     cout << parents_in_Tissue[0][i] << endl;
    // }

    if (gen_Phase == 1 || gen_Phase == 2)
    {
        int new_Parent_Count = -1;

        if (current_Generation == 0)
        {
            parents_Prev_generation[tissue] = num_Viral_particles;
        }

        if (gen_Phase == 1)
        {
            /**
             * If stationary population is maintained as is, equal to previous generations parent population.
             * *It can vary around a given variance.
             **/
            cout << "Tissue phase: Stationary phase\n";
            if (num_Viral_particles >= parents_Prev_generation[tissue])
            {
                normal_distribution<float> distribution(parents_Prev_generation[tissue], variable_1);
                new_Parent_Count = distribution(gen);
                // cout << "parents_Prev_generation: " << parents_Prev_generation << endl;
                // cout << new_Parent_Count << endl;
                if (new_Parent_Count < num_Viral_particles && new_Parent_Count >= 0)
                {
                    cout << "Parent population maintained at: " << new_Parent_Count << endl;
                }
                else
                {
                    new_Parent_Count = -1;
                }
            }
        }
        else if (gen_Phase == 2)
        {
            /**
             * If in depriciation the population is reduced relative to the previous generation.
             **/
            cout << "Tissue phase: Depriciation phase\n";
            if (num_Viral_particles >= parents_Prev_generation[tissue])
            {
                new_Parent_Count = functions.beta_Distribution(variable_1, variable_2, gen) * parents_Prev_generation[tissue];
                new_Parent_Count = parents_Prev_generation[tissue] - new_Parent_Count;
                cout << "Parent population reduced to: " << new_Parent_Count << endl;
            }
        }
        if (new_Parent_Count != -1)
        {
            // if (new_Parent_Count > 0)
            // {
            //     int **temp = functions.create_INT_2D_arrays(2, new_Parent_Count);

            //     for (int parent = 0; parent < new_Parent_Count; parent++)
            //     {
            //         temp[0][parent] = parents_in_Tissue[0][parent];
            //         temp[1][parent] = parents_in_Tissue[1][parent];
            //     }

            //     functions.clear_Array_int_CPU(parents_in_Tissue, 2);

            //     parents_in_Tissue = functions.create_INT_2D_arrays(2, new_Parent_Count);

            //     for (int parent = 0; parent < new_Parent_Count; parent++)
            //     {
            //         parents_in_Tissue[0][parent] = temp[0][parent];
            //         parents_in_Tissue[1][parent] = temp[1][parent];
            //     }

            //     functions.clear_Array_int_CPU(temp, 2);

            //     for (int i = 0; i < new_Parent_Count; i++)
            //     {
            //         cout << parents_in_Tissue[0][i] << endl;
            //     }
            // }

            cout << "Resizing parent cell array\n";

            /**
             * The new population is controlled by the stride array's accesing of the parents_in_Tisse array.
             **/

            vector<int> temp_Cells;
            if (new_Parent_Count > 0)
            {
                cout << "Configuring cell index\n";
                int max_Cell = parents_in_Tissue[1][new_Parent_Count - 1] + 1;
                for (int cell = 0; cell < max_Cell; cell++)
                {
                    temp_Cells.push_back(start_Stop_cells[cell]);
                }
                temp_Cells.push_back(start_Stop_cells[max_Cell]);

                if (temp_Cells[temp_Cells.size() - 1] > new_Parent_Count)
                {
                    temp_Cells[temp_Cells.size() - 1] = new_Parent_Count;
                }
            }
            else
            {
                temp_Cells.push_back(0);
            }

            start_Stop_cells.clear();
            start_Stop_cells = temp_Cells;

            parents_Prev_generation[tissue] = new_Parent_Count;
        }
        else
        {
            parents_Prev_generation[tissue] = num_Viral_particles;
        }
    }
    else
    {
        cout << "Tissue phase: Neutral phase, all particles are viable\n";
        parents_Prev_generation[tissue] = num_Viral_particles;
    }

    /**
     * If the occupied cell count exceeds the set cell limit for the tissue it has to be resized.
     **/

    if ((start_Stop_cells.size() - 1) > cell_Limit[tissue])
    {
        vector<int> temp;
        for (int cell = 0; cell < cell_Limit[tissue] + 1; cell++)
        {
            temp.push_back(start_Stop_cells[cell]);
        }
        start_Stop_cells.clear();
        start_Stop_cells = temp;
        parents_Prev_generation[tissue] = start_Stop_cells[start_Stop_cells.size() - 1];
    }

    return start_Stop_cells;
}

void node_within_host::intialize_Tissues(string &host_Folder, vector<vector<string>> &tissue_Sequences, functions_library &functions)
{
    current_Viral_load_per_Tissue = (int *)malloc(sizeof(int) * num_Tissues);
    dead_Particle_count = (int *)malloc(sizeof(int) * num_Tissues);
    parents_Prev_generation = (int *)malloc(sizeof(int) * num_Tissues);
    current_Generation = 0;

    set<int> init_removed_by_Transfer_Indexes;

    for (int tissue = 0; tissue < num_Tissues; tissue++)
    {
        vector<string> tissue_Sequence;
        functions.config_Folder(host_Folder + "/" + to_string(tissue), to_string(cave_ID) + "_" + to_string(host_ID) + " Tissue " + to_string(tissue + 1));

        tissue_Sequences.push_back(tissue_Sequence);
        current_Viral_load_per_Tissue[tissue] = 0;
        dead_Particle_count[tissue] = 0;

        parents_Prev_generation[tissue] = 0;

        removed_by_Transfer_Indexes.push_back(init_removed_by_Transfer_Indexes);
    }
}

void node_within_host::clear_Arrays_end()
{
    cout << "\nClearing all arrays for host " << get_Name() << endl;
    free(current_Viral_load_per_Tissue);
    free(dead_Particle_count);
    free(parents_Prev_generation);

    removed_by_Transfer_Indexes.clear();
}

int node_within_host::get_Load(int &num_tissues_Calc, int *tissue_array)
{
    int sum = 0;

    for (int tissue = 0; tissue < num_tissues_Calc; tissue++)
    {
        sum = sum + (current_Viral_load_per_Tissue[tissue_array[tissue]] - dead_Particle_count[tissue_array[tissue]] - removed_by_Transfer_Indexes[tissue_array[tissue]].size());
    }

    return sum;
}

int node_within_host::infectious_status(int &num_tissues, int *tissue_array)
{
    if (get_Load(num_tissues, tissue_array) >= infectious_Load)
    {
        set_Infectious();
        return 1;
    }
    else
    {
        return 0;
    }
}
int node_within_host::terminal_status(int &num_tissues, int *tissue_array)
{
    if (get_Load(num_tissues, tissue_array) >= terminal_Load)
    {
        set_Dead();
        // if (enable_Folder_management == "YES")
        // {
        //     compress_Folder(source_sequence_Data_folder, enable_Compression);
        // }
        return 1;
    }
    else
    {
        return 0;
    }
}

int node_within_host::terminal_status(int &num_tissues, int *tissue_array, string source_sequence_Data_folder,
                                      string enable_Folder_management, string enable_Compression)
{
    if (get_Load(num_tissues, tissue_array) >= terminal_Load)
    {
        set_Dead();
        if (enable_Folder_management == "YES")
        {
            compress_Folder(source_sequence_Data_folder, enable_Compression);
        }
        return 1;
    }
    else
    {
        return 0;
    }
}

string node_within_host::get_Name()
{
    return to_string(cave_ID) + "_" + to_string(host_ID);
}

string node_within_host::get_Status()
{
    return this->status;
}

int node_within_host::get_Profile()
{
    return profile_ID;
}

int node_within_host::get_host_Index()
{
    return host_Index;
}

int node_within_host::get_Generation()
{
    return current_Generation;
}

int *node_within_host::get_current_Viral_load_per_Tissue()
{
    return current_Viral_load_per_Tissue;
}

float node_within_host::get_infection_probability()
{
    return infection_probability;
}

void node_within_host::set_Infected()
{
    this->status = "Infected";
}
void node_within_host::set_Infectious()
{
    this->status = "Infectious";
}
void node_within_host::set_Removed()
{
    this->status = "Removed";
    // // ! Compress the hosts intermediate folder
}
void node_within_host::set_Dead()
{
    this->status = "Dead";
    clear_Arrays_end();
    // // ! Compress the hosts intermediate folder
}

void node_within_host::set_Susceptible()
{
    this->status = "Susceptible";
    // // ! Compress the hosts intermediate folder
}