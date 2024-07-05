#include "simulator_Master.cuh"

simulator_Master::simulator_Master(string parameter_Master_Location)
{
    cout << "Intializing Apollo\n";

    parameter_load Parameters = parameter_load();
    functions_library function = functions_library();

    vector<string> parameters_List = {
        "\"CUDA Device IDs\"",
        "\"CPU cores\"",
        "\"GPU max units\"",
        "\"Intermediate folders\"",
        "\"Output folders\"",
        "\"Multi read\"",
        "\"Network profile\"",
        "\"Nodes master profile\"",
        "\"Sequence master profile\"",
        "\"Intermediate Sequences per file\"",
        "\"Process cell rate\"",
        "\"Start date\"",
        "\"Stop after generations\"",
        "\"Enable folder management\"",
        "\"First infection\""};

    vector<string> found_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List);

    cout << "\nConfiguring folders:\n";

    output_Folder_location = Parameters.get_STRING(found_Parameters[4]);
    intermediate_Folder_location = Parameters.get_STRING(found_Parameters[3]);
    max_sequences_per_File = Parameters.get_INT(found_Parameters[9]);
    max_Cells_at_a_time = Parameters.get_INT(found_Parameters[10]);
    start_Date = Parameters.get_STRING(found_Parameters[11]);
    stop_after_generations = function.to_Upper_Case(Parameters.get_STRING(found_Parameters[12]));
    first_Infection = function.to_Upper_Case(Parameters.get_STRING(found_Parameters[14]));

    if (max_sequences_per_File <= 0)
    {
        cout << "ERROR: Intermediate Sequences per file PARAMETER MUST BE GREATER THAN ZERO.\n";
        exit(-1);
    }

    if (max_Cells_at_a_time <= 0)
    {
        cout << "ERROR: Process cell rate PARAMETER MUST BE GREATER THAN ZERO.\n";
        exit(-1);
    }

    function.config_Folder(intermediate_Folder_location, "Intermediate");
    cout << "\nFolder management: ";
    if (function.to_Upper_Case(Parameters.get_STRING(found_Parameters[13])) == "YES")
    {
        enable_Folder_management = "YES";

        vector<string> folder_Management = {"\"Compress folders\""};
        vector<string> folder_management_Parameters = Parameters.get_parameters(parameter_Master_Location, folder_Management);

        if (function.to_Upper_Case(Parameters.get_STRING(folder_Management[0])) == "YES")
        {
            enable_Compression = "YES";
        }
    }
    cout << enable_Folder_management << endl;
    cout << "Folder compression: " << enable_Compression << endl
         << endl;
    // exit(-1);

    function.config_Folder(output_Folder_location, "Output");

    output_Network_location = this->output_Folder_location + "/network_Data";
    function.config_Folder(output_Network_location, "Network");
    output_Node_location = this->output_Folder_location + "/node_Data";
    function.config_Folder(output_Node_location, "Node");
    network_File_location = output_Network_location + "/node_node_Relationships.csv";
    function.create_File(network_File_location, "Source\tTarget");
    Host_source_target_network_location = output_Network_location + "/hosts_source_target_Relationships.csv";
    function.create_File(Host_source_target_network_location, "Source\tTarget\tTime_Infection\tDate_Infection");

    cout << "\nConfiguring node master profiles:\n";
    node_Master_location = Parameters.get_STRING(found_Parameters[7]);

    cout << "\nConfiguring sequence master profiles:\n";
    sequence_Master_location = Parameters.get_STRING(found_Parameters[8]);

    if (stop_after_generations == "YES")
    {
        cout << "\nConfiguring simulation termination by overall generations run\n";

        vector<string> parameters_List_stop_Gen_Mode = {"\"Mode to stop\""};
        vector<string> stop_Gen_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List_stop_Gen_Mode);

        if (function.to_Upper_Case(Parameters.get_STRING(stop_Gen_Parameters[0])) == "GENERATIONS")
        {
            stop_gen_Mode = 0;
            cout << "Stop after given number of generations: ";
            parameters_List_stop_Gen_Mode.clear();
            stop_Gen_Parameters.clear();

            vector<string> parameters_List_stop_Gen_Mode = {"\"Number of generations\""};
            vector<string> stop_Gen_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List_stop_Gen_Mode);

            stop_generations_Count = Parameters.get_INT(stop_Gen_Parameters[0]);
            cout << stop_generations_Count << endl;
        }
        else
        {
            cout << "Stop after given date: ";
            stop_gen_Mode = 1;
            parameters_List_stop_Gen_Mode.clear();
            stop_Gen_Parameters.clear();

            vector<string> parameters_List_stop_Gen_Mode = {"\"End date\""};
            vector<string> stop_Gen_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List_stop_Gen_Mode);

            string stop_Date_String = Parameters.get_STRING(stop_Gen_Parameters[0]);

            cout << stop_Date_String << endl;
            vector<string> split_Date;
            function.split(split_Date, stop_Date_String, '-');

            stop_Date = function.date_to_Decimal(stoi(split_Date[0]), stoi(split_Date[1]), stoi(split_Date[2]));
            cout << "Decimal date: " << stop_Date << endl;
        }
    }
    else
    {
        cout << "\nSimulation termination by overall generations run is deactivated\n";
    }

    cout << "\nConfiguring hardware resources:\n\n";
    this->CPU_cores = Parameters.get_INT(found_Parameters[1]);
    cout << "Available CPU cores: " << this->CPU_cores << endl
         << endl;

    cout << "Maximum cells processed at a time: " << this->max_Cells_at_a_time << endl
         << endl;

    cout << "Maximum sequences per file (Intermediary): " << this->max_sequences_per_File << endl
         << endl;
    this->multi_Read = Parameters.get_STRING(found_Parameters[5]);
    transform(multi_Read.begin(), multi_Read.end(), multi_Read.begin(), ::toupper);
    cout << "Multiple read and write: " << this->multi_Read << endl
         << endl;

    // this->CUDA_device_number = Parameters.get_INT(found_Parameters[0]);
    string cuda_IDs_String = Parameters.get_STRING(found_Parameters[0]);
    vector<string> cuda_IDs;
    function.split(cuda_IDs, cuda_IDs_String, ',');
    if (cuda_IDs.size() > 0)
    {
        this->num_Cuda_devices = cuda_IDs.size();
        CUDA_device_IDs = (int *)malloc(sizeof(int) * num_Cuda_devices);
        tot_Blocks = (int *)malloc(sizeof(int) * num_Cuda_devices);
        tot_ThreadsperBlock = (int *)malloc(sizeof(int) * num_Cuda_devices);
        function.print_Cuda_devices(cuda_IDs, this->CUDA_device_IDs, num_Cuda_devices, this->tot_Blocks, this->tot_ThreadsperBlock);
    }
    else
    {
        cout << "ERROR: THERE HAS TO BE AT LEAST ONE CUDA DEVICE SELECTED\n";
        exit(-1);
    }
    // function.print_Cuda_device(this->CUDA_device_number, this->tot_Blocks, this->tot_ThreadsperBlock);

    this->gpu_Limit = Parameters.get_INT(found_Parameters[2]);

    cout << "Per round GPU max unit: " << this->gpu_Limit << endl
         << endl;

    // exit(-1);

    configure_Network_Profile(Parameters.get_STRING(found_Parameters[6]), Parameters);
    cout << "\n";

    // exit(-1);
}

void simulator_Master::configure_Network_Profile(string network_Profile_File, parameter_load &Parameters)
{
    cout << "Configuring network profile: " << network_Profile_File << endl;

    vector<string> parameters_List = {"\"Network type\""};
    vector<string> found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

    transform(found_Parameters[0].begin(), found_Parameters[0].end(), found_Parameters[0].begin(), ::toupper);

    // parameters_List.clear();
    // found_Parameters.clear();

    if (Parameters.get_STRING(found_Parameters[0]) == "BA MODEL")
    {
        cout << "\nBarabsi Albert model selected: \n";
        network_Model = "BA";

        parameters_List = {"\"BA model number of nodes\"",
                           "\"BA model standard new connections\""};
        found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

        Total_number_of_Nodes = Parameters.get_INT(found_Parameters[0]);
        cout << "Number of nodes: " << Total_number_of_Nodes << endl;
        connection_Model = Parameters.get_STRING(found_Parameters[1]);
        transform(connection_Model.begin(), connection_Model.end(), connection_Model.begin(), ::toupper);

        cout << "Node connection type: " << connection_Model << endl;

        if (connection_Model == "FIXED")
        {
            parameters_List = {"\"BA model fixed new connections\""};
            found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

            BA_FIXED = Parameters.get_INT(found_Parameters[0]);
            cout << "Fixed new connections: " << BA_FIXED << endl;
        }
        else if (connection_Model == "NEGATIVE BINOMIAL")
        {
            parameters_List = {"\"BA model Negative binomial sucesses\"",
                               "\"BA model Negative binomial probability\""};
            found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

            BA_NB_sucesses = Parameters.get_INT(found_Parameters[0]);
            BA_NB_probability = Parameters.get_FLOAT(found_Parameters[1]);

            cout << "Negative Binomial sucesses: " << BA_NB_sucesses << endl;
            cout << "Negative Binomial probability: " << BA_NB_probability << endl;
        }
        else if (connection_Model == "POISSON")
        {
            parameters_List = {"\"BA model Poisson mean\""};
            found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

            BA_Poisson_mean = Parameters.get_FLOAT(found_Parameters[0]);

            cout << "Poisson mean: " << BA_Poisson_mean << endl;
        }
    }
    else if (Parameters.get_STRING(found_Parameters[0]) == "SC MODEL")
    {
        cout << "\nStandard Caveman model selected: \n";
        network_Model = "SCM";

        parameters_List = {"\"SC_model number of caves\"",
                           "\"SC_model number of nodes per caves\""};
        found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

        SCM_number_of_caves = Parameters.get_INT(found_Parameters[0]);
        SCM_number_of_nodes_per_cave = Parameters.get_INT(found_Parameters[1]);

        cout << "Number of caves: " << SCM_number_of_caves;
        cout << "\nNumber of nodes per cave: " << SCM_number_of_nodes_per_cave;
        Total_number_of_Nodes = SCM_number_of_caves * SCM_number_of_nodes_per_cave;
        cout << "\nTotal nodes in network: " << Total_number_of_Nodes << endl;
    }
    else if (Parameters.get_STRING(found_Parameters[0]) == "DC MODEL")
    {
        cout << "\nDynamic Caveman model selected: \n";
        network_Model = "DCM";

        parameters_List = {"\"DC_model number of caves\"",
                           "\"DC_model node distribution\"",
                           "\"DC_model neighbouring nodes percentage\"",
                           "\"DC_model global nodes percentage\""};

        found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

        DCM_number_of_caves = Parameters.get_INT(found_Parameters[0]);
        cout << "Number of caves: " << DCM_number_of_caves << endl;

        connection_Model = Parameters.get_STRING(found_Parameters[1]);
        transform(connection_Model.begin(), connection_Model.end(), connection_Model.begin(), ::toupper);

        DC_percent_Neighbouring = Parameters.get_FLOAT(found_Parameters[2]);
        DC_percent_Global_freedom = Parameters.get_FLOAT(found_Parameters[3]);

        cout << "Node connection type: " << connection_Model << endl;

        if (connection_Model == "NEGATIVE BINOMIAL")
        {
            parameters_List = {"\"DC_model node Negative binomial sucesses\"",
                               "\"DC_model node Negative binomial probability\""};

            found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

            DC_ND_sucesses = Parameters.get_INT(found_Parameters[0]);
            DC_ND_probability = Parameters.get_FLOAT(found_Parameters[1]);

            cout << "Negative Binomial sucesses: " << DC_ND_sucesses << endl;
            cout << "Negative Binomial probability: " << DC_ND_probability << endl;
        }
        else if (connection_Model == "POISSON")
        {
            parameters_List = {"\"DC_model node Poisson mean\""};
            found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

            DC_Poisson_mean = Parameters.get_FLOAT(found_Parameters[0]);

            cout << "Poisson mean: " << DC_Poisson_mean << endl;
        }
    }
    else if (Parameters.get_STRING(found_Parameters[0]) == "RANDOM MODEL")
    {
        cout << "\nRandom model selected: \n";
        network_Model = "RANDOM";

        parameters_List = {"\"Random model number of nodes\""};
        found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

        Total_number_of_Nodes = Parameters.get_INT(found_Parameters[0]);
    }
    else if (Parameters.get_STRING(found_Parameters[0]) == "ER MODEL")
    {
        cout << "\nErdos and Renyi Random model selected: \n";
        network_Model = "ER_RANDOM";

        parameters_List = {"\"ER Random model total number of nodes\"",
                           "\"ER Random model probability of linkage\""};

        found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

        Total_number_of_Nodes = Parameters.get_INT(found_Parameters[0]);
        ER_link_probability = Parameters.get_FLOAT(found_Parameters[1]);

        // exit(-1);
    }
    else
    {
        cout << "ERROR Incorrect network selected. Please check \"Network type\" in the network parameter file: \"" << network_Profile_File << "\"";
        exit(-1);
    }
}

void simulator_Master::ingress()
{
    functions_library functions = functions_library(tot_Blocks, tot_ThreadsperBlock, CUDA_device_IDs, num_Cuda_devices, gpu_Limit, CPU_cores);

    cout
        << "STEP 1: Configuring population network\n\n";

    //// Compatible for both BA and Caveman
    // INT, INT = Cave_ID and Node, for BA CaveID is 0 for all.
    // vector<vector<pair<int, int>>> each_Nodes_Connection;
    //// Only for CAVEMAN models, to keep track of each nodes Caves

    network_Manager(functions);

    // exit(-1);

    cout << "STEP 2: Configuring node profiles and within host mechanics\n\n";
    node_Master_Manager(functions);

    // exit(-1);

    cout << "STEP 3: Configuring sequence profiles\n\n";
    sequence_Master_Manager(functions);

    cout << "STEP 4: Assigning profiles to network nodes\n\n";
    vector<node_within_host> Hosts = node_Profile_assignment_Manager(functions);

    // for (int host = 0; host < Total_number_of_Nodes; host++)
    // {
    //     Hosts[host].print_All();
    // }

    // exit(-1);

    cout << "STEP 5: Infection begins\n\n";
    apollo(functions, Hosts);
}

void simulator_Master::thread_compress_Folders(int new_dead_Host, vector<node_within_host> &Hosts)
{
    //  for (int host = start; host < stop; host++)
    //{
    Hosts[new_dead_Host].compress_Folder(intermediary_Sequence_location + "/" + to_string(Hosts[new_dead_Host].get_host_Index()), enable_Compression, 1);
    //}
}

void simulator_Master::apollo(functions_library &functions, vector<node_within_host> &Hosts)
{
    // // TODO: host infection times
    // // TODO: create an overall generational summary
    // // TODO: Per individual generational summary
    // // TODO: SIRS
    // // TODO: Terminate after given number of generations or time
    // // TODO: CHECK IF ALL ARRAYS ARE CLEARED

    // Source\tTarget\tInfection_time
    // Generation_Network\tGeneration_time_decimal\tDate\tSusceptible_population\tInfected\tInfectious\tDead\tRemoved
    // Generation\tGeneration_time_decimal\tDate\tTissue\tParents\tParticles_produced\tDead\tTransferred\tRemaining

    cout << "Configuring infection start date: " << start_Date << endl;

    vector<string> split_Date;
    functions.split(split_Date, start_Date, '-');

    float decimal_Date = functions.date_to_Decimal(stoi(split_Date[0]), stoi(split_Date[1]), stoi(split_Date[2]));
    float date_Increment = generation_Time / (float)365.25;

    cout << "Decimal date: " << decimal_Date << endl;
    cout << "Date incerement by generation: " << date_Increment << endl;
    // cout << "Decimal date: " << (decimal_Date + date_Increment) << endl;
    // int year, month, day;
    // functions.decimal_to_Date(decimal_Date + date_Increment, year, month, day);
    // cout << "Conversion check:\n"
    //      << year << endl
    //      << month << endl
    //      << day << endl;

    // exit(-1);

    cout << "\nConfiguring susceptible population\n\n";

    vector<int> susceptible_Population;

    for (size_t i = 0; i < Total_number_of_Nodes; i++)
    {
        susceptible_Population.push_back(i);
        // removed_Population.push_back(-1);
    }

    cout << "Starting infection\n\n";

    vector<int> infected_Population;

    int first_Infected = get_first_Infected(susceptible_Population, infected_Population, functions);

    Hosts[first_Infected].begin_Infection(functions, intermediary_Sequence_location, entry_tissues, entry_array, max_sequences_per_File, output_Node_location, tissue_Names, first_Infection);

    cout << endl;
    functions.folder_Delete(intermediate_Folder_location + "/sequence_Data/reference_Sequences");

    // cout << Total_number_of_Nodes << endl;
    // cout << susceptible_Population.size() << endl;
    // cout << infected_Population.size() << endl;
    // cout << to_string(susceptible_Population.size() + infected_Population.size()) << endl;

    int stop = 0;

    // test
    // vector<pair<int, int>> host_Connections_TEST;
    // host_Connections_TEST.push_back(make_pair(0, 88));
    // host_Connections_TEST.push_back(make_pair(0, 77));
    // host_Connections_TEST.push_back(make_pair(0, 1));
    // Node_search(host_Connections_TEST);
    // exit(-1);

    random_device rd;
    mt19937 gen(rd());

    int overall_Generations = 0;
    int removed_Count = 0;
    int dead_Count = 0;

    string overall_Generational_Summary = output_Network_location + "/overall_generational_Summary.csv";
    functions.create_File(overall_Generational_Summary, "overall_Generation\tdecimal_Date\tDate\tsusceptible_Population\tinfected_Population\tinfectious_Population\tremoved_Population\tdead_Population");

    // for (int test = 0; test < all_node_IDs.size(); test++)
    // {
    //     for (size_t i = 0; i < each_Nodes_Connection_INT[test].size(); i++)
    //     {
    //         cout << all_node_IDs[each_Nodes_Connection_INT[test][i]].first << "_" << all_node_IDs[each_Nodes_Connection_INT[test][i]].second << ", ";
    //     }
    //     cout << endl;
    // }

    cout << "\nWriting node index information\n";
    functions.create_File(intermediary_Index_location + "/node_Index.csv", "node_Index\tnode_ID");
    fstream index_Data_file;
    index_Data_file.open(intermediary_Index_location + "/node_Index.csv", ios::app);
    if (index_Data_file.is_open())
    {
        for (int node = 0; node < all_node_IDs.size(); node++)
        {
            index_Data_file << node << "\t" << all_node_IDs[node].first << "_" << all_node_IDs[node].second << "\n";
        }
        index_Data_file.close();
    }
    else
    {
        cout << "\nERROR: INDEX DATA FILE CANNOT BE OPENED: " << intermediary_Index_location << "/node_Index.csv\n";
    }

    // exit(-1);

    cout << "\nCreating processing time file\n";
    fstream time_Track;
    functions.create_File(output_Folder_location + "/run_Time_track.csv", "Generation\tTime");
    time_Track.open(output_Folder_location + "/run_Time_track.csv", ios::app);

    // exit(-1);

    do
    {
        auto start_Time = chrono::high_resolution_clock::now();
        // check dead and remove from infected

        /**
         * Structuring the epidemiological compartments
         **/

        cout << "\nDetermining removed/ dead and newly infectious hosts\n";
        vector<int> temp;
        vector<int> infectious_Population;
        vector<int> new_dead_Population;
        for (int host = 0; host < infected_Population.size(); host++)
        {
            if (Hosts[infected_Population[host]].get_Status() == "Dead")
            {
                cout << "Node " << Hosts[infected_Population[host]].get_Name() << " is dead\n";
                dead_Count++;
                // removed_Population.push_back(infected_Population[host]);
            }
            else if (Hosts[infected_Population[host]].get_Status() == "Removed")
            {
                cout << "Node " << Hosts[infected_Population[host]].get_Name() << " removed\n";
                removed_Count++;
            }
            else if (Hosts[infected_Population[host]].terminal_status(terminal_tissues, terminal_array) == 1)
            {
                cout << "Node " << Hosts[infected_Population[host]].get_Name() << " is dead\n";
                new_dead_Population.push_back(infected_Population[host]);
                dead_Count++;
            }
            else if (Hosts[infected_Population[host]].get_Status() != "Susceptible")
            {
                temp.push_back(infected_Population[host]);
                if (Hosts[infected_Population[host]].get_Status() == "Infectious" || Hosts[infected_Population[host]].infectious_status(infectious_tissues, infectious_array) == 1)
                {
                    cout << "Node " << Hosts[infected_Population[host]].get_Name() << " is infectious\n";
                    infectious_Population.push_back(infected_Population[host]);
                }
            }
        }

        infected_Population = temp;
        temp.clear();

        /**
         * Processing the dead population data
         **/
        // // TODO: Compress all dead at once if multiread is available or else one by one
        if (new_dead_Population.size() > 0 && enable_Folder_management == "YES")
        {
            cout << "\nCompressing " << new_dead_Population.size() << " dead host(s)\n";
            //  intermediary_Sequence_location + "/" + to_string(Hosts[infected_Population[host]].get_host_Index()), enable_Folder_management, enable_Compression
            if (multi_Read == "NO")
            {
                cout << "Single thread compression\n";
                for (int host = 0; host < new_dead_Population.size(); host++)
                {
                    Hosts[new_dead_Population[host]].compress_Folder(intermediary_Sequence_location + "/" + to_string(Hosts[new_dead_Population[host]].get_host_Index()), enable_Compression);
                }
            }
            else
            {
                cout << "Multi threaded compression\n";

                // int num_per_Core = new_dead_Population.size() / CPU_cores;
                // int remainder = new_dead_Population.size() % CPU_cores;

                vector<thread> threads_vec;

                int current_thread_Usage = 0;

                for (int host = 0; host < new_dead_Population.size(); host++)
                {
                    threads_vec.push_back(thread{&simulator_Master::thread_compress_Folders, this, new_dead_Population[host], ref(Hosts)});

                    current_thread_Usage++;

                    if (current_thread_Usage >= CPU_cores)
                    {
                        for (thread &t : threads_vec)
                        {
                            if (t.joinable())
                            {
                                t.join();
                            }
                        }

                        threads_vec.clear();
                        current_thread_Usage = 0;
                    }
                }

                for (thread &t : threads_vec)
                {
                    if (t.joinable())
                    {
                        t.join();
                    }
                }

                // for (int core_ID = 0; core_ID < CPU_cores; core_ID++)
                // {
                //     int start_Node = core_ID * num_per_Core;
                //     int stop_Node = start_Node + num_per_Core;

                //     threads_vec.push_back(thread{&simulator_Master::thread_compress_Folders, this, start_Node, stop_Node, new_dead_Population, ref(Hosts)});
                // }

                // if (remainder != 0)
                // {
                //     int start_Node = new_dead_Population.size() - remainder;
                //     int stop_Node = new_dead_Population.size();

                //     threads_vec.push_back(thread{&simulator_Master::thread_compress_Folders, this, start_Node, stop_Node, new_dead_Population, ref(Hosts)});
                // }

                // for (thread &t : threads_vec)
                // {
                //     if (t.joinable())
                //     {
                //         t.join();
                //     }
                // }

                threads_vec.clear();
            }
            cout << endl;
            // exit(-1);
        }

        new_dead_Population.clear();

        // exit(-1);

        // // TODO WRITE to overall generational summary
        fstream overall_Generational_summary_File;
        overall_Generational_summary_File.open(overall_Generational_Summary, ios::app);
        if (overall_Generational_summary_File.is_open())
        {
            cout << "Writing generation " << overall_Generations << " summary: " << overall_Generational_Summary << endl;

            int year, month, day;
            functions.decimal_to_Date(decimal_Date, year, month, day);

            overall_Generational_summary_File << to_string(overall_Generations)
                                              << "\t" << to_string(decimal_Date)
                                              << "\t" << to_string(year) << "-" << to_string(month) << "-" << to_string(day)
                                              << "\t" << to_string(susceptible_Population.size() - infected_Population.size() - removed_Count - dead_Count)
                                              << "\t" << to_string(infected_Population.size())
                                              << "\t" << to_string(infectious_Population.size())
                                              << "\t" << to_string(removed_Count)
                                              << "\t" << to_string(dead_Count) << endl;
            overall_Generational_summary_File.close();
        }
        else
        {
            cout << "ERROR: UNABLE TO OPEN OVERALL GENERATIONAL SUMMARY FILE: " << overall_Generational_Summary << endl;
            exit(-1);
        }

        // TEST BLOCKS
        // if (overall_Generations == 1)
        // {
        //     cout << "\n******DONE UPTO HERE********\n";
        //     exit(-1);
        // }

        // infect
        /**
         * * INFECTING THE SUSCEPTIBLE POPULATION FROM THE INFECTIOUS POPULATION
         **/
        if (infectious_Population.size() > 0)
        {
            cout << "\nNew host to host infections\n";
            for (int host = 0; host < infectious_Population.size(); host++)
            {
                /**
                 * Get the probability of an infectious individual being able to infect a susceptible host
                 * attached to the infectious host in the contact network.
                 **/

                float infection_probability = Hosts[infectious_Population[host]].get_infection_probability();
                int go_for_Infection = 0;

                if (infection_probability == 1)
                {
                    go_for_Infection = 1;
                }
                else if (infection_probability > 0)
                {
                    /**
                     * Using the probability determine if the individual being infectious or not.
                     **/
                    bernoulli_distribution infection_draw(infection_probability);
                    bool result = infection_draw(gen);

                    if (result)
                    {
                        go_for_Infection = 1;
                    }
                }

                if (go_for_Infection == 1)
                {
                    /**
                     * If the individual is infectious we determine whom they infect using the contact network.
                     **/

                    cout << "\nNode: " << Hosts[infectious_Population[host]].get_Name() << " is infecting new nodes\n";
                    int node_Profile = Hosts[infectious_Population[host]].get_Profile();

                    /**
                     * Contact network connections for the infectious host.
                     **/
                    vector<int> host_Connections = each_Nodes_Connection_INT[infectious_Population[host]];

                    // Node_search(host_Connections);

                    // test
                    // cout << Hosts[infectious_Population[host]].get_Name() << endl
                    //      << endl;
                    // for (int test = 0; test < host_Connections.size(); test++)
                    // {
                    //     cout << host_Connections[test].first << "_" << host_Connections[test].second << endl;
                    // }

                    // cout << "\ncheck\n";

                    // for (int test = 0; test < host_Connections.size(); test++)
                    // {
                    //     cout << Hosts[search_Indexes[test]].get_Name() << endl;
                    // }

                    // exit(-1);

                    cout << "Checking index search\n";
                    vector<int> possible_Infections;

                    vector<int> new_Hosts_Indexes;

                    /**
                     * If reinfection is prevented then we make sure that even though hosts are connected to the infectious hosts
                     * if they are already infected then they will not be infected again.
                     * * ELSE they too can be infected again with newer strains by the reinfection.
                     **/
                    if (reinfection_Availability == 0)
                    {
                        cout << "Preventing reinfection of infeceted hosts\n";
                        // index and position in host_Connections

                        /**
                         * Identify connected hosts that are susceptible and not already infected
                         **/
                        for (int host = 0; host < host_Connections.size(); host++)
                        {
                            // if (search_Indexes[host] == -1)
                            // {
                            //     cout << "ERROR: ID MISSING\n";
                            //     exit(-1);
                            // }
                            // if (host_Connections[host] != -1)
                            // {
                            if (Hosts[host_Connections[host]].get_Status() == "Susceptible")
                            {
                                possible_Infections.push_back(host_Connections[host]);
                            }
                            // }
                        }
                        // search_Indexes.clear();
                        // overall_Found = 0;
                        if (possible_Infections.size() > 0)
                        {
                            new_Hosts_Indexes = get_new_Hosts_Indexes(node_Profile, gen, possible_Infections);
                        }
                        else
                        {
                            cout << "No new infection from node as all surrounding nodes have already been infected.\n";
                        }
                    }
                    else
                    {
                        cout << "Reinfection of infected hosts can occur\n";
                        // vector<pair<int, int>> host_Connections = each_Nodes_Connection[infectious_Population[host]];
                        // Node_search(host_Connections);

                        // cout << "Checking index search\n";
                        // vector<int> possible_Infections;
                        for (int host = 0; host < host_Connections.size(); host++)
                        {
                            // if (search_Indexes[host] == -1)
                            // {
                            //     cout << "ERROR: ID MISSING\n";
                            //     exit(-1);
                            // }
                            // if (host_Connections[host] != -1)
                            // {
                            /**
                             * * Only hosts that are susceptible or infected or infectious can be infected.
                             **/

                            if (Hosts[host_Connections[host]].get_Status() != "Dead" && Hosts[host_Connections[host]].get_Status() != "Removed")
                            {
                                possible_Infections.push_back(host_Connections[host]);
                            }
                            //}
                        }

                        // search_Indexes.clear();
                        // overall_Found = 0;

                        if (possible_Infections.size() > 0)
                        {
                            new_Hosts_Indexes = get_new_Hosts_Indexes(node_Profile, gen, possible_Infections);
                        }
                        else
                        {
                            cout << "No new infection from node as all surrounding nodes are dead or removed.\n";
                        }
                    }

                    // exit(-1);

                    /**
                     * Infecting target hosts
                     **/

                    if (new_Hosts_Indexes.size() > 0)
                    {
                        cout << "\nNovel host(s) recognised\n\n";

                        //// DO: 1 Make document to host to target, then sequence/ sequences of infection, TEST above code
                        int num_viruses_to_transfer = (int)transmission_parameters[1];

                        string source_Target_file_Location = intermediary_Sequence_location + "/" + to_string(Hosts[infectious_Population[host]].get_host_Index());
                        int source_Index = Hosts[infectious_Population[host]].get_host_Index();
                        int source_Generation = Hosts[infectious_Population[host]].get_Generation();
                        string source_Name = Hosts[infectious_Population[host]].get_Name();

                        //// INDEX FOLDER of Source HOST
                        vector<vector<pair<int, int>>> indexed_Source_Folders = functions.index_sequence_Folders(source_Target_file_Location, num_tissues_per_Node, source_Generation, multi_Read);
                        // for (int test = 0; test < num_tissues_per_Node; test++)
                        // {
                        //     cout << "tissue: " << test << endl;
                        //     for (size_t i = 0; i < indexed_Source_Folders[test].size(); i++)
                        //     {
                        //         cout << indexed_Source_Folders[test][i].first << "_" << indexed_Source_Folders[test][i].second << endl;
                        //     }
                        //     cout << endl;
                        // }

                        // exit(-1);

                        /**
                         * Infect the hosts
                         **/

                        for (int target_Host = 0; target_Host < new_Hosts_Indexes.size(); target_Host++)
                        {
                            if (transmission_parameters[0] == 1)
                            {
                                binomial_distribution<int> viruses_to_transfer_distribution(transmission_parameters[1], transmission_parameters[2]);
                                num_viruses_to_transfer = viruses_to_transfer_distribution(gen);
                            }

                            string status = Hosts[new_Hosts_Indexes[target_Host]].transfer_Infection(functions, intermediary_Sequence_location, source_Target_file_Location,
                                                                                                     source_Index, source_Generation, source_Name, Hosts[infectious_Population[host]].get_current_Viral_load_per_Tissue(),
                                                                                                     num_viruses_to_transfer,
                                                                                                     entry_tissues, entry_array, Hosts[infectious_Population[host]].get_Load(exit_tissues, exit_array), exit_tissues, exit_array,
                                                                                                     Hosts[infectious_Population[host]].removed_by_Transfer_Indexes,
                                                                                                     max_sequences_per_File,
                                                                                                     indexed_Source_Folders,
                                                                                                     decimal_Date,
                                                                                                     Host_source_target_network_location,
                                                                                                     output_Node_location, tissue_Names,
                                                                                                     gen);
                            // exit(-1);
                            /**
                             * UPDATE the INFECTED population
                             **/
                            if (status == "Infected")
                            {
                                infected_Population.push_back(new_Hosts_Indexes[target_Host]);
                            }
                        }
                    }
                }
            }

            // // TEST BLOCKS
            // if (overall_Generations == 1)
            // {
            //     cout << "\n******DONE UPTO HERE********\n";
            //     exit(-1);
            // }
        }

        // exit(-1);

        if (infected_Population.size() > 0)
        {
            /**
             * Each host in the infected population is processed at a time.
             **/

            cout << "\nAttempting to simulate " << infected_Population.size() << " hosts\n";
            for (int host = 0; host < infected_Population.size(); host++)
            {
                // string source_Target_file_Location = intermediary_Sequence_location + "/" + to_string(Hosts[infectious_Population[host]].get_host_Index());
                Hosts[infected_Population[host]].run_Generation(functions, this->multi_Read, this->max_Cells_at_a_time, this->gpu_Limit, CUDA_device_IDs, this->num_Cuda_devices, this->genome_Length,
                                                                CPU_cores, max_sequences_per_File,
                                                                intermediary_Sequence_location + "/" + to_string(Hosts[infected_Population[host]].get_host_Index()), output_Node_location,
                                                                tissue_Names,
                                                                num_replication_phases, tissue_replication_data, tissue_param_profile_Stride,
                                                                terminal_tissues, terminal_array,
                                                                cell_Distribution_Type, ALL_profiles_Tissue_cell_disribution[Hosts[infected_Population[host]].get_Profile()],
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
                                                                viral_Migration,
                                                                viral_Migration_Values,
                                                                migration_start_Generation,
                                                                overall_Generations,
                                                                infected_to_Recovered,
                                                                enable_Folder_management,
                                                                enable_Compression,
                                                                gen);
            }

            cout << "\nCompleted simulating hosts\n";
            decimal_Date = decimal_Date + date_Increment;

            auto stop_Time = chrono::high_resolution_clock::now();

            // Calculate the duration
            auto duration_Time = chrono::duration_cast<chrono::milliseconds>(stop_Time - start_Time);

            time_Track << to_string(overall_Generations) << "\t" << to_string(duration_Time.count()) << "\n";

            overall_Generations++;
            cout << "\nMoved generation forward\n";

            // // TEST BLOCKS
            // if (overall_Generations == 1)
            // {
            //     cout << "\n******DONE UPTO HERE********\n";
            //     exit(-1);
            // }

            if (trials_Sampling != -1)
            {
                cout << "\nSampling infected hosts\n";
                binomial_distribution<int> sequences_to_Sample(sampling_trials, sampling_probability);
                int num_Nodes_to_Sample = sequences_to_Sample(gen);

                uniform_int_distribution<> distribution_Indexes(0, infected_Population.size() - 1);
                set<int> nodes_Indexes;

                for (int host = 0; host < num_Nodes_to_Sample; host++)
                {
                    int node_Index = distribution_Indexes(gen);

                    if (resampling == -1)
                    {
                        auto it = sampled_Nodes.find(infected_Population[node_Index]);

                        if (it == sampled_Nodes.end())
                        {
                            nodes_Indexes.insert(node_Index);
                            sampled_Nodes.insert(infected_Population[node_Index]);
                        }
                    }
                    else
                    {
                        nodes_Indexes.insert(node_Index);
                    }
                }

                if (nodes_Indexes.size() > 0)
                {
                    vector<int> indexes_of_Sampling_Nodes(nodes_Indexes.begin(), nodes_Indexes.end());
                    nodes_Indexes.clear();

                    cout << "Attempting to sample " << indexes_of_Sampling_Nodes.size() << " hosts\n";

                    int num_Samples = per_Node_sampling[1];

                    // // TODO: FIX
                    int success_Sampling = 0;

                    for (int host = 0; host < indexes_of_Sampling_Nodes.size(); host++)
                    {
                        if (Hosts[infected_Population[indexes_of_Sampling_Nodes[host]]].get_Status() != "Dead" && Hosts[infected_Population[indexes_of_Sampling_Nodes[host]]].get_Status() != "Removed")
                        {
                            cout << "\nSampling host " << Hosts[infected_Population[indexes_of_Sampling_Nodes[host]]].get_Name() << endl;
                            for (int tissue = 0; tissue < sampling_tissues; tissue++)
                            {
                                cout << "Attempting to sample " << tissue_Names[sampling_array[tissue]] << " tissue\n";
                                if (per_Node_sampling[0] == 1)
                                {
                                    binomial_distribution<int> num_samples_Obtained(per_Node_sampling[1], per_Node_sampling[2]);
                                    num_Samples = num_samples_Obtained(gen);
                                }

                                if (num_Samples > 0)
                                {
                                    success_Sampling = Hosts[infected_Population[indexes_of_Sampling_Nodes[host]]].sample_Host(functions, decimal_Date,
                                                                                                                               tissue_Names,
                                                                                                                               intermediary_Sequence_location + "/" + to_string(Hosts[infected_Population[indexes_of_Sampling_Nodes[host]]].get_host_Index()), sampling_array[tissue], num_Samples,
                                                                                                                               sampled_sequences_Folder,
                                                                                                                               gen);
                                }
                            }
                            // REMOVE
                            // exit(-1);

                            if (Hosts[infected_Population[indexes_of_Sampling_Nodes[host]]].get_infection_probability() <= 0)
                            {
                                Hosts[infected_Population[indexes_of_Sampling_Nodes[host]]].set_Infection_prob_Zero(intermediary_Sequence_location + "/" + to_string(Hosts[infected_Population[indexes_of_Sampling_Nodes[host]]].get_host_Index()), enable_Folder_management, enable_Compression);
                            }
                        }
                    }

                    if (success_Sampling == 1)
                    {
                        count_Sampling_instances++;
                    }
                }
                else
                {
                    cout << "No new hosts to sample\n";
                }

                if (limit_Sampled != -1)
                {
                    if (count_Sampling_instances >= limit_Sampled)
                    {
                        stop = 3;
                    }
                }
            }

            // REMOVE
            // exit(-1);

            // TEST BLOCKS
            // if (overall_Generations == 2)
            // {
            //     cout << "\n******DONE UPTO HERE********\n";
            //     exit(-1);
            // }
        }
        else
        {
            stop = 2;
        }

        if (stop_after_generations == "YES")
        {
            if (stop_gen_Mode == 0)
            {
                if (overall_Generations >= stop_generations_Count)
                {
                    stop = 4;
                }
            }
            else
            {
                if (decimal_Date >= stop_Date)
                {
                    stop = 5;
                }
            }
        }

        // stop = 1;

        // cudaError_t err = cudaGetLastError();

        // if (err != cudaSuccess)
        // {
        //     printf("CUDA Error: %s\n", cudaGetErrorString(err));
        //     exit(-1);
        //     // Possibly: exit(-1) if program cannot continue....
        // }
        // else
        // {
        //     cout << "OK\n";
        // }

        // exit(-1);

    } while (stop == 0);

    cout << "\nSimulation has concluded: ";

    time_Track.close();

    functions.clear_Array_int_CPU(cell_Distribution_Type, number_of_node_Profiles);

    if (stop == 1)
    {
        cout << "GOD Mode STOP\n";
    }
    else if (stop == 2)
    {
        cout << "No infected population remaining\n";
    }
    else if (stop == 3)
    {
        cout << "Maximum sucessfull sequencing instances of " << limit_Sampled << " have been reached\n";
    }
    else if (stop == 4)
    {
        cout << "Maximum number of " << overall_Generations << " generations has been reached\n";
    }
    else if (stop == 5)
    {
        cout << "Maximum date of " << stop_Date << " has been reached\n";
    }

    // uniform_int_distribution<int> distribution(0, Total_number_of_Nodes - 1);
}

vector<int> simulator_Master::get_new_Hosts_Indexes(int &node_Profile, mt19937 &gen, vector<int> &possible_Infections)
{
    /**
     * * Determine who will get infected.
     **/
    vector<int> new_Hosts_Indexes;

    cout << "Novel hosts infected: ";
    int num_New_hosts = -1;
    if (infection_parameters[node_Profile][0] == 0)
    {
        num_New_hosts = infection_parameters[node_Profile][1];
    }
    else
    {
        binomial_distribution<> infections_distribution(infection_parameters[node_Profile][1], infection_parameters[node_Profile][2]);
        num_New_hosts = infections_distribution(gen);
    }

    if (possible_Infections.size() <= num_New_hosts)
    {
        new_Hosts_Indexes = possible_Infections;
        cout << possible_Infections.size() << endl;
        // cout << "x" << endl;
    }
    else
    {
        int new_Hosts = 0;
        cout << num_New_hosts << endl;
        // // CHANGE so changes in infection rate can be accounted for, by sampleing effects
        uniform_int_distribution<> distribution_Hosts(0, possible_Infections.size() - 1);
        do
        {
            int host_index = distribution_Hosts(gen);

            if (find(new_Hosts_Indexes.begin(), new_Hosts_Indexes.end(), possible_Infections[host_index]) == new_Hosts_Indexes.end())
            {
                new_Hosts_Indexes.push_back(possible_Infections[host_index]);
                new_Hosts++;
            }

        } while (new_Hosts < num_New_hosts);
    }

    return new_Hosts_Indexes;
}

// void simulator_Master::Node_search(vector<pair<int, int>> &host_Connections)
// {
//     cout << "Intiating node index search\n";
//     search_Indexes.clear();
//     overall_Found = 0;
//     for (int host = 0; host < host_Connections.size(); host++)
//     {
//         search_Indexes.push_back(-1);
//     }

//     int num_per_Core = Total_number_of_Nodes / this->CPU_cores;
//     int remainder = Total_number_of_Nodes % this->CPU_cores;

//     vector<thread> threads_vec;

//     for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
//     {
//         int start_Node = core_ID * num_per_Core;
//         int stop_Node = start_Node + num_per_Core;

//         threads_vec.push_back(thread{&simulator_Master::thread_Node_search, this, start_Node, stop_Node, host_Connections});
//     }

//     if (remainder != 0)
//     {
//         int start_Node = Total_number_of_Nodes - remainder;
//         int stop_Node = Total_number_of_Nodes;

//         threads_vec.push_back(thread{&simulator_Master::thread_Node_search, this, start_Node, stop_Node, host_Connections});
//     }

//     for (thread &t : threads_vec)
//     {
//         if (t.joinable())
//         {
//             t.join();
//         }
//     }

//     threads_vec.clear();

//     // for (int host = 0; host < host_Connections.size(); host++)
//     // {
//     //     cout << "Index: " << search_Indexes[host] << endl;
//     // }
// }

// void simulator_Master::thread_Node_search(int start_Node, int stop_Node, vector<pair<int, int>> host_Connections)
// {
//     vector<int> found_Indexes;

//     for (int host = 0; host < host_Connections.size(); host++)
//     {
//         found_Indexes.push_back(-1);
//     }

//     int found = 0;

//     for (int node = start_Node; node < stop_Node; node++)
//     {
//         if (overall_Found == host_Connections.size())
//         {
//             break;
//         }
//         for (int check = 0; check < host_Connections.size(); check++)
//         {
//             if (all_node_IDs[node].first == host_Connections[check].first && all_node_IDs[node].second == host_Connections[check].second)
//             {
//                 found_Indexes[check] = node;
//                 found++;
//                 break;
//             }
//         }
//         if (found == host_Connections.size())
//         {
//             break;
//         }
//     }

//     unique_lock<shared_mutex> ul(g_mutex);
//     for (int host = 0; host < host_Connections.size(); host++)
//     {
//         if (found_Indexes[host] != -1)
//         {
//             search_Indexes[host] = found_Indexes[host];
//             overall_Found++;
//         }
//     }
// }

int simulator_Master::get_first_Infected(vector<int> &susceptible_Population,
                                         vector<int> &infected_Population, functions_library &functions)
{
    int node_infected = -1;

    random_device rd;
    mt19937 gen(rd());

    uniform_int_distribution<int> distribution(0, susceptible_Population.size() - 1);

    int susceptible_Index = distribution(gen);
    node_infected = susceptible_Population[susceptible_Index];

    cout << "Node " << all_node_IDs[node_infected].first << "_" << all_node_IDs[node_infected].second << " infected first\n\n";

    infected_Population.push_back(node_infected);
    sort(infected_Population.begin(), infected_Population.end());

    // CHECK change
    // vector<int> susceptible_Population_Temp;

    // for (int fill_Temp = 0; fill_Temp < susceptible_Population.size(); fill_Temp++)
    // {
    //     if (fill_Temp != susceptible_Index)
    //     {
    //         susceptible_Population_Temp.push_back(susceptible_Population[fill_Temp]);
    //     }
    // }

    // susceptible_Population = susceptible_Population_Temp;

    cout << "Infecting node with reference genomes\n";

    intermediary_Sequence_location = intermediate_Folder_location + "/sequence_Data";
    functions.config_Folder(intermediary_Sequence_location, "Intermediary sequence data");

    intermediary_Index_location = intermediate_Folder_location + "/index_Data";
    functions.config_Folder(intermediary_Index_location, "Node index data");

    string reference_Sequences = intermediary_Sequence_location + "/reference_Sequences";
    functions.config_Folder(reference_Sequences, "Converted reference sequence(s)");

    // functions.process_Reference_Sequences(read_Reference_Sequences(node_infected),genome_Length);
    vector<int> tissue_Sequence_Count;
    vector<string> collect_Sequences = read_Reference_Sequences(tissue_Sequence_Count);
    // exit(-1);

    int total_Sequences = collect_Sequences.size();

    cout << "\nProcessing " << total_Sequences << " collected sequence(s)\n";

    // exit(-1);

    if (first_Infection == "RANDOM")
    {

        int full_Rounds = total_Sequences / this->gpu_Limit;
        int partial_Rounds = total_Sequences % this->gpu_Limit;

        vector<pair<int, int>> start_stops;

        for (int full = 0; full < full_Rounds; full++)
        {
            int start = full * this->gpu_Limit;
            int stop = start + this->gpu_Limit;
            start_stops.push_back(make_pair(start, stop));
        }

        if (partial_Rounds != 0)
        {
            int start = total_Sequences - partial_Rounds;
            start_stops.push_back(make_pair(start, total_Sequences));
        }

        vector<string> sequence_Write_Store_All;
        int last_seq_Num = 0;

        vector<char> seq_Status;

        for (int round = 0; round < start_stops.size(); round++)
        {
            cout << "\nExecuting " << round + 1 << " of " << start_stops.size() << " rounds\n";

            int num_of_Sequences_current = start_stops[round].second - start_stops[round].first;
            // vector<string> collect_Sequences, int &genome_Length, int &round, vector<pair<int, int>> &start_stops, int num_of_Sequences_current

            int **sequences = functions.process_Reference_Sequences(collect_Sequences, genome_Length, num_of_Sequences_current);

            vector<string> sequence_Write_Store = functions.convert_Sequences_Master(sequences, genome_Length, num_of_Sequences_current);

            functions.clear_Array_int_CPU(sequences, num_of_Sequences_current);

            functions.sequence_Write_Configurator(sequence_Write_Store_All, sequence_Write_Store,
                                                  max_sequences_per_File, reference_Sequences, last_seq_Num, seq_Status);

            // for (int row = 0; row < num_of_Sequences_current; row++)
            // {
            //     for (size_t c = 0; c < genome_Length; c++)
            //     {
            //         cout << sequence[row][c];
            //     }
            //     cout << "\n\n";
            // }

            // functions.clear_Array_int_CPU(sequences, num_of_Sequences_current);
        }

        functions.partial_Write_Check(sequence_Write_Store_All, reference_Sequences, last_seq_Num, seq_Status);
        collect_Sequences.clear();
    }
    else
    {
        vector<vector<string>> collect_Sequences_Tissue;
        for (int tissue = 0; tissue < tissue_Names.size(); tissue++)
        {
            vector<string> tissue_Sequences;
            for (int start = tissue_Sequence_Count[tissue]; start < tissue_Sequence_Count[tissue + 1]; start++)
            {
                tissue_Sequences.push_back(collect_Sequences[start]);
            }
            collect_Sequences_Tissue.push_back(tissue_Sequences);
        }
        collect_Sequences.clear();
        tissue_Sequence_Count.clear();

        for (int tissue = 0; tissue < tissue_Names.size(); tissue++)
        {
            cout << "\nConverting sequences of tissue: " << tissue_Names[tissue] << endl;

            string reference_Sequences_Tissue = reference_Sequences + "/" + tissue_Names[tissue];
            functions.config_Folder(reference_Sequences_Tissue, tissue_Names[tissue] + " converted reference sequence(s)");

            int full_Rounds = collect_Sequences_Tissue[tissue].size() / this->gpu_Limit;
            int partial_Rounds = collect_Sequences_Tissue[tissue].size() % this->gpu_Limit;

            vector<pair<int, int>> start_stops;

            for (int full = 0; full < full_Rounds; full++)
            {
                int start = full * this->gpu_Limit;
                int stop = start + this->gpu_Limit;
                start_stops.push_back(make_pair(start, stop));
            }

            if (partial_Rounds != 0)
            {
                int start = collect_Sequences_Tissue[tissue].size() - partial_Rounds;
                start_stops.push_back(make_pair(start, collect_Sequences_Tissue[tissue].size()));
            }

            vector<string> sequence_Write_Store_All;
            int last_seq_Num = 0;

            vector<char> seq_Status;

            for (int round = 0; round < start_stops.size(); round++)
            {
                cout << "\nExecuting " << round + 1 << " of " << start_stops.size() << " rounds\n";

                int num_of_Sequences_current = start_stops[round].second - start_stops[round].first;
                // vector<string> collect_Sequences, int &genome_Length, int &round, vector<pair<int, int>> &start_stops, int num_of_Sequences_current

                int **sequences = functions.process_Reference_Sequences(collect_Sequences_Tissue[tissue], genome_Length, num_of_Sequences_current);

                vector<string> sequence_Write_Store = functions.convert_Sequences_Master(sequences, genome_Length, num_of_Sequences_current);

                functions.clear_Array_int_CPU(sequences, num_of_Sequences_current);

                functions.sequence_Write_Configurator(sequence_Write_Store_All, sequence_Write_Store,
                                                      max_sequences_per_File, reference_Sequences_Tissue, last_seq_Num, seq_Status);
            }
            functions.partial_Write_Check(sequence_Write_Store_All, reference_Sequences_Tissue, last_seq_Num, seq_Status);
        }
    }

    cout << endl;
    // exit(-1);

    return node_infected;
}

vector<string> simulator_Master::read_Reference_Sequences(vector<int> &tissue_Sequence_Count)
{
    // parent_Sequence_Folder
    cout << "Reading parent sequence folder: " << parent_Sequence_Folder << endl;

    if (filesystem::exists(parent_Sequence_Folder) && filesystem::is_directory(parent_Sequence_Folder))
    {
        if (first_Infection == "RANDOM")
        {
            vector<string> sequences_Paths;
            for (const auto &entry : filesystem::directory_iterator(parent_Sequence_Folder))
            {
                if (filesystem::is_regular_file(entry))
                {
                    string check_Extenstion = entry.path().extension();
                    if (check_Extenstion == ".fasta" || check_Extenstion == ".fa" || check_Extenstion == ".nfa" || check_Extenstion == ".nfasta")
                    {
                        // cout << "Found sequence: " << entry.path() << endl;
                        sequences_Paths.push_back(entry.path().string());
                        // cout << "Found sequence: " << sequences_Paths[sequences_Paths.size() - 1] << endl;
                    }
                }
            }

            if (sequences_Paths.size() > 0)
            {

                cout << "Identified " << sequences_Paths.size() << " parent sequence file(s)\n\n";

                vector<string> collect_Sequences = read_Reference_Sequence_Files(sequences_Paths);
                if (collect_Sequences.size() > 0)
                {
                    cout << "\nIdentified " << collect_Sequences.size() << " parent sequences\n";

                    cout << "Validating collected parent sequences\n";
                    genome_Length = collect_Sequences[0].size();

                    for (int genome = 1; genome < collect_Sequences.size(); genome++)
                    {
                        if (collect_Sequences[genome].size() != genome_Length)
                        {
                            cout << "ERROR ALL GENOMES MUST BE OF EQUAL LENGTH\n";
                            exit(-1);
                        }
                    }
                    cout << "All sequences are of valid lenth: " << genome_Length << endl;
                    return collect_Sequences;
                }
                else
                {
                    cout << "ERROR: THERE SHOULD BE MORE THAN ZERO PARENT SEQUENCES.\n";
                    exit(-1);
                }
            }
            else
            {
                cout << "ERROR: THERE SHOULD BE MORE THAN ZERO PARENT SEQUENCE FILES.\n";
                exit(-1);
            }
        }
        else
        {
            tissue_Sequence_Count.push_back(0);
            vector<string> collect_Sequences_Full;
            for (int tissue = 0; tissue < tissue_Names.size(); tissue++)
            {
                vector<string> sequences_Paths;
                vector<string> collect_Sequences;

                if (filesystem::exists(parent_Sequence_Folder + "/" + tissue_Names[tissue]) && filesystem::is_directory(parent_Sequence_Folder + "/" + tissue_Names[tissue]))
                {
                    cout << "\nReading reference tissue: " << tissue_Names[tissue] << endl;
                    for (const auto &entry : filesystem::directory_iterator(parent_Sequence_Folder + "/" + tissue_Names[tissue]))
                    {
                        if (filesystem::is_regular_file(entry))
                        {
                            string check_Extenstion = entry.path().extension();
                            if (check_Extenstion == ".fasta" || check_Extenstion == ".fa" || check_Extenstion == ".nfa" || check_Extenstion == ".nfasta")
                            {
                                // cout << "Found sequence: " << entry.path() << endl;
                                sequences_Paths.push_back(entry.path().string());
                                // cout << "Found sequence: " << sequences_Paths[sequences_Paths.size() - 1] << endl;
                            }
                        }
                    }
                    if (sequences_Paths.size() > 0)
                    {
                        cout << "Identified " << sequences_Paths.size() << " parent sequence file(s)\n\n";
                        collect_Sequences = read_Reference_Sequence_Files(sequences_Paths);

                        if (collect_Sequences.size() > 0)
                        {
                            cout << "\nIdentified " << collect_Sequences.size() << " parent sequences\n";
                            cout << "Validating collected parent sequences\n";
                            if (genome_Length == 0)
                            {
                                genome_Length = collect_Sequences[0].size();
                            }
                            for (int genome = 0; genome < collect_Sequences.size(); genome++)
                            {
                                if (collect_Sequences[genome].size() != genome_Length)
                                {
                                    cout << "ERROR ALL GENOMES MUST BE OF EQUAL LENGTH\n";
                                    exit(-1);
                                }
                                else
                                {
                                    collect_Sequences_Full.push_back(collect_Sequences[genome]);
                                }
                            }
                        }
                        else
                        {
                            cout << "ERROR: THERE SHOULD BE MORE THAN ZERO PARENT SEQUENCES.\n";
                            exit(-1);
                        }
                    }
                    else
                    {
                        cout << "ERROR: THERE SHOULD BE MORE THAN ZERO PARENT SEQUENCE FILES.\n";
                        exit(-1);
                    }
                }
                tissue_Sequence_Count.push_back(collect_Sequences_Full.size());
            }
            cout << "\nAll sequences are of valid lenth: " << genome_Length << endl;
            return collect_Sequences_Full;
        }
    }
    else
    {
        cout << "ERROR: PARENT SEQUENCE FOLDER DOES NOT EXIST AT THE GIVEN PATH: " << parent_Sequence_Folder << endl;
        exit(-1);
    }
}

vector<string> simulator_Master::read_Reference_Sequence_Files(vector<string> &reference_Files)
{
    vector<string> collect_Sequences;
    // cout << "Reading FASTA file: \n";

    for (int file_Index = 0; file_Index < reference_Files.size(); file_Index++)
    {
        // cout << "Loop\n";
        // cout << "Reading FASTA file: " << reference_Files[file_Index];
        fstream fasta_Read;
        fasta_Read.open(reference_Files[file_Index], ios::in);
        if (fasta_Read.is_open())
        {
            cout << "Reading FASTA file: " << reference_Files[file_Index] << "\n";

            string line;
            string sequence = "";

            while (getline(fasta_Read, line))
            {
                if (line.at(0) == '>')
                {
                    if (sequence != "")
                    {
                        collect_Sequences.push_back(sequence);
                        sequence = "";
                    }
                }
                else
                {
                    sequence.append(line);
                }
            }

            if (sequence != "")
            {
                collect_Sequences.push_back(sequence);
            }

            fasta_Read.close();
        }
        else
        {
            cout << "ERROR: UNABLE TO OPEN FASTA FILE: " << reference_Files[file_Index] << endl;
            exit(-1);
        }
    }

    return collect_Sequences;
}

vector<node_within_host> simulator_Master::node_Profile_assignment_Manager(functions_library &functions)
{
    vector<node_within_host> Hosts;
    cout << "Configuring Node Profile arrays\n";

    // each_Node_Profile = (int *)malloc(sizeof(int) * Total_number_of_Nodes);
    each_Node_Profile_Configuration = functions.create_FLOAT_2D_arrays(Total_number_of_Nodes, 5 + num_tissues_per_Node);

    // cout << "Configuring Node Profile arrays\n";

    int num_per_Core = Total_number_of_Nodes / this->CPU_cores;
    int remainder = Total_number_of_Nodes % this->CPU_cores;

    vector<thread> threads_vec;

    for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
    {
        int start_Node = core_ID * num_per_Core;
        int stop_Node = start_Node + num_per_Core;

        threads_vec.push_back(thread{&simulator_Master::node_Profile_assignment_thread, this, start_Node, stop_Node});
    }

    if (remainder != 0)
    {
        int start_Node = Total_number_of_Nodes - remainder;
        int stop_Node = Total_number_of_Nodes;

        threads_vec.push_back(thread{&simulator_Master::node_Profile_assignment_thread, this, start_Node, stop_Node});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    cout << "Node Profiles arrays configured: " << Total_number_of_Nodes << "\n\n";

    string standard_Node_File_location = output_Network_location + "/nodes_Configuration.csv";
    string header = "ID\tCave_Index\tProfile_name\tGenerations_projected\tInfectious_load\tTerminal_load\tSampling_effect";

    for (int tissue = 0; tissue < num_tissues_per_Node; tissue++)
    {
        header = header + "\tTissue_" + to_string(tissue + 1) + "_limit";
    }

    functions.create_File(standard_Node_File_location, header);

    fstream node_Configuration;
    node_Configuration.open(standard_Node_File_location, ios::app);

    if (node_Configuration.is_open())
    {
        cout << "Writing node(s) configurations: " << standard_Node_File_location << endl;
        for (int node = 0; node < Total_number_of_Nodes; node++)
        {
            Hosts.push_back(node_within_host());

            node_Configuration << all_node_IDs[node].first << "_" << all_node_IDs[node].second
                               << "\t" << all_node_IDs[node].first
                               << "\t" << profile_names[each_Node_Profile_Configuration[node][0]];

            Hosts[Hosts.size() - 1].setHost(node, all_node_IDs[node].first, all_node_IDs[node].second, (int)each_Node_Profile_Configuration[node][0], num_tissues_per_Node);

            for (int col = 1; col < 4; col++)
            {
                node_Configuration << "\t" << each_Node_Profile_Configuration[node][col];
            }

            Hosts[Hosts.size() - 1].setNum_Generation((int)each_Node_Profile_Configuration[node][1]);
            Hosts[Hosts.size() - 1].setInfectious_Load((int)each_Node_Profile_Configuration[node][2]);
            Hosts[Hosts.size() - 1].setTerminal_Load((int)each_Node_Profile_Configuration[node][3]);

            if (each_Node_Profile_Configuration[node][4] == 1)
            {
                node_Configuration << "\tNo Change";
            }
            else if (each_Node_Profile_Configuration[node][4] == 0)
            {
                node_Configuration << "\tRemoved";
            }
            else
            {
                node_Configuration << "\t" << each_Node_Profile_Configuration[node][4];
            }

            Hosts[Hosts.size() - 1].setSampling_Effect(each_Node_Profile_Configuration[node][4]);

            vector<int> tissue_Limits;
            for (int col = 5; col < (5 + num_tissues_per_Node); col++)
            {
                if (each_Node_Profile_Configuration[node][col] == -1)
                {
                    node_Configuration << "\tUnlimited";
                }
                else
                {
                    node_Configuration << "\t" << each_Node_Profile_Configuration[node][col];
                }
                tissue_Limits.push_back((int)each_Node_Profile_Configuration[node][col]);
            }
            Hosts[Hosts.size() - 1].setCell_Limit(tissue_Limits);
            node_Configuration << "\n";
        }
        node_Configuration.close();
    }
    else
    {
        cout << "ERROR CANNOT FIND NODE CONFIGURATION FILE LOCATION: " << standard_Node_File_location << endl;
        exit(-1);
    }

    cout << "\nPurging network memory\n";

    if (network_Model == "DCM")
    {
        free(per_cave_Stride);
    }
    free(node_profile_Distributions);

    functions.clear_Array_float_CPU(node_sampling_effect, number_of_node_Profiles);
    functions.clear_Array_float_CPU(infectious_load_Profiles_param, number_of_node_Profiles);
    functions.clear_Array_float_CPU(terminal_load_Profiles_param, number_of_node_Profiles);
    functions.clear_Array_float_CPU(profile_tissue_Limits, num_tissues_per_Node * number_of_node_Profiles);
    functions.clear_Array_float_CPU(each_Node_Profile_Configuration, Total_number_of_Nodes);

    cout << endl;

    // for (int i = 0; i < number_of_node_Profiles; i++)
    // {
    //     for (size_t t = 0; t < num_tissues_per_Node; t++)
    //     {
    //         for (size_t c = 0; c < 3; c++)
    //         {
    //             cout << profile_tissue_Limits[(i * num_tissues_per_Node) + t][c] << "\t";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }

    // // profile_tissue_Limits[(profile_Assignments[profile_Assignments.size() - 1] * num_tissues_per_Node) + tissue][0]

    // cout << "\n*****************************\n";
    // for (size_t i = 0; i < Total_number_of_Nodes; i++)
    // {
    //     for (size_t c = 0; c < (5 + num_tissues_per_Node); c++)
    //     {
    //         cout << each_Node_Profile_Configuration[i][c] << "\t";
    //     }
    //     cout << endl;
    // }
    return Hosts;
}

void simulator_Master::node_Profile_assignment_thread(int start_Node, int stop_Node)
{
    vector<int> profile_Assignments;
    vector<int> number_of_generations_per_node;
    vector<pair<int, int>> infectious_terminal_loads;
    vector<float> sampling_effect;
    vector<vector<int>> tissue_cell_Limits_per_Node;

    random_device rd;
    mt19937 gen(rd());

    gamma_distribution<float> days_in_Host(shape_days_in_Host, scale_days_in_Host);

    for (int node = start_Node; node < stop_Node; node++)
    {
        number_of_generations_per_node.push_back((int)(days_in_Host(gen) / generation_Time));

        int infectious_Load = 0;
        int terminal_Load = 0;

        // cout << "Configuring node: " << node + 1 << endl;
        float randomValue = (float)rand() / RAND_MAX;
        float cum_Prob = 0;

        for (int profile = 0; profile < number_of_node_Profiles; profile++)
        {
            cum_Prob += node_profile_Distributions[profile];
            if (randomValue <= cum_Prob)
            {
                profile_Assignments.push_back(profile);
                break;
            }
        }

        if (infectious_load_Profiles_param[profile_Assignments[profile_Assignments.size() - 1]][0] == -1)
        {
            infectious_Load = infectious_load_Profiles_param[profile_Assignments[profile_Assignments.size() - 1]][1];
        }
        else
        {
            binomial_distribution<int> infectious_Load_distribution(infectious_load_Profiles_param[profile_Assignments[profile_Assignments.size() - 1]][1], infectious_load_Profiles_param[profile_Assignments[profile_Assignments.size() - 1]][2]);
            infectious_Load = infectious_Load_distribution(gen);
        }

        if (terminal_load_Profiles_param[profile_Assignments[profile_Assignments.size() - 1]][0] == -1)
        {
            terminal_Load = terminal_load_Profiles_param[profile_Assignments[profile_Assignments.size() - 1]][1];
        }
        else
        {
            binomial_distribution<int> terminal_Load_distribution(terminal_load_Profiles_param[profile_Assignments[profile_Assignments.size() - 1]][1], terminal_load_Profiles_param[profile_Assignments[profile_Assignments.size() - 1]][2]);
            terminal_Load = terminal_Load_distribution(gen);
        }

        infectious_terminal_loads.push_back(make_pair(infectious_Load, terminal_Load));

        if (trials_Sampling == 1)
        {
            if (node_sampling_effect[profile_Assignments[profile_Assignments.size() - 1]][0] == 0)
            {
                sampling_effect.push_back(1);
            }
            else if (node_sampling_effect[profile_Assignments[profile_Assignments.size() - 1]][0] == 1)
            {
                sampling_effect.push_back(0);
            }
            else
            {
                functions_library functions = functions_library();
                sampling_effect.push_back(functions.beta_Distribution(node_sampling_effect[profile_Assignments[profile_Assignments.size() - 1]][1], node_sampling_effect[profile_Assignments[profile_Assignments.size() - 1]][2], gen));
            }
        }
        else
        {
            sampling_effect.push_back(-1);
        }

        vector<int> tissue_Limits;

        // profile * num_tissues_per_Node

        for (int tissue = 0; tissue < num_tissues_per_Node; tissue++)
        {
            if (profile_tissue_Limits[(profile_Assignments[profile_Assignments.size() - 1] * num_tissues_per_Node) + tissue][0] == 0)
            {
                tissue_Limits.push_back(-1);
            }
            else
            {
                binomial_distribution<int> tissue_limit_distribution(profile_tissue_Limits[(profile_Assignments[profile_Assignments.size() - 1] * num_tissues_per_Node) + tissue][1], profile_tissue_Limits[(profile_Assignments[profile_Assignments.size() - 1] * num_tissues_per_Node) + tissue][2]);
                tissue_Limits.push_back(tissue_limit_distribution(gen));
            }
        }

        tissue_cell_Limits_per_Node.push_back(tissue_Limits);
    }

    int index = 0;
    unique_lock<shared_mutex> ul(g_mutex);
    for (int node = start_Node; node < stop_Node; node++)
    {
        each_Node_Profile_Configuration[node][0] = profile_Assignments[index];
        each_Node_Profile_Configuration[node][1] = number_of_generations_per_node[index];
        each_Node_Profile_Configuration[node][2] = infectious_terminal_loads[index].first;
        each_Node_Profile_Configuration[node][3] = infectious_terminal_loads[index].second;
        each_Node_Profile_Configuration[node][4] = sampling_effect[index];

        for (int tissue = 0; tissue < num_tissues_per_Node; tissue++)
        {
            each_Node_Profile_Configuration[node][5 + tissue] = tissue_cell_Limits_per_Node[index][tissue];
        }

        index++;
    }
}

void simulator_Master::sequence_Master_Manager(functions_library &functions)
{
    parameter_load Parameters = parameter_load();
    cout << "Loading sequence master profile: " << this->sequence_Master_location << endl;

    vector<string> parameters_List = {
        "\"Parent sequences folder\"",
        "\"Mutation availability\"",
        "\"Recombination availability\"",
        "\"Reference Fitness\"",
        "\"Reference Survivability\"",
        "\"Transmitted sequence distribution\"",
        "\"Reinfection availability\""};

    vector<string> found_Parameters = Parameters.get_parameters(sequence_Master_location, parameters_List);

    parent_Sequence_Folder = Parameters.get_STRING(found_Parameters[0]);
    cout << "\nParent sequences folder: " << parent_Sequence_Folder << endl;

    // cout << "\nConfiguration of host to host infection\n";
    // infection_parameters = (float *)malloc(sizeof(float) * 3);
    // vector<string> infected_Parameters;
    // cout << "Infection distribution: ";
    // if (functions.to_Upper_Case(Parameters.get_STRING(found_Parameters[7])) == "FIXED")
    // {
    //     infection_parameters[0] = 0;
    //     cout << "FIXED\n";
    //     parameters_List = {"\"Infection sequence Fixed\""};
    //     infected_Parameters = Parameters.get_parameters(sequence_Master_location, parameters_List);
    //     infection_parameters[1] = Parameters.get_INT(infected_Parameters[0]);
    //     cout << "Infection sequence Fixed: " << infection_parameters[1] << endl;
    // }
    // else if (functions.to_Upper_Case(Parameters.get_STRING(found_Parameters[7])) == "BINOMIAL")
    // {
    //     infection_parameters[0] = 1;
    //     cout << "BINOMIAL\n";
    //     parameters_List = {"\"Infection sequence Binomial trials\"",
    //                        "\"Infection sequence Binomial probability\""};
    //     infected_Parameters = Parameters.get_parameters(sequence_Master_location, parameters_List);
    //     infection_parameters[1] = Parameters.get_INT(infected_Parameters[0]);
    //     infection_parameters[2] = Parameters.get_FLOAT(infected_Parameters[1]);
    //     cout << "Infection sequence Binomial trials: " << infection_parameters[1] << endl;
    //     cout << "Infection sequence Binomial probability: " << infection_parameters[2] << endl;
    // }
    // else
    // {
    //     cout << "ERROR Infection sequence distribution HAS TO BE EITHER \"FIXED\" OR \"NEGATIVE BINOMIAL\"";
    //     exit(1);
    // }

    // exit(-1);

    cout << "\nReinfection availability: ";
    if (functions.to_Upper_Case(Parameters.get_STRING(found_Parameters[6])) == "YES")
    {
        reinfection_Availability = 1;
        cout << "YES\n";
    }
    else
    {
        cout << "NO\n";
    }

    // exit(-1);

    cout << "\nConfiguration of Viral Tranmssion bottleneck\n";
    transmission_parameters = (float *)malloc(sizeof(float) * 3);
    vector<string> transmitted_Parameters;
    cout << "Transmission distribution: ";
    if (functions.to_Upper_Case(Parameters.get_STRING(found_Parameters[5])) == "FIXED")
    {
        transmission_parameters[0] = 0;
        cout << "FIXED\n";
        parameters_List = {"\"Transmitted sequence Fixed\""};
        transmitted_Parameters = Parameters.get_parameters(sequence_Master_location, parameters_List);
        transmission_parameters[1] = Parameters.get_INT(transmitted_Parameters[0]);
        cout << "Transmitted sequence Fixed: " << transmission_parameters[1] << endl;
    }
    else if (functions.to_Upper_Case(Parameters.get_STRING(found_Parameters[5])) == "BINOMIAL")
    {
        transmission_parameters[0] = 1;
        cout << "BINOMIAL\n";
        parameters_List = {"\"Transmitted sequence Binomial trials\"",
                           "\"Transmitted sequence Binomial probability\""};
        transmitted_Parameters = Parameters.get_parameters(sequence_Master_location, parameters_List);
        transmission_parameters[1] = Parameters.get_INT(transmitted_Parameters[0]);
        transmission_parameters[2] = Parameters.get_FLOAT(transmitted_Parameters[1]);
        cout << "Transmitted sequence Binomial trials: " << transmission_parameters[1] << endl;
        cout << "Transmitted sequence Binomial probability: " << transmission_parameters[2] << endl;
    }
    else
    {
        cout << "ERROR Transmitted sequence distribution HAS TO BE EITHER \"FIXED\" OR \"NEGATIVE BINOMIAL\"";
        exit(1);
    }

    cout << "\nConfiguring reference genome parameters:\n";
    Reference_fitness_survivability_proof_reading = (float *)malloc(sizeof(float) * 3);

    cout << "Reference Fitness: ";
    Reference_fitness_survivability_proof_reading[0] = Parameters.get_FLOAT(found_Parameters[3]);
    cout << Reference_fitness_survivability_proof_reading[0] << endl;

    cout << "Reference Survivability: ";
    Reference_fitness_survivability_proof_reading[1] = Parameters.get_FLOAT(found_Parameters[4]);
    cout << Reference_fitness_survivability_proof_reading[1] << endl;

    mutation_recombination_proof_Reading_availability = (int *)malloc(sizeof(int) * 3);
    cout << "\nSequence mechanisms: \n";

    string status;

    for (int i = 1; i <= 2; i++)
    {
        status = Parameters.get_STRING(found_Parameters[i]);
        transform(status.begin(), status.end(), status.begin(), ::toupper);

        if (i == 1)
        {
            cout << "Mutations: ";
        }
        else
        {
            cout << "Recombinations: ";
        }
        if (status == "YES")
        {
            mutation_recombination_proof_Reading_availability[i - 1] = 1;
            cout << "Active\n";

            if (i == 1)
            {
                vector<string> parameter_Proof_Reading = {"\"Proof reading availability\""};
                vector<string> found_Proof_Reading = Parameters.get_parameters(sequence_Master_location, parameter_Proof_Reading);

                transform(found_Proof_Reading[0].begin(), found_Proof_Reading[0].end(), found_Proof_Reading[0].begin(), ::toupper);

                cout << "Proof reading: ";
                if (Parameters.get_STRING(found_Proof_Reading[0]) == "YES")
                {
                    mutation_recombination_proof_Reading_availability[2] = 1;
                    cout << "Active\n";

                    parameter_Proof_Reading = {"\"Reference Proof Reading\""};
                    found_Proof_Reading = Parameters.get_parameters(sequence_Master_location, parameter_Proof_Reading);
                    cout << "Reference Proof reading: ";
                    Reference_fitness_survivability_proof_reading[2] = Parameters.get_FLOAT(found_Proof_Reading[0]);
                    cout << Reference_fitness_survivability_proof_reading[2] << endl;
                }
                else
                {
                    mutation_recombination_proof_Reading_availability[2] = 0;
                    Reference_fitness_survivability_proof_reading[2] = -1;
                    cout << "Not active\n";
                }
            }
        }
        else
        {
            mutation_recombination_proof_Reading_availability[i - 1] = 0;
            if (i == 1)
            {
                mutation_recombination_proof_Reading_availability[2] = 0;
                Reference_fitness_survivability_proof_reading[2] = -1;
            }
            cout << "Not active\n";
        }
    }

    cout << "\nConfiguring profiles:\n\n";
    tot_prob_selectivity = (int *)malloc(sizeof(int) * 2);
    for (int fill = 0; fill < 2; fill++)
    {
        tot_prob_selectivity[fill] = 0;
    }

    if (mutation_recombination_proof_Reading_availability[0] == 1 || mutation_recombination_proof_Reading_availability[1] == 1)
    {
        parameters_List = {
            "\"Fitness profile file\"",
            "\"Survivability profile file\""};

        found_Parameters = Parameters.get_parameters(sequence_Master_location, parameters_List);

        string fitness_Profile_Location = Parameters.get_STRING(found_Parameters[0]);
        string survivability_Profile_Location = Parameters.get_STRING(found_Parameters[1]);

        num_effect_Segregating_sites = (int *)malloc(sizeof(int) * 3);

        if (fitness_Profile_Location != "NA")
        {
            sequence_Fitness_changes = Parameters.get_Profile_Array(fitness_Profile_Location, num_effect_Segregating_sites[0], functions);

            cout << "Printing Fitness matrix:\n";

            for (int row = 0; row < num_effect_Segregating_sites[0]; row++)
            {
                for (int col = 0; col < 5; col++)
                {
                    cout << sequence_Fitness_changes[row][col] << "\t";
                }
                cout << endl;
            }
            cout << endl;
        }
        else
        {
            cout << "No fitness profile\n\n";
            num_effect_Segregating_sites[0] = 0;
        }

        if (survivability_Profile_Location != "NA")
        {
            sequence_Survivability_changes = Parameters.get_Profile_Array(survivability_Profile_Location, num_effect_Segregating_sites[1], functions);

            cout << "Printing Survivability matrix:\n";

            for (int row = 0; row < num_effect_Segregating_sites[1]; row++)
            {
                for (int col = 0; col < 5; col++)
                {
                    cout << sequence_Survivability_changes[row][col] << "\t";
                }
                cout << endl;
            }
            cout << endl;
        }
        else
        {
            cout << "No survivability profile\n\n";
            num_effect_Segregating_sites[1] = 0;
        }

        if (mutation_recombination_proof_Reading_availability[0] == 1)
        {
            if (mutation_recombination_proof_Reading_availability[2] == 1)
            {
                parameters_List = {
                    "\"Proof reading profile file\""};

                found_Parameters = Parameters.get_parameters(sequence_Master_location, parameters_List);
                string proof_Reading_Profile_Location = Parameters.get_STRING(found_Parameters[0]);

                if (proof_Reading_Profile_Location != "NA")
                {
                    sequence_Proof_reading_changes = Parameters.get_Profile_Array(proof_Reading_Profile_Location, num_effect_Segregating_sites[2], functions);

                    cout << "Printing Proof reading matrix:\n";

                    for (int row = 0; row < num_effect_Segregating_sites[2]; row++)
                    {
                        for (int col = 0; col < 5; col++)
                        {
                            cout << sequence_Proof_reading_changes[row][col] << "\t";
                        }
                        cout << endl;
                    }
                }
                else
                {
                    cout << "No proof reading profile\n\n";
                    num_effect_Segregating_sites[2] = 0;
                }
            }
            else
            {
                num_effect_Segregating_sites[2] = 0;
            }

            // CREATE MUTATION ARRAYs
            vector<pair<string, string>> mutations_Block = Parameters.get_block_from_File(sequence_Master_location, "Mutations");

            mutation_Hotspots = Parameters.get_INT(mutations_Block, "Number of hotspots");

            if (mutation_Hotspots > 0)
            {
                cout << "\nProcessing " << this->mutation_Hotspots << " mutation hotspots: \n";
                // exit(-1);
                A_0_mutation = functions.create_Fill_2D_array_FLOAT(mutation_Hotspots, 4, -1);
                T_1_mutation = functions.create_Fill_2D_array_FLOAT(mutation_Hotspots, 4, -1);
                G_2_mutation = functions.create_Fill_2D_array_FLOAT(mutation_Hotspots, 4, -1);
                C_3_mutation = functions.create_Fill_2D_array_FLOAT(mutation_Hotspots, 4, -1);

                mutation_hotspot_parameters = functions.create_Fill_2D_array_FLOAT(mutation_Hotspots, 5, -1);

                for (int mutation_hotspot = 0; mutation_hotspot < mutation_Hotspots; mutation_hotspot++)
                {
                    string hotspot_ID = "Hotspot " + to_string(mutation_hotspot + 1);
                    cout << "\nProcessing: " << hotspot_ID << endl;
                    vector<pair<string, string>> mutations_hotspot_Block = Parameters.get_block_from_block(mutations_Block, hotspot_ID);

                    string region = Parameters.get_STRING(mutations_hotspot_Block, "Region");
                    string clock_Model = Parameters.get_STRING(mutations_hotspot_Block, "Clock model");
                    transform(clock_Model.begin(), clock_Model.end(), clock_Model.begin(), ::toupper);

                    vector<string> split_Region;
                    functions.split(split_Region, region, '_');

                    mutation_hotspot_parameters[mutation_hotspot][0] = stof(split_Region[0]);
                    mutation_hotspot_parameters[mutation_hotspot][1] = stof(split_Region[1]);
                    cout << "Region: From " << mutation_hotspot_parameters[mutation_hotspot][0] << " to " << mutation_hotspot_parameters[mutation_hotspot][1] << endl;

                    cout << "Clock model: " << clock_Model << endl;

                    if (clock_Model == "POISSON")
                    {
                        mutation_hotspot_parameters[mutation_hotspot][2] = 0;
                        mutation_hotspot_parameters[mutation_hotspot][3] = Parameters.get_FLOAT(mutations_hotspot_Block, "Clock model Poisson mean");

                        cout << "Poisson mean: " << mutation_hotspot_parameters[mutation_hotspot][3] << endl;
                    }
                    else if (clock_Model == "NEGATIVE BINOMIAL")
                    {
                        mutation_hotspot_parameters[mutation_hotspot][2] = 1;
                        mutation_hotspot_parameters[mutation_hotspot][3] = Parameters.get_FLOAT(mutations_hotspot_Block, "Clock model Negative Binomial sucesses");
                        mutation_hotspot_parameters[mutation_hotspot][4] = Parameters.get_FLOAT(mutations_hotspot_Block, "Clock model Negative Binomial probability");

                        cout << "Negative Binomial sucesses: " << mutation_hotspot_parameters[mutation_hotspot][3] << endl;
                        cout << "Negative Binomial probability: " << mutation_hotspot_parameters[mutation_hotspot][4] << endl;
                    }
                    else if (clock_Model == "FIXED")
                    {
                        mutation_hotspot_parameters[mutation_hotspot][2] = 2;
                        mutation_hotspot_parameters[mutation_hotspot][3] = Parameters.get_FLOAT(mutations_hotspot_Block, "Clock model Fixed probability");

                        cout << "Fixed probability: " << mutation_hotspot_parameters[mutation_hotspot][3] << endl;
                    }
                    else
                    {
                        cout << "HOTSPOT " << mutation_hotspot + 1 << "'S CLOCK MODEL HAS TO BE POISSON, NEGATIVE BINOMIAL OR FIXED.\n";
                        exit(-1);
                    }

                    vector<string> bases = {"A", "T", "G", "C"};

                    cout << "\nConfiguring mutation probabilities:\n";

                    for (int reference_Base = 0; reference_Base < bases.size(); reference_Base++)
                    {
                        for (int mutation_Base = 0; mutation_Base < bases.size(); mutation_Base++)
                        {
                            string search_Query = bases[reference_Base] + bases[mutation_Base];
                            cout << search_Query << ": ";

                            if (reference_Base == 0)
                            {
                                A_0_mutation[mutation_hotspot][mutation_Base] = Parameters.get_FLOAT(mutations_hotspot_Block, search_Query);
                                cout << A_0_mutation[mutation_hotspot][mutation_Base];
                            }
                            else if (reference_Base == 1)
                            {
                                T_1_mutation[mutation_hotspot][mutation_Base] = Parameters.get_FLOAT(mutations_hotspot_Block, search_Query);
                                cout << T_1_mutation[mutation_hotspot][mutation_Base];
                            }
                            else if (reference_Base == 2)
                            {
                                G_2_mutation[mutation_hotspot][mutation_Base] = Parameters.get_FLOAT(mutations_hotspot_Block, search_Query);
                                cout << G_2_mutation[mutation_hotspot][mutation_Base];
                            }
                            else
                            {
                                C_3_mutation[mutation_hotspot][mutation_Base] = Parameters.get_FLOAT(mutations_hotspot_Block, search_Query);
                                cout << C_3_mutation[mutation_hotspot][mutation_Base];
                            }

                            if ((mutation_Base + 1) != bases.size())
                            {
                                cout << " | ";
                            }
                        }
                        cout << endl;
                    }
                }
            }
            else
            {
                cout << "No mutational hotspots present to configure.\n";
            }
        }
        else
        {
            // num_effect_Segregating_sites[0] = 0;
            num_effect_Segregating_sites[2] = 0;
        }

        if (mutation_recombination_proof_Reading_availability[1] == 1)
        {
            vector<pair<string, string>> recombinations_Block = Parameters.get_block_from_File(sequence_Master_location, "Recombination");

            recombination_Hotspots = Parameters.get_INT(recombinations_Block, "Number of hotspots");

            if (recombination_Hotspots > 0)
            {
                cout << "\nProcessing " << this->recombination_Hotspots << " recombination hotspots: \n";
                recombination_hotspot_parameters = functions.create_Fill_2D_array_FLOAT(recombination_Hotspots, 4, -1);

                string recombination_Hotspots_profiles_folder = Parameters.get_STRING(recombinations_Block, "Recombination profiles folder");
                cout << "Recombination profiles folder: " << recombination_Hotspots_profiles_folder << endl;

                vector<vector<pair<int, vector<float>>>> recom_probability_Changes;
                vector<vector<pair<int, vector<float>>>> recom_survivability_Changes;

                for (int recombination_hotspot = 0; recombination_hotspot < recombination_Hotspots; recombination_hotspot++)
                {
                    string hotspot_ID = "Hotspot " + to_string(recombination_hotspot + 1);
                    cout << "\nProcessing: " << hotspot_ID << endl;
                    vector<pair<string, string>> recombination_hotspot_Block = Parameters.get_block_from_block(recombinations_Block, hotspot_ID);

                    string region = Parameters.get_STRING(recombination_hotspot_Block, "Region");

                    vector<string> split_Region;
                    functions.split(split_Region, region, '_');

                    recombination_hotspot_parameters[recombination_hotspot][0] = stof(split_Region[0]);
                    recombination_hotspot_parameters[recombination_hotspot][1] = stof(split_Region[1]);
                    cout << "Region: From " << recombination_hotspot_parameters[recombination_hotspot][0] << " to " << recombination_hotspot_parameters[recombination_hotspot][1] << endl;

                    recombination_hotspot_parameters[recombination_hotspot][2] = Parameters.get_FLOAT(recombination_hotspot_Block, "Reference probability of recombination");
                    recombination_hotspot_parameters[recombination_hotspot][3] = Parameters.get_FLOAT(recombination_hotspot_Block, "Reference selectivity");

                    cout << "Reference probability of recombination: " << recombination_hotspot_parameters[recombination_hotspot][2] << endl;
                    cout << "Reference selectivity: " << recombination_hotspot_parameters[recombination_hotspot][3] << endl;

                    cout << endl;
                    if (recombination_Hotspots_profiles_folder != "NA")
                    {
                        recom_probability_Changes.push_back(Parameters.get_recombination_Hotspot_Parameters("probability", to_string(recombination_hotspot + 1), recombination_Hotspots_profiles_folder + "/hotspot_" + to_string(recombination_hotspot + 1) + ".csv", functions, tot_prob_selectivity[0]));
                        recom_survivability_Changes.push_back(Parameters.get_recombination_Hotspot_Parameters("selectivity", to_string(recombination_hotspot + 1), recombination_Hotspots_profiles_folder + "/hotspot_" + to_string(recombination_hotspot + 1) + ".csv", functions, tot_prob_selectivity[1]));
                    }
                    else
                    {
                        cout << "No recombination hotspot mutations present to configure.\n";
                    }

                    // if (recombination_hotspot == 1)
                    // {
                    //     vector<pair<int, vector<float>>> recom_matrix;
                    //     recom_matrix = recom_survivability_Changes[recombination_hotspot];

                    //     for (int check = 0; check < recom_matrix.size(); check++)
                    //     {
                    //         cout << recom_matrix[check].first << "\t";
                    //         for (size_t i = 0; i < recom_matrix[check].second.size(); i++)
                    //         {
                    //             cout << recom_matrix[check].second[i] << " | ";
                    //         }
                    //         cout << endl;
                    //     }
                    // }

                    // REMOVE
                    // exit(-1);
                }

                // cout << tot_prob_selectivity[0] << endl
                //      << tot_prob_selectivity[1] << endl;

                cout << endl;

                cout << "Configuring recombination hotspot mutations\n";

                recombination_Prob_matrix = functions.create_FLOAT_2D_arrays(tot_prob_selectivity[0], 5);
                recombination_Select_matrix = functions.create_FLOAT_2D_arrays(tot_prob_selectivity[1], 5);

                recombination_prob_Stride = (int *)malloc(sizeof(int) * (recombination_Hotspots + 1));
                recombination_select_Stride = (int *)malloc(sizeof(int) * (recombination_Hotspots + 1));

                int pos_curent_prob = 0;
                int pos_curent_select = 0;

                for (int recombination_hotspot = 0; recombination_hotspot < recombination_Hotspots; recombination_hotspot++)
                {
                    recombination_prob_Stride[recombination_hotspot] = pos_curent_prob;

                    vector<pair<int, vector<float>>> recom_matrix;
                    recom_matrix = recom_probability_Changes[recombination_hotspot];

                    sort(recom_matrix.begin(), recom_matrix.end());

                    for (int position = 0; position < recom_matrix.size(); position++)
                    {
                        recombination_Prob_matrix[pos_curent_prob][0] = recom_matrix[position].first;

                        for (int base = 0; base < 4; base++)
                        {
                            recombination_Prob_matrix[pos_curent_prob][base + 1] = recom_matrix[position].second[base];
                        }
                        pos_curent_prob++;
                    }

                    recombination_select_Stride[recombination_hotspot] = pos_curent_select;
                    recom_matrix = recom_survivability_Changes[recombination_hotspot];

                    sort(recom_matrix.begin(), recom_matrix.end());
                    for (int position = 0; position < recom_matrix.size(); position++)
                    {
                        recombination_Select_matrix[pos_curent_select][0] = recom_matrix[position].first;

                        for (int base = 0; base < 4; base++)
                        {
                            recombination_Select_matrix[pos_curent_select][base + 1] = recom_matrix[position].second[base];
                        }
                        pos_curent_select++;
                    }
                }

                recombination_prob_Stride[recombination_Hotspots] = pos_curent_prob;
                recombination_select_Stride[recombination_Hotspots] = pos_curent_select;

                if (tot_prob_selectivity[0] > 0)
                {
                    cout << "\nPrinting recombination mutations probability matrix:\n";
                    for (int recombination_hotspot = 0; recombination_hotspot < recombination_Hotspots; recombination_hotspot++)
                    {
                        cout << "\nHotspot " << recombination_hotspot + 1 << endl;
                        for (int stride = recombination_prob_Stride[recombination_hotspot]; stride < recombination_prob_Stride[recombination_hotspot + 1]; stride++)
                        {
                            for (int column = 0; column < 5; column++)
                            {
                                cout << recombination_Prob_matrix[stride][column] << "\t";
                            }
                            cout << endl;
                        }
                    }
                }
                if (tot_prob_selectivity[1] > 0)
                {
                    cout << "\nPrinting recombination mutations selectivity matrix:\n";
                    for (int recombination_hotspot = 0; recombination_hotspot < recombination_Hotspots; recombination_hotspot++)
                    {
                        cout << "\nHotspot " << recombination_hotspot + 1 << endl;
                        for (int stride = recombination_select_Stride[recombination_hotspot]; stride < recombination_select_Stride[recombination_hotspot + 1]; stride++)
                        {
                            for (int column = 0; column < 5; column++)
                            {
                                cout << recombination_Select_matrix[stride][column] << "\t";
                            }
                            cout << endl;
                        }
                    }
                }
            }
            else
            {
                cout << "No recombination hotspots present to configure.\n";
            }
        }
    }
    else
    {
        num_effect_Segregating_sites = (int *)malloc(sizeof(int) * 3);
        for (int fill = 0; fill < 3; fill++)
        {
            num_effect_Segregating_sites[fill] = 0;
        }
    }
    cout << endl;
    // exit(-1);
}

void simulator_Master::node_Master_Manager(functions_library &functions)
{
    parameter_load Parameters = parameter_load();
    cout << "Loading nodes master profile: " << this->node_Master_location << endl;

    vector<string> parameters_List = {
        "\"Shape replication time\"",
        "\"Scale replication time\"",
        "\"Shape days in host\"",
        "\"Scale days in host\"",
        "\"Sampling present\"",
        "\"Progeny distribution type\"",
        "\"Number of node profiles\"",
        "\"Location of node profiles\"",
        "\"Infected to Recovered\""};

    vector<string> found_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);
    parameters_List.clear();

    random_device rd;
    mt19937 gen(rd());

    gamma_distribution<float> days_in_Host(Parameters.get_FLOAT(found_Parameters[0]), Parameters.get_FLOAT(found_Parameters[1]));
    this->generation_Time = days_in_Host(gen);

    if (this->generation_Time > 0)
    {
        cout << "\n";
        cout << "Time taken for viral replication (from cell attachment to release) (days): " << generation_Time << endl;

        this->shape_days_in_Host = Parameters.get_FLOAT(found_Parameters[2]);
        this->scale_days_in_Host = Parameters.get_FLOAT(found_Parameters[3]);
        cout << "Days in host gamma distribution: Shape: " << this->shape_days_in_Host << " Scale: " << this->scale_days_in_Host << endl;

        if (functions.to_Upper_Case(Parameters.get_STRING(found_Parameters[8])) == "YES")
        {
            infected_to_Recovered = "YES";
        }

        cout << "SIRS status: " << infected_to_Recovered << endl;

        cout << "\nSampling status: ";
        if (functions.to_Upper_Case(Parameters.get_STRING(found_Parameters[4])) == "YES")
        {
            cout << "Active\n";
            parameters_List = {"\"Sampling rate Binomial trials\"",
                               "\"Sampling rate Binomial probability\"",
                               "\"Resampling of nodes\"",
                               "\"Distribution per node sequences sampled\"",
                               "\"Limit samples obtained\""};

            vector<string> sampling_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);

            // sampling_effect = Parameters.get_STRING(sampling_Parameters[0]);
            // transform(sampling_effect.begin(), sampling_effect.end(), sampling_effect.begin(), ::toupper);

            // cout << "Effect of sampling: " << sampling_effect << endl;

            trials_Sampling = 1;
            sampled_sequences_Folder = output_Network_location + "/sampled_Sequences";
            functions.config_Folder(sampled_sequences_Folder, "Sampled sequences");
            functions.create_File(sampled_sequences_Folder + "/sampled_Sequences_FASTA.nfasta");
            functions.create_File(sampled_sequences_Folder + "/sampled_Sequences_summary.csv", "Host_ID\tTissue\tSequence_ID\tSampling_time\t_Sampling_date");
            sampling_trials = Parameters.get_INT(sampling_Parameters[0]);
            sampling_probability = Parameters.get_FLOAT(sampling_Parameters[1]);

            cout << "Sampling rate distribution trials: " << sampling_trials << endl;
            cout << "Sampling rate distribution probability: " << sampling_probability << endl;

            // string resampling_Status = Parameters.get_STRING(sampling_Parameters[2]);
            // transform(resampling_Status.begin(), resampling_Status.end(), resampling_Status.begin(), ::toupper);

            cout << "Resampling of the same nodes: ";
            if (functions.to_Upper_Case(Parameters.get_STRING(sampling_Parameters[2])) == "YES")
            {
                cout << "Active\n";
                resampling = 1;
            }
            else
            {
                cout << "Inactive\n";
            }

            per_Node_sampling = (float *)malloc(sizeof(float) * 3);

            cout << "Sequences sampled per node\nDistribution type: ";
            vector<string> sequence_sampling_Parameters;

            if (functions.to_Upper_Case(Parameters.get_STRING(sampling_Parameters[3])) == "FIXED")
            {
                cout << "Fixed number of sequences sampled\n";
                per_Node_sampling[0] = 0;
                parameters_List = {"\"Per node sampling rate Fixed\""};
                sequence_sampling_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);
                per_Node_sampling[1] = Parameters.get_INT(sequence_sampling_Parameters[0]);
                cout << "Sequences sampled per node: " << per_Node_sampling[1] << endl;
            }
            else if (functions.to_Upper_Case(Parameters.get_STRING(sampling_Parameters[3])) == "BINOMIAL")
            {
                cout << "Binomial distribution\n";
                per_Node_sampling[0] = 1;
                parameters_List = {"\"Per node sampling rate Binomial trials\"",
                                   "\"Per node sampling rate Binomial probability\""};
                sequence_sampling_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);
                per_Node_sampling[1] = Parameters.get_INT(sequence_sampling_Parameters[0]);
                per_Node_sampling[2] = Parameters.get_FLOAT(sequence_sampling_Parameters[1]);

                cout << "Per node sampling rate Binomial trials: " << per_Node_sampling[1] << endl;
                cout << "Per node sampling rate Binomial probability: " << per_Node_sampling[2] << endl;
            }

            cout << "Limited number of sequences to be sampled: ";
            if (functions.to_Upper_Case(Parameters.get_STRING(sampling_Parameters[4])) == "YES")
            {
                cout << "YES\n";
                cout << "Sequenced limit: ";
                parameters_List = {"\"Max samples to obtain\""};
                vector<string> sequence_sampling_limit = Parameters.get_parameters(node_Master_location, parameters_List);

                limit_Sampled = Parameters.get_INT(sequence_sampling_limit[0]);
                cout << limit_Sampled << endl;
            }
            else
            {
                cout << "Unlimited sampling\n";
            }
        }
        else
        {
            cout << "Inactive\n";
        }

        // exit(-1);

        progeny_distribution_Model = Parameters.get_STRING(found_Parameters[5]);
        transform(progeny_distribution_Model.begin(), progeny_distribution_Model.end(), progeny_distribution_Model.begin(), ::toupper);
        progeny_distribution_parameters_Array = (float *)malloc(sizeof(float) * 3);

        vector<string> progeny_distribution_Parameters;

        if (progeny_distribution_Model == "NEGATIVE BINOMIAL")
        {
            progeny_distribution_parameters_Array[0] = 0;
            parameters_List = {"\"Progeny Negative Binomial sucesses\"",
                               "\"Progeny Negative Binomial probability\""};

            progeny_distribution_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);

            progeny_distribution_parameters_Array[1] = Parameters.get_INT(progeny_distribution_Parameters[0]);
            progeny_distribution_parameters_Array[2] = Parameters.get_FLOAT(progeny_distribution_Parameters[1]);

            cout << "\nProgeny generation via a Negative Binomial disribution\n";
            cout << "Negative Binomial sucesses: " << progeny_distribution_parameters_Array[1]
                 << "\nNegative Binomial probability: " << progeny_distribution_parameters_Array[2] << endl;
        }
        else if (progeny_distribution_Model == "GAMMA")
        {
            progeny_distribution_parameters_Array[0] = 1;
            parameters_List = {"\"Progeny Gamma shape\"",
                               "\"Progeny Gamma scale\""};

            progeny_distribution_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);
            progeny_distribution_parameters_Array[1] = Parameters.get_FLOAT(progeny_distribution_Parameters[0]);
            progeny_distribution_parameters_Array[2] = Parameters.get_FLOAT(progeny_distribution_Parameters[1]);

            cout << "\nProgeny generation via a Gamma disribution\n";
            cout << "Gamma shape: " << progeny_distribution_parameters_Array[1]
                 << "\nGamma scale: " << progeny_distribution_parameters_Array[2] << endl;
        }
        else if (progeny_distribution_Model == "POISSON")
        {
            progeny_distribution_parameters_Array[0] = 2;
            progeny_distribution_parameters_Array[2] = -1;

            parameters_List = {"\"Progeny Poisson mean\""};

            progeny_distribution_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);
            progeny_distribution_parameters_Array[1] = Parameters.get_FLOAT(progeny_distribution_Parameters[0]);

            cout << "\nProgeny generation via a Poisson disribution\n";
            cout << "Poisson mean: " << progeny_distribution_parameters_Array[1] << endl;
        }
        else
        {
            cout << "\nERROR: PROGENY GENERATION DISTRIBUTION MODEL TYPE\n";
            exit(-1);
        }
        progeny_distribution_Parameters.clear();

        // exit(-1);

        // Configure tissue profile read
        vector<pair<string, string>> Tissue_profiles_block_Data = Parameters.get_block_from_File(node_Master_location, "Tissue profiles");

        // for (size_t i = 0; i < Tissue_profiles_block_Data.size(); i++)
        // {
        //     cout << Tissue_profiles_block_Data[i].first << " : " << Tissue_profiles_block_Data[i].second << endl;
        // }

        num_tissues_per_Node = Parameters.get_INT(Tissue_profiles_block_Data, "Number of tissues");

        if (num_tissues_per_Node > 0)
        {
            cout << "\nNumber of tissues in a node: " << num_tissues_per_Node << endl;

            for (int tissue = 0; tissue < num_tissues_per_Node; tissue++)
            {
                string check_Tissue = "Tissue " + to_string(tissue + 1) + " Name";
                tissue_Names.push_back(Parameters.get_STRING(Tissue_profiles_block_Data, check_Tissue));
                cout << check_Tissue << ": " << tissue_Names[tissue] << endl;
            }
        }
        else
        {
            cout << "ERROR: TISSUE NUMBER HAS TO BE GREATER THAN ZERO.\n\n";
        }

        cout << endl;

        string viral_Entry = Parameters.get_STRING(Tissue_profiles_block_Data, "Viral entry tissues");
        string viral_Infectious = Parameters.get_STRING(Tissue_profiles_block_Data, "Infectious load tissues");
        string viral_Terminal = Parameters.get_STRING(Tissue_profiles_block_Data, "Terminal load tissues");
        string viral_Exit = Parameters.get_STRING(Tissue_profiles_block_Data, "Viral exit tissues");
        string sampling_Tissues;
        if (trials_Sampling != -1)
        {
            sampling_Tissues = Parameters.get_STRING(Tissue_profiles_block_Data, "Sampling tissues");
        }

        vector<string> viral_tissue_Split;

        functions.split(viral_tissue_Split, viral_Entry, ',');
        this->entry_tissues = viral_tissue_Split.size();
        entry_array = (int *)malloc(entry_tissues * sizeof(int));
        cout << this->entry_tissues << " tissue(s) availble for entry: ";
        for (size_t i = 0; i < entry_tissues; i++)
        {
            entry_array[i] = stoi(viral_tissue_Split[i]) - 1;
            cout << tissue_Names[entry_array[i]];

            if ((i + 1) != entry_tissues)
            {
                cout << ", ";
            }
            else
            {
                cout << endl;
            }
        }

        functions.split(viral_tissue_Split, viral_Infectious, ',');
        this->infectious_tissues = viral_tissue_Split.size();
        this->infectious_array = (int *)malloc(infectious_tissues * sizeof(int));
        cout << this->infectious_tissues << " tissue(s) determine node infectivity: ";
        for (size_t i = 0; i < infectious_tissues; i++)
        {
            infectious_array[i] = stoi(viral_tissue_Split[i]) - 1;
            cout << tissue_Names[infectious_array[i]];

            if ((i + 1) != infectious_tissues)
            {
                cout << ", ";
            }
            else
            {
                cout << endl;
            }
        }

        functions.split(viral_tissue_Split, viral_Terminal, ',');
        this->terminal_tissues = viral_tissue_Split.size();
        this->terminal_array = (int *)malloc(terminal_tissues * sizeof(int));
        cout << this->terminal_tissues << " tissue(s) determine node termination (death): ";
        for (size_t i = 0; i < terminal_tissues; i++)
        {
            terminal_array[i] = stoi(viral_tissue_Split[i]) - 1;
            cout << tissue_Names[terminal_array[i]];

            if ((i + 1) != terminal_tissues)
            {
                cout << ", ";
            }
            else
            {
                cout << endl;
            }
        }

        functions.split(viral_tissue_Split, viral_Exit, ',');
        this->exit_tissues = viral_tissue_Split.size();
        this->exit_array = (int *)malloc(exit_tissues * sizeof(int));
        cout << this->exit_tissues << " tissue(s) available for viral exit: ";
        for (size_t i = 0; i < exit_tissues; i++)
        {
            exit_array[i] = stoi(viral_tissue_Split[i]) - 1;
            cout << tissue_Names[exit_array[i]];

            if ((i + 1) != exit_tissues)
            {
                cout << ", ";
            }
            else
            {
                cout << endl;
            }
        }

        if (trials_Sampling != -1)
        {
            functions.split(viral_tissue_Split, sampling_Tissues, ',');
            this->sampling_tissues = viral_tissue_Split.size();
            this->sampling_array = (int *)malloc(sampling_tissues * sizeof(int));
            cout << this->sampling_tissues << " tissue(s) available for host sampling: ";
            for (size_t i = 0; i < sampling_tissues; i++)
            {
                sampling_array[i] = stoi(viral_tissue_Split[i]) - 1;
                cout << tissue_Names[sampling_array[i]];

                if ((i + 1) != sampling_tissues)
                {
                    cout << ", ";
                }
                else
                {
                    cout << endl;
                }
            }
        }

        cout << endl;
        // exit(-1);

        // cout << "******\n\n";
        viral_Migration = Parameters.get_STRING(Tissue_profiles_block_Data, "Viral tissue migration");
        transform(viral_Migration.begin(), viral_Migration.end(), viral_Migration.begin(), ::toupper);

        if (viral_Migration == "YES")
        {
            cout << "Viral migration: Activated\n";

            viral_Migration_Values = functions.create_Fill_2D_array_FLOAT(num_tissues_per_Node * (num_tissues_per_Node - 1), 2, -1);
            migration_start_Generation = (int *)malloc(num_tissues_per_Node * (num_tissues_per_Node - 1) * sizeof(int));

            for (int fill_mig = 0; fill_mig < num_tissues_per_Node * (num_tissues_per_Node - 1); fill_mig++)
            {
                migration_start_Generation[fill_mig] = -1;
            }

            vector<pair<string, string>> Viral_migration_block_Data = Parameters.get_block_from_block(Tissue_profiles_block_Data, "Viral migration");

            cout << "Configuring tissue to tissue migrations:\n\n";

            for (int migration_Check = 0; migration_Check < (num_tissues_per_Node * (num_tissues_per_Node - 1)); migration_Check++)
            {
                int source = migration_Check / (num_tissues_per_Node - 1);
                int destination = migration_Check % (num_tissues_per_Node - 1);

                if (destination >= source)
                {
                    destination = destination + 1;
                }

                string check_source_destination = to_string(source + 1) + "_" + to_string(destination + 1);

                vector<pair<string, string>> block_Migration = Parameters.check_block_from_block(Viral_migration_block_Data, check_source_destination);

                if (block_Migration.size() > 0)
                {
                    cout << "From " << tissue_Names[source] << " to " << tissue_Names[destination] << endl;
                    for (int i = 0; i < block_Migration.size(); i++)
                    {
                        if (Parameters.get_STRING(block_Migration[i].first) == "Cell migration Binomial trials")
                        {
                            viral_Migration_Values[migration_Check][0] = Parameters.get_INT(block_Migration[i].second);
                        }
                        else if (Parameters.get_STRING(block_Migration[i].first) == "Cell migration Binomial probability")
                        {
                            viral_Migration_Values[migration_Check][1] = Parameters.get_FLOAT(block_Migration[i].second);
                        }
                        else if (Parameters.get_STRING(block_Migration[i].first) == "Start generation")
                        {
                            migration_start_Generation[migration_Check] = Parameters.get_INT(block_Migration[i].second);
                        }
                        else
                        {
                            cout << "ERROR INVALID ENTRY AT " << check_source_destination << endl;
                            exit(-1);
                        }
                    }
                    cout << "Migration start generation: " << migration_start_Generation[migration_Check] << endl;
                    cout << "Cell migration Binomial trials: " << viral_Migration_Values[migration_Check][0] << endl;
                    cout << "Cell migration Binomial probability: " << viral_Migration_Values[migration_Check][1] << endl;
                    cout << endl;
                    // exit(-1);
                }
            }
        }
        else
        {
            cout << "Viral migration: Does not occur\n";
        }
        // for (size_t i = 0; i < Viral_migration_block_Data.size(); i++)
        // {
        //     cout << Viral_migration_block_Data[i].first << " : " << Viral_migration_block_Data[i].second << endl;
        // }

        number_of_node_Profiles = Parameters.get_INT(found_Parameters[6]);
        if (number_of_node_Profiles > 0)
        {
            cout << number_of_node_Profiles << " node profiles present" << endl;

            vector<vector<int>> replication_phases_Profile_tissues;

            vector<float> time_Ratios;
            vector<string> phase_Type;
            vector<pair<float, float>> phase_paramaters;

            // (int *)malloc(infectious_tissues * sizeof(int));
            node_profile_Distributions = (float *)malloc(sizeof(float) * number_of_node_Profiles);
            profile_tissue_Limits = functions.create_Fill_2D_array_FLOAT(num_tissues_per_Node * number_of_node_Profiles, 3, 0);
            node_sampling_effect = functions.create_Fill_2D_array_FLOAT(number_of_node_Profiles, 3, 0);

            tissue_param_profile_Stride = (int *)malloc(sizeof(int) * (number_of_node_Profiles + 1));
            tissue_param_profile_Stride[0] = 0;

            if (infectious_tissues > 0)
            {
                infectious_load_Profiles_param = functions.create_Fill_2D_array_FLOAT(number_of_node_Profiles, 3, -1);
            }

            if (terminal_tissues > 0)
            {
                terminal_load_Profiles_param = functions.create_Fill_2D_array_FLOAT(number_of_node_Profiles, 3, -1);
            }

            node_Profile_folder_Location = Parameters.get_STRING(found_Parameters[7]);
            cout << "Extracting profiles from: " << node_Profile_folder_Location << endl
                 << endl;

            infection_parameters = functions.create_Fill_2D_array_FLOAT(number_of_node_Profiles, 3, -1);

            cell_Distribution_Type = functions.create_INT_2D_arrays(number_of_node_Profiles, num_tissues_per_Node);

            for (int profile = 0; profile < number_of_node_Profiles; profile++)
            {
                parameters_List = {"\"Profile name\"",
                                   "\"Probability of occurence\"",
                                   "\"Infection sequence distribution\""};

                string profile_check = node_Profile_folder_Location + "/profile_" + to_string(profile + 1) + ".json";

                if (filesystem::exists(profile_check))
                {
                    cout << profile + 1 << " Profile found: " << profile_check << endl;

                    vector<string> profile_Parameters = Parameters.get_parameters(profile_check, parameters_List);
                    profile_names.push_back(Parameters.get_STRING(profile_Parameters[0]));
                    node_profile_Distributions[profile] = Parameters.get_FLOAT(profile_Parameters[1]);
                    vector<int> replication_phases_tissues;

                    cout << "Configuring profile: " << profile_names[profile] << endl;
                    cout << "Probability of occurence: " << node_profile_Distributions[profile] << endl;

                    cout << "\nConfiguration of host to host infection\n";
                    vector<string> infected_Parameters;
                    cout << "Infection distribution: ";
                    if (functions.to_Upper_Case(Parameters.get_STRING(profile_Parameters[2])) == "FIXED")
                    {
                        infection_parameters[profile][0] = 0;
                        cout << "FIXED\n";
                        parameters_List = {"\"Infection sequence Fixed\""};
                        infected_Parameters = Parameters.get_parameters(profile_check, parameters_List);
                        infection_parameters[profile][1] = Parameters.get_INT(infected_Parameters[0]);
                        cout << "Infection sequence Fixed: " << infection_parameters[profile][1] << endl;
                    }
                    else if (functions.to_Upper_Case(Parameters.get_STRING(profile_Parameters[2])) == "BINOMIAL")
                    {
                        infection_parameters[profile][0] = 1;
                        cout << "BINOMIAL\n";
                        parameters_List = {"\"Infection sequence Binomial trials\"",
                                           "\"Infection sequence Binomial probability\""};
                        infected_Parameters = Parameters.get_parameters(profile_check, parameters_List);
                        infection_parameters[profile][1] = Parameters.get_INT(infected_Parameters[0]);
                        infection_parameters[profile][2] = Parameters.get_FLOAT(infected_Parameters[1]);
                        cout << "Infection sequence Binomial trials: " << infection_parameters[profile][1] << endl;
                        cout << "Infection sequence Binomial probability: " << infection_parameters[profile][2] << endl;
                    }
                    else
                    {
                        cout << "ERROR Infection sequence distribution HAS TO BE EITHER \"FIXED\" OR \"NEGATIVE BINOMIAL\"";
                        exit(1);
                    }
                    cout << endl;

                    // exit(-1);

                    // Configure tisses and their phases
                    if (trials_Sampling != -1)
                    {
                        parameters_List = {"\"Sampling effect\""};
                        vector<string> sampling_Parameters = Parameters.get_parameters(profile_check, parameters_List);

                        string sampling_effect = Parameters.get_STRING(sampling_Parameters[0]);
                        transform(sampling_effect.begin(), sampling_effect.end(), sampling_effect.begin(), ::toupper);

                        if (sampling_effect == "NO CHANGE" || sampling_effect == "REMOVED" || sampling_effect == "LESS INFECTIOUS")
                        {
                            if (sampling_effect == "REMOVED")
                            {
                                node_sampling_effect[profile][0] = 1;
                            }
                            else if (sampling_effect == "LESS INFECTIOUS")
                            {
                                node_sampling_effect[profile][0] = -1;
                                parameters_List = {"\"Sampling less infectious effect Alpha\"",
                                                   "\"Sampling less infectious effect Beta\""};
                                vector<string> profile_sampling_Less = Parameters.get_parameters(profile_check, parameters_List);
                                node_sampling_effect[profile][1] = Parameters.get_FLOAT(profile_sampling_Less[0]);
                                node_sampling_effect[profile][2] = Parameters.get_FLOAT(profile_sampling_Less[1]);
                            }
                        }
                        else
                        {
                            cout << "ERROR: PROFILE " << profile + 1 << " SAMPLING EFFECT: " << sampling_effect << " IS NOT A VALID ENTRY." << endl;
                            cout << "IT HAS TO BE ONE OF THE FOLLOWING: NO CHANGE, REMOVED OR LESS INFECTIOUS.\n";
                            exit(-1);
                        }

                        cout << "Sampling effect: " << sampling_effect << endl;
                        if (node_sampling_effect[profile][0] == -1)
                        {
                            cout << "Beta distribution alpha value: " << node_sampling_effect[profile][1] << endl;
                            cout << "Beta distribution beta value: " << node_sampling_effect[profile][2] << endl;
                        }
                    }

                    if (infectious_tissues > 0)
                    {
                        parameters_List = {"\"Infectious load distribution type\""};

                        vector<string> infectious_Tissue_distribution = Parameters.get_parameters(profile_check, parameters_List);
                        transform(infectious_Tissue_distribution[0].begin(), infectious_Tissue_distribution[0].end(), infectious_Tissue_distribution[0].begin(), ::toupper);

                        cout << "\nInfectious viral load distribution: " << infectious_Tissue_distribution[0] << endl;

                        if (infectious_Tissue_distribution[0] == "\"BINOMIAL\"")
                        {
                            parameters_List = {"\"Infectious load Binomial trials\"",
                                               "\"Infectious load Binomial probability\""};
                            infectious_Tissue_distribution = Parameters.get_parameters(profile_check, parameters_List);

                            infectious_load_Profiles_param[profile][0] = 0;
                            infectious_load_Profiles_param[profile][1] = Parameters.get_INT(infectious_Tissue_distribution[0]);
                            infectious_load_Profiles_param[profile][2] = Parameters.get_FLOAT(infectious_Tissue_distribution[1]);

                            cout << "Infectious load Binomial trials: " << infectious_load_Profiles_param[profile][1] << endl;
                            cout << "Infectious load Binomial probability: " << infectious_load_Profiles_param[profile][2] << endl;
                        }
                        else if (infectious_Tissue_distribution[0] == "\"FIXED\"")
                        {
                            parameters_List = {"\"Infectious load Fixed\""};
                            infectious_Tissue_distribution = Parameters.get_parameters(profile_check, parameters_List);

                            infectious_load_Profiles_param[profile][0] = -1;
                            infectious_load_Profiles_param[profile][1] = Parameters.get_INT(infectious_Tissue_distribution[0]);

                            cout << "Infectious load Fixed: " << infectious_load_Profiles_param[profile][1] << endl;
                        }
                        else
                        {
                            cout << "PROFILE " << profile + 1 << " ERROR INFECTIOUS TISSUE DISTRIBUTION SHOULD BE FIXED OR BINOMIAL.\n";
                            exit(-1);
                        }
                    }

                    if (terminal_tissues > 0)
                    {
                        parameters_List = {"\"Terminal load distribution type\""};

                        vector<string> terminal_Tissue_distribution = Parameters.get_parameters(profile_check, parameters_List);
                        transform(terminal_Tissue_distribution[0].begin(), terminal_Tissue_distribution[0].end(), terminal_Tissue_distribution[0].begin(), ::toupper);

                        cout << "\nTerminal viral load distribution: " << terminal_Tissue_distribution[0] << endl;

                        if (terminal_Tissue_distribution[0] == "\"BINOMIAL\"")
                        {
                            parameters_List = {"\"Terminal load Binomial trials\"",
                                               "\"Terminal load Binomial probability\""};
                            terminal_Tissue_distribution = Parameters.get_parameters(profile_check, parameters_List);

                            terminal_load_Profiles_param[profile][0] = 0;
                            terminal_load_Profiles_param[profile][1] = Parameters.get_INT(terminal_Tissue_distribution[0]);
                            terminal_load_Profiles_param[profile][2] = Parameters.get_FLOAT(terminal_Tissue_distribution[1]);

                            cout << "Terminal load Binomial trials: " << terminal_load_Profiles_param[profile][1] << endl;
                            cout << "Terminal load Binomial probability: " << terminal_load_Profiles_param[profile][2] << endl;
                        }
                        else if (terminal_Tissue_distribution[0] == "\"FIXED\"")
                        {
                            parameters_List = {"\"Terminal load Fixed\""};
                            terminal_Tissue_distribution = Parameters.get_parameters(profile_check, parameters_List);

                            terminal_load_Profiles_param[profile][0] = -1;
                            terminal_load_Profiles_param[profile][1] = Parameters.get_INT(terminal_Tissue_distribution[0]);

                            cout << "Terminal load Fixed: " << terminal_load_Profiles_param[profile][1] << endl;
                        }
                        else
                        {
                            cout << "PROFILE " << profile + 1 << " ERROR TERMINAL TISSUE DISTRIBUTION SHOULD BE FIXED OR BINOMIAL.\n";
                            exit(-1);
                        }
                    }

                    // exit(-1);

                    cout << "\nCollecting tissue data\n";
                    vector<pair<string, string>> Tissue_profiles_block_Data = Parameters.get_block_from_File(profile_check, "Tissue profiles");

                    // for (int test = 0; test < Tissue_profiles_block_Data.size(); test++)
                    // {
                    //     cout << Tissue_profiles_block_Data[test].first << endl;
                    // }

                    int tissue_Limit_Start = profile * num_tissues_per_Node;

                    vector<pair<float, float>> Tissue_cell_disribution;

                    for (int tissue = 0; tissue < num_tissues_per_Node; tissue++)
                    {
                        string get_Tissue = "Tissue " + to_string(tissue + 1);
                        cout << "Processing: " << get_Tissue << endl;

                        vector<pair<string, string>> current_tissue_Profile_block_Data = Parameters.get_block_from_block(Tissue_profiles_block_Data, get_Tissue);

                        string cell_Limit = Parameters.get_STRING(current_tissue_Profile_block_Data, "Cell limit");
                        transform(cell_Limit.begin(), cell_Limit.end(), cell_Limit.begin(), ::toupper);

                        cout << tissue_Names[tissue] << " tissue cell limit: ";

                        if (cell_Limit == "YES")
                        {
                            cout << "YES\n";
                            profile_tissue_Limits[tissue_Limit_Start + tissue][0] = 1;
                            profile_tissue_Limits[tissue_Limit_Start + tissue][1] = Parameters.get_INT(current_tissue_Profile_block_Data, "Cell limit Binomial trials");
                            profile_tissue_Limits[tissue_Limit_Start + tissue][2] = Parameters.get_FLOAT(current_tissue_Profile_block_Data, "Cell limit Binomial probability");

                            cout << "Cell limit Binomial trials: " << profile_tissue_Limits[tissue_Limit_Start + tissue][1] << endl
                                 << "Cell limit Binomial probability: " << profile_tissue_Limits[tissue_Limit_Start + tissue][2] << endl;
                        }
                        else
                        {
                            cout << " NO\n";
                        }
                        cout << endl;

                        string viral_distribution_Type = Parameters.get_STRING(current_tissue_Profile_block_Data, "Viral distribution type");
                        transform(viral_distribution_Type.begin(), viral_distribution_Type.end(), viral_distribution_Type.begin(), ::toupper);

                        cout << "Viral distribution type: ";
                        if (viral_distribution_Type == "BINOMIAL")
                        {
                            cell_Distribution_Type[profile][tissue] = 0;
                            cout << viral_distribution_Type << endl;

                            int viral_distribution_Binomial_trials = Parameters.get_INT(current_tissue_Profile_block_Data, "Viral distribution Binomial trials");
                            float viral_distribution_Binomial_probability = Parameters.get_FLOAT(current_tissue_Profile_block_Data, "Viral distribution Binomial probability");

                            Tissue_cell_disribution.push_back(make_pair(viral_distribution_Binomial_trials, viral_distribution_Binomial_probability));

                            cout << "Viral distribution Binomial trials: " << viral_distribution_Binomial_trials << endl;
                            cout << "Viral distribution Binomial probability: " << viral_distribution_Binomial_probability << endl;
                        }
                        else if (viral_distribution_Type == "GAMMA")
                        {
                            cell_Distribution_Type[profile][tissue] = 1;
                            cout << viral_distribution_Type << endl;

                            float viral_distribution_Gamma_shape = Parameters.get_FLOAT(current_tissue_Profile_block_Data, "Viral distribution Gamma shape");
                            float viral_distribution_Gamma_scale = Parameters.get_FLOAT(current_tissue_Profile_block_Data, "Viral distribution Gamma scale");

                            Tissue_cell_disribution.push_back(make_pair(viral_distribution_Gamma_shape, viral_distribution_Gamma_scale));

                            cout << "Viral distribution Gamma shape: " << viral_distribution_Gamma_shape << endl;
                            cout << "Viral distribution Gamma scale: " << viral_distribution_Gamma_scale << endl;
                        }
                        else
                        {
                            cout << "ERROR: VIRAL DISTRIBUTION TYPE HAS TO BE BINOMIAL OR GAMMA\n";
                            exit(-1);
                        }

                        cout << endl;

                        // exit(-1);

                        cout << "Configuring Tissue " << tissue + 1 << " replication phases: \n";
                        vector<pair<string, string>> replication_Phases_Block = Parameters.get_block_from_block(current_tissue_Profile_block_Data, "Replication phases");
                        int num_Phases = Parameters.get_INT(replication_Phases_Block, "Number of phases");

                        if (num_Phases > 0)
                        {
                            replication_phases_tissues.push_back(num_Phases);
                            cout << "Number of phases: " << replication_phases_tissues[replication_phases_tissues.size() - 1] << endl;

                            for (int rep_Phase = 0; rep_Phase < num_Phases; rep_Phase++)
                            {
                                string phase_keyword = "Phase " + to_string(rep_Phase + 1);
                                string phase_Mode = phase_keyword + " Mode";
                                string phase_Time_ratio = phase_keyword + " Time ratio";

                                phase_Type.push_back(Parameters.get_STRING(replication_Phases_Block, phase_Mode));
                                time_Ratios.push_back(Parameters.get_FLOAT(replication_Phases_Block, phase_Time_ratio));

                                transform(phase_Type[phase_Type.size() - 1].begin(), phase_Type[phase_Type.size() - 1].end(), phase_Type[phase_Type.size() - 1].begin(), ::toupper);

                                if (phase_Type[phase_Type.size() - 1] == "NEUTRAL")
                                {
                                    phase_paramaters.push_back(make_pair(-1, -1));
                                }
                                else if (phase_Type[phase_Type.size() - 1] == "STATIONARY")
                                {
                                    phase_paramaters.push_back(make_pair(Parameters.get_FLOAT(replication_Phases_Block, phase_keyword + " Variance"), -1));
                                }
                                else if (phase_Type[phase_Type.size() - 1] == "DEPRICIATION")
                                {
                                    phase_paramaters.push_back(make_pair(Parameters.get_FLOAT(replication_Phases_Block, phase_keyword + " Alpha"), Parameters.get_FLOAT(replication_Phases_Block, phase_keyword + " Beta")));
                                }
                                else
                                {
                                    cout << "ERROR " << profile + 1 << " PROFILE TISSUE: " << tissue + 1 << "TISSUE REPLICATION MODE HAS TO BE ONE OF NEUTRAL, STATIONARY OT DEPRICIATION.\n";
                                    exit(-1);
                                }

                                cout << "\nPhase " << rep_Phase + 1 << ": \n";
                                cout << "Mode: " << phase_Type[phase_Type.size() - 1] << endl;
                                // CHECK TIME RATIO ADDITIONS
                                cout << "Time ratio: " << time_Ratios[time_Ratios.size() - 1] << endl;

                                if (phase_paramaters[phase_paramaters.size() - 1].first != -1 && phase_paramaters[phase_paramaters.size() - 1].second == -1)
                                {
                                    cout << "Variance of Stationary: " << phase_paramaters[phase_paramaters.size() - 1].first << endl;
                                }
                                else if (phase_paramaters[phase_paramaters.size() - 1].first != -1 && phase_paramaters[phase_paramaters.size() - 1].second != -1)
                                {
                                    cout << "Alpha of Depriciation: " << phase_paramaters[phase_paramaters.size() - 1].first;
                                    cout << " ; Beta of Depriciation: " << phase_paramaters[phase_paramaters.size() - 1].second << endl;
                                }
                            }
                        }
                        else
                        {
                            cout << "ERROR: " << profile + 1 << " PROFILE TISSUE: " << tissue + 1 << " CANNOT HAVE LESS THAN ONE PHASES." << endl;
                            exit(-1);
                        }
                        cout << endl;
                    }

                    ALL_profiles_Tissue_cell_disribution.push_back(Tissue_cell_disribution);

                    // exit(-1);

                    replication_phases_Profile_tissues.push_back(replication_phases_tissues);
                    tissue_param_profile_Stride[profile + 1] = phase_Type.size();

                    cout << endl;
                }
                else
                {
                    cout << "ERROR: " << profile + 1 << " PROFILE NOT FOUND, PLEASE CHECK LOCATION: " << profile_check << endl;
                    exit(-1);
                }
            }

            // Finally check profiles distributions
            float distribution_Check = 0;
            cout << "Final configuration of profiles\n";
            cout << "Performing profile distribution check: ";
            for (int profile = 0; profile < number_of_node_Profiles; profile++)
            {
                distribution_Check = distribution_Check + node_profile_Distributions[profile];
            }
            if (distribution_Check != 1)
            {
                cout << "ERROR: PROFILE DISTRIBUTIONS MUST SUM UPTO 1, NOW THEY SUM UPTO: " << distribution_Check << endl;
                exit(-1);
            }
            else
            {
                cout << "\nDistribution check passed\n";
            }

            cout << endl;

            cout << "Performing profile tissue configurations and checks: \n";
            num_replication_phases = (int *)malloc(sizeof(int) * number_of_node_Profiles * num_tissues_per_Node);
            tissue_replication_data = functions.create_Fill_2D_array_FLOAT(time_Ratios.size(), 4, -1);

            for (int profile = 0; profile < number_of_node_Profiles; profile++)
            {
                cout << "\nConfiguring Profile " << profile + 1 << endl;
                vector<int> phases_per_Tissue = replication_phases_Profile_tissues[profile];

                for (int tissue = 0; tissue < phases_per_Tissue.size(); tissue++)
                {
                    num_replication_phases[(profile * num_tissues_per_Node) + tissue] = phases_per_Tissue[tissue];
                }

                int tissue = 0;
                int num_phases_per_tissue = 0;

                float time_Check = 0;

                for (int param_Index = tissue_param_profile_Stride[profile]; param_Index < tissue_param_profile_Stride[profile + 1]; param_Index++)
                {
                    tissue_replication_data[param_Index][0] = time_Ratios[param_Index];

                    time_Check = time_Check + tissue_replication_data[param_Index][0];

                    num_phases_per_tissue++;

                    if (num_phases_per_tissue == num_replication_phases[(profile * num_tissues_per_Node) + tissue])
                    {
                        num_phases_per_tissue = 0;
                        tissue++;
                        if (time_Check == 1)
                        {
                            cout << "Tissue " << tissue << " passed\n";
                        }
                        else
                        {
                            cout << "ERROR TISSUE " << tissue << " FAILED, PHASES SHOULD ADD TO 1 NOT " << time_Check << "\n";
                            exit(-1);
                        }
                        time_Check = 0;
                    }

                    if (phase_Type[param_Index] == "NEUTRAL")
                    {
                        tissue_replication_data[param_Index][1] = 0;
                    }
                    else if (phase_Type[param_Index] == "STATIONARY")
                    {
                        tissue_replication_data[param_Index][1] = 1;
                        tissue_replication_data[param_Index][2] = phase_paramaters[param_Index].first;
                    }
                    else if (phase_Type[param_Index] == "DEPRICIATION")
                    {
                        tissue_replication_data[param_Index][1] = 2;
                        tissue_replication_data[param_Index][2] = phase_paramaters[param_Index].first;
                        tissue_replication_data[param_Index][3] = phase_paramaters[param_Index].second;
                    }
                }
            }

            // for (int profile = 0; profile < number_of_node_Profiles; profile++)
            // {
            //     cout << "Profile :" << profile + 1 << endl;

            //     vector<int> phases_per_Tissue = replication_phases_Profile_tissues[profile];

            //     int tissue = 0;
            //     int phase_Count = 0;

            //     cout << "Tissue: " << tissue + 1 << endl;

            //     for (size_t i = tissue_param_profile_Stride[profile]; i < tissue_param_profile_Stride[profile + 1]; i++)
            //     {
            //         cout << "Time: " << time_Ratios[i] << endl;
            //         phase_Count++;

            //         if (phase_Count == phases_per_Tissue[tissue])
            //         {
            //             cout << "********\n";
            //             tissue++;
            //             cout << "Tissue: " << tissue + 1 << endl;
            //             phase_Count = 0;
            //         }
            //     }
            //     break;
            // }
        }
        else
        {
            cout << "ERROR: NUMBER OF PROFILES MUST BE GREATER THAN ZERO\n\n";
            exit(-1);
        }
    }
    else
    {
        cout << "\nERROR: GENERATION SHOULE BE A NON ZERO POSITIVE VALUE. PLEASE CHECK THE REPLICATION PARAMETERS\n";
        exit(-1);
    }
    cout << endl;
}

void simulator_Master::network_Manager(functions_library &functions)
{
    for (int initialize = 0; initialize < Total_number_of_Nodes; initialize++)
    {
        // vector<pair<int, int>> intialize_Vector;
        // each_Nodes_Connection.push_back(intialize_Vector);

        vector<int> nodes_Intialize_INT;
        each_Nodes_Connection_INT.push_back(nodes_Intialize_INT);
    }

    if (network_Model == "BA")
    {
        BA_Model_Engine(functions);

        // TEST node connections = DONE
        // for (int test = 0; test < each_Nodes_Connection.size(); test++)
        // {
        //     for (size_t i = 0; i < each_Nodes_Connection[test].size(); i++)
        //     {
        //         cout << each_Nodes_Connection[test][i].second << " ";
        //     }
        //     cout << endl;
        // }

        cout << "Completed Barbasi Albert model network engine\n\n";
    }
    else if (network_Model == "SCM")
    {
        SCM_Model_Engine(functions);
    }
    else if (network_Model == "DCM")
    {
        DCM_Model_Engine(functions);
    }
    else if (network_Model == "RANDOM")
    {
        RANDOM_Model_Engine(functions);
    }
    else if (network_Model == "ER_RANDOM")
    {
        ER_RANDOM_Model_Engine(functions);
    }
    else
    {
        cout << "ERROR Incorrect network selected. Please check \"Network type\" in the network parameter file.\n";
        exit(-1);
    }

    // cout << Total_number_of_Nodes << endl
    //      << all_node_IDs.size() << endl;
    // exit(-1);
}

void simulator_Master::ER_RANDOM_Model_Engine(functions_library &functions)
{
    cout << "Intializing Erdos and Renyi Random model network engine\n";

    cout << "Total nodes: " << Total_number_of_Nodes << endl;
    cout << "Probability of paired linakge: " << ER_link_probability << endl
         << endl;

    for (int node = 0; node < Total_number_of_Nodes; node++)
    {
        all_node_IDs.push_back(make_pair(0, node));
    }

    // int number_of_Rounds = (Total_number_of_Nodes * (Total_number_of_Nodes - 1)) / 2;

    random_device rd;
    mt19937 gen(rd());
    // uniform_int_distribution<int> distribution_Neighbour(0, Total_number_of_Nodes - 1);
    bernoulli_distribution link_or_NOT(ER_link_probability);

    fstream network_File;
    network_File.open(network_File_location, ios::app);

    cout << "Forming node to node relationships: \n\n ";

    for (int node_1 = 0; node_1 < Total_number_of_Nodes; node_1++)
    {
        for (int node_2 = node_1 + 1; node_2 < Total_number_of_Nodes; node_2++)
        {
            if (link_or_NOT(gen) == 1)
            {
                cout << "Forming linkage between nodes: ";

                cout << all_node_IDs[node_1].first << "_" << all_node_IDs[node_1].second
                     << " and " << all_node_IDs[node_2].first << "_" << all_node_IDs[node_2].second
                     << endl;

                each_Nodes_Connection_INT[node_1].push_back(node_2);
                each_Nodes_Connection_INT[node_2].push_back(node_1);

                network_File << "0"
                             << "_" << node_1 << "\t"
                             << "0"
                             << "_" << node_2 << "\n";
            }
        }
    }

    network_File.close();
    cout << endl;

    // exit(-1);
}

void simulator_Master::RANDOM_Model_Engine(functions_library &functions)
{
    cout << "Intializing Random model network engine\n";

    random_device rd;
    mt19937 gen(rd());

    fstream network_File;
    network_File.open(network_File_location, ios::app);

    cout << "Forming node to node relationships\n";

    cout << "Configuring node 1 of " << Total_number_of_Nodes << " node(s)" << endl;

    all_node_IDs.push_back(make_pair(0, 0));

    for (int node = 1; node < Total_number_of_Nodes; node++)
    {
        all_node_IDs.push_back(make_pair(0, node));
        cout << "Configuring node " << (node + 1) << " of " << Total_number_of_Nodes << " node(s)" << endl;

        uniform_int_distribution<int> distribution_Neighbour(0, node - 1);

        int attach_Node = -1;
        do
        {
            attach_Node = distribution_Neighbour(gen);
        } while (attach_Node == node);

        // each_Nodes_Connection[node].push_back(make_pair(0, attach_Node));
        // each_Nodes_Connection[attach_Node].push_back(make_pair(0, node));

        each_Nodes_Connection_INT[node].push_back(attach_Node);
        each_Nodes_Connection_INT[attach_Node].push_back(node);

        network_File << "0"
                     << "_" << node << "\t"
                     << "0"
                     << "_" << attach_Node << "\n";
    }

    network_File.close();
    cout << endl;
}

void simulator_Master::DCM_Model_Engine(functions_library &functions)
{
    cout << "Intializing Dynamic Caveman model network engine\n";

    // exit(-1);

    random_device rd;
    mt19937 gen(rd());

    cout << "Determining per cave node counts\n";

    per_cave_Stride = (int *)malloc((DCM_number_of_caves + 1) * sizeof(int));
    per_cave_Stride[0] = 0;

    int **network_Array = functions.create_INT_2D_arrays(DCM_number_of_caves, 3);
    vector<vector<int>> neighbour_Nodes;
    vector<vector<int>> global_Nodes;

    for (int cave_ID = 0; cave_ID < DCM_number_of_caves; cave_ID++)
    {
        int number_of_Nodes;

        if (connection_Model == "NEGATIVE BINOMIAL")
        {
            negative_binomial_distribution<int> distribution(DC_ND_sucesses, DC_ND_probability);
            number_of_Nodes = distribution(gen);
        }
        else if (connection_Model == "POISSON")
        {
            poisson_distribution<int> distribution(BA_Poisson_mean);
            number_of_Nodes = distribution(gen);
        }

        number_of_Nodes = (number_of_Nodes < 1) ? 1 : number_of_Nodes;

        // if (number_of_Nodes < 1)
        // {
        //     number_of_Nodes = 1;
        // }

        network_Array[cave_ID][0] = number_of_Nodes;
        per_cave_Stride[cave_ID + 1] = per_cave_Stride[cave_ID] + network_Array[cave_ID][0];

        vector<int> intialize;
        neighbour_Nodes.push_back(intialize);
        global_Nodes.push_back(intialize);
    }

    Total_number_of_Nodes = per_cave_Stride[DCM_number_of_caves];

    cout << "Total nodes in network: " << Total_number_of_Nodes << endl;

    cout << "Determining relationships for neighour and global nodes\n";

    for (int cave_ID = 0; cave_ID < DCM_number_of_caves; cave_ID++)
    {
        binomial_distribution<int> binomialDist_Neighbour(network_Array[cave_ID][0], DC_percent_Neighbouring);
        network_Array[cave_ID][1] = binomialDist_Neighbour(gen);

        uniform_int_distribution<int> distribution_Neighbour(0, network_Array[cave_ID][0] - 1);

        if (network_Array[cave_ID][1] < 1)
        {
            network_Array[cave_ID][1] = 1;
        }

        vector<int> no_repetition;

        while (neighbour_Nodes[cave_ID].size() < network_Array[cave_ID][1])
        {
            int check_indicator = 0;
            int get_Node = distribution_Neighbour(gen);
            for (int check = 0; check < no_repetition.size(); check++)
            {
                if (get_Node == no_repetition[check])
                {
                    check_indicator = 1;
                    break;
                }
            }
            if (check_indicator == 0)
            {
                neighbour_Nodes[cave_ID].push_back(get_Node);
                no_repetition.push_back(get_Node);
            }
        }

        binomial_distribution<int> binomialDist_Global(network_Array[cave_ID][1], DC_percent_Global_freedom);
        network_Array[cave_ID][2] = binomialDist_Global(gen);

        uniform_int_distribution<int> distribution_Global(0, network_Array[cave_ID][1] - 1);

        no_repetition.clear();

        while (global_Nodes[cave_ID].size() < network_Array[cave_ID][2])
        {
            int check_indicator = 0;
            int get_Index = distribution_Global(gen);

            for (int check = 0; check < no_repetition.size(); check++)
            {
                if (get_Index == no_repetition[check])
                {
                    check_indicator = 1;
                    break;
                }
            }

            if (check_indicator == 0)
            {
                global_Nodes[cave_ID].push_back(neighbour_Nodes[cave_ID][get_Index]);
                no_repetition.push_back(get_Index);
            }
        }
    }

    // vector<vector<pair<int, int>>> each_Nodes_Connections;

    for (size_t i = 0; i < Total_number_of_Nodes; i++)
    {
        // vector<pair<int, int>> nodes_Intialize;
        // each_Nodes_Connection.push_back(nodes_Intialize);

        vector<int> nodes_Intialize_INT;
        each_Nodes_Connection_INT.push_back(nodes_Intialize_INT);
    }

    cout << "Configuring cave cohort relationships\n";

    fstream network_File;
    network_File.open(network_File_location, ios::app);

    for (int cave_ID = 0; cave_ID < DCM_number_of_caves; cave_ID++)
    {
        int nodes_Total = network_Array[cave_ID][0];

        // cout << cave_ID << endl;

        for (int node_parent = 0; node_parent < nodes_Total; node_parent++)
        {
            all_node_IDs.push_back(make_pair(cave_ID, node_parent));

            int node_Count = per_cave_Stride[cave_ID] + node_parent;
            // cout << node_Count << endl;

            // node_Summary << cave_ID << "_" << node_parent << "\t" << cave_ID << "\n";
            for (int node_child = node_parent + 1; node_child < nodes_Total; node_child++)
            {
                int node_child_Main = node_Count + (node_child - node_parent);
                network_File << cave_ID << "_" << node_parent << "\t" << cave_ID << "_" << node_child << "\n";
                // each_Nodes_Connection[node_Count].push_back(make_pair(cave_ID, node_child));
                // each_Nodes_Connection[node_child_Main].push_back(make_pair(cave_ID, node_parent));
                // cout << node_Count << "\t" << node_child_Main << endl;

                each_Nodes_Connection_INT[node_Count].push_back(node_child_Main);
                each_Nodes_Connection_INT[node_child_Main].push_back(node_Count);
            }
        }
        // cout << endl;
    }

    cout << "Configuring neigbouring relationships\n";

    for (int cave_ID = 0; cave_ID < DCM_number_of_caves; cave_ID++)
    {
        uniform_int_distribution<int> distribution_Neighbour_parent(0, neighbour_Nodes[cave_ID].size() - 1);

        int parent_Node_Index = neighbour_Nodes[cave_ID][distribution_Neighbour_parent(gen)];

        int attach_Cave = cave_ID + 1;
        int attach_Node_Index;
        if (attach_Cave != DCM_number_of_caves)
        {
            uniform_int_distribution<int> distribution_Neighbour_attach(0, neighbour_Nodes[attach_Cave].size() - 1);
            attach_Node_Index = neighbour_Nodes[attach_Cave][distribution_Neighbour_attach(gen)];
        }
        else
        {
            attach_Cave = 0;
            uniform_int_distribution<int> distribution_Neighbour_attach(0, neighbour_Nodes[attach_Cave].size() - 1);
            attach_Node_Index = neighbour_Nodes[attach_Cave][distribution_Neighbour_attach(gen)];
        }

        network_File << cave_ID << "_" << parent_Node_Index << "\t" << attach_Cave << "_" << attach_Node_Index << "\n";

        int parent_Node_all_Index = per_cave_Stride[cave_ID] + parent_Node_Index;
        int attach_Node_all_Index = per_cave_Stride[attach_Cave] + attach_Node_Index;

        // each_Nodes_Connection[parent_Node_all_Index].push_back(make_pair(attach_Cave, attach_Node_Index));
        // each_Nodes_Connection[attach_Node_all_Index].push_back(make_pair(cave_ID, parent_Node_Index));

        each_Nodes_Connection_INT[parent_Node_all_Index].push_back(attach_Node_all_Index);
        each_Nodes_Connection_INT[attach_Node_all_Index].push_back(parent_Node_all_Index);
    }

    cout << "Configuring Global networks attachments\n";

    vector<string> repeat_Check;

    for (int cave_ID = 0; cave_ID < DCM_number_of_caves; cave_ID++)
    {
        // cout << cave_ID;
        // if (cave_ID + 1 != DCM_number_of_caves)
        // {
        //     cout << ", ";
        // }
        for (size_t i = 0; i < global_Nodes[cave_ID].size(); i++)
        {
            int cave_Global;
            do
            {
                uniform_int_distribution<int> distribution_global_CAVEs(0, DCM_number_of_caves - 1);
                cave_Global = distribution_global_CAVEs(gen);
                // cout << cave_ID << ":\t" << cave_Global << endl;
            } while (cave_Global == cave_ID);

            if (global_Nodes[cave_Global].size() > 0)
            {
                // cout << "Nodes global: " << global_Nodes[cave_Global].size() << endl;
                uniform_int_distribution<int> distribution_Parent_attach(0, global_Nodes[cave_ID].size() - 1);
                uniform_int_distribution<int> distribution_Global_attach(0, global_Nodes[cave_Global].size() - 1);

                int parent_Cave_Node = global_Nodes[cave_ID][distribution_Parent_attach(gen)];
                int Global_Cave_Node = global_Nodes[cave_Global][distribution_Global_attach(gen)];

                string current = to_string(cave_ID) + "_" + to_string(parent_Cave_Node) + "\t" + to_string(cave_Global) + "_" + to_string(Global_Cave_Node);
                string repeat = to_string(cave_Global) + "_" + to_string(Global_Cave_Node) + "\t" + to_string(cave_ID) + "_" + to_string(parent_Cave_Node);

                int repeat_Index = -1;

                for (int check = 0; check < repeat_Check.size(); check++)
                {
                    if (repeat_Check[check] == repeat || repeat_Check[check] == current)
                    {
                        repeat_Index = 0;
                        break;
                    }
                }

                if (repeat_Index == -1)
                {
                    repeat_Check.push_back(current);
                    repeat_Check.push_back(repeat);
                    network_File << current << "\n";

                    int parent_Node_all_Index = per_cave_Stride[cave_ID] + parent_Cave_Node;
                    int attach_Node_all_Index = per_cave_Stride[cave_Global] + Global_Cave_Node;

                    // each_Nodes_Connection[parent_Node_all_Index].push_back(make_pair(cave_Global, Global_Cave_Node));
                    // each_Nodes_Connection[attach_Node_all_Index].push_back(make_pair(cave_ID, parent_Cave_Node));

                    each_Nodes_Connection_INT[parent_Node_all_Index].push_back(attach_Node_all_Index);
                    each_Nodes_Connection_INT[attach_Node_all_Index].push_back(parent_Node_all_Index);
                }
            }
        }
    }

    functions.clear_Array_int_CPU(network_Array, DCM_number_of_caves);
    network_File.close();

    cout << endl;
}

void simulator_Master::SCM_Model_Engine(functions_library &functions)
{
    cout << "Intializing Standard Caveman model network engine\n";

    random_device rd;
    mt19937 gen(rd());

    cout << "Configuring cave cohort relationships\n";

    vector<int> neighbour_Node;

    uniform_int_distribution<int> distribution(0, SCM_number_of_nodes_per_cave - 1);

    cout << "Configuring cave neighbour nodes" << endl;
    for (int cave_ID = 0; cave_ID < SCM_number_of_caves; cave_ID++)
    {
        neighbour_Node.push_back(distribution(gen));
        // cout << cave_ID << ": " << neighbour_Node[cave_ID] << endl;
    }

    cout << "Configuring cave cohort nodes" << endl;
    for (int cave_ID = 0; cave_ID < SCM_number_of_caves; cave_ID++)
    {
        int node_start_Index = cave_ID * SCM_number_of_nodes_per_cave;

        for (int node_Index_Parent = 0; node_Index_Parent < SCM_number_of_nodes_per_cave; node_Index_Parent++)
        {
            all_node_IDs.push_back(make_pair(cave_ID, node_Index_Parent));

            for (int node_Index_child = node_Index_Parent + 1; node_Index_child < SCM_number_of_nodes_per_cave; node_Index_child++)
            {
                // each_Nodes_Connection[node_start_Index + node_Index_Parent].push_back(make_pair(cave_ID, node_Index_child));
                // each_Nodes_Connection[node_start_Index + node_Index_child].push_back(make_pair(cave_ID, node_Index_Parent));

                each_Nodes_Connection_INT[node_start_Index + node_Index_Parent].push_back(node_start_Index + node_Index_child);
                each_Nodes_Connection_INT[node_start_Index + node_Index_child].push_back(node_start_Index + node_Index_Parent);
            }
        }
    }

    cout << "Re-routing cave cohort nodes to neighbours" << endl;

    vector<int> index_of_neighbour_CAVES;

    for (int cave_ID = 0; cave_ID < SCM_number_of_caves; cave_ID++)
    {
        int node_start_Index = cave_ID * SCM_number_of_nodes_per_cave;
        int rerouting_Node = neighbour_Node[cave_ID] + node_start_Index;

        uniform_int_distribution<int> rerouting_Connection_Draw(0, each_Nodes_Connection_INT[rerouting_Node].size() - 1);

        int rerouting_Connection = rerouting_Connection_Draw(gen);
        int rerouting_Connection_Node = each_Nodes_Connection_INT[rerouting_Node][rerouting_Connection];

        int index = -1;
        for (int find = 0; find < each_Nodes_Connection_INT[rerouting_Connection_Node].size(); find++)
        {
            if (each_Nodes_Connection_INT[rerouting_Connection_Node][find] == rerouting_Node)
            {
                index = find;
                break;
            }
        }

        if (index == -1)
        {
            cout << "ERROR in node rerouting removal\n";
        }
        else
        {
            each_Nodes_Connection_INT[rerouting_Node].erase(each_Nodes_Connection_INT[rerouting_Node].begin() + rerouting_Connection);
            each_Nodes_Connection_INT[rerouting_Connection_Node].erase(each_Nodes_Connection_INT[rerouting_Connection_Node].begin() + index);
        }
    }

    for (int cave_ID = 0; cave_ID < SCM_number_of_caves; cave_ID++)
    {
        int node_start_Index = cave_ID * SCM_number_of_nodes_per_cave;
        int rerouting_Node = neighbour_Node[cave_ID] + node_start_Index;

        if ((cave_ID + 1) != SCM_number_of_caves)
        {
            int node_start_Index_neighboour = (cave_ID + 1) * SCM_number_of_nodes_per_cave;
            int rerouting_Node_neighboour = neighbour_Node[cave_ID + 1] + node_start_Index_neighboour;
            // each_Nodes_Connection[rerouting_Node][rerouting_Connection].first = cave_ID + 1;
            // each_Nodes_Connection[rerouting_Node][rerouting_Connection].second = neighbour_Node[cave_ID + 1];

            each_Nodes_Connection_INT[rerouting_Node].push_back(rerouting_Node_neighboour);
            each_Nodes_Connection_INT[rerouting_Node_neighboour].push_back(rerouting_Node);
        }
        else
        {
            int rerouting_Node_neighboour = neighbour_Node[0];

            each_Nodes_Connection_INT[rerouting_Node].push_back(rerouting_Node_neighboour);
            each_Nodes_Connection_INT[rerouting_Node_neighboour].push_back(rerouting_Node);

            // each_Nodes_Connection[rerouting_Node][rerouting_Connection].first = 0;
            // each_Nodes_Connection[rerouting_Node][rerouting_Connection].second = neighbour_Node[0];
        }
    }

    // for (int cave_ID = 0; cave_ID < SCM_number_of_caves; cave_ID++)
    // {
    //     int node_start_Index = cave_ID * SCM_number_of_nodes_per_cave;
    //     int rerouting_Node = neighbour_Node[cave_ID] + node_start_Index;

    //     uniform_int_distribution<int> rerouting_Connection_Draw(0, each_Nodes_Connection_INT[rerouting_Node].size() - 1);

    //     int rerouting_Connection = rerouting_Connection_Draw(gen);
    //     int rerouting_Connection_Node = each_Nodes_Connection_INT[rerouting_Node][rerouting_Connection];

    //     int index = -1;
    //     for (int find = 0; find < each_Nodes_Connection[rerouting_Connection_Node].size(); find++)
    //     {
    //         if (each_Nodes_Connection[rerouting_Connection_Node][find] == rerouting_Node)
    //         {
    //             index = find;
    //             break;
    //         }
    //     }

    //     if (index == -1)
    //     {
    //         cout << "ERROR in node rerouting removal\n";
    //     }
    //     else
    //     {
    //         each_Nodes_Connection_INT[rerouting_Connection_Node].erase(each_Nodes_Connection[rerouting_Connection_Node].begin() + index);
    //     }

    //     if ((cave_ID + 1) != SCM_number_of_caves)
    //     {
    //         // each_Nodes_Connection[rerouting_Node][rerouting_Connection].first = cave_ID + 1;
    //         // each_Nodes_Connection[rerouting_Node][rerouting_Connection].second = neighbour_Node[cave_ID + 1];

    //         each_Nodes_Connection_INT[rerouting_Node][rerouting_Connection] = neighbour_Node[cave_ID + 1] + ((cave_ID + 1) * SCM_number_of_nodes_per_cave);
    //     }
    //     else
    //     {
    //         // each_Nodes_Connection[rerouting_Node][rerouting_Connection].first = 0;
    //         // each_Nodes_Connection[rerouting_Node][rerouting_Connection].second = neighbour_Node[0];
    //     }
    // }

    fstream network_File;
    network_File.open(network_File_location, ios::app);

    for (int cave_ID = 0; cave_ID < SCM_number_of_caves; cave_ID++)
    {
        cout << "Finalizing cave " << cave_ID + 1 << " of " << SCM_number_of_caves << endl;
        int node_start_Index = cave_ID * SCM_number_of_nodes_per_cave;

        for (int node = 0; node < SCM_number_of_nodes_per_cave; node++)
        {
            for (int connections = 0; connections < each_Nodes_Connection_INT[node_start_Index + node].size(); connections++)
            {
                if (all_node_IDs[each_Nodes_Connection_INT[node_start_Index + node][connections]].first == cave_ID)
                {
                    if (all_node_IDs[each_Nodes_Connection_INT[node_start_Index + node][connections]].second > node)
                    {
                        network_File << all_node_IDs[node_start_Index + node].first << "_" << all_node_IDs[node_start_Index + node].second
                                     << "\t" << all_node_IDs[each_Nodes_Connection_INT[node_start_Index + node][connections]].first << "_" << all_node_IDs[each_Nodes_Connection_INT[node_start_Index + node][connections]].second << endl;

                        // network_File << cave_ID << "_" << to_string(node) << "\t"
                        //              << cave_ID << "_" << to_string(each_Nodes_Connection[node_start_Index + node][connections].second) << "\n";
                    }
                }
                else if (all_node_IDs[each_Nodes_Connection_INT[node_start_Index + node][connections]].first > cave_ID)
                {
                    network_File << all_node_IDs[node_start_Index + node].first << "_" << all_node_IDs[node_start_Index + node].second
                                 << "\t" << all_node_IDs[each_Nodes_Connection_INT[node_start_Index + node][connections]].first << "_" << all_node_IDs[each_Nodes_Connection_INT[node_start_Index + node][connections]].second << endl;
                }
                // if (each_Nodes_Connection[node_start_Index + node][connections].first == cave_ID)
                // {
                //     if (each_Nodes_Connection[node_start_Index + node][connections].second > node)
                //     {
                //         network_File << cave_ID << "_" << to_string(node) << "\t"
                //                      << cave_ID << "_" << to_string(each_Nodes_Connection[node_start_Index + node][connections].second) << "\n";
                //     }
                // }
                // else
                // {
                //     network_File << to_string(cave_ID) << "_" << to_string(node) << "\t"
                //                  << to_string(each_Nodes_Connection[node_start_Index + node][connections].first) << "_" << to_string(each_Nodes_Connection[node_start_Index + node][connections].second) << "\n";
                // }
            }
        }
    }

    network_File.close();

    // for (int cave_ID = 0; cave_ID < SCM_number_of_caves; cave_ID++)
    // {
    //     int node_start_Index = cave_ID * SCM_number_of_nodes_per_cave;
    //     for (size_t i = node_start_Index; i < (node_start_Index + SCM_number_of_nodes_per_cave); i++)
    //     {
    //         for (int x = 0; x < each_Nodes_Connection[i].size(); x++)
    //         {
    //             cout << each_Nodes_Connection[i][x].first << "_" << each_Nodes_Connection[i][x].second << ", ";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }

    // exit(-1);

    cout << endl;
}

void simulator_Master::BA_Model_Engine(functions_library &functions)
{
    cout << "Intializing Barbasi Albert model network engine\n";

    // string output_Network_location = this->output_Folder_location + "/network_Data";
    // functions.config_Folder(output_Network_location, "Network");
    // network_File_location = output_Network_location + "/node_node_Relationships.csv";
    // functions.create_File(network_File_location, "Source\tTarget");

    // for (int initialize = 0; initialize < Total_number_of_Nodes; initialize++)
    // {
    //     vector<pair<int, int>> intialize_Vector;
    //     each_Nodes_Connection.push_back(intialize_Vector);
    // }

    cout << "Forming node to node relationships\n";

    vector<int> connections_per_Node;
    connections_per_Node.push_back(1);
    all_node_IDs.push_back(make_pair(0, 0));
    cout << "Configuring node 1 of " << Total_number_of_Nodes << " node(s)" << endl;

    int tot_connections = 2;

    // int iterative_Attach = BA_FIXED;

    fstream network_File;
    network_File.open(network_File_location, ios::app);

    random_device rd;
    mt19937 gen(rd());

    for (int node = 1; node < Total_number_of_Nodes; node++)
    {
        cout << "Configuring node " << (node + 1) << " of " << Total_number_of_Nodes << " node(s)" << endl;
        connections_per_Node.push_back(0);

        all_node_IDs.push_back(make_pair(0, node));

        int iterative_Attach;
        if (connection_Model == "FIXED")
        {
            iterative_Attach = BA_FIXED;
        }
        else if (connection_Model == "NEGATIVE BINOMIAL")
        {
            negative_binomial_distribution<int> distribution(BA_NB_sucesses, BA_NB_probability);
            iterative_Attach = distribution(gen);
        }
        else if (connection_Model == "POISSON")
        {
            poisson_distribution<int> distribution(BA_Poisson_mean);
            iterative_Attach = distribution(gen);
        }
        else
        {
            cout << "ERROR in network parameter \"BA model fixed new connections\". CHECK PARAMETER\n";
            exit(-1);
        }

        iterative_Attach = (iterative_Attach < 1) ? 1 : iterative_Attach;

        for (int interative = 0; interative < iterative_Attach; interative++)
        {
            tot_connections++;
            int attach_Node = node;

            while (attach_Node == node)
            {
                float randomNum = static_cast<float>(std::rand()) / RAND_MAX;
                float cum_Prob = 0;

                for (int check_Node = 0; check_Node < connections_per_Node.size(); check_Node++)
                {
                    cum_Prob += ((float)(connections_per_Node[check_Node] + 1) / (float)tot_connections);
                    if (randomNum < cum_Prob)
                    {
                        attach_Node = check_Node;
                        break;
                    }
                }
            }

            if (attach_Node != -1)
            {
                int check_Present = 0;

                for (int check = 0; check < each_Nodes_Connection_INT[node].size(); check++)
                {
                    if (attach_Node == each_Nodes_Connection_INT[node][check])
                    {
                        check_Present = 1;
                        break;
                    }
                }

                if (check_Present == 0)
                {
                    // cout << "Node " << node + 1 << " attached to " << attach_Node + 1 << endl;
                    network_File << "0_" << to_string(attach_Node) << "\t"
                                 << "0_" << to_string(node) << "\n";

                    connections_per_Node[attach_Node] = connections_per_Node[attach_Node] + 1;

                    // each_Nodes_Connection[attach_Node].push_back(make_pair(0, node));
                    // each_Nodes_Connection[node].push_back(make_pair(0, attach_Node));

                    each_Nodes_Connection_INT[attach_Node].push_back(node);
                    each_Nodes_Connection_INT[node].push_back(attach_Node);

                    tot_connections++;
                }
            }
            else
            {
                cout << "ERROR in node configuration\n";
                exit(-1);
            }
        }
    }

    network_File.close();
}
