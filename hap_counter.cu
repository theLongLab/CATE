#include "hap_counter.cuh"

hap_counter::hap_counter(string parameter_Master_Location)
{
    cout << "Intializing Haplotype retriever\n";

    //exit(-1);

    parameter_load Parameters = parameter_load();
    functions_library function = functions_library();

    vector<string> parameters_List = {
        "\"Intermediate folders\"",
        "\"Output folders\"",
        "\"CPU cores\"",
        "\"Multi read\"",
        "\"Haplotype Node IDs\"",
        "\"Nodes master profile\""};

    vector<string> found_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List);

    node_Master_location = Parameters.get_STRING(found_Parameters[5]);

    cout << "\nNode master profile file: " << node_Master_location << endl;

    cout << "\nConfiguring folders:\n";

    this->intermediate_Folder_location = Parameters.get_STRING(found_Parameters[0]);
    this->output_Folder_location = Parameters.get_STRING(found_Parameters[1]);

    if (filesystem::exists(intermediate_Folder_location) && filesystem::is_directory(intermediate_Folder_location))
    {
        cout << "Intermediary folder configured: " << intermediate_Folder_location << "\n";
        intermediary_Sequence_location = intermediate_Folder_location + "/sequence_Data";
        intermediary_Index_location = intermediate_Folder_location + "/index_Data";

        if (filesystem::exists(intermediary_Sequence_location) && filesystem::is_directory(intermediary_Sequence_location))
        {
            cout << "Intermediary sequence data folder configured: " << intermediary_Sequence_location << "\n";
            if (filesystem::exists(intermediary_Index_location) && filesystem::is_directory(intermediary_Index_location))
            {
                cout << "Intermediary node index folder configured: " << intermediary_Index_location << "\n";
            }
            else
            {
                cout << "INTERMEDIARY NODE INDEX FOLDER NOT FOUND AT: " << this->intermediary_Index_location << endl;
            }
        }
        else
        {
            cout << "INTERMEDIARY SEQUENCE DATA FOLDER NOT FOUND AT: " << this->intermediary_Sequence_location << endl;
        }

        if (filesystem::exists(output_Folder_location) && filesystem::is_directory(output_Folder_location))
        {
            cout << "Output folder configured: " << output_Folder_location << "\n";
        }
        else
        {
            cout << "ERROR: OUTPUT FOLDER NOT FOUND: " << this->output_Folder_location << endl;
        }
    }
    else
    {
        cout << "INTERMEDIARY FOLDER NOT FOUND AT: " << this->intermediate_Folder_location << endl;
    }

    cout << "\nConfiguring hardware resources:\n\n";

    this->CPU_cores = Parameters.get_INT(found_Parameters[2]);
    this->multi_Read = Parameters.get_STRING(found_Parameters[3]);

    cout << "Available CPU cores: " << this->CPU_cores << endl
         << endl;

    transform(multi_Read.begin(), multi_Read.end(), multi_Read.begin(), ::toupper);
    cout << "Multiple read and write: " << this->multi_Read << endl
         << endl;

    nodes_to_Analyse = Parameters.get_STRING(found_Parameters[4]);
    function.to_Upper_Case(nodes_to_Analyse);
    cout << "Node(s) to analyse: " << nodes_to_Analyse << endl
         << endl;
}

void hap_counter::ingress()
{
    functions_library functions = functions_library();
    parameter_load Parameters = parameter_load();

    cout << "Loading nodes master profile: " << this->node_Master_location << endl;

    vector<string> tissue_Names;
    vector<pair<string, string>> Tissue_profiles_block_Data = Parameters.get_block_from_File(node_Master_location, "Tissue profiles");
    int num_tissues_per_Node = Parameters.get_INT(Tissue_profiles_block_Data, "Number of tissues");

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

    cout << "\nReading node index data: " << intermediary_Index_location + "/node_Index.csv\n";

    fstream index_Data_file;
    index_Data_file.open(intermediary_Index_location + "/node_Index.csv", ios::in);

    vector<pair<int, string>> index_Information;

    if (index_Data_file.is_open())
    {
        string line;
        getline(index_Data_file, line);

        vector<string> split_Data;

        while (getline(index_Data_file, line))
        {
            functions.split(split_Data, line, '\t');
            index_Information.push_back(make_pair(stoi(split_Data[0]), split_Data[1]));
        }
        index_Data_file.close();
    }
    else
    {
        cout << "ERROR: UNABLE TO OPEN NODE INDEX FILE: " << intermediary_Index_location + "/node_Index.csv\n";
    }

    cout << "Node index data read, " << index_Information.size() << " node(s) indexed\n";

    if (nodes_to_Analyse != "ALL")
    {
        vector<string> nodes;
        functions.split(nodes, nodes_to_Analyse, ',');

        for (int node = 0; node < nodes.size(); node++)
        {
            int index = -1;
            cout << "\nProcessing node: " << nodes[node] << endl;

            for (int find = 0; find < index_Information.size(); find++)
            {
                if (index_Information[find].second == nodes[node])
                {
                    index = find;
                    break;
                }
            }
            if (index != -1)
            {
                cout << "Node found, index : " << index << endl;

                if (filesystem::exists(output_Folder_location + "/node_Data/" + nodes[node]) && filesystem::is_directory(output_Folder_location + "/node_Data/" + nodes[node]))
                {
                    cout << "Node was infected: " << output_Folder_location << "/node_Data/" << nodes[node];

                    if (filesystem::exists(intermediary_Sequence_location + "/" + to_string(index) + ".tar") && !filesystem::exists(intermediary_Sequence_location + "/" + to_string(index)))
                    {
                        cout << "Extracting tar directory: " << intermediary_Sequence_location << "/" << to_string(index) << ".tar\n";
                        string command = "tar -xf" + intermediary_Sequence_location + "/" + to_string(index) + ".tar -C .";

                        int result = system(command.c_str());

                        if (result == 0)
                        {
                            // The command executed successfully
                            cout << "Successfully untarred the folder." << endl;
                        }
                        else
                        {
                            // An error occurred during the execution of the command
                            cout << "Failed to untar the folder." << endl;
                            exit(-1);
                        }
                    }

                    if (filesystem::exists(intermediary_Sequence_location + "/" + to_string(index)))
                    {
                        string intermediary_Folder = intermediary_Sequence_location + "/" + to_string(index);
                        cout << "\nIntermediary directory found: " << intermediary_Folder << endl;

                        cout << "\nLoading Parent sequence data: ";

                        fstream parent_Sequence_file;
                        parent_Sequence_file.open(output_Folder_location + "/node_Data/" + nodes[node] + "/cells_of_Parents.csv", ios::in);

                        // 3D vector
                        vector<vector<vector<int>>> tissue_generation_Sequence;

                        for (int tissue = 0; tissue < num_tissues_per_Node; tissue++)
                        {
                            vector<vector<int>> tissue_Vector;
                            tissue_generation_Sequence.push_back(tissue_Vector);
                        }

                        if (parent_Sequence_file.is_open())
                        {
                            cout << output_Folder_location << "/node_Data/" << nodes[node] << "/cells_of_Parents.csv\n";
                            string line;
                            vector<string> line_Data;

                            getline(parent_Sequence_file, line);

                            while (getline(parent_Sequence_file, line))
                            {
                                functions.split(line_Data, line, '\t');
                                vector<string> sequence_Data;
                                functions.split(sequence_Data, line_Data[0], '_');

                                int tissue_Index = -1;
                                for (int tissue = 0; tissue < num_tissues_per_Node; tissue++)
                                {
                                    if (tissue_Names[tissue] == sequence_Data[2])
                                    {
                                        tissue_Index = tissue;
                                        break;
                                    }
                                }

                                if (tissue_Index != -1)
                                {
                                    int generation = stoi(sequence_Data[3]);
                                    while (tissue_generation_Sequence[tissue_Index].size() <= generation)
                                    {
                                        vector<int> generation_Vector;
                                        tissue_generation_Sequence[tissue_Index].push_back(generation_Vector);
                                    }
                                    tissue_generation_Sequence[tissue_Index][generation].push_back(stoi(sequence_Data[4]));
                                }
                                else
                                {
                                    cout << "ERROR UNKNOWN TISSUE FOUND: " << sequence_Data[2] << endl;
                                    exit(-1);
                                }
                            }

                            parent_Sequence_file.close();
                        }
                        else
                        {
                            cout << "ERROR: UNABLE TO OPEN PARENT SEQUENCE FILE: " << output_Folder_location << "/node_Data/" << nodes[node] << "/cells_of_Parents.csv\n";
                            exit(-1);
                        }

                        for (int tissue = 0; tissue < num_tissues_per_Node; tissue++)
                        {
                            for (int gen = 0; gen < tissue_generation_Sequence[tissue].size(); gen++)
                            {
                                sort(tissue_generation_Sequence[tissue][gen].begin(), tissue_generation_Sequence[tissue][gen].end());
                            }
                        }

                        // for (int test = 0; test < tissue_generation_Sequence[1].size(); test++)
                        // {
                        //     for (int seq = 0; seq < tissue_generation_Sequence[1][test].size(); seq++)
                        //     {
                        //         cout << tissue_generation_Sequence[1][test][seq] << endl;
                        //     }
                        //     cout << endl;
                        // }

                        // exit(-1);

                        functions.create_File(output_Folder_location + "/node_Data/" + nodes[node] + "/all_Haplotype_Frequencies.csv", "Tissue\tGeneration\tSequence\tCount\tFrequency");
                        functions.create_File(output_Folder_location + "/node_Data/" + nodes[node] + "/alive_Haplotype_Frequencies.csv", "Tissue\tGeneration\tSequence\tCount\tFrequency");
                        functions.create_File(output_Folder_location + "/node_Data/" + nodes[node] + "/parent_Haplotype_Frequencies.csv", "Tissue\tGeneration\tSequence\tCount\tFrequency");

                        functions.create_File(output_Folder_location + "/node_Data/" + nodes[node] + "/all_Haplotype_Count.csv", "Tissue\tGeneration\thaplotype_Count");
                        functions.create_File(output_Folder_location + "/node_Data/" + nodes[node] + "/alive_Haplotype_Count.csv", "Tissue\tGeneration\thaplotype_Count");
                        functions.create_File(output_Folder_location + "/node_Data/" + nodes[node] + "/parent_Haplotype_Count.csv", "Tissue\tGeneration\thaplotype_Count");

                        for (int tissue = 0; tissue < num_tissues_per_Node; tissue++)
                        {
                            cout << "\nConfiguring tissue: " << tissue_Names[tissue] << endl;
                            string tissue_Folder = intermediary_Folder + "/" + to_string(tissue);

                            vector<pair<int, string>> generation_Folder_Path;
                            vector<string> split_Data;

                            for (const auto &entry : filesystem::directory_iterator(tissue_Folder))
                            {
                                if (filesystem::is_directory(entry.path()) || entry.path().extension().string() == ".tar")
                                {
                                    string generation_Directory = entry.path().stem().string();
                                    functions.split(split_Data, generation_Directory, '_');
                                    generation_Folder_Path.push_back(make_pair(stoi(split_Data[1]), entry.path().string()));
                                }
                            }

                            sort(generation_Folder_Path.begin(), generation_Folder_Path.end());
                            cout << generation_Folder_Path.size() << " generation(s) found\n";

                            for (int generation = 0; generation < generation_Folder_Path.size(); generation++)
                            {
                                int track_Seq = 0;

                                cout << "\nProcessing generation: " << generation << endl;
                                if (filesystem::path(generation_Folder_Path[generation].second).extension().string() == ".tar")
                                {
                                    cout << "Extracting tar directory: " << generation_Folder_Path[generation].second << endl;
                                    string command = "tar -xf" + generation_Folder_Path[generation].second + " -C .";

                                    int result = system(command.c_str());

                                    if (result == 0)
                                    {
                                        // The command executed successfully
                                        cout << "Successfully untarred the folder." << endl;
                                    }
                                    else
                                    {
                                        // An error occurred during the execution of the command
                                        cout << "Failed to untar the folder." << endl;
                                        exit(-1);
                                    }
                                }
                                // cout << filesystem::path(generation_Folder_Path[generation].second).stem().string() << endl;
                                if (filesystem::exists(tissue_Folder + "/" + filesystem::path(generation_Folder_Path[generation].second).stem().string()))
                                {
                                    cout << "Reading folder: " << tissue_Folder + "/" + filesystem::path(generation_Folder_Path[generation].second).stem().string() << endl
                                         << endl;

                                    vector<pair<int, int>> nFASTA_files = functions.index_Source_folder(intermediary_Folder, tissue, generation);

                                    // Multithread
                                    // vector<pair<int, string>> all_Hap_Count;
                                    // vector<pair<int, string>> all_Hap_Alive_Count;
                                    // vector<pair<int, string>> all_Hap_Parent_Count;

                                    cout << "Detecting haplotypes and their counts: ";

                                    for (int nFASTA_file = 0; nFASTA_file < nFASTA_files.size(); nFASTA_file++)
                                    {
                                        string file = tissue_Folder + "/" + (filesystem::path(generation_Folder_Path[generation].second).stem().string()) + "/" + to_string(nFASTA_files[nFASTA_file].first) + "_" + to_string(nFASTA_files[nFASTA_file].second) + ".nfasta";

                                        fstream read_nFASTA;
                                        read_nFASTA.open(file, ios::in);

                                        if (read_nFASTA.is_open())
                                        {
                                            vector<pair<string, string>> line_Data;
                                            string line_Name;
                                            while (getline(read_nFASTA, line_Name))
                                            {
                                                string line_Sequence;
                                                getline(read_nFASTA, line_Sequence);
                                                line_Data.push_back(make_pair(line_Name.substr(1, line_Name.length()), line_Sequence));
                                            }
                                            read_nFASTA.close();

                                            // for (int test = 0; test < line_Data.size(); test++)
                                            // {
                                            //     cout << line_Data[test].first << endl;
                                            //     //      << endl
                                            //     //      << line_Data[test].second << endl;
                                            // }

                                            // exit(-1);

                                            vector<thread> threads_vec;

                                            threads_vec.push_back(thread{&hap_counter::all_Haplotype_Counter, this, ref(line_Data)});
                                            threads_vec.push_back(thread{&hap_counter::all_Haplotype_Alive_Counter, this, ref(line_Data), ref(functions)});
                                            if (generation < tissue_generation_Sequence[tissue].size())
                                            {
                                                threads_vec.push_back(thread{&hap_counter::all_Haplotype_Parent_Counter, this, ref(line_Data), ref(functions), ref(tissue_generation_Sequence), tissue, generation, ref(track_Seq)});
                                            }

                                            for (thread &t : threads_vec)
                                            {
                                                if (t.joinable())
                                                {
                                                    t.join();
                                                }
                                            }

                                            threads_vec.clear();

                                            // for (int test = 0; test < all_Hap_Count.size(); test++)
                                            // {
                                            //     // if (all_Hap_Count[test].first != 1)
                                            //     //{
                                            //     cout << all_Hap_Count[test].first << endl;
                                            //     cout << all_Hap_Count[test].second << endl;
                                            //     //}
                                            // }

                                            // exit(-1);
                                        }
                                        else
                                        {
                                            cout << "ERROR: UNABLE TO OPEN NFASTA FILE: " << file << endl;
                                            exit(-1);
                                        }
                                    }
                                    // exit(-1);
                                    cout << "Done\n";

                                    if (filesystem::path(generation_Folder_Path[generation].second).extension().string() == ".tar")
                                    {
                                        string command = "rm -r " + tissue_Folder + "/" + filesystem::path(generation_Folder_Path[generation].second).stem().string();
                                        int result = system(command.c_str());

                                        if (result == 0)
                                        {
                                            // The command executed successfully
                                            cout << "Successfully deleted the folder." << endl;
                                        }
                                        else
                                        {
                                            // An error occurred during the execution of the command
                                            cout << "Failed to delete the folder." << endl;
                                            exit(-1);
                                        }
                                    }
                                }
                                // exit(-1);

                                // Calculate totals and write the frequencies
                                write_Files(tissue_Names[tissue], generation, all_Hap_Count,
                                            output_Folder_location + "/node_Data/" + nodes[node] + "/all_Haplotype_Frequencies.csv",
                                            output_Folder_location + "/node_Data/" + nodes[node] + "/all_Haplotype_Count.csv", all_Hap_Total);

                                write_Files(tissue_Names[tissue], generation, all_Hap_Alive_Count,
                                            output_Folder_location + "/node_Data/" + nodes[node] + "/alive_Haplotype_Frequencies.csv",
                                            output_Folder_location + "/node_Data/" + nodes[node] + "/alive_Haplotype_Count.csv", all_Hap_Alive_Total);

                                write_Files(tissue_Names[tissue], generation, all_Hap_Parent_Count,
                                            output_Folder_location + "/node_Data/" + nodes[node] + "/parent_Haplotype_Frequencies.csv",
                                            output_Folder_location + "/node_Data/" + nodes[node] + "/parent_Haplotype_Count.csv", all_Hap_Parent_Total);

                                all_Hap_Count.clear();
                                all_Hap_Alive_Count.clear();
                                all_Hap_Parent_Count.clear();

                                all_Hap_Total = 0;
                                all_Hap_Alive_Total = 0;
                                all_Hap_Parent_Total = 0;

                                // track_Seq = 0;
                            }
                            // exit(-1);
                        }
                    }
                }
                else
                {
                    cout << "Node was not infected, no folder at: " << output_Folder_location << "/node_Data/" << nodes[node] << "\n";
                }
            }
            else
            {
                cout << "Node not found: " << nodes[node] << endl;
            }
        }
    }
}

void hap_counter::write_Files(string tissue_Name, int generation, vector<pair<int, string>> &Hap_count, string location_Frequencies, string location_Summaries, int total)
{
    fstream all_hap_File;
    all_hap_File.open(location_Frequencies, ios::app);
    if (all_hap_File.is_open())
    {
        for (int haplotype = 0; haplotype < Hap_count.size(); haplotype++)
        {
            all_hap_File << tissue_Name << "\t"
                         << to_string(generation) << "\t"
                         << "\"" << Hap_count[haplotype].second << "\""
                         << "\t"
                         << Hap_count[haplotype].first << "\t";
            float frequency = (float)Hap_count[haplotype].first / (float)total;
            all_hap_File << frequency << "\n";
        }
        all_hap_File.close();
    }

    fstream all_hap_Summary;
    all_hap_Summary.open(location_Summaries, ios::app);
    all_hap_Summary << tissue_Name << "\t"
                    << to_string(generation) << "\t"
                    << to_string(Hap_count.size()) << "\n";
    all_hap_Summary.close();
}

void hap_counter::all_Haplotype_Counter(vector<pair<string, string>> &line_Data)
{
    for (int check = 0; check < line_Data.size(); check++)
    {
        int found = 0;
        for (int hap = 0; hap < all_Hap_Count.size(); hap++)
        {
            if (line_Data[check].second == all_Hap_Count[hap].second)
            {
                all_Hap_Count[hap].first = all_Hap_Count[hap].first + 1;
                found = 1;
                break;
            }
        }
        if (found == 0)
        {
            all_Hap_Count.push_back(make_pair(1, line_Data[check].second));
        }
        all_Hap_Total++;
    }
}

void hap_counter::all_Haplotype_Alive_Counter(vector<pair<string, string>> &line_Data, functions_library &functions)
{
    vector<string> sequence_Name_data;

    for (int check = 0; check < line_Data.size(); check++)
    {
        functions.split(sequence_Name_data, line_Data[check].first, '_');
        if (sequence_Name_data[1] == "A")
        {
            int found = 0;
            for (int hap = 0; hap < all_Hap_Alive_Count.size(); hap++)
            {
                if (line_Data[check].second == all_Hap_Alive_Count[hap].second)
                {
                    all_Hap_Alive_Count[hap].first = all_Hap_Alive_Count[hap].first + 1;
                    found = 1;
                    break;
                }
            }
            if (found == 0)
            {
                all_Hap_Alive_Count.push_back(make_pair(1, line_Data[check].second));
            }
            all_Hap_Alive_Total++;
        }
    }
}

void hap_counter::all_Haplotype_Parent_Counter(vector<pair<string, string>> &line_Data, functions_library &functions,
                                               vector<vector<vector<int>>> &tissue_generation_Sequence,
                                               int tissue, int generation, int &track_Seq)
{
    vector<int> sequences = tissue_generation_Sequence[tissue][generation];

    if (sequences.size() != 0 && track_Seq < sequences.size())
    {
        vector<string> sequence_Name_data;
        for (int check = 0; check < line_Data.size(); check++)
        {
            functions.split(sequence_Name_data, line_Data[check].first, '_');
            if (stoi(sequence_Name_data[0]) == sequences[track_Seq])
            {
                int found = 0;
                for (int hap = 0; hap < all_Hap_Parent_Count.size(); hap++)
                {
                    if (line_Data[check].second == all_Hap_Parent_Count[hap].second)
                    {
                        all_Hap_Parent_Count[hap].first = all_Hap_Parent_Count[hap].first + 1;
                        found = 1;
                        break;
                    }
                }
                if (found == 0)
                {
                    all_Hap_Parent_Count.push_back(make_pair(1, line_Data[check].second));
                }

                track_Seq++;
                all_Hap_Parent_Total++;

                if (track_Seq >= sequences.size())
                {
                    break;
                }
            }
        }
    }
}