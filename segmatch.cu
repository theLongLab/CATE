#include "segmatch.cuh"

segmatch::segmatch(string parameter_Master_Location)
{
    cout << "Intializing seg matching\n";

   // exit(-1);

    parameter_load Parameters = parameter_load();
    functions_library function = functions_library();

    vector<string> parameters_List = {
        "\"Intermediate folders\"",
        "\"Output folders\"",
        "\"Nodes master profile\"",
        "\"Sequence master profile\"",
        "\"CPU cores\"",
        "\"Segregating match sequences\"",
        "\"Segregating match Node ID\"",
        "\"Segregating match tissue\"",
        "\"Segregating match cutoff\"",
        "\"Multi read\""};

    vector<string> found_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List);

    intermediate_Folder_location = Parameters.get_STRING(found_Parameters[0]);
    output_Folder_location = Parameters.get_STRING(found_Parameters[1]);
    CPU_cores = Parameters.get_INT(found_Parameters[4]);

    if (function.to_Upper_Case(found_Parameters[9]) == "NO")
    {
        CPU_cores = 1;
    }

    cout << "\nConfiguring segregating sites: \n";
    string sequence_Master_location = Parameters.get_STRING(found_Parameters[3]);

    vector<pair<string, string>> mutations_Block = Parameters.get_block_from_File(sequence_Master_location, "Mutations");
    int num_mutation_Hotspots = Parameters.get_INT(mutations_Block, "Number of hotspots");

    if (num_mutation_Hotspots > 0)
    {
        cout << "\nProcessing " << num_mutation_Hotspots << " mutation hotspots: \n";
        vector<string> split_Region;
        for (int hotspot = 0; hotspot < num_mutation_Hotspots; hotspot++)
        {
            string hotspot_ID = "Hotspot " + to_string(hotspot + 1);
            cout << "\nProcessing: " << hotspot_ID << endl;
            vector<pair<string, string>> mutations_hotspot_Block = Parameters.get_block_from_block(mutations_Block, hotspot_ID);

            string region = Parameters.get_STRING(mutations_hotspot_Block, "Region");
            function.split(split_Region, region, '_');

            int start_Pos = stoi(split_Region[0]) - 1;
            int stop_Pos = stoi(split_Region[1]) - 1;
            positions_start_end.push_back(make_pair(start_Pos, stop_Pos));

            cout << "Start position: " << (start_Pos + 1) << " | Stop position: " << (stop_Pos + 1) << endl;
            int bases = stop_Pos - start_Pos + 1;
            cout << "Number of bases: " << bases << endl;
            total_Bases = total_Bases + bases;
        }

        cout << "\nTotal number of segregating sites to check: " << total_Bases << endl;

        tissue_Name = Parameters.get_STRING(found_Parameters[7]);
        cout << "\nIdentifying tissue index: " << tissue_Name << endl;

        string node_Master_location = Parameters.get_STRING(found_Parameters[2]);
        vector<pair<string, string>> Tissues_Block = Parameters.get_block_from_File(node_Master_location, "Tissue profiles");
        int num_Tissues = Parameters.get_INT(Tissues_Block, "Number of tissues");

        if (num_Tissues > 0)
        {
            for (int tissue = 0; tissue < num_Tissues; tissue++)
            {
                string check = "Tissue " + to_string(tissue + 1) + " Name";
                string Name = Parameters.get_STRING(Tissues_Block, check);

                if (Name == tissue_Name)
                {
                    tissue_Index = tissue;
                    break;
                }
            }

            if (tissue_Index != -1)
            {
                cout << "Index: " << tissue_Index << endl;

                cout << "\nGetting node index\n";
                string index_Node_File = intermediate_Folder_location + "/index_Data/node_Index.csv";
                node_ID = Parameters.get_STRING(found_Parameters[6]);
                cout << "Node: " << node_ID << endl;

                fstream index_File;
                index_File.open(index_Node_File, ios::in);

                if (index_File.is_open())
                {
                    string line;
                    getline(index_File, line);

                    vector<string> line_Data;

                    while (getline(index_File, line))
                    {
                        function.split(line_Data, line, '\t');
                        if (line_Data[1] == node_ID)
                        {
                            node_Index = stoi(line_Data[0]);
                            break;
                        }
                    }

                    index_File.close();

                    if (node_Index != -1)
                    {
                        cout << "Node index: " << node_Index << endl;
                        node_intermediary_location = intermediate_Folder_location + "/sequence_Data/" + to_string(node_Index) + "/" + to_string(tissue_Index);
                        node_results_Location = output_Folder_location + "/node_Data/" + node_ID;

                        cout << "\nProcessing query sequences\n";
                        cutoff = Parameters.get_FLOAT(found_Parameters[8]);
                        cout << "cutoff Percentage: " << cutoff << endl;
                        cutoff_Count = total_Bases * cutoff;
                        cout << "cutoff Count: " << cutoff_Count << endl;

                        string query_Sequence_Location = Parameters.get_STRING(found_Parameters[5]);
                        cout << "Reading sequences file: " << query_Sequence_Location << endl;

                        fstream sequences_File;
                        sequences_File.open(query_Sequence_Location, ios::in);

                        if (sequences_File.is_open())
                        {
                            string line;
                            while (getline(sequences_File, line))
                            {
                                if (line.at(0) == '>')
                                {
                                    string line_Fix = function.clean_Invisible(line.substr(1, line.length()));
                                    sequences_to_Check_ID.push_back(function.clean_Line(line_Fix));
                                }
                                else
                                {
                                    string n_Line = line;
                                    for (int base = 0; base < line.size(); base++)
                                    {
                                        if (n_Line.at(base) == 'A' || n_Line.at(base) == 'a')
                                        {
                                            n_Line.at(base) = '0';
                                        }
                                        else if (n_Line.at(base) == 'T' || n_Line.at(base) == 't')
                                        {
                                            n_Line.at(base) = '1';
                                        }
                                        else if (n_Line.at(base) == 'G' || n_Line.at(base) == 'g')
                                        {
                                            n_Line.at(base) = '2';
                                        }
                                        else if (n_Line.at(base) == 'C' || n_Line.at(base) == 'c')
                                        {
                                            n_Line.at(base) = '3';
                                        }
                                    }

                                    sequences_to_Check.push_back(n_Line);
                                }
                            }
                            sequences_File.close();

                            if (sequences_to_Check.size() == sequences_to_Check_ID.size())
                            {
                                cout << "Sequences to check: " << sequences_to_Check.size() << endl;
                                cout << "\nMulti read: " << Multi_Read << endl;
                                cout << "CPU cores being used: " << CPU_cores << endl;

                                cout << "\nConfiguring output folders: " << endl;
                                function.config_Folder(node_results_Location + "/seg_Match", "Master output");

                                for (int sequence = 0; sequence < sequences_to_Check.size(); sequence++)
                                {
                                    function.config_Folder(node_results_Location + "/seg_Match/" + sequences_to_Check_ID[sequence], sequences_to_Check_ID[sequence]);
                                    function.create_File(node_results_Location + "/seg_Match/" + sequences_to_Check_ID[sequence] + "/" + tissue_Name + "_" + to_string(cutoff) + ".csv", "Target_sequence\tTissue\tGeneration\tquery_Sequence_ID\tmatching_Percentage\tmatch_Count\tmismatch_Count\tMismatch_bases");
                                }
                            }
                            else
                            {
                                cout << "IDs and sequences do not match in count:" << sequences_to_Check_ID.size() << "\t" << sequences_to_Check.size() << endl;
                                exit(-1);
                            }
                        }
                        else
                        {
                            cout << "UNABLE TO OPEN QUERY SEQUENCE FILE: " << query_Sequence_Location << endl;
                            exit(-1);
                        }
                    }
                    else
                    {
                        cout << "UNABLE TO FIND NODE: " << node_ID << endl;
                        exit(-1);
                    }
                }
                else
                {
                    cout << "UNABLE TO OPEN NODE INDEX FILE: " << index_Node_File << endl;
                    exit(-1);
                }
            }
            else
            {
                cout << "TISSUE " << tissue_Name << " NOT FOUND\n";
                exit(-1);
            }
        }
        else
        {
            cout << "HAS TO HAVE AT LEAST 1 TISSUE\n";
            exit(-1);
        }
    }
    else
    {
        cout << "HAS TO HAVE AT LEAST 1 MUTATION HOTSPOT REGION\n";
        exit(-1);
    }
}

void segmatch::ingress()
{
    functions_library functions = functions_library();
    vector<string> line_Data;

    cout << "\nIntializing generation detection\n";
    vector<pair<int, string>> generations_Paths;

    for (const auto &entry : filesystem::directory_iterator(node_intermediary_location))
    {
        if (filesystem::is_directory(entry))
        {
            // cout << entry.path().string() << endl;
            string gen_Name = entry.path().filename().string();
            functions.split(line_Data, gen_Name, '_');
            generations_Paths.push_back(make_pair(stoi(line_Data[1]), entry.path().string()));
        }
    }

    sort(generations_Paths.begin(), generations_Paths.end());

    int num_per_Core = sequences_to_Check_ID.size() / CPU_cores;
    int remainder = sequences_to_Check_ID.size() % CPU_cores;

    for (int generation = 0; generation < generations_Paths.size(); generation++)
    {
        cout << "\nProcessing generation " << generation + 1 << " of " << generations_Paths.size() << endl;
        cout << "Current generation: " << generations_Paths[generation].first << endl;

        vector<pair<int, int>> nFASTA_files = functions.index_Source_folder(generations_Paths[generation].second);

        cout << "Detecting haplotypes and their counts: \n";

        for (int nFASTA_file = 0; nFASTA_file < nFASTA_files.size(); nFASTA_file++)
        {
            cout << "Processing sequence file " << nFASTA_file + 1 << " of " << nFASTA_files.size() << endl;

            string n_FASTA_location = generations_Paths[generation].second + "/" + to_string(nFASTA_files[nFASTA_file].first) + "_" + to_string(nFASTA_files[nFASTA_file].second) + ".nfasta";
            cout << "Reading file: " << n_FASTA_location << endl;

            fstream nFASTA;
            nFASTA.open(n_FASTA_location, ios::in);

            vector<string> headers;
            vector<string> lines;

            if (nFASTA.is_open())
            {
                string line;
                string header = "";
                while (getline(nFASTA, line))
                {
                    if (line.at(0) == '>')
                    {
                        string header = functions.clean_Invisible(line.substr(1, line.length()));
                        header = functions.clean_Line(header);
                        functions.split(line_Data, header, '_');
                        headers.push_back(functions.clean_Line(line_Data[0]));
                    }
                    else
                    {
                        lines.push_back(line);
                    }
                }
                nFASTA.close();
            }
            else
            {
                cout << "UNABLE TO OPEN nFASTA FILE: " << n_FASTA_location << endl;
                exit(-1);
            }

            cout << "Processing file:" << endl;

            for (int line_Num = 0; line_Num < lines.size(); line_Num++)
            {
                string line = lines[line_Num];
                vector<thread> threads_vec;

                for (int core_ID = 0; core_ID < CPU_cores; core_ID++)
                {
                    int start_Node = core_ID * num_per_Core;
                    int stop_Node = start_Node + num_per_Core;

                    threads_vec.push_back(thread{&segmatch::get_Match, this, start_Node, stop_Node, ref(line), headers[line_Num], to_string(generations_Paths[generation].first)});
                }

                if (remainder != 0)
                {
                    int start_Node = sequences_to_Check_ID.size() - remainder;
                    int stop_Node = sequences_to_Check_ID.size();

                    threads_vec.push_back(thread{&segmatch::get_Match, this, start_Node, stop_Node, ref(line), headers[line_Num], to_string(generations_Paths[generation].first)});
                }

                for (thread &t : threads_vec)
                {
                    if (t.joinable())
                    {
                        t.join();
                    }
                }

                threads_vec.clear();
            }
            // exit(-1);
        }
    }
    cout << "\nDone\n";
}

void segmatch::get_Match(int start, int stop, string &sequence_Query, string header, string generation)
{
    for (int sequence = start; sequence < stop; sequence++)
    {
        int match_Count = 0;
        vector<pair<int, string>> mismatch_Bases;

        for (int site = 0; site < positions_start_end.size(); site++)
        {
            for (int base = positions_start_end[site].first; base <= positions_start_end[site].second; base++)
            {
                // cout << sequences_to_Check[sequence].at(base) << "\t" << sequence_Query.at(base) << endl;
                if (sequences_to_Check[sequence].at(base) == sequence_Query.at(base))
                {
                    match_Count++;
                }
                else
                {
                    // cout << sequences_to_Check[sequence].at(base) << "\t" << sequence_Query.at(base) << "\t" << to_string(sequence_Query.at(base)) << endl;
                    if (sequence_Query.at(base) == '0')
                    {
                        mismatch_Bases.push_back(make_pair(base, "A"));
                    }
                    else if (sequence_Query.at(base) == '1')
                    {
                        mismatch_Bases.push_back(make_pair(base, "T"));
                    }
                    else if (sequence_Query.at(base) == '2')
                    {
                        mismatch_Bases.push_back(make_pair(base, "G"));
                    }
                    else if (sequence_Query.at(base) == '3')
                    {
                        mismatch_Bases.push_back(make_pair(base, "C"));
                    }
                }
            }
        }

        float match_Percentage = (float)match_Count / (float)total_Bases;

        if (match_Percentage >= cutoff)
        {
            fstream write_File;
            write_File.open(node_results_Location + "/seg_Match/" + sequences_to_Check_ID[sequence] + "/" + tissue_Name + "_" + to_string(cutoff) + ".csv", ios::app);
            cout << "Writing to file: " << node_results_Location + "/seg_Match/" + sequences_to_Check_ID[sequence] + "/" + tissue_Name + "_" + to_string(cutoff) + ".csv" << endl;

            if (write_File.is_open())
            {
                write_File << sequences_to_Check_ID[sequence] << "\t" << tissue_Name << "\t" << generation << "\t"
                           << node_ID << "_" << tissue_Name << "_" << generation << "_" << header << "\t" << match_Percentage << "\t" << match_Count << "\t" << to_string(mismatch_Bases.size()) << "\t";
                if (mismatch_Bases.size() > 0)
                {
                    for (int mismatch = 0; mismatch < mismatch_Bases.size(); mismatch++)
                    {
                        if (mismatch != 0)
                        {
                            write_File << "|";
                        }
                        write_File << to_string(mismatch_Bases[mismatch].first + 1) << ":" << mismatch_Bases[mismatch].second;
                    }
                }
                else
                {
                    write_File << "NA";
                }

                write_File << endl;

                write_File.close();
            }
            // cout << sequences_to_Check_ID[sequence] << endl;
            // cout << match_Count << "/" << total_Bases << endl;
            // cout << match_Percentage << endl;
        }
        // else
        // {
        //     cout << match_Count << endl;
        // }
    }
}