#include "bfs.cuh"

bfs::bfs(string parameter_Master_Location)
{
    cout << "Initiating Breath First Search to identify pedigree of a sequence\n";

    parameter_load Parameters = parameter_load();
    functions_library function = functions_library();

    // exit(-1);

    vector<string> parameters_List = {
        "\"Intermediate folders\"",
        "\"Output folders\"",
        "\"Pedigree Node ID\"",
        "\"Pedigree Tissue\"",
        "\"Pedigree Generation\"",
        "\"Pedigree Sequence\"",
        "\"Nodes master profile\"",
        "\"Pedigree get Sequences\""};

    vector<string> found_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List);

    string nsequence = "";
    int num_tissues_per_Node = 0;

    this->intermediate_Folder_location = Parameters.get_STRING(found_Parameters[0]);
    this->output_Folder_location = Parameters.get_STRING(found_Parameters[1]);
    this->print_Sequences = function.to_Upper_Case(Parameters.get_STRING(found_Parameters[7]));

    cout << "\nReading target sequence\n";
    string pedigree_Sequence_loation = Parameters.get_STRING(found_Parameters[5]);
    fstream pedigree_File;
    pedigree_File.open(pedigree_Sequence_loation, ios::in);

    if (pedigree_File.is_open())
    {
        string line;
        // skip first line;
        getline(pedigree_File, line);
        // sequence line
        getline(pedigree_File, line);
        // cout << line << endl;
        nsequence = line;
        pedigree_File.close();
    }
    else
    {
        cout << "ERROR: UNABLE TO OPEN TARGET SEQUENCE FILE: " << pedigree_Sequence_loation << endl;
        exit(-1);
    }

    for (int base = 0; base < nsequence.size(); base++)
    {
        if (nsequence[base] == 'A' || nsequence[base] == 'a')
        {
            nsequence[base] = '0';
        }
        else if (nsequence[base] == 'T' || nsequence[base] == 't')
        {
            nsequence[base] = '1';
        }
        else if (nsequence[base] == 'G' || nsequence[base] == 'g')
        {
            nsequence[base] = '2';
        }
        else if (nsequence[base] == 'C' || nsequence[base] == 'c')
        {
            nsequence[base] = '3';
        }
        else if (nsequence[base] != '0' && nsequence[base] != '1' && nsequence[base] != '2' && nsequence[base] != '3')
        {
            cout << "ERROR: UNRECOGNISED BASE IN SEQUENCE: " << nsequence[base];
            exit(-1);
        }
    }

    cout << "Target sequence loaded\n";
    // cout << nsequence << endl;

    string node_ID = Parameters.get_STRING(found_Parameters[2]);
    cout << "\nGetting index of target node: " << node_ID << "\n";

    string node_Index_file_location = intermediate_Folder_location + "/index_Data/node_Index.csv";
    fstream node_index_File;
    node_index_File.open(node_Index_file_location, ios::in);

    if (node_index_File.is_open())
    {
        string line;
        vector<string> line_Data;

        // skip first header line
        getline(node_index_File, line);

        while (getline(node_index_File, line))
        {
            function.split(line_Data, line, '\t');
            node_Indexes.push_back(make_pair(stoi(line_Data[0]), line_Data[1]));
            if (line_Data[1] == node_ID)
            {
                node_main_Index = line_Data[0];
                // break;
            }
        }
        node_index_File.close();
    }
    else
    {
        cout << "ERROR: UNABLE TO OPEN NODE INDEX FILE: " << node_Index_file_location;
        exit(-1);
    }

    if (node_main_Index != "")
    {
        cout << "Target node's index: " << node_main_Index << endl;
    }
    else
    {
        cout << "ERROR: UNABLE TO FIND THE NODE ID: " << node_ID << endl;
        exit(-1);
    }

    sort(node_Indexes.begin(), node_Indexes.end());

    string tissue_Name = Parameters.get_STRING(found_Parameters[3]);
    cout << "\nGetting index of tissue: " << tissue_Name << endl;

    string node_Master_location = Parameters.get_STRING(found_Parameters[6]);
    cout << "Node master file: " << node_Master_location << endl;
    vector<pair<string, string>> Tissue_profiles_block_Data = Parameters.get_block_from_File(node_Master_location, "Tissue profiles");
    num_tissues_per_Node = Parameters.get_INT(Tissue_profiles_block_Data, "Number of tissues");

    if (num_tissues_per_Node > 0)
    {
        cout << "\nNumber of tissues in a node: " << num_tissues_per_Node << endl;

        for (int tissue = 0; tissue < num_tissues_per_Node; tissue++)
        {
            string check_Tissue = "Tissue " + to_string(tissue + 1) + " Name";
            tissue_Names.push_back(Parameters.get_STRING(Tissue_profiles_block_Data, check_Tissue));
            cout << check_Tissue << ": " << tissue_Names[tissue] << endl;

            if (tissue_Names[tissue] == tissue_Name)
            {
                tissue_main_Index = tissue;
            }
        }
    }
    else
    {
        cout << "ERROR: TISSUE NUMBER HAS TO BE GREATER THAN ZERO.\n\n";
    }

    if (tissue_main_Index != -1)
    {
        cout << "\nTissue index of " << tissue_Name << " : " << tissue_main_Index << endl;
    }
    else
    {
        cout << "UNABLE TO FIND TISSUE INDEX: " << tissue_Name << endl;
        exit(-1);
    }

    generation = Parameters.get_INT(found_Parameters[4]);

    cout << "\nIndentifying matching sequences from target\n";
    string sequence_Search_folder = intermediate_Folder_location + "/sequence_Data/" + node_main_Index;

    re_tar_sequence_Folder = check_Tar_Folder(sequence_Search_folder);

    sequence_Search_folder = sequence_Search_folder + "/" + to_string(tissue_main_Index);

    re_tar_Tissue_folder = check_Tar_Folder(sequence_Search_folder);

    sequence_Search_folder = sequence_Search_folder + "/generation_" + to_string(generation);

    re_tar_Generation = check_Tar_Folder(sequence_Search_folder);

    for (const auto &entry : filesystem::directory_iterator(sequence_Search_folder))
    {
        if (filesystem::is_regular_file(entry) && entry.path().extension() == ".nfasta")
        {
            string check_File_location = entry.path().string();
            cout << "Checking nfasta file: " << check_File_location << endl;

            fstream nfasta_File;
            nfasta_File.open(check_File_location, ios::in);
            if (nfasta_File.is_open())
            {
                string line;
                string header;
                while (getline(nfasta_File, line))
                {
                    if (line.at(0) == '>')
                    {
                        header = line;
                    }
                    else
                    {
                        if (line == nsequence)
                        {
                            header = header.substr(1);
                            header = header.substr(0, header.find('_'));
                            search_sequence_IDs.push_back(node_ID + "_" + tissue_Name + "_" + to_string(generation) + "_" + header);
                        }
                    }
                }

                nfasta_File.close();
            }
            else
            {
                cout << "ERROR UNABLE TO OPEN NFASTA FILE: " << check_File_location << endl;
                exit(-1);
            }
        }
    }

    if (search_sequence_IDs.size() > 0)
    {
        cout << "\nFound " << search_sequence_IDs.size() << " matching sequence(s)\n";
        // this->current_node_ID = node_ID;
        //  cout << search_sequence_IDs[0] << endl;

        if (!filesystem::exists(this->output_Folder_location + "/node_Data/" + node_ID + "/pedigree_Information"))
        {
            filesystem::create_directory(this->output_Folder_location + "/node_Data/" + node_ID + "/pedigree_Information");
        }

        if (!filesystem::exists(this->output_Folder_location + "/node_Data/" + node_ID + "/pedigree_Information/" + filesystem::path(pedigree_Sequence_loation).stem().string()))
        {
            filesystem::create_directory(this->output_Folder_location + "/node_Data/" + node_ID + "/pedigree_Information/" + filesystem::path(pedigree_Sequence_loation).stem().string());
        }

        pedigree_Folder_location = this->output_Folder_location + "/node_Data/" + node_ID + "/pedigree_Information/" + filesystem::path(pedigree_Sequence_loation).stem().string();
    }
    else
    {
        cout << "ERROR SEQUENCE NOT FOUND\n";
        exit(-1);
    }
}

void bfs::ingress()
{
    functions_library functions = functions_library();

    for (int sequence = 0; sequence < search_sequence_IDs.size(); sequence++)
    {
        string ID_Sequence_Original = search_sequence_IDs[sequence];
        cout << "\nIdentifying pedigree of sequence: " << ID_Sequence_Original << endl;

        vector<string> sequence_Information;
        functions.split(sequence_Information, ID_Sequence_Original, '_');
        string current_node_ID = sequence_Information[0] + "_" + sequence_Information[1];

        vector<pair<string, string>> queue;
        queue.push_back(make_pair(sequence_Information[0] + "_" + sequence_Information[1], ID_Sequence_Original));
        int queue_Count = 0;

        functions.create_File(pedigree_Folder_location + "/" + ID_Sequence_Original + "_pedigree_Relationships.csv", "Source\tTarget\tType");
        functions.create_File(pedigree_Folder_location + "/" + ID_Sequence_Original + "_sequence_Information.csv", "ID\tHost\tTissue\tGeneration");

        set<string> sequence_IDs;
        sequence_IDs.insert(ID_Sequence_Original);

        vector<string> line_Data;

        do
        {
            string ID_Sequence = queue[queue_Count].second;
            string line;

            fstream node_File;
            node_File.open(this->output_Folder_location + "/node_Data/" + queue[queue_Count].first + "/sequence_parent_Progeny_relationships.csv",
                           ios::in);

            // initialize this using node generational_Summary
            vector<pair<string, string>> progeny_Parent;
            vector<pair<string, int>> generation_Line;
            vector<string> Type;

            int collect = -1;

            int line_Count = 0;
            int hit = 0;
            // int gen_Track = 0;
            string generation_Current = "-1";

            vector<string> sequence_temp_Info;

            if (node_File.is_open())
            {
                cout << "\nReading node parent progeny file: " << queue[queue_Count].first << endl;

                // skip header
                getline(node_File, line);

                while (getline(node_File, line))
                {
                    functions.split(line_Data, line, '\t');

                    functions.split(sequence_temp_Info, line_Data[0], '_');
                    if (generation_Current == "-1")
                    {
                        generation_Current = sequence_temp_Info[3];
                        hit++;
                    }
                    if (generation_Current != sequence_temp_Info[3])
                    {
                        hit++;
                        functions.split(sequence_temp_Info, progeny_Parent[progeny_Parent.size() - 1].second, '_');
                        generation_Current = sequence_temp_Info[3];
                        if (hit == 2)
                        {
                            generation_Line.push_back(make_pair(generation_Current, line_Count - 1));
                            cout << "Generation index: " << generation_Current << " line: " << line_Count - 1 << endl;
                            hit = 0;
                        }
                    }

                    if (collect == 1 && line_Data[1] != ID_Sequence)
                    {
                        generation_Line.push_back(make_pair(generation_Current, line_Count - 1));
                        cout << "Generation index: " << generation_Current << " line: " << line_Count - 1 << endl;
                        cout << "Captured target sequence: " << ID_Sequence << endl;
                        break;
                    }
                    progeny_Parent.push_back(make_pair(line_Data[1], line_Data[0]));
                    Type.push_back(line_Data[2]);
                    if (line_Data[1] == ID_Sequence)
                    {
                        collect = 1;
                    }
                    line_Count++;
                }
                node_File.close();
            }
            else
            {
                cout << "ERROR UNABLE TO OPEN PARENT PROGENY FILE: " << this->output_Folder_location << "/node_Data/" << queue[queue_Count].first << "/sequence_parent_Progeny_relationships.csv" << endl;
                exit(-1);
            }

            cout << "Collecting parent progeny pedigree\n";

            int get_Parents = 0;
            vector<string> progeny_List;
            progeny_List.push_back(ID_Sequence);
            int track_progeny_Parent = progeny_Parent.size() - 1;

            fstream source_Target_pedigree;
            source_Target_pedigree.open(pedigree_Folder_location + "/" + ID_Sequence_Original + "_pedigree_Relationships.csv", ios::app);

            functions.split(sequence_temp_Info, ID_Sequence, '_');
           // generation_Current = sequence_temp_Info[3];

            if (source_Target_pedigree.is_open())
            {
                int catch_Check = -1;
                do
                {
                    if (progeny_Parent[track_progeny_Parent].first == progeny_List[get_Parents])
                    {
                        catch_Check = 1;
                        source_Target_pedigree << progeny_Parent[track_progeny_Parent].second << "\t"
                                               << progeny_Parent[track_progeny_Parent].first << "\t"
                                               << Type[track_progeny_Parent] << "\n";

                        // if parent not present in progeny list add it to progeny list
                        auto it = sequence_IDs.find(progeny_Parent[track_progeny_Parent].second);
                        if (it == sequence_IDs.end())
                        {
                            sequence_IDs.insert(progeny_Parent[track_progeny_Parent].second);
                            // check if new host
                            functions.split(line_Data, progeny_Parent[track_progeny_Parent].second, '_');
                            string check_Node = line_Data[0] + "_" + line_Data[1];
                            if (check_Node != queue[queue_Count].first)
                            {
                                queue.push_back(make_pair(check_Node, progeny_Parent[track_progeny_Parent].second));
                            }
                            else
                            {
                                progeny_List.push_back(progeny_Parent[track_progeny_Parent].second);
                            }
                        }
                    }
                    track_progeny_Parent--;

                    if (catch_Check == 1 && progeny_Parent[track_progeny_Parent].first != progeny_List[get_Parents])
                    {
                        get_Parents++;
                        catch_Check = -1;
                        cout << "New progeny: " << progeny_List[get_Parents] << endl;
                        functions.split(sequence_temp_Info, progeny_List[get_Parents], '_');
                        string generation_to_Check = sequence_temp_Info[3];

                        for (int gen = 0; gen < generation_Line.size(); gen++)
                        {
                            if (generation_Line[gen].first == generation_to_Check)
                            {
                                track_progeny_Parent = generation_Line[gen].second;
                                break;
                            }
                        }
                    }

                } while (get_Parents < progeny_List.size() && track_progeny_Parent >= 0);

                source_Target_pedigree.close();

                //exit(-1);
            }
            else
            {
                cout << "ERROR: UNABLE TO OPEN FILE: " << pedigree_Folder_location << "/" << ID_Sequence_Original << "_pedigree_Relationships.csv" << endl;
                exit(-1);
            }

            queue_Count++;
        } while (queue_Count < queue.size());

        cout << "Writing sequence information\n";

        vector<string> IDs_complete(sequence_IDs.begin(), sequence_IDs.end());
        fstream sequence_Information_File;
        sequence_Information_File.open(pedigree_Folder_location + "/" + ID_Sequence_Original + "_sequence_Information.csv", ios::app);
        if (sequence_Information_File.is_open())
        {
            for (string write_Line : IDs_complete)
            {
                sequence_Information_File << write_Line << "\t";
                functions.split(line_Data, write_Line, '_');
                sequence_Information_File << line_Data[0] << "_" << line_Data[1] << "\t"
                                          << line_Data[2] << "\t" << line_Data[3] << endl;
            }
            sequence_Information_File.close();
        }
        else
        {
            cout << "ERROR: UNABLE TO OPEN FILE: " << pedigree_Folder_location << "/" << ID_Sequence_Original << "_sequence_Information.csv\n";
            exit(-1);
        }

        if (print_Sequences != "NO")
        {
            cout << "\nPrinting sequences:\n";
            fstream sequence_File;
            sequence_File.open(pedigree_Folder_location + "/" + ID_Sequence_Original + "_pedigree_Sequences.fasta", ios::out);

            vector<string> node_locations_Unique;
            vector<vector<pair<int, string>>> sequence_IDs_per_Location;

            if (sequence_File.is_open())
            {
                vector<string> untar_List;
                for (int seq = 0; seq < IDs_complete.size(); seq++)
                {
                    string sequence_Label = IDs_complete[seq];
                    cout << "Processing sequence: " << sequence_Label << endl;

                    functions.split(line_Data, sequence_Label, '_');

                    int tissue_Index = -1;
                    for (int tissue = 0; tissue < tissue_Names.size(); tissue++)
                    {
                        if (line_Data[2] == tissue_Names[tissue])
                        {
                            tissue_Index = tissue;
                            break;
                        }
                    }
                    if (tissue_Index == -1)
                    {
                        cout << "ERROR: UNRECOGNISED TISSUE LABEL: " << line_Data[2] << endl;
                        exit(-1);
                    }

                    string node_ID = line_Data[0] + "_" + line_Data[1];
                    int node_Index = -1;

                    for (int node = 0; node < node_Indexes.size(); node++)
                    {
                        if (node_Indexes[node].second == node_ID)
                        {
                            node_Index = node_Indexes[node].first;
                            break;
                        }
                    }
                    if (node_Index == -1)
                    {
                        cout << "ERROR: UNABLE TO FIND NODE: " << node_ID << endl;
                        exit(-1);
                    }

                    string folder_Location = this->intermediate_Folder_location + "/sequence_Data/" + to_string(node_Index) + "/" + to_string(tissue_Index) + "/generation_" + line_Data[3];

                    int location_Check = -1;

                    for (int location = 0; location < node_locations_Unique.size(); location++)
                    {
                        if (node_locations_Unique[location] == folder_Location)
                        {
                            location_Check = location;
                            break;
                        }
                    }

                    if (location_Check == -1)
                    {
                        node_locations_Unique.push_back(folder_Location);
                        vector<pair<int, string>> init_Seq_IDs;
                        sequence_IDs_per_Location.push_back(init_Seq_IDs);
                        sequence_IDs_per_Location[sequence_IDs_per_Location.size() - 1].push_back(make_pair(stoi(line_Data[4]), sequence_Label));

                        if (check_Tar_Folder(this->intermediate_Folder_location + "/sequence_Data/" + to_string(node_Index)) == 1)
                        {
                            untar_List.push_back(this->intermediate_Folder_location + "/sequence_Data/" + to_string(node_Index));
                        }

                        if (check_Tar_Folder(this->intermediate_Folder_location + "/sequence_Data/" + to_string(node_Index) + "/" + to_string(tissue_Index)) == 1)
                        {
                            untar_List.push_back(this->intermediate_Folder_location + "/sequence_Data/" + to_string(node_Index) + "/" + to_string(tissue_Index));
                        }

                        if (check_Tar_Folder(this->intermediate_Folder_location + "/sequence_Data/" + to_string(node_Index) + "/" + to_string(tissue_Index) + "/generation_" + line_Data[3]) == 1)
                        {
                            untar_List.push_back(this->intermediate_Folder_location + "/sequence_Data/" + to_string(node_Index) + "/" + to_string(tissue_Index) + "/generation_" + line_Data[3]);
                        }
                    }
                    else
                    {
                        sequence_IDs_per_Location[location_Check].push_back(make_pair(stoi(line_Data[4]), sequence_Label));
                    }
                }

                cout << "\nExtracting sequences: \n\n";

                for (int location = 0; location < node_locations_Unique.size(); location++)
                {
                    vector<pair<int, int>> index_Folder = functions.index_Source_folder(node_locations_Unique[location]);
                    vector<pair<int, string>> sequences_to_Extract = sequence_IDs_per_Location[location];
                    sort(sequences_to_Extract.begin(), sequences_to_Extract.end());

                    int seq_Track = 0;

                    for (int folder = 0; folder < index_Folder.size(); folder++)
                    {
                        if (seq_Track >= sequences_to_Extract.size())
                        {
                            break;
                        }

                        fstream seq_nFASTA_File;
                        seq_nFASTA_File.open(node_locations_Unique[location] + "/" + to_string(index_Folder[folder].first) + "_" + to_string(index_Folder[folder].second) + ".nfasta", ios::in);
                        if (seq_nFASTA_File.is_open())
                        {
                            string line;
                            int current_Line = 0;
                            int line_t0_check = (sequences_to_Extract[seq_Track].first - index_Folder[folder].first) * 2;

                            while (getline(seq_nFASTA_File, line))
                            {
                                if ((sequences_to_Extract[seq_Track].first >= index_Folder[folder].first) && (sequences_to_Extract[seq_Track].first <= index_Folder[folder].second))
                                {
                                    if (line_t0_check == current_Line)
                                    {
                                        cout << "Writing sequence: " << sequences_to_Extract[seq_Track].second << ": ";
                                        vector<string> line_Data;
                                        functions.split(line_Data, line, '_');

                                        if (stoi(line_Data[0].substr(1)) == sequences_to_Extract[seq_Track].first)
                                        {
                                            getline(seq_nFASTA_File, line);
                                            current_Line++;
                                            sequence_File << ">" << sequences_to_Extract[seq_Track].second << endl;
                                            // sequence_File << line << endl;
                                            for (int base = 0; base < line.size(); base++)
                                            {
                                                if (line.at(base) == '0')
                                                {
                                                    sequence_File << "A";
                                                }
                                                else if (line.at(base) == '1')
                                                {
                                                    sequence_File << "T";
                                                }
                                                else if (line.at(base) == '2')
                                                {
                                                    sequence_File << "G";
                                                }
                                                else if (line.at(base) == '3')
                                                {
                                                    sequence_File << "C";
                                                }
                                            }
                                            sequence_File << endl;
                                            sequence_File.flush();

                                            cout << "DONE\n";

                                            seq_Track++;
                                            if (seq_Track >= sequences_to_Extract.size())
                                            {
                                                break;
                                            }
                                            else
                                            {
                                                cout << "Looking for sequence: " << sequences_to_Extract[seq_Track].second << "\n";
                                                line_t0_check = (sequences_to_Extract[seq_Track].first - index_Folder[folder].first) * 2;
                                            }
                                        }
                                        else
                                        {
                                            cout << "ERROR: CORRECT SEQUENCE NOT FOUND AT INDEX\n";
                                            cout << "Looking for: " << sequences_to_Extract[seq_Track].first << endl
                                                 << "Sequence ID at location: " << line << endl
                                                 << "File: " << node_locations_Unique[location] << "/" << index_Folder[folder].first << "_" << index_Folder[folder].second << ".nfasta\n";
                                            exit(-1);
                                        }
                                    }
                                }
                                else
                                {
                                    break;
                                }
                                current_Line++;
                            }
                            seq_nFASTA_File.close();
                        }
                        else
                        {
                            cout << "ERROR: UNABLE TO OPEN FILE: " << node_locations_Unique[location] << "/" << index_Folder[folder].first << "_" << index_Folder[folder].second << ".nfasta\n";
                            exit(-1);
                        }
                    }
                }

                sequence_File.close();
                cout << "\nCompleted extraction of sequences\n";

                cout << "\nPurging untar folders\n";
                for (int untar = untar_List.size() - 1; untar >= 0; untar--)
                {
                    if (filesystem::exists(untar_List[untar]))
                    {
                        filesystem::remove_all(untar_List[untar]);
                        cout << "Folder purged: " << untar_List[untar] << endl;
                    }
                }
            }
            else
            {
                cout << "ERROR: UNABLE TO CREATE SEQUENCE FILE: " << pedigree_Folder_location << "/" << ID_Sequence_Original << "_pedigree_Sequences.fasta\n";
                exit(-1);
            }
        }
        //break;
    }

    // retar code;
    string sequence_Search_folder;
    if (re_tar_Generation == 1)
    {
        sequence_Search_folder = intermediate_Folder_location + "/sequence_Data/" + node_main_Index;
        sequence_Search_folder = sequence_Search_folder + "/" + to_string(tissue_main_Index);
        sequence_Search_folder = sequence_Search_folder + "/generation_" + to_string(generation);
        functions.folder_Delete(sequence_Search_folder);
    }
    if (re_tar_Tissue_folder == 1)
    {
        sequence_Search_folder = intermediate_Folder_location + "/sequence_Data/" + node_main_Index;
        sequence_Search_folder = sequence_Search_folder + "/" + to_string(tissue_main_Index);
        functions.folder_Delete(sequence_Search_folder);
    }
    if (re_tar_sequence_Folder == 1)
    {
        sequence_Search_folder = intermediate_Folder_location + "/sequence_Data/" + node_main_Index;
        functions.folder_Delete(sequence_Search_folder);
    }

    cout << "\nPedigree extraction complete\n";
}

int bfs::check_Tar_Folder(string location)
{
    int re_Tar = -1;

    if (filesystem::exists(location + ".tar"))
    {
        cout << "Extracting folder: " << location + ".tar\n";
        string command = "tar -xf" + location + ".tar -C .";

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

        re_Tar = 1;
    }

    return re_Tar;
}
