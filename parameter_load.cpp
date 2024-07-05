#include "parameter_load.h"

parameter_load::parameter_load()
{
}

parameter_load::parameter_load(string file_location)
{
    this->file_location = file_location;
}

vector<string> parameter_load::get_parameters(string file_Location, vector<string> &parameters_List)
{
    vector<pair<string, int>> parameter_List_to_Find;
    vector<string> found_Parameters;
    int number_to_Find = parameters_List.size();

    for (size_t i = 0; i < number_to_Find; i++)
    {
        parameter_List_to_Find.push_back(make_pair(parameters_List[i], i));
        found_Parameters.push_back("");
    }

    functions_library function = functions_library();

    fstream parameter_File;
    parameter_File.open(file_Location, ios::in);

    if (parameter_File.is_open())
    {
        string line;
        getline(parameter_File, line);

        vector<string> line_Data;

        while (getline(parameter_File, line))
        {
            if (line != "}" && line != "")
            {
                string remove = line;
                int i = 0;
                while (remove[i] == ' ')
                {
                    i++; // Skip leading spaces
                }
                remove.erase(0, i);

                if (remove.at(0) != '#')
                {
                    if (remove.at(remove.size() - 1) == ',')
                    {
                        remove = remove.substr(0, remove.length() - 1);
                    }

                    // cout << remove << endl;
                    function.split(line_Data, remove, ':');

                    for (int check_Pos = 0; check_Pos < parameter_List_to_Find.size(); check_Pos++)
                    {
                        if (line_Data[0] == parameter_List_to_Find[check_Pos].first)
                        {
                            found_Parameters[parameter_List_to_Find[check_Pos].second] = line_Data[1];
                            parameter_List_to_Find.erase(parameter_List_to_Find.begin() + check_Pos);

                            break;
                        }
                    }
                }
            }
            if (parameter_List_to_Find.size() == 0)
            {
                break;
            }
        }
        parameter_File.close();
    }

    // check if all were found
    for (size_t i = 0; i < found_Parameters.size(); i++)
    {
        // cout << found_Parameters[i] << endl;
        if (found_Parameters[i] == "")
        {
            cout << "ERROR:\nCHECK FILE " << file_Location << "\nThe following parameter is missing: " << parameters_List[i] << endl;
            exit(-1);
        }
    }

    return found_Parameters;
}

void parameter_load::get_parameters(int &CUDA_device_ID, string &parent_SEQ_folder,
                                    float &mean_rep_time, float &standard_deviation_rep_time,
                                    float &mean_days_host, float &standard_deviation_host_time)
{
    functions_library function = functions_library();

    fstream parameter_File;
    parameter_File.open(this->file_location, ios::in);

    if (parameter_File.is_open())
    {
        string line;
        getline(parameter_File, line);

        vector<string> line_Data;

        while (getline(parameter_File, line))
        {
            if (line != "}" && line != "")
            {
                string remove = line;
                int i = 0;
                while (remove[i] == ' ')
                {
                    i++; // Skip leading spaces
                }
                remove.erase(0, i);

                if (remove.at(0) != '#')
                {
                    if (remove.at(remove.size() - 1) == ',')
                    {
                        remove = remove.substr(0, remove.length() - 1);
                    }

                    // cout << remove << endl;
                    function.split(line_Data, remove, ':');

                    if (line_Data[0] == "\"CUDA Device ID\"")
                    {
                        // cout << line_Data[1] << endl;
                        CUDA_device_ID = get_INT(line_Data[1]);
                        // cout << CUDA_device_ID << endl;
                    }
                    else if (line_Data[0] == "\"Parent sequences folder\"")
                    {
                        // cout << line_Data[1] << endl;
                        parent_SEQ_folder = get_STRING(line_Data[1]);
                        // cout << parent_SEQ_folder << endl;
                    }
                    else if (line_Data[0] == "\"Mean replication time\"")
                    {
                        mean_rep_time = get_FLOAT(line_Data[1]);
                    }
                    else if (line_Data[0] == "\"Standard_deviation replication time\"")
                    {
                        standard_deviation_rep_time = get_FLOAT(line_Data[1]);
                    }
                    else if (line_Data[0] == "\"Mean days in host\"")
                    {
                        mean_days_host = get_FLOAT(line_Data[1]);
                    }
                    else if (line_Data[0] == "\"Standard_deviation in host\"")
                    {
                        standard_deviation_host_time = get_FLOAT(line_Data[1]);
                    }
                }
            }
        }

        parameter_File.close();
    }
}

float parameter_load::get_FLOAT(string value)
{
    return stof(value.substr(1, value.length() - 2));
}

int parameter_load::get_INT(string value)
{
    return stoi(value);
}

string parameter_load::get_STRING(string value)
{
    return (value.substr(1, value.length() - 2));
}

float **parameter_load::get_Profile_Array(string profile_File_location, int &num_Sites, functions_library &function)
{
    float **profile_Array;

    num_Sites = 0;

    fstream profile_File;
    profile_File.open(profile_File_location, ios::in);

    if (profile_File.is_open())
    {
        cout << "Reading profile file: " << profile_File_location << endl;
        string line;
        getline(profile_File, line);

        vector<string> line_Data_Positions;

        char split_Char = ',';

        // getline(profile_File, line);
        // cout << line << endl;

        function.split(line_Data_Positions, line, split_Char);

        if (line_Data_Positions.size() == 1)
        {
            split_Char = '\t';
            function.split(line_Data_Positions, line, split_Char);
        }

        if (line_Data_Positions.size() == 1)
        {
            cout << "ERROR: PROFILE FILE SHOULD BE COMMA(,) OR TAB DELIMITED.\n";
            exit(-1);
        }

        vector<int> positions;
        num_Sites = line_Data_Positions.size() - 1;

        for (int column = 1; column < line_Data_Positions.size(); column++)
        {
            positions.push_back(stoi(line_Data_Positions[column]));
        }

        sort(positions.begin(), positions.end());

        // for (const auto &pair : positions)
        // {
        //     cout << pair.first << ", " << pair.second << endl;
        // }

        profile_Array = function.create_FLOAT_2D_arrays(num_Sites, 5);

        for (int position = 0; position < positions.size(); position++)
        {
            profile_Array[position][0] = positions[position];
        }

        vector<string> line_Data;

        while (getline(profile_File, line))
        {
            function.split(line_Data, line, split_Char);

            int base_Index = function.get_base_Index(line_Data[0]) + 1;

            for (int column = 1; column < line_Data.size(); column++)
            {
                profile_Array[function.binary_Search(positions, stoi(line_Data_Positions[column]))][base_Index] = stof(line_Data[column]);
            }
        }

        profile_File.close();
    }
    else
    {
        cout << "ERROR: UNABLE TO OPEN PROFILE FILE: " << profile_File_location << endl;
        exit(-1);
    }

    return profile_Array;
}

vector<pair<int, vector<float>>> parameter_load::get_recombination_Hotspot_Parameters(string parameter, string hotspot, string path_to_file, functions_library &function, int &total)
{
    vector<pair<int, vector<float>>> recom_matrix;

    fstream parameter_File;
    parameter_File.open(path_to_file, ios::in);

    char split_Char = ',';

    if (parameter_File.is_open())
    {
        string line;
        getline(parameter_File, line);

        string line_clean = function.clean_Line(line);

        if (line_clean.at(0) == '#')
        {
            vector<string> line_Data;
            function.split(line_Data, line_clean, split_Char);

            vector<string> get_Hotspot_ID;

            if (line_Data.size() > 0)
            {
                function.split(get_Hotspot_ID, line_Data[0], ':');
            }
            else
            {
                split_Char = '\t';
                function.split(line_Data, line_clean, split_Char);
                if (line_Data.size() > 0)
                {
                    function.split(get_Hotspot_ID, line_Data[0], ':');
                }
                else
                {
                    cout << "ERROR IN RECOMBINATION PARAMETER FILE'S, HAS TO BE TAB OR COMMA DELIMITED: " << path_to_file << endl;
                    exit(-1);
                }
            }

            if (get_Hotspot_ID[1] == hotspot)
            {
                transform(parameter.begin(), parameter.end(), parameter.begin(), ::toupper);
                cout << "Extracting " << parameter << " mutation data\nFrom: " << path_to_file << endl;

                int catch_Parameter = 0;

                while (getline(parameter_File, line))
                {
                    line_clean = function.clean_Line(line);
                    function.split(line_Data, line_clean, split_Char);
                    transform(line_Data[0].begin(), line_Data[0].end(), line_Data[0].begin(), ::toupper);
                    if (line_Data[0] == parameter)
                    {
                        catch_Parameter = 1;
                        break;
                    }
                }

                if (catch_Parameter == 1)
                {
                    // int sites_Total = 0;
                    vector<int> positions;
                    for (int check_Data = 1; check_Data < line_Data.size(); check_Data++)
                    {
                        if (line_Data[check_Data].size() > 0)
                        {
                            positions.push_back(stoi(line_Data[check_Data]));
                            // sites_Total++;
                        }
                        else
                        {
                            break;
                        }
                    }

                    cout << positions.size() << " mutation effect sites found" << endl;
                    // if (positions.size() > max)
                    // {
                    total = total + positions.size();
                    //  }
                    for (int site = 0; site < positions.size(); site++)
                    {
                        vector<float> bases;
                        for (int base = 0; base < 4; base++)
                        {
                            bases.push_back(0);
                        }

                        recom_matrix.push_back(make_pair(positions[site], bases));
                    }

                    for (int base = 0; base < 4; base++)
                    {
                        getline(parameter_File, line);
                        line_clean = function.clean_Line(line);

                        function.split(line_Data, line_clean, split_Char);

                        int base_Current = -1;
                        if (line_Data[0] == "A" || line_Data[0] == "a")
                        {
                            base_Current = 0;
                        }
                        else if (line_Data[0] == "T" || line_Data[0] == "t")
                        {
                            base_Current = 1;
                        }
                        else if (line_Data[0] == "G" || line_Data[0] == "g")
                        {
                            base_Current = 2;
                        }
                        else if (line_Data[0] == "C" || line_Data[0] == "c")
                        {
                            base_Current = 3;
                        }
                        else
                        {
                            cout << "ERROR IN RECOMBINATION PARAMETER FILE, BASES SHOULD BE A, T. G OR C. NOT \"" << line_Data[0] << "\": " << path_to_file << endl;
                            exit(-1);
                        }

                        for (int position = 0; position < positions.size(); position++)
                        {
                            recom_matrix[position].second[base_Current] = stof(line_Data[position + 1]);
                        }
                    }
                }
                else
                {
                    cout << "ERROR IN RECOMBINATION PARAMETER FILE, " << parameter << " NOT FOUND : " << path_to_file << endl;
                    exit(-1);
                }
            }
        }
        else
        {
            cout << "ERROR IN RECOMBINATION PARAMETER FILE, LINE ONE HAS TO START WITH A # AND STATE THE REOCMBINATION HOTSPOT NUMBER: " << path_to_file << endl;
            exit(-1);
        }
        parameter_File.close();
    }
    else
    {
        cout << "ERROR IN RECOMBINATION PARAMETER FILE, CHECK IF FILE EXISTS: " << path_to_file << endl;
        exit(-1);
    }

    return recom_matrix;
}

vector<pair<string, string>> parameter_load::get_block_from_File(string &file_parameter_Location, string block_Header)
{
    vector<pair<string, string>> block_Data;
    int activate_Collection = 0;

    functions_library function = functions_library();

    fstream parameter_File;
    parameter_File.open(file_parameter_Location, ios::in);

    if (parameter_File.is_open())
    {
        string line;
        getline(parameter_File, line);

        vector<string> line_Data;

        while (getline(parameter_File, line))
        {
            line_Data = clean_Line(line, function);

            if (line_Data.size() != 0)
            {
                if (line_Data[0] == "\"" + block_Header + "\"")
                {
                    activate_Collection = 1;
                    // cout << line_Data[0] << endl;
                    break;
                }
            }
        }

        if (activate_Collection == 1)
        {
            int count_bracket = 0;

            while (getline(parameter_File, line))
            {
                line_Data = clean_Line(line, function);

                if (line_Data.size() != 0)
                {
                    if (line_Data.size() > 1)
                    {
                        // cout << line << endl;
                        // cout << line_Data[1] << endl;
                        // exit(-1);
                        if (line_Data[1] == "{")
                        {
                            count_bracket++;
                        }

                        if (count_bracket >= 0)
                        {
                            block_Data.push_back(make_pair(line_Data[0], line_Data[1]));
                        }
                        else
                        {
                            break;
                        }
                    }
                    else if (line_Data[0] == "}")
                    {
                        count_bracket--;

                        if (count_bracket >= 0)
                        {
                            block_Data.push_back(make_pair(line_Data[0], ""));
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }
        }
        else
        {
            cout << "SYSTEM ERROR: " << block_Header << " DOES NOT EXIST IN THE CURRENT BLOCK.\n";
            exit(-1);
        }

        parameter_File.close();
    }

    return block_Data;
}

vector<pair<string, string>> parameter_load::get_block_from_block(vector<pair<string, string>> &block, string block_Header)
{
    vector<pair<string, string>> block_Data;

    int catch_Index = -1;

    for (int check = 0; check < block.size(); check++)
    {
        if (block[check].first == "\"" + block_Header + "\"")
        {
            catch_Index = check + 1;
            break;
        }
    }

    if (catch_Index != -1)
    {
        int count_bracket = 0;

        for (int check = catch_Index; check < block.size(); check++)
        {
            if (block[check].second == "{")
            {
                count_bracket++;
            }
            else if (block[check].first == "}")
            {
                count_bracket--;
            }

            if (count_bracket >= 0)
            {
                block_Data.push_back(make_pair(block[check].first, block[check].second));
            }
            else
            {
                break;
            }
        }
    }
    else
    {
        cout << "SYSTEM ERROR: " << block_Header << " DOES NOT EXIST IN THE CURRENT BLOCK.\n";
        exit(-1);
    }

    return block_Data;
}

vector<pair<string, string>> parameter_load::check_block_from_block(vector<pair<string, string>> &block, string block_Header)
{
    vector<pair<string, string>> block_Data;

    int catch_Index = -1;

    for (int check = 0; check < block.size(); check++)
    {
        if (block[check].first == "\"" + block_Header + "\"")
        {
            catch_Index = check + 1;
            break;
        }
    }

    if (catch_Index != -1)
    {
        int count_bracket = 0;

        for (int check = catch_Index; check < block.size(); check++)
        {
            if (block[check].second == "{")
            {
                count_bracket++;
            }
            else if (block[check].first == "}")
            {
                count_bracket--;
            }

            if (count_bracket >= 0)
            {
                block_Data.push_back(make_pair(block[check].first, block[check].second));
            }
            else
            {
                break;
            }
        }
    }

    return block_Data;
}

int parameter_load::get_INT(vector<pair<string, string>> block, string value)
{
    int integer_Value;
    int found = -1;

    for (int i = 0; i < block.size(); i++)
    {
        if (block[i].first == "\"" + value + "\"")
        {
            integer_Value = stoi(block[i].second);
            found = 1;
            break;
        }
    }

    if (found == 1)
    {
        return integer_Value;
    }
    else
    {
        cout << "SYSTEM ERROR: " << value << " DOES NOT EXIST IN THE CURRENT BLOCK.\n";
        exit(-1);
    }
}

float parameter_load::get_FLOAT(vector<pair<string, string>> block, string value)
{
    float return_Value;
    int found = -1;

    for (int i = 0; i < block.size(); i++)
    {
        if (block[i].first == "\"" + value + "\"")
        {
            return_Value = stof(block[i].second.substr(1, block[i].second.length() - 2));
            found = 1;
            break;
        }
    }

    if (found == 1)
    {
        return return_Value;
    }
    else
    {
        cout << "SYSTEM ERROR: " << value << " DOES NOT EXIST IN THE CURRENT BLOCK.\n";
        exit(-1);
    }
}

string parameter_load::get_STRING(vector<pair<string, string>> block, string value)
{
    string return_Value;
    int found = -1;

    for (int i = 0; i < block.size(); i++)
    {
        if (block[i].first == "\"" + value + "\"")
        {
            return_Value = block[i].second.substr(1, block[i].second.length() - 2);
            found = 1;
            break;
        }
    }

    if (found == 1)
    {
        return return_Value;
    }
    else
    {
        cout << "SYSTEM ERROR: " << value << " DOES NOT EXIST IN THE CURRENT BLOCK.\n";
        exit(-1);
    }
}

vector<string> parameter_load::clean_Line(string line, functions_library &function)
{
    vector<string> line_Data;

    if (line != "")
    {
        // string trim_Line = line;
        int i = 0;
        while (line[i] == ' ')
        {
            i++; // Skip leading spaces
        }

        line.erase(0, i);

        if (line.at(0) != '#')
        {
            if (line.at(line.size() - 1) == ',')
            {
                line = line.substr(0, line.length() - 1);
            }
            function.split(line_Data, line, ':');
        }
    }

    return line_Data;
}