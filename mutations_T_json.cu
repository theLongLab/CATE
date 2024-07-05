#include "mutations_T_json.cuh"

mutations_T_json::mutations_T_json(string parameter_Master_Location)
{
    cout << "Intializing conversion to JSON format\n";

    parameter_load Parameters = parameter_load();
    functions_library function = functions_library();
    
    //exit(-1);

    vector<string> parameters_List = {
        "\"Convert file location\"",
        "\"Output folders\""};

    vector<string> found_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List);

    site_model_Location = Parameters.get_STRING(found_Parameters[0]);
    output_Folder = Parameters.get_STRING(found_Parameters[1]);
    function.config_Folder(output_Folder, "Output");
}

void mutations_T_json::recombinations_Convert()
{
    cout << "\nReading recombination file: " << site_model_Location << endl;

    functions_library functions = functions_library();

    string file_Name = filesystem::path(site_model_Location).filename().stem().string();

    string output_File_Location = output_Folder + "/" + file_Name + ".json";
    cout << "Writing recombination JSON file to: " << output_File_Location << endl;

    fstream site_model_File;
    site_model_File.open(site_model_Location, ios::in);

    if (site_model_File.is_open())
    {
        string line;
        getline(site_model_File, line);

        vector<string> line_Data;
        char delim = '\t';
        functions.split(line_Data, line, delim);

        if (line_Data.size() == 1)
        {
            delim = ',';
            functions.split(line_Data, line, delim);
        }

        if (line_Data.size() == 5)
        {
            fstream json_File;
            json_File.open(output_File_Location, ios::out);

            json_File << "    \"Recombination\":{\n\n";

            vector<string> write_Lines;
            cout << "Configuring lines\n";

            while (getline(site_model_File, line))
            {
                functions.split(line_Data, line, delim);

                string write_Line = "        \"Hotspot " + line_Data[0] + "\":{\n\n" + "            \"Region\":\"" + line_Data[1] + "_" + line_Data[2] + "\",\n\n";
                write_Line = write_Line + "            \"Reference probability of recombination\":\"" + line_Data[3] + "\",\n";
                write_Line = write_Line + "            \"Reference selectivity\":\"" + line_Data[4] + "\"\n";
                write_Line = write_Line + "        },\n\n";
                write_Lines.push_back(write_Line);
            }

            cout << "Writing lines\n";

            json_File << "        \"Number of hotspots\":" << write_Lines.size() << ",\n\n";

            for (int line = 0; line < write_Lines.size(); line++)
            {
                json_File << write_Lines[line];
            }

            json_File << "    },\n";

            json_File.close();
        }
        else
        {
            cout << "ERROR: RECOMBINATION FILE DOES NOT HAVE THE CORRECT NUMBER OF COLUMNS\n\n";
            exit(-1);
        }

        site_model_File.close();
    }
}

void mutations_T_json::mutations_Convert()
{
    cout << "\nReading site model file: " << site_model_Location << endl;
    functions_library functions = functions_library();

    string file_Name = filesystem::path(site_model_Location).filename().stem().string();

    string output_File_Location = output_Folder + "/" + file_Name + ".json";
    cout << "Writing site model JSON file to: " << output_File_Location << endl;

    // exit(-1);

    fstream site_model_File;
    site_model_File.open(site_model_Location, ios::in);

    if (site_model_File.is_open())
    {
        string line;
        getline(site_model_File, line);

        vector<string> line_Data;
        char delim = '\t';
        functions.split(line_Data, line, delim);

        if (line_Data.size() == 1)
        {
            delim = ',';
            functions.split(line_Data, line, delim);
        }

        if (line_Data.size() == 24)
        {
            vector<string> base_mutations;
            cout << "\nCollecting base change headers\n\n";
            for (int col = 8; col < line_Data.size(); col++)
            {
                base_mutations.push_back(line_Data[col]);
            }

            fstream json_File;
            json_File.open(output_File_Location, ios::out);

            json_File << "    \"Mutations\":{\n\n";

            vector<string> write_Lines;
            cout << "Configuring lines\n";
            while (getline(site_model_File, line))
            {
                functions.split(line_Data, line, delim);
                string write_Line = "        \"Hotspot " + line_Data[0] + "\":{\n\n" + "            \"Region\":\"" + line_Data[1] + "_" + line_Data[2] + "\",\n\n";

                write_Line = write_Line + "            \"Clock model\":\"" + line_Data[3] + "\",\n";

                if (functions.to_Upper_Case(line_Data[3]) == "POISSON")
                {
                    write_Line = write_Line + "            \"Clock model Poisson mean\":\"" + line_Data[4] + "\",\n";
                }
                else if (functions.to_Upper_Case(line_Data[3]) == "NEGATIVE_BINOMIAL")
                {
                    write_Line = write_Line + "            \"Clock model Negative Binomial sucesses\":\"" + line_Data[5] + "\",\n";
                    write_Line = write_Line + "            \"Clock model Negative Binomial probability\":\"" + line_Data[6] + "\",\n";
                }
                else
                {
                    write_Line = write_Line + "            \"Clock model Fixed probability\":\"" + line_Data[7] + "\",\n";
                }

                write_Line = write_Line + "\n";

                int base_pair = 0;
                for (int col = 8; col < line_Data.size(); col++)
                {
                    write_Line = write_Line + "            \"" + base_mutations[base_pair] + "\":\"" + line_Data[col];
                    if (col != line_Data.size() - 1)
                    {
                        write_Line = write_Line + "\",\n";
                    }
                    else
                    {
                        write_Line = write_Line + "\"\n";
                    }
                    base_pair++;
                }
                write_Line = write_Line + "        },\n\n";
                write_Lines.push_back(write_Line);
            }

            cout << "Writing lines\n";

            json_File << "        \"Number of hotspots\":" << write_Lines.size() << ",\n\n";

            for (int line = 0; line < write_Lines.size(); line++)
            {
                json_File << write_Lines[line];
            }

            json_File << "    },\n";

            json_File.close();
        }
        else
        {
            cout << "ERROR: BASE SITE MODEL FILE DOES NOT HAVE THE CORRECT NUMBER OF COLUMNS\n\n";
            exit(-1);
        }

        site_model_File.close();
    }
}

void mutations_T_json::ingress(string convert_Type)
{
    if (convert_Type == "mutations")
    {
        mutations_Convert();
    }
    else
    {
        recombinations_Convert();
    }
}