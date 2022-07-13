#include "gff2gene.cuh"
#include "functions.cuh"

gff2gene::gff2gene(string input_File, string output_Path)
{
    cout << "Starting up GFF to Gene list" << endl
         << endl;
    this->input_File = input_File;
    this->output_Path = output_Path;
}

void gff2gene::ingress()
{
    functions function = functions();

    cout << "Processing GFF file: " << this->input_File << endl
         << endl;

    fstream gff_File;
    gff_File.open(input_File, ios::in);

    if (gff_File.is_open())
    {
        string line;

        while (getline(gff_File, line))
        {
            if (line.at(0) != '#')
            {
                break;
            }
        }

        string file_Name = this->output_Path + "/" + filesystem::path(input_File).stem().string() + ".txt";
        fstream gene_List;
        gene_List.open(file_Name, ios::out);

        while (getline(gff_File, line))
        {
            // cout << line << endl;

            vector<string> split_Data;
            function.split(split_Data, line, '\t');

            if (split_Data.size() == 9)
            {
               // cout << line << endl;
                string feature = split_Data[2];
                transform(feature.begin(), feature.end(), feature.begin(), ::toupper);

                if (feature == "GENE")
                {
                    string attribute = split_Data[8];
                    vector<string> split_attributes;
                    function.split(split_attributes, attribute, ';');

                    string gene_Name = "NA";

                    for (string split_attribute : split_attributes)
                    {
                        // string ID_check;
                        vector<string> ID_check_split;
                        function.split(ID_check_split, split_attribute, '=');

                        transform(ID_check_split[0].begin(), ID_check_split[0].end(), ID_check_split[0].begin(), ::toupper);
                        if (ID_check_split[0] == "ID")
                        {
                            gene_Name = ID_check_split[1];
                            break;
                        }
                    }

                    string write_Line = gene_Name + "\t" + split_Data[0] + ":" + split_Data[3] + ":" + split_Data[4] + "\n";
                    cout << write_Line;
                    gene_List << write_Line;
                }

                // REMOVE AFTER TESTING
               // break;
            }
        }

        gff_File.close();
        gene_List.close();
    }
}