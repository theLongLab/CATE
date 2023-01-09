#include "functions.cuh"
#include "map2gene.cuh"

map2gene::map2gene(string input_map_File_path, string ouput_Path, string pre_fix)
{
    cout << "Starting up MAP to Gene list" << endl
         << endl;

    this->input_map_File_path = input_map_File_path;
    this->ouput_Path = ouput_Path;

    this->pre_fix = pre_fix;
}

void map2gene::ingress()
{
    functions function = functions();

    fstream map_File;
    map_File.open(this->input_map_File_path, ios::in);

    string line;
    vector<string> line_Data;

    if (map_File.is_open())
    {
        cout << "Processing map file: " << this->input_map_File_path << endl
             << endl;

        // int pre_index = -1;
        int index;
        string pre_CHR = "";
        fstream chr_gene_File;

        while (getline(map_File, line))
        {
            function.split(line_Data, line, '\t');

            if (pre_CHR != line_Data[0])
            {
                cout << "Writing chromosome: " << line_Data[0] << endl;
                index = chr_Index(line_Data[0]);
                pre_CHR = line_Data[0];

                if (chr_gene_File.is_open())
                {
                    chr_gene_File.close();
                }
                chr_gene_File.open(this->ouput_Path + "/" + filesystem::path(input_map_File_path).stem().string() + "_" + line_Data[0] + ".txt", ios::app);
            }
            // pre_index = index;

            chr_gene_File << this->pre_fix + "_" + to_string(count[index]) + "\t" + line_Data[1] << "\n";
            count[index] = count[index] + 1;
        }
        chr_gene_File.close();
        map_File.close();
    }
}

int map2gene::chr_Index(string chr)
{
    functions function = functions();
    int index = -1;

    for (size_t i = 0; i < chromosomes.size(); i++)
    {
        if (chromosomes[i] == chr)
        {
            index = i;
            break;
        }
    }

    if (index == -1)
    {
        chromosomes.push_back(chr);
        count.push_back(1);
        index = chromosomes.size() - 1;

        function.createFile(this->ouput_Path + "/" + filesystem::path(input_map_File_path).stem().string() + "_" + chr + ".txt");
    }

    return index;
}