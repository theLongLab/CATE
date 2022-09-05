#include "functions.cuh"
#include "gene_extract.h"

gene_extract::gene_extract(string gene_List, string reference_Path, string output_Path, string intermediate_Path)
{
    cout << "Starting up Gene extractor" << endl;
    this->gene_List = gene_List;
    this->reference_Path = reference_Path;
    this->output_Path = output_Path;
    this->intermediate_Path = intermediate_Path;
}

void gene_extract::ingress()
{
    functions function = functions();

    fstream reference;
    reference.open(this->reference_Path, ios::in);

    string full_Reference = "";

    if (reference.is_open())
    {
        cout << endl
             << "Loading reference file: " << this->reference_Path << endl;
        string line;
        // skip header
        getline(reference, line);
        while (getline(reference, line))
        {
            full_Reference.append(line);
        }
        reference.close();
    }

    transform(full_Reference.begin(), full_Reference.end(), full_Reference.begin(), ::toupper);

    cout << "Reference file has been loaded" << endl
         << endl;

    // cout << 10505 << "\t" << full_Reference.at(10505 - 1) << endl;
    // cout << 10506 << "\t" << full_Reference.at(10506 - 1) << endl;

    fstream gene_File;
    gene_File.open(gene_List, ios::in);

    string intermediate_File = intermediate_Path + "/" + filesystem::path(gene_List).stem().string() + ".log_ex";

    cout << "Extracting genes to: " << this->output_Path << endl
         << endl;

    if (gene_File.is_open())
    {
        string gene_Combo;

        if (filesystem::exists(intermediate_File) == 0)
        {
            function.createFile(intermediate_File);
        }
        else
        {
            fstream intermediate;
            intermediate.open(intermediate_File, ios::in);
            string get_finished;
            while (getline(intermediate, get_finished))
            {
                getline(gene_File, gene_Combo);
                if (gene_Combo != get_finished)
                {
                    break;
                }
            }
            intermediate.close();
        }

        fstream intermediate;
        intermediate.open(intermediate_File, ios::app);

        while (getline(gene_File, gene_Combo))
        {
            vector<string> split_Data;
            function.split(split_Data, gene_Combo, '\t');
            string gene_Name = split_Data[0];
            cout << "Gene name\t: " << gene_Name << endl;
            vector<string> coordinates;
            function.split(coordinates, split_Data[1], ':');
            int start_Co = stoi(coordinates[1]);
            int end_Co = stoi(coordinates[2]);
            cout << "Coordinates\t: Chromosome: " << coordinates[0] << " Start: " << start_Co << " End: " << end_Co << endl
                 << endl;

            string output_File = output_Path + "/" + gene_Name + ".fasta";

            if (filesystem::exists(output_File) != 0)
            {
                cout << "Deleting existing file at: " << output_File << endl;
                filesystem::remove(output_File);
            }

            cout << "Writing gene to file: " << output_File << endl;
            function.createFile(output_File, ">" + gene_Name + ":" + split_Data[1]);

            fstream output;
            output.open(output_File, ios::app);

            for (int i = start_Co - 1; i < end_Co; i++)
            {
                output << full_Reference.at(i);
            }

            output.flush();
            intermediate << gene_Combo << "\n";
            output.close();
            cout << endl;
        }
        intermediate.close();
        gene_File.close();
    }
}