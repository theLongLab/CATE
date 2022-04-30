#include "fasta_splitter.h"
#include "functions.cuh"

fasta_splitter::fasta_splitter(string fasta_File, string output_Path, string fasta_Name)
{
    cout << "Starting up FASTA SPLITTER" << endl
         << endl;
    this->fasta_File = fasta_File;
    cout << "Raw FASTA file\t : " << this->fasta_File << endl;
    this->output_Folder = output_Path;

    string check = fasta_Name;
    transform(check.begin(), check.end(), check.begin(), ::tolower);
    if (check == "all")
    {
        this->fasta_Name = check;
    }
    else
    {
        if (fasta_Name.at(0) == '>')
        {
            this->fasta_Name = fasta_Name.substr(1, fasta_Name.length());
        }
        else
        {
            this->fasta_Name = fasta_Name;
        }
    }
    cout << endl;
}

void fasta_splitter::ingress()
{

    if (filesystem::exists(this->fasta_File) == 0)
    {
        cout << "INVALID FASTA FILE." << endl
             << "FILE \"" << this->fasta_File << "\" WAS NOT FOUND AT THE LOCATION." << endl;
    }
    else
    {
        if (fasta_Name == "all")
        {
            split_all();
        }
        else
        {
            split_select();
        }
    }
}

void fasta_splitter::split_select()
{

    functions function = functions();
    fstream fasta_File;
    string split_Fasta;
    string check = ">" + fasta_Name;

    fasta_File.open(this->fasta_File, ios::in);

    if (fasta_File.is_open())
    {
        int found = 0;
        cout << "Reading FASTA file and extracting sequence: " << fasta_Name << endl
             << endl;
        split_Fasta = this->output_Folder + "/" + fasta_Name + ".fasta";
        function.createFile(split_Fasta, check);
        fstream output;
        output.open(split_Fasta, ios::app);

        string line;
        while (getline(fasta_File, line))
        {
            if (line.at(0) == '>')
            {
                if (found == 1)
                {
                    break;
                }

                if (line == check)
                {
                    found = 1;
                    cout << "WRITING sequence: " << line.substr(1, line.length()) << endl;
                }
                else
                {
                    cout << "Skipping sequence: " << line.substr(1, line.length()) << endl;
                }
            }
            else
            {
                if (found == 1)
                {
                    output << line << "\n";
                }
            }
        }
        output.close();
        fasta_File.close();
    }
}

void fasta_splitter::split_all()
{
    functions function = functions();
    fstream fasta_File;
    string split_Fasta;

    fasta_File.open(this->fasta_File, ios::in);

    if (fasta_File.is_open())
    {
        cout << "Reading and splitting entire FASTA file" << endl
             << endl;

        string line;
        fstream output;
        int first = 0;
        while (getline(fasta_File, line))
        {
            if (line.at(0) == '>')
            {
                if (first == 1)
                {
                    output.close();
                }
                first = 1;
                cout << "Writing sequence: " << line.substr(1, line.length()) << endl;
                split_Fasta = this->output_Folder + "/" + line.substr(1, line.length()) + ".fasta";
                function.createFile(split_Fasta, line);
                output.open(split_Fasta, ios::app);
            }
            else
            {
                //fstream output;
                //output.open(split_Fasta, ios::app);
                output << line << "\n";
                //output.close();
            }
        }
        output.close();
        fasta_File.close();
    }
}