#include "fasta_merge.h"
#include "functions.cuh"

fasta_merge::fasta_merge(string fasta_Folder, string output_FASTA)
{
    cout << "Starting up FASTA MERGER" << endl
         << endl;
    this->fasta_Folder = fasta_Folder;
    cout << "Raw FASTA folder: " << this->fasta_Folder << endl;
    this->output_FASTA = output_FASTA;
}

void fasta_merge::ingress()
{
    functions function = functions();
    vector<string> files = index_FASTA_folder();

    cout << endl
         << "Intitating MERGE:" << endl;

    function.createFile(output_FASTA);
    fstream FASTA_file;
    FASTA_file.open(this->output_FASTA, ios::app);

    for (string file : files)
    {
        fstream raw_File;
        raw_File.open(file, ios::in);
        if (raw_File.is_open())
        {
            cout << "Merging file: " << file << endl;
            string line;
            while (getline(raw_File, line))
            {
                FASTA_file << line << "\n";
            }
            raw_File.close();
        }
    }

    FASTA_file.close();
    cout << endl
         << "All files MERGED to \"" << output_FASTA << "\" Fasta File" << endl;
}

vector<string> fasta_merge::index_FASTA_folder()
{
    cout << "Indexing FASTA Folder" << endl;
    vector<string> files;

    for (const auto &entry : filesystem::directory_iterator(this->fasta_Folder))
    {
        string file = entry.path().string();
        string ext = file.substr(file.find_last_of("."), file.length());
        transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext == ".fasta" || ext == ".fna" || ext == ".ffn" || ext == ".faa" || ext == ".frn" || ext == "fa")
        {
            files.push_back(file);
        }
    }

    cout << files.size() << " FASTA files found" << endl;

    return files;
}