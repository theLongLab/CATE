#include "fasta_splitter.h"
#include "functions.cuh"

fasta_splitter::fasta_splitter(string fasta_File, string output_Path, string fasta_Name)
{
    /**
     * * Constructor Function
     * Assigns passed variables to the classes' private variable.
     **/

    cout << "Starting up FASTA SPLITTER" << endl
         << endl;
    this->fasta_File = fasta_File;
    cout << "Raw FASTA file\t : " << this->fasta_File << endl;
    this->output_Folder = output_Path;

    /**
     * First we convert the FASTA ID to lowercase to prevent user error if they have selected the option to extract all the
     * sequences from the merged FASTA file.
     **/

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

    /**
     * Execution function.
     **/

    /**
     * Ensures first and foremost that the entered file does exist in the said location.
     **/

    if (filesystem::exists(this->fasta_File) == 0)
    {

        /**
         *  ! Initialised if the merged FASTA file does not exists in the said location.
         **/

        cout << "INVALID FASTA FILE." << endl
             << "FILE \"" << this->fasta_File << "\" WAS NOT FOUND AT THE LOCATION." << endl;
    }
    else
    {
        /**
         * Executes respective function based on the split requirement of the user.
         **/
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
    /**
     * This function is used to extract a specific sequence from the merged file.
     **/

    /**
     * Call the "functions" class. Bespoke functions commonly used by CATE.
     **/
    functions function = functions();

    fstream fasta_File;

    /**
     * @param split_Fasta defines the FASTA file names of the output files.
     **/
    string split_Fasta;

    /**
     * @param check defines the sequence ID that needs to be extracted.
     **/
    string check = ">" + fasta_Name;

    fasta_File.open(this->fasta_File, ios::in);

    if (fasta_File.is_open())
    {

        /**
         * @param found acts as a boolean variable.
         * Once the query sequence is found it is used to trigger the writing of the sequence to the output file.
         **/
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
                /**
                 * If found is 1 and another line starts with ">" it means that the query sequence is complete.
                 * Therefore the file read loop will be broken and the program will end
                 **/
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
                /**
                 * All lines will be written till the next sequence.
                 **/
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
    /**
     * This function is used to extract all the sequences from the merged file.
     **/

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

        /**
         * @param first acts as a boolean variable.
         * Ensures close is not triggered on the first sequence itself.
         **/
        int first = 0;

        while (getline(fasta_File, line))
        {
            if (line.at(0) == '>')
            {
                if (first == 1)
                {
                    /**
                     * Ensures close is not triggered on the first sequence itself.
                     * Ensures currently written sequence is completed and the next sequence has begun.
                     * The current sequence's file will be closed.
                     **/
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
                // fstream output;
                // output.open(split_Fasta, ios::app);
                output << line << "\n";
                // output.close();
            }
        }
        output.close();
        fasta_File.close();
    }
}