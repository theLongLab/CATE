#include "functions.cuh"
#include "gene_extract.h"

gene_extract::gene_extract(string gene_List, string reference_Path, string output_Path, string intermediate_Path)
{
    /**
     * * Constructor Function
     * Assigns passed variables to the classes' private variable.
     **/

    cout << "Starting up Gene extractor" << endl;
    this->gene_List = gene_List;
    this->reference_Path = reference_Path;
    this->output_Path = output_Path;
    this->intermediate_Path = intermediate_Path;
}

void gene_extract::ingress()
{
    /**
     * Execution function.
     **/

    /**
     * Call the "functions" class. Bespoke functions commonly used by CATE.
     **/
    functions function = functions();

    /**
     * Load the reference file to the RAM.
     * This is a double edged sword.
     * Immensely increases the speed of finding and writing the genes.
     * But,
     * Requires RAM to be large at least to match the size of the reference sequence.
     * Does take time to load the reference sequence to the RAM.
     **/
    fstream reference;
    reference.open(this->reference_Path, ios::in);

    string full_Reference = "";

    if (reference.is_open())
    {
        cout << endl
             << "Loading reference file: " << this->reference_Path << endl;

        /**
         * @param line captures the sequences line read from the FASTA file.
         **/
        string line;

        /**
         * Skip the first line of the FASTA sequence.
         * This is the ID of the sequence (begins with a ">" symbol).
         **/
        getline(reference, line);

        while (getline(reference, line))
        {
            full_Reference.append(line);
        }
        reference.close();
    }

    /**
     * REFERENCE sequence will be converted to UPPERCASE.
     * Helps create uniformity.
     **/
    transform(full_Reference.begin(), full_Reference.end(), full_Reference.begin(), ::toupper);

    cout << "Reference file has been loaded" << endl
         << endl;

    // cout << 10505 << "\t" << full_Reference.at(10505 - 1) << endl;
    // cout << 10506 << "\t" << full_Reference.at(10506 - 1) << endl;

    /**
     * Initiate the reading of the gene file.
     **/
    fstream gene_File;
    gene_File.open(gene_List, ios::in);

    /**
     * Log file created in the intermediate folder.
     * ! This helps with the resume function. Automatically resumes from the last completely processed gene in the event of a program crash.
     **/
    string intermediate_File = intermediate_Path + "/" + filesystem::path(gene_List).stem().string() + ".log_ex";

    cout << "Extracting genes to: " << this->output_Path << endl
         << endl;

    if (gene_File.is_open())
    {
        /**
         * @param gene_Combo used to capture and extract info of each gene combination.
         **/
        string gene_Combo;

        /**
         * If the intermediate log file is absent this run will be considered as a brand new run of this query and,
         * the intermediate log file will be created.
         **/
        if (filesystem::exists(intermediate_File) == 0)
        {
            function.createFile(intermediate_File);
        }
        else
        {
            /**
             * If the intermediate log file present then the resume process will initiated.
             * This is a unintelligent resume. Essentially it matches the each read line written with the lines read from the gene file.
             * The break will occur as soon as their is a mismatch.
             * To counter any errors it is advised to have a new gene file name or a new intermediate folder per new run.
             **/
            fstream intermediate;
            intermediate.open(intermediate_File, ios::in);

            /**
             * @param get_finished comparison variable. Used o compare the intermediate file data with that of the gene file.
             **/
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
            /**
             * @param split_Data vector captures split function's outputs on the genes information.
             **/
            vector<string> split_Data;
            function.split(split_Data, gene_Combo, '\t');

            /**
             * @param gene_Name captures the gene's name.
             **/
            string gene_Name = split_Data[0];
            cout << "Gene name\t: " << gene_Name << endl;

            /**
             * @param coordinates vector captures split function's outputs on gene coordinates.
             * [0] = chromosome
             * [1] = start position
             * [2] = end position
             **/
            vector<string> coordinates;
            function.split(coordinates, split_Data[1], ':');

            /**
             * @param start_Co captures query gene's start position as an integer.
             **/
            int start_Co = stoi(coordinates[1]);
            /**
             * @param end_Co captures query gene's end position as an integer.
             **/
            int end_Co = stoi(coordinates[2]);
            cout << "Coordinates\t: Chromosome: " << coordinates[0] << " Start: " << start_Co << " End: " << end_Co << endl
                 << endl;

            /**
             * @param output_File defines the output file for the query gene sequence.
             **/
            string output_File = output_Path + "/" + gene_Name + ".fasta";

            /**
             * If the file already exists it will be deleted without prejudice.
             * So be advised about this.
             **/
            if (filesystem::exists(output_File) != 0)
            {
                cout << "Deleting existing file at: " << output_File << endl;
                filesystem::remove(output_File);
            }

            /**
             * Writing of the sequence to the FILE. The sequence ID will be a combination of the gene's name and coordinates.
             **/
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