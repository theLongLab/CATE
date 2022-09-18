#include "gff2gene.cuh"
#include "functions.cuh"

gff2gene::gff2gene(string input_File, string output_Path)
{
    /**
     * * Constructor Function
     * Assigns passed variables to the classes' private variable.
     **/

    cout << "Starting up GFF to Gene list" << endl
         << endl;
    this->input_File = input_File;
    this->output_Path = output_Path;
}

void gff2gene::ingress()
{
    /**
     * Execution function.
     **/

    /**
     * Call the "functions" class. Bespoke functions commonly used by CATE.
     **/
    functions function = functions();

    cout << "Processing GFF file: " << this->input_File << endl
         << endl;

    /**
     * Open the GFF file from the location.
     **/
    fstream gff_File;
    gff_File.open(input_File, ios::in);

    if (gff_File.is_open())
    {
        /**
         * @param line captures the sequences line read from the FASTA file.
         **/
        string line;

        while (getline(gff_File, line))
        {
            /**
             * Skip the header lines beginning with a "#" or "##" at the start of a GFF file.
             **/
            if (line.at(0) != '#')
            {
                break;
            }
        }

        /**
         * @param file_Name defines the output gene file, it will be a derivation of the GFF file's name.
         **/
        string file_Name = this->output_Path + "/" + filesystem::path(input_File).stem().string() + ".txt";
        fstream gene_List;
        gene_List.open(file_Name, ios::out);

        while (getline(gff_File, line))
        {
            // cout << line << endl;

            /**
             * @param split_Data vector captures split function's outputs.
             * [0] = sequence ID
             * [1] = source
             * [2] = feature type
             * [3] = start
             * [4] = end
             * [5] = score
             * [6] = strand
             * [7] = phase
             * [8] = atributes
             **/
            vector<string> split_Data;
            function.split(split_Data, line, '\t');

            /**
             * Ensures that it is a target data row. They all have 9 columns as standard.
             **/
            if (split_Data.size() == 9)
            {
                // cout << line << endl;

                /**
                 * @param feature captures the feature type which is the third column.
                 **/
                string feature = split_Data[2];
                /**
                 * Capitalizes the string to normalize the data and prevent mismatch errors.
                 **/
                transform(feature.begin(), feature.end(), feature.begin(), ::toupper);

                /**
                 * If the feature is equal to being a GENE it's data will be extracted.
                 **/
                if (feature == "GENE")
                {
                    /**
                     * @param attribute Get attribute information to get the gene's NAME if available else it will be "NA".
                     **/
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
                        /**
                         * Gene names are followed by the "ID" tag.
                         **/
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