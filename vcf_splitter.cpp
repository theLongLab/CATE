#include "vcf_splitter.h"

vcf_splitter::vcf_splitter(char folder[], string vcf_Folder, string population_File, string output_Path, int REF, int ALT, int snp_Count)
{
    cout << "Starting up VCF SPLITTER" << endl
         << endl;
    strcpy(this->folder, folder);
    //cout << "Execution folder\t: " << folder << endl;
    this->vcf_Folder =  vcf_Folder;
    cout << "VCF folder path\t: " << this->vcf_Folder << endl;
    this->population_File = population_File;
    cout << "Population file path\t: " << this->population_File << endl;
    this->output_Path = output_Path;
    this->REF = REF;
    this->ALT = ALT;
    cout << endl
         << "Parameters:" << endl;
    cout << "Reference Allele(s)\t: " << REF << endl;
    cout << "Alternate Allele(s)\t: " << ALT << endl;
    this->snp_Count = snp_Count;
    cout << "Maximum number of SNPs per file\t: " << snp_Count << endl;
}

void vcf_splitter::index_population()
{
    cout << endl
         << "Indexing Populations and Patients" << endl;
    fstream pop_File;
    pop_File.open(this->population_File, ios::in);
    if (pop_File.is_open())
    {

        string line;
        //skip heading
        getline(pop_File, line);

        while (getline(pop_File, line))
        {
            vector<string> line_Data;
            split(line_Data, line);
            super_pops.insert(line_Data[5]);
            //cout<<line_Data[0]<<endl<<line_Data[5]<<endl;
            pop_Index.push_back(make_pair(line_Data[0], line_Data[5]));
            //break;
        }
        pop_File.close();
    }

    int pop_Num = super_pops.size();
    cout << pop_Num << " populations found\t: ";
    int count = 0;
    for (auto pop : super_pops)
    {
        cout << pop;
        if (count < pop_Num - 1)
        {
            cout << ", ";
        }

        count++;
    }
    cout << endl;
}

void vcf_splitter::read_File()
{
    list<string> vcf_Files;
    find_VCF(vcf_Folder, vcf_Files);

    for (string file : vcf_Files)
    {
        cout << endl
             << "Processing File\t: " << file << endl
             << endl;

        fstream myFile;
        string file_path = vcf_Folder + "/" + file;

        for (auto pop : super_pops)
        {
            string pop_Folder = output_Path + "/" + pop;

            if (filesystem::exists(pop_Folder) == 0)
            {
                cout << "Creating population folder\t: " << pop_Folder << endl;
                filesystem::create_directory(pop_Folder);
            }

            // fstream pop_File;
            // pop_File.open(pop_Folder + "/" + file.substr(0, file.find_last_of(".")) + "_" + pop + ".vcf", ios::out);
            // pop_File.close();
        }

        myFile.open(file_path, ios::in);
        if (myFile.is_open())
        {

            string line;
            while (getline(myFile, line))
            {
                if (line.substr(0, 2).compare("##") != 0)
                {
                    break;
                }
            }
            //cout << line << endl;
            vector<string> line_Data;
            split(line_Data, line);
            const int size = line_Data.size();
            vector<string> patient_Coutry;
            //cout << line_Data[9];

            for (int c = 9; c < line_Data.size(); c++)
            {
                string country = get_Country(line_Data[c]);
                patient_Coutry.push_back(country);
                //cout << line_Data[c] << "\t" << country<<endl;
            }

            myFile.close();
            write_File_SNP_only(patient_Coutry, file, file_path);
        }
    }
}

void vcf_splitter::write_Header(vector<string> &patient_Coutry, string &file, vector<string> &line_Data_header, string &first_SNP_size)
{
    for (auto pop : super_pops)
    {
        fstream pop_File;
        string pop_Folder = output_Path + "/" + pop;
        pop_File.open(pop_Folder + "/" + file.substr(0, file.find_last_of(".")) + "_" + pop + "_" + first_SNP_size + ".vcf", ios::app);

        for (int c = 0; c < 9; c++)
        {
            pop_File << line_Data_header[c] << "\t";
        }

        for (int c = 9; c < line_Data_header.size(); c++)
        {
            string check = "no";
            if (patient_Coutry[c - 9] == pop)
            {
                //cout<<line_Data[c]<<endl;
                check = "yes";
                pop_File << line_Data_header[c];
            }

            if (check == "yes" && c != (line_Data_header.size() - 1))
            {
                pop_File << "\t";
            }
        }
        pop_File << "\n";
        pop_File.close();
    }
}

void vcf_splitter::write_File_SNP_only(vector<string> &patient_Coutry, string &file, string &file_path)
{

    cout << endl
         << "Splitting VCF and writing to file(s): " << endl;
    fstream myFile;
    myFile.open(file_path, ios::in);

    if (myFile.is_open())
    {

        string line;
        while (getline(myFile, line))
        {
            if (line.substr(0, 2).compare("##") != 0)
            {
                break;
            }
        }

        vector<string> line_Data_header;
        split(line_Data_header, line);

        //here
        //write_Header(patient_Coutry, file, line_Data_header);

        int count = 0;
        string first_SNP_size;
        vector<string> line_Data;
        while (getline(myFile, line))
        {
            //cout<<line;
            string check = split_check(line_Data, line);
            if (check != "no")
            {
                //if count = 0 then write_Header(patient_Coutry,file,line_Data_header) with intial snp size;
                if (count == 0)
                {
                    first_SNP_size = line_Data[1];
                    write_Header(patient_Coutry, file, line_Data_header, first_SNP_size);
                }
                for (auto pop : super_pops)
                {
                    fstream pop_File;
                    string pop_Folder = output_Path + "/" + pop;
                    pop_File.open(pop_Folder + "/" + file.substr(0, file.find_last_of(".")) + "_" + pop + "_" + first_SNP_size + ".vcf", ios::app);

                    for (int c = 0; c < 9; c++)
                    {
                        pop_File << line_Data[c] << "\t";
                    }

                    for (int c = 9; c < line_Data.size(); c++)
                    {
                        string check = "no";
                        if (patient_Coutry[c - 9] == pop)
                        {
                            //cout<<line_Data[c]<<endl;
                            check = "yes";
                            pop_File << line_Data[c];
                        }

                        if (check == "yes" && c != (line_Data.size() - 1))
                        {
                            pop_File << "\t";
                        }
                    }
                    pop_File << "\n";
                    pop_File.close();
                }
                count++;
            }
            //if count = snp_count then rename all pop files with the first and last size value and count = 0;
            if (count == snp_Count)
            {
                string last_SNP_size = line_Data[1];
                for (auto pop : super_pops)
                {
                    string pop_Folder = output_Path + "/" + pop + "/";
                    string old_Name, new_Name;
                    old_Name = pop_Folder + file.substr(0, file.find_last_of(".")) + "_" + pop + "_" + first_SNP_size + ".vcf";
                    new_Name = pop_Folder + file.substr(0, file.find_last_of(".")) + "_" + pop + "_" + first_SNP_size + "_" + last_SNP_size + ".vcf";
                    int check = rename(old_Name.c_str(), new_Name.c_str());
                    if (!check)
                    {
                        cout << "Generated split VCF\t: " << file.substr(0, file.find_last_of(".")) + "_" + pop + "_" + first_SNP_size + "_" + last_SNP_size + ".vcf" << endl;
                    }
                    else
                    {
                        cout << "FAILED TO GENERATE VCF\t: " << file.substr(0, file.find_last_of(".")) + "_" + pop + "_" + first_SNP_size + "_" + last_SNP_size + ".vcf" << endl;
                    }
                }
                count = 0;
            }
        }
        //if count less than 100 then rename all pop files with remaining last snp and finish
        if (count < snp_Count)
        {
            string last_SNP_size = line_Data[1];
            for (auto pop : super_pops)
            {
                string pop_Folder = output_Path + "/" + pop + "/";
                string old_Name, new_Name;
                old_Name = pop_Folder + file.substr(0, file.find_last_of(".")) + "_" + pop + "_" + first_SNP_size + ".vcf";
                new_Name = pop_Folder + file.substr(0, file.find_last_of(".")) + "_" + pop + "_" + first_SNP_size + "_" + last_SNP_size + ".vcf";
                rename(old_Name.c_str(), new_Name.c_str());
            }
        }
        myFile.close();
    }
}

string vcf_splitter::split_check(vector<string> &line_Data, string &line)
{
    string check = "yes";

    vector<string>().swap(line_Data);
    char *convert;
    string capture(line);
    convert = &capture[0];
    //cout<<convert;

    char *split_data;
    split_data = strtok(convert, "\t");
    while (split_data != NULL)
    {
        //cout<<split_data<<endl;
        string char2string;
        char2string.append(split_data);
        //cout << char2string << endl;
        line_Data.push_back(char2string);

        if (line_Data.size() == 5)
        {
            if ((line_Data[3].length()) != REF || line_Data[4].length() != ALT)
            {
                //cout << line_Data[3] << "\t" << line_Data[4] << endl;
                //cout << line_Data[3].length() << "\t" << line_Data[4].length() << endl;
                check = "no";
                break;
            }
        }

        split_data = strtok(NULL, "\t");
    }

    return check;
}

string vcf_splitter::get_Country(string &ID)
{
    string country = "";
    for (auto check : pop_Index)
    {
        if (check.first == ID)
        {
            //cout << check.first << "\t" << check.second << endl;
            country = check.second;
            break;
        }
    }

    return country;
}

void vcf_splitter::split(vector<string> &line_Data, string &line)
{
    vector<string>().swap(line_Data);
    char *convert;
    string capture(line);
    convert = &capture[0];
    //cout<<convert;

    char *split_data;
    split_data = strtok(convert, "\t");

    while (split_data != NULL)
    {
        //cout<<split_data<<endl;
        string char2string;
        char2string.append(split_data);
        //cout << char2string << endl;
        line_Data.push_back(char2string);
        split_data = strtok(NULL, "\t");
    }

    //delete convert;
    //delete split_data;
}

void vcf_splitter::find_VCF(string &vcf_Folder_Path, list<string> &vcf_Files)
{

    for (const auto &entry : filesystem::directory_iterator(vcf_Folder_Path))
    {
        //cout << entry.path().filename().extension() <<endl;
        string extension = entry.path().filename().extension().string();
        if (extension == ".vcf")
        {
            vcf_Files.push_back(entry.path().filename().string());
        }
    }

    //cout<<vcf_Folder_Path<<endl;
}
