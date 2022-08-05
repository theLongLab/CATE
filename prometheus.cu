#include "prometheus.cuh"
#include "functions.cuh"

prometheus::prometheus(vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, long int combinations, float an, float e1, float e2, int N, int CPU_cores, int SNPs_per_Run, int number_of_genes)
{
    // TAJIMA CONSTRUCTOR

    set_Values(folder_Index, tot_Blocks, tot_ThreadsperBlock, Multi_read, CPU_cores, SNPs_per_Run, number_of_genes);

    this->an = an;
    this->e1 = e1;
    this->e2 = e2;
    this->N = N;
    this->combinations = combinations;

    // this->gene_Size = gene_Size;
}

prometheus::prometheus(string output_File, int window_Size, int step_Size, vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, long int combinations, float an, float e1, float e2, int N, int CPU_cores, int SNPs_per_Run, int number_of_genes)
{
    // TAJIMA WINDOW CONSTRUCTOR

    set_Values(folder_Index, tot_Blocks, tot_ThreadsperBlock, Multi_read, CPU_cores, SNPs_per_Run, number_of_genes);
    set_Values_Window(output_File, window_Size, step_Size, number_of_genes);

    this->an = an;
    this->e1 = e1;
    this->e2 = e2;
    this->N = N;
    this->combinations = combinations;
}

prometheus::prometheus(vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star)
{
    // FU LI Constructor

    set_Values(folder_Index, tot_Blocks, tot_ThreadsperBlock, Multi_read, CPU_cores, SNPs_per_Run, number_of_genes);

    this->N = N;
    this->N_float = (float)N;
    this->combinations = combinations;

    this->an = an;
    this->vd = vd;
    this->ud = ud;
    this->vd_star = vd_star;
    this->ud_star = ud_star;
    this->uf = uf;
    this->vf = vf;
    this->uf_star = uf_star;
    this->vf_star = vf_star;
}

prometheus::prometheus(string output_File, int window_Size, int step_Size, vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star)
{
    // FU LI CONSTRUCTOR WINDOW

    set_Values(folder_Index, tot_Blocks, tot_ThreadsperBlock, Multi_read, CPU_cores, SNPs_per_Run, number_of_genes);
    set_Values_Window(output_File, window_Size, step_Size, number_of_genes);

    this->N = N;
    this->N_float = (float)N;
    this->combinations = combinations;

    this->an = an;
    this->vd = vd;
    this->ud = ud;
    this->vd_star = vd_star;
    this->ud_star = ud_star;
    this->uf = uf;
    this->vf = vf;
    this->uf_star = uf_star;
    this->vf_star = vf_star;
}

prometheus::prometheus(vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float bn, float bn_plus1)
{
    // FAY WU Constructor

    set_Values(folder_Index, tot_Blocks, tot_ThreadsperBlock, Multi_read, CPU_cores, SNPs_per_Run, number_of_genes);

    this->N = N;
    this->N_float = (float)N;
    this->combinations = combinations;

    this->an = an;
    this->bn = bn;
    this->bn_plus1 = bn_plus1;
}

prometheus::prometheus(string output_File, int window_Size, int step_Size, vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float bn, float bn_plus1)
{
    // FAY WU Constructor WINDOW

    set_Values(folder_Index, tot_Blocks, tot_ThreadsperBlock, Multi_read, CPU_cores, SNPs_per_Run, number_of_genes);
    set_Values_Window(output_File, window_Size, step_Size, number_of_genes);

    this->N = N;
    this->N_float = (float)N;
    this->combinations = combinations;

    this->an = an;
    this->bn = bn;
    this->bn_plus1 = bn_plus1;
}

prometheus::prometheus(vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float e1, float e2, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star, float bn, float bn_plus1)
{
    // NEUTRALITY
    set_Values(folder_Index, tot_Blocks, tot_ThreadsperBlock, Multi_read, CPU_cores, SNPs_per_Run, number_of_genes);

    this->N = N;
    this->N_float = (float)N;
    this->combinations = combinations;

    this->an = an;
    this->e1 = e1;
    this->e2 = e2;

    this->an = an;
    this->vd = vd;
    this->ud = ud;
    this->vd_star = vd_star;
    this->ud_star = ud_star;
    this->uf = uf;
    this->vf = vf;
    this->uf_star = uf_star;
    this->vf_star = vf_star;

    this->bn = bn;
    this->bn_plus1 = bn_plus1;
}

prometheus::prometheus(string output_File, int window_Size, int step_Size, vector<pair<string, string>> folder_Index, string Multi_read, int tot_Blocks, int tot_ThreadsperBlock, int CPU_cores, int SNPs_per_Run, int number_of_genes, int N, long int combinations, float an, float e1, float e2, float vd, float ud, float vd_star, float ud_star, float uf, float vf, float uf_star, float vf_star, float bn, float bn_plus1)
{
    // NEUTRALITY CONSTRUCTOR WINDOW

    set_Values(folder_Index, tot_Blocks, tot_ThreadsperBlock, Multi_read, CPU_cores, SNPs_per_Run, number_of_genes);
    set_Values_Window(output_File, window_Size, step_Size, number_of_genes);

    this->N = N;
    this->N_float = (float)N;
    this->combinations = combinations;

    this->an = an;
    this->e1 = e1;
    this->e2 = e2;

    this->an = an;
    this->vd = vd;
    this->ud = ud;
    this->vd_star = vd_star;
    this->ud_star = ud_star;
    this->uf = uf;
    this->vf = vf;
    this->uf_star = uf_star;
    this->vf_star = vf_star;

    this->bn = bn;
    this->bn_plus1 = bn_plus1;
}

void prometheus::set_Values(vector<pair<string, string>> folder_Index, int tot_Blocks, int tot_ThreadsperBlock, string Multi_read, int CPU_cores, int SNPs_per_Run, int number_of_genes)
{
    this->folder_Index = folder_Index;

    this->tot_Blocks = tot_Blocks;
    this->tot_ThreadsperBlock = tot_ThreadsperBlock;

    this->Multi_read = Multi_read;
    this->CPU_cores = CPU_cores - 1;
    this->SNPs_per_Run = SNPs_per_Run;

    cout << "Processing on " << this->CPU_cores + 1 << " CPU cores" << endl;
    cout << "Processing " << number_of_genes << " genes at a time" << endl;
    cout << "Processing " << this->SNPs_per_Run << " SNPs at a time" << endl;
    if (this->Multi_read == "YES")
    {
        cout << "Multi read: Available" << endl;
    }
    else
    {
        cout << "Multi read: Unavailable" << endl;
    }
    cout << endl;
}

void prometheus::set_Values_Window(string output_File, int window_Size, int step_Size, int number_of_genes)
{
    this->output_File = output_File;

    this->window_Size = window_Size;
    this->step_Size = step_Size;
    this->number_of_genes_Window = number_of_genes;

    this->calc_Mode = "WINDOW";

    cout << "Window size: " << window_Size << endl;
    if (step_Size != 0)
    {
        cout << "Step size: " << step_Size << endl;
    }
    else
    {
        cout << "Sliding Window" << endl;
    }
    cout << endl;
}

void prometheus::process_C_sliding_Window(string test)
{
    functions function = functions();
    this->sliding_Mode = "YES";

    int pre_Lock = 1;
    int file_Count_Start = 0;
    int line_Num = 0;

    if (filesystem::exists(output_File) != 0)
    {
        int found = 0;

        fstream output_Check;
        output_Check.open(output_File, ios::in);
        if (output_Check.is_open())
        {
            string line_Check;
            getline(output_Check, line_Check); // skip first header line

            for (int file_Count = 0; file_Count < folder_Index.size(); file_Count++)
            {
                string file_Path = folder_Index[file_Count].second;
                fstream file;
                file.open(file_Path, ios::in);
                int line_Current = 0;

                if (file.is_open())
                {
                    string line;
                    getline(file, line); // skip first header line
                    while (getline(file, line))
                    {
                        line_Current++;
                        int VALID = function.get_Valid(line);
                        if (VALID != -1)
                        {
                            getline(output_Check, line_Check);
                            string trim = line_Check.substr(0, line_Check.find('\t'));

                            vector<string> positions;
                            function.split_getPos_ONLY(positions, line, '\t');
                            string pos = positions[1] + ":" + to_string((stoi(positions[1]) + window_Size));

                            if (pos != trim)
                            {
                                found = 1;
                                file_Count_Start = file_Count;
                                line_Num = line_Current;
                                break;
                            }
                        }
                    }
                    file.close();
                }
                if (found == 1)
                {
                    break;
                }
            }
            output_Check.close();
        }
    }
    else
    {
        if (test == "T")
        {
            function.createFile(output_File, "Coordinates\tPi\tS\tTajimas_D");
        }
        else if (test == "FU")
        {
            function.createFile(output_File, "Coordinates\tPi\tS\tne\tns\tD\tD_star\tF\tF_star");
        }
        else if (test == "FA")
        {
            function.createFile(output_File, "Coordinates\tPi\tS\tTotal_iEi\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E");
        }
        else if (test == "N")
        {
            function.createFile(output_File, "Coordinates\tPi\tS\tne\tns\tTotal_iEi\tTajimas_D\tD\tD_star\tF\tF_star\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E");
        }
    }

    vector<int> row_vec;

    fstream output;
    output.open(output_File, ios::app);

    cout << "Writing to file\t: " << output_File << endl;
    cout << endl;

    int line_Current = 0;

    for (int file_Count = file_Count_Start; file_Count < folder_Index.size(); file_Count++)
    {
        string file_Path = folder_Index[file_Count].second;
        fstream file;
        file.open(file_Path, ios::in);

        if (file.is_open())
        {
            string line;
            getline(file, line); // skip first header line
            while (getline(file, line))
            {
                if (line_Current < line_Num)
                {
                    line_Current++;
                }
                else
                {
                    vector<string> positions;
                    function.split_getPos_ONLY(positions, line, '\t');
                    int start_Co = stoi(positions[1]);
                    int end_Co = start_Co + window_Size;

                    // vector<string> file_List = function.forward_Search_Only(file_Count, this->folder_Index, start_Co, end_Co);
                    all_start_Co.push_back(start_Co);
                    all_end_Co.push_back(end_Co);

                    string write = to_string(start_Co) + ":" + to_string(end_Co);
                    write_Lines.push_back(write);

                    seg_catch_points_ALL.push_back(-1);
                    seg_catch_index_ALL.push_back(-1);
                    seg_backward_index_ALL.push_back(row_vec);
                    seg_forward_index_ALL.push_back(row_vec);

                    catch_Point.push_back(1);

                    if (write_Lines.size() == number_of_genes_Window)
                    {
                        // Process
                        // intialize(coordinates);
                        int folder_Start = all_start_Co[0];
                        int folder_End = all_end_Co[all_end_Co.size() - 1];
                        this->gene_Size = write_Lines.size();

                        cout << "Processing windows from " << folder_Start << " to " << folder_End << endl;

                        // Collect FILEs
                        vector<string> file_List;
                        if (folder_Index.size() > 1)
                        {
                            file_List = function.compound_interpolationSearch(folder_Index, folder_Start, folder_End);
                        }
                        else
                        {
                            file_List.push_back(folder_Index[0].second);
                        }

                        if (this->prev_file_List != file_List)
                        {
                            if (pre_Lock == 1)
                            {
                                pre_Lock = 0;
                            }
                            else
                            {
                                if (test == "T")
                                {
                                    free(pre_MA);
                                }
                                else if (test == "FU")
                                {
                                    free(pre_MA);
                                    free(pre_ne);
                                    free(pre_ns);
                                }
                                else if (test == "FA")
                                {
                                    free(pre_MA);
                                    free(pre_Theta_partials);
                                }
                                else if (test == "N")
                                {
                                    free(pre_MA);

                                    free(pre_Theta_partials);

                                    free(pre_ne);
                                    free(pre_ns);
                                }

                                pre_Lock = 1;
                            }

                            same_Files = "NO";

                            tot_Segs = 0;

                            all_Lines.clear();
                            position_index_Segs.clear();

                            concat_Segs.clear();
                            all_site_Index.clear();

                            start_stop.clear();

                            if (Multi_read == "YES")
                            {
                                cout << "Intitating multi read based segregating site search" << endl;

                                vector<thread> threads_Multi_read;

                                for (string files : file_List)
                                {
                                    // cout << files << endl;
                                    threads_Multi_read.push_back(thread{&prometheus::file_Reader_multi, this, files});
                                }

                                for (thread &t : threads_Multi_read)
                                {
                                    if (t.joinable())
                                    {
                                        t.join();
                                    }
                                }

                                threads_Multi_read.clear();
                            }
                            else
                            {
                                cout << "Intitating single read based segregating site search" << endl;
                                for (string files : file_List)
                                {
                                    fstream file;
                                    file.open(files, ios::in);

                                    if (file.is_open())
                                    {
                                        string line;
                                        getline(file, line); // skip first header line
                                        while (getline(file, line))
                                        {
                                            all_Lines.push_back(line);
                                        }
                                        file.close();
                                    }
                                }
                            }

                            // Process Segs
                            tot_Segs = all_Lines.size();
                            cout << "System is processing and filtering " << tot_Segs << " segregating site(s)" << endl;
                            // all_Files_index.clear();

                            int GPU_rounds_full = tot_Segs / SNPs_per_Run;
                            int GPU_rounds_partial = tot_Segs % SNPs_per_Run;

                            // vector<pair<int, int>> start_stop;
                            for (int i = 0; i < GPU_rounds_full; i++)
                            {
                                int start = i * SNPs_per_Run;
                                int stop = start + SNPs_per_Run;
                                start_stop.push_back(make_pair(start, stop));
                            }

                            if (GPU_rounds_partial != 0)
                            {
                                int start = tot_Segs - GPU_rounds_partial;
                                int stop = tot_Segs;
                                start_stop.push_back(make_pair(start, stop));
                            }

                            vector<thread> threads_vec;

                            for (int rounds = 0; rounds < start_stop.size(); rounds++)
                            {
                                concat_Segs.push_back("");
                                int *row;
                                all_site_Index.push_back(row);
                            }

                            for (int rounds = 0; rounds < start_stop.size(); rounds++)
                            {
                                threads_vec.push_back(thread{&prometheus::seg_Concat, this, rounds, start_stop[rounds].first, start_stop[rounds].second});
                            }

                            for (thread &t : threads_vec)
                            {
                                if (t.joinable())
                                {
                                    t.join();
                                }
                            }

                            threads_vec.clear();

                            prev_file_List.clear();
                            this->prev_file_List = file_List;
                        }
                        else
                        {
                            same_Files = "YES";
                        }

                        // Process test
                        if (test == "T")
                        {
                            process_Tajima();
                            cout << "System has completed Tajima's D for the gene(s)" << endl;
                        }
                        else if (test == "FU")
                        {
                            process_Fu_Li();
                            cout << "System has completed Fu and Li for the gene(s)" << endl;
                        }
                        else if (test == "FA")
                        {
                            process_Fay_Wu();
                            cout << "System has completed Fay and Wu for the gene(s)" << endl;
                        }
                        else if (test == "N")
                        {
                            process_Neutrality();
                            cout << "System has completed Neutrality tests for the gene(s)" << endl;
                        }

                        for (int gene_Count = 0; gene_Count < write_Lines.size(); gene_Count++)
                        {
                            if (seg_catch_points_ALL[gene_Count] != -1)
                            {
                                output << write_Lines[gene_Count] << "\n";
                            }
                        }

                        // for (string line : write_Lines)
                        // {
                        //     output << line << "\n";
                        // }

                        output.flush();

                        cout << "System has written the results for the batch" << endl
                             << endl;

                        // coordinates.clear();
                        all_start_Co.clear();
                        all_end_Co.clear();
                        write_Lines.clear();

                        catch_Point.clear();

                        seg_catch_points_ALL.clear();
                        seg_catch_index_ALL.clear();
                        seg_backward_index_ALL.clear();
                        seg_forward_index_ALL.clear();

                        gene_Size = 0;
                    }
                }
            }

            file.close();
        }
        // REMOVE AFTER testing
        break;
    }

    if (write_Lines.size() != 0)
    {
    }
}

void prometheus::process_Window(string test)
{
    functions function = functions();

    // CREATE the output files
    // OUTPUT FILE columns changes based on test type

    int pre_Lock = 1;

    int start_Value = stoi(folder_Index[0].first.substr(0, folder_Index[0].first.find('_')));
    int end_Value = stoi(folder_Index[folder_Index.size() - 1].first.substr(folder_Index[folder_Index.size() - 1].first.find('_') + 1));

    // cout << start_Value << endl;
    // cout << end_Value << endl;

    int start_Co = 0;
    int end_Co = start_Co + window_Size;

    while (start_Value > end_Co)
    {
        start_Co = start_Co + step_Size;
        end_Co = start_Co + window_Size;
    }

    // vector<pair<int, int>> coordinates;

    // RESUME FUNCTION
    int caught = 0;

    if (filesystem::exists(output_File) != 0)
    {
        // skip
        fstream output;
        output.open(output_File, ios::in);

        string output_Line;

        while (start_Co <= end_Value)
        {

            // skip header
            getline(output, output_Line);

            while (getline(output, output_Line))
            {
                string trim = output_Line.substr(0, output_Line.find('\t'));
                string check = to_string(start_Co) + ":" + to_string(end_Co);
                if (trim != check)
                {
                    caught = 1;
                    break;
                }
            }

            if (caught == 1)
            {
                break;
            }

            start_Co = start_Co + step_Size;
            end_Co = start_Co + window_Size;
        }
        output.close();
    }
    else
    {
        if (test == "T")
        {
            function.createFile(output_File, "Coordinates\tPi\tS\tTajimas_D");
        }
        else if (test == "FU")
        {
            function.createFile(output_File, "Coordinates\tPi\tS\tne\tns\tD\tD_star\tF\tF_star");
        }
        else if (test == "FA")
        {
            function.createFile(output_File, "Coordinates\tPi\tS\tTotal_iEi\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E");
        }
        else if (test == "N")
        {
            function.createFile(output_File, "Coordinates\tPi\tS\tne\tns\tTotal_iEi\tTajimas_D\tD\tD_star\tF\tF_star\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E");
        }
    }

    vector<int> row_vec;

    fstream output;
    output.open(output_File, ios::app);

    cout << "Writing to file\t: " << output_File << endl;
    cout << endl;

    while (start_Co <= end_Value)
    {
        // coordinates.push_back(make_pair(start_Co, end_Co));
        all_start_Co.push_back(start_Co);
        all_end_Co.push_back(end_Co);

        string write = to_string(start_Co) + ":" + to_string(end_Co);
        write_Lines.push_back(write);

        seg_catch_points_ALL.push_back(-1);
        seg_catch_index_ALL.push_back(-1);
        seg_backward_index_ALL.push_back(row_vec);
        seg_forward_index_ALL.push_back(row_vec);

        catch_Point.push_back(1);

        // cout << end_Co << endl;

        if (write_Lines.size() == number_of_genes_Window)
        {
            // Process
            // intialize(coordinates);
            int folder_Start = all_start_Co[0];
            int folder_End = all_end_Co[all_end_Co.size() - 1];
            this->gene_Size = write_Lines.size();

            cout << "Processing windows from " << folder_Start << " to " << folder_End << endl;

            // Collect FILEs
            vector<string> file_List;
            if (folder_Index.size() > 1)
            {
                file_List = function.compound_interpolationSearch(folder_Index, folder_Start, folder_End);
            }
            else
            {
                file_List.push_back(folder_Index[0].second);
            }

            if (this->prev_file_List != file_List)
            {
                if (pre_Lock == 1)
                {
                    pre_Lock = 0;
                }
                else
                {
                    if (test == "T")
                    {
                        free(pre_MA);
                    }
                    else if (test == "FU")
                    {
                        free(pre_MA);
                        free(pre_ne);
                        free(pre_ns);
                    }
                    else if (test == "FA")
                    {
                        free(pre_MA);
                        free(pre_Theta_partials);
                    }
                    else if (test == "N")
                    {
                        free(pre_MA);

                        free(pre_Theta_partials);

                        free(pre_ne);
                        free(pre_ns);
                    }

                    pre_Lock = 1;
                }

                same_Files = "NO";

                tot_Segs = 0;

                all_Lines.clear();
                position_index_Segs.clear();

                concat_Segs.clear();
                all_site_Index.clear();

                start_stop.clear();

                if (Multi_read == "YES")
                {
                    cout << "Intitating multi read based segregating site search" << endl;

                    vector<thread> threads_Multi_read;

                    for (string files : file_List)
                    {
                        // cout << files << endl;
                        threads_Multi_read.push_back(thread{&prometheus::file_Reader_multi, this, files});
                    }

                    for (thread &t : threads_Multi_read)
                    {
                        if (t.joinable())
                        {
                            t.join();
                        }
                    }

                    threads_Multi_read.clear();
                }
                else
                {
                    cout << "Intitating single read based segregating site search" << endl;
                    for (string files : file_List)
                    {
                        fstream file;
                        file.open(files, ios::in);

                        if (file.is_open())
                        {
                            string line;
                            getline(file, line); // skip first header line
                            while (getline(file, line))
                            {
                                all_Lines.push_back(line);
                            }
                            file.close();
                        }
                    }
                }

                // Process Segs
                tot_Segs = all_Lines.size();
                cout << "System is processing and filtering " << tot_Segs << " segregating site(s)" << endl;
                // all_Files_index.clear();

                int GPU_rounds_full = tot_Segs / SNPs_per_Run;
                int GPU_rounds_partial = tot_Segs % SNPs_per_Run;

                // vector<pair<int, int>> start_stop;
                for (int i = 0; i < GPU_rounds_full; i++)
                {
                    int start = i * SNPs_per_Run;
                    int stop = start + SNPs_per_Run;
                    start_stop.push_back(make_pair(start, stop));
                }

                if (GPU_rounds_partial != 0)
                {
                    int start = tot_Segs - GPU_rounds_partial;
                    int stop = tot_Segs;
                    start_stop.push_back(make_pair(start, stop));
                }

                vector<thread> threads_vec;

                for (int rounds = 0; rounds < start_stop.size(); rounds++)
                {
                    concat_Segs.push_back("");
                    int *row;
                    all_site_Index.push_back(row);
                }

                for (int rounds = 0; rounds < start_stop.size(); rounds++)
                {
                    threads_vec.push_back(thread{&prometheus::seg_Concat, this, rounds, start_stop[rounds].first, start_stop[rounds].second});
                }

                for (thread &t : threads_vec)
                {
                    if (t.joinable())
                    {
                        t.join();
                    }
                }

                threads_vec.clear();

                prev_file_List.clear();
                this->prev_file_List = file_List;
            }
            else
            {
                same_Files = "YES";
            }

            // Process test
            if (test == "T")
            {
                process_Tajima();
                cout << "System has completed Tajima's D for the gene(s)" << endl;
            }
            else if (test == "FU")
            {
                process_Fu_Li();
                cout << "System has completed Fu and Li for the gene(s)" << endl;
            }
            else if (test == "FA")
            {
                process_Fay_Wu();
                cout << "System has completed Fay and Wu for the gene(s)" << endl;
            }
            else if (test == "N")
            {
                process_Neutrality();
                cout << "System has completed Neutrality tests for the gene(s)" << endl;
            }

            for (string line : write_Lines)
            {
                output << line << "\n";
            }

            output.flush();

            cout << "System has written the results for the batch" << endl
                 << endl;

            // coordinates.clear();
            all_start_Co.clear();
            all_end_Co.clear();
            write_Lines.clear();

            catch_Point.clear();

            seg_catch_points_ALL.clear();
            seg_catch_index_ALL.clear();
            seg_backward_index_ALL.clear();
            seg_forward_index_ALL.clear();

            gene_Size = 0;
        }

        start_Co = start_Co + step_Size;
        end_Co = start_Co + window_Size;
        // cout << end_Co;
    }

    // check if any remains
    if (write_Lines.size() != 0)
    {
        // intialize(coordinates);
        int folder_Start = all_start_Co[0];
        int folder_End = all_end_Co[all_end_Co.size() - 1];
        this->gene_Size = write_Lines.size();

        cout << "Processing windows from " << folder_Start << " to " << folder_End << endl;
        vector<string> file_List;
        if (folder_Index.size() > 1)
        {
            file_List = function.compound_interpolationSearch(folder_Index, folder_Start, folder_End);
        }
        else
        {
            file_List.push_back(folder_Index[0].second);
        }
        if (this->prev_file_List != file_List)
        {
            if (pre_Lock == 1)
            {
                pre_Lock = 0;
            }
            else
            {
                if (test == "T")
                {
                    free(pre_MA);
                }
                else if (test == "FU")
                {
                    free(pre_MA);
                    free(pre_ne);
                    free(pre_ns);
                }
                else if (test == "FA")
                {
                    free(pre_MA);
                    free(pre_Theta_partials);
                }
                else if (test == "N")
                {
                    free(pre_MA);

                    free(pre_Theta_partials);

                    free(pre_ne);
                    free(pre_ns);
                }

                pre_Lock = 1;
            }

            same_Files = "NO";

            tot_Segs = 0;

            all_Lines.clear();
            position_index_Segs.clear();

            concat_Segs.clear();
            all_site_Index.clear();

            start_stop.clear();

            if (Multi_read == "YES")
            {
                cout << "Intitating multi read based segregating site search" << endl;

                vector<thread> threads_Multi_read;

                for (string files : file_List)
                {
                    // cout << files << endl;
                    threads_Multi_read.push_back(thread{&prometheus::file_Reader_multi, this, files});
                }

                for (thread &t : threads_Multi_read)
                {
                    if (t.joinable())
                    {
                        t.join();
                    }
                }

                threads_Multi_read.clear();
            }
            else
            {
                cout << "Intitating single read based segregating site search" << endl;
                for (string files : file_List)
                {
                    fstream file;
                    file.open(files, ios::in);

                    if (file.is_open())
                    {
                        string line;
                        getline(file, line); // skip first header line
                        while (getline(file, line))
                        {
                            all_Lines.push_back(line);
                        }
                        file.close();
                    }
                }
            }

            // Process Segs
            tot_Segs = all_Lines.size();
            cout << "System is processing and filtering " << tot_Segs << " segregating site(s)" << endl;
            // all_Files_index.clear();

            int GPU_rounds_full = tot_Segs / SNPs_per_Run;
            int GPU_rounds_partial = tot_Segs % SNPs_per_Run;

            // vector<pair<int, int>> start_stop;
            for (int i = 0; i < GPU_rounds_full; i++)
            {
                int start = i * SNPs_per_Run;
                int stop = start + SNPs_per_Run;
                start_stop.push_back(make_pair(start, stop));
            }

            if (GPU_rounds_partial != 0)
            {
                int start = tot_Segs - GPU_rounds_partial;
                int stop = tot_Segs;
                start_stop.push_back(make_pair(start, stop));
            }

            vector<thread> threads_vec;

            for (int rounds = 0; rounds < start_stop.size(); rounds++)
            {
                concat_Segs.push_back("");
                int *row;
                all_site_Index.push_back(row);
            }

            for (int rounds = 0; rounds < start_stop.size(); rounds++)
            {
                threads_vec.push_back(thread{&prometheus::seg_Concat, this, rounds, start_stop[rounds].first, start_stop[rounds].second});
            }

            for (thread &t : threads_vec)
            {
                if (t.joinable())
                {
                    t.join();
                }
            }

            threads_vec.clear();

            prev_file_List.clear();
            this->prev_file_List = file_List;
        }
        else
        {
            same_Files = "YES";
        }

        // Process test
        if (test == "T")
        {
            process_Tajima();
            cout << "System has completed Tajima's D for the gene(s)" << endl;
        }
        else if (test == "FU")
        {
            process_Fu_Li();
            cout << "System has completed Fu and Li for the gene(s)" << endl;
        }
        else if (test == "FA")
        {
            process_Fay_Wu();
            cout << "System has completed Fay and Wu for the gene(s)" << endl;
        }
        else if (test == "N")
        {
            process_Neutrality();
            cout << "System has completed Neutrality tests for the gene(s)" << endl;
        }

        for (string line : write_Lines)
        {
            output << line << "\n";
        }

        output.flush();

        cout << "System has written the results for the batch" << endl
             << endl;

        // coordinates.clear();
        all_start_Co.clear();
        all_end_Co.clear();
        write_Lines.clear();

        catch_Point.clear();

        seg_catch_points_ALL.clear();
        seg_catch_index_ALL.clear();
        seg_backward_index_ALL.clear();
        seg_forward_index_ALL.clear();

        gene_Size = 0;

        // COMPLETE THIS
    }

    output.close();
}

// void prometheus::intialize(vector<pair<int, int>> &coordinates)
// {
//     for (int i = 0; i < coordinates.size(); i++)
//     {
//         // vector<int> row;
//         // Coordinates\tPi\tS\tTajimas_D
//         all_start_Co.push_back(coordinates[i].first);
//         all_end_Co.push_back(coordinates[i].second);
//         string write = coordinates[i].first + ":" + coordinates[i].second;
//         write_Lines.push_back(write);
//     }
// }

vector<string> prometheus::collection_Engine(vector<string> &gene_Collect, string test_Type)
{
    // FIX all_lines

    this->gene_Size = gene_Collect.size();
    // MAKE HDD version : DONE
    for (size_t i = 0; i < gene_Size; i++)
    {
        // multithread
        // vector<string> row;
        vector<int> row;
        // all_Lines.push_back("");
        // all_Files.push_back(row);

        all_start_Co.push_back(-1);
        all_end_Co.push_back(-1);
        write_Lines.push_back("");

        catch_Point.push_back(-1);
        // catch_files.push_back("");
        catch_forward_index.push_back(row);
        catch_back_index.push_back(row);

        seg_catch_points_ALL.push_back(-1);
        seg_catch_index_ALL.push_back(-1);
        seg_backward_index_ALL.push_back(row);
        seg_forward_index_ALL.push_back(row);
    }

    cout << "Intitating file index collection" << endl;
    if (folder_Index.size() > 1)
    {
        vector<thread> threads_vec;
        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            threads_vec.push_back(thread{&prometheus::get_Gene_info_and_catch, this, gene_Collect[gene_ID], gene_ID});
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }
        threads_vec.clear();

        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            // int pos, int start_Co, int end_Co, int gene_ID
            if (catch_Point[gene_ID] != -1)
            {
                threads_vec.push_back(thread{&prometheus::forward_Search, this, catch_Point[gene_ID], all_start_Co[gene_ID], all_end_Co[gene_ID], gene_ID});
                threads_vec.push_back(thread{&prometheus::backward_Search, this, catch_Point[gene_ID], all_start_Co[gene_ID], all_end_Co[gene_ID], gene_ID});
            }
        }
        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }
        threads_vec.clear();

        // clear all catch forward and back after compile

        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            if (catch_Point[gene_ID] != -1)
            {
                threads_vec.push_back(thread{&prometheus::compile, this, gene_ID});
            }
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        // catch_Point.clear();
        // catch_files.clear();
        catch_forward_index.clear();
        catch_back_index.clear();
    }
    else
    {
        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            // vector<string> file_List;
            all_Files_index.insert(0);
        }
    }

    // sort(all_Files_index.begin(), all_Files_index.end());

    if (Multi_read == "YES")
    {
        cout << "Intitating multi read based segregating site search" << endl;

        vector<thread> threads_Multi_read;
        // vector<int> files_READ;

        // for (int gene_ID = 0; gene_ID < gene_Collect.size(); gene_ID++)
        // {
        // vector<int> gene_Files = all_Files_index[gene_ID];

        for (int index : all_Files_index)
        {
            // if (binary_search(files_READ.begin(), files_READ.end(), index))
            // {
            string files = folder_Index[index].second;
            threads_Multi_read.push_back(thread{&prometheus::file_Reader_multi, this, files});
            // files_READ.push_back(index);
            // sort(files_READ.begin(), files_READ.end());
            //  }
        }
        //}

        for (thread &t : threads_Multi_read)
        {
            if (t.joinable())
            {
                t.join();
            }
        }
        // files_READ.clear();
        threads_Multi_read.clear();
    }
    else
    {
        cout << "Intitating single read based segregating site search" << endl;
        file_Reader_single();
    }

    tot_Segs = all_Lines.size();
    cout << "System is processing and filtering " << tot_Segs << " segregating site(s)" << endl;
    all_Files_index.clear();

    int GPU_rounds_full = tot_Segs / SNPs_per_Run;
    int GPU_rounds_partial = tot_Segs % SNPs_per_Run;

    // vector<pair<int, int>> start_stop;
    for (int i = 0; i < GPU_rounds_full; i++)
    {
        int start = i * SNPs_per_Run;
        int stop = start + SNPs_per_Run;
        start_stop.push_back(make_pair(start, stop));
    }

    if (GPU_rounds_partial != 0)
    {
        int start = tot_Segs - GPU_rounds_partial;
        int stop = tot_Segs;
        start_stop.push_back(make_pair(start, stop));
    }

    vector<thread> threads_vec;

    for (int rounds = 0; rounds < start_stop.size(); rounds++)
    {
        concat_Segs.push_back("");
        int *row;
        all_site_Index.push_back(row);
    }

    for (int rounds = 0; rounds < start_stop.size(); rounds++)
    {
        // all_end_Index.push_back(row);
        // char **matrix;
        // all_segs_Matrix.push_back(matrix);
        // max_Track.push_back(0);

        threads_vec.push_back(thread{&prometheus::seg_Concat, this, rounds, start_stop[rounds].first, start_stop[rounds].second});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    if (test_Type == "T")
    {
        process_Tajima();
        cout << "System has completed Tajima's D for the gene(s)" << endl;
    }
    else if (test_Type == "FU")
    {
        process_Fu_Li();
        cout << "System has completed Fu and Li for the gene(s)" << endl;
    }
    else if (test_Type == "FA")
    {
        process_Fay_Wu();
        cout << "System has completed Fay and Wu for the gene(s)" << endl;
    }
    else if (test_Type == "N")
    {
        process_Neutrality();
        cout << "System has completed Neutrality tests for the gene(s)" << endl;
    }

    return write_Lines;
}

__global__ void cuda_neutrality_Prometheus(char *sites, int *index, int num_Segregrating_sites, int *theta_Partials, int *VALID_or_NOT, int *MA_count, int *ne, int *ns, int *cuda_pos_start_Index, int *cuda_pos_end_Index, int start)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Segregrating_sites)
    {
        int column = 0;
        int site_Start = index[tid];
        int site_End = index[tid + 1];

        int i = site_Start;

        while (column < 1)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // POS column
        cuda_pos_start_Index[tid] = i;
        while (column < 2)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // will point to the tab but makes < easier later
        cuda_pos_end_Index[tid] = i - 1;

        char REF = 'N';
        while (column < 3)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        if (sites[i] >= 97)
        {
            REF = sites[i] - 32;
        }
        else
        {
            REF = sites[i];
        }

        char ALT = 'N';
        while (column < 4)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        if (sites[i] >= 97)
        {
            ALT = sites[i] - 32;
        }
        else
        {
            ALT = sites[i];
        }

        // printf("REF: %c", REF);
        // printf("\n");
        // printf("ALT: %c", ALT);
        // printf("\n");

        while (column < 7)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        char AA = 'N';

        while (column < 8)
        {
            if ((sites[i] == 'A' && sites[i + 1] == 'A' && sites[i + 2] == '=') || (sites[i] == 'a' && sites[i + 1] == 'a' && sites[i + 2] == '='))
            {
                char CHECK = 'N';

                if (sites[i + 3] >= 97)
                {
                    CHECK = sites[i + 3] - 32;
                }
                else
                {
                    CHECK = sites[i + 3];
                }

                if (CHECK == 'A' || CHECK == 'T' || CHECK == 'G' || CHECK == 'C')
                {
                    AA = CHECK;
                }
            }

            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        if (AA == 'N')
        {
            AA = REF;
        }

        while (column < 9)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        int ALT_count = 0;
        int REF_count = 0;
        while (i < site_End)
        {
            if (sites[i] == '1')
            {
                ALT_count = ALT_count + 1;
            }
            else if (sites[i] == '0')
            {
                REF_count = REF_count + 1;
            }

            i++;
        }

        if (ALT_count == 0 || REF_count == 0)
        {
            VALID_or_NOT[tid] = 0;
            MA_count[tid + start] = 0;
            theta_Partials[tid + start] = 0;
        }
        else
        {
            VALID_or_NOT[tid] = 1;
            char MA = 'N';

            if (ALT_count < REF_count)
            {
                MA = ALT;
                MA_count[tid + start] = ALT_count;
            }
            else
            {
                MA = REF;
                MA_count[tid + start] = REF_count;
            }

            if (MA_count[tid + start] == 1)
            {
                ns[tid + start] = 1;

                if (MA == AA)
                {
                    ne[tid + start] = 0;
                }
                else
                {
                    ne[tid + start] = 1;
                }
            }
            else
            {
                ne[tid + start] = 0;
                ns[tid + start] = 0;
            }

            int theta_Partial = 0;
            if (MA != AA)
            {
                theta_Partial = MA_count[tid + start];
            }

            theta_Partials[tid + start] = theta_Partial;
        }

        tid += blockDim.x * gridDim.x;
    }
}

void prometheus::process_Neutrality()
{
    vector<thread> threads_vec;

    int *MA_Count, *ne, *ns, *Theta_partials;
    MA_Count = (int *)malloc(tot_Segs * sizeof(int));
    ne = (int *)malloc(tot_Segs * sizeof(int));
    ns = (int *)malloc(tot_Segs * sizeof(int));
    Theta_partials = (int *)malloc(tot_Segs * sizeof(int));

    if (same_Files == "NO")
    {

        int *cuda_MA_Count;
        cudaMallocManaged(&cuda_MA_Count, tot_Segs * sizeof(int));

        int *ne_CUDA, *ns_CUDA;
        cudaMallocManaged(&ne_CUDA, tot_Segs * sizeof(int));
        cudaMallocManaged(&ns_CUDA, tot_Segs * sizeof(int));

        int *cuda_Theta_partials;
        cudaMallocManaged(&cuda_Theta_partials, tot_Segs * sizeof(int));

        for (int rounds = 0; rounds < start_stop.size(); rounds++)
        {
            string Seg_sites = concat_Segs[rounds];
            int *site_Index = all_site_Index[rounds];

            int start = start_stop[rounds].first;
            int stop = start_stop[rounds].second;

            int total_Segs = stop - start;

            char *full_Char;
            full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
            strcpy(full_Char, Seg_sites.c_str());

            char *cuda_full_Char;
            cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
            int *cuda_site_Index;
            cudaMallocManaged(&cuda_site_Index, (total_Segs + 1) * sizeof(int));
            cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
            cudaMemcpy(cuda_site_Index, site_Index, (total_Segs + 1) * sizeof(int), cudaMemcpyHostToDevice);

            int *cuda_VALID_or_NOT, *VALID_or_NOT;
            cudaMallocManaged(&cuda_VALID_or_NOT, total_Segs * sizeof(int));
            VALID_or_NOT = (int *)malloc(total_Segs * sizeof(int));

            int *pos_start_Index, *pos_end_Index, *cuda_pos_start_Index, *cuda_pos_end_Index;
            cudaMallocManaged(&cuda_pos_start_Index, total_Segs * sizeof(int));
            cudaMallocManaged(&cuda_pos_end_Index, total_Segs * sizeof(int));
            pos_start_Index = (int *)malloc(total_Segs * sizeof(int));
            pos_end_Index = (int *)malloc(total_Segs * sizeof(int));

            cuda_neutrality_Prometheus<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_full_Char, cuda_site_Index, total_Segs, cuda_Theta_partials, cuda_VALID_or_NOT, cuda_MA_Count, ne_CUDA, ns_CUDA, cuda_pos_start_Index, cuda_pos_end_Index, start);
            cudaError_t err = cudaGetLastError();

            if (err != cudaSuccess)
            {
                printf("CUDA Error: %s\n", cudaGetErrorString(err));

                // Possibly: exit(-1) if program cannot continue....
            }
            cudaDeviceSynchronize();

            cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(pos_start_Index, cuda_pos_start_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(pos_end_Index, cuda_pos_end_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

            cudaFree(cuda_site_Index);
            cudaFree(cuda_full_Char);
            cudaFree(cuda_VALID_or_NOT);
            cudaFree(cuda_pos_start_Index);
            cudaFree(cuda_pos_end_Index);

            cout << "Round " << rounds + 1 << ": System is indexing " << total_Segs << " processed segregating site(s)" << endl;

            int segs_per_Thread = total_Segs / CPU_cores;
            int remainder = total_Segs % CPU_cores;

            for (int core_ID = 0; core_ID < CPU_cores; core_ID++)
            {
                int start_Seg = core_ID * segs_per_Thread;
                int stop_Seg = start_Seg + segs_per_Thread;

                threads_vec.push_back(thread{&prometheus::seg_Indexer, this, start_Seg, stop_Seg, full_Char, VALID_or_NOT, pos_start_Index, pos_end_Index, start});
            }

            if (remainder != 0)
            {
                int start_Seg = total_Segs - remainder;
                int stop_Seg = total_Segs;

                threads_vec.push_back(thread{&prometheus::seg_Indexer, this, start_Seg, stop_Seg, full_Char, VALID_or_NOT, pos_start_Index, pos_end_Index, start});
            }

            for (thread &t : threads_vec)
            {
                if (t.joinable())
                {
                    t.join();
                }
            }

            threads_vec.clear();

            free(full_Char);
            free(site_Index);
            free(VALID_or_NOT);
            free(pos_start_Index);
            free(pos_end_Index);
        }

        cudaMemcpy(MA_Count, cuda_MA_Count, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(ne, ne_CUDA, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(ns, ns_CUDA, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(Theta_partials, cuda_Theta_partials, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        if (calc_Mode != "FILE")
        {
            pre_MA = (int *)malloc(tot_Segs * sizeof(int));
            cudaMemcpy(pre_MA, cuda_MA_Count, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);

            pre_Theta_partials = (int *)malloc(tot_Segs * sizeof(int));
            cudaMemcpy(pre_Theta_partials, cuda_Theta_partials, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);

            pre_ne = (int *)malloc(tot_Segs * sizeof(int));
            pre_ns = (int *)malloc(tot_Segs * sizeof(int));
            cudaMemcpy(pre_ne, ne_CUDA, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(pre_ns, ns_CUDA, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        }

        cudaFree(cuda_MA_Count);
        cudaFree(ne_CUDA);
        cudaFree(ns_CUDA);
        cudaFree(cuda_Theta_partials);

        sort(position_index_Segs.begin(), position_index_Segs.end());
    }
    else
    {
        memcpy(MA_Count, pre_MA, tot_Segs * sizeof(int));

        memcpy(Theta_partials, pre_Theta_partials, tot_Segs * sizeof(int));

        memcpy(ne, pre_ne, tot_Segs * sizeof(int));
        memcpy(ns, pre_ns, tot_Segs * sizeof(int));
    }

    cout << "Filteration and selection of segregating site(s)" << endl;

    if (this->sliding_Mode == "YES")
    {
        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            threads_vec.push_back(thread{&prometheus::get_POS_VALID, this, all_start_Co[gene_ID], gene_ID});
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();

        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            if (seg_catch_points_ALL[gene_ID] != -1)
            {
                threads_vec.push_back(thread{&prometheus::seg_Search_forward, this, gene_ID});
            }
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }
    else
    {
        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            if (catch_Point[gene_ID] != -1)
            {
                threads_vec.push_back(thread{&prometheus::seg_Search_catch_point, this, gene_ID});
            }
            else
            {
                // add to write_lines the blank
                // Pi\tS\tne\tns\tTotal_iEi\tTajimas_D\tD\tD_star\tF\tF_star\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E
                string write = "0\t0\t0\t0\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
                write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;
            }
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();

        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            if (catch_Point[gene_ID] != -1)
            {
                if (seg_catch_index_ALL[gene_ID] != -1)
                {
                    threads_vec.push_back(thread{&prometheus::seg_Search_forward, this, gene_ID});
                    threads_vec.push_back(thread{&prometheus::seg_Search_backward, this, gene_ID});
                }
                else
                {
                    // add to write_lines the blank
                    string write = "0\t0\t0\t0\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
                    write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;
                }
            }
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }

    // PROCESS NEUTRALITY ALL
    for (size_t gene_ID = 0; gene_ID < gene_Size; gene_ID++)
    {
        if ((catch_Point[gene_ID] != -1) && (seg_catch_index_ALL[gene_ID] != -1))
        {
            threads_vec.push_back(thread{&prometheus::calc_Neutrality_Segs, this, gene_ID, MA_Count, ne, ns, Theta_partials});
        }
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    free(MA_Count);
    free(ne);
    free(ns);
    free(Theta_partials);

    cout << "System has processed all segregating site(s)" << endl;
}

void prometheus::calc_Neutrality_Segs(int gene_ID, int *MA_Count, int *ne, int *ns, int *Theta_partials)
{
    int real_segregrating_Sites = seg_backward_index_ALL[gene_ID].size() + 1 + seg_forward_index_ALL[gene_ID].size();
    float tot_pairwise_Differences = 0;

    int total_iTheta = 0;

    int singletons_ne = 0;
    int singletons_ns = 0;

    for (int index : seg_backward_index_ALL[gene_ID])
    {
        float MAF = (float)MA_Count[index] / (float)N;
        tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));

        total_iTheta = total_iTheta + Theta_partials[index];

        singletons_ne = singletons_ne + ne[index];
        singletons_ns = singletons_ns + ns[index];
    }

    float MAF = (float)MA_Count[seg_catch_index_ALL[gene_ID]] / (float)N;
    tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));

    total_iTheta = total_iTheta + Theta_partials[seg_catch_index_ALL[gene_ID]];

    singletons_ne = singletons_ne + ne[seg_catch_index_ALL[gene_ID]];
    singletons_ns = singletons_ns + ns[seg_catch_index_ALL[gene_ID]];

    for (int index : seg_forward_index_ALL[gene_ID])
    {
        float MAF = (float)MA_Count[index] / (float)N;
        tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));

        total_iTheta = total_iTheta + Theta_partials[index];

        singletons_ne = singletons_ne + ne[index];
        singletons_ns = singletons_ns + ns[index];
    }

    float pi = (float)tot_pairwise_Differences / combinations;
    // TAJIMA
    float Tajima_D = 0;

    // pi = (float)tot_pairwise_Differences / combinations;
    Tajima_D = (float)(pi - (real_segregrating_Sites / an)) / sqrt(((e1 * real_segregrating_Sites) + (e2 * real_segregrating_Sites * (real_segregrating_Sites - 1))));

    // FU LI

    float D = (float)(real_segregrating_Sites - (an * singletons_ne)) / sqrt(((ud * real_segregrating_Sites) + (vd * (pow(real_segregrating_Sites, 2)))));
    float D_star = (float)(((N_float / (N_float - 1)) * real_segregrating_Sites) - (an * singletons_ns)) / sqrt(((ud_star * real_segregrating_Sites) + (vd_star * (pow(real_segregrating_Sites, 2.0)))));

    // float pi = (float)tot_pairwise_Differences / combinations;

    float F = (float)(pi - singletons_ne) / sqrt(((uf * real_segregrating_Sites) + (vf * (pow(real_segregrating_Sites, 2)))));
    float F_star = (float)(pi - (((N_float - 1) / N_float) * singletons_ns)) / sqrt(((uf_star * real_segregrating_Sites) + (vf_star * (pow(real_segregrating_Sites, 2.0)))));

    // FAY WU

    float theta_L = (float)(1 / (N_float - 1)) * (float)total_iTheta;

    float theta_squared = (float)(real_segregrating_Sites * (real_segregrating_Sites - 1)) / (pow(an, 2) + bn);
    float theta_W = (float)real_segregrating_Sites / an;

    // float pi = (float)tot_pairwise_Differences / combinations;

    float VAR_pi_MINUS_theta_L = (float)(((N_float - 2.0) / (6.0 * (N_float - 1.0))) * theta_W) + ((((18.0 * pow(N_float, 2) * ((3.0 * N_float) + 2.0) * bn_plus1) - ((88.0 * pow(N_float, 3)) + (9.0 * pow(N_float, 2)) - (13.0 * N_float) + 6.0)) / (9.0 * N_float * pow(N_float - 1, 2))) * theta_squared);
    float VAR_theta_L_MINUS_theta_W = (float)(((N_float / (2.0 * (N_float - 1.0))) - (1.0 / an)) * theta_W) + (((bn / (pow(an, 2))) + (2.0 * pow((N_float / (N_float - 1.0)), 2) * bn) - ((2.0 * ((N_float * bn) - N_float + 1.0)) / ((N_float - 1.0) * an)) - (((3.0 * N_float) + 1) / (N_float - 1.0))) * theta_squared);

    float H = (float)(pi - theta_L) / (sqrt(VAR_pi_MINUS_theta_L));
    float E = (float)(theta_L - theta_W) / (sqrt(VAR_theta_L_MINUS_theta_W));

    // Pi\tS\tne\tns\tTotal_iEi\tTajimas_D\tD\tD_star\tF\tF_star\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E"
    string write = to_string(pi) + "\t" + to_string(real_segregrating_Sites) + "\t" + to_string(singletons_ne) + "\t" + to_string(singletons_ns) + "\t" + to_string(total_iTheta) + "\t" + to_string(Tajima_D) + "\t" + to_string(D) + "\t" + to_string(D_star) + "\t" + to_string(F) + "\t" + to_string(F_star) + "\t" + to_string(H) + "\t" + to_string(E);
    unique_lock<shared_mutex> ul(g_mutex);
    write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;
}

__global__ void cuda_fay_wu_Prometheus(char *sites, int *index, int num_Segregrating_sites, int *theta_Partials, int *VALID_or_NOT, int *MA_count, int *cuda_pos_start_Index, int *cuda_pos_end_Index, int start)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Segregrating_sites)
    {
        int column = 0;
        int site_Start = index[tid];
        int site_End = index[tid + 1];

        int i = site_Start;

        while (column < 1)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // POS column
        cuda_pos_start_Index[tid] = i;
        while (column < 2)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // will point to the tab but makes < easier later
        cuda_pos_end_Index[tid] = i - 1;

        char REF = 'N';
        while (column < 3)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        if (sites[i] >= 97)
        {
            REF = sites[i] - 32;
        }
        else
        {
            REF = sites[i];
        }

        char ALT = 'N';
        while (column < 4)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        if (sites[i] >= 97)
        {
            ALT = sites[i] - 32;
        }
        else
        {
            ALT = sites[i];
        }

        while (column < 7)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        char AA = 'N';

        while (column < 8)
        {
            if ((sites[i] == 'A' && sites[i + 1] == 'A' && sites[i + 2] == '=') || (sites[i] == 'a' && sites[i + 1] == 'a' && sites[i + 2] == '='))
            {
                char CHECK = 'N';

                if (sites[i + 3] >= 97)
                {
                    CHECK = sites[i + 3] - 32;
                }
                else
                {
                    CHECK = sites[i + 3];
                }

                if (CHECK == 'A' || CHECK == 'T' || CHECK == 'G' || CHECK == 'C')
                {
                    AA = CHECK;
                }
            }

            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        if (AA == 'N')
        {
            AA = REF;
        }

        while (column < 9)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        int ALT_count = 0;
        int REF_count = 0;
        while (i < site_End)
        {
            if (sites[i] == '1')
            {
                ALT_count = ALT_count + 1;
            }
            else if (sites[i] == '0')
            {
                REF_count = REF_count + 1;
            }

            i++;
        }

        if (ALT_count == 0 || REF_count == 0)
        {
            VALID_or_NOT[tid] = 0;
            MA_count[tid + start] = 0;
            theta_Partials[tid + start] = 0;
        }
        else
        {
            VALID_or_NOT[tid] = 1;
            char MA = 'N';

            if (ALT_count < REF_count)
            {
                MA = ALT;
                MA_count[tid + start] = ALT_count;
            }
            else
            {
                MA = REF;
                MA_count[tid + start] = REF_count;
            }

            int theta_Partial = 0;
            if (MA != AA)
            {
                theta_Partial = MA_count[tid + start];
            }

            theta_Partials[tid + start] = theta_Partial;
        }

        tid += blockDim.x * gridDim.x;
    }
}

void prometheus::process_Fay_Wu()
{
    vector<thread> threads_vec;

    int *MA_Count;
    MA_Count = (int *)malloc(tot_Segs * sizeof(int));

    int *Theta_partials;
    Theta_partials = (int *)malloc(tot_Segs * sizeof(int));

    if (same_Files == "NO")
    {
        int *cuda_MA_Count;
        cudaMallocManaged(&cuda_MA_Count, tot_Segs * sizeof(int));

        int *cuda_Theta_partials;
        cudaMallocManaged(&cuda_Theta_partials, tot_Segs * sizeof(int));

        for (int rounds = 0; rounds < start_stop.size(); rounds++)
        {
            string Seg_sites = concat_Segs[rounds];
            int *site_Index = all_site_Index[rounds];

            int start = start_stop[rounds].first;
            int stop = start_stop[rounds].second;

            int total_Segs = stop - start;

            char *full_Char;
            full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
            strcpy(full_Char, Seg_sites.c_str());

            char *cuda_full_Char;
            cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
            int *cuda_site_Index;
            cudaMallocManaged(&cuda_site_Index, (total_Segs + 1) * sizeof(int));
            cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
            cudaMemcpy(cuda_site_Index, site_Index, (total_Segs + 1) * sizeof(int), cudaMemcpyHostToDevice);

            int *cuda_VALID_or_NOT, *VALID_or_NOT;
            cudaMallocManaged(&cuda_VALID_or_NOT, total_Segs * sizeof(int));
            VALID_or_NOT = (int *)malloc(total_Segs * sizeof(int));

            int *pos_start_Index, *pos_end_Index, *cuda_pos_start_Index, *cuda_pos_end_Index;
            cudaMallocManaged(&cuda_pos_start_Index, total_Segs * sizeof(int));
            cudaMallocManaged(&cuda_pos_end_Index, total_Segs * sizeof(int));
            pos_start_Index = (int *)malloc(total_Segs * sizeof(int));
            pos_end_Index = (int *)malloc(total_Segs * sizeof(int));

            cuda_fay_wu_Prometheus<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_full_Char, cuda_site_Index, total_Segs, cuda_Theta_partials, cuda_VALID_or_NOT, cuda_MA_Count, cuda_pos_start_Index, cuda_pos_end_Index, start);
            cudaError_t err = cudaGetLastError();

            if (err != cudaSuccess)
            {
                printf("CUDA Error: %s\n", cudaGetErrorString(err));

                // Possibly: exit(-1) if program cannot continue....
            }
            cudaDeviceSynchronize();

            cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(pos_start_Index, cuda_pos_start_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(pos_end_Index, cuda_pos_end_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

            cudaFree(cuda_site_Index);
            cudaFree(cuda_full_Char);
            cudaFree(cuda_VALID_or_NOT);
            cudaFree(cuda_pos_start_Index);
            cudaFree(cuda_pos_end_Index);

            cout << "Round " << rounds + 1 << ": System is indexing " << total_Segs << " processed segregating site(s)" << endl;

            int segs_per_Thread = total_Segs / CPU_cores;
            int remainder = total_Segs % CPU_cores;

            for (int core_ID = 0; core_ID < CPU_cores; core_ID++)
            {
                int start_Seg = core_ID * segs_per_Thread;
                int stop_Seg = start_Seg + segs_per_Thread;

                threads_vec.push_back(thread{&prometheus::seg_Indexer, this, start_Seg, stop_Seg, full_Char, VALID_or_NOT, pos_start_Index, pos_end_Index, start});
            }

            if (remainder != 0)
            {
                int start_Seg = total_Segs - remainder;
                int stop_Seg = total_Segs;

                threads_vec.push_back(thread{&prometheus::seg_Indexer, this, start_Seg, stop_Seg, full_Char, VALID_or_NOT, pos_start_Index, pos_end_Index, start});
            }

            for (thread &t : threads_vec)
            {
                if (t.joinable())
                {
                    t.join();
                }
            }

            threads_vec.clear();

            free(full_Char);
            free(site_Index);
            free(VALID_or_NOT);
            free(pos_start_Index);
            free(pos_end_Index);
        }

        cudaMemcpy(MA_Count, cuda_MA_Count, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(Theta_partials, cuda_Theta_partials, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        if (calc_Mode != "FILE")
        {
            pre_MA = (int *)malloc(tot_Segs * sizeof(int));
            cudaMemcpy(pre_MA, cuda_MA_Count, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
            pre_Theta_partials = (int *)malloc(tot_Segs * sizeof(int));
            cudaMemcpy(pre_Theta_partials, cuda_Theta_partials, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        }

        cudaFree(cuda_MA_Count);
        cudaFree(cuda_Theta_partials);

        sort(position_index_Segs.begin(), position_index_Segs.end());
    }
    else
    {
        memcpy(MA_Count, pre_MA, tot_Segs * sizeof(int));
        memcpy(Theta_partials, pre_Theta_partials, tot_Segs * sizeof(int));
    }

    cout << "Filteration and selection of segregating site(s)" << endl;

    if (this->sliding_Mode == "YES")
    {
        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            threads_vec.push_back(thread{&prometheus::get_POS_VALID, this, all_start_Co[gene_ID], gene_ID});
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();

        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            if (seg_catch_points_ALL[gene_ID] != -1)
            {
                threads_vec.push_back(thread{&prometheus::seg_Search_forward, this, gene_ID});
            }
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }
    else
    {
        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            if (catch_Point[gene_ID] != -1)
            {
                threads_vec.push_back(thread{&prometheus::seg_Search_catch_point, this, gene_ID});
            }
            else
            {
                // add to write_lines the blank
                // Pi\tS\tTotal_iEi\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E
                string write = "0\t0\t0\tNA\tNA";
                write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;
            }
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();

        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            if (catch_Point[gene_ID] != -1)
            {
                if (seg_catch_index_ALL[gene_ID] != -1)
                {
                    threads_vec.push_back(thread{&prometheus::seg_Search_forward, this, gene_ID});
                    threads_vec.push_back(thread{&prometheus::seg_Search_backward, this, gene_ID});
                }
                else
                {
                    // add to write_lines the blank
                    string write = "0\t0\t0\tNA\tNA";
                    write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;
                }
            }
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }

    // Process Fay Wu
    for (size_t gene_ID = 0; gene_ID < gene_Size; gene_ID++)
    {
        if ((catch_Point[gene_ID] != -1) && (seg_catch_index_ALL[gene_ID] != -1))
        {
            threads_vec.push_back(thread{&prometheus::calc_Fay_Wu_Segs, this, gene_ID, MA_Count, Theta_partials});
        }
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    free(MA_Count);
    free(Theta_partials);

    cout << "System has processed all segregating site(s)" << endl;
}

void prometheus::calc_Fay_Wu_Segs(int gene_ID, int *MA_Count, int *Theta_partials)
{
    int real_segregrating_Sites = seg_backward_index_ALL[gene_ID].size() + 1 + seg_forward_index_ALL[gene_ID].size();
    float tot_pairwise_Differences = 0;
    int total_iTheta = 0;

    for (int index : seg_backward_index_ALL[gene_ID])
    {
        float MAF = (float)MA_Count[index] / (float)N;
        tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));

        total_iTheta = total_iTheta + Theta_partials[index];
    }

    float MAF = (float)MA_Count[seg_catch_index_ALL[gene_ID]] / (float)N;
    tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));

    total_iTheta = total_iTheta + Theta_partials[seg_catch_index_ALL[gene_ID]];

    for (int index : seg_forward_index_ALL[gene_ID])
    {
        float MAF = (float)MA_Count[index] / (float)N;
        tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));

        total_iTheta = total_iTheta + Theta_partials[index];
    }

    float theta_L = (float)(1 / (N_float - 1)) * (float)total_iTheta;

    float theta_squared = (float)(real_segregrating_Sites * (real_segregrating_Sites - 1)) / (pow(an, 2) + bn);
    float theta_W = (float)real_segregrating_Sites / an;

    float pi = (float)tot_pairwise_Differences / combinations;

    float VAR_pi_MINUS_theta_L = (float)(((N_float - 2.0) / (6.0 * (N_float - 1.0))) * theta_W) + ((((18.0 * pow(N_float, 2) * ((3.0 * N_float) + 2.0) * bn_plus1) - ((88.0 * pow(N_float, 3)) + (9.0 * pow(N_float, 2)) - (13.0 * N_float) + 6.0)) / (9.0 * N_float * pow(N_float - 1, 2))) * theta_squared);
    float VAR_theta_L_MINUS_theta_W = (float)(((N_float / (2.0 * (N_float - 1.0))) - (1.0 / an)) * theta_W) + (((bn / (pow(an, 2))) + (2.0 * pow((N_float / (N_float - 1.0)), 2) * bn) - ((2.0 * ((N_float * bn) - N_float + 1.0)) / ((N_float - 1.0) * an)) - (((3.0 * N_float) + 1) / (N_float - 1.0))) * theta_squared);

    float H = (float)(pi - theta_L) / (sqrt(VAR_pi_MINUS_theta_L));
    float E = (float)(theta_L - theta_W) / (sqrt(VAR_theta_L_MINUS_theta_W));

    // Pi\tS\tTotal_iEi\tFay_Wu_Normalized_H\tFay_Wu_Normalized_E
    string write = to_string(pi) + "\t" + to_string(real_segregrating_Sites) + "\t" + to_string(total_iTheta) + "\t" + to_string(H) + "\t" + to_string(E);
    unique_lock<shared_mutex> ul(g_mutex);
    write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;
}

__global__ void cuda_fu_li_Prometheus(char *sites, int *index, int tot_Segregrating_sites, int *VALID_or_NOT, int *MA_count, int *ne, int *ns, int *cuda_pos_start_Index, int *cuda_pos_end_Index, int start)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < tot_Segregrating_sites)
    {
        int column = 0;
        int site_Start = index[tid];
        int site_End = index[tid + 1];

        int i = site_Start;

        while (column < 1)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // POS column
        cuda_pos_start_Index[tid] = i;
        while (column < 2)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // will point to the tab but makes < easier later
        cuda_pos_end_Index[tid] = i - 1;

        char REF = 'N';
        while (column < 3)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        if (sites[i] >= 97)
        {
            REF = sites[i] - 32;
        }
        else
        {
            REF = sites[i];
        }

        char ALT = 'N';
        while (column < 4)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        if (sites[i] >= 97)
        {
            ALT = sites[i] - 32;
        }
        else
        {
            ALT = sites[i];
        }

        while (column < 7)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        char AA = 'N';

        while (column < 8)
        {
            if ((sites[i] == 'A' && sites[i + 1] == 'A' && sites[i + 2] == '=') || (sites[i] == 'a' && sites[i + 1] == 'a' && sites[i + 2] == '='))
            {
                char CHECK = 'N';

                if (sites[i + 3] >= 97)
                {
                    CHECK = sites[i + 3] - 32;
                }
                else
                {
                    CHECK = sites[i + 3];
                }

                if (CHECK == 'A' || CHECK == 'T' || CHECK == 'G' || CHECK == 'C')
                {
                    AA = CHECK;
                }
            }

            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        if (AA == 'N')
        {
            AA = REF;
        }

        while (column < 9)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        int ALT_count = 0;
        int REF_count = 0;
        while (i < site_End)
        {
            if (sites[i] == '1')
            {
                ALT_count = ALT_count + 1;
            }
            else if (sites[i] == '0')
            {
                REF_count = REF_count + 1;
            }

            i++;
        }

        if (ALT_count == 0 || REF_count == 0)
        {
            VALID_or_NOT[tid] = 0;
            MA_count[tid + start] = 0;
            ne[tid + start] = 0;
            ns[tid + start] = 0;
        }
        else
        {
            VALID_or_NOT[tid] = 1;

            char MA = 'N';

            if (ALT_count < REF_count)
            {
                MA_count[tid + start] = ALT_count;
                MA = ALT;
            }
            else
            {
                MA = REF;
                MA_count[tid + start] = REF_count;
            }

            if (MA_count[tid + start] == 1)
            {
                ns[tid + start] = 1;

                if (MA == AA)
                {
                    ne[tid + start] = 0;
                }
                else
                {
                    ne[tid + start] = 1;
                }
            }
            else
            {
                ne[tid + start] = 0;
                ns[tid + start] = 0;
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

void prometheus::process_Fu_Li()
{
    vector<thread> threads_vec;

    int *MA_Count;
    MA_Count = (int *)malloc(tot_Segs * sizeof(int));

    int *ne, *ns;
    ne = (int *)malloc(tot_Segs * sizeof(int));
    ns = (int *)malloc(tot_Segs * sizeof(int));

    if (same_Files == "NO")
    {
        int *cuda_MA_Count;
        cudaMallocManaged(&cuda_MA_Count, tot_Segs * sizeof(int));

        int *ne_CUDA, *ns_CUDA;
        cudaMallocManaged(&ne_CUDA, tot_Segs * sizeof(int));
        cudaMallocManaged(&ns_CUDA, tot_Segs * sizeof(int));

        for (int rounds = 0; rounds < start_stop.size(); rounds++)
        {
            string Seg_sites = concat_Segs[rounds];
            // cout << Seg_sites << endl;
            int *site_Index = all_site_Index[rounds];

            int start = start_stop[rounds].first;
            int stop = start_stop[rounds].second;

            int total_Segs = stop - start;

            char *full_Char;
            full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
            strcpy(full_Char, Seg_sites.c_str());

            char *cuda_full_Char;
            cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
            int *cuda_site_Index;
            cudaMallocManaged(&cuda_site_Index, (total_Segs + 1) * sizeof(int));
            cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
            cudaMemcpy(cuda_site_Index, site_Index, (total_Segs + 1) * sizeof(int), cudaMemcpyHostToDevice);

            int *cuda_VALID_or_NOT, *VALID_or_NOT;
            cudaMallocManaged(&cuda_VALID_or_NOT, total_Segs * sizeof(int));
            VALID_or_NOT = (int *)malloc(total_Segs * sizeof(int));

            int *pos_start_Index, *pos_end_Index, *cuda_pos_start_Index, *cuda_pos_end_Index;
            cudaMallocManaged(&cuda_pos_start_Index, total_Segs * sizeof(int));
            cudaMallocManaged(&cuda_pos_end_Index, total_Segs * sizeof(int));
            pos_start_Index = (int *)malloc(total_Segs * sizeof(int));
            pos_end_Index = (int *)malloc(total_Segs * sizeof(int));

            cuda_fu_li_Prometheus<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_full_Char, cuda_site_Index, total_Segs, cuda_VALID_or_NOT, cuda_MA_Count, ne_CUDA, ns_CUDA, cuda_pos_start_Index, cuda_pos_end_Index, start);
            cudaError_t err = cudaGetLastError();

            if (err != cudaSuccess)
            {
                printf("CUDA Error: %s\n", cudaGetErrorString(err));

                // Possibly: exit(-1) if program cannot continue....
            }
            cudaDeviceSynchronize();

            cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(pos_start_Index, cuda_pos_start_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(pos_end_Index, cuda_pos_end_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

            cudaFree(cuda_site_Index);
            cudaFree(cuda_full_Char);
            cudaFree(cuda_VALID_or_NOT);
            cudaFree(cuda_pos_start_Index);
            cudaFree(cuda_pos_end_Index);

            cout << "Round " << rounds + 1 << ": System is indexing " << total_Segs << " processed segregating site(s)" << endl;

            int segs_per_Thread = total_Segs / CPU_cores;
            int remainder = total_Segs % CPU_cores;

            for (int core_ID = 0; core_ID < CPU_cores; core_ID++)
            {
                int start_Seg = core_ID * segs_per_Thread;
                int stop_Seg = start_Seg + segs_per_Thread;

                threads_vec.push_back(thread{&prometheus::seg_Indexer, this, start_Seg, stop_Seg, full_Char, VALID_or_NOT, pos_start_Index, pos_end_Index, start});
            }

            if (remainder != 0)
            {
                int start_Seg = total_Segs - remainder;
                int stop_Seg = total_Segs;

                threads_vec.push_back(thread{&prometheus::seg_Indexer, this, start_Seg, stop_Seg, full_Char, VALID_or_NOT, pos_start_Index, pos_end_Index, start});
            }

            for (thread &t : threads_vec)
            {
                if (t.joinable())
                {
                    t.join();
                }
            }

            threads_vec.clear();

            free(full_Char);
            free(site_Index);
            free(VALID_or_NOT);
            free(pos_start_Index);
            free(pos_end_Index);
        }

        cudaMemcpy(MA_Count, cuda_MA_Count, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(ne, ne_CUDA, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(ns, ns_CUDA, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        if (calc_Mode != "FILE")
        {
            pre_MA = (int *)malloc(tot_Segs * sizeof(int));
            cudaMemcpy(pre_MA, cuda_MA_Count, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
            pre_ne = (int *)malloc(tot_Segs * sizeof(int));
            pre_ns = (int *)malloc(tot_Segs * sizeof(int));
            cudaMemcpy(pre_ne, ne_CUDA, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(pre_ns, ns_CUDA, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        }

        cudaFree(cuda_MA_Count);
        cudaFree(ne_CUDA);
        cudaFree(ns_CUDA);

        sort(position_index_Segs.begin(), position_index_Segs.end());
    }
    else
    {
        memcpy(MA_Count, pre_MA, tot_Segs * sizeof(int));
        memcpy(ne, pre_ne, tot_Segs * sizeof(int));
        memcpy(ns, pre_ns, tot_Segs * sizeof(int));
    }

    cout << "Filteration and selection of segregating site(s)" << endl;

    if (this->sliding_Mode == "YES")
    {
        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            threads_vec.push_back(thread{&prometheus::get_POS_VALID, this, all_start_Co[gene_ID], gene_ID});
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();

        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            if (seg_catch_points_ALL[gene_ID] != -1)
            {
                threads_vec.push_back(thread{&prometheus::seg_Search_forward, this, gene_ID});
            }
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }
    else
    {

        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            if (catch_Point[gene_ID] != -1)
            {
                threads_vec.push_back(thread{&prometheus::seg_Search_catch_point, this, gene_ID});
            }
            else
            {
                // add to write_lines the blank
                string write = "0\t0\t0\t0\tNA\tNA\tNA\tNA";
                write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;
            }
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();

        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            if (catch_Point[gene_ID] != -1)
            {
                if (seg_catch_index_ALL[gene_ID] != -1)
                {
                    threads_vec.push_back(thread{&prometheus::seg_Search_forward, this, gene_ID});
                    threads_vec.push_back(thread{&prometheus::seg_Search_backward, this, gene_ID});
                }
                else
                {
                    // add to write_lines the blank
                    string write = "0\t0\t0\t0\tNA\tNA\tNA\tNA";
                    write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;
                }
            }
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }

    // Process Fu Li
    for (size_t gene_ID = 0; gene_ID < gene_Size; gene_ID++)
    {
        if ((catch_Point[gene_ID] != -1) && (seg_catch_index_ALL[gene_ID] != -1))
        {
            threads_vec.push_back(thread{&prometheus::calc_Fu_Li_Segs, this, gene_ID, MA_Count, ne, ns});
        }
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    free(MA_Count);
    free(ne);
    free(ns);

    cout << "System has processed all segregating site(s)" << endl;
}

void prometheus::calc_Fu_Li_Segs(int gene_ID, int *MA_Count, int *ne, int *ns)
{
    // int real_segregrating_Sites = 0;
    float tot_pairwise_Differences = 0;
    int singletons_ne = 0;
    int singletons_ns = 0;

    int real_segregrating_Sites = seg_backward_index_ALL[gene_ID].size() + 1 + seg_forward_index_ALL[gene_ID].size();

    for (int index : seg_backward_index_ALL[gene_ID])
    {
        float MAF = (float)MA_Count[index] / N;
        tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));

        singletons_ne = singletons_ne + ne[index];
        singletons_ns = singletons_ns + ns[index];
    }

    float MAF = (float)MA_Count[seg_catch_index_ALL[gene_ID]] / (float)N;
    tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));

    singletons_ne = singletons_ne + ne[seg_catch_index_ALL[gene_ID]];
    singletons_ns = singletons_ns + ns[seg_catch_index_ALL[gene_ID]];

    for (int index : seg_forward_index_ALL[gene_ID])
    {
        float MAF = (float)MA_Count[index] / (float)N;
        tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));

        singletons_ne = singletons_ne + ne[index];
        singletons_ns = singletons_ns + ns[index];
    }

    float D = (float)(real_segregrating_Sites - (an * singletons_ne)) / sqrt(((ud * real_segregrating_Sites) + (vd * (pow(real_segregrating_Sites, 2)))));
    float D_star = (float)(((N_float / (N_float - 1)) * real_segregrating_Sites) - (an * singletons_ns)) / sqrt(((ud_star * real_segregrating_Sites) + (vd_star * (pow(real_segregrating_Sites, 2.0)))));

    float pi = (float)tot_pairwise_Differences / combinations;

    float F = (float)(pi - singletons_ne) / sqrt(((uf * real_segregrating_Sites) + (vf * (pow(real_segregrating_Sites, 2)))));
    float F_star = (float)(pi - (((N_float - 1) / N_float) * singletons_ns)) / sqrt(((uf_star * real_segregrating_Sites) + (vf_star * (pow(real_segregrating_Sites, 2.0)))));

    string write = to_string(pi) + "\t" + to_string(real_segregrating_Sites) + "\t" + to_string(singletons_ne) + "\t" + to_string(singletons_ns) + "\t" + to_string(D) + "\t" + to_string(D_star) + "\t" + to_string(F) + "\t" + to_string(F_star);
    unique_lock<shared_mutex> ul(g_mutex);
    write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;
}

void prometheus::compile(int gene_ID)
{
    unique_lock<shared_mutex> ul(g_mutex);

    for (int file_index : catch_back_index[gene_ID])
    {
        all_Files_index.insert(file_index);
    }

    all_Files_index.insert(catch_Point[gene_ID]);

    for (int file_index : catch_forward_index[gene_ID])
    {
        all_Files_index.insert(file_index);
    }
}

void prometheus::file_Reader_multi(string files)
{
    vector<string> lines_Collected;

    fstream file;
    file.open(files, ios::in);
    if (file.is_open())
    {
        string line;
        getline(file, line); // skip first header line
        while (getline(file, line))
        {
            lines_Collected.push_back(line);
        }
        file.close();
    }

    unique_lock<shared_mutex> ul(g_mutex);
    for (string lines_Read : lines_Collected)
    {
        all_Lines.push_back(lines_Read);
    }
}

void prometheus::file_Reader_single()
{
    // vector<int> files_READ;

    // for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
    // {
    // cout << gene_ID<<endl;
    // vector<int> gene_Files_Indexes = all_Files_index[gene_ID];

    for (int index : all_Files_index)
    {
        // cout << files << endl;

        // if (binary_search(files_READ.begin(), files_READ.end(), index))
        // {
        string files = folder_Index[index].second;
        // cout << files << endl;
        //  files_READ.push_back(index);

        fstream file;
        file.open(files, ios::in);

        if (file.is_open())
        {
            string line;
            getline(file, line); // skip first header line
            while (getline(file, line))
            {
                all_Lines.push_back(line);
            }
            file.close();
        }
        // sort(files_READ.begin(), files_READ.end());
        //}
    }
    //}
}

void prometheus::erase()
{
    // all_Files.clear();
    all_Lines.clear();
    all_start_Co.clear();
    all_end_Co.clear();
    tot_Segs = 0;
    write_Lines.clear();

    gene_Size = 0;

    catch_Point.clear();
    catch_forward_index.clear();
    catch_back_index.clear();

    seg_catch_points_ALL.clear();
    seg_catch_index_ALL.clear();
    seg_backward_index_ALL.clear();
    seg_forward_index_ALL.clear();

    all_Files_index.clear();
    position_index_Segs.clear();

    concat_Segs.clear();
    all_site_Index.clear();
    start_stop.clear();

    // all_end_Index.clear();
    // all_segs_Matrix.clear();
    // max_Track.clear();
}

__global__ void cuda_tajima_Prometheus(char *sites, int *index, int tot_Segregrating_sites, int *VALID_or_NOT, int *MA_count, int *cuda_pos_start_Index, int *cuda_pos_end_Index, int start)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < tot_Segregrating_sites)
    {
        // printf("run\n");
        int column = 0;
        int site_Start = index[tid];
        int site_End = index[tid + 1];

        int i = site_Start;

        while (column < 1)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // POS column
        cuda_pos_start_Index[tid] = i;
        while (column < 2)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // will point to the tab but makes < easier later
        cuda_pos_end_Index[tid] = i - 1;

        while (column < 9)
        {
            if (sites[i] == '\t')
            {
                column++;
            }
            i++;
        }

        // printf("Column 1: %c\n", sites[i]);

        int ALT_count = 0;
        int REF_count = 0;
        while (i < site_End)
        {
            if (sites[i] == '1')
            {
                ALT_count = ALT_count + 1;
            }
            else if (sites[i] == '0')
            {
                REF_count = REF_count + 1;
            }

            i++;
        }

        if (ALT_count == 0 || REF_count == 0)
        {
            VALID_or_NOT[tid] = 0;
            MA_count[tid + start] = 0;
        }
        else
        {
            VALID_or_NOT[tid] = 1;

            if (ALT_count < REF_count)
            {
                MA_count[tid + start] = ALT_count;
            }
            else
            {
                MA_count[tid + start] = REF_count;
            }
        }
        // printf("run\n");
        tid += blockDim.x * gridDim.x;
    }
}

// void prometheus::seg_Concat_New(int round_ID, int start_Seg, int stop_Seg)
// {
//     int *end_Index;
//     char **seg_Matrix;
//     int max = 0;

//     int total_Segs = stop_Seg - start_Seg;
//     end_Index = (int *)malloc(total_Segs * sizeof(int));

//     // rows = number of Segs
//     seg_Matrix = (char **)malloc(total_Segs * sizeof(char *));

//     // columns = each char
//     for (int i = start_Seg; i < stop_Seg; i++)
//     {
//         seg_Matrix[i - start_Seg] = (char *)malloc((all_Lines[i].size() + 1) * sizeof(char));
//         end_Index[i - start_Seg] = all_Lines[i].size();

//         char *full_Char;
//         full_Char = (char *)malloc((all_Lines[i].size() + 1) * sizeof(char));
//         strcpy(full_Char, all_Lines[i].c_str());

//         seg_Matrix[i - start_Seg] = full_Char;

//         if (all_Lines[i].size() > max)
//         {
//             max = all_Lines[i].size();
//         }
//     }

//     unique_lock<shared_mutex> ul(g_mutex);
//     all_end_Index[round_ID] = end_Index;
//     all_segs_Matrix[round_ID] = seg_Matrix;
//     if (max > max_Track[round_ID])
//     {
//         max_Track[round_ID] = max;
//     }
// }

void prometheus::seg_Concat(int round_ID, int start_Seg, int stop_Seg)
{
    int *site_Index;
    string Seg_sites = "";

    int total_Segs = stop_Seg - start_Seg;

    site_Index = (int *)malloc((total_Segs + 1) * sizeof(int));
    site_Index[0] = 0;

    // vector<string> total_Segregrating_sites = all_Lines;
    for (size_t i = start_Seg; i < stop_Seg; i++)
    {
        Seg_sites.append(all_Lines[i]);
        site_Index[i + 1 - start_Seg] = site_Index[i - start_Seg] + all_Lines[i].size();
        // if (round_ID == 1)
        // {
        //     cout << site_Index[i] << endl;
        // }
    }

    unique_lock<shared_mutex> ul(g_mutex);
    all_site_Index[round_ID] = site_Index;
    concat_Segs[round_ID] = Seg_sites;
}

void prometheus::process_Tajima()
{
    vector<thread> threads_vec;

    int *MA_Count = (int *)malloc(tot_Segs * sizeof(int));

    if (same_Files == "NO")
    {
        int *cuda_MA_Count;
        cudaMallocManaged(&cuda_MA_Count, tot_Segs * sizeof(int));
        // MA_Count = (int *)malloc(tot_Segs * sizeof(int));

        for (int rounds = 0; rounds < start_stop.size(); rounds++)
        {
            string Seg_sites = concat_Segs[rounds];
            // cout << Seg_sites << endl;
            int *site_Index = all_site_Index[rounds];

            int start = start_stop[rounds].first;
            int stop = start_stop[rounds].second;

            int total_Segs = stop - start;

            char *full_Char;
            full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
            strcpy(full_Char, Seg_sites.c_str());

            char *cuda_full_Char;
            cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
            int *cuda_site_Index;
            cudaMallocManaged(&cuda_site_Index, (total_Segs + 1) * sizeof(int));
            cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
            cudaMemcpy(cuda_site_Index, site_Index, (total_Segs + 1) * sizeof(int), cudaMemcpyHostToDevice);

            int *cuda_VALID_or_NOT, *VALID_or_NOT;
            cudaMallocManaged(&cuda_VALID_or_NOT, total_Segs * sizeof(int));
            VALID_or_NOT = (int *)malloc(total_Segs * sizeof(int));

            int *pos_start_Index, *pos_end_Index, *cuda_pos_start_Index, *cuda_pos_end_Index;
            cudaMallocManaged(&cuda_pos_start_Index, total_Segs * sizeof(int));
            cudaMallocManaged(&cuda_pos_end_Index, total_Segs * sizeof(int));
            pos_start_Index = (int *)malloc(total_Segs * sizeof(int));
            pos_end_Index = (int *)malloc(total_Segs * sizeof(int));

            cuda_tajima_Prometheus<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_full_Char, cuda_site_Index, total_Segs, cuda_VALID_or_NOT, cuda_MA_Count, cuda_pos_start_Index, cuda_pos_end_Index, start);

            cudaError_t err = cudaGetLastError();

            if (err != cudaSuccess)
            {
                printf("CUDA Error: %s\n", cudaGetErrorString(err));

                // Possibly: exit(-1) if program cannot continue....
            }
            cudaDeviceSynchronize();

            cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(pos_start_Index, cuda_pos_start_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(pos_end_Index, cuda_pos_end_Index, total_Segs * sizeof(int), cudaMemcpyDeviceToHost);

            cudaFree(cuda_site_Index);
            cudaFree(cuda_full_Char);
            cudaFree(cuda_VALID_or_NOT);
            cudaFree(cuda_pos_start_Index);
            cudaFree(cuda_pos_end_Index);

            cout << "Round " << rounds + 1 << ": System is indexing " << total_Segs << " processed segregating site(s)" << endl;
            // cout << Seg_sites.size() << endl;

            int segs_per_Thread = total_Segs / CPU_cores;
            int remainder = total_Segs % CPU_cores;

            for (int core_ID = 0; core_ID < CPU_cores; core_ID++)
            {
                int start_Seg = core_ID * segs_per_Thread;
                // cout << start_Seg << endl;
                int stop_Seg = start_Seg + segs_per_Thread;

                threads_vec.push_back(thread{&prometheus::seg_Indexer, this, start_Seg, stop_Seg, full_Char, VALID_or_NOT, pos_start_Index, pos_end_Index, start});
            }

            if (remainder != 0)
            {
                int start_Seg = total_Segs - remainder;
                // cout << start_Seg << endl;
                int stop_Seg = total_Segs;

                threads_vec.push_back(thread{&prometheus::seg_Indexer, this, start_Seg, stop_Seg, full_Char, VALID_or_NOT, pos_start_Index, pos_end_Index, start});
            }

            for (thread &t : threads_vec)
            {
                if (t.joinable())
                {
                    t.join();
                }
            }

            threads_vec.clear();

            // cout << position_index_Segs.size() << endl;
            // if ((rounds + 1) == 1)
            // {
            //     for (int c = 0; c < total_Segs+1; c++)
            //     {
            //         cout << site_Index[c] << endl;
            //     }
            // }

            free(full_Char);
            free(site_Index);
            free(VALID_or_NOT);
            free(pos_start_Index);
            free(pos_end_Index);

            // delete full_Char, site_Index, VALID_or_NOT, pos_start_Index, pos_end_Index;
        }

        cudaMemcpy(MA_Count, cuda_MA_Count, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        // USED in WINDOW MODE
        if (calc_Mode != "FILE")
        {
            pre_MA = (int *)malloc(tot_Segs * sizeof(int));
            cudaMemcpy(pre_MA, cuda_MA_Count, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        }

        // string Seg_sites = "";

        // int *site_Index;
        // site_Index = (int *)malloc((tot_Segs + 1) * sizeof(int));
        // site_Index[0] = 0;

        // // vector<string> total_Segregrating_sites = all_Lines;
        // for (size_t i = 0; i < all_Lines.size(); i++)
        // {
        //     Seg_sites.append(all_Lines[i]);
        //     site_Index[i + 1] = site_Index[i] + all_Lines[i].size();
        // }

        // char *full_Char;
        // full_Char = (char *)malloc((Seg_sites.size() + 1) * sizeof(char));
        // strcpy(full_Char, Seg_sites.c_str());
        // all_Lines.clear();

        // char *cuda_full_Char;
        // cudaMallocManaged(&cuda_full_Char, (Seg_sites.size() + 1) * sizeof(char));
        // int *cuda_site_Index;
        // cudaMallocManaged(&cuda_site_Index, (tot_Segs + 1) * sizeof(int));
        // cudaMemcpy(cuda_full_Char, full_Char, (Seg_sites.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
        // cudaMemcpy(cuda_site_Index, site_Index, (tot_Segs + 1) * sizeof(int), cudaMemcpyHostToDevice);

        // int *cuda_VALID_or_NOT, *VALID_or_NOT;
        // cudaMallocManaged(&cuda_VALID_or_NOT, tot_Segs * sizeof(int));
        // VALID_or_NOT = (int *)malloc(tot_Segs * sizeof(int));

        // int *cuda_MA_Count, *MA_Count;
        // cudaMallocManaged(&cuda_MA_Count, tot_Segs * sizeof(int));
        // MA_Count = (int *)malloc(tot_Segs * sizeof(int));

        // int *pos_start_Index, *pos_end_Index, *cuda_pos_start_Index, *cuda_pos_end_Index;
        // cudaMallocManaged(&cuda_pos_start_Index, tot_Segs * sizeof(int));
        // cudaMallocManaged(&cuda_pos_end_Index, tot_Segs * sizeof(int));
        // pos_start_Index = (int *)malloc(tot_Segs * sizeof(int));
        // pos_end_Index = (int *)malloc(tot_Segs * sizeof(int));

        // cuda_tajima_Prometheus<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_full_Char, cuda_site_Index, tot_Segs, cuda_VALID_or_NOT, cuda_MA_Count, cuda_pos_start_Index, cuda_pos_end_Index);
        // cudaDeviceSynchronize();

        // cudaMemcpy(VALID_or_NOT, cuda_VALID_or_NOT, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        // cudaMemcpy(MA_Count, cuda_MA_Count, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        // cudaMemcpy(pos_start_Index, cuda_pos_start_Index, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);
        // cudaMemcpy(pos_end_Index, cuda_pos_end_Index, tot_Segs * sizeof(int), cudaMemcpyDeviceToHost);

        // cudaFree(cuda_site_Index);
        cudaFree(cuda_MA_Count);
        // cudaFree(cuda_full_Char);

        // sort the index and positions of segs to global seg handler
        // cout << "System is indexing segregating site(s)" << endl;

        // int segs_per_Thread = tot_Segs / CPU_cores;
        // int remainder = tot_Segs % CPU_cores;

        // vector<thread> threads_vec;

        // for (int core_ID = 0; core_ID < CPU_cores; core_ID++)
        // {
        //     int start_Seg = core_ID * segs_per_Thread;
        //     int stop_Seg = start_Seg + segs_per_Thread;

        //     threads_vec.push_back(thread{&prometheus::seg_Indexer, this, start_Seg, stop_Seg, full_Char, VALID_or_NOT, pos_start_Index, pos_end_Index});
        // }

        // if (remainder != 0)
        // {
        //     int start_Seg = tot_Segs - remainder;
        //     int stop_Seg = tot_Segs;

        //     threads_vec.push_back(thread{&prometheus::seg_Indexer, this, start_Seg, stop_Seg, full_Char, VALID_or_NOT, pos_start_Index, pos_end_Index});
        // }

        // for (thread &t : threads_vec)
        // {
        //     if (t.joinable())
        //     {
        //         t.join();
        //     }
        // }

        // threads_vec.clear();

        // for (int seg_No = 0; seg_No < tot_Segs; seg_No++)
        // {
        //     if (VALID_or_NOT[seg_No] == 1)
        //     {
        //         string POS_string = "";

        //         for (int i = pos_start_Index[seg_No]; i < pos_end_Index[seg_No]; i++)
        //         {
        //             POS_string = POS_string + full_Char[i];
        //         }

        //         int POS = stoi(POS_string);

        //         position_index_Segs.push_back(make_pair(POS, seg_No));
        //     }
        // }

        // free(full_Char);

        // free(site_Index);

        sort(position_index_Segs.begin(), position_index_Segs.end());

        // TESTING purposes
        // cout << position_index_Segs[0].first << endl;
        // cout << all_Lines[position_index_Segs[0].second] << endl;
        // cout << position_index_Segs.size() << endl;
        // cout << endl;

        //  SEG_POS, array_pos/index
        // ENSURE CATCH POINT IS NOT -1 before submission
    }
    else
    {
        memcpy(MA_Count, pre_MA, tot_Segs * sizeof(int));
        // MA_Count = pre_MA;
    }

    cout << "Filteration and selection of segregating site(s)" << endl;

    if (this->sliding_Mode == "YES")
    {
        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            threads_vec.push_back(thread{&prometheus::get_POS_VALID, this, all_start_Co[gene_ID], gene_ID});
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();

        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            if (seg_catch_points_ALL[gene_ID] != -1)
            {
                threads_vec.push_back(thread{&prometheus::seg_Search_forward, this, gene_ID});
            }
            // else
            // {
            //     string write = "0\t0\tNA";
            //     write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;
            // }
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }
    else
    {

        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            if (catch_Point[gene_ID] != -1)
            {
                threads_vec.push_back(thread{&prometheus::seg_Search_catch_point, this, gene_ID});
            }
            else
            {
                // add to write_lines the blank
                string write = "0\t0\tNA";
                write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;
            }
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();

        for (int gene_ID = 0; gene_ID < gene_Size; gene_ID++)
        {
            if (catch_Point[gene_ID] != -1)
            {
                if (seg_catch_points_ALL[gene_ID] != -1)
                {
                    threads_vec.push_back(thread{&prometheus::seg_Search_forward, this, gene_ID});
                    threads_vec.push_back(thread{&prometheus::seg_Search_backward, this, gene_ID});
                }
                else
                {
                    // add to write_lines the blank
                    string write = "0\t0\tNA";
                    write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;
                }
            }
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }

    // process Tajima

    for (size_t gene_ID = 0; gene_ID < gene_Size; gene_ID++)
    {
        if ((catch_Point[gene_ID] != -1) && (seg_catch_points_ALL[gene_ID] != -1))
        {
            threads_vec.push_back(thread{&prometheus::calc_Tajima_Segs, this, gene_ID, MA_Count});
        }
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    free(MA_Count);
    cout << "System has processed all segregating site(s)" << endl;
}

void prometheus::get_POS_VALID(int start_Co, int gene_ID)
{
    int top = 0;
    int bottom = position_index_Segs.size() - 1;
    int middle = top + ((bottom - top) / 2);

    int seg_catch_point = -1;

    while (top <= bottom)
    {
        if (position_index_Segs[middle].first == start_Co)
        {
            seg_catch_point = middle;
            // found = 'Y';
            break;
        }
        else if (position_index_Segs[middle].first < start_Co)
        {
            top = middle + 1;
        }
        else
        {
            bottom = middle - 1;
        }
        middle = top + ((bottom - top) / 2);
    }

    unique_lock<shared_mutex> ul(g_mutex);
    seg_catch_points_ALL[gene_ID] = seg_catch_point;
    if (seg_catch_point != -1)
    {
        seg_catch_index_ALL[gene_ID] = position_index_Segs[seg_catch_point].second;
    }
}

void prometheus::calc_Tajima_Segs(int gene_ID, int *MA_Count)
{
    float tot_pairwise_Differences = 0;

    for (int index : seg_backward_index_ALL[gene_ID])
    {
        float MAF = (float)MA_Count[index] / (float)N;
        tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));
    }

    float MAF = (float)MA_Count[seg_catch_index_ALL[gene_ID]] / (float)N;
    tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));

    for (int index : seg_forward_index_ALL[gene_ID])
    {
        float MAF = (float)MA_Count[index] / (float)N;
        tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));
    }

    float pi = 0;
    float D = 0;
    string Tajima_D = "NA";

    int real_segregrating_Sites = seg_backward_index_ALL[gene_ID].size() + 1 + seg_forward_index_ALL[gene_ID].size();

    pi = (float)tot_pairwise_Differences / combinations;
    D = (float)(pi - (real_segregrating_Sites / an)) / sqrt(((e1 * real_segregrating_Sites) + (e2 * real_segregrating_Sites * (real_segregrating_Sites - 1))));
    Tajima_D = to_string(D);

    string write = to_string(pi) + "\t" + to_string(real_segregrating_Sites) + "\t" + Tajima_D;
    unique_lock<shared_mutex> ul(g_mutex);
    write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;

    //     // shared_lock<shared_mutex> sl(g_mutex);
    //    // int start = gene_Index[gene_ID];
    //     // int stop = gene_Index[gene_ID + 1];

    //     for (size_t i = start; i < stop; i++)
    //     {
    //         if (VALID_or_NOT[i] == 1)
    //         {
    //             real_segregrating_Sites = real_segregrating_Sites + 1;
    //             float MAF = (float)MA_Count[i] / (float)N;
    //             tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N, 2));
    //         }
    //     }

    //     // "Gene_name\tCoordinates\tPi\tS\tTajimas_D"
    //     float pi = 0;
    //     float D = 0;
    //     string Tajima_D;

    //     if (real_segregrating_Sites != 0)
    //     {
    //         pi = (float)tot_pairwise_Differences / combinations;
    //         D = (float)(pi - (real_segregrating_Sites / an)) / sqrt(((e1 * real_segregrating_Sites) + (e2 * real_segregrating_Sites * (real_segregrating_Sites - 1))));
    //         Tajima_D = to_string(D);
    //     }
    //     else
    //     {
    //         Tajima_D = "NA";
    //     }

    //     string write = to_string(pi) + "\t" + to_string(real_segregrating_Sites) + "\t" + Tajima_D;
    //     unique_lock<shared_mutex> ul(g_mutex);
    //     write_Lines[gene_ID] = write_Lines[gene_ID] + "\t" + write;
}

void prometheus::seg_Indexer(int start_Seg, int stop_Seg, char *full_Char, int *VALID_or_NOT, int *pos_start_Index, int *pos_end_Index, int start)
{

    // int start = core_ID * segs_per_Thread;
    // int stop = start + segs_per_Thread;

    vector<pair<int, int>> position_index_Segs_partial;

    for (int seg_No = start_Seg; seg_No < stop_Seg; seg_No++)
    {
        if (VALID_or_NOT[seg_No] == 1)
        {
            string POS_string = "";

            for (int i = pos_start_Index[seg_No]; i < pos_end_Index[seg_No]; i++)
            {
                POS_string = POS_string + full_Char[i];
            }

            int POS = stoi(POS_string);
            position_index_Segs_partial.push_back(make_pair(POS, (seg_No + start)));
        }
    }

    unique_lock<shared_mutex> ul(g_mutex);
    for (int i = 0; i < position_index_Segs_partial.size(); i++)
    {
        position_index_Segs.push_back(position_index_Segs_partial[i]);
    }
}

void prometheus::seg_Search_backward(int gene_ID)
{
    int pos = seg_catch_points_ALL[gene_ID] - 1;

    int start_Co = all_start_Co[gene_ID];
    int end_Co = all_end_Co[gene_ID];

    vector<int> seg_backward_index;

    while (pos >= 0)
    {
        int value_at_POS = position_index_Segs[pos].first;

        if (start_Co > value_at_POS)
        {
            break;
        }

        if ((value_at_POS >= start_Co) && (value_at_POS <= end_Co))
        {
            seg_backward_index.push_back(position_index_Segs[pos].second);
        }

        pos--;
    }

    unique_lock<shared_mutex> ul(g_mutex);
    seg_backward_index_ALL[gene_ID] = seg_backward_index;
}

void prometheus::seg_Search_forward(int gene_ID)
{
    int pos = seg_catch_points_ALL[gene_ID] + 1;

    int start_Co = all_start_Co[gene_ID];
    int end_Co = all_end_Co[gene_ID];

    vector<int> seg_forward_index;

    while (pos < position_index_Segs.size())
    {
        int value_at_POS = position_index_Segs[pos].first;

        if (end_Co < value_at_POS)
        {
            break;
        }

        if ((value_at_POS >= start_Co) && (value_at_POS <= end_Co))
        {
            seg_forward_index.push_back(position_index_Segs[pos].second);
        }

        pos++;
    }

    unique_lock<shared_mutex> ul(g_mutex);
    seg_forward_index_ALL[gene_ID] = seg_forward_index;
}

void prometheus::seg_Search_catch_point(int gene_ID)
{
    // binary search
    // SAVES time from interpolation since we cannot guarantee even distribution

    int top = 0;
    int bottom = position_index_Segs.size() - 1;
    int middle = top + ((bottom - top) / 2);

    int seg_catch_point = -1;

    int start_Co = all_start_Co[gene_ID];
    int end_Co = all_end_Co[gene_ID];

    while (top <= bottom)
    {
        if ((position_index_Segs[middle].first >= start_Co) && (position_index_Segs[middle].first <= end_Co))
        {
            seg_catch_point = middle;
            // found = 'Y';
            break;
        }
        else if (position_index_Segs[middle].first < start_Co)
        {
            top = middle + 1;
        }
        else
        {
            bottom = middle - 1;
        }
        middle = top + ((bottom - top) / 2);
    }

    // int start = 0;
    // int end = position_index_Segs.size() - 1;

    // int low_Value = position_index_Segs[0].first;
    // int high_Value = position_index_Segs[end].first;

    // int start_Co = all_start_Co[gene_ID];
    // int end_Co = all_end_Co[gene_ID];

    // int seg_catch_point = -1;

    // while (start <= end && start_Co >= low_Value && start_Co <= high_Value)
    // {
    //     // vector<string> line_Data_get;

    //     int pos = start + ((double)(end - start) / ((high_Value - low_Value)) * (start_Co - low_Value));

    //     int value_at_POS = position_index_Segs[pos].first;

    //     if ((value_at_POS >= start_Co) && (value_at_POS <= end_Co))
    //     {

    //         // catch point
    //         seg_catch_point = pos;
    //         // catch_file = folder_Index[pos].second;
    //         // file_List.push_back(folder_Index[pos].second);

    //         break;
    //     }
    //     else if (start_Co > value_at_POS)
    //     {
    //         start = pos + 1;
    //     }
    //     else
    //     {
    //         end = pos - 1;
    //     }

    //     low_Value = position_index_Segs[start].first;
    //     high_Value = position_index_Segs[end].first;
    // }

    unique_lock<shared_mutex> ul(g_mutex);
    seg_catch_points_ALL[gene_ID] = seg_catch_point;
    if (seg_catch_point != -1)
    {
        seg_catch_index_ALL[gene_ID] = position_index_Segs[seg_catch_point].second;
    }
}

void prometheus::get_Gene_info_and_catch(string gene_Combo, int gene_ID)
{
    // functions function = functions();
    vector<string> collect_Segregrating_sites;

    vector<string> split_Data;
    string line_trim = gene_Combo.substr(0, gene_Combo.find('\t'));
    string gene_Name = line_trim;

    // function.split(split_Data, gene_Combo, "\t");

    // cout << "Gene name\t: " << gene_Name << endl;
    vector<string> coordinates;

    string coordinates_String = gene_Combo.substr(gene_Combo.find('\t') + 1);
    while (coordinates_String.find(':') != string::npos)
    {
        coordinates.push_back(coordinates_String.substr(0, coordinates_String.find(':')));
        coordinates_String = coordinates_String.substr(coordinates_String.find(':') + 1);
    }

    coordinates.push_back(coordinates_String);

    // function.split(coordinates, split_Data[1], ":");
    int start_Co = stoi(coordinates[1]);
    int end_Co = stoi(coordinates[2]);

    // vector<string> file_List;
    // cout << endl;
    // cout << "System is retrieving file(s)" << endl;

    // file_List = compound_interpolationSearch(start_Co, end_Co);
    int catch_pos = -1;

    if (start_Co <= end_Co)
    {

        vector<string> line_Data;
        string query_Line = folder_Index[0].first;
        while (query_Line.find('_') != string::npos)
        {
            line_Data.push_back(query_Line.substr(0, query_Line.find('_')));
            query_Line = query_Line.substr(query_Line.find('_') + 1);
        }

        line_Data.push_back(query_Line);

        // function.split(line_Data, folder_Index[0].first, "_");

        int low_Value = stoi(line_Data[0]);
        line_Data.clear();

        query_Line = folder_Index[folder_Index.size() - 1].first;
        while (query_Line.find('_') != string::npos)
        {
            line_Data.push_back(query_Line.substr(0, query_Line.find('_')));
            query_Line = query_Line.substr(query_Line.find('_') + 1);
        }

        line_Data.push_back(query_Line);

        // function.split(line_Data, folder_Index[folder_Index.size() - 1].first, "_");
        int high_Value = stoi(line_Data[1]);
        line_Data.clear();

        int start = 0;
        int end = folder_Index.size() - 1;

        // string catch_file = "";

        while (start <= end && start_Co <= high_Value)
        {
            vector<string> line_Data_get;

            int pos = start + ((double)(end - start) / ((high_Value - low_Value)) * (start_Co - low_Value));

            query_Line = folder_Index[pos].first;
            while (query_Line.find('_') != string::npos)
            {
                line_Data_get.push_back(query_Line.substr(0, query_Line.find('_')));
                query_Line = query_Line.substr(query_Line.find('_') + 1);
            }

            line_Data_get.push_back(query_Line);

            // function.split(line_Data_get, folder_Index[pos].first, "_");
            int low_Value_atpos = stoi(line_Data_get[0]);
            int high_Value_atpos = stoi(line_Data_get[1]);

            if ((start_Co >= low_Value_atpos && start_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co >= high_Value_atpos))
            {

                // catch point
                catch_pos = pos;
                // catch_file = folder_Index[pos].second;
                // file_List.push_back(folder_Index[pos].second);

                break;
            }
            else if (start_Co > low_Value_atpos)
            {
                start = pos + 1;
            }
            else
            {
                end = pos - 1;
            }

            line_Data_get.clear();
            query_Line = folder_Index[start].first;
            while (query_Line.find('_') != string::npos)
            {
                line_Data_get.push_back(query_Line.substr(0, query_Line.find('_')));
                query_Line = query_Line.substr(query_Line.find('_') + 1);
            }

            line_Data_get.push_back(query_Line);

            // function.split(line_Data_get, folder_Index[start].first, "_");
            low_Value = stoi(line_Data_get[0]);

            line_Data_get.clear();
            query_Line = folder_Index[end].first;
            while (query_Line.find('_') != string::npos)
            {
                line_Data_get.push_back(query_Line.substr(0, query_Line.find('_')));
                query_Line = query_Line.substr(query_Line.find('_') + 1);
            }

            line_Data_get.push_back(query_Line);

            // function.split(line_Data_get, folder_Index[end].first, "_");
            high_Value = stoi(line_Data_get[1]);
        }

        // cout << "System has retrieved all file(s)" << endl;

        // cout << "System is collecting segregrating site(s)" << endl;

        // cout << file_List.size()<<endl;

        // cout << "Coordinates\t: Chromosome: " + coordinates[0] + " Start: " + to_string(start_Co) + " End: " + to_string(end_Co) + "\n";
    }
    unique_lock<shared_mutex> ul(g_mutex);
    // gi_mutex.lock();
    // all_Files[gene_ID] = file_List;
    // gi_mutex.unlock();

    // gi_mutex.lock();
    all_start_Co[gene_ID] = start_Co;
    // gi_mutex.unlock();

    // gi_mutex.lock();
    all_end_Co[gene_ID] = end_Co;
    // gi_mutex.unlock();

    // gi_mutex.lock();
    write_Lines[gene_ID] = gene_Name + "\t" + coordinates[0] + ":" + to_string(start_Co) + ":" + to_string(end_Co);

    catch_Point[gene_ID] = catch_pos;
    // catch_files[gene_ID] = catch_file;
    // gi_mutex.unlock();

    // cout << file_List.size() << endl;

    // vector<thread> threads_vec;

    // for (string files : file_List)
    // {
    //     threads_vec.push_back(thread{&prometheus::file_Reader, this, files, start_Co, end_Co, gene_ID});
    // }

    // for (thread &t : threads_vec)
    // {
    //     if (t.joinable())
    //     {
    //         t.join();
    //     }
    // }

    // // write_Lines[gene_ID] = gene_Name + "\t" + coordinates[0] + ":" + to_string(start_Co) + ":" + to_string(end_Co);
    // set_values(gene_ID, gene_Name, coordinates[0], to_string(start_Co), to_string(end_Co));
}

void prometheus::backward_Search(int pos, int start_Co, int end_Co, int gene_ID)
{
    // shared_lock<shared_mutex> sl(g_mutex);
    // functions function = functions();
    vector<int> backward_get;
    pos = pos - 1;

    while (pos >= 0)
    {
        vector<string> line_Data_get;
        string query_Line = folder_Index[pos].first;
        while (query_Line.find('_') != string::npos)
        {
            line_Data_get.push_back(query_Line.substr(0, query_Line.find('_')));
            query_Line = query_Line.substr(query_Line.find('_') + 1);
        }

        line_Data_get.push_back(query_Line);

        // function.split(line_Data_get, folder_Index[pos].first, "_");
        int low_Value_atpos = stoi(line_Data_get[0]);
        int high_Value_atpos = stoi(line_Data_get[1]);

        line_Data_get.clear();

        // cout << low_Value_atpos << endl;

        if (start_Co > high_Value_atpos)
        {
            break;
        }

        if ((start_Co >= low_Value_atpos && start_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co >= high_Value_atpos))
        {
            backward_get.push_back(pos);
            // backward_get.push_back(folder_Index[pos].second);
        }

        pos = pos - 1;
    }

    unique_lock<shared_mutex> ul(g_mutex);
    catch_back_index[gene_ID] = backward_get;
}

void prometheus::forward_Search(int pos, int start_Co, int end_Co, int gene_ID)
{
    // shared_lock<shared_mutex> sl(g_mutex);
    vector<int> forward_get;
    // functions function = functions();
    pos = pos + 1;

    while (pos < folder_Index.size())
    {
        vector<string> line_Data_get;
        string query_Line = folder_Index[pos].first;
        while (query_Line.find('_') != string::npos)
        {
            line_Data_get.push_back(query_Line.substr(0, query_Line.find('_')));
            query_Line = query_Line.substr(query_Line.find('_') + 1);
        }

        line_Data_get.push_back(query_Line);

        // function.split(line_Data_get, folder_Index[pos].first, "_");
        int low_Value_atpos = stoi(line_Data_get[0]);
        int high_Value_atpos = stoi(line_Data_get[1]);

        line_Data_get.clear();
        // cout << low_Value_atpos << endl;
        if (end_Co < low_Value_atpos)
        {
            break;
        }

        if ((start_Co >= low_Value_atpos && start_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co <= high_Value_atpos) || (start_Co <= low_Value_atpos && end_Co >= high_Value_atpos))
        {

            forward_get.push_back(pos);
            // forward_get.push_back(folder_Index[pos].second);
        }

        pos = pos + 1;
    }

    unique_lock<shared_mutex> ul(g_mutex);
    catch_forward_index[gene_ID] = forward_get;
}
