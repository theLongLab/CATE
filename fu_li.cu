#include "fu_li.cuh"
#include "functions.cuh"

fu_li::fu_li(string gene_List, string input_Folder, string ouput_Path, int cuda_ID, string intermediate_Path, int ploidy)
{
    cout << "Initiating CUDA powered Fu and Li's D, D*, F and F* calculator" << endl
         << endl;
    this->gene_List = gene_List;
    cout << "Gene list file path\t: " << gene_List << endl
         << endl;
    this->input_Folder = input_Folder;
    this->ouput_Path = ouput_Path;
    this->intermediate_Path = intermediate_Path;
    this->ploidy = ploidy;

    cudaSetDevice(cuda_ID);
    cout << "Properties of selected CUDA GPU:" << endl;
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, cuda_ID);
    cout << "GPU number\t: " << cuda_ID << endl;
    cout << "GPU name\t: " << prop.name << endl;
    size_t l_free = 0;
    size_t l_Total = 0;
    cudaError_t error_id = cudaMemGetInfo(&l_free, &l_Total);
    cout << "GPU memory (GB)\t: " << l_Total / (1000 * 1000 * 1000) << endl;
    cout << "GPU number of multiprocessor(s)\t: " << prop.multiProcessorCount << endl;
    cout << "GPU block(s) per multiprocessor\t: " << prop.maxBlocksPerMultiProcessor << endl;
    this->tot_Blocks = prop.maxBlocksPerMultiProcessor;
    this->tot_ThreadsperBlock = prop.maxThreadsPerBlock;
    cout << "GPU thread(s) per block\t: " << tot_ThreadsperBlock << endl
         << endl;
}

void fu_li::ingress()
{
    functions function = functions();

    vector<string> countries = function.get_Countries(this->input_Folder);
    cout << countries.size() << " populations were found: ";
    for (int count = 0; count < countries.size(); count++)
    {
        string folder = countries[count];
        cout << folder.substr(folder.find_last_of("/") + 1, folder.length());
        if (count < countries.size() - 1)
        {
            cout << ", ";
        }
    }
    cout << endl
         << endl;

    for (string country : countries)
    {
        cout << "Processing country\t: " << country.substr(country.find_last_of("/") + 1, country.length()) << endl
             << endl;

        vector<pair<string, string>> folder_Index = function.index_Folder(country);
        cout << "Completed indexing folder\t: " << country << endl;

        cout << endl;
        int samples = function.getN_Split(folder_Index[0].second);
        cout << "Number of samples in " << country.substr(country.find_last_of("/") + 1, country.length()) << " population\t: " << samples << endl;

        int N = samples * ploidy;
        float N_float = (float)N;
        cout << "Number of sequences in " << country.substr(country.find_last_of("/") + 1, country.length()) << " population [ " << samples << " x " << ploidy << " ] (N)\t: " << N << endl;
        long int combinations = function.combos_N(N);
        cout << "Pairwise combinations\t: " << combinations << endl;
        //float soft_Singl = 1 / N_float;
        //string SOFT_singleton_MAF = function.roundoff(soft_Singl, 4);
        //cout << "Singleton MAF\t: " << SOFT_singleton_MAF << endl;
        cout << endl;

        //calculate prerequisites
        float an, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star;
        calc_Pre(N, an, vd, ud, vd_star, ud_star, uf, vf, uf_star, vf_star);

        fstream gene_File;
        gene_File.open(gene_List, ios::in);
        cout << "Processing gene list:" << endl;
        string output_File = ouput_Path + "/" +
                             country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                             filesystem::path(gene_List).stem().string() +
                             ".fl";
        string intermediate_File = intermediate_Path + "/" +
                                   country.substr(country.find_last_of("/") + 1, country.length()) + "_" +
                                   filesystem::path(gene_List).stem().string() +
                                   ".log_fl";
        cout << endl;
        cout << "Writing to file\t: " << output_File << endl;
        cout << endl;

        if (gene_File.is_open())
        {
            string gene_Combo;

            if (filesystem::exists(output_File) == 0)
            {
                function.createFile(output_File, "Gene_name\tCoordinates\tPi\tS\tne\tns\tD\tD_star\tF\tF_star");
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

            fstream output;
            fstream intermediate;
            output.open(output_File, ios::app);
            intermediate.open(intermediate_File, ios::app);

            while (getline(gene_File, gene_Combo))
            {
                vector<string> split_Data;
                function.split(split_Data, gene_Combo, "\t");
                string gene_Name = split_Data[0];
                cout << "Gene name\t: " << gene_Name << endl;
                vector<string> coordinates;
                function.split(coordinates, split_Data[1], ":");
                int start_Co = stoi(coordinates[1]);
                int end_Co = stoi(coordinates[2]);
                cout << "Coordinates\t: Chromosome: " << coordinates[0] << " Start: " << start_Co << " End: " << end_Co << endl;

                float tot_pairwise_Differences = 0;
                int segregating_Sites = 0;
                int singletons_ns = 0;
                int singletons_ne = 0;

                vector<string> file_List;
                cout << endl;
                cout << "System is retrieving file(s)" << endl;
                if (folder_Index.size() > 1)
                {
                    file_List = function.compound_interpolationSearch(folder_Index, start_Co, end_Co);
                }
                else
                {
                    file_List.push_back(folder_Index[0].second);
                }
                cout << "System has retrieved all file(s)" << endl;
                cout << "System is collecting segregrating site(s)" << endl;
                vector<string> collect_Segregrating_sites;

                for (string files : file_List)
                {
                    //cout << files << endl;
                    fstream file;
                    file.open(files, ios::in);
                    if (file.is_open())
                    {
                        string line;
                        getline(file, line); //skip first header line
                        while (getline(file, line))
                        {
                            vector<string> positions;
                            function.split_getPos_ONLY(positions, line, "\t");
                            int pos = stoi(positions[1]);

                            if (pos >= start_Co && pos <= end_Co)
                            {
                                collect_Segregrating_sites.push_back(line);
                                // string check_0 = country.substr(country.find_last_of("/") + 1, country.length()) + "_AF=0";
                                // string seg_Check = "GO";
                                // vector<string> info;
                                // function.split(info, positions[7], ";");
                                // for (string AF_check : info)
                                // {
                                //     if (AF_check == check_0)
                                //     {
                                //         seg_Check = "NO";
                                //         break;
                                //     }
                                // }
                                // if (seg_Check == "GO")
                                // {
                                //     string check_AF_country = country.substr(country.find_last_of("/") + 1, country.length()) + "_AF";
                                //     float MAF_float = 0.0000;
                                //     segregating_Sites = segregating_Sites + 1;
                                //     for (string AF_check : info)
                                //     {
                                //         vector<string> split_info;
                                //         function.split(split_info, AF_check, "=");
                                //         if (split_info[0] == check_AF_country)
                                //         {
                                //             MAF_float = stof(split_info[1]);
                                //             string MAF = split_info[1];
                                //             string singleton_MAF;
                                //             string check;
                                //             if (MAF_float > 0.5)
                                //             {
                                //                 MAF_float = 1 - MAF_float;
                                //                 //singleton_MAF = function.roundoff(soft_Singl, MAF.length() - 2);
                                //                 string new_MAF = function.roundoff(MAF_float, MAF.length() - 2);
                                //                 check = check_AF_country + "=" + new_MAF;
                                //             }
                                //             else
                                //             {
                                //                 //MAF = split_info[1];
                                //                 check = AF_check;
                                //             }

                                //             singleton_MAF = function.roundoff(soft_Singl, MAF.length() - 2);
                                //             string check_singleton = check_AF_country + "=" + singleton_MAF;
                                //             if (check == check_singleton)
                                //             {
                                //                 //cout << "Singleton MAF check\t: " << check_singleton << "\t" << AF_check << endl;
                                //                 singletons_ns = singletons_ns + 1;
                                //                 singletons_ne = singletons_ne + outgroup_Singleton(info, positions);
                                //             }
                                //             break;
                                //         }
                                //     }
                                //     //int pairwise_Differences = function.calc_Pairwise(line, N, this->tot_Blocks, this->tot_ThreadsperBlock);
                                //     tot_pairwise_Differences = tot_pairwise_Differences + (MAF_float * (1 - MAF_float) * pow(N_float, 2));
                                //     //break;
                                // }
                            }
                            else if (pos > end_Co)
                            {
                                break;
                            }
                        }
                        file.close();
                    }
                }

                function.process_Seg_sites_fu_li(collect_Segregrating_sites, N_float, segregating_Sites, tot_pairwise_Differences, singletons_ne, singletons_ns, this->tot_Blocks, this->tot_ThreadsperBlock);

                cout << endl;

                //test
                // segregating_Sites = 18;
                // singletons_ne = 9;
                // singletons_ns = 10;
                //"Gene_name\tCoordinates\tPi\tS\tne\tns\tD\tD_star\tF\tF_star"
                cout << "Total segregating sites (S)\t: " << segregating_Sites << endl;

                string Fu_Li_D;
                string Fu_Li_D_star;
                string Fu_Li_F;
                string Fu_Li_F_star;
                float pi = 0;

                if (segregating_Sites != 0)
                {
                    //test
                    //N_float = 24.0;
                    float D = (float)(segregating_Sites - (an * singletons_ne)) / sqrt(((ud * segregating_Sites) + (vd * (pow(segregating_Sites, 2)))));
                    float D_star = (float)(((N_float / (N_float - 1)) * segregating_Sites) - (an * singletons_ns)) / sqrt(((ud_star * segregating_Sites) + (vd_star * (pow(segregating_Sites, 2.0)))));

                    pi = (float)tot_pairwise_Differences / combinations;
                    //test
                    //pi = 3.16;

                    cout << "Average pairwise polymorphisms (pi)\t: " << pi << endl;
                    cout << "Total ns singletons\t: " << singletons_ns << endl;
                    cout << "Total ne singletons\t: " << singletons_ne << endl;
                    cout << endl;

                    float F = (float)(pi - singletons_ne) / sqrt(((uf * segregating_Sites) + (vf * (pow(segregating_Sites, 2)))));
                    float F_star = (float)(pi - (((N_float - 1) / N_float) * singletons_ns)) / sqrt(((uf_star * segregating_Sites) + (vf_star * (pow(segregating_Sites, 2.0)))));

                    Fu_Li_D = to_string(D);
                    Fu_Li_D_star = to_string(D_star);
                    Fu_Li_F = to_string(F);
                    Fu_Li_F_star = to_string(F_star);

                    cout << "Fu and Li's D\t: " << Fu_Li_D << endl;
                    cout << "Fu and Li's D*\t: " << Fu_Li_D_star << endl;
                    cout << "Fu and Li's F\t: " << Fu_Li_F << endl;
                    cout << "Fu and Li's F*\t: " << Fu_Li_F_star << endl;
                }
                else
                {
                    cout << endl;
                    Fu_Li_D = "NA";
                    Fu_Li_D_star = "NA";
                    Fu_Li_F = "NA";
                    Fu_Li_F_star = "NA";
                    cout << "Fu and Li's D, D*, F and F*\t: "
                         << "Not Available" << endl;
                }

                cout << endl;
                //"Gene_name\tCoordinates\tPi\tS\tne\tns\tD\tD_star\tF\tF_star"
                output << gene_Name << "\t"
                       << coordinates[0] << ":" << to_string(start_Co) << ":" << to_string(end_Co)
                       << "\t" << to_string(pi)
                       << "\t" << to_string(segregating_Sites)

                       << "\t" << to_string(singletons_ne)
                       << "\t" << to_string(singletons_ns)

                       << "\t" << Fu_Li_D
                       << "\t" << Fu_Li_D_star
                       << "\t" << Fu_Li_F
                       << "\t" << Fu_Li_F_star << "\n";

                intermediate << gene_Combo << "\n";
                output.flush();
                intermediate.flush();
            }
            output.close();
            intermediate.close();
            gene_File.close();
        }
    }
}

int fu_li::outgroup_Singleton(vector<string> &info, vector<string> &positions)
{
    functions function = functions();
    string present = "false";
    string MA;
    string AA;
    int return_Val = 0;
    vector<string> line_partial;

    for (string info_Data : info)
    {
        vector<string> AA_check;
        function.split(AA_check, info_Data, "=");
        if (AA_check[0] == "AA")
        {
            //cout << "AA: " << AA_check[1].at(0) << endl;
            if ((toupper(AA_check[1].at(0)) == 'A') || (toupper(AA_check[1].at(0)) == 'T') || (toupper(AA_check[1].at(0)) == 'G') || (toupper(AA_check[1].at(0)) == 'C'))
            {
                //cout << "AA Caught: " << AA_check[1].at(0) << endl;
                AA = AA_check[1].at(0);
                present = "true";
            }
            break;
        }
    }

    if (present == "false")
    {
        AA = positions[3];
    }

    int REF_0 = 0;
    int ALT_1 = 0;

    vector<string> sample_1;
    vector<string> sample_2;

    function.split(sample_1, positions[9], "|");
    function.split(sample_2, positions[10], "|");

    for (string check : sample_1)
    {
        if (check == "0")
        {
            REF_0 = REF_0 + 1;
        }
        else
        {
            ALT_1 = ALT_1 + 1;
        }
    }

    for (string check : sample_2)
    {
        if (check == "0")
        {
            REF_0 = REF_0 + 1;
        }
        else
        {
            ALT_1 = ALT_1 + 1;
        }
    }

    if (REF_0 > ALT_1)
    {
        MA = positions[4];
    }
    else
    {
        MA = positions[3];
    }

    if (toupper(AA.at(0)) == toupper(MA.at(0)))
    {
        return_Val = 0;
    }
    else
    {
        return_Val = 1;
    }

    return return_Val;
}

__global__ void fuli_Calculation(int N, float *a1_CUDA, float *a2_CUDA)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < N)
    {
        a1_CUDA[tid] = (float)1 / (tid + 1);
        a2_CUDA[tid] = (float)1 / ((tid + 1) * (tid + 1));
        tid += blockDim.x * gridDim.x;
    }
}

void fu_li::calc_Pre(int N_tot, float &an, float &vd, float &ud, float &vd_star, float &ud_star, float &uf, float &vf, float &uf_star, float &vf_star)
{
    //functions function = functions();
    //test
    //N_tot = 24;
    int N = N_tot - 1;
    float *a1_CUDA, *a2_CUDA;
    float *a1_partial, *a2_partial, bn;

    a1_partial = (float *)malloc(N * sizeof(float));
    a2_partial = (float *)malloc(N * sizeof(float));

    cudaMallocManaged(&a1_CUDA, N * sizeof(int));
    cudaMallocManaged(&a2_CUDA, N * sizeof(int));

    fuli_Calculation<<<tot_Blocks, tot_ThreadsperBlock>>>(N, a1_CUDA, a2_CUDA);
    cudaDeviceSynchronize();

    cudaMemcpy(a1_partial, a1_CUDA, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(a2_partial, a2_CUDA, N * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(a1_CUDA);
    cudaFree(a2_CUDA);

    // an = function.add(N, a1_partial);
    // bn = function.add(N, a2_partial);

    an = 0;
    bn = 0;
    for (size_t i = 0; i < N; i++)
    {
        an += a1_partial[i];
        bn += a2_partial[i];
    }

    free(a1_partial);
    free(a2_partial);

    float N_float_tot = (float)N_tot;
    float N_float = (float)N;

    float an_plus_1 = (float)an + (1 / N_float_tot);
    //cout << "an+1: " << an_plus_1 << endl;
    float cn = (float)2 * (((N_float_tot * an) - (2 * N_float)) / (N_float * (N_float_tot - 2)));
    float an_square = (float)pow(an, 2.0);

    //float test = (N_float_tot_float + 1) / N_float_float;
    //cout << "test: " << test << endl;

    vd = (float)1 + ((an_square / (bn + an_square)) * (cn - (N_float_tot + 1) / N_float));
    ud = (float)an - 1 - vd;

    float dn = (float)cn + ((N_float_tot - 2) / pow(N_float, 2.0)) + ((2 / N_float) * (1.5 - (((2 * an_plus_1) - 3) / (N_float_tot - 2)) - (1 / N_float_tot)));
    //cout << "dn: " << dn << endl;
    vd_star = (float)((pow((N_float_tot / N_float), 2.0) * bn) + (pow(an, 2.0) * dn) - (2 * ((N_float_tot * an * (an + 1)) / (pow(N_float, 2.0))))) / (pow(an, 2.0) + bn);
    ud_star = (float)N_float_tot / N_float * (an - (N_float_tot / N_float)) - vd_star;

    vf = (float)(cn + ((2 * (pow(N_float_tot, 2) + N_float_tot + 3)) / (9 * N_float_tot * N_float)) - (2 / N_float)) / (pow(an, 2) + bn);
    uf = (float)((1 + ((N_float_tot + 1) / (3 * N_float)) - ((4 * ((N_float_tot + 1) / pow(N_float, 2))) * (an_plus_1 - ((2 * N_float_tot) / (N_float_tot + 1))))) / an) - vf;

    //Original papers incorrect equations due to paper's typos
    //vf_star = (float)(dn + ((2 * (pow(N_float_tot, 2) + N_float_tot + 3)) / (9 * N_float_tot * N_float)) - ((2 * (1 / N_float)) * ((4 * bn) - 6 + (8 / N_float_tot)))) / (pow(an, 2) + bn);
    //uf_star = (float)(((N_float_tot / N_float) + ((N_float_tot + 1) / (3 * N_float)) - (2 * (2 / (N_float_tot * N_float))) + ((2 * ((N_float_tot + 1) / pow(N_float, 2))) * (an_plus_1 - ((2 * N_float_tot) / (N_float_tot + 1))))) / an) - vf_star;

    //corrected equations Simonsen et al 1995
    vf_star = (float)((((2 * pow(N_float_tot, 3.0)) + (110 * pow(N_float_tot, 2.0)) - (225 * N_float_tot) + 153) / (9 * pow(N_float_tot, 2.0) * N_float)) + ((2 * N_float * an) / (pow(N_float_tot, 2.0))) - ((8 * bn) / N_float_tot)) / (pow(an, 2.0) + bn);
    uf_star = (float)((((4 * pow(N_float_tot, 2.0)) + (19 * N_float_tot) + 3 - (12 * (N_float_tot + 1) * an_plus_1)) / (3 * N_float_tot * (N_float))) / an) - vf_star;

    cout << "an: " << an << "\t"
         << "bn: " << bn << "\t"
         << "cn: " << cn << endl;

    cout << endl;

    cout << "ud\t: " << ud << "\t"
         << "vd\t: " << vd << endl;

    cout << "ud*\t: " << ud_star << "\t"
         << "vd*\t: " << vd_star << endl;

    cout << endl;

    cout << "an+1\t: " << an_plus_1 << "\t"
         << "dn\t: " << dn << endl;

    cout << endl;

    cout << "uf\t: " << uf << "\t"
         << "vf\t: " << vf << endl;

    cout << "uf*\t: " << uf_star << "\t"
         << "vf*\t: " << vf_star << endl;

    cout << endl;
}