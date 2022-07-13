#include "functions.cuh"
#include "fst_test_pop.cuh"

fst_test_pop::fst_test_pop() {}

fst_test_pop::fst_test_pop(string Pop_ID, vector<string> super_Pops, vector<string> sample_IDs, vector<string> countries_Folder, int ploidy)
{
    this->Pop_ID = Pop_ID;
    cout << "Intializing: " << this->Pop_ID << endl
         << endl;
    this->super_Pops = super_Pops;
    this->sample_IDs = sample_IDs;
    this->countries_Folder = countries_Folder;
    this->ploidy = ploidy;

    index_Population();
    index_Samples();
}

void fst_test_pop::clear_final_Seg_Collection()
{
    final_Seg_Collection.clear();
}

vector<pair<int, string>> fst_test_pop::get_final_Seg_Collection()
{
    return final_Seg_Collection;
}

void fst_test_pop::clear_sample_Location()
{
    sample_Location.clear();
}

vector<int> fst_test_pop::get_sample_Location()
{
    return sample_Location;
}

void fst_test_pop::combine_Segs()
{
    if (super_Pops.size() > 1)
    {
        cout << "\t : ";
        vector<pair<int, string>> Super_pop_One = collect_Segregrating_sites_All[0];
        // cout << "TEST " << collect_Segregrating_sites_All.size();
        sort(Super_pop_One.begin(), Super_pop_One.end());

        vector<int> pos_Found_index;
        vector<pair<int, string>> concat_Sites;

        for (size_t i = 0; i < Super_pop_One.size(); i++)
        {
            pos_Found_index.push_back(1);
        }

        cout << "Merging and concatanating segregating sites: " << super_Pops[0] << ", ";
        for (size_t i = 1; i < super_Pops.size(); i++)
        {
            cout << super_Pops[i];

            vector<pair<int, string>> query_Super_Pop = collect_Segregrating_sites_All[i];
            sort(query_Super_Pop.begin(), query_Super_Pop.end());

            for (size_t j = 0; j < Super_pop_One.size(); j++)
            {
                // binary search
                int top = 0;
                int bottom = query_Super_Pop.size() - 1;
                int middle = top + ((bottom - top) / 2);

                while (top <= bottom)
                {
                    if (query_Super_Pop[middle].first == Super_pop_One[j].first)
                    {
                        int val = pos_Found_index[j];
                        val = val + 1;
                        pos_Found_index[j] = val;

                        // concat
                        string one_Line = Super_pop_One[j].second;
                        string query_Line = query_Super_Pop[middle].second;

                        for (int t = 0; t < 9; t++)
                        {
                            int tenth_Column = query_Line.find('\t') + 1;
                            query_Line = query_Line.substr(tenth_Column, query_Line.length());
                        }

                        one_Line = one_Line + query_Line;
                        Super_pop_One[j].second = one_Line;

                        break;
                    }
                    else if (query_Super_Pop[middle].first < Super_pop_One[j].first)
                    {
                        top = middle + 1;
                    }
                    else
                    {
                        bottom = middle - 1;
                    }
                    middle = top + ((bottom - top) / 2);
                }
            }

            if (i != (super_Pops.size() - 1))
            {
                cout << ", ";
            }
        }

        // cout << endl;

        // fstream test;
        // test.open("test_concat", ios::out);

        for (size_t i = 0; i < Super_pop_One.size(); i++)
        {
            if (pos_Found_index[i] == super_Pops.size())
            {
                // cout << pos_Found_index[i];
                // test << Super_pop_One[i].second << "\n";
                final_Seg_Collection.push_back(make_pair(Super_pop_One[i].first, Super_pop_One[i].second));
            }
        }
        // test.close();
        // for (int test : pos_Found_index)
        // {
        //     if (test == super_Pops.size())
        //     {
        //         cout << test << endl;
        //     }
        // }

        // clear collect_Segregrating_sites_All after this
    }
    else
    {
        final_Seg_Collection = collect_Segregrating_sites_All[0];
    }
    collect_Segregrating_sites_All.clear();
    sort(final_Seg_Collection.begin(), final_Seg_Collection.end());
    cout << endl;
}

vector<pair<int, string>> fst_test_pop::return_Seg_site(string super_Pop)
{
    vector<pair<int, string>> retrieve;
    for (size_t i = 0; i < super_Pops.size(); i++)
    {
        if (super_Pops[i] == super_Pop)
        {
            cout << "Bringing " << super_Pop << " forward" << endl;
            retrieve = collect_Segregrating_sites_All[i];
            break;
        }
    }
    return retrieve;
}

void fst_test_pop::seg_Retrival_with_Found(int start_Co, int end_Co, vector<string> super_Pops_FOUND, vector<vector<pair<int, string>>> collect_Segregrating_sites_FOUND)
{
    // We do not clear cause we need the previous ones
    // collect_Segregrating_sites_All.clear();
    functions function = functions();

    for (int i = 0; i < super_Pops.size(); i++)
    {
        cout << "Collecting seg sites: " << super_Pops[i] << endl;
        int found = 0;
        vector<pair<int, string>> collect_Segregrating_sites;

        for (size_t check = 0; check < super_Pops_FOUND.size(); check++)
        {
            if (super_Pops[i] == super_Pops_FOUND[check])
            {
                collect_Segregrating_sites = collect_Segregrating_sites_FOUND[check];
                cout << "Carried " << super_Pops[i] << " forward\n";
                collect_Segregrating_sites_All.push_back(collect_Segregrating_sites);
                found = 1;
                break;
            }
        }

        if (found == 0)
        {
            vector<string> file_List = file_List_All[i];
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
                        vector<string> positions;
                        function.split_getPos_ONLY(positions, line, '\t');
                        int pos = stoi(positions[1]);

                        if (pos >= start_Co && pos <= end_Co)
                        {
                            collect_Segregrating_sites.push_back(make_pair(pos, line));
                        }
                        else if (pos > end_Co)
                        {
                            break;
                        }
                    }
                    file.close();
                }
            }
            collect_Segregrating_sites_All.push_back(collect_Segregrating_sites);
        }
    }
    cout << "System has finished collecting segregrating site(s) for " << this->Pop_ID << endl;
}

void fst_test_pop::seg_Retrival(int start_Co, int end_Co)
{
    // Since this is the FIRST run we clear the list
    collect_Segregrating_sites_All.clear();
    functions function = functions();

    for (int i = 0; i < super_Pops.size(); i++)
    {
        cout << "Collecting seg sites: " << super_Pops[i] << endl;
        vector<string> file_List = file_List_All[i];
        vector<pair<int, string>> collect_Segregrating_sites;

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
                    vector<string> positions;
                    function.split_getPos_ONLY(positions, line, '\t');
                    int pos = stoi(positions[1]);

                    if (pos >= start_Co && pos <= end_Co)
                    {
                        collect_Segregrating_sites.push_back(make_pair(pos, line));
                    }
                    else if (pos > end_Co)
                    {
                        break;
                    }
                }
                file.close();
            }
        }
        collect_Segregrating_sites_All.push_back(collect_Segregrating_sites);
    }
    cout << "System has finished collecting segregrating site(s) for " << this->Pop_ID << endl;
}

void fst_test_pop::folder_Search(int start_Co, int end_Co)
{
    cout << "File retrival for: " << this->Pop_ID << endl;

    functions function = functions();

    this->file_List_All.clear();

    for (int i = 0; i < super_Pops.size(); i++)
    {
        cout << "System is retrieving file(s) for " << super_Pops[i] << endl;
        vector<string> file_List;
        vector<pair<string, string>> folder_Index = this->folder_Index_Super_Pops[i];

        if (folder_Index.size() > 1)
        {
            file_List = function.compound_interpolationSearch(folder_Index, start_Co, end_Co);
        }
        else
        {
            file_List.push_back(folder_Index[0].second);
        }

        file_List_All.push_back(file_List);
        cout << "System has retrieved all file(s) for " << super_Pops[i] << endl;
    }
}

void fst_test_pop::index_Samples()
{
    functions function = functions();

    string complete_Line = "";

    cout << "Mapping sample IDs:" << endl;
    for (size_t i = 0; i < super_Pops.size(); i++)
    {
        vector<pair<string, string>> folder_Index = folder_Index_Super_Pops[i];

        string first_file = folder_Index[0].second;
        fstream super_Pop_file;
        super_Pop_file.open(first_file, ios::in);
        string line;

        // get first line
        if (super_Pop_file.is_open())
        {
            getline(super_Pop_file, line);
            super_Pop_file.close();
        }

        if (i > 0)
        {
            for (int t = 0; t < 9; t++)
            {
                int tenth_Column = line.find('\t') + 1;
                line = line.substr(tenth_Column, line.length());
            }

            // cout << line << endl;
        }

        complete_Line = complete_Line + line;
    }

    // cout << complete_Line << endl;
    //  fstream test;
    //  test.open("test.csv", ios::out);
    //  test << complete_Line;
    //  test.close();

    vector<string> split_Header;
    function.split(split_Header, complete_Line, '\t');

    for (size_t i = 9; i < split_Header.size(); i++)
    {
        string query_ID = split_Header[i];

        for (string sample_ID : this->sample_IDs)
        {
            if (query_ID == sample_ID)
            {
                // cout << sample_ID << endl;
                this->sample_Location.push_back(i);
                break;
            }
        }
    }

    this->samples = sample_Location.size();
    int N = (ploidy * samples);
    this->N = N;
    this->N_float = (float)this->N;
    this->combinations = function.combos_N(this->N);

    cout << "Number of samples in " << this->Pop_ID << " population: " << samples << endl;
    cout << "Number of sequences in " << this->Pop_ID << " population [ " << samples << " x " << ploidy << " ] (N): " << N << endl;
    cout << "Pairwise combinations: " << combinations << endl;
}

int fst_test_pop::get_Sample_Size()
{
    return this->samples;
}

int fst_test_pop::get_Sequence_Size()
{
    return this->N;
}

void fst_test_pop::index_Population()
{
    functions function = functions();
    for (string super_Pop : this->super_Pops)
    {
        cout << "Processing super population: " << super_Pop << endl;
        string folder_Path;

        for (string folder : this->countries_Folder)
        {
            string country = folder.substr(folder.find_last_of("/") + 1, folder.length());
            if (country == super_Pop)
            {
                folder_Path = folder;
                break;
            }
        }

        vector<pair<string, string>> folder_Index = function.index_Folder(folder_Path);
        folder_Index_Super_Pops.push_back(folder_Index);
        cout << "Completed indexing folder\t: " << folder_Path << endl;
        cout << endl;
    }
}
