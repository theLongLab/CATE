#include <iostream>
#include <filesystem>
#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>

using namespace std;

class parameter{

    private:
    string path;

    public:
    parameter(string &path);
    string where(string word);
    int where_Int(string word);
    
    void split(vector<string> &line_Data, string &line);
};