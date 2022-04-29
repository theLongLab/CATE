#include "parameter.h"

parameter::parameter(string &path)
{
    this->path = path;
}

int parameter::where_Int(string word)
{
    string parametric = where(word);
    return stoi(parametric);
}

string parameter::where(string word)
{
    string paramteric = "NA";
    string quotes = "\"" + word + "\"";

    fstream parameter_File;
    parameter_File.open(path, ios::in);
    if (parameter_File.is_open())
    {
        string line;
        vector<string> line_Data;
        getline(parameter_File, line);
        while (getline(parameter_File, line))
        {
            if (line.find("\"") > line.length() && line.find("#") < line.length())
            {
                // cout << line<< endl;
            }
            else
            {
                line.erase(remove(line.begin(), line.end(), '{'), line.end());
                line.erase(remove(line.begin(), line.end(), '}'), line.end());
                if (line.length() > 0)
                {
                    split(line_Data, line);
                    if (line_Data.size() > 0)
                    {

                        if (line_Data[0] == quotes)
                        {
                            line_Data[1];

                            if (line_Data[1].find_last_of(',') == line_Data[1].length() - 1)
                            {
                                // line_Data[1].erase(remove(line_Data[1].begin(), line_Data[1].end(), ','), line_Data[1].end());
                                line_Data[1] = line_Data[1].substr(0, line_Data[1].length() - 1);
                            }

                            if (line_Data[1].find("\"") < line_Data[1].length())
                            {
                                line_Data[1].erase(remove(line_Data[1].begin(), line_Data[1].end(), '\"'), line_Data[1].end());
                            }
                            paramteric = line_Data[1];
                            break;
                        }
                    }
                }
            }
        }
        parameter_File.close();
    }
    if (paramteric != "NA")
    {
        return paramteric;
    }
    else
    {
        cout << "SPECIFIED PARAMETER WAS NOT FOUND IN THE JSON FILE. PLEASE CHECK THE FOLLOWING SECTION\t: " << word << endl;
        cout << "PATH OF ERRORED PARAMTER JSON FILE\t: " << path << endl;
        exit(4);
    }
}

void parameter::split(vector<string> &line_Data, string &line)
{
    vector<string>().swap(line_Data);
    char *convert;
    string capture(line);
    convert = &capture[0];
    // cout<<convert;

    char *split_data;
    split_data = strtok(convert, ":");

    while (split_data != NULL)
    {
        // cout<<split_data<<endl;
        string char2string;
        char2string.append(split_data);
        string complete;
        if (char2string.find("\"") < char2string.length())
        {
            string begin = char2string.substr(0, char2string.find("\""));
            begin.erase(remove(begin.begin(), begin.end(), ' '), begin.end());
            complete = begin + char2string.substr(char2string.find("\""), char2string.length());
        }
        else
        {
            complete = char2string;
        }
        // cout << char2string << endl;
        line_Data.push_back(complete);
        split_data = strtok(NULL, "\t");
    }
}
