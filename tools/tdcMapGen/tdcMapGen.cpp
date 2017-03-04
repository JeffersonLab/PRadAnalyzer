#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <unordered_map>

using namespace std;

struct DSC_ADDR
{
    int id;
    int channel;
    DSC_ADDR() : id(0), channel(0) {};
    DSC_ADDR(const int &i, const int &c): id(i), channel(c) {};
    bool operator==(const DSC_ADDR &rhs) const
    {
        return ((id == rhs.id)&&(channel == rhs.channel));
    }
};

namespace std {
    template <>
    struct hash<DSC_ADDR>
    {
        size_t operator()(const DSC_ADDR& addr) const
        {
            return (addr.id << 24 | addr.channel);
        }
    };
}

string trim(const string& str,
            const string& whitespace = " \t")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

string comment_out(const string &str,
                   const string &comment_mark = "#")
{
    if(comment_mark.empty()) {
        return str;
    }

    const auto cmtEnd = str.find(comment_mark);
    return str.substr(0, cmtEnd);
}

int main()
{
    int crate, slot;
    bool read_crate = false, read_slot = false;
    unordered_map <DSC_ADDR, int> tdc_map;

    ifstream infile("tdc_map.txt");
    string line, trim_line;

    // read tdc map
    while(getline(infile, line))
    {
        trim_line = trim(comment_out(line));
        if(trim_line.empty())
            continue;
        istringstream iss(trim_line);
        if(read_crate == false) {
            iss >> crate;
            read_crate = true;
            cout << "read tdc crate: " << crate << endl;
        } else if(read_slot == false) {
            iss >> slot;
            read_slot = true;
            cout << "read tdc slot: " << slot << endl;
        } else {
            string tdc_module;
            int dsc_id, dsc_channel, channel;
            iss >> tdc_module >> dsc_id >> dsc_channel;
            channel = stoi(tdc_module.substr(1));
            if(tdc_module.at(0) >= 'A' && tdc_module.at(0) <= 'D') {
                channel += 32*int(tdc_module.at(0) - 'A');
            } else {
                cout << "unrecognized module: " << tdc_module << endl;
                continue;
            }
            tdc_map[DSC_ADDR(dsc_id, dsc_channel)] = channel;
        }
    }
    cout << "read tdc map, " << tdc_map.size() << " channels in total." << endl;
    infile.close();

    infile.open("dsc_map.txt");
    ofstream outfile("tdc_group_list.txt");
    outfile << "#"
            << setw(11) << "name"
            << setw(12) << "crate"
            << setw(12) << "slot"
            << setw(12) << "channel"
            << endl;

    int count = 0;
    while(getline(infile, line))
    {
        trim_line = trim(line);
        if(trim_line.at(0) == '#' || trim_line.empty())
            continue;
        istringstream iss(trim_line);
        string tdc_name;
        int dsc_id, dsc_channel;
        iss >> dsc_id >> dsc_channel >> tdc_name;
        if(tdc_name == "NONE") 
            continue;
        const auto &it = tdc_map.find(DSC_ADDR(dsc_id, dsc_channel));
        if(it == tdc_map.end()) {
            cout << "cannot find tdc channel for discriminator id " << dsc_id << " channel " << dsc_channel << endl;
            continue;
        } else {
            ++ count;
            outfile << setw(12) << tdc_name
                    << setw(12) << crate
                    << setw(12) << slot
                    << setw(12) << it->second
                    << endl;
        }
    }
    cout << "read dsc map, found " << count << " matched with tdc map." << endl;
    cout << "tdc group list is written in tdc_group_list.txt" << endl;
    outfile.close();
    infile.close();
    return 0;
}
