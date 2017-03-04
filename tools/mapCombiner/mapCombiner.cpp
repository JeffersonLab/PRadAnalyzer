#include <string>
#include <vector>
#include <queue>
#include <list>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include "ConfigParser.h"

using namespace std;

int main(int argc, char *argv[])
{
    char *ptr;
    string output = "combined_map.txt", input;
    queue<string> input_list;

    for(int i = 1; i < argc; ++i)
    {
        ptr = argv[i];
        if(*(ptr++) == '-') {
            switch(*(ptr++))
            {
            case 'o':
                output = argv[++i];
                break;
            case 'i':
                for(int j = i + 1; j < argc; ++j)
                {
                    if(*argv[j] == '-') {
                        i = j-1;
                        break;
                    }
                    input = argv[j];
                    if(!input.empty())
                        input_list.push(input);
                }
                break;
            default:
                printf("Unkown option!\n");
                exit(1);
            }
        }
    }

    ConfigParser parser;
    unordered_map<string, list<ConfigValue>> eles_map;

    if(input_list.size() < 2) {
        cerr << "Less than 2 input files, no need to combine" << endl;
    }

    bool first_map = true;
    while(input_list.size())
    {
        int line = 0;
        parser.ReadFile(input_list.front());
        input_list.pop();
        while(parser.ParseLine())
        {
            if(!parser.NbofElements())
                continue;

            ++line;
            string key;
            if(parser.NbofElements() == 1) { // PrimEx map, line number = primex id
                if(line < 1000)
                    key = "G" + to_string(line);
                else
                    key = "W" + to_string(line - 1000);
            } else {
                key = parser.TakeFirst();
            }

            list<ConfigValue> eles = parser.TakeAll<list>();

            auto it = eles_map.find(key);
            if(it == eles_map.end()) {
                if(!first_map) { // do not create map element if it is not the first map
                    cout << "Cannot find key " << key << " in original map" << endl;
                    continue;
                }
                eles_map[key] = eles;
            } else {
                it->second.splice(it->second.end(), eles);
            }
        }
        first_map = false;
    }

    size_t nb_ele = 0;
    ofstream outfile(output);
    for(auto ele : eles_map)
    {
        // sanity check
        if(nb_ele && nb_ele != ele.second.size())
        {
            cout << "Skipped " << ele.first
                 << ", it has " << ele.second.size()
                 << " elements, while the previous one has " << nb_ele
                 << ", one of the map might miss this element."
                 << endl;
            continue;
        }

        nb_ele = ele.second.size();

        outfile << setw(6) << ele.first;
        list<ConfigValue> eles = ele.second;
        for(auto e : eles)
        {
            outfile << setw(13) << e;
        }
        outfile << endl;
    }

    outfile.close();
    return 0;
}
