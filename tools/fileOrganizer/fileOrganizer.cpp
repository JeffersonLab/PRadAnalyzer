#include <string>
#include <vector>
#include <queue>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include "ConfigParser.h"

using namespace std;

int main(int argc, char *argv[])
{
    char *ptr;
    string order, input;

    for(int i = 1; i < argc; ++i)
    {
        ptr = argv[i];
        if(*(ptr++) == '-') {
            switch(*(ptr++))
            {
            case 'o':
                order = argv[++i];
                break;
            case 'i':
                input = argv[++i];
                break;
            default:
                printf("Unkown option!\n");
                exit(1);
            }
        }
    }

    ConfigParser parser;
    auto order_list = ConfigParser::split(order, ",");

    vector<int> orders;

    while(order_list.size())
    {
        orders.push_back(stoi(order_list.front()));
        order_list.pop();
    }

    int line = 0;
    if(!parser.ReadFile(input)) {
        cout << "Cannot open file " << input
             << endl;
        return -1;
    }

    ofstream out("re_ordered.txt");

    while(parser.ParseLine())
    {
        if(!parser.NbofElements())
            continue;

        vector<ConfigValue> eles = parser.TakeAll<vector>();

        for(auto &ord : orders)
            out << setw(12) << eles.at(ord);

        out << endl;
    }

    out.close();
    return 0;
}
