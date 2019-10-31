//============================================================================//
// An example showing how to train the neural network for cosmic rejection.   //
//                                                                            //
// Chao Peng                                                                  //
// 01/27/2017                                                                 //
//============================================================================//

#include "CNeuralNetwork.h"
#include "cosmicEval.h"
#include "PRadHyCalSystem.h"
#include "PRadDataHandler.h"
#include "PRadDSTParser.h"
#include "PRadBenchMark.h"
#include "ConfigOption.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <vector>

#define PROGRESS_COUNT 1000
#define NN_INPUT 6
#define NN_OUTPUT 1

using namespace std;

void NeuralTrain(CNeuralNetwork &net,
                 PRadHyCalSystem &sys,
                 const string &path,
                 const string &path2,
                 int train_times,
                 unsigned int cap);

void FillParams(PRadHyCalSystem &sys,
                const string &path,
                vector<vector<double>> &params,
                unsigned int cap);

int main(int argc, char *argv[])
{
    string net_path, save_path, layer_str;
    double learn_factor = 0.1;
    int learn_times = 5000, cap = 1000;

    save_path = "saved.net";

    ConfigOption conf_opt;
    conf_opt.AddOpts(ConfigOption::arg_require, 'n', "net-path");
    conf_opt.AddOpts(ConfigOption::arg_require, 's', "save-path");
    conf_opt.AddOpts(ConfigOption::arg_require, 'l', "layer");
    conf_opt.AddOpts(ConfigOption::arg_require, 'f', "learning-factor");
    conf_opt.AddOpts(ConfigOption::arg_require, 't', "training-times");
    conf_opt.AddOpts(ConfigOption::arg_require, 'c', "bank-capacity");
    conf_opt.AddOpts(ConfigOption::help_message, 'h', "help");

    conf_opt.SetDesc("usage: %0 <cosmic_data> <good_data>");
    conf_opt.SetDesc('n', "create network from file, a new network will be created by default.");
    conf_opt.SetDesc('l', "set the hidden layers for neural network in the format of \"val1, val2, ...\".");
    conf_opt.SetDesc('s', "set the path to save the trained network, save to \"saved.net\" by default.");
    conf_opt.SetDesc('f', "define learning factor, 0.1 is the default value.");
    conf_opt.SetDesc('t', "set the training times (1,000 as the unit), 5,000k is the default value.");
    conf_opt.SetDesc('c', "set the training bank capacity (1,000 as the unit), 1,000k is the default value.");
    conf_opt.SetDesc('h', "show instruction.");

    if(!conf_opt.ParseArgs(argc, argv) || conf_opt.NbofArgs() != 2) {
        std::cout << conf_opt.GetInstruction() << std::endl;
        return -1;
    }

    for(auto &opt : conf_opt.GetOptions())
    {
        switch(opt.mark)
        {
        case 'n':
            net_path = opt.var.String();
            break;
        case 's':
            save_path = opt.var.String();
            break;
        case 'f':
            learn_factor = opt.var.Double();
            break;
        case 't':
            learn_times = opt.var.Int();
            break;
        case 'c':
            cap = opt.var.Int();
            break;
        case 'l':
            layer_str = opt.var.String();
            break;
        default:
            std::cout << conf_opt.GetInstruction() << std::endl;
            return -1;
        }
    }

    string cosmic_file = conf_opt.GetArgument(0).String();
    string good_file = conf_opt.GetArgument(1).String();
    CNeuralNetwork my_net(learn_factor);

    vector<unsigned int> hidden;
    if(layer_str.empty())
    {
        // 4 hidden layers with 20, 10, 5, 3 neurons
        hidden = {20, 10, 5, 3};
    }
    else
    {
        auto vals = ConfigParser::split(layer_str, ",");
        while(vals.size())
        {
            string val = move(vals.front());
            vals.pop_front();
            hidden.push_back(stoi(val));
        }
    }

    if(net_path.empty())
    {
        // create net with dimensions
        // inputs and outputs are hard coded
        my_net.CreateNet(NN_INPUT, NN_OUTPUT, hidden);
        // initialize the weights with random values
        my_net.InitializeWeights();
    }
    else {
        // or create net from saved network data
        // exit if fail to create
        if(my_net.CreateNet(net_path.c_str()) == 0) {
            cout << "Failed to create the network." << endl;
            return -1;
        }
    }

    PRadHyCalSystem sys;
    sys.Configure("config/hycal.conf");

    NeuralTrain(my_net, sys, cosmic_file, good_file, learn_times*1000, cap*1000);

    my_net.SaveNet(save_path.c_str());
    return 0;
}

void NeuralTrain(CNeuralNetwork &net,
                 PRadHyCalSystem &sys,
                 const string &path,
                 const string &path2,
                 int number,
                 unsigned int cap)
{

    vector<vector<double>> cosmic_params;
    vector<double> cosmic_expect(NN_OUTPUT, 1.0);
    vector<vector<double>> good_params;
    vector<double> good_expect(NN_OUTPUT, 0.0);

    // read events and fill in the params array
    FillParams(sys, path, cosmic_params, cap);
    FillParams(sys, path2, good_params, cap);

    cout << "Start to train the network with files "
         << endl
         << "\"" << path << "\" (cosmic)."
         << endl
         << "\"" << path2 << "\"(good events)."
         << endl;

    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::uniform_int_distribution<unsigned int> cosmic_dist(0, cosmic_params.size()-1);
    std::uniform_int_distribution<unsigned int> good_dist(0, good_params.size()-1);

    int count = 0;
    PRadBenchMark timer;
    while(++count <= number)
    {
        if(count%PROGRESS_COUNT == 0) {
            cout << "----------training " << count
                 << "-------[ " << timer.GetElapsedTimeStr() << " ]------"
                 << "\r" << flush;
        }
        // randomly pick param sets
        auto &cosmic_input = cosmic_params.at(cosmic_dist(rng));
        net.BP_Train(cosmic_input, cosmic_expect);

        auto &good_input = good_params.at(good_dist(rng));
        net.BP_Train(good_input, good_expect);
    }

    cout << "----------training " << count - 1
         << "-------[ " << timer.GetElapsedTimeStr() << " ]------"
         << endl;
    cout << "Training finished, now test the result." << endl;

    // check the training result
    count = 0;
    double average1 = 0., average2 = 0.;
    while(++count <= 10000)
    {
        auto &cosmic_input = cosmic_params.at(cosmic_dist(rng));
        net.Update(cosmic_input);
        average1 += net.GetOutput().at(0);

        auto &good_input = good_params.at(good_dist(rng));
        net.Update(good_input);
        average2 += net.GetOutput().at(0);
    }

    cout << "The average output for cosmic is " << average1/(double)count
         << endl
         << "The average output for good event is " << average2/(double)count
         << endl;
}

void FillParams(PRadHyCalSystem &sys,
                const string &path,
                vector<vector<double>> &params,
                unsigned int cap)
{

    PRadDSTParser dst;
    dst.OpenInput(path);

    int count = 0;
    PRadBenchMark timer;

    cout << "Prepare inputs from file "
         << "\"" << path << "\" (Capacity: " << cap << ")"
         << endl;

    params.clear();
    params.reserve(cap);

    while(dst.Read())
    {
        if(dst.EventType() == PRadDSTParser::Type::event)
        {
            auto event = dst.GetEvent();

            if(!event.is_physics_event())
                continue;

            if((++count)%PROGRESS_COUNT == 0) {
                cout << "----------event " << count
                     << "-------[ " << timer.GetElapsedTimeStr() << " ]------"
                     << "\r" << flush;
            }

            auto param = AnalyzeEvent(&sys, event);
            if(param.group_size == 0)
                continue;

            params.emplace_back(move(param.GetParamList()));

            if(params.size() >= cap)
                break;
        }
    }

    cout << "----------event " << count
         << "-------[ " << timer.GetElapsedTimeStr() << " ]------"
         << endl;
}
