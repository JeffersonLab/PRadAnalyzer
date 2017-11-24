#include "ConfigParser.h"
#include "cneural/CNeuralNetwork.h"
#include "canalib.h"

using namespace std;
const string prad_root = getenv("PRAD_PATH");

#define MAX_PROFILE_ENERGY 3000
#define MAX_PROFILE_DIST 5
#define LG_SIZE 38.15
#define PWO_SIZE 20.75

void train_profile(string net = "", bool pwo = true, int Ntrain = 50000)
{
    ConfigParser cp;
    struct DataPoint {double dist, ratio, error;};
    struct DataSet {double energy; vector<DataPoint> data;};
    vector<DataSet> data;

    string type = (pwo) ? "pwo" : "lg";
    double size = (pwo) ? PWO_SIZE : LG_SIZE;

    // read-in data
    for(int energy = 200; energy < MAX_PROFILE_ENERGY; energy += 100)
    {
        if(!cp.ReadFile(prad_root + "/database/cluster_profiles/" + type + "/prof_" + to_string(energy) + "MeV.txt"))
            continue;

        DataSet newset;
        newset.energy = energy/double(MAX_PROFILE_ENERGY);
        DataPoint newp;
        while(cp.ParseLine())
        {
            cp >> newp.dist >> newp.ratio >> newp.error;
            // quantized by size
            newp.dist /= size;
            // normailze
            newp.dist /= double(MAX_PROFILE_DIST);
            newset.data.push_back(newp);
        }

        data.emplace_back(newset);
    }

    CNeuralNetwork my_net(0.1);
    if(net.empty()) {
        vector<unsigned int> hidden = {10, 20, 5};
        my_net.CreateNet(2, 1, hidden);
        my_net.InitializeWeights();
        net = "prof_" + type + ".net";
    } else {
        my_net.CreateNet(net.c_str());
    }

    int iter = 1;
    while(iter++ <= Ntrain)
    {
        for(auto &dset : data)
        {
            for(auto &p : dset.data)
            {
                my_net.BP_Train({dset.energy, p.dist}, {p.ratio});
            }
        }
        cout << "Training : " << iter << "/" << Ntrain << "\r" << flush;
    }
    cout << "Training : " << Ntrain << "/" << Ntrain << "\r" << endl;

    my_net.SaveNet(net.c_str());
}

void test_profile(string net, bool pwo = true, int energy = 1000)
{
    CNeuralNetwork my_net;
    my_net.CreateNet(net.c_str());
    ConfigParser cp;
    string type = (pwo) ? "pwo" : "lg";
    if(!cp.ReadFile(prad_root + "/database/cluster_profiles/" + type + "/prof_" + to_string(energy) + "MeV.txt"))
        return;

    TGraphErrors *gre = new TGraphErrors;
    TGraph *gr = new TGraph();
    double dist, ratio, error;
    double e = double(energy)/double(MAX_PROFILE_ENERGY);
    double size = (pwo) ? PWO_SIZE : LG_SIZE;
    while(cp.ParseLine())
    {
        cp >> dist >> ratio >> error;
        auto N = gre->GetN();
        gre->SetPoint(N, dist, ratio);
        gre->SetPointError(N, 0., error);
        my_net.Update({e, dist/(size*MAX_PROFILE_DIST)});
        gr->SetPoint(gr->GetN(), dist, my_net.GetOutput()[0]);
    }

    gre->SetMarkerStyle(20);
    gre->Draw("AP");
    gr->SetLineColor(2);
    gr->Draw("C");

}

