#include "ConfigParser.h"
#include "CNeuralNetwork.h"
#include "canalib.h"

using namespace std;
const string prad_root = getenv("PRAD_PATH");

// profile energy and distance ranges and the step sizes
#define CLPROF_MIN_ENE 200      // max energy in the profile, MeV
#define CLPROF_MAX_ENE 2100     // min energy in the profile, MeV
#define CLPROF_MAX_DIST 5       // max distance in the profile, module size
#define CLPROF_STEP_ENE 100     // each step in energy
#define CLPROF_STEP_DIST 0.001  // each step in distance

// for training normailzation
#define ENERGY_NORM CLPROF_MAX_ENE + 500
#define DIST_NORM CLPROF_MAX_DIST + 1
#define LG_SIZE 38.15
#define PWO_SIZE 20.75

void train_profile(bool pwo = true, int Ntrain = 50000, string net = "")
{
    ConfigParser cp;
    struct DataPoint {double dist, ratio, error;};
    struct DataSet {double energy; vector<DataPoint> data;};
    vector<DataSet> data;

    string type = (pwo) ? "pwo" : "lg";
    double enorm = ENERGY_NORM;
    double dnorm = ((pwo) ? PWO_SIZE : LG_SIZE)*DIST_NORM;

    // read-in data
    for(int energy = CLPROF_MIN_ENE; energy < CLPROF_MAX_ENE; energy += CLPROF_STEP_ENE)
    {
        if(!cp.ReadFile("train_sets/" + type + "/prof_" + to_string(energy) + "MeV.txt")) {
            cout << "Cannot open file for energy = " << energy << endl;
            continue;
        }

        DataSet newset;
        newset.energy = energy/enorm;
        DataPoint newp;
        while(cp.ParseLine())
        {
            cp >> newp.dist >> newp.ratio >> newp.error;
            // quantized by step
            newp.dist /= dnorm;
            newset.data.push_back(newp);
        }

        data.emplace_back(newset);
    }

    CNeuralNetwork my_net(0.1);
    if(net.empty()) {
        vector<unsigned int> hidden = {10, 20};
        my_net.CreateNet(2, 2, hidden);
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
                my_net.BP_Train({dset.energy, p.dist}, {p.ratio, p.error});
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
    if(!cp.ReadFile("train_sets/" + type + "/prof_" + to_string(energy) + "MeV.txt")) {
        cout << "Cannot open file for energy = " << energy << endl;
        return;
    }

    TGraph *grd = new TGraph;
    TGraph *gre = new TGraph;
    TGraph *gr1 = new TGraph();
    TGraph *gr2 = new TGraph();
    double dist, ratio, error;
    double e = double(energy)/double(ENERGY_NORM);
    double dnorm = ((pwo) ? PWO_SIZE : LG_SIZE)*DIST_NORM;
    while(cp.ParseLine())
    {
        cp >> dist >> ratio >> error;
        grd->SetPoint(grd->GetN(), dist, ratio);
        gre->SetPoint(gre->GetN(), dist, error);
        my_net.Update({e, dist/dnorm});
        gr1->SetPoint(gr1->GetN(), dist, my_net.GetOutput()[0]);
        gr2->SetPoint(gr2->GetN(), dist, my_net.GetOutput()[1]);
    }

    grd->SetMarkerStyle(20);
    grd->Draw("AP");
    gre->SetMarkerStyle(20);
    gre->Draw("P");
    gr1->SetLineColor(2);
    gr1->Draw("C");
    gr2->SetLineColor(4);
    gr2->Draw("C");

}

void gen_profile(string net, string output)
{
    ofstream outf(output);
    CNeuralNetwork my_net;
    my_net.CreateNet(net.c_str());

    int Ne = (CLPROF_MAX_ENE - CLPROF_MIN_ENE)/CLPROF_STEP_ENE + 1;
    int Nd = CLPROF_MAX_DIST/CLPROF_STEP_DIST + 1;
    outf << "# min energy, max energy, energy step, max distance, distance step\n"
         << CLPROF_MIN_ENE << ", " << CLPROF_MAX_ENE << ", " << CLPROF_STEP_ENE
         << ", " << CLPROF_MAX_DIST << ", " << CLPROF_STEP_DIST
         << endl;
    for(int ie = 0; ie < Ne; ++ie)
    {
        double e = double(CLPROF_MIN_ENE + CLPROF_STEP_ENE*ie)/double(ENERGY_NORM);

        for(int id = 0; id < Nd; ++id)
        {
            double d = id*CLPROF_STEP_DIST/double(DIST_NORM);
            my_net.Update({e, d});
            outf << setw(8) << ie
                 << setw(8) << id
                 << setw(15) << my_net.GetOutput().at(0)
                 << setw(15) << my_net.GetOutput().at(1)
                 << endl;
        }
    }

    outf.close();
}
