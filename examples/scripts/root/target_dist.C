#include "ConfigParser.h"
#include "canalib.h"
#include "time.h"

using namespace std;

struct TarDistPoint
{
    double z, density, cdf;

    // define operators for binary search of cdf
    bool operator <(const double &rhs) const {return cdf < rhs;}
    bool operator >(const double &rhs) const {return cdf > rhs;}
    bool operator ==(const double &rhs) const {return cdf == rhs;}
    bool operator !=(const double &rhs) const {return cdf != rhs;}
};

// read-in points from the file
// hear points must be sorted in z transcendent order, other wise the probability
// calculation will fail
vector<TarDistPoint> fill_points(const string &path)
{
    vector<TarDistPoint> res;

    ConfigParser c_parser;
    c_parser.ReadFile(path);

    double integral = 0.;
    while(c_parser.ParseLine())
    {
        TarDistPoint p;
        c_parser >> p.z >> p.density;

        // convert z to mm, and convert density to H_atom/cm^3
        p.z *= 10.;
        p.density *= 2.;

        // calculate the cumulative probability function
        // first point
        if(res.empty()) {
            p.cdf = 0;
        } else {
            // linear approximation between two points
            double mean_density = (res.back().density + p.density)/2.;
            double bin_size = res.back().z - p.z;
            p.cdf = res.back().cdf + mean_density*bin_size;
        }

        res.push_back(p);
    }

    // normalize the probability
    for(auto &p : res)
    {
        p.cdf /= res.back().cdf;
    }

    return res;
}

// plot target density distribution
void show_target_dist()
{
    string prad_path = getenv("PRAD_PATH");
    string empty_file = prad_path + "/database/target_comsol/empty_target.dat";
    string full_file = prad_path + "/database/target_comsol/full_target.dat";

    vector<TarDistPoint> full_points = fill_points(full_file);
    vector<TarDistPoint> empty_points = fill_points(empty_file);

    TGraph *g1 = new TGraph();
    TGraph *g2 = new TGraph();
    for(auto &p : full_points)
    {
        g1->SetPoint(g1->GetN(), p.z, p.density);
    }

    for(auto &p: empty_points)
    {
        g2->SetPoint(g2->GetN(), p.z, p.density);
    }

    TCanvas *c1 = new TCanvas("tar_dist", "target distribution", 200, 10, 700, 500);
    c1->SetGrid();
    c1->SetLogy();

    g1->Draw("LAP");
    g2->SetLineColor(2);
    g2->SetMarkerColor(2);
    g2->Draw("LP");
}

void monte_carlo_generator(size_t event_number)
{
    string prad_path = getenv("PRAD_PATH");
    string target_file = prad_path + "/database/target_comsol/empty_target.dat";
    vector<TarDistPoint> data = fill_points(target_file);

    TFile f("z_target.root", "RECREATE");
    TH1D z_hist("Z_TAR", "Target Z", 1600, -800, 800);
    TRandom3 rnd_gen(time(NULL));

    while(event_number--)
    {
        double z_target;
        double rnd_val = rnd_gen.Rndm();
        auto res = cana::binary_search_interval(data.begin(), data.end(), rnd_val);

        // should not happen since it is in [0, 1]
        if(res.first == data.end() && res.second == data.end()) {
            cout << "Cannot find interval for random number " << rnd_val << endl;
            continue;
        }

        // it happenly sit exactly on one point
        if(res.first == res.second) {
            z_target = res.first->z;
        } else {
            // linear interpolation
            z_target = (res.first->z * (rnd_val - res.first->cdf)
                        + res.second->z * (res.second->cdf - rnd_val));
            z_target /= res.second->cdf - res.first->cdf;
        }

        z_hist.Fill(z_target);
    }

    z_hist.Write();
    f.Close();
}
