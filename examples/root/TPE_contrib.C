#include "ConfigParser.h"
#include "canalib.h"

using namespace std;

struct TPE_params
{
    double D, m, mu;
    double a[4][6];
};

TPE_params read_TPE_params(std::string path)
{
    ConfigParser parser;
    parser.ReadFile("TPE_params.dat");

    TPE_params param;

    if (parser.ParseLine()) {
        parser >> param.D >> param.m >> param.mu;
    }

    int l = 0;
    while (parser.ParseLine())
    {
        for (int i = 0; i < 6; ++i)
        {
            parser >> param.a[l][i];
        }
        l++;
    }
    return param;
}

// GeV
double calc_TPE_contrib(double Q2, double eps, const TPE_params &param)
{
    double eps_x[4] = {std::sqrt(1. - eps),
                       std::sqrt(eps*(1. - eps)),
                       1. - eps,
                       (1. - eps)*(1. - eps)};

    double Q2_x[6] = {1., Q2, std::log(Q2 + param.D), Q2*std::log(Q2 + param.D),
                      1./(1. + Q2/(param.m*param.m)), 1./(1. + Q2/(param.mu*param.mu))};

    double result = 0.;
    for (int i = 0; i < 4; ++i)
    {
        double a = 0.;
        for (int j = 0; j < 6; ++j)
        {
            a += param.a[i][j]*Q2_x[j];
        }
        result += a*eps_x[i];
    }

    return result;
}

void show_TPE_contrib()
{
    std::vector<double> energies = {1.097, 2.142};
    std::vector<TGraph*> plots;
    std::vector<int> colors = {2, 4};
    double angle_min = 0.3, angle_max = 8.0;
    int angle_bins = 1000;
    double angle_step = (angle_max - angle_min)/(double)angle_bins;

    auto param = read_TPE_params("TPE_params.dat");

    TCanvas *c1 = new TCanvas("unpol tail","unpol tail", 10, 10, 700, 500);
    c1->SetGrid();
    auto fr = c1->DrawFrame(1e-4, -0.1, 1.0, 0.1);
    fr->GetYaxis()->SetTitle("#deltaG_{E}/G_{E} (%)");
    fr->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
    c1->SetLogx();
    auto legend = new TLegend(0.7,0.7,0.9,0.9);
    for (size_t k = 0; k < energies.size(); ++k)
    {
        double E = energies[k];
        TGraph *gr = new TGraph;
        gr->SetLineColor(colors[k]);
        gr->SetMarkerColor(colors[k]);
        for (int i = 0; i < angle_bins; ++i)
        {
            double angle = angle_min + angle_step*i;
            double sin2 = cana::pow2(std::sin(angle*cana::deg2rad/2.));
            double nu = E - E/(1. + 2.*E/(cana::proton_mass/1e3)*sin2);
            double Q2 = 4.*E*(E - nu)*sin2;
            double eps = 1./(1. + 2.*(nu*nu/Q2 + 1.)*cana::pow2(tan(angle/2.*cana::deg2rad)));

            double delta_GE = calc_TPE_contrib(Q2, eps, param)*cana::mup;
            gr->SetPoint(gr->GetN(), Q2, delta_GE*100.);
        }
        legend->AddEntry(gr, (std::to_string(int(E*1000. + 0.5)) + " MeV").c_str(), "lp");
        gr->Draw("l");
    }
    legend->Draw();
}
