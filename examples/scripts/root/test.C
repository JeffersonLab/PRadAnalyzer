#include "ConfigParser.h"
#include "ConfigObject.h"
#include "PRadMollerGen.h"
#include "canalib.h"

using namespace std;

void test_find_pair()
{
    string tstr = "{|}asvs{asdasdss}asdasd}";
    string op = "{";
    string cl = "}";
    auto p = ConfigParser::find_pair(tstr, op, cl);

    if(p.first == string::npos || p.second == string::npos) {
        cout << "not found" << endl;
    } else {
        cout << tstr.substr(p.first + op.size(), p.second - op.size() - p.first) << endl;
    }
}

void test_reform()
{
    ConfigObject conf_obj;

    conf_obj.ReadConfigString("Var1 = 2 \n"
             "Var2 = 1{Var1}3456 \n"
             "Var3 = {Var{Var1}}789 \n"
             "Var4 = {Var{Var1}}{Var1}{Var1} \n"
             "Path = variable_test.conf \n"
             "INCLUDE({Path})");

    string var1 = conf_obj.GetConfig<string>("Var1");
    string var2 = conf_obj.GetConfig<string>("Var2");
    string var3 = conf_obj.GetConfig<string>("Var3");
    string var4 = conf_obj.GetConfig<string>("Var4");
    string var5 = conf_obj.GetConfig<string>("Var5");

    cout << var1 << endl
         << var2 << endl
         << var3 << endl
         << var4 << endl
         << var5 << endl;
}

void test_block_read(const string &path = "block_test.conf")
{
    // read in file
    string buffer = ConfigParser::file_to_string(path);

    // remove comments
    ConfigParser::comment_between(buffer, "/*", "*/");
    ConfigParser::comment_line(buffer, "//", "\n");
    ConfigParser::comment_line(buffer, "#", "\n");

    // get content blocks
    auto blocks = ConfigParser::break_into_blocks(buffer, "{", "}");

    for(auto &b : blocks)
    {
        cout << b.label << endl;
        cout << b.content << endl;
    }
}

// Mandelstam variables for Moller process
inline void get_moller_stu(double Es, double angle, double &s, double &t, double &u)
{

    double m = cana::ele_mass, m2 = m*m;
    double theta = angle*cana::deg2rad;
    double k1p = sqrt(Es*Es - m2);
    double cosE = (Es - m)*std::pow(cos(theta), 2);
    double Ep = m*(Es + m + cosE)/(Es + m - cosE);
    double k2p = sqrt(Ep*Ep - m*m);

    s = 2.*m*(Es + m);
    t = std::pow(Es - Ep, 2) - std::pow(k2p*sin(theta), 2) - std::pow(k2p*cos(theta) - k1p, 2);
    u = 4.*m2 - s - t;
}

void moller_test(double v_min = 0.5)
{
    TGraph *g1a = new TGraph();
    TGraph *g2a = new TGraph();
    TGraph *g3a = new TGraph();
    TGraph *g1b = new TGraph();
    TGraph *g2b = new TGraph();
    TGraph *g3b = new TGraph();
    PRadMollerGen moller1(1, 1000);
    PRadMollerGen moller2(10, 1000);
    PRadMollerGen moller3(100, 1000);
    double s1 = 2.*cana::ele_mass*(1100. + cana::ele_mass);
    double s2 = 2.*cana::ele_mass*(2200. + cana::ele_mass);
    for(double logq2 = -6; logq2 < -2; logq2 += 0.01)
    {
        double q2 = std::pow(10., logq2);
        double t = -q2*1e6;
        double born = 0., non_rad = 0., rad = 0.;
        if(4.*cana::ele_mass*cana::ele_mass - s1 - t < -120.) {
            moller1.GetXSdQsq(s1, t, born, non_rad, rad);
            g1a->SetPoint(g1a->GetN(), q2, ((non_rad)/born - 1.)*100.);
            moller2.GetXSdQsq(s1, t, born, non_rad, rad);
            g2a->SetPoint(g2a->GetN(), q2, ((non_rad)/born - 1.)*100.);
            moller3.GetXSdQsq(s1, t, born, non_rad, rad);
            g3a->SetPoint(g3a->GetN(), q2, ((non_rad)/born - 1.)*100.);
        }
        if(4.*cana::ele_mass*cana::ele_mass - s2 - t < -120.) {
            moller1.GetXSdQsq(s2, t, born, non_rad, rad);
            g1b->SetPoint(g1b->GetN(), q2, ((non_rad)/born - 1.)*100.);
            moller2.GetXSdQsq(s2, t, born, non_rad, rad);
            g2b->SetPoint(g2b->GetN(), q2, ((non_rad)/born - 1.)*100.);
            moller3.GetXSdQsq(s2, t, born, non_rad, rad);
            g3b->SetPoint(g3b->GetN(), q2, ((non_rad)/born - 1.)*100.);
        }
        cout << logq2 << endl;
    }

    g1b->SetLineStyle(9);
    g2b->SetLineStyle(9);
    g3b->SetLineStyle(9);
    g1a->SetLineColor(2);
    g2a->SetLineColor(1);
    g3a->SetLineColor(8);
    g1b->SetLineColor(2);
    g2b->SetLineColor(1);
    g3b->SetLineColor(8);
    TCanvas *c1 = new TCanvas("Moller XS", "MollerXS", 200, 10, 700, 500);
    c1->SetGrid();
    c1->DrawFrame(1e-6, -35, 1e-2, 5);
    c1->SetLogx();
    g1a->Draw("C");
    g2a->Draw("C");
    g3a->Draw("C");
    g1b->Draw("C");
    g2b->Draw("C");
    g3b->Draw("C");
}

void moller_gen_test(int Nevents, const char *path = "moller_test.dat")
{
    PRadMollerGen moller(1, 3000);
    moller.Generate(2142, 0.3, 3.0, Nevents, path);
}

void show_moller_gen(const char *path)
{
    ConfigParser c_parser;
    c_parser.OpenFile(path);

    double p1, p2, p, th1, th2, th, ph1, ph2, ph;
    TH1F *hist_th = new TH1F("theta dist", "theta dist", 200, 0, 3);
    TH2F *hist_2d = new TH2F("2d dist", "theta vs. energy", 200, 0, 3, 200, 0, 2200);
    double m = cana::ele_mass;
    while(c_parser.ParseLine())
    {
        c_parser >> p1 >> th1 >> ph1 >> p2 >> th2 >> ph2 >> p >> th >> ph;
        hist_2d->Fill(th1*cana::rad2deg, sqrt(p1*p1 + m*m));
        hist_2d->Fill(th2*cana::rad2deg, sqrt(p2*p2 + m*m));
        hist_th->Fill(th1*cana::rad2deg);
        hist_th->Fill(th2*cana::rad2deg);
        if(p > 0.) {
            hist_2d->Fill(th*cana::rad2deg, sqrt(p*p + m*m));
        }
    }

    TCanvas *c1 = new TCanvas("Moller XS", "Moller XS", 200, 10, 1200, 500);
    c1->Divide(2, 1);
    c1->cd(1);
    hist_2d->Draw("colz");

    c1->cd(2);
    hist_th->Draw("colz");
}

void moller_vmin_test(double energy = 2142, double v_max = 1000)
{
    TGraph *g1 = new TGraph();
    TGraph *g2 = new TGraph();
    TGraph *g3 = new TGraph();
    TGraph *g4 = new TGraph();
    TGraph *g5 = new TGraph();
    TGraph *g6 = new TGraph();

    PRadMollerGen moller1(1., v_max);
    PRadMollerGen moller2(5., v_max);
    PRadMollerGen moller3(20., v_max);

    for(double angle = 0.3; angle < 3.0; angle += 0.01)
    {
        double born, non_rad, rad, xs1, xs2, xs3;
        moller1.GetXS(energy, angle, born, non_rad, rad);
        xs1 = non_rad + rad;
        moller2.GetXS(energy, angle, born, non_rad, rad);
        xs2 = non_rad + rad;
        moller3.GetXS(energy, angle, born, non_rad, rad);
        xs3 = non_rad + rad;

        g1->SetPoint(g1->GetN(), angle, xs1);
        g2->SetPoint(g2->GetN(), angle, xs2);
        g3->SetPoint(g3->GetN(), angle, xs3);

        g4->SetPoint(g4->GetN(), angle, 0.);
        g5->SetPoint(g5->GetN(), angle, (xs2 - xs1)/xs1*100.);
        g6->SetPoint(g6->GetN(), angle, (xs3 - xs1)/xs1*100.);
    }

    TCanvas *c1 = new TCanvas("v_min test", "v_min test", 200, 10, 1200, 500);
    c1->Divide(2, 1);
    c1->SetGrid();

    g2->SetLineColor(2);
    g5->SetLineColor(2);
    g3->SetLineColor(4);
    g6->SetLineColor(4);

    c1->cd(1);
    g1->Draw("AC");
    g2->Draw("C");
    g3->Draw("C");

    c1->cd(2);
    g6->Draw("AC");
    g5->Draw("C");
    g4->Draw("C");
}

