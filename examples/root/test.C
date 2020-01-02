#include "ConfigParser.h"
#include "ConfigObject.h"
#include "PRadMollerGen.h"
#include "PRadEpElasGen.h"
#include "canalib.h"

using namespace std;

// helper function
inline bool compare_str(const char *buf1, const char *buf2, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        if (buf1[i] != buf2[i]) {
            return false;
        }
    }
    return true;
}

// helper function
inline std::string break_line(const std::string &buf, const std::string &delim, size_t &i, const std::string &glue)
{
    size_t beg = i;
    if (delim.empty()) {
        i = buf.size();
        return buf.substr(beg);
    }

    for (; i <= buf.size() - delim.size(); ++i) {
        if (compare_str(&buf[i], delim.c_str(), delim.size())) {
            std::string res = buf.substr(beg, i - beg);
            i += delim.size();
            if (glue.size() && (res.size() > glue.size())) {
                size_t pos = res.size() - glue.size();
                if (compare_str(&res[pos], glue.c_str(), glue.size())) {
                    return res.substr(0, pos) + break_line(buf, delim, i, glue);
                }
            }
            return res;
        }
    }
    i = buf.size();
    return buf.substr(beg);
}

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

    string var1 = conf_obj.Value<string>("Var1");
    string var2 = conf_obj.Value<string>("Var2");
    string var3 = conf_obj.Value<string>("Var3");
    string var4 = conf_obj.Value<string>("Var4");
    string var5 = conf_obj.Value<string>("Var5");

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
    auto text = ConfigParser::break_into_blocks(buffer, "{", "}");

    cout << "residual: " << text.residual << endl;

    for(auto &b : text.blocks)
    {
        cout << "label: " << b.label << endl;
        cout << "content: " << b.content << endl;
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

void moller_test(double v_max = 1000)
{
    TGraph *g1a = new TGraph();
    TGraph *g2a = new TGraph();
    TGraph *g3a = new TGraph();
    TGraph *g1b = new TGraph();
    TGraph *g2b = new TGraph();
    TGraph *g3b = new TGraph();
    PRadMollerGen moller1(1, v_max);
    PRadMollerGen moller2(10, v_max);
    PRadMollerGen moller3(100, v_max);
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
    PRadMollerGen moller(1, 3000, 100, 1e-4, 1e-4);
    moller.Generate(2142, 0.3, 3.0, Nevents, path);
}

void show_moller_gen(const char *path)
{
    ConfigParser c_parser;
    c_parser.ReadFile(path);

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

void show_ep_gen(const char *path)
{
    ConfigParser c_parser;
    c_parser.ReadFile(path);

    double p1, p2, p, th1, th2, th, ph1, ph2, ph;
    TH1F *hist_th = new TH1F("theta dist", "theta dist", 200, 5, 13);
    TH1F *hist_E = new TH1F("energy dist", "energy dist", 200, 0, 2200);
    double m = cana::ele_mass;
    while(c_parser.ParseLine())
    {
        c_parser >> p1 >> th1 >> ph1 >> p2 >> th2 >> ph2 >> p >> th >> ph;
        hist_E->Fill(sqrt(p1*p1 + m*m));
        hist_th->Fill(th1*cana::rad2deg);
    }

    TCanvas *c1 = new TCanvas("Elastic XS", "Elastic XS", 200, 10, 1200, 500);
    c1->Divide(2, 1);
    c1->cd(1);
    hist_E->Draw("colz");

    c1->cd(2);
    hist_th->Draw("colz");
}

void interpolate_F1F2(const char *path, const char *outf, double step = 1.0)
{
    ConfigParser c_parser;
    c_parser.ReadFile(path);

    double nu, F1, F2;
    std::vector<double> vnu, vF1, vF2;
    TGraph *f1a, *f2a, *f1b, *f2b;
    f1a = new TGraph();
    f1b = new TGraph();
    f2a = new TGraph();
    f2b = new TGraph();

    // read table
    while(c_parser.ParseLine())
    {
        c_parser >> nu >> F1 >> F2;
        vnu.push_back(nu);
        vF1.push_back(F1);
        vF2.push_back(F2);
        f1a->SetPoint(f1a->GetN(), nu, F1);
        f2a->SetPoint(f2a->GetN(), nu, F2);
    }

    // spline fit
    TSpline3 F1_spl("F1 Spline", &vnu[0], &vF1[0], vnu.size(), "b3e3");
    TSpline3 F2_spl("F2 Spline", &vnu[0], &vF2[0], vnu.size(), "b3e3");

    // output
    ofstream output(outf);
    for(double nu_val = vnu.front(); nu_val <= vnu.back(); nu_val += 1.0)
    {
        double F1_interp = F1_spl.Eval(nu_val);
        double F2_interp = F2_spl.Eval(nu_val);
        output << setw(8) << nu_val
               << setw(15) << F1_interp
               << setw(15) << F2_interp
               << endl;
        f1b->SetPoint(f1b->GetN(), nu_val, F1_interp);
        f2b->SetPoint(f2b->GetN(), nu_val, F2_interp);
    }

    f1a->SetMarkerStyle(20);
    f1a->SetMarkerColor(2);
    f1b->SetMarkerStyle(22);
    f2a->SetMarkerStyle(20);
    f2a->SetMarkerColor(2);
    f2b->SetMarkerStyle(22);

    // show the interpolate result
    TCanvas *c1 = new TCanvas("Interp_F1F2", "F1 F2 Interpolation", 200, 10, 1200, 500);
    c1->Divide(2, 1);
    c1->cd(1);
    f1a->Draw("AP");
    f1b->Draw("P");
    c1->cd(2);
    f2a->Draw("AP");
    f2b->Draw("P");
}

void ep_test(double energy)
{
    double S = 2.*energy*cana::proton_mass;
    TCanvas *c1 = new TCanvas("Elastic Ep XS", "EpXS", 200, 10, 700, 500);
    TGraph *g1 = new TGraph();

    PRadEpElasGen ep_gen;

    for(double logq2 = -6; logq2 < -2; logq2 += 0.01)
    {
        double q2 = std::pow(10., logq2)*1e6;
        double born, non_rad, rad;
        ep_gen.GetXSdQsq(S, q2, born, non_rad, rad);
        std::cout << q2 << ", " << born << ", " << non_rad << std::endl;
        g1->SetPoint(g1->GetN(), q2, (non_rad/born - 1.)*100.);
    }

    c1->SetLogx();
    g1->Draw("AP");
}

void ep_vmin_test(double energy = 2142, double v_max = 1000)
{
    TGraph *g1 = new TGraph();
    TGraph *g2 = new TGraph();
    TGraph *g3 = new TGraph();
    TGraph *g4 = new TGraph();
    TGraph *g5 = new TGraph();
    TGraph *g6 = new TGraph();

    PRadEpElasGen ep1(1., v_max);
    PRadEpElasGen ep2(50., v_max);
    PRadEpElasGen ep3(200., v_max);

    double S = 2.*energy*cana::proton_mass;
    for(double logq2 = -4; logq2 < -2; logq2 += 0.03)
    {
        double Q2 = std::pow(10., logq2)*1e6;
        double born, non_rad, rad, xs1, xs2, xs3;

        ep1.GetXSdQsq(S, Q2, born, non_rad, rad);
        std::cout << non_rad << ", " << rad << std::endl;
        xs1 = non_rad + rad;
        ep2.GetXSdQsq(S, Q2, born, non_rad, rad);
        xs2 = non_rad + rad;
        ep3.GetXSdQsq(S, Q2, born, non_rad, rad);
        xs3 = non_rad + rad;

        g1->SetPoint(g1->GetN(), Q2, xs1);
        g2->SetPoint(g2->GetN(), Q2, xs2);
        g3->SetPoint(g3->GetN(), Q2, xs3);

        g4->SetPoint(g4->GetN(), Q2, 0.);
        g5->SetPoint(g5->GetN(), Q2, (xs2 - xs1)/xs1*100.);
        g6->SetPoint(g6->GetN(), Q2, (xs3 - xs1)/xs1*100.);

        std::cout << Q2 << ", " << (xs2 - xs1)/xs1*100. << ", " << (xs3 - xs1)/xs1*100. << std::endl;
    }

    TCanvas *c1 = new TCanvas("v_min test", "v_min test", 200, 10, 1200, 500);
    c1->Divide(2, 1);
    c1->SetGrid();
    c1->SetLogx();

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


