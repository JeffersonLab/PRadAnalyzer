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

void moller_test()
{
    double energy = 1100;
    TGraph *g1 = new TGraph();
    TGraph *g2 = new TGraph();
    TGraph *g3 = new TGraph();
    PRadMollerGen moller(4.69, 100);
    for(double angle = 0.7; angle < 3.0; angle += 0.01)
    {
        double born, non_rad, rad;
        moller.GetXS(energy, angle, born, non_rad, rad);
        g1->SetPoint(g1->GetN(), angle, born);
        g2->SetPoint(g2->GetN(), angle, non_rad + rad);
        g3->SetPoint(g3->GetN(), angle, ((non_rad + rad)/born - 1.)*100.);
        cout << angle << ", " << born << ", " << non_rad + rad << ", " << non_rad << ", " << rad << endl;
    }

    TCanvas *c1 = new TCanvas("Moller XS", "Moller XS", 200, 10, 1200, 500);
    c1->Divide(2, 1);
    c1->SetGrid();

    c1->cd(1);
    g1->Draw("AC");
    g2->SetLineColor(2);
    g2->Draw("C");
    c1->cd(2);
    g3->Draw("AC");
}

void moller_gen_test(int Nevents, const char *path = "moller_test.dat")
{
    PRadMollerGen moller(5, 3000);
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
