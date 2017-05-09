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
    double energy = 2142;
    TGraph *g1 = new TGraph();
    TGraph *g2 = new TGraph();
    TGraph *g3 = new TGraph();
    PRadMollerGen moller(1e-6, 10);
    for(double angle = 0.8; angle < 3.8; angle += 0.01)
    {
        double born, non_rad, rad;
        moller.GetXS(energy, angle, born, non_rad, rad);
        g1->SetPoint(g1->GetN(), angle, born);
        g2->SetPoint(g2->GetN(), angle, non_rad);
        g3->SetPoint(g3->GetN(), angle, (non_rad/born - 1.)*100.);
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
