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
    for(double angle = 0.05; angle < 3.8; angle += 0.01)
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

void moller_gen_test()
{
    PRadMollerGen moller;
    moller.Generate(2142, 0.3, 3.0, 10000000);
}

void four_momentum_boost_z(double *p, double *pp, double beta)
{
    double gamma = sqrt(1./(1. - beta*beta));
    p[0] = pp[0];
    p[1] = pp[1];
    p[2] = gamma*(pp[2] - beta*pp[3]);
    p[3] = gamma*(pp[3] - beta*pp[2]);
}

void Lorentz_test(double energy = 1100., double angle_CM = 90)
{
    double m = cana::ele_mass;
    double theta_CM = angle_CM*cana::deg2rad;
    double beta_CM = sqrt((energy - m)/(energy + m));

    double p1[4] = {0., 0., sqrt(energy*energy - m*m), energy};
    double k1[4] = {0., 0., 0., m};
    double p1_CM[4], k1_CM[4];
    four_momentum_boost_z(p1_CM, p1, beta_CM);
    four_momentum_boost_z(k1_CM, k1, beta_CM);

    double p_CM = sqrt(p1_CM[0]*p1_CM[0] + p1_CM[1]*p1_CM[1] + p1_CM[2]*p1_CM[2]);
    double E_CM = p1_CM[3];

    cout << p_CM << ", " << sqrt(m*(energy - m)/2.) << ", "
         << E_CM << ", " << sqrt(m*(energy + m)/2.)
         << endl;

    double p2_CM[4] = {-sin(theta_CM)*p_CM, 0., -cos(theta_CM)*p_CM, E_CM};
    double k2_CM[4] = {sin(theta_CM)*p_CM, 0., cos(theta_CM)*p_CM, E_CM};
    double p2[4], k2[4];

    four_momentum_boost_z(k2, k2_CM, -beta_CM);
    four_momentum_boost_z(p2, p2_CM, -beta_CM);

    double tangent1 = sqrt((k2[0]*k2[0] + k2[1]*k2[1])/k2[2]/k2[2]);
    double tangent2 = sqrt((p2[0]*p2[0] + p2[1]*p2[1])/p2[2]/p2[2]);
    double angle1 = atan(tangent1)*cana::rad2deg;
    double angle2 = atan(tangent2)*cana::rad2deg;
    cout << angle1 << ", " << angle2 << ", "
         << acos(sqrt((energy + m)/(energy + 3.*m)))*cana::rad2deg
         << endl;
}

    // equation (A.14) - (A.15)
double S_phi(const double &s_1, const double &s_2, const double &s_3)
{
    double m = cana::ele_mass, m2 = m*m;
    double lamda_1 = s_1*s_1 - 16.*m2*m2, slamda_1 = sqrt(lamda_1);
    double lamda_2 = s_2*s_2 - 16.*m2*m2, slamda_2 = sqrt(lamda_2);
    double lamda_3 = s_3*s_3 - 16.*m2*m2, slamda_3 = sqrt(lamda_3);
    // z_u and z_d
    double z_ud[2] = {slamda_1/slamda_2 - 1., (s_1*s_2 - 4.*m2*s_3)/lamda_2 - 1.};
    // z_1, z_2, z_3, z_4
    double z[4] = {1./slamda_2*(4.*m2*(s_3 - slamda_3)/(s_2 - slamda_2) - s_1 - slamda_2),
                   1./slamda_2*(4.*m2*(s_3 + slamda_3)/(s_2 - slamda_2) - s_1 - slamda_2),
                   1./slamda_2*(s_1 - slamda_2 - 4.*m2*(s_3 + slamda_3)/(s_2 + slamda_2)),
                   1./slamda_2*(s_1 - slamda_2 - 4.*m2*(s_3 - slamda_3)/(s_2 + slamda_2))};

    // Sj
    double Sj[4] = {1, 1, -1, -1};
    // (-1)^(i + 1), i from 1 to 4 but index is from 0 to 3
    double Si[4] = {1, -1, 1, -1};
    // z_u term - z_d term
    double Sk[2] = {1, -1};

    double result = 0.;
    for(int k = 0; k < 2; ++k)
    {
        // first term
        double term = log((s_2 - slamda_2)/(s_2 + slamda_2))
                      * log((z_ud[k] - z[0])*(z_ud[k] - z[2])/(z_ud[k] - z[1])/(z_ud[k] - z[3]));

        // second term
        double sum_term = 0.;
        for(int i = 0; i < 4; ++i)
        {
            for(int j = 0; j < 4; ++j)
            {
                double inner_term;
                double logzi = log(fabs(z_ud[k] - z[i]));
                if(i == j) {
                    inner_term = logzi;
                } else {
                    // the input to log may be negative and thus
                    // the result may be a complex number
                    // TODO check with authors for this part
                    double log_term = logzi*log(fabs(z[j] - z[i]));
                    double spence_term = cana::spence((z_ud[k] - z[i])/(z[j] - z[i]));
                    inner_term = log_term - spence_term;
                }
                sum_term += Sj[j]*Si[i]*inner_term;
                cout << i+1 << ", " << j+1 << ", " << Sj[j]*Si[i]*inner_term << ", " << sum_term << endl;
            }
        }
        cout << s_3/2./slamda_3*sum_term << endl;

        result += s_3/2./slamda_3*(term + sum_term)*Sk[k];
    }
    return result;
}

void S_phi_test()
{
    double s1 = 961.632;
    double s2 = 2248.4;
    double s3 = 1287.81;

    cout << S_phi(s1, s2, s3) << ", " << S_phi(s2, s1, s3) << endl;
}

