#include "cana_interp.h"
#include <iostream>



bool cana::interp_1d(const std::vector<PointErr> &data, PointErr &p, double res)
{
    // safety checks
    if(data.size() < 2) return false;
    if(p.x < data.front().x) {
        if(std::abs(p.x - data.front().x) < std::abs(res*data.front().x)) {
            p = data.front();
            return true;
        }
        return false;
    }
    if(p.x > data.back().x) {
        if(std::abs(p.x - data.back().x) < std::abs(res*data.back().x)) {
            p = data.back();
            return true;
        }
        return false;
    }
    auto itp = binary_search_interval(data.begin(), data.end(), p.x);

    // exactly matched
    if(itp.first == itp.second) {
        p = *itp.first;
        return true;
    }

    // count how many points can be used
    // we will use mostly 4 points, 2 before and 2 after
    if(itp.first - data.begin() > 0) itp.first--;
    if(data.end() - itp.second > 1) itp.second++;

    int Np = itp.second - itp.first + 1;

    // linear interpolation
    if(Np == 2) {
        double c[2];
        auto &p1 = *itp.first;
        auto &p2 = *(itp.first + 1);
        linear_coeff(p1.x, p2.x, p.x, &c[0]);
        p.y = c[0]*p1.y + c[1]*p2.y;
        p.syst = c[0]*p1.syst + c[1]*p2.syst;
        p.stat = std::sqrt(pow2(c[0]*p1.stat) + pow2(c[1]*p2.stat));
    // quadratic interpolation
    } else if(Np == 3) {
        double c[3];
        auto &p1 = *itp.first;
        auto &p2 = *(itp.first + 1);
        auto &p3 = *(itp.first + 2);
        quad_coeff(p1.x, p2.x, p3.x, p.x, &c[0]);
        p.y = c[0]*p1.y + c[1]*p2.y + c[2]*p3.y;
        p.syst = c[0]*p1.syst + c[1]*p2.syst + c[2]*p3.syst;
        p.stat = std::sqrt(pow2(c[0]*p1.stat) + pow2(c[1]*p2.stat) + pow2(c[2]*p3.stat));
    } else if(Np == 4) {
        double c[4];
        auto &p1 = *itp.first;
        auto &p2 = *(itp.first + 1);
        auto &p3 = *(itp.first + 2);
        auto &p4 = *(itp.first + 3);
        cubic_coeff(p1.x, p2.x, p3.x, p4.x, p.x, &c[0]);
        p.y = c[0]*p1.y + c[1]*p2.y + c[2]*p3.y + c[3]*p4.y;
        p.syst = c[0]*p1.syst + c[1]*p2.syst + c[2]*p3.syst + c[3]*p4.syst;
        p.stat = std::sqrt(pow2(c[0]*p1.stat) + pow2(c[1]*p2.stat) + pow2(c[2]*p3.stat) + pow2(c[3]*p4.stat));
    } else {
        std::cerr << "interpolation error, find " << Np << " points around x = "
                 << p.x << ", while expected 2 ~ 4 points." << std::endl;
        return false;
    }

    return true;
}
