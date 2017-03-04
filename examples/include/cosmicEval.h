#include "PRadHyCalSystem.h"
#include <vector>

#define DEFAULT_ENERGY_THRES 1.0

std::vector<PRadHyCalModule*> __modules;
std::vector<std::vector<PRadHyCalModule*>> __groups;

struct LineParam
{
    bool success;
    double k, b, rsq, chisq;

    LineParam()
    : success(false), rsq(0.), chisq(10.) {};
};

struct EnergyParam
{
    double total;
    double maximum;
    double uniform;

    EnergyParam()
    : total(0.), maximum(0.), uniform(10.)
    {};
};

struct EventParam
{
    EnergyParam group_energy;
    LineParam group_line;
    unsigned int group_size;
    unsigned int max_group_size;
    EnergyParam max_group_energy;
    LineParam max_group_line;

    EventParam()
    : group_size(0), max_group_size(0) {};

    std::vector<double> GetParamList()
    const
    {
        std::vector<double> params = {(double)group_size,
                                      (double)max_group_size,
                                      group_energy.maximum,
                                      max_group_energy.uniform,
                                      max_group_line.rsq,
                                      max_group_line.chisq};

        return params;
    }
};

void GetFiredModules(const PRadHyCalDetector *det, std::vector<PRadHyCalModule*> &m, double thres);
void GroupHits(std::vector<PRadHyCalModule*> &hits, std::vector<std::vector<PRadHyCalModule*>> &groups);
LineParam LinearRegression(const std::vector<PRadHyCalModule*> &group);
EnergyParam EvalEnergy(const std::vector<PRadHyCalModule*> &group);

EventParam AnalyzeEvent(PRadHyCalSystem *sys, const EventData &event, double thres = DEFAULT_ENERGY_THRES)
{
    EventParam res;
    sys->ChooseEvent(event);
    GetFiredModules(sys->GetDetector(), __modules, thres);
    res.group_line = LinearRegression(__modules);
    res.group_energy = EvalEnergy(__modules);

    GroupHits(__modules, __groups);

    if(__groups.empty())
        return res;

    res.group_size = __groups.size();

    // find the group that has the most hits
    unsigned int max = 0;
    unsigned int i_max = 0;
    for(unsigned int i = 0; i < __groups.size(); ++i)
    {
        if(__groups.at(i).size() > max)
        {
            max = __groups.at(i).size();
            i_max = i;
        }
    }

    auto &m_group = __groups.at(i_max);
    res.max_group_size = m_group.size();

    res.max_group_line = LinearRegression(m_group);
    res.max_group_energy = EvalEnergy(m_group);

    return res;
}

// fill fired modules with energy larger than threshold
void GetFiredModules(const PRadHyCalDetector *det, std::vector<PRadHyCalModule*> &m, double thres)
{
    m.clear();

    for(auto &module : det->GetModuleList())
    {
        double e = module->GetEnergy();
        if(e > thres)
        {
            m.push_back(module);
        }
    }
}

// get the hit distance normalized by module size
inline double HitDistance(const PRadHyCalModule *m1, const PRadHyCalModule *m2)
{
    double dx = (m1->GetX() - m2->GetX())/(m1->GetSizeX() + m2->GetSizeX());
    double dy = (m1->GetY() - m2->GetY())/(m1->GetSizeY() + m2->GetSizeY());

    return sqrt(dx*dx + dy*dy)*2.;
}

// check if two hits are adjacent
inline bool IsAdjacent(const PRadHyCalModule *m1, const PRadHyCalModule *m2)
{
    return HitDistance(m1, m2) < 1.5;
}

// check if two clusters are adjacent
inline bool IsAdjacent(const std::vector<PRadHyCalModule*> &g1, const std::vector<PRadHyCalModule*> &g2)
{
    for(auto &m1 : g1)
    {
        for(auto &m2 : g2)
        {
            if(IsAdjacent(m1, m2)) {
                return true;
            }
        }
    }
    return false;
}

// fill hit into hits groups
bool FillGroup(PRadHyCalModule *hit, std::vector<std::vector<PRadHyCalModule*>> &groups)
{
    for(auto &group : groups)
    {
        for(auto &prev_hit : group)
        {
            // it belongs to a existing cluster
            if(IsAdjacent(hit, prev_hit)) {
                group.push_back(hit);
                return true;
            }
        }
    }
    return false;
}

// group the adjacent hits
void GroupHits(std::vector<PRadHyCalModule*> &hits, std::vector<std::vector<PRadHyCalModule*>> &groups)
{
    groups.clear();

    // roughly combine all adjacent hits
    for(auto &hit : hits)
    {
        // not belong to any existing cluster
        if(!FillGroup(hit, groups)) {
            std::vector<PRadHyCalModule*> new_group;
            new_group.reserve(50);
            new_group.push_back(hit);
            groups.emplace_back(std::move(new_group));
        }
    }

    // merge adjacent groups
    for(auto it = groups.begin(); it != groups.end(); ++it)
    {
        auto it_next = it;
        while(++it_next != groups.end())
        {
            if(IsAdjacent(*it, *it_next)) {
                it_next->insert(it_next->end(), it->begin(), it->end());
                groups.erase(it--);
                break;
            }
        }
    }
}

LineParam LinearRegression(const std::vector<PRadHyCalModule*> &group)
{
    LineParam res;

    if(group.size() < 3)
        return res;

    // get average x and y of these hits
    double ave_x = 0., ave_y = 0.;

    for(auto &hit : group)
    {
        ave_x += hit->GetX();
        ave_y += hit->GetY();
    }

    ave_x /= (double)group.size();
    ave_y /= (double)group.size();

    // fit the line
    double denom = 0., numer = 0.;

    for(auto &hit : group)
    {
        denom += (hit->GetX() - ave_x)*(hit->GetX() - ave_x);
        numer += (hit->GetX() - ave_x)*(hit->GetY() - ave_y);
    }

    if(denom == 0.)
        return res;

    res.k = numer/denom;
    res.b = ave_y - res.k*ave_x;

    // R-square and chi-square test
    denom = 0.;
    numer = 0.;
    res.chisq = 0.;
    for(auto &hit : group)
    {
        double yp = hit->GetX()*res.k + res.b;
        numer += (yp - ave_y)*(yp - ave_y);
        denom += (hit->GetY() - ave_y)*(hit->GetY() - ave_y);
        res.chisq += ((hit->GetY() - yp)*(hit->GetY() - yp))/hit->GetSizeY()/hit->GetSizeY();
    }

    res.chisq /= (double)(group.size() - 2);

    if(denom > 0.)
        res.rsq = numer/denom;

    return res;
}

EnergyParam EvalEnergy(const std::vector<PRadHyCalModule*> &group)
{
    EnergyParam res;

    res.total = 0.;
    res.maximum = 0.;
    for(auto &hit : group)
    {
        double e = hit->GetEnergy();
        if(e > res.maximum)
            res.maximum = e;

        res.total += e;
    }

    double average = res.total/(double)group.size();


    if(group.size() < 2)
        return res;


    res.uniform = 0.;
    for(auto &hit : group)
    {
        res.uniform += fabs(hit->GetEnergy() - average)/average;
    }
    res.uniform /= (double)(group.size() - 1);

    return res;
}


