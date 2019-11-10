//============================================================================//
// Basic PRad Cluster Reconstruction Class For GEM                            //
// GEM Planes send hits infromation and container for the GEM Clusters to this//
// class reconstruct the clusters and send it back to GEM Planes.             //
// Thus the clustering algorithm can be adjusted in this class.               //
//                                                                            //
// Xinzhan Bai & Kondo Gnanvo, first version coding of the algorithm          //
// Chao Peng, adapted to PRad analysis software package                       //
// Xinzhan Bai, add cross talk removal                                        //
// 10/21/2016                                                                 //
//============================================================================//

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <deque>
#include <cmath>
#include "PRadGEMCluster.h"
#include "PRadGEMDetector.h"

// constructor
PRadGEMCluster::PRadGEMCluster(const std::string &config_path)
{
    Configure(config_path);
}

// destructor
PRadGEMCluster::~PRadGEMCluster()
{
    // place holder
}

// configure the cluster method
void PRadGEMCluster::Configure(const std::string &path)
{
    // if no configuration file specified, load the default value quietly
    bool verbose = false;

    if(!path.empty()) {
        ConfigObject::Configure(path);
        verbose = true;
    }

    CONF_CONN(min_cluster_hits, "Min Cluster Hits", 1, verbose);
    CONF_CONN(max_cluster_hits, "Max Cluster Hits", 20, verbose);
    CONF_CONN(split_cluster_diff, "Split Threshold", 14, verbose);
    CONF_CONN(cross_talk_width, "Cross Talk Width", 2, verbose);
    CONF_CONN(consecutive_thres, "Consecutive Threshold", 1, verbose);

    // get cross talk characteristic distance
    charac_dists.clear();
    std::string dist_str = GetConfigValue<std::string>("Characteristic Distance");
    charac_dists = ConfigParser::stofs(dist_str, ",", " \t");
}

// group hits into clusters
void PRadGEMCluster::FormClusters(std::vector<StripHit> &hits,
                                  std::vector<StripCluster> &clusters)
const
{
    // clean container first
    clusters.clear();

    // group consecutive hits as the preliminary clusters
    groupHits(hits, clusters);

    // reconstruct the cluster position
    for(auto &cluster : clusters)
    {
        reconstructCluster(cluster);
    }

    // set cross talk flag
    setCrossTalk(clusters);
}

// a helper function to further separate hits at minimum
template<class Iter>
void split_cluster(Iter beg, Iter end, double thres, std::vector<StripCluster> &clusters)
{
    auto size = end - beg;
    if(size < 3) {
        clusters.emplace_back(std::vector<StripHit>(beg, end));
        return;
    }

    // find the first local minimum
    bool descending = false, extremum = false;
    auto minimum = beg;
    for(auto it = beg, it_n = beg + 1; it_n != end; ++it, ++it_n)
    {
        if(descending) {
            // update minimum
            if(it->charge < minimum->charge)
                minimum = it;

            // transcending trend, confirm a local minimum (valley)
            if(it_n->charge - it->charge > thres) {
                extremum = true;
                // only needs the first local minimum, thus exit the loop
                break;
            }
        } else {
            // descending trend, expect a local minimum
            if(it->charge - it_n->charge > thres) {
                descending = true;
                minimum = it_n;
            }
        }
    }

    if(extremum) {
        // half the charge of overlap strip
        minimum->charge /= 2.;

        // new split cluster
        clusters.emplace_back(std::vector<StripHit>(beg, minimum));

        // check the leftover strips
        split_cluster(minimum, end, thres, clusters);
    } else {
        clusters.emplace_back(std::vector<StripHit>(beg, end));
    }
}

// cluster consecutive hits
template<class Iter>
inline void cluster_hits(Iter beg, Iter end, int con_thres, double diff_thres, std::vector<StripCluster> &clusters)
{
    auto cbeg = beg;
    for(auto it = beg; it != end; ++it)
    {
        auto it_n = it + 1;
        if((it_n == end) || (it_n->strip - it->strip > con_thres)) {
            split_cluster(cbeg, it_n, diff_thres, clusters);
            cbeg = it_n;
        }
    }
}

// used for separate hits of overlapped APVs
#define IS_FROM_APV_SET1(hit) ( ((hit).apv_addr == APVAddress(1, 8)) || \
                                ((hit).apv_addr == APVAddress(6, 8)) )

template<class Iter>
inline bool is_duplicated_strip(Iter it)
{
    // strip number range for duplicates
    if( ((it->strip < 1392) && (it->strip > 1263)) &&
    // possible apv that connect with duplicated strips
        ((it->apv_addr == APVAddress(1, 8)) ||
         (it->apv_addr == APVAddress(0, 8)) ||
         (it->apv_addr == APVAddress(3, 8)) ||
         (it->apv_addr == APVAddress(6, 8)) ||
         (it->apv_addr == APVAddress(5, 8)) ||
         (it->apv_addr == APVAddress(7, 8))) )
        return true;
    return false;
}

template<class Iter, template<class, class> class Container>
inline void separate_duplicates(Iter beg, Iter end,
                                Container<StripHit, std::allocator<StripHit>> &dup1,
                                Container<StripHit, std::allocator<StripHit>> &dup2)
{
    for(auto it = beg; it != end; ++it)
    {
        if(IS_FROM_APV_SET1(*it))
            dup1.emplace_back(*it);
        else
            dup2.emplace_back(*it);
    }
}

// group consecutive hits
void PRadGEMCluster::groupHits(std::vector<StripHit> &hits,
                               std::vector<StripCluster> &clusters)
const
{
    // sort the hits by its strip number
    std::sort(hits.begin(), hits.end(),
              // lamda expr, compare hit by their strip numbers
              [](const StripHit &h1, const StripHit &h2)
              {
                  return h1.strip < h2.strip;
              });

    // for X-plane, we have strips at the same x-position (same strip number),
    // but they are segmented due to the central hole, this needs special treatment
    // separate hits into normal group, upper group, lower group
    auto normal_end = hits.begin();
    // find the end of the normal group
    while((normal_end != hits.end()) && !is_duplicated_strip(normal_end)) {normal_end ++;}
    // separate the duplicated strips into two groups (upper and lower)
    std::deque<StripHit> dup1, dup2;
    separate_duplicates(normal_end, hits.end(), dup1, dup2);

    // no need the special treatment
    if(dup1.empty() || dup2.empty()) {
        cluster_hits(hits.begin(), hits.end(), consecutive_thres, split_cluster_diff, clusters);
    // no normal hits
    } else if(normal_end == hits.begin()) {
        cluster_hits(dup1.begin(), dup1.end(), consecutive_thres, split_cluster_diff, clusters);
        cluster_hits(dup2.begin(), dup2.end(), consecutive_thres, split_cluster_diff, clusters);
    } else {
        // group normal first
        cluster_hits(hits.begin(), normal_end, consecutive_thres, split_cluster_diff, clusters);
        // check if the last hit is consecutive to the duplicated groups
        auto lastn = normal_end - 1;
        bool neighbor1 = (dup1.begin()->strip - lastn->strip > (int)consecutive_thres);
        bool neighbor2 = (dup2.begin()->strip - lastn->strip > (int)consecutive_thres);
        // last hits
        auto &last_hits = clusters.back().hits;
        // share last hits
        if(neighbor1 & neighbor2) {
            // determine share by the closest hit charge
            double factor1 = 1./(dup2.begin()->charge/dup1.begin()->charge + 1.);
            for(auto it = last_hits.rbegin(); it != last_hits.rend(); ++it)
            {
                auto share_hit = *it;
                share_hit.charge *= factor1;
                dup1.emplace_front(share_hit);
            }
            double factor2 = 1. - factor1;
            for(auto it = last_hits.rbegin(); it != last_hits.rend(); ++it)
            {
                auto share_hit = *it;
                share_hit.charge *= factor2;
                dup2.emplace_front(share_hit);
            }
        // merge to dup1
        } else if (neighbor1) {
            dup1.insert(dup1.begin(), last_hits.begin(), last_hits.end());
        // merge to dup2
        } else if(neighbor2) {
            dup2.insert(dup2.begin(), last_hits.begin(), last_hits.end());
        }

        // discard last cluster, since it is merged into duplicates groups
        if(neighbor1 | neighbor2)
            clusters.pop_back();

        // cluster the duplicates groups
        cluster_hits(dup1.begin(), dup1.end(), consecutive_thres, split_cluster_diff, clusters);
        cluster_hits(dup2.begin(), dup2.end(), consecutive_thres, split_cluster_diff, clusters);
    }
}

// helper function to check cross talk strips
inline bool is_pure_ct(const StripCluster &cl)
{
    for(auto &hit : cl.hits)
    {
        // still has non-cross-talk strips
        if(!hit.cross_talk)
            return false;
    }

    // pure cross talk strips
    return true;
}

// helper function to check cross talk characteristic distance
typedef std::vector<StripCluster>::iterator SCit;
inline bool ct_distance(SCit it, SCit end, float width, const std::vector<float> &charac)
{
    if(it == end)
        return false;

    for(auto itn = it + 1; itn != end; ++itn)
    {
        float delta = std::abs(it->position - itn->position);

        for(auto &dist : charac)
        {
            if((delta > dist - width) && (delta < dist + width))
                return true;
        }
    }

    return false;
}

void PRadGEMCluster::setCrossTalk(std::vector<StripCluster> &clusters)
const
{
    // sort by peak charge
    std::sort(clusters.begin(), clusters.end(),
              [](const StripCluster &c1, const StripCluster &c2)
              {
                  return c1.peak_charge < c2.peak_charge;
              });

    for(auto it = clusters.begin(); it != clusters.end(); ++it)
    {
        // only remove cross talk clusters that is a composite of cross talk strips
        if(!is_pure_ct(*it))
            continue;

        it->cross_talk = ct_distance(it, clusters.end(), cross_talk_width, charac_dists);
    }
}

// calculate the cluster position
// it reconstruct the position of cluster using linear weight of charge fraction
void PRadGEMCluster::reconstructCluster(StripCluster &cluster)
const
{
    // no hits
    if(cluster.hits.empty())
        return;

    // determine position, peak charge and total charge of the cluster
    cluster.total_charge = 0.;
    cluster.peak_charge = 0.;
    float weight_pos = 0.;

    for(auto &hit : cluster.hits)
    {
        if(cluster.peak_charge < hit.charge)
            cluster.peak_charge = hit.charge;

        cluster.total_charge += hit.charge;
        weight_pos +=  hit.position*hit.charge;
    }

    cluster.position = weight_pos/cluster.total_charge;
}


// is it a good cluster
bool PRadGEMCluster::IsGoodCluster(const StripCluster &cluster)
const
{
    // bad size
    if((cluster.hits.size() < min_cluster_hits) ||
       (cluster.hits.size() > max_cluster_hits))
        return false;

    // not a cross talk cluster
    return !cluster.cross_talk;
}

// this function accepts x, y clusters from detectors and then form GEM Cluster
// it return the number of clusters
void PRadGEMCluster::CartesianReconstruct(const std::vector<StripCluster> &x_cluster,
                                          const std::vector<StripCluster> &y_cluster,
                                          std::vector<GEMHit> &container,
                                          int det_id,
                                          float res)
const
{
    // empty first
    container.clear();

    // TODO, probably add some criteria here to filter out some bad clusters
    // fill possible clusters in
    for(auto &xc : x_cluster)
    {
        if(!IsGoodCluster(xc))
            continue;
        for(auto &yc : y_cluster)
        {
            if(!IsGoodCluster(yc))
                continue;

            container.emplace_back(xc.position, yc.position, 0.,        // by default z = 0
                                   det_id,                              // detector id
                                   xc.total_charge, yc.total_charge,    // fill in total charge
                                   xc.peak_charge, yc.peak_charge,      // fill in peak charge
                                   xc.hits.size(), yc.hits.size(),      // number of hits
                                   res);                                // position resolution
        }
    }
}
