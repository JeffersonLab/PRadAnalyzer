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
    std::string dist_str = GetConfig<std::string>("Characteristic Distance");
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

typedef std::vector<StripHit>::iterator SHit;
void cluster_hits(SHit beg, SHit end, double thres, std::vector<StripCluster> &clusters)
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
        cluster_hits(minimum, end, thres, clusters);
    } else {
        clusters.emplace_back(std::vector<StripHit>(beg, end));
    }
}


// used for separate hits of overlapped APVs
#define IS_FROM_APV_SET1(hit) ( ((hit).apv_addr == APVAddress(1, 8)) || \
                                ((hit).apv_addr == APVAddress(6, 8)) )

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

    // group the hits that have consecutive strip number
    auto cluster_begin = hits.begin();
    for(auto it = hits.begin(); it != hits.end(); ++it)
    {
        auto next = it + 1;

        // end of list, group the last cluster
        if(next == hits.end()) {
            cluster_hits(cluster_begin, next, split_cluster_diff, clusters);
            break;
        }

        // recursively group hits for overlapped APVs
        if(next->strip == it->strip) {
            std::vector<StripHit> dup1, dup2;
            dup1.reserve(hits.size());
            dup2.reserve(hits.size());

            // duplicate shared strips but halve the charge
            for(auto itd = cluster_begin; itd != it; ++itd)
            {
                StripHit shared_hit(*itd);
                shared_hit.charge /= 2.;

                dup1.push_back(shared_hit);
                dup2.push_back(shared_hit);
            }

            // dispatch the rest strips by APV address
            // duplicating always occurs at the end
            for(auto itd = it; itd != hits.end(); ++itd)
            {
                if(IS_FROM_APV_SET1(*itd))
                    dup1.push_back(*itd);
                else
                    dup2.push_back(*itd);
            }

            // group the remaining strips into clusters
            groupHits(dup1, clusters);
            groupHits(dup2, clusters);

            break; // finished all grouping, break out the entire loop

        // not consecutive, create a new cluster
        } else if(next->strip - it->strip > (int)consecutive_thres) {
            cluster_hits(cluster_begin, next, split_cluster_diff, clusters);
            cluster_begin = next;
        }
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
