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

    min_cluster_hits = getDefConfig<unsigned int>("Min Cluster Hits", 1, verbose);
    max_cluster_hits = getDefConfig<unsigned int>("Max Cluster Hits", 20, verbose);
    split_cluster_diff = getDefConfig<float>("Split Threshold", 14, verbose);
    cross_talk_width = getDefConfig<float>("Cross Talk Width", 2, verbose);

    // get cross talk characteristic distance
    charac_distance.clear();
    std::string dist_str = GetConfig<std::string>("Characteristic Distance");
    auto dists = ConfigParser::split(dist_str, ",");  // split input string
    while(dists.size())
    {
	    std::string dist = ConfigParser::trim(dists.front(), " \t"); // trim off white spaces
	    dists.pop_front();
	    charac_distance.push_back(std::stod(dist)); // convert string to double and save it
    }
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

    // split the clusters that may contain multiple physical hits
    splitCluster(clusters);

    // reconstruct the cluster position
    reconstructCluster(clusters);

    // remove the clusters that does not pass certain criteria
    filterCluster(clusters);
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
            clusters.emplace_back(std::vector<StripHit>(cluster_begin, next));
            break;
        }

        // recursively group hits for overlapped APVs
        if(next->strip == it->strip)
        {
            std::vector<StripHit> dup1, dup2;

            // duplicate shared strips but halve the charge
            for(auto itd = cluster_begin; itd != it; ++itd)
            {
                StripHit shared_hit(*itd);
                shared_hit.charge /= 2.;

                dup1.push_back(shared_hit);
                dup2.push_back(shared_hit);
            }

            // dispatch the rest strips by APV address
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
        } else if(next->strip - it->strip > 1) {
            clusters.emplace_back(std::vector<StripHit>(cluster_begin, next));
            cluster_begin = next;
        }
    }
}

// split cluster at valley
void PRadGEMCluster::splitCluster(std::vector<StripCluster> &clusters)
const
{
    // We are trying to find the valley that satisfies certain critieria,
    // i.e., less than a sigma comparing to neighbor strips on both sides.
    // Such valley probably means two clusters are grouped together, so it
    // will be separated, and each gets 1/2 of the charge from the overlap
    // strip.

    // loop over the cluster list
    // since the vector size is changing, it cannot use iterator
    for(size_t i = 0; i < clusters.size(); ++i)
    {
        // no need to do separation if less than 3 hits
        if(clusters[i].hits.size() < 3)
            continue;

        // new cluster for the latter part after split
        StripCluster split_cluster;

        // insert the splited cluster if there is one
        if(splitCluster_sub(clusters[i], split_cluster)) {
            // insert keeps the original order, but has a worse performance
            clusters.insert(clusters.begin()+i+1, split_cluster);
        }
    }
}

// This function helps splitCluster
// It finds the FIRST local minimum, separate the cluster at its position
// The charge at local minimum strip will be halved, and kept for both original
// and split clusters.
// It returns true if a cluster is split, and vice versa
// The split part of the original cluster c will be removed, and filled in c1
bool PRadGEMCluster::splitCluster_sub(StripCluster &c, StripCluster &c1)
const
{
    // we use 2 consecutive iterator
    auto it = c.hits.begin();
    auto it_next = it + 1;

    // loop to find the local minimum
    bool descending = false, extremum = false;
    auto minimum = it;
    for(; it_next != c.hits.end(); ++it, ++it_next)
    {
        if(descending) {
            // update minimum
            if(it->charge < minimum->charge)
                minimum = it;

            // transcending trend, confirm a local minimum (valley)
            if(it_next->charge - it->charge > split_cluster_diff) {
                extremum = true;
                // only needs the first local minimum, thus exit the loop
                break;
            }
        } else {
            // descending trend, expect a local minimum
            if(it->charge - it_next->charge > split_cluster_diff) {
                descending = true;
                minimum = it_next;
            }
        }
    }

    if(extremum) {
        // half the charge of overlap strip
        minimum->charge /= 2.;

        // new split cluster
        c1 = StripCluster(std::vector<StripHit>(minimum, c.hits.end()));

        // remove the hits that are moved into new cluster, but keep the minimum
        c.hits.erase(std::next(minimum, 1), c.hits.end());
    }

    return extremum;
}

// filter out bad clusters
#define MAX_CLUSTER_WIDTH 2.0
void PRadGEMCluster::filterCluster(std::vector<StripCluster> &clusters)
const
{
    // remove cluster that has too less/many hits
    for(auto it = clusters.begin(); it != clusters.end(); ++it)
    {
        if((it->hits.size() < min_cluster_hits) ||
           (it->hits.size() > max_cluster_hits))
            clusters.erase(it--);
    }

    // remove cross talk cluster
    for(auto it = clusters.begin(); it != clusters.end(); ++it)
    {
	    if(filterCrossTalk(*it, clusters)) {
            clusters.erase(it--);
        }
    }

}

bool PRadGEMCluster::filterCrossTalk(const StripCluster &cluster,
                                     const std::vector<StripCluster> &clusters)
const
{
    if(!cluster.IsCrossTalk())
        return false;

    for(auto it = clusters.begin(); it != clusters.end(); ++it)
    {
        double delta = fabs(it->position - cluster.position);

        for(auto &dist : charac_distance)
        {
            if(delta > dist - cross_talk_width &&
               delta < dist + cross_talk_width)
                return true;
        }
    }

    return false;
}

// calculate the cluster position
// it reconstruct the position of cluster using linear weight of charge portion
void PRadGEMCluster::reconstructCluster(std::vector<StripCluster> &clusters)
const
{
    for(auto &c : clusters)
    {
        // here determine position, peak charge and total charge of the cluster
        c.total_charge = 0.;
        c.peak_charge = 0.;
        float weight_pos = 0.;

        // no hits
        if(!c.hits.size())
            continue;

        for(auto &hit : c.hits)
        {
            if(c.peak_charge < hit.charge)
                c.peak_charge = hit.charge;

            c.total_charge += hit.charge;

            weight_pos +=  hit.position*hit.charge;
        }

        c.position = weight_pos/c.total_charge;
    }
}

// this function accepts x, y clusters from detectors and then form GEM Cluster
// it return the number of clusters
void PRadGEMCluster::CartesianReconstruct(const std::vector<StripCluster> &x_cluster,
                                          const std::vector<StripCluster> &y_cluster,
                                          std::vector<GEMHit> &container,
                                          int det_id)
const
{
    // empty first
    container.clear();

    // TODO, probably add some criteria here to filter out some bad clusters
    // fill possible clusters in
    for(auto &xc : x_cluster)
    {
        for(auto &yc : y_cluster)
        {
            container.emplace_back(xc.position, yc.position, 0.,        // by default z = 0
                                   det_id,                              // detector id
                                   xc.total_charge, yc.total_charge,    // fill in total charge
                                   xc.peak_charge, yc.peak_charge,      // fill in peak charge
                                   xc.hits.size(), yc.hits.size());     // number of hits
        }
    }
}
