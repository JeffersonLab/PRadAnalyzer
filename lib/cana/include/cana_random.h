#ifndef CANA_RANDOM_H
#define CANA_RANDOM_H

#include <random>
#include <iostream>
#include "cana_interp.h"

namespace cana
{
    // seeding random engine
    template<class T = std::mt19937, std::size_t N = T::state_size>
    auto seeded_random_engine() -> typename std::enable_if<!!N, T>::type
    {
        typename T::result_type random_data[N];
        std::random_device rd;
        std::generate(std::begin(random_data), std::end(random_data), std::ref(rd));
        std::seed_seq seeds(std::begin(random_data), std::end(random_data));
        T re(seeds);
        return re;
    }

    // random number generator
    template<typename T = double, class Engine = std::mt19937>
    class rand_gen
    {
    public:
        // constructor
        rand_gen()
        : engine(seeded_random_engine<Engine>())
        {
            typename Engine::result_type range = engine.max() - engine.min();
            divisor = static_cast<T>(range) + 1;

            if(divisor <= 0) {
                std::cerr << "CRandom: Output type cannot handle engine range, "
                          << "undefined value may be returned."
                          << std::endl;
            }
        }

        // get a random number in (0, 1)
        T operator() ()
        {
            return static_cast<T>(engine() - engine.min())/divisor;
        }

        // get a random number in (min, max)
        T operator() (T min, T max)
        {
            return static_cast<T>(engine() - engine.min())/divisor*(max - min) + min;
        }

    private:
        Engine engine;
        T divisor;
    };

    // simple structure for converting uniform distribution to any distribution
    struct val_cdf
    {
        double val, cdf;

        // constructors
        val_cdf() {}
        val_cdf(double v, double c) : val(v), cdf(c) {}

        // for binary search
        bool operator ==(double v) const {return cdf == v;}
        bool operator <(double v) const {return cdf < v;}
        bool operator >(double v) const {return cdf > v;}
        bool operator != (double v) const {return cdf != v;}
    };

    // convert uniform distribution to another distribution according to the
    // val-cdf distribution provided
    template<class RdmIt, typename T>
    inline T uni2dist(RdmIt begin, RdmIt end, T val)
    {
        auto itv = cana::binary_search_interval(begin, end, val);

        // should not happen
        if(itv.first == end || itv.second == end)
        {
            std::cerr << "Fail to convert distribution from value = " << val << std::endl;
            return 0.;
        }

        if(itv.first == itv.second) {
            return itv.first->val;
        } else {
            return cana::linear_interp(itv.first->cdf, itv.first->val,
                                       itv.second->cdf, itv.second->val,
                                       val);
        }
    }
} // namespace cana

#endif // CANA_RANDOM_H
