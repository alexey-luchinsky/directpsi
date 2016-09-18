/* 
 * File:   Random.h
 * Author: luchinsky
 *
 * Created on March 24, 2016, 9:16 AM
 */

#ifndef GGGPSIPSI_RANDOM_H
#define GGGPSIPSI_RANDOM_H

#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stdlib.h>


class Random {
private:
    unsigned int seed;
    boost::mt11213b inner_random;
    boost::uniform_real<double> flat;

public:
    Random() { }

    Random(long _seed_) {
        seed = _seed_;
        if(seed==0) {
            inner_random.seed(time(NULL));
            seed=(int)rand(0,INT_MAX);
        }
        inner_random.seed(seed);      
    }

    unsigned int get_seed() { return seed; };

    double rand(double min, double max) {
        return min + (max - min) * flat(inner_random);
    }
};

extern Random random_generator;


#endif //GGGPSIPSI_RANDOM_H