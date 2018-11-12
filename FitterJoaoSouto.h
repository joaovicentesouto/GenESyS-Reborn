/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FitterJoaoSouto.h
 * Author: joaovicentesouto
 *
 * Created on October 13, 2018, 9:10 PM
 */

#ifndef FITTERJOAOSOUTO_H
#define FITTERJOAOSOUTO_H

#include "Fitter_if.h"
#include "Statistics_if.h"
#include "CollectorDatafile_if.h"

class FitterJoaoSouto : public Fitter_if {
public:
    FitterJoaoSouto();
    FitterJoaoSouto(const FitterJoaoSouto& orig) = default;
    ~FitterJoaoSouto();
public:
    bool isNormalDistributed(double confidencelevel);
    void fitUniform(double *sqrerror, double *min, double *max);
    void fitTriangular(double *sqrerror, double *min, double *mo, double *max);
    void fitNormal(double *sqrerror, double *avg, double *stddev);
    void fitExpo(double *sqrerror, double *avg1);
    void fitErlang(double *sqrerror, double *avg, int *m);
    void fitBeta(double *sqrerror, double *alpha, double *beta, double *infLimit, double *supLimit);
    void fitWeibull(double *sqrerror, double *alpha, double *scale);
    void fitAll(double *sqrerror, std::string *name);
public:
    void setDataFilename(std::string dataFilename);
    std::string getDataFilename();

private:
    template<typename Integrator, typename F, typename ... Args>
    double square_error(F f, Args ... args);

private:
    Statistics_if * _stats{nullptr};
    CollectorDatafile_if * _collector{nullptr};
};

#include <cmath>

template<typename Integrator, typename F, typename ... Args>
double FitterJoaoSouto::square_error(F f, Args ... args)
{
    double error = 0;
    double range = (_stats->max() - _stats->min()) / _stats->histogramNumClasses();
    double i = _stats->min();
    double j = i + range;
    double total_elements = _stats->numElements();

    Integrator integrator;

    for (auto c = 0; c < _stats->histogramNumClasses(); c++, i = j, j += range)
    {
        double observed_frequency = _stats->histogramClassFrequency(c);
        double expected_frequency = integrator.integrate(i, j, f, args...) * total_elements;
        
        error += pow(observed_frequency - expected_frequency, 2) / expected_frequency;
    }

    return error;
}

#endif /* FITTERJOAOSOUTO_H */
