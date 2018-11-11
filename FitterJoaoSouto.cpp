/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FitterJoaoSouto.cpp
 * Author: João Vicente Souto
 * 
 * Created on October 13, 2018, 9:10 PM
 */

#include "FitterJoaoSouto.h"

#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <cmath>
#include <vector>

#include "Traits.h"
#include "ProbDistrib.h"

template<typename F, typename ... Args>
double FitterJoaoSouto::square_error(F f, Args ... args)
{
    double error = 0;
    double range = (_stats->max() - _stats->min()) / _stats->histogramNumClasses();
    double i = _stats->min();
    double j = i + range;
    double total_elements = _stats->numElements();

    Traits<Integrator_if>::Implementation integrator;

    for (auto c = 0; c < _stats->histogramNumClasses(); c++, i = j, j += range)
    {
        double area = integrator.integrate(i, j, f, args...);
        error += pow(_stats->histogramClassFrequency(c) - area * total_elements, 2);
    }

    return error;
}

FitterJoaoSouto::FitterJoaoSouto()
{
    _collector = new Traits<Collector_if>::Implementation();
    _stats = new Traits<Statistics_if>::Implementation();
    _stats->setCollector(_collector);
}

FitterJoaoSouto::~FitterJoaoSouto()
{
    //! Não é possível deletar o ponteiro porque a interface não
    //! possui um destrutor virtual implementado!!!
//    delete _stats;
//    delete _collector;
}

bool FitterJoaoSouto::isNormalDistributed(double confidencelevel)
{
    double sqrerror, avg, stddev;
    
    fitNormal(&sqrerror, &avg, &stddev);

    Traits<HypothesisTester_if>::Implementation tester;
    tester.setDataFilename(_collector->getDataFilename());
    
    return tester.testAverage(confidencelevel, avg, HypothesisTester_if::EQUAL)
      && tester.testVariance(confidencelevel, stddev, HypothesisTester_if::EQUAL);
}

void FitterJoaoSouto::fitUniform(double *sqrerror, double *min, double *max)
{
    *min = _stats->min();
    *max = _stats->max();
    *sqrerror = square_error(ProbDistrib::uniform, *min, *max);
}

void FitterJoaoSouto::fitTriangular(double *sqrerror, double *min, double *mo, double *max)
{
    *min = _stats->min();
    *mo  = _stats->mode();
    *max = _stats->max();
    *sqrerror = square_error(ProbDistrib::triangular, *min, *mo, *max);
}

void FitterJoaoSouto::fitNormal(double *sqrerror, double *avg, double *stddev)
{
    *avg      = _stats->average();
    *stddev   = _stats->stddeviation();
    *sqrerror = square_error(ProbDistrib::normal, *avg, *stddev);
}

void FitterJoaoSouto::fitExpo(double *sqrerror, double *avg1)
{
    *avg1 = _stats->average();
    *sqrerror = square_error(ProbDistrib::exponential, *avg1);
}

void FitterJoaoSouto::fitErlang(double *sqrerror, double *avg, int *m)
{
    *avg = _stats->average();
    *m = *avg * std::pow(*avg / _stats->stddeviation(), 2);   
    *sqrerror = square_error(ProbDistrib::erlang, *avg, *m);
}

void FitterJoaoSouto::fitBeta(double *sqrerror, double *alpha, double *beta, double *infLimit, double *supLimit)
{
    *infLimit = _stats->min();
    *supLimit = _stats->max();

    double avg = _stats->average();
    double dev = _stats->stddeviation();

    *alpha = ((1 - avg) / std::pow(dev, 2) - 1 / avg) * std::pow(avg, 2);
    *beta = *alpha * (1 / avg - 1);

    *sqrerror = square_error(ProbDistrib::beta, *alpha, *beta);
}

void FitterJoaoSouto::fitWeibull(double *sqrerror, double *alpha, double *scale)
{
    double avg = _stats->average();
    double dev = _stats->stddeviation();
    
    *alpha = std::pow(dev / avg, -1.086);
    *scale = avg / std::tgamma(1 + 1 / *alpha);
    
    *sqrerror = square_error(ProbDistrib::weibull, *alpha, *scale);
}

void FitterJoaoSouto::fitAll(double *sqrerror, std::string *name)
{
    double error, aux1, aux2, aux3, aux4;
    std::pair<double, std::string> fit;

    fitUniform(&fit.first, &aux1, &aux2);
    fit.second = "Uniform";

    fitTriangular(&error, &aux1, &aux2, &aux3);
    if (error < fit.first) {
        fit.first = error;
        fit.second = "Triangular";
    }

    fitNormal(&error, &aux1, &aux2);
    if (error < fit.first) {
        fit.first = error;
        fit.second = "Normal";
    }

    fitExpo(&error, &aux1);
    if (error < fit.first) {
        fit.first = error;
        fit.second = "Exponencial";
    }

    int auxInt;
    fitErlang(&error, &aux1, &auxInt);
    if (error < fit.first) {
        fit.first = error;
        fit.second = "Erlang";
    }

    fitBeta(&error, &aux1, &aux2, &aux3, &aux4);
    if (error < fit.first) {
        fit.first = error;
        fit.second = "Beta";
    }

    fitWeibull(&error, &aux1, &aux2);
    if (error < fit.first) {
        fit.first = error;
        fit.second = "Weibull";
    }

    *sqrerror = fit.first;
    *name = fit.second;
}

void FitterJoaoSouto::setDataFilename(std::string dataFilename)
{
    _collector->setDataFilename(dataFilename);
}

std::string FitterJoaoSouto::getDataFilename()
{
    return _collector->getDataFilename();
}
