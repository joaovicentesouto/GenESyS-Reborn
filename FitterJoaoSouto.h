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
    template<typename F, typename ... Args>
    double square_error(F f, Args ... args);

private:
    Statistics_if * _stats{nullptr};
    CollectorDatafile_if * _collector{nullptr};
};

#endif /* FITTERJOAOSOUTO_H */
