/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FitterDefault.cpp
 * Author: João Vicente Souto
 * 
 * Created on 19 de Setembro de 2018, 16:00
 */

#include "FitterDefault.h"

void FitterDefault::readFile() {
    delete _stats;
    _stats = new Statistics(_dataFilename);
}

LocalStatistics::LocalStatistics(std::string dataFilename)
{
    std::ifstream input(dataFilename, std::ios::in);

    if (!input.is_open())
        throw std::out_of_range("Could not open file!");

    bool init = false;
    double num = 0;

    if (input >> num) {
        _min = _max = num;
        _values.push_back(num);
        init = true;
    }

    while (input >> num) {
        if (num < _min)
            _min = num;
        else if (_max < num)
            _max = num;

        _values.push_back(num);
    }

    input.close();

    if (!init)
        throw std::out_of_range("Values does not exist!");

    _values_amount = _values.size();
    _class_amount = std::ceil(1 + 3.322 * (log10(_values_amount)));
    _range = (_max - _min) / _class_amount;
}

unsigned int LocalStatistics::histogramClassFrequency(unsigned short classNum)
{
    double min = histogramClassLowerLimit(classNum);
    double max = min + _range;
    unsigned int amount = 0;

    for (double value : _values)
        if (min <= value && value < max)
            amount++;
        
    return amount;
}

double LocalStatistics::average() {
    if (!_values_amount)
        return 0;

    double total = 0;
    for (double value : _values)
            total += value;
    
    return total / _values_amount;
}

double LocalStatistics::stddeviation() {
    if (_values_amount < 2)
        throw std::out_of_range("StdDev: It takes at least two values!");

    double value_mean = average();

    double total = 0;
    for (double value : _values)
            total += pow(value - value_mean, 2);

    return sqrt(total / (_values_amount - 1));
}

double LocalStatistics::mode() {
    int clss = 0, amount = histogramClassFrequency(0);

    for (int i = 0; i < _class_amount; i++)
        if (amount < histogramClassFrequency(i))
        {
            amount = histogramClassFrequency(i);
            clss = i;
        }
    
    double limit_inf = histogramClassLowerLimit(clss);
    return (limit_inf + (limit_inf + _range)) / 2;
}

bool FitterDefault::isNormalDistributed(double confidencelevel) {
    std::ofstream output("./NormalDistributionSample.txt", std::ios::out);

    if (!output.is_open())
        throw std::out_of_range("Could not open file!");

	Traits<Sampler_if>::Implementation sampler;

    for (int i = 0; i < 1000; i++)
		output << sampler.sampleUniform(10, 20) << std::endl;

	output.close();

    Traits<HypothesisTester_if>::Implementation tester;
    tester.setDataFilename("./NormalDistributionSample.txt");

    return tester.testAverage(confidencelevel, _dataFilename, HypothesisTester_if::EQUAL)
        && tester.testVariance(confidencelevel, _dataFilename, HypothesisTester_if::EQUAL);
}

void FitterDefault::fitUniform(double *sqrerror, double *min, double *max) {
	if (!_stats)
        readFile();

    *min = _stats->min();
    *max = _stats->max();
    *sqrerror = square_error(ProbDistrib::uniform, *min, *max);
}

void FitterDefault::fitTriangular(double *sqrerror, double *min, double *mo, double *max) {
	if (!_stats)
        readFile();

    *min = _stats->min();
    *mo  = _stats->mode();
    *max = _stats->max();
    *sqrerror = square_error(ProbDistrib::triangular, *min, *mo, *max);
}

void FitterDefault::fitNormal(double *sqrerror, double *avg, double *stddev) {
	if (!_stats)
        readFile();

    *avg      = _stats->average();
    *stddev   = _stats->stddeviation();
    *sqrerror = square_error(ProbDistrib::normal, *avg, *stddev);
}

void FitterDefault::fitExpo(double *sqrerror, double *avg1) {
    if (!_stats)
        readFile();

    *avg1 = _stats->average();
    *sqrerror = square_error(ProbDistrib::exponential, *avg1);
}

void FitterDefault::fitErlang(double *sqrerror, double *a, double *b, double *offset, double *mult) {
    if (!_stats)
        readFile();

    *a = _stats->min();;
    *b = _stats->max();
    *offset = _stats->stddeviation();
    *mult = _stats->average();
    *sqrerror = square_error(FitterDefault::erlang, *offset, *mult);
}

void FitterDefault::fitBeta(double *sqrerror, double *a, double *b, double *offset, double *mult) {
    if (!_stats)
        readFile();

    *a = _stats->min();
    *b = _stats->max();

    double mean = _stats->average();
    double std = _stats->stddeviation();

    //! https://en.wikipedia.org/wiki/Beta_distribution
    *offset = pow(mean, 2) * (((1 - mean) / pow(std, 2)) - (1 / mean));
    *mult = *offset * ( 1 / mean - 2);

    *sqrerror = square_error(ProbDistrib::beta, *offset, *mult);
}

void FitterDefault::fitWeibull(double *sqrerror, double *a, double *b, double *offset, double *mult) {
    if (!_stats)
        readFile();

    *a = _stats->min();
    *b = _stats->max();
    *offset = _stats->stddeviation();
    *mult = _stats->average();
    *sqrerror = square_error(ProbDistrib::weibull, *offset, *mult);
}

void FitterDefault::fitAll(double *sqrerror, std::string *name) {
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

    fitErlang(&error, &aux1, &aux2, &aux3, &aux4);
    if (error < fit.first) {
        fit.first = error;
        fit.second = "Erlang";
    }

    fitBeta(&error, &aux1, &aux2, &aux3, &aux4);
    if (error < fit.first) {
        fit.first = error;
        fit.second = "Beta";
    }

    fitWeibull(&error, &aux1, &aux2, &aux3, &aux4);
    if (error < fit.first) {
        fit.first = error;
        fit.second = "Weibull";
    }

    *sqrerror = fit.first;
    *name = fit.second;
}

int FitterDefault::factorial(int n) {
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double FitterDefault::erlang(double x, double mean, double M) {
    return pow(x, M - 1) * exp(- x / mean) / (pow(mean, M) * factorial(M - 1));
}

void FitterDefault::setDataFilename(std::string dataFilename) {
    delete _stats;
    _stats = nullptr;
    _dataFilename = dataFilename;
}

std::string FitterDefault::getDataFilename() {
    return _dataFilename;
}
