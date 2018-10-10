/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FitterDefault.h
 * Author: João Vicente Souto
 *
 * Created on 23 de Agosto de 2018, 15:36
 */

#ifndef FITTERDEFAULT_h
#define FITTERDEFAULT_h

#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <cmath>
#include <vector>

#include "Fitter_if.h"
#include "Traits.h"
#include "ProbDistrib.h"

class LocalStatistics
{
public:
	LocalStatistics(std::string dataFilename);
	~LocalStatistics() = default;

	unsigned int numElements() { return _values_amount; }
	double min() { return _min; }
	double max() { return _max; }
	double average();
	double mode(); //! Forced mode
	double stddeviation();
	unsigned short histogramNumClasses() { return _class_amount; }
	double histogramClassLowerLimit(unsigned short classNum) { return _min + classNum * _range; }
	unsigned int histogramClassFrequency(unsigned short classNum);

private:
	double _min, _max, _range;
	size_t _values_amount;
	size_t _class_amount;
	std::vector<double> _values;
};

class FitterDefault: public Fitter_if
{
	//! Pode ser alterado para Statistics_if se ela ler de um arquivo!
	typedef LocalStatistics Statistics;

public:
	FitterDefault() = default;
	FitterDefault(const FitterDefault& orig) = default;
	~FitterDefault() = default;
public:
	bool isNormalDistributed(double confidencelevel);
	void fitUniform (double *sqrerror, double *min, double *max);
	void fitTriangular (double *sqrerror, double *min, double *mo, double *max);
	void fitNormal (double *sqrerror, double *avg, double *stddev);
	void fitExpo (double *sqrerror, double *avg1);
	void fitErlang (double *sqrerror, double *a, double *b, double *offset,double *mult);
	void fitBeta (double *sqrerror, double *a, double *b, double *offset,double *mult);
	void fitWeibull (double *sqrerror, double *a, double *b, double *offset, double *mult);
	void fitAll (double *sqrerror, std::string *name);
public:
	void setDataFilename(std::string dataFilename);
	std::string getDataFilename();

private:
	void readFile();

	template<typename F, typename ... Args>
	double square_error(F f, Args ... args);
	static int factorial(int n);
	static double erlang(double x, double mean, double M);
	
private:
	std::string _dataFilename = "";
	Statistics * _stats{nullptr};
};

template<typename F, typename ... Args>
double FitterDefault::square_error(F f, Args ... args)
{
	if (!_stats)
		return std::numeric_limits<double>::max();

	double error = 0;
	double range = (_stats->max() - _stats->min()) / _stats->histogramNumClasses();
	double i = _stats->min();
	double j = i + range;
	double total_elements = _stats->numElements();

	Traits<Integrator_if>::Implementation integrator;

	for (auto c = 0; c < _stats->histogramNumClasses(); c++, i = j, j += range)
	{
		double area = integrator.integrate(i, j, f, args...);
		error += pow(_stats->histogramClassFrequency(c) / total_elements - area, 2);
	}

	return error;
}

#endif /* FITTERDEFAULT_h */

