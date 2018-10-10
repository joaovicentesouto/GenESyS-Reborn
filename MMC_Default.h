/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MMC_Default.h
 * Author: Bruno Bonotto, Fabíola Kretzer e João Vicente Souto
 */

#ifndef MMC_IF_H
#define MMC_IF_H

#include <cstdlib>
#include <cmath>
#include <stdarg.h>

#include "Sampler_if.h"

class MMC_Default : public Sampler_if {
public:

	struct MyRNG_Parameters : public RNG_Parameters {
		unsigned int seed;
		unsigned int module;
		unsigned int multiplier;
	};
	
public:
	MMC_Default() = default;
	MMC_Default(const MMC_Default& orig) = delete;
	~MMC_Default() = default;
public: // probability distributions
	double random();
	double sampleUniform(double min, double max);
	double sampleExponential(double mean);
	double sampleErlang(double mean, int M);
	double sampleNormal(double mean, double stddev);
	double sampleGamma(double mean, double alpha);
	double sampleBeta(double alpha, double beta, double infLimit, double supLimit);
	double sampleWeibull(double alpha, double scale);
	double sampleLogNormal(double mean, double stddev);
	double sampleTriangular(double min, double mode, double max);
	double sampleDiscrete(double value, double acumProb, ...);

  public:
	void setRNGparameters(RNG_Parameters* param);
	RNG_Parameters* getRNGparameters() const;

private:
	double sampleGammaJonk(double alpha);

	unsigned int seed{12345};
	unsigned int module{3};
	unsigned int multiplier{5};
	double _normal_value{0};
	bool _normal{true};

};

#endif /* MMC_Default_H */
