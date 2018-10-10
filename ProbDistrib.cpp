/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ProbDistrib.cpp
 * Author: João Vicente Souto
 * 
 * Created on 23 de Agosto de 2018, 17:25
 */

#include "ProbDistrib.h"

double ProbDistrib::uniform(double x, double min, double max) {
    return 1 / (max - min);
}

double ProbDistrib::exponential(double x, double mean) {
    return (1 / mean) * exp(-(1 / mean) * x);
}

int factorial(int n) {
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double ProbDistrib::erlang(double x, double mean, int M) {
    return pow(x, M - 1) * exp(- x / mean) / (pow(mean, M) * factorial(M - 1));
}

double ProbDistrib::normal(double x, double mean, double stddev) {
    return (1 / (stddev * sqrt(2 * M_PI))) * exp(-0.5 * pow((x - mean), 2) / pow(stddev, 2));
}

double ProbDistrib::gamma(double x, double mean, double alpha) {
}

double gamma_function(double n) {
	return sqrt(2 * M_PI * n) * pow((n/exp(1)), n);
}

double beta_auxiliar(double x, double y) {
	return gamma_function(x) * gamma_function(y) / gamma_function(x + y);
}

double ProbDistrib::beta(double x, double alpha, double beta) {
    return (1 / beta_auxiliar(alpha, beta)) * pow(x, alpha - 1) * pow(1 - x, beta - 1);
}

double ProbDistrib::weibull(double x, double alpha, double scale) {
    if (x < 0)
		return 0;
	
	return alpha / scale * pow(x / scale, alpha - 1) * exp(-pow(x/scale, alpha));
}

double ProbDistrib::logNormal(double x, double mean, double stddev) {
}

double ProbDistrib::triangular(double x, double min, double mode, double max) {
    if (min <= x && x < mode)
		return 2 * (x - min) / ((max - min) * (mode - min));

	if (x == mode)
		return 2 / (max - min);
	
	if (mode < x && x <= max)
		return 2 * (max - x) / ((max - min) * (max - mode));
	
	return 0;
}

