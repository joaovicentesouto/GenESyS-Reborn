/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ProcessAnalyser_if.h
 * Author: cancian
 *
 * Created on 10 de Outubro de 2018, 14:26
 */

#ifndef PROCESSANALYSER_IF_H
#define PROCESSANALYSER_IF_H

#include "List.h"
#include "SimulationScenario_if.h"
#include "SimulationControl.h"
#include "SimulationResponse.h"
#include "Listener.h"

class ProcessAnalyser_if {
public:
	virtual List<SimulationScenario_if*>* getScenarios() const = 0;
	virtual List<SimulationControl_if*>* getControls() const = 0;
	virtual List<SimulationResponse*>* getResponses() const = 0;
	virtual List<SimulationControl_if*>* extractControlsFromModel(std::string modelFilename) const = 0;
	virtual List<SimulationResponse*>* extractResponsesFromModel(std::string modelFilename) const = 0;
	virtual bool attachTogether(SimulationScenario_if* scenario, SimulationControl_if* control, SimulationResponse* response);
	virtual void startSimulation() = 0;
	virtual void stopSimulation() = 0;
	virtual void addTraceSimulationListener(traceSimulationProcessListener traceSimulationProcessListener) = 0;
};

#endif /* PROCESSANALYSER_IF_H */

