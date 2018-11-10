/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Queue.h
 * Author: cancian
 *
 * Created on 21 de Agosto de 2018, 17:12
 */

#ifndef QUEUE_H
#define QUEUE_H

#include "ModelInfrastructure.h"
#include "LinkedBy.h"

#include "List.h"
#include "Entity.h"
#include "Waiting.h"

#include "StatisticsCollector.h"

class Queue : public ModelInfrastructure, public LinkedBy {
public:
	Queue(Model* model);
	Queue(const Queue& orig);
	virtual ~Queue();
public:
	virtual std::string show();
public:
	void insertElement(Waiting* element);
	void removeElement(Waiting* element, double tnow);
	unsigned int size();
	Waiting* first();
	//List<Waiting*>* getList() const; // can't give direct access so Queue can collect statistics
public: //g&s
	
	
protected: 
	virtual void _loadInstance(std::list<std::string> words);
	virtual std::list<std::string>* _saveInstance();
	virtual bool _verifySymbols(std::string* errorMessage);
	

private: //1::n
	List<Waiting*>* _list = new List<Waiting*>();
private: //1::1
	StatisticsCollector* _cstatNumberInQueue; // = new StatisticsCollector("Number In Queue");
	StatisticsCollector* _cstatTimeInQueue; // = new StatisticsCollector("Time In Queue ");
};

#endif /* QUEUE_H */

