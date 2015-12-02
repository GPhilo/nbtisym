#pragma once
#include "utils.hpp"

using namespace std;

class ProcessorState
{
public:
	string id;
	double dutycycle;
	double freq;
	double vdd;
	int index=-1;

	ProcessorState() 
	{
	}

	ProcessorState(string id, double vdd, double freq, double dutycycle) :
		id(id), vdd(vdd), freq(freq), dutycycle(dutycycle)
	{
	}
};