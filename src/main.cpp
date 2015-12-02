#include <stdio.h>
#include <iostream>
#include "ezOptionParser.hpp"
#include "ProcessorState.hpp"
#include "parser\driver.hh"
#include "easylogging++.h"

INITIALIZE_EASYLOGGINGPP

using namespace ez;
using namespace std;

void Usage(ezOptionParser& opt) {
	std::string usage;
	opt.getUsage(usage);
	std::cout << usage;
};

struct temp_record {
	long timestamp = 0;
	long duration = 0;
	double temp = 0;
	double prob = 0;
};

/*****************CONSTANTS DEFINITIONS****************/
//Physical constants
const double q = 1.602e-19; //C - charge of the electron
const double k = 8.6173324e-5; //ev/°K - Boltzmann constant
/******************************************************/

void runStaticSym(
	string output_filename,
	double Vdd,
	double alpha,
	double Tclk,
	double Vth,
	bool do_delay,
	string path_descr_filename,
	string delay_filename,
	int start_time,
	int end_time,
	int timestep,
	double tox,
	double eta_ox,
	double Ea,
	double n,
	double Tref,
	double T0,
	double K,
	double csi1,
	double csi2,
	double E0
	);

void runDynamicSym(
	string state_def_filename,
	string trace_filename,
	string output_filename,
	bool do_delay,
	string path_descr_filename,
	string delay_filename,
	int start_time,
	int end_time,
	int timestep,
	double Vth,
	double tox,
	double eta_ox,
	double Ea,
	double n,
	double Tref,
	double T0,
	double K,
	double csi1,
	double csi2,
	double E0
	);

int main(int argc, const char * argv[]) {
	
	ezOptionParser opt;
	opt.overview = "Simulator for NBTI-induced aging in single-core microprocessors.";
	opt.syntax = "nbtisym MODEL_TYPE [OPTIONS] ARGUMENTS\n\nMODEL_TYPE = static|dynamic\n\nARGUMENTS for dynamic use: <simulation_stop_time> <state_description_filename> <trace_filename>\nARGUMENTS for static use:  <simulation_stop_time>"; 
	opt.example = "nbtisym dynamic -vth -0.18 -tox 1.4 state_def.txt trace.txt 2000\nnbtisym static -vdd 1.2 -vth -0.18 -alpha 0.5 -freq 1e9 -T 350 -tox 1.4 2000\n\n";
	opt.footer = "NBTISym 1.0.0 - Filippo Garolla\nThis program is free and without warranty.\n";

#pragma region command line parameters definitions

	// Show the usage of the tool
	opt.add(
		"", // Default.
		0, // Required?
		0, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Display usage instructions.", // Help description.
		"-h",     // Flag token. 
		"-help",  // Flag token.
		"--help", // Flag token.
		"--usage" // Flag token.
		);

	// Set the output file name
	opt.add(
		"output.txt", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Custom output file name", // Help description.
		"--outname"     // Flag token. 
		);

	// Activate delays calculation and set the pcp path definition file name
	opt.add(
		"", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Probable Critical Path definition file name", // Help description.
		"--pcp-desc"     // Flag token. 
		);

	// Set the path delay output file name
	opt.add(
		"delay.txt", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Custom critical path delay output file name", // Help description.
		"--delayname"     // Flag token. 
		);

	//Set the simulation starting time
	opt.add(
		"0", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Simulation start time", // Help description.
		"--start"     // Flag token. 
		);

	//Set the simulation time step
	opt.add(
		"100", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Simulation time step (seconds)", // Help description.
		"--step"     // Flag token. 
		);

	//Set the Vth0
	opt.add(
		"", // Default.
		1, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Threshold voltage at t=0 (volts)", // Help description.
		"-vth"     // Flag token. 
		);

	//Set the tox
	opt.add(
		"", // Default.
		1, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Oxide layer thickness in nanometers", // Help description.
		"-tox"     // Flag token. 
		);

	//Set the eta ox
	opt.add(
		"3.45313e-11", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Permittivity of the SiO2 oxide", // Help description.
		"-eox"     // Flag token. 
		);

	//Set the activation energy for the reaction
	opt.add(
		"0.49", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Activation energy of the reaction (elettronVolts)", // Help description.
		"-ea"     // Flag token. 
		);

	//Set n
	opt.add(
		"0.166666666666666666", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Exponent of the power law", // Help description.
		"-n"     // Flag token. 
		);

	//Set Tref
	opt.add(
		"350", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Temperature to which all traced temperatures are normalized", // Help description.
		"-tref"     // Flag token. 
		);

	//Set T0
	opt.add(
		"1e-8", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Fitting parameter", // Help description.
		"-T0"     // Flag token. 
		);

	//Set K
	opt.add(
		"8e4", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Fitting parameter", // Help description.
		"-K"     // Flag token. 
		);

	//Set csi1
	opt.add(
		"0.9", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Fitting parameter", // Help description.
		"-csi1"     // Flag token. 
		);

	//Set csi2
	opt.add(
		"0.5", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Fitting parameter", // Help description.
		"-csi2"     // Flag token. 
		);

	//Set E0
	opt.add(
		"0.335", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Fitting parameter (volts/nanometer)", // Help description.
		"-E0"     // Flag token. 
		);

	//Set Vdd
	opt.add(
		"1.2", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Vdd of the processor (volt)", // Help description.
		"-vdd"     // Flag token. 
		);

	//Set alpha
	opt.add(
		"0.5", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Duty cycle of the processor", // Help description.
		"-alpha"     // Flag token. 
		);

	//Set core frequency
	opt.add(
		"1e9", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Frequency of the processor", // Help description.
		"-freq"     // Flag token. 
		);

	//Set core temperature
	opt.add(
		"350", // Default.
		0, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Temperature of the processor", // Help description.
		"-temp"     // Flag token. 
		);

#pragma endregion

	opt.parse(argc, argv);

#pragma region command line parameters validation and evaluation

	if (opt.isSet("-h")) {
		Usage(opt);
		return 1;
	}

	if (opt.firstArgs.size() != 2) {
		std::cerr << "ERROR: Expected model type but received none.\n\n";
		Usage(opt);
		return 1;
	}

	if ( opt.firstArgs[1]->compare("dynamic") == 0 && opt.lastArgs.size() != 3) {
		std::cerr << "ERROR: Expected 3 arguments for dynamic model simulation.\n\n";
		Usage(opt);
		return 1;
	}

	if (opt.firstArgs[1]->compare("static") == 0 && opt.lastArgs.size() != 1) {
		std::cerr << "ERROR: Expected 1 argument for static model simulation.\n\n";
		Usage(opt);
		return 1;
	}

	//Special checks for static-only options
	if (opt.firstArgs[1]->compare("static") == 0)
	{
		bool fail = false;
		if (!opt.isSet("-temp"))
		{
			std::cerr << "ERROR: Option -temp is required for static simulation.\n\n";
			fail = true;
		}
		if (!opt.isSet("-alpha"))
		{
			std::cerr << "ERROR: Option -alpha is required for static simulation.\n\n";
			fail = true;
		}
		if (!opt.isSet("-freq"))
		{
			std::cerr << "ERROR: Option -freq is required for static simulation.\n\n";
			fail = true;
		}
		if (!opt.isSet("-vdd"))
		{
			std::cerr << "ERROR: Option -vdd is required for static simulation.\n\n";
			fail = true;
		}
		if (fail)
		{
			Usage(opt);
			return 1;
		}
	}

	std::vector<std::string> badOptions;
	int i;
	if (!opt.gotRequired(badOptions)) {
		for (i = 0; i < badOptions.size(); ++i)
			std::cerr << "ERROR: Missing required option " << badOptions[i] << ".\n\n";
		Usage(opt);
		return 1;
	}

	if (!opt.gotExpected(badOptions)) {
		for (i = 0; i < badOptions.size(); ++i)
			std::cerr << "ERROR: Got unexpected number of arguments for option " << badOptions[i] << ".\n\n";

		Usage(opt);
		return 1;
	}

#pragma endregion


	string state_def_filename;
	string trace_filename;
	if (opt.firstArgs[1]->compare("dynamic") == 0)
	{
		state_def_filename = *opt.lastArgs[1];
		trace_filename = *opt.lastArgs[2];
	}
	string output_filename;
	opt.get("--outname")->getString(output_filename);
	string delay_out_filename;
	opt.get("--delayname")->getString(delay_out_filename);
	string path_descr_filename;
	opt.get("--pcp-desc")->getString(path_descr_filename);

	int start_time; //s - time from which to generate the output trace
	opt.get("--start")->getInt(start_time);
	long end_time = atoi((*opt.lastArgs[0]).c_str()); //s - time until which the output trace is generated
	int timestep; //s - time between two trace points in the output
	opt.get("--step")->getInt(timestep);

#pragma region Model parameters

	//Processor parameters
	double Vth;  //-0.18; //V - Starting threshold voltage
	opt.get("-vth")->getDouble(Vth);

	//Mosfet technology parameters
	double tox;// = 1.4e-9; //m - oxide thickness
	opt.get("-tox")->getDouble(tox);
	tox *= 1e-9;
	double eta_ox;// = 3.45313e-11; //F/m - SiO2 permittivity
	opt.get("-eox")->getDouble(eta_ox);
	double Cox = eta_ox / tox; //F/m^2 - Capacity of the oxide per area unit

	//Reaction model parameters
	double Ea; // = 0.49; //eV - Activation energy for H2 diffusion
	opt.get("-ea")->getDouble(Ea);
	double n; // = 1.0 / 6.0; //Exponent of the power law for H2 diffusion
	opt.get("-n")->getDouble(n);
	double Tref; // = 350; //Reference temperature to which all temp intervals are mapped back to (measured in Kelvin)
	opt.get("-tref")->getDouble(Tref);

	//Fitting parameters
	double T0;// = 1e-8;
	opt.get("-T0")->getDouble(T0);
	double K;// = 8e4;
	opt.get("-K")->getDouble(K);
	double csi1;// = 0.9;
	opt.get("-csi1")->getDouble(csi1);
	double csi2;// = 0.5;
	opt.get("-csi2")->getDouble(csi2);
	double E0;// = 0.335; //V/nm - Electrical field somewhere...
	opt.get("-E0")->getDouble(E0);

	//Static model specific parameters
	double alpha;
	opt.get("-alpha")->getDouble(alpha);
	double Vdd;
	opt.get("-vdd")->getDouble(Vdd);
	double T;
	opt.get("-temp")->getDouble(T);
	double freq;
	opt.get("-freq")->getDouble(freq);

#pragma endregion

	if (opt.firstArgs[1]->compare("dynamic") == 0)
		runDynamicSym(state_def_filename, trace_filename, output_filename, opt.isSet("--pcp-desc"), path_descr_filename, delay_out_filename,
		start_time, end_time, timestep, Vth, tox, eta_ox, Ea, n, Tref, T0, K, csi1, csi2, E0);
	else if (opt.firstArgs[1]->compare("static") == 0)
		runStaticSym(output_filename, Vdd, alpha, 1/freq, Vth, opt.isSet("--pcp-desc"), path_descr_filename, delay_out_filename,
		start_time, end_time, timestep, tox, eta_ox, Ea, n, Tref, T0, K, csi1, csi2, E0);
	else
	{
		LOG(ERROR) << "Unknown model type " << opt.firstArgs[0] << "." << endl;
		Usage(opt);
		return 1;
	}

	cout << "Press ENTER to close the program." << endl << endl;
	fflush(stdin);
	getchar();
	return 0;
}

void runDynamicSym(
	string state_def_filename,
	string trace_filename,
	string output_filename,
	bool do_delay,
	string path_descr_filename,
	string delay_filename,
	int start_time,
	int end_time,
	int timestep,
	double Vth,
	double tox,
	double eta_ox,
	double Ea,
	double n,
	double Tref,
	double T0,
	double K,
	double csi1,
	double csi2,
	double E0
	)
{
	double Cox = eta_ox / tox;

	ifstream _trace(trace_filename);
	if (!_trace.is_open())
	{
		LOG(ERROR) << "Trace file '" + trace_filename + "' not found.";
		return;
	}

#pragma region _1: Read states specifications and init state objects
	map<string, ProcessorState> states;

	// Multi-state reading from the input using a scanner and parser
	calcxx_driver driver;
	driver.trace_parsing = false;
	driver.trace_scanning = false;
	if (!driver.parse(state_def_filename))
	{
		cout << "State definition parsing terminated successfully." << endl;
	}
	else
	{
		LOG(ERROR) << "Errors detected in the states definitions." << endl;
	}

	states = driver.states;

	if (states.size() == 0)
	{
		LOG(ERROR) << "No defined states found." << endl;
		return;
	}
	else
	{
		cout << "Identified the states: " << endl;
		int idx = 0;
		for (auto i = states.begin(); i != states.end(); i++)
		{
			cout << i->first << endl;
			i->second.index = idx++; //Assign an index to all the valid states to ease process on data structures later
		}
		cout << endl;
	}
#pragma endregion

#pragma region _2: Extract states probabilities from the states trace file

	//Data that will be visible outside this section
	double* state_prob = new double[states.size()];

	//Do the rest in a {..} block so to hide locally-used variables (this should be extracted to a function in the future)
	{
		//0) Prepare data structures
		long total_time = 0;
		long* _state_times = new long[states.size()];
		memset(_state_times, 0, sizeof(long)*states.size());

		//1) Read initial state
		long t0;
		float temp; //Ignore the temperature data at this stage
		string id_s0;
		//_stateTrace >> t0 >> id_s0;
		_trace >> t0 >> temp >> id_s0;
		long last = t0;

		//2) Read each new line and update the total time and the corresponding state time.
		long ti;
		string id_si = id_s0;
		string next_s;
		do
		{
			//_stateTrace >> ti >> id_si;
			_trace >> ti >> temp >> next_s;
			if (states.find(id_si) == states.end())
			{
				//Read a state that was not defined. Signal and abort.
				LOG(ERROR) << "Trace file contains record for undefined state " + id_si + "." << endl;
				return;
			}
			total_time += ti - last;
			_state_times[states[id_si].index] += ti - last;
			last = ti;
			id_si = next_s;
			//} while (!_stateTrace.eof());
		} while (!_trace.eof());

		if (total_time <= 0 || total_time != ti - t0) //Minus t0 is needed to account for initial times different from 0
		{
			LOG(ERROR) << "Error in reading the states trace. Total time <= 0 or not adding up to match the last traced record." << endl;
			return;
		}

		//3) Calculate probabilities.
		for (int i = 0; i < states.size(); i++)
			state_prob[i] = _state_times[i] / ((double)total_time);

		cout << "Finished parsing states trace. State probabilities are:" << endl;
		for (auto i = states.begin(); i != states.end(); i++)
		{
			cout << i->first << " : " << state_prob[i->second.index] << endl;
		}
		cout << endl;
	}
#pragma endregion

#pragma region _3: Extract teperatre traces for each state

	//Data that will be visible ourside this section
	vector<temp_record>* _temp_per_state_trace = new vector<temp_record>[states.size()];

	{
		long* _last_times = new long[states.size()];
		long t0;
		long last_t, new_t;
		double last_T, new_T;
		string last_s, new_s;
		temp_record rec;
		_trace.clear();
		_trace.seekg(0, ios::beg); //Bring the file pointer back to the starting position

		//Init
		_trace >> last_t >> last_T >> last_s;
		_trace >> new_t >> new_T >> new_s;
		for (int i = 0; i < states.size(); i++)
			_last_times[i] = -1;

		while (!_trace.fail())
		{
			int s_idx = states[last_s].index;
			if (_last_times[s_idx] < 0)
			{
				rec.timestamp = last_t;
				_last_times[s_idx] = new_t;
			}
			else
			{
				long delta = last_t - _last_times[s_idx];
				rec.timestamp = _last_times[s_idx];
				_last_times[s_idx] = new_t - delta;
			}
			rec.temp = last_T;
			_temp_per_state_trace[s_idx].push_back(rec);
			last_s = new_s;
			last_t = new_t;
			last_T = new_T;
			_trace >> new_t >> new_T >> new_s;
		}

		//At the end of the trace, apped for each sub-trace a record with a fake temperature, just to have the final time for each state
		for (int i = 0; i < states.size(); i++)
		{
			if (_last_times[i] < 0)
				continue;
			rec.temp = 0;
			rec.timestamp = _last_times[i];
			_temp_per_state_trace[i].push_back(rec);
		}
		delete _last_times;

		cout << "Finished building the per-state temperature traces." << endl << endl;
	}
#pragma endregion

#pragma region _4: For each state calculate the thermal distribution

	//Data visible outside of the section
	map<double, temp_record>* temp_distrib = new map <double, temp_record>[states.size()]; //Use a map instead of a list for the temperatures because the map is ordered
	for (int i = 0; i < states.size(); i++)
	{
		if (_temp_per_state_trace[i].size() <= 1) //Safety check to exclude states with 0 probability
			continue;

		//For each state 
		long t0;
		long last_t, new_t;
		double last_T, new_T;
		auto it = _temp_per_state_trace[i].begin();

		//Init
		last_t = it->timestamp;
		last_T = it->temp;
		it++;

		t0 = last_t;

		while (it != _temp_per_state_trace[i].end())
		{
			new_t = it->timestamp;
			new_T = it->temp;
			temp_distrib[i][last_T].duration += new_t - last_t;
			temp_distrib[i][last_T].temp = last_T; //Redundant for all updates but the first one
			last_t = new_t;
			last_T = new_T;
			it++; //TODO: Check that this works for traces with just 2 temperature records
		}

		for (auto it = temp_distrib[i].begin(); it != temp_distrib[i].end(); it++)
			it->second.prob = it->second.duration / (double)(last_t - t0);
	}
	cout << "Finished parsing temperature trace. Temperature probabilities are:" << endl;
	for (auto it_s = states.begin(); it_s != states.end(); it_s++)
	{
		cout << "State " << it_s->first << ":" << endl;
		for (auto it = temp_distrib[it_s->second.index].begin(); it != temp_distrib[it_s->second.index].end(); it++)
			cout << it->second.temp << " : " << it->second.prob << endl;
		cout << endl;
	}
	cout << endl;
#pragma endregion

#pragma region _5: Calcola gli Ri per ogni stato

	//Scorri la lista di prob temp e per ogni item calcola la formula e somma 
	map<string, double> Ri;
	for (auto i = states.begin(); i != states.end(); i++)
	{
		//Get the temp distribution of this state

		//Do the integral to calculate Ri (go through the distribution and add together the results for each temperature)
		Ri[i->second.id] = 0.0f;
		for (auto tr = temp_distrib[i->second.index].begin(); tr != temp_distrib[i->second.index].end(); tr++)
		{
			Ri[i->second.id] += std::exp((Ea / k)*(1 / Tref - 1 / tr->second.temp))*state_prob[i->second.index] * tr->second.prob;
		}
	}

	cout << "Finished Calculating Ri:" << endl;
	for (auto it = Ri.begin(); it != Ri.end(); it++)
	{
		cout << "Ri[" + it->first + "] : " << it->second << endl;
	}
	cout << endl;

#pragma endregion

#pragma region _6: Calcola C a Tref e Kv_i a Tref per ogni stato.
	double C = (1 / T0)*exp(-Ea / (k*Tref));
	map<string, double> Kv;
	for (auto i = states.begin(); i != states.end(); i++)
	{
		double Eox = (-i->second.vdd - Vth) / (tox*1e9);
		double m1 = pow(q*tox / eta_ox, 3);
		double m2 = K*K*Cox*(-i->second.vdd - Vth)*sqrt(C);
		double m3 = exp(2 * Eox / E0);
		Kv[i->first] = m1*m2*m3*1e27; //Multiplication by 10^27 is to account for the nm->m conversion
	}

	cout << "C at Tref=" + to_string(Tref) + " K : " << C << endl;
	for (auto i = Kv.begin(); i != Kv.end(); i++)
		cout << "Kv[" + i->first + "] : " << i->second << endl;
	cout << endl;
#pragma endregion

#pragma region _7: Calcola i phi_i per ogni stato

	map<string, double> phi;
	for (auto i = states.begin(); i != states.end(); i++)
	{
		double te = tox;
		double Tclk_hat = Ri[i->first] * (1.0 / i->second.freq);
		double alpha_hat = i->second.dutycycle; //Alpha hat = alpha because Duty Cycle is independent from the temperature in a state.
		phi[i->first] = pow((2 * n*sqrt(Kv[i->first] * Kv[i->first] * alpha_hat*C*Tclk_hat)) / (2 * csi1*te + sqrt(csi2*C*(1 - alpha_hat)*Tclk_hat)), 2 * n);
	}

	for (auto i = phi.begin(); i != phi.end(); i++)
		cout << "Phi[" + i->first + "] : " << i->second << endl;
	cout << endl;
#pragma endregion

#pragma region _8: Calcola omega
	double omega = 0.0f;
	for (auto i = states.begin(); i != states.end(); i++)
	{
		omega += pow(phi[i->first], 1.0 / n)*Ri[i->first];
	}
	omega = pow(omega, n);
	cout << "Omega : " << omega << endl << endl;
#pragma endregion

#pragma region _9: Calcola traccia di output delta_Vth

	ofstream out(output_filename);
	for (long t = start_time; t <= end_time; t += timestep)
	{
		out << t << ' ' << omega*pow((double)t, n) * 1000 << endl;
	}
	out.close();
	cout << "Output written to " << output_filename << endl;

#pragma endregion 

#pragma region _10: Calcola delays

	// Do this section only if explicitly requested by the user
	if (do_delay)
	{
		// Load the path description parameters
		list<double> vel_sat_indices;
		list<double> fresh_delays;

		ifstream _path(path_descr_filename);
		if (!_path.is_open())
		{
			LOG(ERROR) << "Requested delay elaboration but path decription file '" + path_descr_filename + "' was not found!";
			return;
		}

		while (!_path.fail())
		{
			double vsi, fd;
			_path >> vsi >> fd;
			vel_sat_indices.push_back(vsi);
			fresh_delays.push_back(fd);
		}

		// Calculate the delays for each timestep
		ofstream out_d(delay_filename);
		for (long t = start_time; t <= end_time; t += timestep)
		{
			double delta_D = 0;
			double Vgs = 0; //We estimate it as a weighted average of the Vgs of all the processor states, weighted by the state probability
			for (auto s = states.begin(); s != states.end(); s++)
				Vgs += -s->second.vdd*state_prob[s->second.index]; //Vgs = -Vdd
			for (auto beta = vel_sat_indices.begin(), fd = fresh_delays.begin(); beta != vel_sat_indices.end() && fd != fresh_delays.end(); fd++, beta++)
			{
				delta_D += ((*beta * omega*pow((double)t, n)) / abs(Vgs - Vth)) * *fd;
			}
			out_d << t << ' ' << delta_D << endl;
		}
		out_d.close();
		cout << "Output written to " << delay_filename << endl;
	}

#pragma endregion
}

void runStaticSym(
	string output_filename,
	double Vdd,
	double alpha,
	double Tclk,
	double Vth,
	bool do_delay,
	string path_descr_filename,
	string delay_filename,
	int start_time,
	int end_time,
	int timestep,
	double tox,
	double eta_ox,
	double Ea,
	double n,
	double Tref,
	double T0,
	double K,
	double csi1,
	double csi2,
	double E0
	)
{
	double Cox = eta_ox / tox;

	//1) Compute C e Kv
	double C = (1 / T0)*exp(-Ea / (k*Tref));
	double Kv;
	double Eox = (-Vdd - Vth) / (tox*1e9);
	double m1 = pow(q*tox / eta_ox, 3);
	double m2 = K*K*Cox*(-Vdd - Vth)*sqrt(C);
	double m3 = exp(2 * Eox / E0);
	Kv = m1*m2*m3*1e27; //Multiplication by 10^27 is to account for the nm->m conversion

	cout << "C = " << C << endl;
	cout << "Kv = " << Kv << endl;
	cout << endl;

	//2) Iterate over t and write output

	ofstream out(output_filename);
	for (long t = start_time; t <= end_time; t += timestep)
	{
		//Compute beta_t
		double beta_t = 1 - ((2*csi1*tox + sqrt(csi2*C*(1-alpha)*Tclk)) / (2*tox+sqrt(C*t)));

		//compute delta_vth;
		double delta_vth = pow(sqrt(Kv*Kv*alpha*Tclk) / (1 - pow(beta_t, 1 / (2 * n))), 2 * n);

		out << t << ' ' << delta_vth * 1000 << endl;
	}
	out.close();
	cout << "Output written to " << output_filename << endl;
}