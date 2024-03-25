/*****************************************************************************************
	
	Written by Nate Lynd and Louise Kuehster (C) 2024
	
	This program will fit a data to a numerically defined model... hopefully
	
	g++ -I /Users/nl7597/Documents/Calculations/Eigen/ -O3 -msse3 -o fit_stochastic fit_stochastic_reversible.cpp

	./fit_stochastic input.in output.out 0.043524 3.386726 5.080089 14.0 3.11 3.897 11.241 0.0 0.54 0.54 0.0 0.0000001 0.01
	
	./fit_stochastic inputPLGA.in outputPLGA.out 0.043524 3.386726 5.080089 16.4504 6.14505 2.14938 8.78906 0.917694 0.917744 0.853395 0.703005 0.00268994 0.0038

	
	Next up:
	
	Dynamic sizing of stochastic model
	
	Fully reversible copolymerization of one cyclic dimer and one monomer- 
		forward parameters PA+AA (kAA)	PA+B (kAB)	PB+B (kBB)	PB+AA (kBA)
		bakward parameters PA-AA (kA-A)	PA-B (kA-B)	PB-B (kB-B)	PB-AA (kB-A)

*****************************************************************************************/
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include "Eigen"
#include "polymer.h"

using namespace std;
using namespace Eigen;

// the objective functions are sum of squared residuals
double		residual(VectorXd x_data,VectorXd y_data,VectorXd x_model,VectorXd y_model);
double		residual_time(VectorXd x_data,VectorXd y1_data,VectorXd y2_data,
					VectorXd x_model,VectorXd y1_model,VectorXd y2_model,double MonomerA, double MonomerB);
double		residual_compositional(VectorXd y1_data,VectorXd y2_data,VectorXd conv_data,
				VectorXd conv_model,VectorXd y1_model,VectorXd y2_model);

VectorXd	conversion(VectorXd concentrations,VectorXd y1y10_data,VectorXd y2y20_data);

// line search - only used on conjugate gradient
double		linear_optim(double x1,double x2,double fx1,double fx2);
double		brents_optim(double x1,double x2,double x3,double fx1,double fx2,double fx3);

// conjugate gradient optimization
VectorXd	cg_optimization8(VectorXd parameters_old,VectorXd concentrations,int number,
				VectorXd x_data,VectorXd y1_data,VectorXd y2_data,double accuracy_goal);
// line searches using Brent's method may converge to a max or min. This is a test
bool		is_it_a_minimum(VectorXd lambda,VectorXd residuals);

// Downhill simplex - brute force but no derivative required and seems tolerant of noise
void 		get_psum(MatrixXd &p, VectorXd &psum, const int ihi);
double 		amotry_time(MatrixXd &p, VectorXd &y, VectorXd &psum, const int ihi, const double fac,
			VectorXd concentrations,int number,VectorXd x_data,VectorXd y1_data,VectorXd y2_data);
double 		amotry_composition(MatrixXd &p, VectorXd &y, VectorXd &psum, const int ihi, const double fac,
			VectorXd concentrations,int number,VectorXd x_data,VectorXd y1_data,VectorXd y2_data);
void		swap(double &a,double &b);

VectorXd	simplex_time(VectorXd parameters_old,VectorXd concentrations,int number,
				VectorXd x_data,VectorXd y1_data,VectorXd y2_data,double accuracy_goal);
VectorXd	simplex_composition(VectorXd parameters_old,VectorXd concentrations,int number,
				VectorXd x_data,VectorXd conv_data_A,VectorXd conv_data_B,double accuracy_goal);

// the model
MatrixXd	stochastic2(VectorXd concentrations,int i_number,VectorXd parameters,double t_final);
MatrixXd	stochastic8(VectorXd concentrations,int i_number,VectorXd parameters,double t_final);

VectorXd	linspace(double low,double high,int points);

// interpolate that connects model to data
double		y_at_x(double at_x,VectorXd x_model,VectorXd y_model);

// filehandling functions
double		doubleFromStringBuffer(string buffer);

int			integerFromStringBuffer(string buffer);

string		stringAfterEqualsSign(string buffer);

string		stringBetweenBrackets(string buffer);

double		doubleFromString(string substring);

bool		existsInString(string test,string buffer);

// this is just to add a function for VectorXd that is like vector< class T>'s method
// vector<class T>.push_back(class T)
namespace Eigen {
	namespace Vector {
	template <class T>
		inline void push_back(Eigen::Matrix<T,Eigen::Dynamic,1>& v, const T d) 
		{
			Eigen::Matrix<T,Eigen::Dynamic,1> tmp = v;
			v.resize(tmp.size() + 1);
			v.head(tmp.size()) = tmp;
			v[v.size()-1] = d;
		}
	} // namespace Vector
} // namespace Eigen


int main(int argc, char* argv[])
{
	string	input_file;
	string	output_file;
	
	ifstream	input;
	ofstream	output;
	string		buffer;
	
	//think abou this
	double		Initiator;	 // command line input for initiator concentration
	double		MonomerA,MonomerB;	 // command line input for monomer concentration
	
	int			n_coordinates = 1; // time
	int			n_observables = 3; // time and two conversions pG and pL
	int			n_parameters = 9;
	
	// make these commandline parameters
	bool		compositional = false;
	bool		temporal = true;
	
	VectorXd	concentrations = VectorXd::Zero(n_observables);
	VectorXd	parameters = VectorXd::Zero(n_parameters);
	VectorXd	ratios = VectorXd::Zero(n_parameters-3);
	
	// for the compositional fit - fix these first, then go after time-parameters
	double		rA,rB,rA_A,rA_B,rB_A,rB_B; // 6
	// for the time-dependent fit
	double		kAA,kAB,kBB,kBA,kAmA,kAmB,kBmA,kBmB,kT; // 9

	double		r;	 // final SSR
	
	double		termination_time;
	double		accuracy_goal;
	
	// Keep this low when far from a solution to keep rate high, then increase when 
	// solution nears to increase accuracy and decrease noise:
	int			number = 50; // higher for lower error/noise
	// stride will decrease the density of the fit output
	int			stride = 1;
	int			datapoints = 8; // default
	
	VectorXd	conv_data_all;
	VectorXd	conv_model_all;
	VectorXd	time_model,conv_model_A,conv_model_B;	// calculated stochastically
	VectorXd	time_data,conv_data_A,conv_data_B;	// read from file
	MatrixXd	model_data_package;		// stochastic will fill this
	
	double		t_final;
	
	// read default values from command line
	
	input_file			=	string(argv[1]);
	output_file			=	string(argv[2]);
	Initiator			=	atof(argv[3]);
	MonomerA			=	atof(argv[4]);
	MonomerB			=	atof(argv[5]);
	kAA					=	atof(argv[6]);
	kAB					=	atof(argv[7]);
	kBB					=	atof(argv[8]);
	kBA					=	atof(argv[9]);
	kAmA				=	atof(argv[10]);
	kAmB				=	atof(argv[11]);
	kBmB				=	atof(argv[12]);
	kBmA				=	atof(argv[13]);
	kT					=	atof(argv[14]);
	accuracy_goal		=	atof(argv[15]);
	
	// package parameters to use functions
	concentrations(0)	= Initiator;
	concentrations(1)	= MonomerA;
	concentrations(2)	= MonomerB;
	parameters(0)		= kAA;
	parameters(1)		= kAB;
	parameters(2)		= kBB;
	parameters(3)		= kBA;
	parameters(4)		= kAmA;
	parameters(5)		= kAmB;
	parameters(6)		= kBmB;
	parameters(7)		= kBmA;
	parameters(8)		= kT;
	// 
	// Stochastic fit to data for a reversible equilibrium polymerization.
	//
	
	// Open file to read data - reads line by line
	input.open(input_file,ios::in);
	if (input)
	{
		do
		{
			getline(input,buffer);
			// extract information from input buffer
			if (existsInString("number_of_points",buffer))
			{
				datapoints = integerFromStringBuffer(buffer);
				// allocate vectors
				time_data.resize(datapoints);
				conv_data_A.resize(datapoints);
				conv_data_B.resize(datapoints);
				
			}
			else if (existsInString("time",buffer))
			{
				// step through until end
				int i=0;
				do {
					getline(input,buffer);
					// get double from buffer
					time_data[i] = doubleFromString(buffer);
					// add to vector
					i += 1;
				} while (i<datapoints); // need error in case number_of_points wrong
			}
			else if (existsInString("conversion A",buffer))
			{
				// step through until end
				int i=0;
				do {
					getline(input,buffer);
					// get double from buffer
					conv_data_A[i] = doubleFromString(buffer);
					// add to vector
					i += 1;
				} while (i<datapoints);
			}
			else if (existsInString("conversion B",buffer))
			{
				// step through until end
				int i=0;
				do {
					getline(input,buffer);
					// get double from buffer
					conv_data_B[i] = doubleFromString(buffer);
					// add to vector
					i += 1;
				} while (i<datapoints);
			}
		} while (buffer.compare("end")!=0);
	}
	else
	{
		// no file, report and quit
		cout << "Fit::File Not Found!" << endl;
		exit(0);
	}
	input.close();

	cout << "********** FITTING TO STOCHASTIC MODEL **********" << endl;
	cout << "input file = " << input_file << endl;
	cout << "output file = " << output_file << endl;
	cout << "********************** INITIAL VALUES *********************" << endl;
	cout << "{c} = " << concentrations.transpose() << endl;
	cout << "{k} = " << parameters.transpose() << endl;
	cout << "kAA = " << parameters(0) << endl;
	cout << "kAB = " << parameters(1) << endl;
	cout << "kBB = " << parameters(2) << endl;
	cout << "kBA = " << parameters(3) << endl;
	cout << "kAmA = " << parameters(4) << endl;
	cout << "kAmB = " << parameters(5) << endl;
	cout << "kBmB = " << parameters(6) << endl;
	cout << "kBmA = " << parameters(7) << endl;
	cout << "kT = " << parameters(8) << endl;
	cout << "rA = " << parameters(0)/parameters(1) << endl;
	cout << "rB = " << parameters(2)/parameters(3) << endl;
	cout << "rA-A = kA-A/kAA = " << parameters(4)/parameters(0) << endl;
	cout << "rA-B = kA-B/kAB = " << parameters(5)/parameters(1) << endl;
	cout << "rB-B = kB-B/kBB = " << parameters(6)/parameters(2) << endl;
	cout << "rB-A = kB-A/kBA = " << parameters(7)/parameters(3) << endl;
	cout << "********************** FILE DATA *********************" << endl;
	cout << "data points = " << datapoints << endl;
	cout << "time array read = " << time_data.size() << " values" << endl;
	cout << "conv A array read = " << conv_data_A.size() << " values" << endl;
	cout << "conv B array read = " << conv_data_B.size() << " values" << endl;
	
	t_final = time_data[time_data.size()-1];
	
	cout << "t_final = " << t_final << endl;
	
	conv_data_all = conversion(concentrations,conv_data_A,conv_data_B);
	conv_model_all = conversion(concentrations,conv_model_A,conv_model_B);

	// begin main fitting routine
	if (compositional==true) // change this to compositional
	{
		cout << "optimizing based on compositional data..." << endl;
		
		// need to fix some of the parameters... this really needs a specialized
		// simplex code for the fewer parameters:
		rA = kAA/kAB; // fix kAA and kBB
		rB = kBB/kBA;
		rA_A = kAmA/kAA; // kAA fixed
		rA_B = kAmB/kAB; // = kAmB/(kAA/rA);
		rB_B = kBmB/kBB; // kBB fixed
		rB_A = kBmA/kBA; // = kBmA/(kBB/rB);
		
		parameters = simplex_composition(parameters,concentrations,number,time_data,conv_data_A,conv_data_B,accuracy_goal);
		
	}
	if (temporal==true)
	{
		cout << "refining solution based on temporal data..." << endl;
		// increase number of chains to improve accuracy
		number = 100;
		//parameters = simplex9(parameters,concentrations,number,time_data,conv_data_A,conv_data_B,accuracy_goal);
		parameters = simplex_time(parameters,concentrations,number,time_data,conv_data_A,conv_data_B,accuracy_goal);
	}
	// use final parameters to create fit curve
	// calculate new r with new parameters and larger number of chains
	cout << "generating output file data..." << endl;
	number = 100;
	model_data_package = stochastic8(concentrations,number,parameters,t_final);
	time_model 		=  model_data_package.col(0);
	conv_model_A	=  model_data_package.col(1);
	conv_model_B	=  model_data_package.col(2);
	conv_model_all = conversion(concentrations,conv_model_A,conv_model_B);
	r =	residual_time( time_data, conv_data_A, conv_data_B, time_model, conv_model_A,conv_model_B,MonomerA, MonomerB);
	// output
	cout << "********************** RESULTS *********************" << endl;
	cout << "Sum of squares = " << r << endl;
	
	output.open(output_file.c_str(),ios::out);
	output << "********************** RESULTS *********************" << endl;
	output << "Sum of squared residuals (SSR) = " << r << endl;
	output << "{c} = " << concentrations.transpose() << endl;
	output << "{k} = " << parameters.transpose() << endl;
	output << "kAA = " << parameters(0) << endl;
	output << "kAB = " << parameters(1) << endl;
	output << "kBB = " << parameters(2) << endl;
	output << "kBA = " << parameters(3) << endl;
	output << "kAmA = " << parameters(4) << endl;
	output << "kAmB = " << parameters(5) << endl;
	output << "kBmB = " << parameters(6) << endl;
	output << "kBmA = " << parameters(7) << endl;
	output << "kT = " << parameters(8) << endl;
	output << "rA = " << parameters(0)/parameters(1) << endl;
	output << "rB = " << parameters(2)/parameters(3) << endl;
	output << "rA-A = kA-A/kAA = " << parameters(4)/parameters(0) << endl;
	output << "rA-B = kA-B/kAB = " << parameters(5)/parameters(1) << endl;
	output << "rB-B = kB-B/kBB = " << parameters(6)/parameters(2) << endl;
	output << "rB-A = kB-A/kBA = " << parameters(7)/parameters(3) << endl;

	output << "********************* FIT DATA *********************" << endl;
	output << "data = {";
	for (int i=0;i<=time_data.size()-2;i++)
	{
		output << "{" << time_data(i) << ",\t" << conv_data_A(i) << ",\t" << conv_data_B(i) << "}," << endl;
	}
	output << "{" << time_data.tail(1) << ",\t" << conv_data_A.tail(1) << ",\t" << conv_data_B.tail(1) << "}};" << endl;
	output << "********************* BEST FIT *********************" << endl;
	output << "fit = {";
	for (int i=0;i<=time_model.size()-2;i+=stride)
	{
		//output << time_model(i,0) << "\t" << fGConv(i,1) << endl;
		output << "{" << time_model(i) << ",\t" << conv_model_A(i) << ",\t" << conv_model_B(i) << "}," << endl;
	}
	output << "{" << time_model.tail(1) << ",\t" << conv_model_A.tail(1) << ",\t" << conv_model_B.tail(1) << "}};" << endl;
	output.close();
	
return 0;
}
/****************************************************************************************/
MatrixXd	stochastic8(VectorXd concentrations,int i_number,VectorXd parameters,double t_final)
{
	VectorXd	time = VectorXd::Zero(1);
	VectorXd	conversionA = VectorXd::Zero(1);
	VectorXd	conversionB = VectorXd::Zero(1);
	MatrixXd	package;
	
	const double	avogadros = 6.022149e+23; // no. per mole.
	
	random_device rd;
	double      rand = rd();
	mt19937		rando(rd());
	double		rando_maximo = 4294967295;
	double		a,b,c;
	int			selected_chain;

	int			last_last_reaction = 0;
	int			last_reaction = 0;
	bool		ester_found = false;
	bool		debug = false;
	int			selected_ester=0;
	int			ester_count=0;
	int			chain1=0; // not used?
	int			chain2=0;
	int			chain3=0;
	int			monomer1=0;
	int			monomer2=0;
	vector<int> chain1sequence;
	vector<int> chain2sequence;
	int			chain1length;
	int			chain2length;
	
	double		I0,A0,B0;
	int			nP,nPA,nPB,nPAA,nPAB,nPBA,nPBB,nPAAA,nPABB,nPABA,nPAAB,nPBAA,nPBBB,nPBBA,nPBAB,nE=0;
	int			nAA, nAB, nBA, nBB, nIM;//nP is total number of polymers. nIM is total number of linear chains that have initiated
	double		kA,kB;
	double		kAA,kAB,kBA,kBB,kAmA,kAmB,kBmA,kBmB,kT,kTGG, kTGL, kTLG, kTLL, kTI; //regular rate constants (units depend on rate equation)
	double		kA_d, kB_d, kAA_d, kAB_d, kBA_d, kBB_d, kAmA_d, kAmB_d, kBmA_d, kBmB_d;
	double		kT_d,kTAA_d, kTAB_d, kTBA_d, kTBB_d, kTI_d; //discretized rate constants (units are all in seconds)
	double		termination_time;
	double		cat_ratio = 1.; //molar ratio of catalyst to initiator ([cat]/[in])
	double		volume;
	double		t = 0.;
	
	int			ester_sum,ester_id;
	int			counter = 0;
	int			chains = i_number; //This is the number of initiators / number of linear chains
	int			I,A,B;
	int			nI0,nA0,nB0;
	int			Atype = 1, Btype = 2;
	
	int			repeats_removed,repeats_added;
	
	vector<int>			degreeOfPolymerization;
	vector<int>			initiated(0);
	//vector<int>			un_initiated(chains); // speed the locating of rare uninitiated chains
	
	int					reactions = 11;
	VectorXd			rho = VectorXd::Zero(reactions);
	VectorXd			lower = VectorXd::Zero(reactions);
	VectorXd			upper = VectorXd::Zero(reactions);
	double				r_tot;
	int					n_tot;
	
	//vector<double>	dyad_fractions(4);
	//vector<double>	triad_fractions(8);
	//vector<double>	tetrad_fractions(16);
	//vector<double>  pentad_fractions(32);
	//cout <<fixed << "rand = " << rand << endl;
	// these don't appear to be used...
	nAA = 0; nAB = 0; nBA = 0; nBB = 0; nIM = 0;

	// read input from command line:
	// kinetics file for structured data output:
	// these all need to be input arguments
	//kinetics_file				=	string(argv[1]);
	// initiator molecules = number of chains
	I0							=	concentrations(0);
	// A monomers
	A0							=	concentrations(1);
	// B monomer
	B0							=	concentrations(2);
	// forward rate constants:
	kAA							=	parameters(0);
	kAB							=	parameters(1);
	kBB							=	parameters(2);
	kBA							=	parameters(3);
	// backward rate constantsparameters
	kAmA						=	parameters(4);
	kAmB						=	parameters(5);
	kBmB						=	parameters(6);
	kBmA						=	parameters(7);
	//transesterification rate constants
	kT							=	parameters(8);
	//kTGG						=	atof(argv[13]);
	//kTGL						=	atof(argv[14]);
	//kTLG						=	atof(argv[15]);
	//kTLL						=	atof(argv[16]);
	//kTI							=	atof(argv[17]);

	termination_time			=	t_final;
	chains						=	i_number; // relate external to internal
	cat_ratio					=	1.0; // I would model this by change of rate constants?

	vector<polymer>		polymers(chains);//linear polymer chains - no cycles because I don't think it matters for this

	// for simplicity, we set kA = kAA and kB = kBB
	kA = kAA;
	kB = kBB;
	
	// The system volume will be based on the initiator concentration (number of chains)
	volume = double(chains)/avogadros/I0;
	// set up initial conditions
	nI0 = I = chains; // I(0)
	nA0 = A = int(A0*volume*avogadros); // A(0)
	nB0 = B = int(B0*volume*avogadros); // B(0)

	// Calculate discretized rate constants
	// Forward
	kA_d	=		(kAA/volume)/avogadros;
	kB_d	=		(kBB/volume)/avogadros;
	kAA_d	=		(kAA/volume)/avogadros;
	kAB_d	=		(kAB/volume)/avogadros;
	kBA_d	=		(kBA/volume)/avogadros;
	kBB_d	=		(kBB/volume)/avogadros;
	// Backwards
	kAmA_d	=		kAmA;
	kAmB_d	=		kAmB;
	kBmA_d	=		kBmA;
	kBmB_d	=		kBmB;
	
	// just a single kT to start with. Set these equal.
	kT_d	=		(kT/volume)/avogadros;
	// set all based on a single kT for now:
	kTAA_d = kTAA_d = kTAB_d = kTBA_d = kTAA_d = kTI_d = kT_d;
	
	time(0) = t;
	conversionA(0) = double(nA0-A)/double(nA0);
	conversionB(0) = double(nB0-B)/double(nB0);
	
	//for (int i=0;i<=un_initiated.size()-1;i++)
	//{
	//	un_initiated[i] = i;
	//}
	//cout << "number chains = " << chains << endl;
	
	while (t < termination_time)
	{
		nP=0; nPA = 0; nPB = 0; nPAA = 0; nPAB = 0; nPBA = 0; nPBB = 0; nPAAA = 0;
		nPABB = 0, nPABA = 0, nPAAB = 0, nPBAA = 0, nPBBB = 0, nPBBA = 0, nPBAB = 0; nE=0;
		//This starts over counting end group dyads and triads every time
		for (int i = 0; i <= polymers.size() - 1; i++)
		{
			nE += polymers[i].getDegreeOfPolymerization(); // NOTE or 2*(nA0-A)+(nB0-B)
			if (polymers[i].hasInitiated())
			{ 
				nIM += 1;
				nP += 1;
			
			if (polymers[i].endGroup() == 1)
			{
				nPA += 1;
				if (polymers[i].endGroupDyad() == 11)
				{
					nPAA += 1;
					if (polymers[i].endGroupTriad() == 111)
					{
						nPAAA += 1;
					}
					else if (polymers[i].endGroupTriad() == 211)
					{
						nPBAA += 1;
					}
				}
				else if (polymers[i].endGroupDyad() == 21)
				{
					nPBA += 1;
					if (polymers[i].endGroupTriad() == 121)
					{
						nPABA += 1;
					}
					else if (polymers[i].endGroupTriad() == 221)
					{
						nPBBA += 1;
					}
				}
			}
			else if (polymers[i].endGroup() == 2)
			{
				nPB += 1;
				if (polymers[i].endGroupDyad() == 22)
				{
					nPBB += 1;
					if (polymers[i].endGroupTriad() == 122)
					{
						nPABB += 1;
					}
					else if (polymers[i].endGroupTriad() == 222)
					{
						nPBBB += 1;
					}
				}
				else if (polymers[i].endGroupDyad() == 12)
				{
					nPAB += 1;
					if (polymers[i].endGroupTriad() == 112)
					{
						nPAAB += 1;
					}
					else if (polymers[i].endGroupTriad() == 212)
					{
						nPBAB += 1;
					}
				}
			}
			}
		}

		if (nE != (2*(nA0-A)+(nB0-B))) // check for conservation
		{
			cout << "Error: Nonconservative repeat unit error" << endl;
			cout << "last reaction = " << last_reaction << endl;
			cout << "reaction before that = " << last_last_reaction << endl;
			cout << "I\tA\tB\tnE\tnE'" << endl;
			cout << I << "\t" << A << "\t" << B << "\t" << nE << "\t" << 2*((nA0-A)+(nB0-B)) << endl;
			cout << "u1 = " << b << endl;
			cout << lower.transpose() << endl;
			cout << upper.transpose() << endl;
			cout << rho.transpose() << endl;
			cout << "A0 = " << nA0 << endl;
			cout << "B0 = " << nB0 << endl;
			cout << "nE0 = " << 2.*double(nA0+nB0) << endl;
			cout << "nPA = " << nPA << endl;
			cout << "nPB = " << nPB << endl;
			cout << "nPAA = " << nPAA << endl; 
			cout << "nPAB = " << nPAB << endl;
			cout << "nPBA = " << nPBA << endl;
			cout << "npBB = " << nPBB << endl;
			cout << "nPAAA = " << nPAAA << endl;
			cout << "nPABB = " << nPABB << endl;
			cout << "nPBAA = " << nPBAA << endl;
			cout << "nPBBB = " << nPBBB << endl;
			exit(0);
		}
		nE = 2*(nA0-A)+(nB0-B);
		//2 esters per A cyclic dimer, one ester per B monomer
		//update nAA, nAB, nBA, nBB with each reaction instead of counting constantly 
		// reaction rates- calculated based on numbers of molecules
		// units of molecules/time
		// initiation
		rho(0) = kA_d * I * cat_ratio * A; // 1
		rho(1) = kB_d * I * cat_ratio * B; // 2
		
		// propagation
		rho(2) = kAA_d * nPA * cat_ratio * A; // 3
		rho(3) = kAB_d * nPA * cat_ratio * B; // 4
		rho(4) = kBA_d * nPB * cat_ratio * A; // 5
		rho(5) = kBB_d * nPB * cat_ratio * B; // 6
		
		// depropagation
		rho(6) = kAmA_d * nPAAA * cat_ratio; // 7
		rho(7) = kAmB_d * nPAB * cat_ratio; // 8
		rho(8) = kBmA_d * nPBAA * cat_ratio; // 9
		rho(9) = kBmB_d * nPBB * cat_ratio; // 9
		
		// transesterification
		rho(10) = kT_d * nE * nI0 *cat_ratio; // 10
		
		//rho(10) = kTAA_d * nAA * nI0 * cat_ratio; //number GG esters * number active ends
		//rho(11) = kTAB_d * nAB * nI0 * cat_ratio; //different reaction for each dyad type, but assume independent of chain end
		////rho(12) = kTBA_d * nBA * nI0 * cat_ratio;
		//rho(13) = kTAA_d * nBB * nI0 * cat_ratio;
		//rho(14) = kTI_d * nIM * nI0 * cat_ratio;
		
		r_tot = rho.sum();

		// calculate limits for reaction decision:
		lower(0) = 0.; 			upper(0) = rho(0)/r_tot;
		lower(1) = upper(0); 	upper(1) = lower(1) + rho(1)/r_tot;
		lower(2) = upper(1);	upper(2) = lower(2) + rho(2)/r_tot;
		lower(3) = upper(2);	upper(3) = lower(3) + rho(3)/r_tot;
		lower(4) = upper(3);	upper(4) = lower(4) + rho(4)/r_tot;
		lower(5) = upper(4);	upper(5) = lower(5) + rho(5)/r_tot;
		lower(6) = upper(5);	upper(6) = lower(6) + rho(6)/r_tot;
		lower(7) = upper(6);	upper(7) = lower(7) + rho(7)/r_tot;
		lower(8) = upper(7);	upper(8) = lower(8) + rho(8)/r_tot;
		lower(9) = upper(8);	upper(9) = lower(9) + rho(9)/r_tot;
		lower(10) = upper(9);	upper(10) = lower(10) + rho(10)/r_tot;

		
		if (debug)
		{
			cout << t << "\t"<< I <<"\t" << A << "\t" << B << "\t" << (double(nA0)-double(A))/double(nA0) << "\t" << (double(nB0)-double(B))/double(nB0) << endl;
		}
		// update time
		a = double(rando()) / rando_maximo; // time
		b = double(rando()) / rando_maximo; // reaction
		c = double(rando()) / rando_maximo; // chain proposal
		selected_chain = int(double(chains) * c);
		//cout << "selected_chain initial is " << selected_chain << endl;
		//cout << "b is " << b << endl;
		t = t - log(a) /r_tot;  
		
		// make decision
		if (b >= lower[0] && b <= upper[0])
		{
			last_last_reaction = last_reaction;
			last_reaction = 1;
			//cout << "reaction 1!" << endl;
			//cout << "\tI = " << I << endl;
			// reaction 1 initiation with A monomer
			// which chain (that has DP = 0)
			if (polymers[selected_chain].hasInitiated() == false)
			{
				// do reaction 
				I = I - 1;
				A = A - 1;
				polymers[selected_chain].addMonomerToEnd(Atype);
				polymers[selected_chain].addMonomerToEnd(Atype); //twice because each A monomer is actually a dimer
				initiated.push_back(selected_chain); // this keeps track of initiated chains
				//un_initiated.erase(selected_chain);
				//polymers[selected_chain].printPolymer();
			}
			else // find new chain
			{
				while (polymers[selected_chain].hasInitiated() == true)
				{
					//selected_chain = selectNewChain(selected_chain, chains);
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "\t1 selected_chain = " << selected_chain << endl;					
					if (polymers[selected_chain].hasInitiated() == false)
					{
						// do reaction
						I = I - 1;
						A = A - 1;
						polymers[selected_chain].addMonomerToEnd(Atype);
						polymers[selected_chain].addMonomerToEnd(Atype);
						initiated.push_back(selected_chain);
						//un_initiated.erase(selected_chain);
						break;
					}
				}
			}
		}
		else if (b > lower[1] && b <= upper[1])
		{
			last_last_reaction = last_reaction;
			last_reaction = 2;
			//cout << "reaction 2!" << endl;
			//cout << "\tI = " << I << endl;
			//reaction 2- Initiation with B
			if (polymers[selected_chain].hasInitiated() == false)
			{
				// do reaction
				I = I - 1;
				B = B - 1;
				polymers[selected_chain].addMonomerToEnd(Btype); //Only one because B is just a regular monomer
				initiated.push_back(selected_chain);
				//un_initiated.erase(selected_chain);
			}
			else // find new chain
			{
				while (polymers[selected_chain].hasInitiated() == true)
				{
					//selected_chain = selectNewChain(selected_chain,chains);
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "\t2 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].hasInitiated() == false)
					{
						// do reaction
						I = I - 1;
						B = B - 1;
						//nBB = nBB + 1;
						polymers[selected_chain].addMonomerToEnd(Btype);
						initiated.push_back(selected_chain);
						//un_initiated.erase(selected_chain);
						break;
					}
				}
			}
		}

		// from here only from previously initiated chains
		else if (b > lower[2] && b <= upper[2])
		{
			last_last_reaction = last_reaction;
			last_reaction = 3;
			// reaction 3  kAA*cPA*cA;
			//cout << "reaction 3!" << endl;
			if (polymers[selected_chain].endGroup() == Atype)
			{
				// do reaction
				A = A - 1;
				//nAA = nAA + 2;
				polymers[selected_chain].addMonomerToEnd(Atype);
				polymers[selected_chain].addMonomerToEnd(Atype);
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroup() != Atype)
				{
					//selected_chain = selectNewChain(selected_chain,chains);
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "\t3 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroup() == Atype)
					{
						// do reaction
						A = A - 1;
						//nAA = nAA + 2;
						polymers[selected_chain].addMonomerToEnd(Atype);
						polymers[selected_chain].addMonomerToEnd(Atype);
						//cout << "\t\t reaction executed!" << endl;
						//polymers[selected_chain].printPolymer();
						break;
					}
				}
			}
		}
		else if (b > lower[3] && b <= upper[3])
		{
			last_last_reaction = last_reaction;
			last_reaction = 4;
			//cout << "reaction 4!" << endl;
			//reaction 4	kAB*cPA*cB;
			if (polymers[selected_chain].endGroup() == Atype)
			{
				// do reaction
				B = B - 1;
				//nAB = nAB + 1; nBB = nBB + 1;
				polymers[selected_chain].addMonomerToEnd(Btype);
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroup() != Atype)
				{
					//selected_chain = selectNewChain(selected_chain,chains);
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "\t4 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroup() == Atype)
					{
						// do reaction
						B = B - 1;
						//nAB = nAB + 1; nBB = nBB + 1;
						polymers[selected_chain].addMonomerToEnd(Btype);
						break;
					}
				}
			}
		}
		else if (b > lower[4] && b <= upper[4])
		{
			last_last_reaction = last_reaction;
			last_reaction = 5;
			// reaction 5	kBA*cPB*cA;
			//cout << "reaction 5!" << endl;
			if (polymers[selected_chain].endGroup() == Btype)
			{
				// do reaction
				A = A - 1;
				//nBA = nBA + 1; nAA = nAA + 1;
				polymers[selected_chain].addMonomerToEnd(Atype);
				polymers[selected_chain].addMonomerToEnd(Atype);
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroup() != Btype)
				{
					//selected_chain = selectNewChain(selected_chain,chains);
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "\t5 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroup() == Btype)
					{
						// do reaction
						A = A - 1;
						//nBA = nBA + 1; nAA = nAA + 1;
						polymers[selected_chain].addMonomerToEnd(Atype);
						polymers[selected_chain].addMonomerToEnd(Atype);
						break;
					}
				}
			}
		}
		else if (b > lower[5] && b <= upper[5])
		{
			last_last_reaction = last_reaction;
			last_reaction = 6;
			// reaction 6	kBB*cPB*cB
			//cout << "reaction 6!" << endl;
			if (polymers[selected_chain].endGroup() == Btype)
			{
				// do reaction
				B = B - 1;
				//nBB = nBB + 2;
				polymers[selected_chain].addMonomerToEnd(Btype);
			}
			else if (polymers[selected_chain].hasInitiated()) // find new chain
			{
				while (polymers[selected_chain].endGroup() != Btype)
				{
					//selected_chain = selectNewChain(selected_chain,chains);
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "\t6 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroup() == Btype)
					{
						// do reaction
						B = B - 1;
						//nBB = nBB + 2;
						polymers[selected_chain].addMonomerToEnd(Btype);
						break;
					}
				}
			}
		}
		else if (b > lower[6] && b <= upper[6])
		{
			last_last_reaction = last_reaction;
			last_reaction = 7;
			// reaction 7	kAmA*cPAAA
			//cout << "reaction 7!" << endl;
			if (polymers[selected_chain].endGroupTriad() == 111)
			{
				// do reaction
				A = A + 1;
				//nAA = nAA - 2;
				polymers[selected_chain].removeMonomerFromEnd();
				polymers[selected_chain].removeMonomerFromEnd();
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroupTriad() != 111) // this could be an issue
				{
					//selected_chain = selectNewChain(selected_chain,chains);
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "7 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroupTriad() == 111)
					{
						// do reaction
						A = A + 1;
						//nAA = nAA - 2;
						polymers[selected_chain].removeMonomerFromEnd();
						polymers[selected_chain].removeMonomerFromEnd();
						break;
					}
				}
			}
			if (polymers[selected_chain].getDegreeOfPolymerization()==0)
			{
				I = I + 1;
				//un_initiated.push_back(selected_chain);
			}
		}
		else if (b > lower[7] && b <= upper[7])
		{
			last_last_reaction = last_reaction;
			last_reaction = 8;
			// reaction 8	kAmB*cPABB
			//cout << "reaction 8!" << endl;
			if (polymers[selected_chain].endGroupDyad() == 12)
			{
				// do reaction
				B = B + 1;
				//nAB = nAB - 1; nBB = nBB - 1;
				polymers[selected_chain].removeMonomerFromEnd();
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroupDyad() != 12)
				{
					//selected_chain = selectNewChain(selected_chain,chains);
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "8 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroupDyad() == 12)
					{
						// do reaction
						B = B + 1;
						//nAB = nAB - 1; nBB = nBB - 1;
						polymers[selected_chain].removeMonomerFromEnd();
						break;
					}
				}
			}
			if (polymers[selected_chain].getDegreeOfPolymerization()==0)
			{
				I = I + 1;
				//un_initiated.push_back(selected_chain);
			}
		}
		else if (b > lower[8] && b <= upper[8])
		{
			last_last_reaction = last_reaction;
			last_reaction = 9;
			// reaction 9	kBmA*cPBAA
			//cout << "reaction 9!" << endl;
			if (polymers[selected_chain].endGroupTriad() == 211)
			{
				// do reaction
				A = A + 1;
				//nBA = nBA - 1; nAA = nAA - 1;
				polymers[selected_chain].removeMonomerFromEnd();
				polymers[selected_chain].removeMonomerFromEnd();
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroupTriad() != 211)
				{
					//selected_chain = selectNewChain(selected_chain, chains);
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "9 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroupTriad() == 211)
					{
						// do reaction
						A = A + 1;
						//nBA = nBA - 1; nAA = nAA - 1;
						polymers[selected_chain].removeMonomerFromEnd();
						polymers[selected_chain].removeMonomerFromEnd();
						break;
					}
				}
			}
		}
		else if (b > lower[9] && b <= upper[9])
		{
			last_last_reaction = last_reaction;
			last_reaction = 10;
			// reaction 10	kBmB*cPBBB
			//cout << "reaction 10!" << endl;
			if (polymers[selected_chain].endGroupDyad() == 22)
			{
				// do reaction
				B = B + 1;
				//nBB = nBB - 2;
				polymers[selected_chain].removeMonomerFromEnd();
			}
			else // find new chain
			{
				while (polymers[selected_chain].endGroupDyad() != 22)
				{
					//selected_chain = selectNewChain(selected_chain, chains);
					selected_chain = int(double(chains) * double(rando()) / rando_maximo);
					//polymers[selected_chain].printPolymer();
					//cout << "10 selected_chain = " << selected_chain << endl;
					if (polymers[selected_chain].endGroupTriad() == 22)
					{
						// do reaction
						B = B + 1;
						//nBB = nBB - 2;
						polymers[selected_chain].removeMonomerFromEnd();
						break;
					}
				}
			}
			if (polymers[selected_chain].getDegreeOfPolymerization()==0)
			{
				I = I + 1;
				//un_initiated.push_back(selected_chain);
			}
		}
		// general and random transesterification:
		else if (b > lower[10] && b <= upper[10])
		{
			last_last_reaction = last_reaction;
			last_reaction = 11;
			// NOTE - if you cut all repeat units off, or totally depolymerize a polymer,
			// you need to release an initiator I += 1;
			//cout << "reaction 11!" << endl;
			//cout << "selected chain (chain 1) is " << selected_chain << endl;
			repeats_removed = 0; repeats_added = 0;
			ester_sum = 0;
			ester_found = false;
			// random transesterification
			selected_ester = int(double(nE) * double(rando()) / rando_maximo);
			//cout << "selected_ester is " << selected_ester << endl;
			chain1 = selected_chain;
			// select chain1 and ester_id (find chain2 with ester_id)
			// find chain with selected_ester
			for (int i=0; i<=polymers.size()-1;i++)
			{
				ester_sum += polymers[i].getDegreeOfPolymerization();
				// we are searching for the polymer with 
				if (ester_sum >= selected_ester)
				{
					//cout << "ester sum before subtraction is" << ester_sum << endl;
					chain2 = i; // done
					// trim off the ith chain DP to reset ester coordinate to start of chain
					ester_sum -= polymers[i].getDegreeOfPolymerization();
					//cout << "degree of polymerization is " << polymers[chain2].getDegreeOfPolymerization() << endl;
					// now get index of ester - I think we can just calculate it
					ester_id = selected_ester - ester_sum; // do I need + 1 or - 1?
					//cout << "selected_ester = " << selected_ester <<  "\t" "ester_id ic = " << ester_id << " on chain " << chain2 << endl;
					ester_found = true;
					//cout << "ester_id is " << ester_id << "and ester_sum is " << ester_sum << endl;
					break;
				}
			}
			if (chain1 == chain2)
			{
				//Just pick a new chain1 because easter
				//cout << "same chain" << endl;
				while (chain1 == chain2)
				{
					chain1 = int(double(chains) * double(rando()) / rando_maximo);
					//cout << "looking for new chain" << endl;
					if (chain1 != chain2)
					{
						break;
					}
				}
			}

			if (polymers[chain1].hasInitiated() == false)
			{
				I = I-1;
				//cout << "initiator consumed" << endl;
			}

			if (ester_found) // it's on a linear chain A
			{
				//cout << "in second loop " << endl;
				for (int i=ester_id;i<=polymers[chain2].getDegreeOfPolymerization()-1;i++)
				{
					//cout << "stepping through adding units " << endl;
					if (polymers[chain2].repeatUnitAtPosition(i)==Atype)
					{
						polymers[chain1].addMonomerToEnd(Atype);
						repeats_added += 1;
					}
					else if (polymers[chain2].repeatUnitAtPosition(i)==Btype)
					{
						polymers[chain1].addMonomerToEnd(Btype);
						repeats_added += 1;
					}
					//cout << "adding to end" << endl;
				}
				//cout << "repeats added = " << repeats_added << endl;
				// trim chain2 from end to length transfered to selected_chain
				ester_sum = polymers[chain2].getDegreeOfPolymerization(); // reusing  variable
				for (int i=ester_id;i<=ester_sum-1;i++) // count up
				{
					// NOTE possibility for non-conservation of repeat units here
					polymers[chain2].removeMonomerFromEnd();
					repeats_removed += 1;
					//cout << "removing from end" << endl;
				}
				//cout << "repeats removed = " << repeats_removed << endl;
				if (polymers[chain2].getDegreeOfPolymerization()==0)
				{
					I = I + 1;
					//cout << "initiator released" << endl;
				}
				if (repeats_added != repeats_removed)
				{
					cout << "error in A " << endl;
				}
			}
			else
			{
				cout << "couldn't find ester" << endl;
			}
		}
		else
		{
			// throw error
			cout << " oh nonononononon shit shit shit!" << endl;
			cout << a << "\t" << b << "\t" << c << endl;
			cout << "selected chain = " << selected_chain << endl;
			cout << polymers[selected_chain].getDegreeOfPolymerization() << endl;
			//cout << lower[0] << lower[1] << lower[2]<<lower[3]<<lower[4]<<lower[5]<<lower[6]<<lower[7]<<lower[8]<<lower[9]<<lower[10]<<upper[10]<<endl;
		}

		if (counter <= 1000)
		{
			// save time
			Eigen::Vector::push_back<double>(time, t);
			// save conversion of A
			Eigen::Vector::push_back(conversionA, (double(nA0) - double(A)) / double(nA0));
			// save conversion of B
			Eigen::Vector::push_back(conversionB, (double(nB0) - double(B)) / double(nB0));
		}
		else if (counter > 1000)
		{
			if (counter % 1000 == 0)
			{
				// save time
				Eigen::Vector::push_back<double>(time, t);
				// save conversion of A
				Eigen::Vector::push_back(conversionA, (double(nA0) - double(A)) / double(nA0));
				// save conversion of B
				Eigen::Vector::push_back(conversionB, (double(nB0) - double(B)) / double(nB0));
				//cout << t << "\t" << I << "\t" << A << "\t" << B << "\t" << (double(nA0) - double(A)) / double(nA0) << "\t" << (double(nB0) - double(B)) / double(nB0) << endl;
			}
		}


		if (t > termination_time)
		{
			cout << "count is " << counter << endl;
			//Print out everything one last time to prevent errors in y_at_x
			// save time
			Eigen::Vector::push_back<double>(time, t);
			// save conversion of A
			Eigen::Vector::push_back(conversionA, (double(nA0) - double(A)) / double(nA0));
			// save conversion of B
			Eigen::Vector::push_back(conversionB, (double(nB0) - double(B)) / double(nB0));
			break;
		}	
		counter += 1;
		// trouble shoot output here
		
		//cout << t << "\t"<< I <<"\t" << A << "\t" << B << "\t" << (double(nA0)-double(A))/double(nA0) << "\t" << (double(nB0)-double(B))/double(nB0) << endl;
		
	}
	
	package = MatrixXd::Zero(time.size(),3);
	package.col(0) = time;
	package.col(1) = conversionA;
	package.col(2) = conversionB;
	
	return package;
}
/****************************************************************************************/
MatrixXd	stochastic2(VectorXd concentrations,int i_number,VectorXd parameters,double t_final)
{
	// 
	double c_initiator = concentrations(0);
	double c_monomer = concentrations(1);
	double kp = parameters(0);
	double kd = parameters(1);
	
	const double	avogadros = 6.022149e+23; // #/mole.
	// convert from concentrations to integer numbers of initiators
	double	volume = double(i_number)/(avogadros*c_initiator); // mol/L
	int		init0 = i_number;
	int		init = init0;
	int		mono0 = int(c_monomer * volume * avogadros);
	int		mono = mono0;
	int		poly = 0;
	double	rho_tot;
	double	cp = parameters(0)/avogadros/volume;
	double	cd = parameters(1);//*avogadros*volume; // not converted correctly??!!
	double	t = 0.;
	double	u1,u2;
	
	VectorXd	time = VectorXd::Zero(1);
	VectorXd	conversion = VectorXd::Zero(1);
	VectorXd	rho = VectorXd::Zero(2); // number of reactions
	VectorXd	prob = VectorXd::Zero(2); // number of parameters
	MatrixXd	package;
	
	mt19937		rando(time(0));
	double		rando_maximo = 4294967295;
	
	// initiation happens at t = 0;
	// model here:
	mono = mono0 - init;
	poly = init;
	init = 0;
	rho(0)	= cp * double(poly) * double(mono);	// second order rate constant
	rho(1)	= cd * double(poly); 				// first order rate constant
	rho_tot = rho(0) + rho(1);
	prob(0) = rho(0)/rho_tot;
	prob(1) = rho(1)/rho_tot;

	time(0) = t;
	conversion(0) = (double(mono0)-double(mono))/double(mono0);
	
	while (t<t_final)
	{
		u1 = double(rando())/rando_maximo; // time
		u2 = double(rando())/rando_maximo; // reaction selection
		// calculate time
		t = t + 1.0/rho_tot*log(1.0/u1);
		// store time
		Eigen::Vector::push_back<double>(time,t);
		if (u2 <= prob(0)) // propagation
		{
			mono = mono - 1;
		}
		else if (u2 > prob(0) && u2 < prob.sum()) // depropagation
		{
			mono = mono + 1;
		}
		// calculate conversion and store
		Eigen::Vector::push_back(conversion,(double(mono0)-double(mono))/double(mono0));
		rho(0)	= cp * double(poly) * double(mono);
		rho(1)	= cd * double(poly);
		rho_tot = rho(0) + rho(1);
		prob(0) = rho(0)/rho_tot;
		prob(1) = rho(1)/rho_tot;
	}
	package = MatrixXd::Zero(time.size(),2);
	package.col(0) = time;
	package.col(1) = conversion;
	return package;	
}
/****************************************************************************************/
double	residual(VectorXd x_data, VectorXd y_data, VectorXd x_model,VectorXd y_model)
{
	// x(i) and y(i) are the data
	// x and out are the model predictions
	// sum (y(t_i)-y_i)^2
	double sum=0.;
	for (int i=0;i<=x_data.size()-1;i++)
	{
		sum += pow(( y_at_x(x_data(i),x_model,y_model) - y_data(i) ),2.);
		//cout << t(i) << "\t" << y_at_x(t(i),x,out) << "\t" << y(i) << endl;
	}
	return sum;
}
/****************************************************************************************/
double	residual_time(VectorXd x_data,VectorXd y1_data,VectorXd y2_data,
				VectorXd x_model,VectorXd y1_model,VectorXd y2_model,double MonomerA, double MonomerB)
{
	// x(i) and y(i) are the data
	// x and out are the model predictions
	// sum (y(t_i)-y_i)^2
	//
	// to regress on conv vs A/A0 and B/B0 x_data -> y1_data then y2_data
	// y1_data -> conv_data_all and y2_data -> conv_data_all
	double sum=0.;
	for (int i=0;i<=x_data.size()-1;i++)
	{
		sum += pow(( y_at_x(x_data(i),x_model,y1_model) - y1_data(i) ),2.);
		cout << "convA is " << y1_data(i) << " " << y_at_x(x_data(i), x_model, y1_model) << endl;
	}
	for (int i=0;i<=x_data.size()-1;i++)
	{
		sum += pow(( y_at_x(x_data(i),x_model,y2_model) - y2_data(i) ),2.);
		cout << "convB is " << y2_data(i) << " " << y_at_x(x_data(i), x_model, y2_model) << endl;
	}
	for (int i = 0; i <= x_data.size() - 1; i++)
	{
		sum += pow((((MonomerA * (1 - y_at_x(x_data(i), x_model, y1_model))) / (MonomerA * (1 - y_at_x(x_data(i), x_model, y1_model)) + MonomerB * (1 - y_at_x(x_data(i), x_model, y2_model)))) - (MonomerA * (1 - y1_data(i)) / (MonomerA * (1 - y1_data(i)) + MonomerB * (1 - y2_data(i))))), 2.);
		cout << "fG data is " << (MonomerA * (1 - y1_data(i)) / (MonomerA * (1 - y1_data(i)) + MonomerB * (1 - y2_data(i)))) << "\t" << (MonomerA * (1 - y_at_x(x_data(i), x_model, y1_model))) / (MonomerA * (1 - y_at_x(x_data(i), x_model, y1_model)) + MonomerB * (1 - y_at_x(x_data(i), x_model, y2_model)))<< endl;
	}
	return sum;
}
/****************************************************************************************/
double	residual_compositional(VectorXd y1_data,VectorXd y2_data,VectorXd conv_data,
				VectorXd conv_model,VectorXd y1_model,VectorXd y2_model)
{
	// x(i) and y(i) are the data
	// x and out are the model predictions
	// sum (y(t_i)-y_i)^2
	//
	// to regress on conv vs A/A0 and B/B0 x_data -> y1_data then y2_data
	// y1_data -> conv_data_all and y2_data -> conv_data_all
	double sum=0.;
	for (int i=0;i<=y1_data.size()-1;i++)
	{
		sum += pow(( y_at_x(y1_data(i),y1_model,conv_model) - y1_data(i) ),2.);
	}
	for (int i=0;i<=y2_data.size()-1;i++)
	{
		sum += pow(( y_at_x(y2_data(i),y2_model,conv_model) - y2_data(i) ),2.);
	}
	return sum;
}
/****************************************************************************************/
VectorXd	conversion(VectorXd concentrations,VectorXd y1y10_data,VectorXd y2y20_data)
{
	double	c_initiator = concentrations(0);
	double	c_monomer_1 = concentrations(1);
	double	c_monomer_2 = concentrations(2);
	VectorXd conv = VectorXd::Zero(y1y10_data.size());
	double	f10 = concentrations(1)/(concentrations(1)+concentrations(2));
	double	f20 = concentrations(2)/(concentrations(1)+concentrations(2));
	conv = f10*y1y10_data+f20*y2y20_data;
	return conv;
}
/****************************************************************************************/
double	y_at_x(double at_x,VectorXd x_model,VectorXd y_model)
{
	// This function takes a numerically defined set of data and uses that to 
	// evaluate y(at_x) where at_x is a double valued input taken from an interpolant
	// modeling x_model vs. y_model. at_x lies within a few x_model points
	int			rel=0;
	MatrixXd	A = MatrixXd::Zero(3,3);
	VectorXd	c = VectorXd::Zero(3);
	VectorXd	b = VectorXd::Zero(3);
	
	
	// scan thru data 	until x(i-1) < at_x < x(i+1)	(increasing x_model)
	// 					or x(i-1) > at_x > x(i+1)		(decreasing x_model)
	for (int i=1;i<x_model.size()-1;i++)
	{
		//cout << at_x << "\t" << i << "\t" <<  x_model(i-1) << "\t" << x_model(i) << "\t" << x_model(i+1) << endl;
		// increasing x_model series, case #1
		if (x_model(i-1) <= at_x && at_x <= x_model(i+1))
		{
			// i is referenced to the above
			rel = i-1;
			// assemble interpolant matrix
			for (int j=0;j<=2;j++)
			{
				for (int k=0;k<=2;k++)
				{
					A(j,k) = pow(x_model(rel+j),double(k));
				}
				b(j) = y_model(rel+j);
			}
			//cout << "\n" << endl;
			//cout << A << endl;
			//cout << "\n" << endl;
			//cout << b << endl;
			// solve
			c = A.fullPivLu().solve(b);
			//cout << "\n" << endl;
			//cout << "coefficients = " << endl;
			//cout << c << endl;
			// return y(at_x)
			//cout << "---------" << endl;
			//cout << "y("<< at_x << ") = " << c(0) + at_x*c(1) + at_x*at_x*c(2) << endl;
			//cout << "from t and y" << endl;
			//cout << x_model(rel) << "\t" << y_model(rel) << endl;
			//cout << x_model(rel+1) << "\t" << y_model(rel+1) << endl;
			//cout << x_model(rel+2) << "\t" << y_model(rel+2) << endl;
			//cout << "---------------------------" << endl;
			//exit(0);
			return (c(0) + at_x*c(1) + at_x*at_x*c(2));
		}
		// or decreasing x_model series, case #2
		else if (x_model(i-1) >= at_x && at_x >= x_model(i+1))
		{
		// i is referenced to the above
			rel = i-1;
			// assemble interpolant matrix
			for (int j=0;j<=2;j++)
			{
				for (int k=0;k<=2;k++)
				{
					A(j,k) = pow(x_model(rel+j),double(k));
				}
				b(j) = y_model(rel+j);
			}
			//cout << "\n" << endl;
			//cout << A << endl;
			//cout << "\n" << endl;
			//cout << b << endl;
			// solve
			c = A.fullPivLu().solve(b);
			//cout << "\n" << endl;
			//cout << "coefficients = " << endl;
			//cout << c << endl;
			// return y(at_x)
			//cout << "---------" << endl;
			//cout << "y("<< at_x << ") = " << c(0) + at_x*c(1) + at_x*at_x*c(2) << endl;
			//cout << "from t and y" << endl;
			//cout << x_model(rel) << "\t" << y_model(rel) << endl;
			//cout << x_model(rel+1) << "\t" << y_model(rel+1) << endl;
			//cout << x_model(rel+2) << "\t" << y_model(rel+2) << endl;
			//cout << "---------------------------" << endl;
			//exit(0);
			return (c(0) + at_x*c(1) + at_x*at_x*c(2));
		}
		// or flat or symmetric
		else if (x_model(i-1) == at_x && at_x == x_model(i+1))
		{
			cout << "case #3 – not implemented yet" << endl;
		}
	}
}
/****************************************************************************************/



/****************************************************************************************/




/****************************************************************************************/
void get_psum(MatrixXd &p, VectorXd &psum, const int ihi) // Utility function.
{
	//should only add up everything but the highest point
	int ndim = psum.size();
	int	mpts = p.rows();
	for (int j=0;j<ndim;j++)
	{
		double sum=0.0;
		for (int i=0;i<mpts;i++)
		{
			if (i != ihi)
			{
				sum += p(i, j);
			}
		}
		psum[j]=sum;
	}
	//cout << "psum = " << psum << "ihi is " << ihi << endl;
}
/****************************************************************************************/
double amotry_time(MatrixXd &p, VectorXd &y, VectorXd &psum, const int ihi, const double fac,
			VectorXd concentrations,int number,VectorXd x_data,VectorXd y1_data,VectorXd y2_data) // T &func for general function
	// Helper function: Extrapolates by a factor fac through the face of the 
	// simplex across from the high point, tries it, and replaces the high point 
	// if the new point is better.
{
	int 		mpts=p.rows();
	int			ndim=p.cols();
	VectorXd	ptry = VectorXd::Zero(ndim);
	
	VectorXd	time_model;
	VectorXd	conv_model_A,conv_model_B;
	//VectorXd	parameters_new; // simplex will return this
	MatrixXd	model_data_package;
	
	double		t_final = x_data(x_data.size()-1);
	double fac1=(1.0-fac)/(mpts-1);
	double ytry;
	double fac2=fac;
	cout << "ptry is ";
	for (int j=0;j<ndim;j++)
	{
		ptry(j)=(psum(j))*fac1+p(ihi,j)*fac2;
		if (ptry(j) < 0) // make any negatives 0 instead
		{
			ptry(j) = 0;
		}
		cout << ptry(j) << " ";
	}
	cout << endl;
	model_data_package = stochastic8(concentrations,number,ptry,t_final);
	time_model =  model_data_package.col(0);
	conv_model_A =  model_data_package.col(1);
	conv_model_B =  model_data_package.col(2);
	// Evaluate the function at the trial point.
	ytry = residual_time(x_data,y1_data,y2_data,time_model,conv_model_A,conv_model_B,concentrations(1), concentrations(2));
	
	if (ytry < y(ihi))
	{	// If it’s better than the highest, then replace the highest.
		y(ihi)=ytry;
		for (int j=0;j<ndim;j++)
		{
			//psum(j) += ptry(j)-p(ihi,j);
			p(ihi,j)=ptry(j);
		}
	}
	return ytry;
}
/****************************************************************************************/
double amotry_composition(MatrixXd &p, VectorXd &y, VectorXd &psum, const int ihi, const double fac,
			VectorXd concentrations,int number,VectorXd x_data,VectorXd conv_data_A,VectorXd conv_data_B) // T &func for general function
	// Helper function: Extrapolates by a factor fac through the face of the 
	// simplex across from the high point, tries it, and replaces the high point 
	// if the new point is better.
{
	int 		mpts=p.rows();
	int			ndim=p.cols();
	VectorXd	ptry = VectorXd::Zero(ndim);
	
	VectorXd	time_model;
	VectorXd	conv_model_A,conv_model_B,conv_model_all;
	static VectorXd	conv_data_all = conversion(concentrations,conv_data_A,conv_data_B);
	MatrixXd	model_data_package;
	
	double		t_final = x_data(x_data.size()-1);
	double fac1=(1.0-fac)/ndim;
	double ytry;
	double fac2=fac1-fac;
	for (int j=0;j<ndim;j++)
	{
		ptry(j)=psum(j)*fac1-p(ihi,j)*fac2;
	}
	model_data_package = stochastic8(concentrations,number,ptry,t_final);
	time_model =  model_data_package.col(0);
	conv_model_A =  model_data_package.col(1);
	conv_model_B =  model_data_package.col(2);
	conv_model_all = conversion(concentrations,conv_model_A,conv_model_B);
	ytry = 	residual_compositional(conv_data_A,conv_data_B,conv_data_all,conv_model_all,conv_model_A,conv_model_B);
	
	if (ytry < y(ihi))
	{	// If it’s better than the highest, then replace the highest.
		y(ihi)=ytry;
		for (int j=0;j<ndim;j++)
		{
			psum(j) += ptry(j)-p(ihi,j);
			p(ihi,j)=ptry(j);
		}
	}
	return ytry;
}
/****************************************************************************************/
void		swap(double &a,double &b)
{
	double dummy = a;
		a = b;
		b = dummy;
	return;
}
/****************************************************************************************/
VectorXd	simplex_time(VectorXd parameters_old,VectorXd concentrations,int number,
				VectorXd x_data,VectorXd y1_data,VectorXd y2_data,double accuracy_goal)
{
	int n_params = parameters_old.size();
	// simplex initial displacements
	VectorXd	dels = VectorXd::Ones(n_params); // set to 0.1 initially
	double		func_at_min;
	// stochastic calculation populates this:
	VectorXd	time_model;
	VectorXd	conv_model_A,conv_model_B;
	VectorXd	parameters_new; // simplex will return this
	MatrixXd	model_data_package;
	
	int			ndim = n_params; // relate outside to inside routine.
	double		t_final = x_data(x_data.size()-1);
	
	// this may be input in the original code
	MatrixXd	pp = MatrixXd::Zero(ndim+1,ndim);
	
//////////////////////////////////////
	//template <class T>
	//VecDoub minimize(MatDoub_I &pp, T &func)
	// Most general interface: initial simplex specified by the matrix 
	// pp[0..ndim][0..ndim-1]. Its ndim+1 rows are ndim-dimensional vectors 
	// that are the vertices of the starting simplex.
	//{
	const int NMAX = 5000; //const Int NMAX=5000; // maximum allowed function evaluations
	const double TINY=1.e-10; //const Doub TINY=1.0e-10;
	double ftol = accuracy_goal; //2.0*abs(accuracy_goal-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY);
	int ihi,ilo,inhi; //Int ihi,ilo,inhi;
	int mpts=pp.rows();
	int nfunc=0;
	//int ndim=pp.cols();
	
	// need to set up the initial simplex
	for (int i=0;i<=pp.rows()-1;i++)
	{
		for (int j=0;j<=pp.cols()-1;j++)
		{
			pp(i,j) = parameters_old(j);
		}		
	}
	//cout << "Successfully made it into simplex time" << endl;
	cout << pp.rows() << ", " << pp.cols() << endl;
	//for (int i=0;i<=parameters_old.size()-1;i++)
	//{
	//	pp(i,i) += 0.025*pow(-1,(i+2))*dels(i);	
	//}
	//pp(0,0)  = ; pp(0,1)  = ; pp(0,2)  = ; pp(0,3)  = ; pp(0,4)  = ; pp(0,5)  = ; pp(0,6)  = ; pp(0,7)  = ; pp(0,8)  = ;
	//Only add to ensure that we start with positive values (although maybe not necessary if we add in checks later)
	pp(1,0) += 2.0; pp(1,1) += 0.0; pp(1,2) += 0.0; pp(1,3) += 0.0; pp(1,4) += 0.0; pp(1,5) += 0.0; pp(1,6) += 0.0; pp(1,7) += 0.0; pp(1,8) += 0.0;
	pp(2,0) += 0.0; pp(2,1) += 0.5; pp(2,2) += 0.0; pp(2,3) += 0.0; pp(2,4) += 0.0; pp(2,5) += 0.0; pp(2,6) += 0.0; pp(2,7) += 0.0; pp(2,8) += 0.0;
	pp(3,0) += 0.0; pp(3,1) += 0.0; pp(3,2) += 0.5; pp(3,3) += 0.0; pp(3,4) += 0.0; pp(3,5) += 0.0; pp(3,6) += 0.0; pp(3,7) += 0.0; pp(3,8) += 0.0;
	pp(4,0) += 0.0; pp(4,1) += 0.0; pp(4,2) += 0.0; pp(4,3) += 1.0; pp(4,4) += 0.0; pp(4,5) += 0.0; pp(4,6) += 0.0; pp(4,7) += 0.0; pp(4,8) += 0.0;
	pp(5,0) += 0.0; pp(5,1) += 0.0; pp(5,2) += 0.0; pp(5,3) += 0.0; pp(5,4) += 0.5; pp(5,5) += 0.0; pp(5,6) += 0.0; pp(5,7) += 0.0; pp(5,8) += 0.0;
	pp(6,0) += 0.0; pp(6,1) += 0.0; pp(6,2) += 0.0; pp(6,3) += 0.0; pp(6,4) += 0.0; pp(6,5) += 0.5; pp(6,6) += 0.0; pp(6,7) += 0.0; pp(6,8) += 0.0;
	pp(7,0) += 0.0; pp(7,1) += 0.0; pp(7,2) += 0.0; pp(7,3) += 0.0; pp(7,4) += 0.0; pp(7,5) += 0.0; pp(7,6) += 0.5; pp(7,7) += 0.0; pp(7,8) += 0.0;
	pp(8,0) += 0.0; pp(8,1) += 0.0; pp(8,2) += 0.0; pp(8,3) += 0.0; pp(8,4) += 0.0; pp(8,5) += 0.0; pp(8,6) += 0.0; pp(8,7) += 0.5; pp(8,8) += 0.0;
	pp(9,0) += 0.0; pp(9,1) += 0.0; pp(9,2) += 0.0; pp(9,3) += 0.0; pp(9,4) += 0.0; pp(9,5) += 0.0; pp(9,6) += 0.0; pp(9,7) += 0.0; pp(9,8) += 0.01;
	
	cout << "starting symplexes: " << endl;
	cout << pp << endl;
	VectorXd	psum = VectorXd::Zero(ndim);
	VectorXd	pmin = VectorXd::Zero(ndim);
	VectorXd	x = VectorXd::Zero(ndim);
	VectorXd	y = VectorXd::Zero(mpts);
	MatrixXd	p=pp; // huh?
	//y.resize(mpts);
	for (int i=0;i<mpts;i++)
	{
		for (int j=0;j<ndim;j++)
		{
			x(j)=p(i,j);
		}
		cout << "now doing parameters: " << x.transpose() << endl;
		// call stochastic
		// calculate residual, set equation to y(i)
		parameters_new = x; // there are many opportunities for simplification... later
		// do stochastic integration of time_model and conv_model
		model_data_package = stochastic8(concentrations,number,parameters_new,t_final);
		time_model =  model_data_package.col(0);
		conv_model_A =  model_data_package.col(1);
		conv_model_B =  model_data_package.col(2);
		y(i) = residual_time(x_data,y1_data,y2_data,time_model,conv_model_A,conv_model_B,concentrations(1),concentrations(2));
		//y(i)=func(x);
		cout << "... yields residual " << y(i) << endl;
	}
	//get_psum(p,psum,ihi); I think we need to get psum after we find ihi
	for (;;)
	{
		ilo=0;
		// First we must determine which point is the highest (worst), 
		// next-highest, and lowest (best), 
		// by looping over the points in the simplex.
		cout << "looping over points in the simplex to find best point" << endl;
		ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
		for (int i=0;i<mpts;i++)
		{
			if (y[i] <= y[ilo]) 
			{
				ilo=i;
			}
			if (y[i] > y[ihi])
			{
				inhi=ihi;
				ihi=i;
			} 
			else if (y[i] > y[inhi] && i != ihi)
			{
				inhi=i;
			}
		}
		get_psum(p, psum, ihi);
		double rtol = 0;//2.0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY); // change to when the actual parameters are close together
		for (int j = 0; j < ndim; j++)
		{
			rtol += abs((p(ihi, j) - p(ilo,j))/(mpts-1)); //issue here if parameters are converging to 0
		}
		rtol = rtol / mpts; // average percent difference between P average and p(ihi)
		cout << "current y(ilo) = " << y(ilo) << endl;
		cout << "current rtol = " << rtol << " is less than precision goal? " << ftol << endl;
		cout << "Current ilo parameters are " << endl;
		for (int ii = 0; ii <= ndim-1;ii++)
		{
			cout << p(ilo, ii) << " ";
		}
		cout << endl;
		// Compute the fractional range from highest to lowest and return if satisfactory.
		if (rtol < ftol) //(y(ilo) < ftol)
		{ // If returning, put best point and value in slot 0.
			swap(y[0],y[ilo]);
			for (int i=0;i<ndim;i++) 
			{
				swap(p(0,i),p(ilo,i));
				pmin(i)=p(0,i);
			}
			func_at_min=y[0];
			cout << "func_at_min = " << func_at_min << endl;
			return pmin;
		}
		if (nfunc >= NMAX)
		{
			cout << "NMAX exceeded" << endl;
			throw("NMAX exceeded");
		}
		nfunc += 2;
		// Begin a new iteration. First extrapolate by a factor –1 through 
		// the face of the simplex across from the high point, 
		// i.e., reflect the simplex from the high point.
			//Doub ytry=amotry(p,y,psum,ihi,-1.0,func);
		double ytry = amotry_time(p,y,psum,ihi,-1.0,concentrations,number,x_data,y1_data,y2_data);
		cout << "reflect ytry = " << ytry << endl;
		if (ytry <= y[ilo])
		{
		// Gives a result better than the best point, so try an additional extrapolation by a factor 2.
			ytry = amotry_time(p,y,psum,ihi,2.0,concentrations,number,x_data,y1_data,y2_data);
			cout << "extrapolation ytry = " << ytry << endl;
		}
		else if (ytry >= y[inhi])
		{
		// The reflected point is worse than the second-highest, so look for an 
		// intermediate lower point, i.e., do a one-dimensional contraction.
			double ysave=y[ihi];
			
			//ytry=amotry(p,y,psum,ihi,0.5,func);
			ytry = amotry_time(p,y,psum,ihi,0.5,concentrations,number,x_data,y1_data,y2_data);
			cout << "contraction ytry = " << ytry << endl;
			if (ytry >= ysave)
			{ // Can’t seem to get rid of that high point.
				for (int i=0;i<mpts;i++)
				{ // Better contract around the lowest
					if (i != ilo)
					{ // (best) point.
						for (int j=0;j<ndim;j++)
						{
							p(i,j)=psum(j)=0.5*(p(i,j)+p(ilo,j)); // not sure why psum is here- feel like that might be a problem
						}
						//y[i]=func(psum);
						parameters_new = psum; // there are many opportunities for simplification... later
						// do stochastic integration of time_model and conv_model
						model_data_package = stochastic8(concentrations,number,parameters_new,t_final);
						time_model =  model_data_package.col(0);
						conv_model_A =  model_data_package.col(1);
						conv_model_B =  model_data_package.col(2);
						y(i) = residual_time(x_data,y1_data,y2_data,time_model,conv_model_A,conv_model_B,concentrations(1),concentrations(2));
						cout << "y(i) = " << y(i) << endl;
					}
				}
				nfunc += ndim;		// Keep track of function evaluations.
				get_psum(p,psum,ihi);	// Recompute psum.
			}
		} else --nfunc;				// Correct the evaluation count.
	}								// Go back for the test of doneness and the next
//}									// iteration.
///////////////////////////
}
/****************************************************************************************/
VectorXd	simplex_composition(VectorXd parameters_old,VectorXd concentrations,int number,
				VectorXd x_data,VectorXd conv_data_A,VectorXd conv_data_B,double accuracy_goal)
{
	// simplex initial displacements
	VectorXd	dels = VectorXd::Ones(parameters_old.size()); // set to 0.1 initially
	double		func_at_min;
	// stochastic calculation populates this:
	VectorXd	time_model;
	VectorXd	conv_model_A,conv_model_B,conv_model_all;
	VectorXd	conv_data_all = conversion(concentrations,conv_data_A,conv_data_B);
	VectorXd	parameters_new; // simplex will return this
	MatrixXd	model_data_package;
	
	int			ndim = parameters_old.size(); // relate outside to inside routine.
	double		t_final = x_data(x_data.size()-1); // this is needed for calculation
	
	// this may be input in the original code
	MatrixXd	pp = MatrixXd::Zero(ndim+1,ndim);
	
//////////////////////////////////////
	//template <class T>
	//VecDoub minimize(MatDoub_I &pp, T &func)
	// Most general interface: initial simplex specified by the matrix 
	// pp[0..ndim][0..ndim-1]. Its ndim+1 rows are ndim-dimensional vectors 
	// that are the vertices of the starting simplex.
	//{
	const int NMAX = 3; //const Int NMAX=5000; // maximum allowed function evaluations
	const double TINY=1.e-10; //const Doub TINY=1.0e-10;
	double ftol = accuracy_goal; //2.0*abs(accuracy_goal-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY);
	int ihi,ilo,inhi; //Int ihi,ilo,inhi;
	int mpts=pp.rows();
	int nfunc=0;
	//int ndim=pp.cols();
	
	// need to set up the initial simplex
	for (int i=0;i<=pp.rows()-1;i++)
	{
		for (int j=0;j<=pp.cols()-1;j++)
		{
			pp(i,j) = parameters_old(j);
		}		
	}
	cout << "In composition for some reason" << endl;
	cout << pp.rows() << ", " << pp.cols() << endl;
	//for (int i=0;i<=parameters_old.size()-1;i++)
	//{
	//	pp(i,i) += 0.025*pow(-1,(i+2))*dels(i);	
	//}
	//pp(0,0)  = ; pp(0,1)  = ; pp(0,2)  = ; pp(0,3)  = ; pp(0,4)  = ; pp(0,5)  = ; pp(0,6)  = ; pp(0,7)  = ; pp(0,8)  = ;
	pp(1,0) += 0.3; pp(1,1) -= 0.0; pp(1,2) += 0.0; pp(1,3) += 0.0; pp(1,4) += 0.0; pp(1,5) += 0.0; pp(1,6) += 0.0; pp(1,7) += 0.0; pp(1,8) += 0.0;
	pp(2,0) -= 0.0; pp(2,1) += 0.3; pp(2,2) += 0.0; pp(2,3) += 0.0; pp(2,4) += 0.0; pp(2,5) += 0.0; pp(2,6) += 0.0; pp(2,7) += 0.0; pp(2,8) += 0.0;
	pp(3,0) += 0.0; pp(3,1) += 0.0; pp(3,2) -= 0.3; pp(3,3) += 0.0; pp(3,4) += 0.0; pp(3,5) += 0.0; pp(3,6) += 0.0; pp(3,7) += 0.0; pp(3,8) += 0.0;
	pp(4,0) += 0.0; pp(4,1) -= 0.0; pp(4,2) += 0.0; pp(4,3) += 0.3; pp(4,4) += 0.0; pp(4,5) += 0.0; pp(4,6) += 0.0; pp(4,7) += 0.0; pp(4,8) += 0.0;
	pp(5,0) += 0.0; pp(5,1) += 0.0; pp(5,2) += 0.0; pp(5,3) -= 0.0; pp(5,4) += 0.3; pp(5,5) += 0.0; pp(5,6) += 0.0; pp(5,7) += 0.0; pp(5,8) += 0.0;
	pp(6,0) += 0.0; pp(6,1) += 0.0; pp(6,2) -= 0.0; pp(6,3) += 0.0; pp(6,4) += 0.0; pp(6,5) += 0.3; pp(6,6) += 0.0; pp(6,7) += 0.0; pp(6,8) += 0.0;
	pp(7,0) += 0.0; pp(7,1) += 0.0; pp(7,2) += 0.0; pp(7,3) += 0.0; pp(7,4) += 0.0; pp(7,5) += 0.0; pp(7,6) += 0.3; pp(7,7) += 0.0; pp(7,8) += 0.0;
	pp(8,0) += 0.0; pp(8,1) += 0.0; pp(8,2) += 0.0; pp(8,3) += 0.0; pp(8,4) += 0.0; pp(8,5) += 0.0; pp(8,6) += 0.0; pp(8,7) += 0.3; pp(8,8) += 0.0;
	pp(9,0) += 0.0; pp(9,1) += 0.0; pp(9,2) += 0.0; pp(9,3) += 0.0; pp(9,4) += 0.0; pp(9,5) -= 0.0; pp(9,6) += 0.0; pp(9,7) += 0.0; pp(9,8) += 0.01;
	
	cout << "starting symplexes: " << endl;
	cout << pp << endl;
	//exit(0);
	
	VectorXd	psum = VectorXd::Zero(ndim);
	VectorXd	pmin = VectorXd::Zero(ndim);
	VectorXd	x = VectorXd::Zero(ndim);
	VectorXd	y = VectorXd::Zero(mpts);
	MatrixXd	p=pp; // huh?
	//y.resize(mpts);
	for (int i=0;i<mpts;i++)
	{
		for (int j=0;j<ndim;j++)
		{
			x(j)=p(i,j);
		}
		cout << "now doing parameters: " << x.transpose() << endl;
		// call stochastic
		// calculate residual, set equation to y(i)
		parameters_new = x; // there are many opportunities for simplification... later
		// do stochastic integration of time_model and conv_model
		model_data_package = stochastic8(concentrations,number,parameters_new,t_final);
		time_model =  model_data_package.col(0);
		conv_model_A =  model_data_package.col(1);
		conv_model_B =  model_data_package.col(2);
		conv_model_all = conversion(concentrations,conv_model_A,conv_model_B);
		y(i) = 	residual_compositional(conv_data_A,conv_data_B,conv_data_all,conv_model_all,conv_model_A,conv_model_B);
		//y(i) = residual2(x_data,y1_data,y2_data,time_model,conv_model_A,conv_model_B);
		
		//y(i)=func(x);
		cout << "... yields residual " << y(i) << endl;
	}
	get_psum(p,psum,ihi);
	for (;;)
	{
		ilo=0;
		// First we must determine which point is the highest (worst), 
		// next-highest, and lowest (best), 
		// by looping over the points in the simplex.
		cout << "looping over points in the simplex to find best point" << endl;
		ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
		for (int i=0;i<mpts;i++)
		{
			if (y[i] <= y[ilo]) 
			{
				ilo=i;
			}
			if (y[i] > y[ihi])
			{
				inhi=ihi;
				ihi=i;
			} 
			else if (y[i] > y[inhi] && i != ihi)
			{
				inhi=i;
			}
		}
		//double rtol=2.0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY);
		double rtol = 0;
		for (int i = 0; i < ndim; i++)
		{
			rtol += abs(p(ihi, i) - p(ilo, i));
			cout << p(ihi, i) << " " << p(ilo, i) << " " << abs(p(ihi, i) - p(ilo, i)) << endl;
		}
		rtol = rtol / (ndim-1);
		cout << ndim - 1 << " " << rtol << endl;
		cout << "current rtol = " << rtol << " is less than accuracy goal? " << ftol << endl;
		cout << "y(ilo) = " << y(ilo) << endl;
		// Compute the fractional range from highest to lowest and return if satisfactory.
		if (rtol < ftol) //(rtol < ftol)
		{ // If returning, put best point and value in slot 0.
			swap(y[0],y[ilo]);
			for (int i=0;i<ndim;i++) 
			{
				swap(p(0,i),p(ilo,i));
				pmin(i)=p(0,i);
			}
			func_at_min=y[0];
			cout << "func_at_min = " << func_at_min << endl;
			return pmin;
		}
		if (nfunc >= NMAX) throw("NMAX exceeded");
		nfunc += 2;
		// Begin a new iteration. First extrapolate by a factor –1 through 
		// the face of the simplex across from the high point, 
		// i.e., reflect the simplex from the high point.
			//Doub ytry=amotry(p,y,psum,ihi,-1.0,func);
		double ytry = amotry_composition(p,y,psum,ihi,-1.0,concentrations,number,x_data,conv_data_A,conv_data_B);
		cout << "reflect ytry = " << ytry << endl;
		if (ytry <= y[ilo])
		{
		// Gives a result better than the best point, so try an additional extrapolation by a factor 2.
			ytry = amotry_composition(p,y,psum,ihi,2.0,concentrations,number,x_data,conv_data_A,conv_data_B);
			cout << "extrapolation ytry = " << ytry << endl;
		}
		else if (ytry >= y[inhi])
		{
		// The reflected point is worse than the second-highest, so look for an 
		// intermediate lower point, i.e., do a one-dimensional contraction.
			double ysave=y[ihi];
			
			//ytry=amotry(p,y,psum,ihi,0.5,func);
			ytry = amotry_composition(p,y,psum,ihi,0.5,concentrations,number,x_data,conv_data_A,conv_data_B);
			cout << "contraction ytry = " << ytry << endl;
			if (ytry >= ysave)
			{ // Can’t seem to get rid of that high point.
				for (int i=0;i<mpts;i++)
				{ // Better contract around the lowest
					if (i != ilo)
					{ // (best) point.
						for (int j=0;j<ndim;j++)
						{
							p(i,j)=psum(j)=0.5*(p(i,j)+p(ilo,j));
						}
						//y[i]=func(psum);
						parameters_new = psum; // there are many opportunities for simplification... later
						// do stochastic integration of time_model and conv_model
						model_data_package = stochastic8(concentrations,number,parameters_new,t_final);
						time_model =  model_data_package.col(0);
						conv_model_A =  model_data_package.col(1);
						conv_model_B =  model_data_package.col(2);
						conv_model_all = conversion(concentrations,conv_model_A,conv_model_B);
						y(i) = 	residual_compositional(conv_data_A,conv_data_B,conv_data_all,conv_model_all,conv_model_A,conv_model_B);
						cout << "y(i) = " << y(i) << endl;
					}
				}
				nfunc += ndim;		// Keep track of function evaluations.
				get_psum(p,psum,ihi);	// Recompute psum.
			}
		} else --nfunc;				// Correct the evaluation count.
	}								// Go back for the test of doneness and the next
//}									// iteration.
///////////////////////////
}

/****************************************************************************************/
double	linear_optim(double x1,double x2,double fx1,double fx2)
{
	// this is not the best because it will just as easily go to a maximum as a minimum
	return (-(fx2*x1) + fx1*x2)/(fx1 - fx2);
}
/****************************************************************************************/
double	brents_optim(double x1,double x2,double x3,double fx1,double fx2,double fx3)
{
	// this is not the best because it will just as easily go to a maximum as a minimum
	return ((fx3*(-pow(x1,2) + pow(x2,2)) + fx2*(pow(x1,2) - pow(x3,2)) 
			+ fx1*(-pow(x2,2) + pow(x3,2)))/(2.*(fx3*(-x1 + x2) + fx2*(x1 - x3) 
			+ fx1*(-x2 + x3))));
}
/****************************************************************************************/
bool	is_it_a_minimum(VectorXd lambda,VectorXd residual)
{
	double fact;
	bool goes_to_min;
	// this is not the best because it will just as easily go to a maximum as a minimum
	fact = (2*((residual(0) - residual(2))/(lambda(0) - lambda(2)) 
	+ (-residual(1) + residual(2))/(lambda(1) - lambda(2))))/(lambda(0) - lambda(1));
	if (fact < 0)
	{
		goes_to_min = false;
	}
	else if (fact > 0)
	{
		goes_to_min = true;
	}
	return goes_to_min;
}
/****************************************************************************************/
VectorXd linspace(double low,double high,int points)
{
	VectorXd x = VectorXd::Zero(points);
	for (int i=0;i<=points-1;i++)
	{
		x(i) = (high-low)/(double(points)-1.) * double(i) + low;
	}
	return x;
}
/****************************************************************************************/
/****************************************************************************************/
double doubleFromStringBuffer(string buffer)
{
	double	number	= 0.;
	int		address = buffer.find("=",0);
	address = buffer.find_first_not_of(" ",address+1);
	stringstream strm(buffer.substr(address,buffer.size()-address));
	strm >> number;
	return number;
}
/****************************************************************************************/
int integerFromStringBuffer(string buffer)
{
	int		number	= 0.;
	int		address = buffer.find("=",0);
	address = buffer.find_first_not_of(" ",address+1);
	stringstream strm(buffer.substr(address,buffer.size()-address));
	strm >> number;
	return number;
}
/****************************************************************************************/
string stringAfterEqualsSign(string buffer)
{
	int address = buffer.find("=",0);
	address = buffer.find_first_not_of(" ",address+1);
	return buffer.substr(address,buffer.size()-address);
}
/****************************************************************************************/
string stringBetweenBrackets(string buffer)
{
	int address1 = buffer.find("{",0);
	int address2 = buffer.find("}",0);
	return buffer.substr(address1+1,address2-address1-1);
}
/****************************************************************************************/
double doubleFromString(string substring)
{
	double number;
	stringstream strm(substring);
	strm >> number;
	return number;
}
/****************************************************************************************/
bool existsInString(string test,string buffer)
{
	if (buffer.find(test,0)!=string::npos)
	{
		return true;
	}
	else
	{
		return false;
	}
}
/****************************************************************************************/
int	selectNewChain(int selected_chain, int chains)
{
	if (selected_chain < (chains-1))
	{
		return (selected_chain + 1);
	}
	else if (selected_chain == (chains-1))
	{
		return 0;
	}
	else
	{
		cout << "selectNewChain error at " << selected_chain << endl;
		return 0;
	}
}
/****************************************************************************************/
double	conversion(int A,int nA0,int B,int nB0)
{
	return ((double(nA0-A)+double(nB0-B))/double(nA0+nB0));
}
/****************************************************************************************/

/****************************************************************************************/