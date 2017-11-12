/*
   Program to solve the two-dimensional Ising model
   with zero external field and no parallelization
   The coupling constant J is set to J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis aolgorithm  is used as well as periodic boundary conditions.
   The code needs an output file on the command line and the variables mcs, nspins,
   initial temp, final temp and temp step.
   Run as
   ./executable Outputfile numberof spins number of MC cycles initial temp final temp tempstep
   ./test.x Lattice 100 10000000 2.1 2.4 0.01
   Compile and link as
   c++ -O3 -std=c++11 -Rpass=loop-vectorize -o Ising.x IsingModel.cpp -larmadillo
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
using namespace  std;
using namespace arma;
// output file
ofstream ofile;

// inline function for PeriodicBoundary boundary conditions
inline int PeriodicBoundary(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}
// Function to initialise energy and magnetization
void InitializeLattice(int, mat &, double&, double&, long *idum);
// The metropolis algorithm including the loop over Monte Carlo cycles
void MetropolisSampling(int, int, double, vec &, long idum);
// prints to file the results of the calculations
void WriteResultstoFile(int, int, double, vec);

double ran2(long *);

// Main program begins here

int main(int argc, char* argv[])
{
  string filename;
  int NSpins, MCcycles;
  double InitialTemp, FinalTemp, TempStep;
  if (argc <= 5) {
    cout << "Bad Usage: " << argv[0] <<
      " read output file, Number of spins, MC cycles, initial and final temperature and tempurate step" << endl;
    exit(1);
  }
  if (argc > 1) {
    filename=argv[1];
    NSpins = atoi(argv[2]);
    MCcycles = atoi(argv[3]);
    InitialTemp = atof(argv[4]);
    FinalTemp = atof(argv[5]);
    TempStep = atof(argv[6]);
  }
  long idum=-1;
  // Declare new file name and add lattice size to file name
  string fileout = filename;
  string argument = to_string(NSpins);
  fileout.append(argument);
  ofile.open(fileout);
  // Start Monte Carlo sampling by looping over the selcted Temperatures
  for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){
    vec ExpectationValues = zeros<mat>(5);
    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(NSpins, MCcycles, Temperature, ExpectationValues, idum);
    //
    //WriteResultstoFile(NSpins, MCcycles, Temperature, ExpectationValues);
  }
  ofile.close();  // close output file
  return 0;
}



// The Monte Carlo part with the Metropolis algo with sweeps over the lattice
void MetropolisSampling(int NSpins, int MCcycles, double Temperature, vec &ExpectationValues, long idum)
{
  // Initialize the seed and call the Mersienne algo
  std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  // Initialize the lattice spin values
  mat SpinMatrix = zeros<mat>(NSpins,NSpins);
  //    initialize energy and magnetization
  double Energy = 0.;     double MagneticMoment = 0.;
  // initialize array for expectation values
  InitializeLattice(NSpins, SpinMatrix, Energy, MagneticMoment, &idum);
  // setup array for possible energy changes
  vec EnergyDifference = zeros<mat>(17);
  for( int de =-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/Temperature);
  // Start Monte Carlo cycles
  for (int cycles = 1; cycles <= MCcycles; cycles++){
      if (cycles%1000==0){
          WriteResultstoFile(NSpins, cycles, Temperature, ExpectationValues);
      }
    // The sweep over the lattice, looping over all spin sites
    for(int x =0; x < NSpins; x++) {
      for (int y= 0; y < NSpins; y++){
    int ix = (int) (RandomNumberGenerator(gen)*(double)NSpins);
    int iy = (int) (RandomNumberGenerator(gen)*(double)NSpins);
    int deltaE =  2*SpinMatrix(ix,iy)*
      (SpinMatrix(ix,PeriodicBoundary(iy,NSpins,-1))+
       SpinMatrix(PeriodicBoundary(ix,NSpins,-1),iy) +
       SpinMatrix(ix,PeriodicBoundary(iy,NSpins,1)) +
       SpinMatrix(PeriodicBoundary(ix,NSpins,1),iy));
    if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) ) {
      SpinMatrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
      MagneticMoment += (double) 2*SpinMatrix(ix,iy);
      Energy += (double) deltaE;
    }
      }
    }
    // update expectation values  for local node
    ExpectationValues(0) += Energy;    ExpectationValues(1) += Energy*Energy;
    ExpectationValues(2) += MagneticMoment;
    ExpectationValues(3) += MagneticMoment*MagneticMoment;
    ExpectationValues(4) += fabs(MagneticMoment);
  }
} // end of Metropolis sampling over spins

// function to initialise energy, spin matrix and magnetization
void InitializeLattice(int NSpins, mat &SpinMatrix,  double& Energy, double& MagneticMoment, long *idum)
{
  // setup spin matrix and initial magnetization
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
        if(ran2(idum)<0.5){
            SpinMatrix(x,y) = 1.0; // spin orientation for the ground state
        }else{
            SpinMatrix(x,y) = -1.0;
        }
        cout <<SpinMatrix(x,y)<< " ";

      MagneticMoment +=  (double) SpinMatrix(x,y);
    }
    cout <<"\n";
  }
  // setup initial energy
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      Energy -=  (double) SpinMatrix(x,y)*
    (SpinMatrix(PeriodicBoundary(x,NSpins,-1),y) +
     SpinMatrix(x,PeriodicBoundary(y,NSpins,-1)));
    }
  }
}// end function initialise




void WriteResultstoFile(int NSpins, int MCcycles, double temperature, vec ExpectationValues)
{
  double norm = 1.0/((double) (MCcycles));  // divided by  number of cycles
  double E_ExpectationValues = ExpectationValues(0)*norm;
  double E2_ExpectationValues = ExpectationValues(1)*norm;
  double M_ExpectationValues = ExpectationValues(2)*norm;
  double M2_ExpectationValues = ExpectationValues(3)*norm;
  double Mabs_ExpectationValues = ExpectationValues(4)*norm;
  // all expectation values are per spin, divide by 1/NSpins/NSpins
  double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/NSpins/NSpins;
  double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/NSpins/NSpins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << E_ExpectationValues/NSpins/NSpins;
  ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
  ofile << setw(15) << setprecision(8) << M_ExpectationValues/NSpins/NSpins;
  ofile << setw(15) << setprecision(8) << Mvariance/temperature;
  ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues/NSpins/NSpins << endl;
} // end output function

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
  int            j;
  long           k;
  static long    idum2 = 123456789;
  static long    iy=0;
  static long    iv[NTAB];
  double         temp;

  if(*idum <= 0) {
    if(-(*idum) < 1) *idum = 1;
    else             *idum = -(*idum);
    idum2 = (*idum);
    for(j = NTAB + 7; j >= 0; j--) {
      k     = (*idum)/IQ1;
      *idum = IA1*(*idum - k*IQ1) - k*IR1;
      if(*idum < 0) *idum +=  IM1;
      if(j < NTAB)  iv[j]  = *idum;
    }
    iy=iv[0];
  }
  k     = (*idum)/IQ1;
  *idum = IA1*(*idum - k*IQ1) - k*IR1;
  if(*idum < 0) *idum += IM1;
  k     = idum2/IQ2;
  idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
  if(idum2 < 0) idum2 += IM2;
  j     = iy/NDIV;
  iy    = iv[j] - idum2;
  iv[j] = *idum;
  if(iy < 1) iy += IMM1;
  if((temp = AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
