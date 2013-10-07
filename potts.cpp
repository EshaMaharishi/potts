#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <utility>
#include <set>
#include <map>

/*
    Note that if you IMPORT something, you best be
    sure that N, numCells, and OFFSETS are correct
*/

#define  IMPORT           0	// 0 or 1
#define  COLLAGEN_OFFSET  1	// 0 or 1
#define  PERIMETER_OFFSET 0	// 0 or 1

/*******************************************************************************/
/*** Global Parameters ***/

  /*** Random Number Generation ***/

  int seed = 			0;

  /*** Simulation Length ( = Loops * Flips ) ***/

  const int numLoops =		10*10;
  const int numFlips =		pow(2,0);
  const int numPrint = 		numLoops/100;

  const int chunkSize =		2;

  /*** Energies ***/

  const double beta =		1.0;

  bool doPrinting = 		true;

  const double J_air =		0.0;
  const double J_cel =		0.0;

  const double J_col =		-50.00000000;

  const double L_vol = 		.05; // heavy penalty for large volume
  /*const*/ double L_ani =	0.0; // negative penalty for large perimeter
  /*const*/ double L_blb =	2.0; // heavy penalty for wiggliness

  /*** System ***/

  const int N =			100;

  const int numCells =		1;
  /*const*/ int numCollagen =	0;

  const double cellSpawn = 	10.0;
  const double cellRadius =	30.0;

  const int collagenWidth =	1;

/*******************************************************************************/
/*** Global Variables ***/

  int lattice[N][N][2] = {0};
  std::map< int, std::set< std::pair<int, int> > > cellVolumeList;
  std::map< int, std::set< std::pair<int, int> > > cellPerimeterList;

  double totalEnergy;
  double avg[3],dev[3];

  const double targetVolume = 3.141593*cellRadius*cellRadius;

/*******************************************************************************/
/*** Functions ***/

  #include "potts_print_.h"
  #include "potts_spawn_.h"
  #include "potts_energy_.h"
  #include "potts_flip_.h"
  #include "potts_analysis_.h"

/*******************************************************************************/
/*** Main ***/

//optional arguments: output directory name, numCollagen, doprinting
int main(int argc, char *argv[])
{
  printf("Running...\n");

  if( argc > 2 )
	  L_blb = atoi(argv[2]);

  if( argc > 3 )
	  doPrinting = false;

  char   dname[100];
  char   fname[100];

  /* Seed random number generator */

  if(seed==0)
    seed=time(0);
  srand(seed);

  /* Create output directory */

  system("rm -rf output");

  if(argv[1]==NULL)
    strcpy(dname,"output");
  else
    strcpy(dname,argv[1]);
  strcpy(fname,"mkdir ");
  strcat(fname,dname);
  strcat(fname," 2>/dev/null");
  system(fname);

  /* Print log file */

  strcpy(fname,dname);
  strcat(fname,"/aaa_log_.txt");
  printLog(fname);

  /* Let there be life */

  #if IMPORT
  if(doPrinting){printf("\n  Warning: importing cancer.\n\n");}
  readCells();
  readCollagen();
  #else
  if(doPrinting){printf("\n  Warning: creating cancer.\n\n");}
  putCells();
  putCollagen();
  #endif

  strcpy(fname,dname);
  strcat(fname,"/lattice_0_.txt");
  printLattice(fname);
  strcpy(fname,dname);
  strcat(fname,"/collagen_.txt");
  printCollagen(fname);

  /* Calculate the initial energy */

  measureCells();

  totalEnergy=Hamiltonian();

  strcpy(fname,dname);
  strcat(fname,"/aaa_energy_.txt");
  FILE* efile=fopen(fname,"w");
  if(doPrinting){
  fprintf(efile,"%5d ",0);
  fprintf(efile,"%8.3lf %8.3lf ",totalEnergy,0.0);
  fprintf(efile,"%8.3lf %8.3lf ",avg[0],dev[0]);
  fprintf(efile,"%8.3lf %8.3lf ",avg[1],dev[1]);
  fprintf(efile,"%8.3lf %8.3lf ",avg[2],dev[2]);
  fprintf(efile,"\n");
  fflush(efile);
  }

  /* Perform flips and conditionally accept the change via the Metropolis algorithm */

  int accepted;

  if(doPrinting){printf("  Mutating: %d x 2^%d spin flips...\n\n",numLoops,(int)log2((double)numFlips));}

  for(int outerCount=1; outerCount<=numLoops; outerCount++){

    printf("\r    count = %d",outerCount);
    fflush(stdout);

    accepted=0;
    for(int count=0; count<numFlips; count++){
      accepted+=flip();
    }
    //measureCells();

    if( doPrinting ){
    fprintf(efile,"%5d ",outerCount);
    fprintf(efile,"%8.3lf %8.3lf ",totalEnergy,(double)accepted/(double)numFlips);
    fprintf(efile,"%8.3lf %8.3lf ",avg[0],dev[0]);
    fprintf(efile,"%8.3lf %8.3lf ",avg[1],dev[1]);
    fprintf(efile,"%8.3lf %8.3lf ",avg[2],dev[2]);
    fprintf(efile,"\n");
    fflush(efile);

 //   if(outerCount%numPrint==0){
      sprintf(fname,"%s/lattice_%d_.txt",dname,outerCount/*/numPrint*/);
      printLattice(fname);
  //  }
    }

  }

  //printf("\r    finished!              \n\n");

  fclose(efile);

  /*** Show some stats ***/
  if( doPrinting ){
  printf("  Final statistics...\n\n");
  printf("    acceptance ratio : %8.3lf\n\n",(double)accepted/(double)(numFlips));
  printf("    cell volume      : %8.3lf +/- %7.3lf\n",avg[0],dev[0]);
  printf("    cell perimeter   : %8.3lf +/- %7.3lf\n",avg[1],dev[1]);
  printf("    cell anisotropy  : %8.3lf +/- %7.3lf\n\n",avg[2],dev[2]);
  printf("  Done.\n\n");
  }

  if( argc < 1 )
	  printf("\n\n Lblb: %f \t\t Anisotropy: %f\n", L_blb, avg[2]);
  else
	  printf("\n\n %f \t\t %f\n", L_blb, avg[2]);
  fflush(stdout);

  return 0;

}

















































