#include <iostream>
#include <vector>
#include <constants.h>
#include <initialconditions.h>
#include <calculations.h>
#include <data.h>
#include <mesh.h>
#include <omp.h>
#include <postprocess.h>

using namespace std;


int main()
{
  //omp_set_dynamic(1);

  omp_set_nested(1);


  omp_set_num_threads(5);

  //cout << omp_get_max_threads();

  cout << "Hello World!" << endl;


  schemeselection();

  cout << "schemeselection done" << endl;


  datafirst();

  cout << "datafirst done" << endl;


  getconstants();

  cout << "getconstants done" << endl;


  mesh();

  cout << "mesh done" << endl;


  datasecond();

  cout << "data second done" << endl;


  initialconditions();

  cout << "initialconditions done" << endl;



  struct timespec start, finish;
  double elapsed;

  clock_gettime(CLOCK_MONOTONIC, &start);


  calculations(10000);


  clock_gettime(CLOCK_MONOTONIC, &finish);


  elapsed = (finish.tv_sec - start.tv_sec);
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

  cout << "Running Time = " << elapsed << " seconds." << endl;

  cout << "with the CFL number of " << C;

  return 0;
}
