#include <iostream>
#include <fstream>
#include <constants.h>

using namespace std;

double R;
double gama;
double Pr;
double cv;
double cp;
double C;
double physicaldeltat;

bool readbackupflag;
int writebackupflag;
bool dimensionalflag;
int wallflag;
bool writeneighborflag;
bool readneighborflag;
bool datoutputflag;
bool vtkoutputflag;
bool wallheatrateflag;

void getconstants()
{
  ifstream read("mesh/constants.txt");


  string trashline;

  getline(read,trashline,':');

  read >> R;

  getline(read,trashline,':');

  read >> gama;

  getline(read,trashline,':');

  read >> Pr;

  cv = R / (gama - 1);

  cp = gama * cv;

  getline(read,trashline,':');

  read >> C;

  getline(read,trashline,':');

  read >> physicaldeltat;




  ifstream optional("mesh/optional.txt");


  getline(optional,trashline,':');

  optional >> readbackupflag;


  getline(optional,trashline,':');

  optional >> wallflag;


  getline(optional,trashline,':');

  optional >> dimensionalflag;


  getline(optional,trashline,':');

  optional >> writebackupflag;



  getline(optional,trashline,':');

  optional >> writeneighborflag;


  getline(optional,trashline,':');

  optional >> readneighborflag;


  getline(optional,trashline,':');

  optional >> datoutputflag;


  getline(optional,trashline,':');

  optional >> vtkoutputflag;
}
