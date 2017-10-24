#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <data.h>
#include <constants.h>
#include <initialconditions.h>
#include <mesh.h>
#include <omp.h>

using namespace std;

extern vector<element> e;

double uinf, vinf, roinf, Tinf, muinf, Reinf, pinf, Machinf, hinf, Vinf, ainf;

void readbackup()
{
  if(fileExists("mesh/backup.txt"))
    {

      cout << "Backup found." << endl;
    }
  else
    {
      cout << "Backup not found." << endl;

      cin.get();
    }

  fstream read ("mesh/backup.txt");

  double p, T, u, v;

  int i = 1;

  while(read >> u >> v >> p >> T)
    {
      e[i].u = u / Vinf;
      e[i].v = v / Vinf;
      e[i].p = p / (roinf * pow(Vinf, 2));
      e[i].T = T / Tinf;


      e[i].ro = e[i].p * (roinf * pow(Vinf, 2)) / (R * e[i].T * Tinf * roinf);

      e[i].mu = (1.461e-6 * pow(e[i].T * Tinf, (3. / 2.)) / (e[i].T * Tinf + 110.3)) / muinf;

      e[i].V = sqrt(pow(e[i].u, 2) + pow(e[i].v, 2));
      e[i].et = (cv * e[i].T * Tinf + pow(e[i].V * Vinf, 2) / 2.) / pow(Vinf, 2);

      i++;
    }


}

void initialconditions()
{
  string line;

  fstream read("mesh/initialconditions.txt");

  getline(read,line,':');

  read >> uinf;

  getline(read,line,':');

  read >> vinf;

  getline(read,line,':');

  read >> pinf;

  getline(read,line,':');

  read >> Tinf;

  Vinf = sqrt(pow(uinf, 2) + pow(vinf, 2));

  ainf = sqrt(gama * R * Tinf);

  muinf = 1.461e-6 * pow(Tinf, (3. / 2.)) / (Tinf + 110.3);

  roinf = pinf / (R * Tinf);

  Reinf = roinf * Vinf / muinf;



  Machinf = Vinf / ainf;

  //hinf = (gama * cv * Tinf + 1. / 2. * (pow(uinf, 2) + pow(vinf, 2))) / pow(uinf,2) ;


  if(readbackupflag == 1)
    readbackup();

  else
    {
#pragma omp parallel for schedule(dynamic,1000)
      for(unsigned int i = 1; i < e.size(); i++)
        {
          e[i].u = uinf / Vinf;
          e[i].v = vinf / Vinf;

          e[i].ro = roinf / roinf;
          e[i].T = Tinf / Tinf;
          e[i].mu = muinf / muinf;

          e[i].V = sqrt(pow(e[i].u, 2) + pow(e[i].v, 2));
          e[i].et = (cv * e[i].T * Tinf + pow(e[i].V * Vinf, 2) / 2.) / pow(Vinf, 2);
          e[i].p = pinf / (roinf * pow(Vinf, 2));


          e[i].unew = e[i].u;
          e[i].vnew =  e[i].v;

          e[i].ronew = e[i].ro;
          e[i].Tnew = e[i].T;
          e[i].munew = e[i].mu;

          e[i].Vnew = e[i].V;
          e[i].etnew = e[i].et;
          e[i].pnew = e[i].p;
        }




#pragma omp parallel for schedule(dynamic,1000)
      for(unsigned int i = 1; i < e.size(); i++)
        {
              e[i].u = Vinf / Vinf;
              e[i].v = 0. / Vinf;

              e[i].ro = roinf / roinf;
              e[i].T = Tinf / Tinf;
              e[i].mu = muinf / muinf;

              e[i].V = sqrt(pow(e[i].u, 2) + pow(e[i].v, 2));
              e[i].et = (cv * e[i].T * Tinf + pow(e[i].V * Vinf, 2) / 2.) / pow(Vinf, 2);
              e[i].p = e[i].ro * roinf * R * e[i].T * Tinf / (roinf * pow(Vinf, 2));



              e[i].unew = e[i].u;
              e[i].vnew =  e[i].v;

              e[i].ronew = e[i].ro;
              e[i].Tnew = e[i].T;
              e[i].munew = e[i].mu;

              e[i].Vnew = e[i].V;
              e[i].etnew = e[i].et;
              e[i].pnew = e[i].p;

        }

    }


}

