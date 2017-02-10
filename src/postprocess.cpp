#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <postprocess.h>
#include <data.h>
#include <constants.h>
#include <initialconditions.h>
#include <mesh.h>
#include <fluxes.h>

using namespace std;

extern vector<node> n;
extern vector<element> e;



void writebackup()
{
  ofstream myfile;


  myfile.open("backup.txt");

  for(int i = 1; i < e.size(); i++)
    myfile << e[i].u * Vinf << "\t" << e[i].v * Vinf << "\t" << e[i].p * roinf * pow(Vinf, 2) << "\t" << e[i].T * Tinf << endl;

  myfile.close();

}

double gradient(int i, string variable)
{
  double grad[2], n[4], var;

  if(variable == "u")
    {
      for(int d = 0; d < 4; d++)
        {
          if(e[i].neigh[d] != 0)
            {
              n[d] = e[e[i].neigh[d]].u;
            }
          else
            {
              n[d] = e[i].u;
            }
        }

      var = e[i].u;
    }
  else if(variable == "v")
    {
      for(int d = 0; d < 4; d++)
        {
          if(e[i].neigh[d] != 0)
            {
              n[d] = e[e[i].neigh[d]].v;
            }
          else
            {
              n[d] = e[i].v;
            }
        }

      var = e[i].v;
    }
  else if(variable == "ro")
    {
      for(int d = 0; d < 4; d++)
        {
          if(e[i].neigh[d] != 0)
            {
              n[d] = e[e[i].neigh[d]].ro;
            }
          else
            {
              n[d] = e[i].ro;
            }
        }

      var = e[i].ro;
    }
  else if(variable == "p")
    {
      for(int d = 0; d < 4; d++)
        {
          if(e[i].neigh[d] != 0)
            {
              n[d] = e[e[i].neigh[d]].p;
            }
          else
            {
              n[d] = e[i].p;
            }
        }

      var = e[i].p;
    }
  else if(variable == "T")
    {
      for(int d = 0; d < 4; d++)
        {
          if(e[i].neigh[d] != 0)
            {
              n[d] = e[e[i].neigh[d]].T;
            }
          else
            {
              n[d] = e[i].T;
            }
        }

      var = e[i].T;
    }
  else if(variable == "et")
    {
      for(int d = 0; d < 4; d++)
        {
          if(e[i].neigh[d] != 0)
            {
              n[d] = e[e[i].neigh[d]].et;
            }
          else
            {
              n[d] = e[i].et;
            }
        }

      var = e[i].et;
    }



  grad[0] = 1. / e[i].area * ((var + n[0]) / 2. * e[i].vn[0].x * e[i].dS[0] + (var + n[1]) / 2. * e[i].vn[1].x * e[i].dS[1] + (var + n[2]) / 2. * e[i].vn[2].x * e[i].dS[2] + (var + n[3]) / 2. * e[i].vn[3].x * e[i].dS[3]);

  grad[1] = 1. / e[i].area * ((var + n[0]) / 2. * e[i].vn[0].y * e[i].dS[0] + (var + n[1]) / 2. * e[i].vn[1].y * e[i].dS[1] + (var + n[2]) / 2. * e[i].vn[2].y * e[i].dS[2] + (var + n[3]) / 2. * e[i].vn[3].y * e[i].dS[3]);



  return sqrt(pow(grad[0], 2) + pow(grad[1], 2));
}

void vtkformatoutput(int it)
{
  ofstream myfile;

  int quadrilateral = 9;

  int numpoints = 4; //number of points needed to present the cell type, 4 points are needed for quadrilateral


  unsigned int i;

  std::ostringstream fileNameStream("output");
  fileNameStream << it << ".vtk";

  string fileName = fileNameStream.str();

  myfile.open(fileName.c_str());


  myfile << "# vtk DataFile Version 2.0 "<< endl << "2D Unstructured" << endl << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " << endl ;

  myfile << "POINTS " << n.size() - 1 << " float" << endl;


  for(i = 1; i < n.size(); i++) //node.x, node.y, node.z (0)
    {
      myfile << n[i].x << " " << n[i].y << " 0" << endl;
    }

  myfile << endl;


  myfile << "CELLS " << e.size() - 1 << " " << (numpoints + 1) * (e.size() - 1) << endl;


  for(i = 1; i < e.size(); i++) // cell type, node0, node1, node2, node3
    {
      myfile << numpoints << " " << e[i].node[0] - 1 << " " << e[i].node[1] - 1 << " " << e[i].node[2] - 1 << " " << e[i].node[3] - 1 << endl;
    }

  myfile << endl;

  myfile << "CELL_TYPES " << e.size() - 1 << endl;

  for(i = 1; i < e.size(); i++) // cell type, node0, node1, node2, node3
    {
      myfile << quadrilateral << endl;
    }

  myfile << endl;

  myfile << "POINT_DATA " << n.size() - 1 << endl << "SCALARS pressure float" << endl << "LOOKUP_TABLE default" << endl;




  double ro, u, v, et, V, T, pres;

  int ctr;


  for (i = 1; i < n.size(); i++)
    {
      double l1a = sqrt(pow((e[n[i].el[0]].x - n[i].x), 2) + pow((e[n[i].el[0]].y - n[i].y), 2));
      double l2a = sqrt(pow((e[n[i].el[1]].x - n[i].x), 2) + pow((e[n[i].el[1]].y - n[i].y), 2));
      double l3a = sqrt(pow((e[n[i].el[2]].x - n[i].x), 2) + pow((e[n[i].el[2]].y - n[i].y), 2));
      double l4a = sqrt(pow((e[n[i].el[3]].x - n[i].x), 2) + pow((e[n[i].el[3]].y - n[i].y), 2));

      if (n[i].el[0] != 0 && n[i].el[1] != 0 && n[i].el[2] != 0 && n[i].el[3] != 0)
        {
          ro = (e[n[i].el[0]].ro / l1a + e[n[i].el[1]].ro / l2a + e[n[i].el[2]].ro / l3a + e[n[i].el[3]].ro / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

          u = (e[n[i].el[0]].u / l1a + e[n[i].el[1]].u / l2a + e[n[i].el[2]].u / l3a + e[n[i].el[3]].u / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

          v = (e[n[i].el[0]].v / l1a + e[n[i].el[1]].v / l2a + e[n[i].el[2]].v / l3a + e[n[i].el[3]].v / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

          et = (e[n[i].el[0]].et / l1a + e[n[i].el[1]].et / l2a + e[n[i].el[2]].et / l3a + e[n[i].el[3]].et / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

          V = sqrt(pow(u, 2) + pow(v, 2));

          T = (et * pow(Vinf, 2) - pow(V * Vinf, 2) / 2.) / (cv * Tinf);

          pres = ro * roinf * R * T * Tinf / (roinf * pow(Vinf, 2));

        }
      else
        {
          ctr = 0;

          for(int j = 0; j < 4; j++)
            {
              if(n[i].el[j] != 0)
                ctr++;
            }

          if (ctr == 1)
            {
              for(int d = 0; d < 4; d++)
                {
                  if(n[i].el[d] != 0)
                    {
                      ro = e[n[i].el[d]].ro;
                      u = e[n[i].el[d]].u;
                      v = e[n[i].el[d]].v;
                      V = e[n[i].el[d]].V;
                      et = e[n[i].el[d]].et;
                      T = e[n[i].el[d]].T;
                      pres = e[n[i].el[d]].p;
                    }
                }
            }
          else
            {
              if(n[i].el[0] != 0 && n[i].el[1] != 0)
                {
                  ro = (e[n[i].el[0]].ro / l1a + e[n[i].el[1]].ro / l2a) / (1. / l1a + 1. / l2a);
                  u = (e[n[i].el[0]].u / l1a + e[n[i].el[1]].u / l2a) / (1. / l1a + 1. / l2a);
                  v = (e[n[i].el[0]].v / l1a + e[n[i].el[1]].v / l2a) / (1. / l1a + 1. / l2a);
                  V = (e[n[i].el[0]].V / l1a + e[n[i].el[1]].V / l2a) / (1. / l1a + 1. / l2a);
                  et = (e[n[i].el[0]].et / l1a + e[n[i].el[1]].et / l2a) / (1. / l1a + 1. / l2a);
                  T = (e[n[i].el[0]].T / l1a + e[n[i].el[1]].T / l2a) / (1. / l1a + 1. / l2a);
                  pres = (e[n[i].el[0]].p / l1a + e[n[i].el[1]].p / l2a) / (1. / l1a + 1. / l2a);
                }

              else if(n[i].el[0] != 0 && n[i].el[2] != 0)
                {
                  ro = (e[n[i].el[0]].ro / l1a + e[n[i].el[2]].ro / l3a) / (1. / l1a + 1. / l3a);
                  u = (e[n[i].el[0]].u / l1a + e[n[i].el[2]].u / l3a) / (1. / l1a + 1. / l3a);
                  v = (e[n[i].el[0]].v / l1a + e[n[i].el[2]].v / l3a) / (1. / l1a + 1. / l3a);
                  V = (e[n[i].el[0]].V / l1a + e[n[i].el[2]].V / l3a) / (1. / l1a + 1. / l3a);
                  et = (e[n[i].el[0]].et / l1a + e[n[i].el[2]].et / l3a) / (1. / l1a + 1. / l3a);
                  T = (e[n[i].el[0]].T / l1a + e[n[i].el[2]].T / l3a) / (1. / l1a + 1. / l3a);
                  pres = (e[n[i].el[0]].p / l1a + e[n[i].el[2]].p / l3a) / (1. / l1a + 1. / l3a);
                }
              else if(n[i].el[0] != 0 && n[i].el[3] != 0)
                {
                  ro = (e[n[i].el[0]].ro / l1a + e[n[i].el[3]].ro / l4a) / (1. / l1a + 1. / l4a);
                  u = (e[n[i].el[0]].u / l1a + e[n[i].el[3]].u / l4a) / (1. / l1a + 1. / l4a);
                  v = (e[n[i].el[0]].v / l1a + e[n[i].el[3]].v / l4a) / (1. / l1a + 1. / l4a);
                  V = (e[n[i].el[0]].V / l1a + e[n[i].el[3]].V / l4a) / (1. / l1a + 1. / l4a);
                  et = (e[n[i].el[0]].et / l1a + e[n[i].el[3]].et / l4a) / (1. / l1a + 1. / l4a);
                  T = (e[n[i].el[0]].T / l1a + e[n[i].el[3]].T / l4a) / (1. / l1a + 1. / l4a);
                  pres = (e[n[i].el[0]].p / l1a + e[n[i].el[3]].p / l4a) / (1. / l1a + 1. / l4a);
                }
              else if(n[i].el[1] != 0 && n[i].el[2] != 0)
                {
                  ro = (e[n[i].el[1]].ro / l2a + e[n[i].el[2]].ro / l3a) / (1. / l2a + 1. / l3a);
                  u = (e[n[i].el[1]].u / l2a + e[n[i].el[2]].u / l3a) / (1. / l2a + 1. / l3a);
                  v = (e[n[i].el[1]].v / l2a + e[n[i].el[2]].v / l3a) / (1. / l2a + 1. / l3a);
                  V = (e[n[i].el[1]].V / l2a + e[n[i].el[2]].V / l3a) / (1. / l2a + 1. / l3a);
                  et = (e[n[i].el[1]].et / l2a + e[n[i].el[2]].et / l3a) / (1. / l2a + 1. / l3a);
                  T = (e[n[i].el[1]].T / l2a + e[n[i].el[2]].T / l3a) / (1. / l2a + 1. / l3a);
                  pres = (e[n[i].el[1]].p / l2a + e[n[i].el[2]].p / l3a) / (1. / l2a + 1. / l3a);
                }
              else if(n[i].el[1] != 0 && n[i].el[3] != 0)
                {
                  ro = (e[n[i].el[1]].ro / l2a + e[n[i].el[3]].ro / l4a) / (1. / l2a + 1. / l4a);
                  u = (e[n[i].el[1]].u / l2a + e[n[i].el[3]].u / l4a) / (1. / l2a + 1. / l4a);
                  v = (e[n[i].el[1]].v / l2a + e[n[i].el[3]].v / l4a) / (1. / l2a + 1. / l4a);
                  V = (e[n[i].el[1]].V / l2a + e[n[i].el[3]].V / l4a) / (1. / l2a + 1. / l4a);
                  et = (e[n[i].el[1]].et / l2a + e[n[i].el[3]].et / l4a) / (1. / l2a + 1. / l4a);
                  T = (e[n[i].el[1]].T / l2a + e[n[i].el[3]].T / l4a) / (1. / l2a + 1. / l4a);
                  pres = (e[n[i].el[1]].p / l2a + e[n[i].el[3]].p / l4a) / (1. / l2a + 1. / l4a);
                }
              else if(n[i].el[3] != 0 && n[i].el[2] != 0)
                {
                  ro = (e[n[i].el[3]].ro / l4a + e[n[i].el[2]].ro / l3a) / (1. / l4a + 1. / l3a);
                  u = (e[n[i].el[3]].u / l4a + e[n[i].el[2]].u / l3a) / (1. / l4a + 1. / l3a);
                  v = (e[n[i].el[3]].v / l4a + e[n[i].el[2]].v / l3a) / (1. / l4a + 1. / l3a);
                  V = (e[n[i].el[3]].V / l4a + e[n[i].el[2]].V / l3a) / (1. / l4a + 1. / l3a);
                  et = (e[n[i].el[3]].et / l4a + e[n[i].el[2]].et / l3a) / (1. / l4a + 1. / l3a);
                  T = (e[n[i].el[3]].T / l4a + e[n[i].el[2]].T / l3a) / (1. / l4a + 1. / l3a);
                  pres = (e[n[i].el[3]].p / l4a + e[n[i].el[2]].p / l3a) / (1. / l4a + 1. / l3a);
                }
            }
        }

      if(dimensionalflag == 1)
        myfile << pres  << endl;
      else
        myfile << pres * roinf * pow(Vinf, 2) << endl;

    }


}


void datformatoutput(int it)
{
  double ro, u, v, et, V, T, pres, grad, gradel[3];

  int ctr;

  unsigned int i;

  ofstream myfile;

  std::ostringstream fileNameStream("output");
  fileNameStream << it << ".dat";

  string fileName = fileNameStream.str();

  myfile.open(fileName.c_str());


  myfile << "VARIABLES = \"X\" , \"Y\" , \"p\" , \"ro\" , \"T\" , \"u\" , \"v\" , \"et\" , \"Mach\" , \"Gradient\""<< endl;
  myfile << "ZONE N = " << n.size() - 1 << ", E = " << e.size() - 1 << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << endl;

  for (i = 1; i < n.size(); i++)
    {
      gradel[0] = gradient(n[i].el[0], "ro");
      gradel[1] = gradient(n[i].el[1], "ro");
      gradel[2] = gradient(n[i].el[2], "ro");
      gradel[3] = gradient(n[i].el[3], "ro");


      double l1a = sqrt(pow((e[n[i].el[0]].x - n[i].x), 2) + pow((e[n[i].el[0]].y - n[i].y), 2));
      double l2a = sqrt(pow((e[n[i].el[1]].x - n[i].x), 2) + pow((e[n[i].el[1]].y - n[i].y), 2));
      double l3a = sqrt(pow((e[n[i].el[2]].x - n[i].x), 2) + pow((e[n[i].el[2]].y - n[i].y), 2));
      double l4a = sqrt(pow((e[n[i].el[3]].x - n[i].x), 2) + pow((e[n[i].el[3]].y - n[i].y), 2));

      if (n[i].el[0] != 0 && n[i].el[1] != 0 && n[i].el[2] != 0 && n[i].el[3] != 0)
        {
          ro = (e[n[i].el[0]].ro / l1a + e[n[i].el[1]].ro / l2a + e[n[i].el[2]].ro / l3a + e[n[i].el[3]].ro / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

          u = (e[n[i].el[0]].u / l1a + e[n[i].el[1]].u / l2a + e[n[i].el[2]].u / l3a + e[n[i].el[3]].u / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

          v = (e[n[i].el[0]].v / l1a + e[n[i].el[1]].v / l2a + e[n[i].el[2]].v / l3a + e[n[i].el[3]].v / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

          et = (e[n[i].el[0]].et / l1a + e[n[i].el[1]].et / l2a + e[n[i].el[2]].et / l3a + e[n[i].el[3]].et / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

          V = sqrt(pow(u, 2) + pow(v, 2));

          T = (et * pow(Vinf, 2) - pow(V * Vinf, 2) / 2.) / (cv * Tinf);

          pres = ro * roinf * R * T * Tinf / (roinf * pow(Vinf, 2));


          grad = (gradel[0] / l1a + gradel[1] / l2a + gradel[2] / l3a + gradel[3] / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

        }
      else
        {
          ctr = 0;

          for(int j = 0; j < 4; j++)
            {
              if(n[i].el[j] != 0)
                ctr++;
            }

          if (ctr == 1)
            {
              for(int d = 0; d < 4; d++)
                {
                  if(n[i].el[d] != 0)
                    {
                      ro = e[n[i].el[d]].ro;
                      u = e[n[i].el[d]].u;
                      v = e[n[i].el[d]].v;
                      V = e[n[i].el[d]].V;
                      et = e[n[i].el[d]].et;
                      T = e[n[i].el[d]].T;
                      pres = e[n[i].el[d]].p;

                      grad = gradel[d];
                    }
                }
            }
          else
            {
              if(n[i].el[0] != 0 && n[i].el[1] != 0)
                {
                  ro = (e[n[i].el[0]].ro / l1a + e[n[i].el[1]].ro / l2a) / (1. / l1a + 1. / l2a);
                  u = (e[n[i].el[0]].u / l1a + e[n[i].el[1]].u / l2a) / (1. / l1a + 1. / l2a);
                  v = (e[n[i].el[0]].v / l1a + e[n[i].el[1]].v / l2a) / (1. / l1a + 1. / l2a);
                  V = (e[n[i].el[0]].V / l1a + e[n[i].el[1]].V / l2a) / (1. / l1a + 1. / l2a);
                  et = (e[n[i].el[0]].et / l1a + e[n[i].el[1]].et / l2a) / (1. / l1a + 1. / l2a);
                  T = (e[n[i].el[0]].T / l1a + e[n[i].el[1]].T / l2a) / (1. / l1a + 1. / l2a);
                  pres = (e[n[i].el[0]].p / l1a + e[n[i].el[1]].p / l2a) / (1. / l1a + 1. / l2a);

                  grad = (gradel[0] / l1a + gradel[1] / l2a) / (1. / l1a + 1. / l2a);
                }

              else if(n[i].el[0] != 0 && n[i].el[2] != 0)
                {
                  ro = (e[n[i].el[0]].ro / l1a + e[n[i].el[2]].ro / l3a) / (1. / l1a + 1. / l3a);
                  u = (e[n[i].el[0]].u / l1a + e[n[i].el[2]].u / l3a) / (1. / l1a + 1. / l3a);
                  v = (e[n[i].el[0]].v / l1a + e[n[i].el[2]].v / l3a) / (1. / l1a + 1. / l3a);
                  V = (e[n[i].el[0]].V / l1a + e[n[i].el[2]].V / l3a) / (1. / l1a + 1. / l3a);
                  et = (e[n[i].el[0]].et / l1a + e[n[i].el[2]].et / l3a) / (1. / l1a + 1. / l3a);
                  T = (e[n[i].el[0]].T / l1a + e[n[i].el[2]].T / l3a) / (1. / l1a + 1. / l3a);
                  pres = (e[n[i].el[0]].p / l1a + e[n[i].el[2]].p / l3a) / (1. / l1a + 1. / l3a);

                  grad = (gradel[0] / l1a + gradel[2] / l3a) / (1. / l1a + 1. / l3a);
                }
              else if(n[i].el[0] != 0 && n[i].el[3] != 0)
                {
                  ro = (e[n[i].el[0]].ro / l1a + e[n[i].el[3]].ro / l4a) / (1. / l1a + 1. / l4a);
                  u = (e[n[i].el[0]].u / l1a + e[n[i].el[3]].u / l4a) / (1. / l1a + 1. / l4a);
                  v = (e[n[i].el[0]].v / l1a + e[n[i].el[3]].v / l4a) / (1. / l1a + 1. / l4a);
                  V = (e[n[i].el[0]].V / l1a + e[n[i].el[3]].V / l4a) / (1. / l1a + 1. / l4a);
                  et = (e[n[i].el[0]].et / l1a + e[n[i].el[3]].et / l4a) / (1. / l1a + 1. / l4a);
                  T = (e[n[i].el[0]].T / l1a + e[n[i].el[3]].T / l4a) / (1. / l1a + 1. / l4a);
                  pres = (e[n[i].el[0]].p / l1a + e[n[i].el[3]].p / l4a) / (1. / l1a + 1. / l4a);

                  grad = (gradel[0] / l1a + gradel[3] / l4a) / (1. / l1a + 1. / l4a);
                }
              else if(n[i].el[1] != 0 && n[i].el[2] != 0)
                {
                  ro = (e[n[i].el[1]].ro / l2a + e[n[i].el[2]].ro / l3a) / (1. / l2a + 1. / l3a);
                  u = (e[n[i].el[1]].u / l2a + e[n[i].el[2]].u / l3a) / (1. / l2a + 1. / l3a);
                  v = (e[n[i].el[1]].v / l2a + e[n[i].el[2]].v / l3a) / (1. / l2a + 1. / l3a);
                  V = (e[n[i].el[1]].V / l2a + e[n[i].el[2]].V / l3a) / (1. / l2a + 1. / l3a);
                  et = (e[n[i].el[1]].et / l2a + e[n[i].el[2]].et / l3a) / (1. / l2a + 1. / l3a);
                  T = (e[n[i].el[1]].T / l2a + e[n[i].el[2]].T / l3a) / (1. / l2a + 1. / l3a);
                  pres = (e[n[i].el[1]].p / l2a + e[n[i].el[2]].p / l3a) / (1. / l2a + 1. / l3a);

                  grad = (gradel[1] / l2a + gradel[2] / l3a) / (1. / l2a + 1. / l3a);
                }
              else if(n[i].el[1] != 0 && n[i].el[3] != 0)
                {
                  ro = (e[n[i].el[1]].ro / l2a + e[n[i].el[3]].ro / l4a) / (1. / l2a + 1. / l4a);
                  u = (e[n[i].el[1]].u / l2a + e[n[i].el[3]].u / l4a) / (1. / l2a + 1. / l4a);
                  v = (e[n[i].el[1]].v / l2a + e[n[i].el[3]].v / l4a) / (1. / l2a + 1. / l4a);
                  V = (e[n[i].el[1]].V / l2a + e[n[i].el[3]].V / l4a) / (1. / l2a + 1. / l4a);
                  et = (e[n[i].el[1]].et / l2a + e[n[i].el[3]].et / l4a) / (1. / l2a + 1. / l4a);
                  T = (e[n[i].el[1]].T / l2a + e[n[i].el[3]].T / l4a) / (1. / l2a + 1. / l4a);
                  pres = (e[n[i].el[1]].p / l2a + e[n[i].el[3]].p / l4a) / (1. / l2a + 1. / l4a);

                  grad = (gradel[1] / l2a + gradel[2] / l4a) / (1. / l1a + 1. / l4a);
                }
              else if(n[i].el[3] != 0 && n[i].el[2] != 0)
                {
                  ro = (e[n[i].el[3]].ro / l4a + e[n[i].el[2]].ro / l3a) / (1. / l4a + 1. / l3a);
                  u = (e[n[i].el[3]].u / l4a + e[n[i].el[2]].u / l3a) / (1. / l4a + 1. / l3a);
                  v = (e[n[i].el[3]].v / l4a + e[n[i].el[2]].v / l3a) / (1. / l4a + 1. / l3a);
                  V = (e[n[i].el[3]].V / l4a + e[n[i].el[2]].V / l3a) / (1. / l4a + 1. / l3a);
                  et = (e[n[i].el[3]].et / l4a + e[n[i].el[2]].et / l3a) / (1. / l4a + 1. / l3a);
                  T = (e[n[i].el[3]].T / l4a + e[n[i].el[2]].T / l3a) / (1. / l4a + 1. / l3a);
                  pres = (e[n[i].el[3]].p / l4a + e[n[i].el[2]].p / l3a) / (1. / l4a + 1. / l3a);

                  grad = (gradel[3] / l4a + gradel[2] / l3a) / (1. / l4a + 1. / l3a);
                }
            }
        }

      if(dimensionalflag == 1)
        myfile << n[i].x << "\t" << n[i].y << "\t" << pres << "\t" << ro << "\t" << T << "\t" << u << "\t" << v << "\t" << et << "\t" << V * Vinf / sqrt(gama * R * T * Tinf) << "\t" << grad << "\t" <<endl;
      else
        myfile << n[i].x << "\t" << n[i].y << "\t" << pres * roinf * pow(Vinf, 2) << "\t" << ro * roinf << "\t" << T * Tinf << "\t" << u * Vinf << "\t" << v * Vinf << "\t" << et * pow(Vinf, 2) << "\t" << V * Vinf / sqrt(gama * R * T * Tinf) << "\t" << grad << "\t" <<endl;

    }

  for (i = 1; i < e.size(); i++)
    {
        myfile << e[i].node[0] << " " << e[i].node[1] << " " << e[i].node[2] << " " << e[i].node[3] << endl;
    }


  myfile.close();

}

void dattransientoutput(int it)
{
  double ro, u, v, et, V, T, pres, grad, gradel[3];

  int ctr;

  unsigned int i;

  ofstream myfile("output.dat", std::ios_base::app | std::ios_base::out);

  myfile << "ZONE N = " << n.size() - 1 << ", E = " << e.size() - 1 << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL, T = \"" << it * physicaldeltat << " seconds \"" <<endl;

  for (i = 1; i < n.size(); i++)
    {
      gradel[0] = gradient(n[i].el[0], "ro");
      gradel[1] = gradient(n[i].el[1], "ro");
      gradel[2] = gradient(n[i].el[2], "ro");
      gradel[3] = gradient(n[i].el[3], "ro");


      double l1a = sqrt(pow((e[n[i].el[0]].x - n[i].x), 2) + pow((e[n[i].el[0]].y - n[i].y), 2));
      double l2a = sqrt(pow((e[n[i].el[1]].x - n[i].x), 2) + pow((e[n[i].el[1]].y - n[i].y), 2));
      double l3a = sqrt(pow((e[n[i].el[2]].x - n[i].x), 2) + pow((e[n[i].el[2]].y - n[i].y), 2));
      double l4a = sqrt(pow((e[n[i].el[3]].x - n[i].x), 2) + pow((e[n[i].el[3]].y - n[i].y), 2));

      if (n[i].el[0] != 0 && n[i].el[1] != 0 && n[i].el[2] != 0 && n[i].el[3] != 0)
        {
          ro = (e[n[i].el[0]].ro / l1a + e[n[i].el[1]].ro / l2a + e[n[i].el[2]].ro / l3a + e[n[i].el[3]].ro / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

          u = (e[n[i].el[0]].u / l1a + e[n[i].el[1]].u / l2a + e[n[i].el[2]].u / l3a + e[n[i].el[3]].u / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

          v = (e[n[i].el[0]].v / l1a + e[n[i].el[1]].v / l2a + e[n[i].el[2]].v / l3a + e[n[i].el[3]].v / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

          et = (e[n[i].el[0]].et / l1a + e[n[i].el[1]].et / l2a + e[n[i].el[2]].et / l3a + e[n[i].el[3]].et / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

          V = sqrt(pow(u, 2) + pow(v, 2));

          T = (et * pow(Vinf, 2) - pow(V * Vinf, 2) / 2.) / (cv * Tinf);

          pres = ro * roinf * R * T * Tinf / (roinf * pow(Vinf, 2));


          grad = (gradel[0] / l1a + gradel[1] / l2a + gradel[2] / l3a + gradel[3] / l4a) / (1. / l1a + 1. / l2a + 1. / l3a + 1. / l4a);

        }
      else
        {
          ctr = 0;

          for(int j = 0; j < 4; j++)
            {
              if(n[i].el[j] != 0)
                ctr++;
            }

          if (ctr == 1)
            {
              for(int d = 0; d < 4; d++)
                {
                  if(n[i].el[d] != 0)
                    {
                      ro = e[n[i].el[d]].ro;
                      u = e[n[i].el[d]].u;
                      v = e[n[i].el[d]].v;
                      V = e[n[i].el[d]].V;
                      et = e[n[i].el[d]].et;
                      T = e[n[i].el[d]].T;
                      pres = e[n[i].el[d]].p;

                      grad = gradel[d];
                    }
                }
            }
          else
            {
              if(n[i].el[0] != 0 && n[i].el[1] != 0)
                {
                  ro = (e[n[i].el[0]].ro / l1a + e[n[i].el[1]].ro / l2a) / (1. / l1a + 1. / l2a);
                  u = (e[n[i].el[0]].u / l1a + e[n[i].el[1]].u / l2a) / (1. / l1a + 1. / l2a);
                  v = (e[n[i].el[0]].v / l1a + e[n[i].el[1]].v / l2a) / (1. / l1a + 1. / l2a);
                  V = (e[n[i].el[0]].V / l1a + e[n[i].el[1]].V / l2a) / (1. / l1a + 1. / l2a);
                  et = (e[n[i].el[0]].et / l1a + e[n[i].el[1]].et / l2a) / (1. / l1a + 1. / l2a);
                  T = (e[n[i].el[0]].T / l1a + e[n[i].el[1]].T / l2a) / (1. / l1a + 1. / l2a);
                  pres = (e[n[i].el[0]].p / l1a + e[n[i].el[1]].p / l2a) / (1. / l1a + 1. / l2a);

                  grad = (gradel[0] / l1a + gradel[1] / l2a) / (1. / l1a + 1. / l2a);
                }

              else if(n[i].el[0] != 0 && n[i].el[2] != 0)
                {
                  ro = (e[n[i].el[0]].ro / l1a + e[n[i].el[2]].ro / l3a) / (1. / l1a + 1. / l3a);
                  u = (e[n[i].el[0]].u / l1a + e[n[i].el[2]].u / l3a) / (1. / l1a + 1. / l3a);
                  v = (e[n[i].el[0]].v / l1a + e[n[i].el[2]].v / l3a) / (1. / l1a + 1. / l3a);
                  V = (e[n[i].el[0]].V / l1a + e[n[i].el[2]].V / l3a) / (1. / l1a + 1. / l3a);
                  et = (e[n[i].el[0]].et / l1a + e[n[i].el[2]].et / l3a) / (1. / l1a + 1. / l3a);
                  T = (e[n[i].el[0]].T / l1a + e[n[i].el[2]].T / l3a) / (1. / l1a + 1. / l3a);
                  pres = (e[n[i].el[0]].p / l1a + e[n[i].el[2]].p / l3a) / (1. / l1a + 1. / l3a);

                  grad = (gradel[0] / l1a + gradel[2] / l3a) / (1. / l1a + 1. / l3a);
                }
              else if(n[i].el[0] != 0 && n[i].el[3] != 0)
                {
                  ro = (e[n[i].el[0]].ro / l1a + e[n[i].el[3]].ro / l4a) / (1. / l1a + 1. / l4a);
                  u = (e[n[i].el[0]].u / l1a + e[n[i].el[3]].u / l4a) / (1. / l1a + 1. / l4a);
                  v = (e[n[i].el[0]].v / l1a + e[n[i].el[3]].v / l4a) / (1. / l1a + 1. / l4a);
                  V = (e[n[i].el[0]].V / l1a + e[n[i].el[3]].V / l4a) / (1. / l1a + 1. / l4a);
                  et = (e[n[i].el[0]].et / l1a + e[n[i].el[3]].et / l4a) / (1. / l1a + 1. / l4a);
                  T = (e[n[i].el[0]].T / l1a + e[n[i].el[3]].T / l4a) / (1. / l1a + 1. / l4a);
                  pres = (e[n[i].el[0]].p / l1a + e[n[i].el[3]].p / l4a) / (1. / l1a + 1. / l4a);

                  grad = (gradel[0] / l1a + gradel[3] / l4a) / (1. / l1a + 1. / l4a);
                }
              else if(n[i].el[1] != 0 && n[i].el[2] != 0)
                {
                  ro = (e[n[i].el[1]].ro / l2a + e[n[i].el[2]].ro / l3a) / (1. / l2a + 1. / l3a);
                  u = (e[n[i].el[1]].u / l2a + e[n[i].el[2]].u / l3a) / (1. / l2a + 1. / l3a);
                  v = (e[n[i].el[1]].v / l2a + e[n[i].el[2]].v / l3a) / (1. / l2a + 1. / l3a);
                  V = (e[n[i].el[1]].V / l2a + e[n[i].el[2]].V / l3a) / (1. / l2a + 1. / l3a);
                  et = (e[n[i].el[1]].et / l2a + e[n[i].el[2]].et / l3a) / (1. / l2a + 1. / l3a);
                  T = (e[n[i].el[1]].T / l2a + e[n[i].el[2]].T / l3a) / (1. / l2a + 1. / l3a);
                  pres = (e[n[i].el[1]].p / l2a + e[n[i].el[2]].p / l3a) / (1. / l2a + 1. / l3a);

                  grad = (gradel[1] / l2a + gradel[2] / l3a) / (1. / l2a + 1. / l3a);
                }
              else if(n[i].el[1] != 0 && n[i].el[3] != 0)
                {
                  ro = (e[n[i].el[1]].ro / l2a + e[n[i].el[3]].ro / l4a) / (1. / l2a + 1. / l4a);
                  u = (e[n[i].el[1]].u / l2a + e[n[i].el[3]].u / l4a) / (1. / l2a + 1. / l4a);
                  v = (e[n[i].el[1]].v / l2a + e[n[i].el[3]].v / l4a) / (1. / l2a + 1. / l4a);
                  V = (e[n[i].el[1]].V / l2a + e[n[i].el[3]].V / l4a) / (1. / l2a + 1. / l4a);
                  et = (e[n[i].el[1]].et / l2a + e[n[i].el[3]].et / l4a) / (1. / l2a + 1. / l4a);
                  T = (e[n[i].el[1]].T / l2a + e[n[i].el[3]].T / l4a) / (1. / l2a + 1. / l4a);
                  pres = (e[n[i].el[1]].p / l2a + e[n[i].el[3]].p / l4a) / (1. / l2a + 1. / l4a);

                  grad = (gradel[1] / l2a + gradel[2] / l4a) / (1. / l1a + 1. / l4a);
                }
              else if(n[i].el[3] != 0 && n[i].el[2] != 0)
                {
                  ro = (e[n[i].el[3]].ro / l4a + e[n[i].el[2]].ro / l3a) / (1. / l4a + 1. / l3a);
                  u = (e[n[i].el[3]].u / l4a + e[n[i].el[2]].u / l3a) / (1. / l4a + 1. / l3a);
                  v = (e[n[i].el[3]].v / l4a + e[n[i].el[2]].v / l3a) / (1. / l4a + 1. / l3a);
                  V = (e[n[i].el[3]].V / l4a + e[n[i].el[2]].V / l3a) / (1. / l4a + 1. / l3a);
                  et = (e[n[i].el[3]].et / l4a + e[n[i].el[2]].et / l3a) / (1. / l4a + 1. / l3a);
                  T = (e[n[i].el[3]].T / l4a + e[n[i].el[2]].T / l3a) / (1. / l4a + 1. / l3a);
                  pres = (e[n[i].el[3]].p / l4a + e[n[i].el[2]].p / l3a) / (1. / l4a + 1. / l3a);

                  grad = (gradel[3] / l4a + gradel[2] / l3a) / (1. / l4a + 1. / l3a);
                }
            }
        }

      if(dimensionalflag == 1)
        myfile << n[i].x << "\t" << n[i].y << "\t" << pres << "\t" << ro << "\t" << T << "\t" << u << "\t" << v << "\t" << et << "\t" << V * Vinf / sqrt(gama * R * T * Tinf) << "\t" << grad << "\t" <<endl;
      else
        myfile << n[i].x << "\t" << n[i].y << "\t" << pres * roinf * pow(Vinf, 2) << "\t" << ro * roinf << "\t" << T * Tinf << "\t" << u * Vinf << "\t" << v * Vinf << "\t" << et * pow(Vinf, 2) << "\t" << V * Vinf / sqrt(gama * R * T * Tinf) << "\t" << grad << "\t" <<endl;

    }

  for (i = 1; i < e.size(); i++)
    {
        myfile << e[i].node[0] << " " << e[i].node[1] << " " << e[i].node[2] << " " << e[i].node[3] << endl;
    }


  myfile.close();
}

void output(int it)
{
  if(schemetype == 1 || schemetype == 2 || schemetype == 4 || schemetype == 5)
    {
      if(vtkoutputflag)
        vtkformatoutput(it);

      if(datoutputflag)
        datformatoutput(it);
    }
  else if(schemetype == 3 || schemetype == 6)
    {
      dattransientoutput(it);
    }

  ofstream myfile;

  unsigned int i;


  switch(wallflag)
    {
    case 2:
      {
        myfile.open("wall-density.txt");

        if(dimensionalflag == 0)
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].ro * roinf << endl;

              }
          }
        else
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].ro << endl;

              }

          }
      }

    case 3:
      {
        myfile.open("wall-uvelocity.txt");

        if(dimensionalflag == 0)
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].u * Vinf << endl;

              }
          }
        else
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].u << endl;

              }

          }
      }

    case 4:
      {
        myfile.open("wall-velocity.txt");

        if(dimensionalflag == 0)
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].v * Vinf << endl;

              }
          }
        else
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].v << endl;

              }

          }
      }

    case 5:
      {
        myfile.open("wall-vvelocity.txt");

        if(dimensionalflag == 0)
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].v * Vinf << endl;

              }
          }
        else
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].v << endl;

              }

          }
      }

    case 6:
      {
        myfile.open("wall-totalenergy.txt");

        if(dimensionalflag == 0)
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].et * pow(Vinf,2) << endl;

              }
          }
        else
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].et << endl;

              }

          }
      }

    case 7:
      {
        myfile.open("wall-temperature.txt");

        if(dimensionalflag == 0)
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].T * Tinf << endl;

              }
          }
        else
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].T << endl;

              }

          }
      }

    case 8:
      {
        myfile.open("wall-pressure.txt");

        if(dimensionalflag == 0)
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].p * roinf * pow(Vinf, 2) << endl;

              }
          }
        else
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].p << endl;

              }

          }
      }

    case 9:
      {
        myfile.open("wall-viscosity.txt");

        if(dimensionalflag == 0)
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].mu * muinf << endl;

              }
          }
        else
          {
            for(i = 1; i < e.size(); i++)
              {
                if(e[i].c == 5 || e[i].c == 4)
                  myfile << e[i].x << "\t" << e[i].y << "\t" << e[i].mu << endl;

              }

          }
      }

    case 10:
      {
       /* myfile.open("wall-heattransferrate.txt");

        int a;

        for(i = 1; i < e.size(); i++)
          {
            if(e[i].c == 5 || e[i].c == 4)
              {

                for(int d = 0; d < 4; d++)
                  {
                    if(e[i].neigh[d] == 0 && e[e[i].neigh[(d + 2) % 4]].c == 0)
                      a = (d + 2) % 4;
                  }


                myfile << e[i].x << "\t" << e[i].y << "\t" << gama * R / (gama - 1.) * e[i].mu * muinf / Pr * (dTdx(i) * e[i].vn[a].x + dTdy(i) * e[i].vn[a].y) * Tinf << endl;
              }
          }*/
      }


      myfile.close();

    }


}


