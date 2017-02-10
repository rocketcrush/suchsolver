#include <iostream>
#include <vector>
#include <data.h>
#include <fstream>
#include <omp.h>

using namespace std;


vector<node> n; // node information
vector<element> e; // element information
vector<type1> tp1; // convective fluxes
vector<type2> tp2; // diffusive fluxes
vector<type3> tp3; // generated for unsteady calculations, can be split more efficiently

int schemetype;

void schemeselection()
{
  string line;

  ifstream inputflux("mesh/solution.txt");


  getline(inputflux,line,':');

  inputflux >> schemetype;
}

void datafirst()
{
  n.push_back(node()); // first initialization of all vectors, for i = 0, empty vector

  e.push_back(element());

  tp1.push_back(type1());

  tp2.push_back(type2());

  tp3.push_back(type3());
}

void datasecond()
{
  if(schemetype == 1 || schemetype == 2) // steady,explicit,inviscid
    {
      for (unsigned int i = 1; i < e.size(); i++)
        {
          tp1.push_back(type1());
        }
    }
  else if(schemetype == 3) // unsteady, explicit,inviscid
    {
      for (unsigned int i = 1; i < e.size(); i++)
        {
          tp1.push_back(type1());
          tp3.push_back(type3());
        }
    }
  else if(schemetype == 4 || schemetype == 5) //steady, explicit,viscous
    {
      for (unsigned int i = 1; i < e.size(); i++)
        {
          tp1.push_back(type1());
          tp2.push_back(type2());
        }
    }
  else if(schemetype == 6)  // unsteady, explicit, viscous
    {
      for (unsigned int i = 1; i < e.size(); i++)
        {
          tp1.push_back(type1());
          tp2.push_back(type2());
          tp3.push_back(type3());
        }
    }

}
