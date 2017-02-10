#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <data.h>
#include <mesh.h>
#include <constants.h>
#include <omp.h>


using namespace std;

extern vector<node> n;
extern vector<element> e;


extern bool writeneighborflag;
extern bool readneighborflag;

bool fileExists(const std::string& file) {
  struct stat buf;
  return (stat(file.c_str(), &buf) == 0);
}

void readgambitmesh()
{
  double x, y, a, b;


  int ctr = 1;

  int nodesize, elementsize, boundarysize;

  string line;

  int i;

  ifstream input("mesh/mesh.neu");


  if (fileExists("mesh/mesh.neu"))
    cout << "mesh.neu found." << endl;
  else
    {
      cout << "mesh.neu does not exist." << endl;
      cin.get();
    }



  if (input.is_open())
    cout << "mesh.neu opened succesfully." << endl;
  else
    {
      cout << "mesh.neu cannot be opened." << endl;
      cin.get();
    }



  for(int i = 0; i < 6; i++)
    getline(input,line,'\n');


  input >> nodesize >> elementsize >> a >> boundarysize;


  for(int i = 0; i < 3; i++)
    getline(input,line,'\n');


  while (input >> i >> x >> y)
    {
      n.push_back(node());

      n[i].x = x;
      n[i].y = y;
    }

  if(nodesize != n.size() - 1) //checking if the vector size is the same as nodesize
    {
      cout << "nodesize in the beginning of the text file does not equal to vector size." << endl;
      cin.get();
    }


  ifstream inputelement("mesh/mesh.neu");

  for(int i = 0; i <= 10 + nodesize; i++)
    getline(inputelement,line,'\n');


  int n1, n2, n3, n4;


  while (inputelement >> i >> a >> b >> n1 >> n2 >> n3 >> n4)
    {
      e.push_back(element());

      e[i].node[0] = n1;
      e[i].node[1] = n2;
      e[i].node[2] = n3;
      e[i].node[3] = n4;
    }

  if(elementsize != e.size() - 1) //checking if the vector size is the same as nodesize
    {
      cout << "elementsize in the beginning of the text file does not equal to vector size." << endl;
      cin.get();
    }


  //first boundary inflow, case 1

  ifstream inputboundary("mesh/mesh.neu");

  if(ctr == 1)
    {
      while(line != "BOUNDARY CONDITIONS  1.2.1")
        {
          getline(inputboundary,line,'\n');
        }
    }

  getline(inputboundary,line,'\n');

  while(inputboundary >> i >> a >> b)
    {
      if(ctr == 1)
        {
          e[i].side[int(b) - 1].ca = ctr;
          e[i].c = ctr;
        }
    }



  ctr += 1;


  //second boundary outflow, case 2

  ifstream inputboundary2("mesh/mesh.neu");

  if(ctr == 2)
    {
      while(line != "BOUNDARY CONDITIONS  1.2.1")
        {
          getline(inputboundary2,line,'\n');
        }

      getline(inputboundary2,line,'\n');

      while(line != "BOUNDARY CONDITIONS  1.2.1")
        {
          getline(inputboundary2,line,'\n');
        }
    }



  getline(inputboundary2,line,'\n');

  while(inputboundary2 >> i >> a >> b)
    {
      if(ctr == 2)
        {
          e[i].side[int(b) - 1].ca = ctr;
          e[i].c = ctr;
        }
    }

  ctr += 1;


  //third boundary slip wall, case 4

  ifstream inputboundary3("mesh/mesh.neu");

  if(ctr == 3)
    {
      while(line != "BOUNDARY CONDITIONS  1.2.1")
        {
          getline(inputboundary3,line,'\n');
        }

      getline(inputboundary3,line,'\n');

      while(line != "BOUNDARY CONDITIONS  1.2.1")
        {
          getline(inputboundary3,line,'\n');
        }

      getline(inputboundary3,line,'\n');

      while(line != "BOUNDARY CONDITIONS  1.2.1")
        {
          getline(inputboundary3,line,'\n');
        }
    }


  getline(inputboundary3,line,'\n');

  while(inputboundary3 >> i >> a >> b)
    {
      if(ctr == 3)
        {
          e[i].side[int(b) - 1].ca = ctr + 1;
          e[i].c = ctr + 1;
        }
    }

}

void readmesh()
{
  double x, y;

  int nodesize, elementsize;

  string line;

  int i = 1;
  int j = 1;


  ifstream inputnode("mesh/nodes.txt");


  if (fileExists("mesh/nodes.txt"))
    cout << "nodes.txt found." << endl;
  else
    {
      cout << "nodex.txt does not exist." << endl;
      cin.get();
    }



  if (inputnode.is_open())
    cout << "nodes.txt opened succesfully." << endl;
  else
    {
      cout << "nodes.txt cannot be opened." << endl;
      cin.get();
    }



  getline(inputnode,line,':'); // getting the number of nodes

  inputnode >> nodesize;

  getline(inputnode,line,'\n');


  while (inputnode >> x >> y)
    {
      n.push_back(node());

      n[i].x = x;
      n[i].y = y;

      i++;
    }

  if(nodesize != n.size() - 1) //checking if the vector size is the same as nodesize
    {
      cout << "nodesize in the beginning of the text file does not equal to vector size." << endl;
      cin.get();
    }


  ifstream inputelement("mesh/elements.txt");


  if (fileExists("mesh/elements.txt"))
    cout << "elements.txt found." << endl;
  else
    {
      cout << "elements.txt does not exist." << endl;
      cin.get();
    }



  if (inputelement.is_open())
    cout << "elements.txt opened succesfully." << endl;
  else
    {
      cout << "elements.txt cannot be opened." << endl;
      cin.get();
    }

  getline(inputelement,line,':'); // getting the number of elements

  inputelement >> elementsize;

  getline(inputelement,line,'\n');

  int n1, n2, n3, n4;

  i = 1;

  while (inputelement >> n1 >> n2 >> n3 >> n4)
    {
      e.push_back(element());

      e[i].node[0] = n1;
      e[i].node[1] = n2;
      e[i].node[2] = n3;
      e[i].node[3] = n4;

      i++;
    }

  if(elementsize != e.size() - 1) //checking if the vector size is the same as nodesize
    {
      cout << "elementsize in the beginning of the text file does not equal to vector size." << endl;
      cin.get();
    }



  int A, B, C;

  int k;

  for(j = 1; j < 6; j++)
    {
      ifstream inputboundary("mesh/boundary.txt");

      for(k = 1; k <= j; k++)
        getline(inputboundary,line,':');

      inputboundary >> A;

      if (A == 1 || A == 2 || A == 3 || A == 4 || A == 5)
        {
          while (inputboundary >> B >> C)
            {
              e[B].c = A;
              e[B].side[C].ca = A;
            }
        }
    }

}

void normalvectors()
{
  /*        neighbour 2
            |---------|
            |         |
neighbour 3 |         | neighbour 1
            |         |
            |---------|
            neighbour 0           */


#pragma omp parallel for schedule(dynamic,1000) collapse(2)
  for (unsigned int i = 1; i < e.size(); i++)
    {
      for(int d = 0; d < 4; d++)
        {
          e[i].Sx[d] = n[e[i].node[(d + 1) % 4]].y - n[e[i].node[d % 4]].y;

          e[i].Sy[d] = n[e[i].node[d % 4]].x - n[e[i].node[(d + 1) % 4]].x;


          e[i].dS[d] = sqrt(pow(e[i].Sx[d], 2) + pow(e[i].Sy[d], 2));

          e[i].vn[d].x = e[i].Sx[d] / e[i].dS[d];

          e[i].vn[d].y = e[i].Sy[d] / e[i].dS[d];

          e[i].r[d].x = e[e[i].neigh[d]].x - e[i].x;

          e[i].r[d].y = e[e[i].neigh[d]].y - e[i].y;

          e[i].rm[d] = sqrt(pow(e[i].r[d].x, 2) + pow(e[i].r[d].y, 2)); //r magnitude
        }
    }





}

void mesh()
{
  readgambitmesh();

#pragma omp parallel for schedule(dynamic,1000)
  for (unsigned int i = 1; i < e.size(); i++)  //coordinates of the elements
    {
      e[i].x = (n[e[i].node[0]].x + n[e[i].node[1]].x + n[e[i].node[2]].x + n[e[i].node[3]].x) / 4.;
      e[i].y = (n[e[i].node[0]].y + n[e[i].node[1]].y + n[e[i].node[2]].y + n[e[i].node[3]].y) / 4.;

      e[i].area = 1. / 2. * ((n[e[i].node[0]].x - n[e[i].node[2]].x) * (n[e[i].node[1]].y - n[e[i].node[3]].y) + (n[e[i].node[3]].x - n[e[i].node[1]].x) * (n[e[i].node[0]].y - n[e[i].node[2]].y));
    }



#pragma omp parallel for schedule(dynamic,1000)
  for (unsigned int i = 1; i < n.size(); i++) //which elements have the same node?
    {
      int ctr = 0;

      for (unsigned int j = 1; j < e.size(); j++)
        {
          if (e[j].node[0] == i || e[j].node[1] == i || e[j].node[2] == i || e[j].node[3] == i)
            {
              n[i].el[ctr] = j; // node is present in which elements?

              ctr += 1;
            }
        }
    }



  if(readneighborflag)
    {
      ifstream readneighbor("mesh/neighbors.txt");


      if (fileExists("mesh/neighbors.txt"))
        cout << "neighbors.txt found." << endl;
      else
        {
          cout << "neighbors.txt does not exist." << endl;
          cin.get();
        }


      if (readneighbor.is_open())
        cout << "neighbors.txt opened succesfully." << endl;
      else
        {
          cout << "neighbors.txt cannot be opened." << endl;
          cin.get();
        }

      int neigh0, neigh1, neigh2, neigh3;

      int i = 1;

      while (readneighbor >> neigh0 >> neigh1 >> neigh2 >> neigh3)
        {
          e[i].neigh[0] = neigh0;

          e[i].neigh[1] = neigh1;

          e[i].neigh[2] = neigh2;

          e[i].neigh[3] = neigh3;

          i++;
        }
    }

  else
    {
#pragma omp parallel for schedule(dynamic,1000) collapse(2)
      for (unsigned int i = 1; i < e.size(); i++) //neighboring elements
        {
          for (unsigned int j = 1; j < e.size(); j++)
            {
              if (e[i].node[0] == e[j].node[3] && e[i].node[1] == e[j].node[2])
                e[i].neigh[0] = j;

              else if (e[i].node[1] == e[j].node[0] && e[i].node[2] == e[j].node[3])
                e[i].neigh[1] = j;

              else if (e[i].node[2] == e[j].node[1] && e[i].node[3] == e[j].node[0])
                e[i].neigh[2] = j;

              else if (e[i].node[3] == e[j].node[2] && e[i].node[0] == e[j].node[1])
                e[i].neigh[3] = j;
            }
        }
    }

  normalvectors();





  if(writeneighborflag)
    {
      ofstream writeneighbors;

      writeneighbors.open("neighbors.txt");

      int n0, n1, n2, n3;

      for (unsigned int i = 1; i <= e.size(); i++)
        {
          if(e[i].c == 0)
            {
              if(e[e[i].neigh[0]].c != 0)
                n0 = 0;
              else
                n0 = e[i].neigh[0];

              if(e[e[i].neigh[1]].c != 0)
                n1 = 0;
              else
                n1 = e[i].neigh[1];

              if(e[e[i].neigh[2]].c != 0)
                n2 = 0;
              else
                n2 = e[i].neigh[2];

              if(e[e[i].neigh[3]].c != 0)
                n3 = 0;
              else
                n3 = e[i].neigh[3];
            }
          writeneighbors << n0 << " " << n1 << " " << n2 << " " << n3 << endl;
        }

    }


  if(schemetype == 4 || schemetype == 5 || schemetype == 6) // wall boundary is read as slip wall for now, for viscous calculations, case 4 slip wall is shifted to case 5 no slip wall
    {
#pragma omp parallel for schedule(dynamic,1000)
      for (unsigned int i = 1; i < e.size(); i++)
        {
          if(e[i].c == 4)
            {
              e[i].c = 5;

              for(int d = 0; d < 4; d++)
                {
                  if(e[i].neigh[d] == 0)
                    {
                      e[i].side[d].ca = e[i].c;
                    }

                }
            }
        }
    }
}

