extern int schemetype;

void schemeselection();

void datafirst();

void datasecond();



struct node
{
  double x, y;

  int el[4];
};

struct element
{
  double x, y;

  double area;
  
  struct normalvector
  {
    double x;
    
    double y;
    
  }vn[4]; //vector n, only the direction of vector S
  
  double Sx[4], Sy[4];
  
  double dS[4];
  
  int node[4], neigh[4];
  
  double ro, u, v, et, T, p, V;
  
  double ronew, unew, vnew, etnew, Tnew, pnew, Vnew;
  
  double mu;
  
  double munew;
  
  double pprev; //for checking convergence, previous pressure
  double roprev; //for checking convergence,
  double uprev; //for checking convergence,
  double vprev; //for checking convergence, 
  double etprev; //for checking convergence,
  double Tprev;
  double muprev;
  double Vprev;

  struct side
  {
    int ca;

  }side[4];

  int c; //case of the boundary
  
  
  struct rvector
  {
    double x;
    double y;
    
  }r[4];
  
  double rm[4]; // magnitude of r vector
};



struct type1
{
  double W[4];
  
  double Fc[4][4];
};

struct type2
{
  double Fv[4][4];
};

struct type3
{
  double W[4][3];
  
  double Q[4];
};
  
  
struct unsteady1
{
  double ro[3], u[3], v[3], et[3], T[3], p[3], V[3];
  
  double mu[3];
  
  
  double W[4][3];
  
  double Fc[4][4][3]; // location, matrix, time
  
  double Q[4];
  
};
  
  
  
  
  
  