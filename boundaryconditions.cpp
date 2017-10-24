#include <iostream>
#include <vector>
#include <cmath>
#include <data.h>
#include <constants.h>
#include <initialconditions.h>
#include <mesh.h>
#include <fluxes.h>
#include <omp.h>

using namespace std;

extern vector<node> n;
extern vector<element> e;
extern vector<type1> tp1;
extern vector<type2> tp2;
extern vector<type3> tp3;



double extrapolate(int i, string variable)
{


  //first order extrapolation from the first interior element

  int a;

  for(int d = 0; d < 4; d++)
    {
      if(e[i].neigh[d] == 0 && e[e[i].neigh[(d + 2) % 4]].c == 0)
        a = (d + 2) % 4;
    }

  if(variable == "u")
    {
      return e[e[i].neigh[a]].u;
    }
  else if(variable == "v")
    {
      return e[e[i].neigh[a]].v;
    }
  else if(variable == "ro")
    {
      return e[e[i].neigh[a]].ro;
    }
  else if(variable == "p")
    {
      return e[e[i].neigh[a]].p;
    }
  else if(variable == "T")
    {
      return e[e[i].neigh[a]].T;
    }
  else if(variable == "et")
    {
      return e[e[i].neigh[a]].et;
    }

}

void boundarycondition(int i)
{

  for(int d = 0; d < 4; d++)
    {
      if (e[i].side[d].ca == 1)
        {
          tp1[i].Fc[d][0] = roinf / roinf * (uinf / Vinf * e[i].vn[d].x + vinf / Vinf * e[i].vn[d].y);

          tp1[i].Fc[d][1] = roinf / roinf * (uinf / Vinf) * (uinf / Vinf * e[i].vn[d].x + vinf / Vinf * e[i].vn[d].y) + e[i].vn[d].x * pinf / (roinf * pow(Vinf, 2));

          tp1[i].Fc[d][2] = roinf / roinf * (vinf / Vinf) * (uinf / Vinf * e[i].vn[d].x + vinf / Vinf * e[i].vn[d].y) + e[i].vn[d].y * pinf / (roinf * pow(Vinf, 2));

          tp1[i].Fc[d][3] = roinf / roinf * (((cv * Tinf + pow(Vinf, 2) / 2.) / pow(Vinf, 2)) + pinf / (roinf * pow(Vinf, 2))) * (uinf / Vinf * e[i].vn[d].x + vinf / Vinf * e[i].vn[d].y);

        }

      else if (e[i].side[d].ca == 2)
        {
          tp1[i].Fc[d][0] = e[i].ro * (e[i].u * e[i].vn[d].x + e[i].v * e[i].vn[d].y);

          tp1[i].Fc[d][1] = e[i].ro * e[i].u * (e[i].u * e[i].vn[d].x + e[i].v * e[i].vn[d].y) + e[i].vn[d].x * e[i].p;

          tp1[i].Fc[d][2] = e[i].ro * e[i].v * (e[i].u * e[i].vn[d].x + e[i].v * e[i].vn[d].y) + e[i].vn[d].y * e[i].p;

          tp1[i].Fc[d][3] = e[i].ro * (e[i].et + e[i].p / e[i].ro) * (e[i].u * e[i].vn[d].x + e[i].v * e[i].vn[d].y);
        }

      else if (e[i].side[d].ca == 4 || e[i].side[d].ca == 5)
        {
          tp1[i].Fc[d][0] = 0;

          tp1[i].Fc[d][1] = e[i].vn[d].x * e[i].p;

          tp1[i].Fc[d][2] = e[i].vn[d].y * e[i].p;

          tp1[i].Fc[d][3] = 0;

        }




    }

  //for viscous

  if(schemetype == 4 || schemetype == 5 || schemetype == 6)
    {
      double gradTx, gradTy, gradux, graduy, gradvx, gradvy, qx, qy;

      double taoxx, taoxy, taoyy, mu;

      double rx,ry, rm;


      for(int d = 0; d < 4; d++)
        {
          rx = (n[e[i].node[d]].x + n[e[i].node[(d + 1) % 4]].x) / 2. - e[i].x;

          ry = (n[e[i].node[d]].y + n[e[i].node[(d + 1) % 4]].y) / 2. - e[i].y;

          rm = sqrt(pow(rx, 2) + pow(ry, 2));



          if (e[i].side[d].ca == 1)
            {
              gradux = (uinf / Vinf - e[i].u) / rm * rx / rm;

              graduy = (uinf / Vinf - e[i].u) / rm * ry / rm;

              gradvx = (vinf / Vinf - e[i].v) / rm * rx / rm;

              gradvy = (vinf / Vinf - e[i].v) / rm * ry / rm;

              gradTx = (Tinf / Tinf - e[i].T) / rm * rx / rm;

              gradTy = (Tinf / Tinf - e[i].T) / rm * ry / rm;



              taoxx = (muinf / muinf) / Reinf * (4. / 3. * gradux - 2. / 3. * gradvy);

              taoxy = (muinf / muinf) / Reinf * (graduy + gradvx);

              taoyy = (muinf / muinf) / Reinf * (4. / 3. * gradvy - 2. / 3. * gradux);

              qx = (muinf / muinf) / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTx;

              qy = (muinf / muinf) / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTy;


              tp2[i].Fv[d][0] = 0;

              tp2[i].Fv[d][1] = e[i].vn[d].x * taoxx + e[i].vn[d].y * taoxy;

              tp2[i].Fv[d][2] = e[i].vn[d].x * taoxy + e[i].vn[d].y * taoyy;

              tp2[i].Fv[d][3] = e[i].vn[d].x * (uinf / Vinf * taoxx + vinf / Vinf * taoxy + qx) + e[i].vn[d].y * (uinf / Vinf * taoxy + vinf / Vinf * taoyy + qy);
            }

          else if (e[i].side[d].ca == 2)
            {
              gradux = (e[i].u - e[i].u) / rm * rx / rm;

              graduy = (e[i].u - e[i].u) / rm * ry / rm;

              gradvx = (e[i].v - e[i].v) / rm * rx / rm;

              gradvy = (e[i].v - e[i].v) / rm * ry / rm;

              gradTx = (e[i].T - e[i].T) / rm * rx / rm;

              gradTy = (e[i].T - e[i].T) / rm * ry / rm;


              mu = (1.461e-6 * pow(e[i].T * Tinf, (3. / 2.)) / (e[i].T * Tinf + 110.3)) / muinf;


              taoxx = mu / Reinf * (4. / 3. * gradux - 2. / 3. * gradvy);

              taoxy = mu / Reinf * (graduy + gradvx);

              taoyy = mu / Reinf * (4. / 3. * gradvy - 2. / 3. * gradux);

              qx = mu / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTx;

              qy = mu / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTy;


              tp2[i].Fv[d][0] = 0;

              tp2[i].Fv[d][1] = e[i].vn[d].x * taoxx + e[i].vn[d].y * taoxy;

              tp2[i].Fv[d][2] = e[i].vn[d].x * taoxy + e[i].vn[d].y * taoyy;

              tp2[i].Fv[d][3] = e[i].vn[d].x * (e[i].u * taoxx +  e[i].v * taoxy + qx) + e[i].vn[d].y * (e[i].u * taoxy + e[i].v * taoyy + qy);
            }

          else if (e[i].side[d].ca == 5)
            {
              gradux = (0. - e[i].u) / rm * rx / rm;

              graduy = (0. - e[i].u) / rm * ry / rm;

              gradvx = (0. - e[i].v) / rm * rx / rm;

              gradvy = (0. - e[i].v) / rm * ry / rm;

              gradTx = (e[i].T - e[i].T) / rm * rx / rm;

              gradTy = (e[i].T - e[i].T) / rm * ry / rm;



              taoxx = e[i].mu / Reinf * (4. / 3. * gradux - 2. / 3. * gradvy);

              taoxy = e[i].mu / Reinf * (graduy + gradvx);

              taoyy = e[i].mu / Reinf * (4. / 3. * gradvy - 2. / 3. * gradux);

              qx = e[i].mu / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTx;

              qy = e[i].mu / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTy;


              tp2[i].Fv[d][0] = 0;

              tp2[i].Fv[d][1] = e[i].vn[d].x * taoxx + e[i].vn[d].y * taoxy;

              tp2[i].Fv[d][2] = e[i].vn[d].x * taoxy + e[i].vn[d].y * taoyy;

              tp2[i].Fv[d][3] = e[i].vn[d].x * qx + e[i].vn[d].y * qy;

            }

        }
    }




}


void boundarycondition3(int i)
{

  for(int d = 0; d < 4; d++)
    {
      if (e[i].side[d].ca == 1)
        {
          tp1[i].Fc[d][0] = roinf / roinf * (uinf / Vinf * e[i].vn[d].x + vinf / Vinf * e[i].vn[d].y);

          tp1[i].Fc[d][1] = roinf / roinf * (uinf / Vinf) * (uinf / Vinf * e[i].vn[d].x + vinf / Vinf * e[i].vn[d].y) + e[i].vn[d].x * pinf / (roinf * pow(Vinf, 2));

          tp1[i].Fc[d][2] = roinf / roinf * (vinf / Vinf) * (uinf / Vinf * e[i].vn[d].x + vinf / Vinf * e[i].vn[d].y) + e[i].vn[d].y * pinf / (roinf * pow(Vinf, 2));

          tp1[i].Fc[d][3] = roinf / roinf * (((cv * Tinf + pow(Vinf, 2) / 2.) / pow(Vinf, 2)) + pinf / (roinf * pow(Vinf, 2))) * (uinf / Vinf * e[i].vn[d].x + vinf / Vinf * e[i].vn[d].y);

        }

      else if (e[i].side[d].ca == 2)
        {
          tp1[i].Fc[d][0] = e[i].ronew * (e[i].unew * e[i].vn[d].x + e[i].vnew * e[i].vn[d].y);

          tp1[i].Fc[d][1] = e[i].ronew * e[i].unew * (e[i].unew * e[i].vn[d].x + e[i].vnew * e[i].vn[d].y) + e[i].vn[d].x * e[i].pnew;

          tp1[i].Fc[d][2] = e[i].ronew * e[i].vnew * (e[i].unew * e[i].vn[d].x + e[i].vnew * e[i].vn[d].y) + e[i].vn[d].y * e[i].pnew;

          tp1[i].Fc[d][3] = e[i].ronew * (e[i].etnew + e[i].pnew / e[i].ronew) * (e[i].unew * e[i].vn[d].x + e[i].vnew * e[i].vn[d].y);
        }

      else if (e[i].side[d].ca == 4 || e[i].side[d].ca == 5)
        {
          tp1[i].Fc[d][0] = 0;

          tp1[i].Fc[d][1] = e[i].vn[d].x * e[i].pnew;

          tp1[i].Fc[d][2] = e[i].vn[d].y * e[i].pnew;

          tp1[i].Fc[d][3] = 0;

        }




    }

  //for viscous

  if(schemetype == 4 || schemetype == 5 || schemetype == 6)
    {
      double gradTx, gradTy, gradux, graduy, gradvx, gradvy, qx, qy;

      double taoxx, taoxy, taoyy, mu;

      double rx,ry, rm;


      for(int d = 0; d < 4; d++)
        {
          rx = (n[e[i].node[d]].x + n[e[i].node[(d + 1) % 4]].x) / 2. - e[i].x;

          ry = (n[e[i].node[d]].y + n[e[i].node[(d + 1) % 4]].y) / 2. - e[i].y;

          rm = sqrt(pow(rx, 2) + pow(ry, 2));



          if (e[i].side[d].ca == 1)
            {
              gradux = (uinf / Vinf - e[i].unew) / rm * rx / rm;

              graduy = (uinf / Vinf - e[i].unew) / rm * ry / rm;

              gradvx = (vinf / Vinf - e[i].vnew) / rm * rx / rm;

              gradvy = (vinf / Vinf - e[i].vnew) / rm * ry / rm;

              gradTx = (Tinf / Tinf - e[i].Tnew) / rm * rx / rm;

              gradTy = (Tinf / Tinf - e[i].Tnew) / rm * ry / rm;



              taoxx = (muinf / muinf) / Reinf * (4. / 3. * gradux - 2. / 3. * gradvy);

              taoxy = (muinf / muinf) / Reinf * (graduy + gradvx);

              taoyy = (muinf / muinf) / Reinf * (4. / 3. * gradvy - 2. / 3. * gradux);

              qx = (muinf / muinf) / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTx;

              qy = (muinf / muinf) / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTy;


              tp2[i].Fv[d][0] = 0;

              tp2[i].Fv[d][1] = e[i].vn[d].x * taoxx + e[i].vn[d].y * taoxy;

              tp2[i].Fv[d][2] = e[i].vn[d].x * taoxy + e[i].vn[d].y * taoyy;

              tp2[i].Fv[d][3] = e[i].vn[d].x * (uinf / Vinf * taoxx + vinf / Vinf * taoxy + qx) + e[i].vn[d].y * (uinf / Vinf * taoxy + vinf / Vinf * taoyy + qy);
            }

          else if (e[i].side[d].ca == 2)
            {
              gradux = (e[i].u - e[i].u) / rm * rx / rm;

              graduy = (e[i].u - e[i].u) / rm * ry / rm;

              gradvx = (e[i].v - e[i].v) / rm * rx / rm;

              gradvy = (e[i].v - e[i].v) / rm * ry / rm;

              gradTx = (e[i].T - e[i].T) / rm * rx / rm;

              gradTy = (e[i].T - e[i].T) / rm * ry / rm;


              mu = (1.461e-6 * pow(e[i].T * Tinf, (3. / 2.)) / (e[i].T * Tinf + 110.3)) / muinf;


              taoxx = mu / Reinf * (4. / 3. * gradux - 2. / 3. * gradvy);

              taoxy = mu / Reinf * (graduy + gradvx);

              taoyy = mu / Reinf * (4. / 3. * gradvy - 2. / 3. * gradux);

              qx = mu / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTx;

              qy = mu / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTy;


              tp2[i].Fv[d][0] = 0;

              tp2[i].Fv[d][1] = e[i].vn[d].x * taoxx + e[i].vn[d].y * taoxy;

              tp2[i].Fv[d][2] = e[i].vn[d].x * taoxy + e[i].vn[d].y * taoyy;

              tp2[i].Fv[d][3] = e[i].vn[d].x * (e[i].u * taoxx +  e[i].v * taoxy + qx) + e[i].vn[d].y * (e[i].u * taoxy + e[i].v * taoyy + qy);
            }

          else if (e[i].side[d].ca == 5)
            {
              gradux = (0. - e[i].unew) / rm * rx / rm;

              graduy = (0. - e[i].unew) / rm * ry / rm;

              gradvx = (0. - e[i].vnew) / rm * rx / rm;

              gradvy = (0. - e[i].vnew) / rm * ry / rm;

              gradTx = (e[i].Tnew - e[i].Tnew) / rm * rx / rm;

              gradTy = (e[i].Tnew - e[i].Tnew) / rm * ry / rm;



              taoxx = e[i].munew / Reinf * (4. / 3. * gradux - 2. / 3. * gradvy);

              taoxy = e[i].munew / Reinf * (graduy + gradvx);

              taoyy = e[i].munew / Reinf * (4. / 3. * gradvy - 2. / 3. * gradux);

              qx = e[i].munew / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTx;

              qy = e[i].munew / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTy;


              tp2[i].Fv[d][0] = 0;

              tp2[i].Fv[d][1] = e[i].vn[d].x * taoxx + e[i].vn[d].y * taoxy;

              tp2[i].Fv[d][2] = e[i].vn[d].x * taoxy + e[i].vn[d].y * taoyy;

              tp2[i].Fv[d][3] = e[i].vn[d].x * qx + e[i].vn[d].y * qy;

            }

        }
    }



}
