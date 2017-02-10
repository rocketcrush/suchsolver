#include <iostream>
#include <vector>
#include <cmath>
#include <data.h>
#include <constants.h>
#include <mesh.h>
#include <initialconditions.h>


using namespace std;

extern vector<node> n;
extern vector<element> e;
extern vector<type1> tp1;
extern vector<type2> tp2;
extern vector<type3> tp3;


void W(int i)
{
  tp1[i].W[0] = e[i].ro;


  tp1[i].W[1] = e[i].ro * e[i].u;

  tp1[i].W[2] = e[i].ro * e[i].v;

  tp1[i].W[3] = e[i].ro * e[i].et;
}


void Fc1(int i)
{
  double fenergyp, fenergym, fmassp, fmassm;

  double ML,MR;

  double VL, VR; //contravariant velocity of right state and left state

  double cL, cR; //speed of sound

  double roL, uL, vL, pL, etL;

  double roR, uR, vR, pR, etR;

  double vnx, vny;

  double Fp[4], Fm[4];

  roL = e[i].ro;
  uL = e[i].u;
  vL = e[i].v;
  pL = e[i].p;
  etL = e[i].et;

  for(int d = 0; d < 4; d++)
    {
      roR = e[e[i].neigh[d]].ro;
      uR = e[e[i].neigh[d]].u;
      vR = e[e[i].neigh[d]].v;
      pR = e[e[i].neigh[d]].p;
      etR = e[e[i].neigh[d]].et;


      vnx = e[i].vn[d].x;
      vny = e[i].vn[d].y;


      VL = vnx * uL + vny * vL;
      VR = vnx * uR + vny * vR;


      cL = sqrt(gama * pL / roL);
      cR = sqrt(gama * pR / roR);


      ML = VL / cL;
      MR = VR / cR;


      if(ML >= 1.)
        {
          Fp[0] = roL * VL;

          Fp[1] = roL * uL * VL + vnx * pL;

          Fp[2] = roL * vL * VL + vny * pL;

          Fp[3] = roL * (etL + pL / roL) * VL;

        }

      else if(abs(ML) < 1)
        {
          fmassp = roL * cL * pow((ML + 1), 2) / 4.;

          fenergyp = fmassp * ((pow(((gama - 1.) * VL + 2. * cL),2) / (2. * (pow(gama, 2) - 1.))) + (pow(uL, 2) + pow(vL, 2) - pow(VL, 2)) / 2.);

          Fp[0] = fmassp;

          Fp[1] = fmassp * (vnx * (-VL + 2. * cL) / gama + uL);

          Fp[2] = fmassp * (vny * (-VL + 2. * cL) / gama + vL);

          Fp[3] = fenergyp;
        }

      else
        {
          Fp[0] = 0;

          Fp[1] = 0;

          Fp[2] = 0;

          Fp[3] = 0;
        }



      if(MR >= 1.)
        {
          Fm[0] = 0;

          Fm[1] = 0;

          Fm[2] = 0;

          Fm[3] = 0;
        }

      else if(abs(MR) < 1)
        {
          fmassm = -roR * cR * pow((MR - 1), 2) / 4.;

          fenergym = fmassm * ((pow(((gama - 1.) * VR - 2. * cR),2) / (2. * (pow(gama, 2) - 1.))) + (pow(uR, 2) + pow(vR, 2) - pow(VR, 2)) / 2.);

          Fm[0] = fmassm;

          Fm[1] = fmassm * (vnx * (-VR - 2. * cR) / gama + uR);

          Fm[2] = fmassm * (vny * (-VR - 2. * cR) / gama + vR);

          Fm[3] = fenergym;
        }

      else
        {
          Fm[0] = roR * VR;

          Fm[1] = roR * uR * VR + vnx * pR;

          Fm[2] = roR * vR * VR + vny * pR;

          Fm[3] = roR * (etR + pR / roR) * VR;
        }


      for(int a = 0; a < 4; a++)
        tp1[i].Fc[d][a] = Fp[a] + Fm[a];
    }

}


void W3(int i, int t)
{
  if(t == 0)
    {
      tp3[i].W[0][t] = e[i].roprev;

      tp3[i].W[1][t] = e[i].roprev * e[i].uprev;

      tp3[i].W[2][t] = e[i].roprev * e[i].vprev;

      tp3[i].W[3][t] = e[i].roprev * e[i].etprev;
    }
  else if(t == 1)
    {
      tp3[i].W[0][t] = e[i].ro;

      tp3[i].W[1][t] = e[i].ro * e[i].u;

      tp3[i].W[2][t] = e[i].ro * e[i].v;

      tp3[i].W[3][t] = e[i].ro * e[i].et;
    }
  else if(t == 2)
    {
      tp3[i].W[0][t] = e[i].ronew;

      tp3[i].W[1][t] = e[i].ronew * e[i].unew;

      tp3[i].W[2][t] = e[i].ronew * e[i].vnew;

      tp3[i].W[3][t] = e[i].ronew * e[i].etnew;
    }

}


void Fc1new(int i)
{
  double fenergyp, fenergym, fmassp, fmassm;

  double ML,MR;

  double VL, VR; //contravariant velocity of right state and left state

  double cL, cR; //speed of sound

  double roL, uL, vL, pL, etL;

  double roR, uR, vR, pR, etR;

  double vnx, vny;

  double Fp[4], Fm[4];

  roL = e[i].ronew;
  uL = e[i].unew;
  vL = e[i].vnew;
  pL = e[i].pnew;
  etL = e[i].etnew;

  for(int d = 0; d < 4; d++)
    {
      roR = e[e[i].neigh[d]].ro;
      uR = e[e[i].neigh[d]].u;
      vR = e[e[i].neigh[d]].v;
      pR = e[e[i].neigh[d]].p;
      etR = e[e[i].neigh[d]].et;


      vnx = e[i].vn[d].x;
      vny = e[i].vn[d].y;


      VL = vnx * uL + vny * vL;
      VR = vnx * uR + vny * vR;


      cL = sqrt(gama * pL / roL);
      cR = sqrt(gama * pR / roR);


      ML = VL / cL;
      MR = VR / cR;


      if(ML >= 1.)
        {
          Fp[0] = roL * VL;

          Fp[1] = roL * uL * VL + vnx * pL;

          Fp[2] = roL * vL * VL + vny * pL;

          Fp[3] = roL * (etL + pL / roL) * VL;

        }

      else if(abs(ML) < 1)
        {
          fmassp = roL * cL * pow((ML + 1), 2) / 4.;

          fenergyp = fmassp * ((pow(((gama - 1.) * VL + 2. * cL),2) / (2. * (pow(gama, 2) - 1.))) + (pow(uL, 2) + pow(vL, 2) - pow(VL, 2)) / 2.);

          Fp[0] = fmassp;

          Fp[1] = fmassp * (vnx * (-VL + 2. * cL) / gama + uL);

          Fp[2] = fmassp * (vny * (-VL + 2. * cL) / gama + vL);

          Fp[3] = fenergyp;
        }

      else
        {
          Fp[0] = 0;

          Fp[1] = 0;

          Fp[2] = 0;

          Fp[3] = 0;
        }



      if(MR >= 1.)
        {
          Fm[0] = 0;

          Fm[1] = 0;

          Fm[2] = 0;

          Fm[3] = 0;
        }

      else if(abs(MR) < 1)
        {
          fmassm = -roR * cR * pow((MR - 1), 2) / 4.;

          fenergym = fmassm * ((pow(((gama - 1.) * VR - 2. * cR),2) / (2. * (pow(gama, 2) - 1.))) + (pow(uR, 2) + pow(vR, 2) - pow(VR, 2)) / 2.);

          Fm[0] = fmassm;

          Fm[1] = fmassm * (vnx * (-VR - 2. * cR) / gama + uR);

          Fm[2] = fmassm * (vny * (-VR - 2. * cR) / gama + vR);

          Fm[3] = fenergym;
        }

      else
        {
          Fm[0] = roR * VR;

          Fm[1] = roR * uR * VR + vnx * pR;

          Fm[2] = roR * vR * VR + vny * pR;

          Fm[3] = roR * (etR + pR / roR) * VR;
        }


      for(int a = 0; a < 4; a++)
        tp1[i].Fc[d][a] = Fp[a] + Fm[a];
    }


}


void Fc1new1(int i)
{

  double fenergyp, fenergym, fmassp, fmassm;

  double ML,MR;

  double VL, VR; //contravariant velocity of right state and left state

  double cL, cR; //speed of sound

  double roL, uL, vL, pL, etL;

  double roR, uR, vR, pR, etR;

  double vnx, vny;

  double Fp[4], Fm[4];

  roL = e[i].ronew;
  uL = e[i].unew;
  vL = e[i].vnew;
  pL = e[i].pnew;
  etL = e[i].etnew;

  for(int d = 0; d < 4; d++)
    {
      roR = e[e[i].neigh[d]].ronew;
      uR = e[e[i].neigh[d]].unew;
      vR = e[e[i].neigh[d]].vnew;
      pR = e[e[i].neigh[d]].pnew;
      etR = e[e[i].neigh[d]].etnew;


      vnx = e[i].vn[d].x;
      vny = e[i].vn[d].y;


      VL = vnx * uL + vny * vL;
      VR = vnx * uR + vny * vR;


      cL = sqrt(gama * pL / roL);
      cR = sqrt(gama * pR / roR);


      ML = VL / cL;
      MR = VR / cR;


      if(ML >= 1.)
        {
          Fp[0] = roL * VL;

          Fp[1] = roL * uL * VL + vnx * pL;

          Fp[2] = roL * vL * VL + vny * pL;

          Fp[3] = roL * (etL + pL / roL) * VL;

        }

      else if(abs(ML) < 1)
        {
          fmassp = roL * cL * pow((ML + 1), 2) / 4.;

          fenergyp = fmassp * ((pow(((gama - 1.) * VL + 2. * cL),2) / (2. * (pow(gama, 2) - 1.))) + (pow(uL, 2) + pow(vL, 2) - pow(VL, 2)) / 2.);

          Fp[0] = fmassp;

          Fp[1] = fmassp * (vnx * (-VL + 2. * cL) / gama + uL);

          Fp[2] = fmassp * (vny * (-VL + 2. * cL) / gama + vL);

          Fp[3] = fenergyp;
        }

      else
        {
          Fp[0] = 0;

          Fp[1] = 0;

          Fp[2] = 0;

          Fp[3] = 0;
        }



      if(MR >= 1.)
        {
          Fm[0] = 0;

          Fm[1] = 0;

          Fm[2] = 0;

          Fm[3] = 0;
        }

      else if(abs(MR) < 1)
        {
          fmassm = -roR * cR * pow((MR - 1), 2) / 4.;

          fenergym = fmassm * ((pow(((gama - 1.) * VR - 2. * cR),2) / (2. * (pow(gama, 2) - 1.))) + (pow(uR, 2) + pow(vR, 2) - pow(VR, 2)) / 2.);

          Fm[0] = fmassm;

          Fm[1] = fmassm * (vnx * (-VR - 2. * cR) / gama + uR);

          Fm[2] = fmassm * (vny * (-VR - 2. * cR) / gama + vR);

          Fm[3] = fenergym;
        }

      else
        {
          Fm[0] = roR * VR;

          Fm[1] = roR * uR * VR + vnx * pR;

          Fm[2] = roR * vR * VR + vny * pR;

          Fm[3] = roR * (etR + pR / roR) * VR;
        }


      for(int a = 0; a < 4; a++)
        tp1[i].Fc[d][a] = Fp[a] + Fm[a];
    }


}


void Fv1(int i)
{
  double ro, u, v, p, et, mu;

  double vnx, vny;

  double gradTx, gradTy, gradux, graduy, gradvx, gradvy, qx, qy;

  double taoxx, taoxy, taoyy;

  for(int d = 0; d < 4; d++)
    {
      vnx = e[i].vn[d].x;
      vny = e[i].vn[d].y;

      ro = (e[i].ro + e[e[i].neigh[d]].ro) / 2.;
      u = (e[i].u + e[e[i].neigh[d]].u) / 2.;
      v = (e[i].v + e[e[i].neigh[d]].v) / 2.;
      p = (e[i].ro + e[e[i].neigh[d]].p) / 2.;
      et = (e[i].et + e[e[i].neigh[d]].et) / 2.;
      mu = (e[i].mu + e[e[i].neigh[d]].mu) / 2.;


      gradux = (e[e[i].neigh[d]].u - e[i].u) / e[i].rm[d] * e[i].r[d].x / e[i].rm[d];

      graduy = (e[e[i].neigh[d]].u - e[i].u) / e[i].rm[d] * e[i].r[d].y / e[i].rm[d];

      gradvx = (e[e[i].neigh[d]].v - e[i].v) / e[i].rm[d] * e[i].r[d].x / e[i].rm[d];

      gradvy = (e[e[i].neigh[d]].v - e[i].v) / e[i].rm[d] * e[i].r[d].y / e[i].rm[d];

      gradTx = (e[e[i].neigh[d]].T - e[i].T) / e[i].rm[d] * e[i].r[d].x / e[i].rm[d];

      gradTy = (e[e[i].neigh[d]].T - e[i].T) / e[i].rm[d] * e[i].r[d].y / e[i].rm[d];


      taoxx = mu / Reinf * (4. / 3. * gradux - 2. / 3. * gradvy);

      taoxy = mu / Reinf * (graduy + gradvx);

      taoyy = mu / Reinf * (4. / 3. * gradvy - 2. / 3. * gradux);

      qx = mu / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTx;

      qy = mu / (Reinf * Pr * (gama - 1.) * pow(Machinf,2)) * gradTy;


      tp2[i].Fv[d][0] = 0;

      tp2[i].Fv[d][1] = vnx * taoxx + vny * taoxy;

      tp2[i].Fv[d][2] = vnx * taoxy + vny * taoyy;

      tp2[i].Fv[d][3] = vnx * (u * taoxx + v * taoxy + qx) + vny * (u * taoxy + v * taoyy + qy);

    }
}
