#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <data.h>
#include <mesh.h>
#include <constants.h>
#include <initialconditions.h>
#include <fluxes.h>
#include <boundaryconditions.h>
#include <calculations.h>
#include <postprocess.h>
#include <omp.h>



using namespace std;

extern vector<node> n;
extern vector<element> e;
extern vector<type1> tp1;
extern vector<type2> tp2;
extern vector<type3> tp3;

int it;

vector <double> t1;


double L2norm(string variable)
{
  double var, varprev;


  double sum = 0;

  for(unsigned int i = 1; i < e.size(); i++)
    {
      if(variable == "u")
        {
          var = e[i].u;
          varprev = e[i].uprev;

        }
      else if(variable == "v")
        {
          var = e[i].v;
          varprev = e[i].vprev;
        }
      else if(variable == "ro")
        {
          var = e[i].ro;
          varprev = e[i].roprev;
        }
      else if(variable == "p")
        {
          var = e[i].p;
          varprev = e[i].pprev;
        }
      else if(variable == "et")
        {
          var = e[i].et;
          varprev = e[i].etprev;
        }


      if(e[i].c == 0)
        sum += pow((var - varprev) / var, 2);


    }

  sum = sqrt(sum);


  return sum;
}


void tstepinviscid()
{
  t1.resize((e.size()), 10);

  double c, lambdaci, lambdacj;

#pragma omp parallel for schedule(dynamic,1000) private(c, lambdaci, lambdacj)
  for (unsigned int i = 1; i < e.size(); i++)
    {
      c = sqrt(gama * e[i].p / e[i].ro);

      /* lambdaci = (abs(e[i].u * (e[i].vn[1].x - e[i].vn[3].x) / 2. + e[i].v * (e[i].vn[1].y - e[i].vn[3].y) / 2.) + c) * (e[i].dS[1] + e[i].dS[3]) / 2.;

          lambdacj = (abs(e[i].u * (e[i].vn[2].x - e[i].vn[0].x) / 2. + e[i].v * (e[i].vn[2].y - e[i].vn[0].y) / 2.) + c) * (e[i].dS[2] + e[i].dS[0]) / 2.;


          t1[i] = C * e[i].area / (lambdaci + lambdacj); */

      t1[i] = C * e[i].area / ((abs(e[i].u) + c) * (abs(e[i].Sx[0]) + abs(e[i].Sx[1]) + abs(e[i].Sx[2]) + abs(e[i].Sx[3])) / 2.
          + (abs(e[i].v) + c) * (abs(e[i].Sy[0]) + abs(e[i].Sy[1]) + abs(e[i].Sy[2]) + abs(e[i].Sy[3])) / 2.);


      // t1[i] = 1e-3;
    }


}


void tstepviscous()
{
  t1.resize((e.size()), 10);

  double c;

#pragma omp parallel for schedule(dynamic,1000) private(c)
  for (unsigned int i = 1; i < e.size(); i++)
    {
          c = sqrt(gama * e[i].p / e[i].ro);

          /*   Si = (abs(e[i].Sx[0]) + abs(e[i].Sx[1]) + abs(e[i].Sx[2]) + abs(e[i].Sx[3])) / 2.;

           Sj = (abs(e[i].Sy[0]) + abs(e[i].Sy[1]) + abs(e[i].Sy[2]) + abs(e[i].Sy[3])) / 2.;

          lambdaci = (abs(e[i].u) + c) * Si;

          lambdacj = (abs(e[i].v) + c) * Sj;


          lambdavi = max((4. / (3. * e[i].ro)), (gama / e[i].ro)) * e[i].mu / Pr * pow(Si, 2) / e[i].area;

          lambdavj = max((4. / (3. * e[i].ro)), (gama / e[i].ro)) * e[i].ro / Pr * pow(Sj, 2) / e[i].area;


          t1[i] = C * e[i].area / (lambdaci + lambdacj + 2. * (lambdavi + lambdavj)); */


          //t1[i] = 1e-4;

          t1[i] = C * e[i].area / ((abs(e[i].u) + 1. * e[i].mu / e[i].ro) * (abs(e[i].Sx[0]) + abs(e[i].Sx[1]) + abs(e[i].Sx[2]) + abs(e[i].Sx[3])) / 2.
              + (abs(e[i].v) + 1. * e[i].mu / e[i].ro) * (abs(e[i].Sy[0]) + abs(e[i].Sy[1]) + abs(e[i].Sy[2]) + abs(e[i].Sy[3])) / 2.);


    }

}


bool converge()
{
  double sum = 0;

#pragma omp parallel for schedule(dynamic,1000) reduction(+:sum)
  for(unsigned int i = 1; i < e.size(); i++)
    {
      if(e[i].c == 0)
        sum += pow((e[i].p - e[i].pprev) / e[i].p, 2);
    }

  sum = sqrt(sum);

  cout << sum << endl;

  if(sum < 1.e-8)
    {
      cout << "Converged?" << endl;

      return 1;

    }

  return 0;
}


void primitiveinviscid(double A[], int i)
{
  e[i].ro = A[0];
  e[i].u = A[1] / e[i].ro;
  e[i].v = A[2] / e[i].ro;
  e[i].et = A[3] / e[i].ro;

  e[i].V = sqrt(pow(e[i].u, 2) + pow(e[i].v, 2));
  e[i].T = (e[i].et * pow(Vinf, 2) - pow(e[i].V * Vinf, 2) / 2.) / (cv * Tinf);
  e[i].p = e[i].ro * roinf * R * e[i].T * Tinf / (roinf * pow(Vinf, 2));
}


void primitiveviscous(double A[], int i)
{
  e[i].ro = A[0];
  e[i].u = A[1] / e[i].ro;
  e[i].v = A[2] / e[i].ro;
  e[i].et = A[3] / e[i].ro;

  e[i].V = sqrt(pow(e[i].u, 2) + pow(e[i].v, 2));
  e[i].T = (e[i].et * pow(Vinf, 2) - pow(e[i].V * Vinf, 2) / 2.) / (cv * Tinf);
  e[i].p = e[i].ro * roinf * R * e[i].T * Tinf / (roinf * pow(Vinf, 2));

  e[i].mu = (1.461e-6 * pow(e[i].T * Tinf, (3. / 2.)) / (e[i].T * Tinf + 110.3)) / muinf;
}


void iterateinviscid()
{
  tstepinviscid();

#pragma omp parallel for schedule(dynamic,1000)
  for (unsigned int i = 1; i < e.size(); i++)
    {
      e[i].pprev = e[i].p;

      e[i].roprev = e[i].ro;

      e[i].uprev = e[i].u;

      e[i].vprev = e[i].v;

      e[i].etprev = e[i].et;


      W(i);
      Fc1(i);

      boundarycondition(i);


      for (int d = 0; d < 4; d++)
        {
          tp1[i].W[d] = tp1[i].W[d] - t1[i] * (tp1[i].Fc[0][d] * e[i].dS[0] + tp1[i].Fc[1][d] * e[i].dS[1] + tp1[i].Fc[2][d] * e[i].dS[2] + tp1[i].Fc[3][d] * e[i].dS[3]) / e[i].area;
        }

      primitiveinviscid(tp1[i].W,i);


    }



}


void iterateviscous()
{
  tstepviscous();

#pragma omp parallel for schedule(dynamic,1000)
  for (unsigned int i = 1; i < e.size(); i++)
    {


      e[i].pprev = e[i].p;

      e[i].roprev = e[i].ro;

      e[i].uprev = e[i].u;

      e[i].vprev = e[i].v;

      e[i].etprev = e[i].et;


      W(i);
      Fc1(i);
      Fv1(i);


      boundarycondition(i);


      for (int d = 0; d < 4; d++)
        {
          tp1[i].W[d] = tp1[i].W[d] - t1[i] / e[i].area * ((tp1[i].Fc[0][d] - tp2[i].Fv[0][d]) * e[i].dS[0] + (tp1[i].Fc[1][d] - tp2[i].Fv[1][d]) * e[i].dS[1] + (tp1[i].Fc[2][d] - tp2[i].Fv[2][d]) * e[i].dS[2] + (tp1[i].Fc[3][d] - tp2[i].Fv[3][d]) * e[i].dS[3]);
        }


      primitiveviscous(tp1[i].W,i);
    }


}


void iterate3stageinviscid()
{
  tstepinviscid();

  double alpha1,alpha2,alpha3;

  alpha1 = 0.1481;

  alpha2 = 0.4;

  alpha3 = 1.;


  double W0[4];

#pragma omp parallel for schedule(dynamic,1000) private(W0)
  for (unsigned int i = 1; i < e.size(); i++)
    {

      e[i].pprev = e[i].p;

      e[i].roprev = e[i].ro;

      e[i].uprev = e[i].u;

      e[i].vprev = e[i].v;

      e[i].etprev = e[i].et;


      W(i);
      Fc1(i);

      if(e[i].c != 0)
        boundarycondition(i);


      for(int d = 0; d < 4; d++)
        {
          W0[d] = tp1[i].W[d];
        }



      for (int d = 0; d < 4; d++) //first stage
        {
          tp1[i].W[d] = W0[d] - alpha1 * t1[i] / e[i].area * (tp1[i].Fc[0][d] * e[i].dS[0] + tp1[i].Fc[1][d] * e[i].dS[1] + tp1[i].Fc[2][d] * e[i].dS[2] + tp1[i].Fc[3][d] * e[i].dS[3]);
        }

      primitiveinviscid(tp1[i].W,i);





      W(i);
      Fc1(i);

      if(e[i].c != 0)
        boundarycondition(i);


      for (int d = 0; d < 4; d++)  //second stage
        {
          tp1[i].W[d] = W0[d] - alpha2 * t1[i] / e[i].area * (tp1[i].Fc[0][d] * e[i].dS[0] + tp1[i].Fc[1][d] * e[i].dS[1] + tp1[i].Fc[2][d] * e[i].dS[2] + tp1[i].Fc[3][d] * e[i].dS[3]);
        }



      primitiveinviscid(tp1[i].W,i);




      W(i);
      Fc1(i);

      if(e[i].c != 0)
        boundarycondition(i);


      for (int d = 0; d < 4; d++)  //third stage
        {
          tp1[i].W[d] = W0[d] - alpha3 * t1[i] / e[i].area * (tp1[i].Fc[0][d] * e[i].dS[0] + tp1[i].Fc[1][d] * e[i].dS[1] + tp1[i].Fc[2][d] * e[i].dS[2] + tp1[i].Fc[3][d] * e[i].dS[3]);
        }


      primitiveinviscid(tp1[i].W,i);

    }

}


void iterate3stageviscous()
{
  tstepviscous();

#pragma omp parallel for schedule(dynamic,1000)
  for (unsigned int i = 1; i < e.size(); i++)
    {

      e[i].pprev = e[i].p;

      e[i].roprev = e[i].ro;

      e[i].uprev = e[i].u;

      e[i].vprev = e[i].v;

      e[i].etprev = e[i].et;



      double alpha1,alpha2,alpha3;

      alpha1 = 0.1481;

      alpha2 = 0.4;

      alpha3 = 1.;


      double W0[4];


      W(i);
      Fc1(i);
      Fv1(i);

      if(e[i].c != 0)
        boundarycondition(i);


      for(int d = 0; d < 4; d++)
        {
          W0[d] = tp1[i].W[d];
        }



      for (int d = 0; d < 4; d++) //first stage
        {
          tp1[i].W[d] = W0[d] - alpha1 * t1[i] / e[i].area * ((tp1[i].Fc[0][d] - tp2[i].Fv[0][d]) * e[i].dS[0] + (tp1[i].Fc[1][d] - tp2[i].Fv[1][d]) * e[i].dS[1] + (tp1[i].Fc[2][d] - tp2[i].Fv[2][d]) * e[i].dS[2] + (tp1[i].Fc[3][d] - tp2[i].Fv[3][d]) * e[i].dS[3]);
        }

      e[i].ro = tp1[i].W[0];
      e[i].u = tp1[i].W[1] / e[i].ro;
      e[i].v = tp1[i].W[2] / e[i].ro;
      e[i].et = tp1[i].W[3] / e[i].ro;

      e[i].V = sqrt(pow(e[i].u, 2) + pow(e[i].v, 2));
      e[i].T = (e[i].et * pow(uinf, 2) - pow(e[i].V * uinf, 2) / 2.) / (cv * Tinf);
      e[i].p = e[i].ro * roinf * R * e[i].T * Tinf / (roinf * pow(uinf, 2));

      e[i].mu = pow(e[i].T, 3. / 2.) * (Tinf + 110) / (e[i].T * Tinf + 110);





      W(i);
      Fc1(i);

      if(e[i].c != 0)
        boundarycondition(i);


      for (int d = 0; d < 4; d++)  //second stage
        {
          tp1[i].W[d] = W0[d] - alpha2 * t1[i] / e[i].area * (tp1[i].Fc[0][d] * e[i].dS[0] + tp1[i].Fc[1][d] * e[i].dS[1] + tp1[i].Fc[2][d] * e[i].dS[2] + tp1[i].Fc[3][d] * e[i].dS[3]);
        }



      e[i].ro = tp1[i].W[0];
      e[i].u = tp1[i].W[1] / e[i].ro;
      e[i].v = tp1[i].W[2] / e[i].ro;
      e[i].et = tp1[i].W[3] / e[i].ro;

      e[i].V = sqrt(pow(e[i].u, 2) + pow(e[i].v, 2));
      e[i].T = (e[i].et * pow(uinf, 2) - pow(e[i].V * uinf, 2) / 2.) / (cv * Tinf);
      e[i].p = e[i].ro * roinf * R * e[i].T * Tinf / (roinf * pow(uinf, 2));
 e[i].mu = pow(e[i].T, 3. / 2.) * (Tinf + 110) / (e[i].T * Tinf + 110);

      W(i);
      Fc1(i);
      Fv1(i);


      if(e[i].c != 0)
        boundarycondition(i);


      for (int d = 0; d < 4; d++)  //third stage
        {
          tp1[i].W[d] = W0[d] - alpha3 * t1[i] / e[i].area * ((tp1[i].Fc[0][d] - tp2[i].Fv[0][d]) * e[i].dS[0] + (tp1[i].Fc[1][d] - tp2[i].Fv[1][d]) * e[i].dS[1] + (tp1[i].Fc[2][d] - tp2[i].Fv[2][d]) * e[i].dS[2] + (tp1[i].Fc[3][d] - tp2[i].Fv[3][d]) * e[i].dS[3]);
        }


      e[i].ro = tp1[i].W[0];
      e[i].u = tp1[i].W[1] / e[i].ro;
      e[i].v = tp1[i].W[2] / e[i].ro;
      e[i].et = tp1[i].W[3] / e[i].ro;

      e[i].V = sqrt(pow(e[i].u, 2) + pow(e[i].v, 2));
      e[i].T = (e[i].et * pow(uinf, 2) - pow(e[i].V * uinf, 2) / 2.) / (cv * Tinf);
      e[i].p = e[i].ro * roinf * R * e[i].T * Tinf / (roinf * pow(uinf, 2));

      e[i].mu = pow(e[i].T, 3. / 2.) * (Tinf + 110) / (e[i].T * Tinf + 110);
    }

}


void calculations(int maxit)
{
  switch (schemetype)
    {
    case 1:
      {
        ofstream myfile;

        myfile.open("L2norms.txt");

        for (it = 1; it <= maxit; it++)
          {
            iterateinviscid();

            /*    if (it % 1000 == 0)
                  {
                    output(it);

                    myfile << it << " " << L2norm("p") << " " << L2norm("ro") << " " << L2norm("u") << " " << L2norm("v") << " " << L2norm("et") << endl;
                  } */

            if(it % writebackupflag == 0)
              writebackup();


            cout << "Iteration " << it << endl;

            if (converge())
              {
                output(it);
                break;
              }


          }
      }
      break;
    case 2:
      {
        for (it = 1; it <= maxit; it++)
          {
            iterate3stageinviscid();

            if (it % 1000 == 0)
              {
                output(it);
              }

            if(it % writebackupflag == 0)
              writebackup();


          //  cout << "Iteration " << it << endl;

            if (converge())
              {
                output(it);
                break;
              }


          }

      }
      break;
    case 3:
      {
        ofstream myfile;

        myfile.open("output.dat");

        myfile << "VARIABLES = \"X\" , \"Y\" , \"p\" , \"ro\" , \"T\" , \"u\" , \"v\" , \"et\" , \"Mach\" , \"Gradient\""<< endl;


        double alpha1,alpha2,alpha3;

        alpha1 = 0.1481;

        alpha2 = 0.4;

        alpha3 = 1.;


        double dt = physicaldeltat * Vinf;


        for (it = 1; it <= maxit; it++)
          {
            double ctr;

            double a;

            tstepinviscid();

#pragma omp parallel for schedule(dynamic,1000) private(ctr, a)
            for (unsigned int i = 1; i < e.size(); i++)
              {
                ctr = 0;

                W3(i, 2);

                W3(i, 1);

                W3(i, 0);

               for(int d = 0; d < 4; d++)
                  {
                  tp3[i].Q[d] = 2. / dt * e[i].area * tp3[i].W[d][1] - 1. / (2. * dt) * e[i].area * tp3[i].W[d][0];
                  }


               do{
                    Fc1new(i);

                    if(e[i].c != 0)
                      boundarycondition3(i);


                    for(int d = 0; d < 4; d++)
                      tp3[i].W[d][1] = tp3[i].W[d][2];


                    for(int d = 0; d < 4; d++) //first stage
                      {
                        tp3[i].W[d][2] = tp3[i].W[d][1] - alpha1 * t1[i] / e[i].area * pow((1. + 3. / (2. * dt) * alpha1 * t1[i]), -1) * (tp1[i].Fc[0][d] * e[i].dS[0] + tp1[i].Fc[1][d] * e[i].dS[1] + tp1[i].Fc[2][d] * e[i].dS[2] + tp1[i].Fc[3][d] * e[i].dS[3] + 3. / (2. * dt) * e[i].area * tp3[i].W[d][1] - tp3[i].Q[d]);
                      }


                    e[i].ronew = tp3[i].W[0][2];
                    e[i].unew = tp3[i].W[1][2] / e[i].ronew;
                    e[i].vnew = tp3[i].W[2][2] / e[i].ronew;
                    e[i].etnew = tp3[i].W[3][2] / e[i].ronew;

                    e[i].Vnew = sqrt(pow(e[i].unew, 2) + pow(e[i].vnew, 2));
                    e[i].Tnew = (e[i].etnew * pow(Vinf, 2) - pow(e[i].Vnew * Vinf, 2) / 2.) / (cv * Tinf);
                    e[i].pnew = e[i].ronew * roinf * R * e[i].Tnew * Tinf / (roinf * pow(Vinf, 2));



                    Fc1new(i);

                    if(e[i].c != 0)
                      boundarycondition3(i);




                    for(int d = 0; d < 4; d++) //second stage
                      {

                        tp3[i].W[d][2] = tp3[i].W[d][1] - alpha2 * t1[i] / e[i].area * pow((1. + 3. / (2. * dt) * alpha2 * t1[i]), -1) * (tp1[i].Fc[0][d] * e[i].dS[0] + tp1[i].Fc[1][d] * e[i].dS[1] + tp1[i].Fc[2][d] * e[i].dS[2] + tp1[i].Fc[3][d] * e[i].dS[3] + 3. / (2. * dt) * e[i].area * tp3[i].W[d][1] - tp3[i].Q[d]);

                      }



                    e[i].ronew = tp3[i].W[0][2];
                    e[i].unew = tp3[i].W[1][2] / e[i].ronew;
                    e[i].vnew = tp3[i].W[2][2] / e[i].ronew;
                    e[i].etnew = tp3[i].W[3][2] / e[i].ronew;

                    e[i].Vnew = sqrt(pow(e[i].unew, 2) + pow(e[i].vnew, 2));
                    e[i].Tnew = (e[i].etnew * pow(Vinf, 2) - pow(e[i].Vnew * Vinf, 2) / 2.) / (cv * Tinf);
                    e[i].pnew = e[i].ronew * roinf * R * e[i].Tnew * Tinf / (roinf * pow(Vinf, 2));



                    Fc1new(i);

                    if(e[i].c != 0)
                      boundarycondition3(i);



                    for(int d = 0; d < 4; d++) //third stage
                      {
                        tp3[i].W[d][2] = tp3[i].W[d][1] - alpha3 * t1[i] / e[i].area * pow((1. + 3. / (2. * dt) * alpha3 * t1[i]), -1) * (tp1[i].Fc[0][d] * e[i].dS[0] + tp1[i].Fc[1][d] * e[i].dS[1] + tp1[i].Fc[2][d] * e[i].dS[2] + tp1[i].Fc[3][d] * e[i].dS[3] + 3. / (2. * dt) * e[i].area * tp3[i].W[d][1] - tp3[i].Q[d]);
                      }

                    e[i].ronew = tp3[i].W[0][2];
                    e[i].unew = tp3[i].W[1][2] / e[i].ronew;
                    e[i].vnew = tp3[i].W[2][2] / e[i].ronew;
                    e[i].etnew = tp3[i].W[3][2] / e[i].ronew;

                    e[i].Vnew = sqrt(pow(e[i].unew, 2) + pow(e[i].vnew, 2));
                    e[i].Tnew = (e[i].etnew * pow(Vinf, 2) - pow(e[i].Vnew * Vinf, 2) / 2.) / (cv * Tinf);
                    e[i].pnew = e[i].ronew * roinf * R * e[i].Tnew * Tinf / (roinf * pow(Vinf, 2));

                    ctr++;

                    a = abs(tp3[i].W[0][2] - tp3[i].W[0][1]) / tp3[i].W[0][1];

                  }
                while(a > 1e-8);

              }


#pragma omp parallel for schedule(dynamic,1000)
            for (unsigned int i = 1; i < e.size(); i++)
              {
                Fc1new1(i);

                if(e[i].c != 0)
                  boundarycondition3(i);


                for (int d = 0; d < 4; d++)
                  {
                    tp3[i].W[d][2] = tp3[i].W[d][2] - dt / e[i].area * (tp1[i].Fc[0][d] * e[i].dS[0] + tp1[i].Fc[1][d] * e[i].dS[1] + tp1[i].Fc[2][d] * e[i].dS[2] + tp1[i].Fc[3][d] * e[i].dS[3]);
                  }

                e[i].ro = tp3[i].W[0][2];
                e[i].u = tp3[i].W[1][2] / e[i].ro;
                e[i].v = tp3[i].W[2][2] / e[i].ro;
                e[i].et = tp3[i].W[3][2] / e[i].ro;

                e[i].V = sqrt(pow(e[i].u, 2) + pow(e[i].v, 2));
                e[i].T = (e[i].et * pow(Vinf, 2) - pow(e[i].V * Vinf, 2) / 2.) / (cv * Tinf);
                e[i].p = e[i].ro * roinf * R * e[i].T * Tinf / (roinf * pow(Vinf, 2));

              }



#pragma omp parallel for schedule(dynamic,1000)
            for (unsigned int i = 1; i < e.size(); i++)
              {
                e[i].roprev = e[i].ro;
                e[i].uprev = e[i].u;
                e[i].vprev = e[i].v;
                e[i].etprev = e[i].et;

                e[i].Vprev = e[i].V;
                e[i].Tprev = e[i].T;
                e[i].pprev = e[i].p;


                e[i].ro = e[i].ronew;
                e[i].u = e[i].unew;
                e[i].v = e[i].vnew;
                e[i].et = e[i].etnew;

                e[i].V = e[i].Vnew;
                e[i].T = e[i].Tnew;
                e[i].p = e[i].pnew;
              }


            if (it % 1000 == 0)
              {
                output(it);
              }

            if(it % writebackupflag == 0)
              writebackup();

            cout << "Iteration " << it << endl;

            if (converge())
              {
                //output(it);
                //break;
              }

          }

        myfile.close();
      }
      break;
    case 4:
      {
        for (it = 1; it <= maxit; it++)
          {
            iterateviscous();

            if (it % 1000 == 0)
              {
                output(it);
              }

            if(it % writebackupflag == 0)
              writebackup();



            cout << "Iteration " << it << endl;

            if (converge())
              {
                output(it);
                break;
              }


          }


      }
      break;
    case 5:
      {
        for (it = 1; it <= maxit; it++)
          {
            iterate3stageviscous();

           if (it % 100000 == 0)
              {
                output(it);
              }

            if(it % writebackupflag == 0)
              writebackup();


            cout << "Iteration " << it << endl;

            if (converge())
              {
                output(it);
                break;
              }


          }


      }
      break;
    case 6:
      {

      }
      break;
    default:
      cout << "This scheme is not available.";
      break;
    }
}



