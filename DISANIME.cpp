/********************************************************************************
*
*   Copyright (C) 2018 Culham Centre for Fusion Energy,
*   United Kingdom Atomic Energy Authority, Oxfordshire OX14 3DB, UK
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*
********************************************************************************
*
*   Program: DISANIME
*            Dislocation segment strain field in anisotropic medium
*   Version: 1.0
*   Date:    14 May 2019
*   Author:  Pui-Wai (Leo) MA
*   Contact: leo.ma@ukaea.uk
*   Address: Culham Centre for Fusion Energy, OX14 3DB, United Kingdom
*
********************************************************************************
* 
*   It only works for small dislocation segment. 
*
*   It is based on the formula in:
*   Micromechanics of Defects in Solids, Second, Revised Edition. Page 47. 
*   by Toshio Mura
*
*   Input units
*   Length: Angstrom
*   Elastic constants: GPa
*
*   Output units
*   Strain tensor: Unitless
*
********************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

using namespace std;

/** Prototypes ******************************************************************/

void Matrix_Inversion(double M[3][3], double Inv_M[3][3]);
int Levi_Civita(int i, int j, int k);

void make_Cijkl(double Cijkl[3][3][3][3], double Cij[6][6]);
void make_Gik(double Cijkl[3][3][3][3], double r_vec[3], double Gir[3][3]);
void make_Gik_j(double Cijkl[3][3][3][3], double r_vec[3], double Gir_s[3][3][3]);
void make_Gik_jl(double Cijkl[3][3][3][3], double r_vec[3], double Gir_sm[3][3][3][3]);

void read_elastic(double Cij[6][6]);

void calculate_strain(double Cijkl[3][3][3][3], double r_vec[3], double Burgers[3], double L_vec[3], double epsilon[3][3]);

/** Main program ****************************************************************/

int main(int argc, char **argv){


    double Cij[6][6]; //in unit of GPa
    double Cijkl[3][3][3][3]; //in unit of eV per Angstrom^3

    read_elastic(Cij);

    //convert Cij into Cijkl
    make_Cijkl(Cijkl, Cij); 

    //input variables
    double L_vec[3]; //Line segment vector
    double Burgers[3] = {1.0, 0.0, 0.0}; //Burgers vector
    double r_vec[3]; // vector from the position of a segment to a spatial point. 

    //output variable
    double epsilon[3][3]; // strain field at a spatial point induce by a segment


    double Total_L = 100000e0; // total dislocation line length

    double L1 = 0.1e0; //each segment length
    double L2 = 1e0;
    double L3 = 10e0;
  
    int N1 = int(Total_L / L1 + 1e-8);
    int N2 = int(Total_L / L2 + 1e-8);
    int N3 = int(Total_L / L3 + 1e-8);
    
    double total_epsilon[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
            total_epsilon[m][n] = 0.0;

    for (int i = 0; i < N1; ++i){

         r_vec[0] = (i - N1/2 + 0.5)*L1; 
         r_vec[1] = 0.0;
         r_vec[2] = 50e0;

         L_vec[0] = L1;
         L_vec[1] = 0.0;
         L_vec[2] = 0.0;

         calculate_strain(Cijkl, r_vec, Burgers, L_vec, epsilon);

         for (int m = 0; m < 3; ++m){
             for (int n = 0; n < 3; ++n){
                 total_epsilon[m][n] += epsilon[m][n];
             }
         }
    }   
    cout << setiosflags(ios::scientific) << setprecision(8);
    cout << total_epsilon[0][0] << " " << total_epsilon[0][1] << " " << total_epsilon[0][2] << " " 
         << total_epsilon[1][0] << " " << total_epsilon[1][1] << " " << total_epsilon[1][2] << " " 
         << total_epsilon[2][0] << " " << total_epsilon[2][1] << " " << total_epsilon[2][2] << "\n";

    for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
            total_epsilon[m][n] = 0.0;

    for (int i = 0; i < N2; ++i){

         r_vec[0] = (i - N2/2 + 0.5)*L2;
         r_vec[1] = 0.0;
         r_vec[2] = 50e0;

         L_vec[0] = L2;
         L_vec[1] = 0.0;
         L_vec[2] = 0.0;

         calculate_strain(Cijkl, r_vec, Burgers, L_vec, epsilon);

         for (int m = 0; m < 3; ++m){
             for (int n = 0; n < 3; ++n){
                 total_epsilon[m][n] += epsilon[m][n];
             }
         }
    }   
    cout << total_epsilon[0][0] << " " << total_epsilon[0][1] << " " << total_epsilon[0][2] << " " 
         << total_epsilon[1][0] << " " << total_epsilon[1][1] << " " << total_epsilon[1][2] << " " 
         << total_epsilon[2][0] << " " << total_epsilon[2][1] << " " << total_epsilon[2][2] << "\n";

    for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
            total_epsilon[m][n] = 0.0;

    for (int i = 0; i < N3; ++i){

         r_vec[0] = (i - N3/2 + 0.5)*L3;
         r_vec[1] = 0.0;
         r_vec[2] = 50e0;

         L_vec[0] = L3;
         L_vec[1] = 0.0;
         L_vec[2] = 0.0;

         calculate_strain(Cijkl, r_vec, Burgers, L_vec, epsilon);

         for (int m = 0; m < 3; ++m){
             for (int n = 0; n < 3; ++n){
                 total_epsilon[m][n] += epsilon[m][n];
             }
         }
    }   
    cout << total_epsilon[0][0] << " " << total_epsilon[0][1] << " " << total_epsilon[0][2] << " " 
         << total_epsilon[1][0] << " " << total_epsilon[1][1] << " " << total_epsilon[1][2] << " " 
         << total_epsilon[2][0] << " " << total_epsilon[2][1] << " " << total_epsilon[2][2] << "\n";


    
    //isotropic analytic solution
    double epsilon_analytic[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    epsilon_analytic[2][0] = -Burgers[0] / (4*3.1415926)*r_vec[1]/(r_vec[1]*r_vec[1] + r_vec[2]*r_vec[2]);
    epsilon_analytic[1][0] = -Burgers[0] / (4*3.1415926)*r_vec[2]/(r_vec[1]*r_vec[1] + r_vec[2]*r_vec[2]);
    epsilon_analytic[0][2] = epsilon_analytic[2][0];
    epsilon_analytic[0][1] = epsilon_analytic[1][0];

    cout << epsilon_analytic[0][0] << " " << epsilon_analytic[0][1] << " " << epsilon_analytic[0][2] << " " 
         << epsilon_analytic[1][0] << " " << epsilon_analytic[1][1] << " " << epsilon_analytic[1][2] << " " 
         << epsilon_analytic[2][0] << " " << epsilon_analytic[2][1] << " " << epsilon_analytic[2][2] << "\n";


}            

/* functions *******************************************************************/

void calculate_strain(double Cijkl[3][3][3][3], double r_vec[3], double Burgers[3], double L_vec[3], double epsilon[3][3]){

    //Toshio Mura -- Micromechanics of Defects in Solids, Second, Revised Edition. Page 47. 
    //Assume short straight dislocation segment, apply mean value theorem.
    //Assume Gip,q(r-r') is the same for the whole segment.

    double beta[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            beta[i][j] = 0e0;

    double Gip_q[3][3][3];

    make_Gik_j(Cijkl, r_vec, Gip_q);

    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            for (int h = 0; h < 3; ++h){
                for (int p = 0; p < 3; ++p){
                    for (int q = 0; q < 3; ++q){
                        for (int m = 0; m < 3; ++m){
                            for (int n = 0; n < 3; ++n){
                                int epsilon_jnh = Levi_Civita(j,n,h);
                                if (epsilon_jnh != 0){
                                    beta[j][i] += epsilon_jnh*Cijkl[p][q][m][n]*Gip_q[i][p][q]*Burgers[m]*L_vec[h];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    //symmetrize epsilon
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            epsilon[i][j] = (beta[i][j] + beta[j][i])/2e0;

}
    




               
void make_Cijkl(double C[3][3][3][3], double C_Voigt[6][6]){

    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            for (int k = 0; k < 3; ++k){
                for (int l = 0; l < 3; ++l){
                    C[i][j][k][l] = 0e0;
                }
            }
        }
    }
    int p, q, r, s;
    for (int i = 0; i < 6; ++i){
        if (i == 0){
            p = 0;
            q = 0;
        } else if (i == 1){
            p = 1;
            q = 1;
        } else if (i == 2){
            p = 2;
            q = 2;
        } else if (i == 3){
            p = 1;
            q = 2;
        } else if (i == 4){
            p = 2;
            q = 0;
        } else if (i == 5){
            p = 0;
            q = 1;
        }
        for (int j = 0; j < 6; ++j){
            if (j == 0){
                r = 0;
                s = 0;
            } else if (j == 1){
                r = 1;
                s = 1;
            } else if (j == 2){
                r = 2;
                s = 2;
            } else if (j == 3){
                r = 1;
                s = 2;
            } else if (j == 4){
                r = 2;
                s = 0;
            } else if (j == 5){
                r = 0;
                s = 1;
            }
            if (i < 3 && j < 3){
                C[p][q][r][s] = C_Voigt[i][j];
                C[q][p][r][s] = C_Voigt[i][j];
                C[p][q][s][r] = C_Voigt[i][j];
                C[q][p][s][r] = C_Voigt[i][j];
            } else if ((i < 3 && j < 6) || (i < 6 && j < 3)){
                C[p][q][r][s] = C_Voigt[i][j];
                C[q][p][r][s] = C_Voigt[i][j];
                C[p][q][s][r] = C_Voigt[i][j];
                C[q][p][s][r] = C_Voigt[i][j];
            } else {
                C[p][q][r][s] = C_Voigt[i][j];
                C[q][p][r][s] = C_Voigt[i][j];
                C[p][q][s][r] = C_Voigt[i][j];
                C[q][p][s][r] = C_Voigt[i][j];
            }
        }
    }
}

void make_Gik(double Cijkl[3][3][3][3], double r_vec[3], double Gir[3][3]){

    //Gik is calculated according to D. M. Barnett Phys. Stat. Sol. (b) 49, 741 (1972)

    double Pi = 3.141592653589793;

    double r0 = sqrt(pow(r_vec[0],2) + pow(r_vec[1],2) + pow(r_vec[2],2));

    double T[3];
    for (int i = 0; i < 3; ++i) T[i] = r_vec[i]/r0;

    //Note: We followed the spherical coordinate system that is used in Barnett's paper. 

    double cosPhi, sinPhi, cosTheta, sinTheta;
    if (fabs(T[0]) < 1e-8 && fabs(T[1]) < 1e-8){
        cosPhi = 1e0;
        sinPhi = 0e0;
        cosTheta = 1e0;
        sinTheta = 0e0;
    } else {
        cosPhi = T[2];
        sinPhi = sin(acos(T[2]));
        cosTheta=T[0]/sinPhi;
        sinTheta=T[1]/sinPhi;
    }

    double a[3], b[3];
    a[0] = sinTheta;
    a[1] = -cosTheta;
    a[2] = 0e0;
    b[0] = cosPhi*cosTheta;
    b[1] = cosPhi*sinTheta;
    b[2] = -sinPhi;

    int NPsi = 20;
    //Note: Using NPsi = 20 should give answer up to 4 significant figures correct.
    double dPsi = Pi/NPsi;

    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
            Gir[i][r] = 0e0;

    for (int nn = 0; nn < NPsi; ++nn){ //perform integration

        double Psi = Pi/NPsi*(nn+0.5);
        double cosPsi = cos(Psi);
        double sinPsi = sin(Psi);

        double z[3];
        for (int i = 0; i < 3; ++i) z[i] = a[i]*cosPsi + b[i]*sinPsi;

        double M[3][3], Inv_M[3][3];

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                M[i][j] = 0e0;

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int r = 0; r < 3; ++r)
                    for (int s = 0; s < 3; ++s)
                        if (Cijkl[i][j][r][s] > 0e0) M[i][r] += Cijkl[i][j][r][s]*z[j]*z[s];

        Matrix_Inversion(M, Inv_M);

        for (int i = 0; i < 3; ++i)
            for (int r = 0; r < 3; ++r)
                Gir[i][r] += Inv_M[i][r];
    }

    double factor = dPsi/(4e0*pow(Pi,2)*r0);
    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
            Gir[i][r] *= factor;
    
}

void make_Gik_j(double Cijkl[3][3][3][3], double r_vec[3], double Gir_s[3][3][3]){

    //Gik_j is calculated according to D. M. Barnett Phys. Stat. Sol. (b) 49, 741 (1972)

    double Pi = 3.141592653589793;

    double rsq = pow(r_vec[0],2) + pow(r_vec[1],2) + pow(r_vec[2],2);
    double r0 = sqrt(rsq);

    double T[3];
    for (int i = 0; i < 3; ++i) T[i] = r_vec[i]/r0;

    //Note: We followed the spherical coordinate system that is used in Barnett's paper. 

    double cosPhi, sinPhi, cosTheta, sinTheta;
    if (fabs(T[0]) < 1e-8 && fabs(T[1]) < 1e-8){
        cosPhi = 1e0;
        sinPhi = 0e0;
        cosTheta = 1e0;
        sinTheta = 0e0;
    } else {
        cosPhi = T[2];
        sinPhi = sin(acos(T[2]));
        cosTheta=T[0]/sinPhi;
        sinTheta=T[1]/sinPhi;
    }

    double a[3], b[3];
    a[0] = sinTheta;
    a[1] = -cosTheta;
    a[2] = 0e0;
    b[0] = cosPhi*cosTheta;
    b[1] = cosPhi*sinTheta;
    b[2] = -sinPhi;

    int NPsi = 20;
    //Note: Using NPsi = 20 should give answer up to 4 significant figures correct.
    double dPsi = Pi/NPsi;

    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
            for (int s = 0; s < 3; ++s)
                    Gir_s[i][r][s] = 0e0;

    for (int nn = 0; nn < NPsi; ++nn){ //perform integration

        double Psi = Pi/NPsi*(nn+0.5);
        double cosPsi = cos(Psi);
        double sinPsi = sin(Psi);

        double z[3];
        for (int i = 0; i < 3; ++i) z[i] = a[i]*cosPsi + b[i]*sinPsi;

        double M[3][3], Inv_M[3][3], F[3][3];

        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                M[i][j] = 0e0;
                F[i][j] = 0e0;
            }
        }

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int r = 0; r < 3; ++r)
                    for (int s = 0; s < 3; ++s)
                        if (Cijkl[i][j][r][s] > 0e0) M[i][r] += Cijkl[i][j][r][s]*z[j]*z[s];

        Matrix_Inversion(M, Inv_M);

        for (int i = 0; i < 3; ++i)
            for (int r = 0; r < 3; ++r)
                for (int j = 0; j < 3; ++j)
                    for (int p = 0; p < 3; ++p)
                        for (int n = 0; n < 3; ++n)
                            for (int w = 0; w < 3; ++w)
                                if (Cijkl[j][p][n][w] > 0e0) 
                                    F[i][r] += Cijkl[j][p][n][w]*Inv_M[i][j]*Inv_M[n][r]*(z[p]*T[w] + z[w]*T[p]);

        for (int i = 0; i < 3; ++i)
            for (int r = 0; r < 3; ++r)
                for (int s = 0; s < 3; ++s)
                        Gir_s[i][r][s] += -T[s]*Inv_M[i][r] + z[s]*F[i][r];
    }

    double factor = dPsi/(4e0*pow(Pi,2)*rsq);
    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
            for (int s = 0; s < 3; ++s)
                    Gir_s[i][r][s] *= factor;
    
}

void make_Gik_jl(double Cijkl[3][3][3][3], double r_vec[3], double Gir_sm[3][3][3][3]){

    //Gik_jl is calculated according to D. M. Barnett Phys. Stat. Sol. (b) 49, 741 (1972)

    double Pi = 3.141592653589793;

    double r0 = sqrt(pow(r_vec[0],2) + pow(r_vec[1],2) + pow(r_vec[2],2));

    double T[3];
    for (int i = 0; i < 3; ++i) T[i] = r_vec[i]/r0;

    //Note: We followed the spherical coordinate system that is used in Barnett's paper. 

    double cosPhi, sinPhi, cosTheta, sinTheta;
    if (fabs(T[0]) < 1e-8 && fabs(T[1]) < 1e-8){
        cosPhi = 1e0;
        sinPhi = 0e0;
        cosTheta = 1e0;
        sinTheta = 0e0;
    } else {
        cosPhi = T[2];
        sinPhi = sin(acos(T[2]));
        cosTheta=T[0]/sinPhi;
        sinTheta=T[1]/sinPhi;
    }

    double a[3], b[3];
    a[0] = sinTheta;
    a[1] = -cosTheta;
    a[2] = 0e0;
    b[0] = cosPhi*cosTheta;
    b[1] = cosPhi*sinTheta;
    b[2] = -sinPhi;

    int NPsi = 20; 
    //Note: Using NPsi = 20 should give answer up to 4 significant figures correct.
    double dPsi = Pi/NPsi;

    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
            for (int s = 0; s < 3; ++s)
                for (int m = 0; m < 3; ++m)
                    Gir_sm[i][r][s][m] = 0e0;

    for (int nn = 0; nn < NPsi; ++nn){ //perform integration

        double Psi = Pi/NPsi*(nn+0.5);
        double cosPsi = cos(Psi);
        double sinPsi = sin(Psi);

        double z[3];
        for (int i = 0; i < 3; ++i) z[i] = a[i]*cosPsi + b[i]*sinPsi;

        double M[3][3], Inv_M[3][3], A[3][3], F[3][3];

        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                M[i][j] = 0e0;
                A[i][j] = 0e0;
                F[i][j] = 0e0;
            }
        }

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int r = 0; r < 3; ++r)
                    for (int s = 0; s < 3; ++s)
                        if (Cijkl[i][j][r][s] > 0e0) M[i][r] += Cijkl[i][j][r][s]*z[j]*z[s];

        Matrix_Inversion(M, Inv_M);

        for (int i = 0; i < 3; ++i)
            for (int r = 0; r < 3; ++r)
                for (int j = 0; j < 3; ++j)
                    for (int p = 0; p < 3; ++p)
                        for (int n = 0; n < 3; ++n)
                            for (int w = 0; w < 3; ++w)
                                if (Cijkl[j][p][n][w] > 0e0) 
                                    F[i][r] += Cijkl[j][p][n][w]*Inv_M[i][j]*Inv_M[n][r]*(z[p]*T[w] + z[w]*T[p]);


        for (int i = 0; i < 3; ++i)
            for (int r = 0; r < 3; ++r)
                for (int j = 0; j < 3; ++j)
                    for (int p = 0; p < 3; ++p)
                        for (int n = 0; n < 3; ++n)
                            for (int w = 0; w < 3; ++w)
                                if (Cijkl[j][p][n][w] > 0e0) 
                                    A[i][r] += Cijkl[j][p][n][w]*((z[p]*T[w]+z[w]*T[p])*(F[i][j]*Inv_M[n][r]+Inv_M[i][j]*F[n][r])
                                               -2e0*Inv_M[i][j]*Inv_M[n][r]*T[p]*T[w]);

        for (int i = 0; i < 3; ++i)
            for (int r = 0; r < 3; ++r)
                for (int s = 0; s < 3; ++s)
                    for (int m = 0; m < 3; ++m)
                        Gir_sm[i][r][s][m] += 2e0*T[s]*T[m]*Inv_M[i][r] - 2e0*(z[s]*T[m]+z[m]*T[s])*F[i][r] + z[s]*z[m]*A[i][r];
    }

    double factor = dPsi/(4e0*pow(Pi,2)*pow(r0,3));
    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
            for (int s = 0; s < 3; ++s)
                for (int m = 0; m < 3; ++m)
                    Gir_sm[i][r][s][m] *= factor;
    
}

void Matrix_Inversion(double M[3][3], double Inv_M[3][3]){

    Inv_M[0][0] = M[1][1]*M[2][2] - M[2][1]*M[1][2];
    Inv_M[0][1] = M[0][2]*M[2][1] - M[0][1]*M[2][2];
    Inv_M[0][2] = M[0][1]*M[1][2] - M[0][2]*M[1][1];
    Inv_M[1][0] = M[1][2]*M[2][0] - M[1][0]*M[2][2];
    Inv_M[1][1] = M[0][0]*M[2][2] - M[0][2]*M[2][0];
    Inv_M[1][2] = M[1][0]*M[0][2] - M[0][0]*M[1][2];
    Inv_M[2][0] = M[1][0]*M[2][1] - M[2][0]*M[1][1];
    Inv_M[2][1] = M[2][0]*M[0][1] - M[0][0]*M[2][1];
    Inv_M[2][2] = M[0][0]*M[1][1] - M[1][0]*M[0][1];

    double Det = M[0][0]*Inv_M[0][0]
               + M[0][1]*Inv_M[1][0]
               + M[0][2]*Inv_M[2][0];

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            Inv_M[i][j] /= Det;
}

int Levi_Civita(int i, int j, int k){

    int output;
    if (((i == 0) && (j == 1) && (k == 2)) ||
        ((i == 1) && (j == 2) && (k == 0)) ||
        ((i == 2) && (j == 0) && (k == 1))){
        output = 1;
    } else if (((i == 2) && (j == 1) && (k == 0)) ||
               ((i == 0) && (j == 2) && (k == 1)) ||
               ((i == 1) && (j == 0) && (k == 2))){
        output = -1;
    } else {
        output = 0;
    }
    return output;
}

void read_elastic(double Cij[6][6]) {

    // read in Cij in Voigt notation 
    ifstream in_file_1("input_elastic");

    if (in_file_1) {
        cout << "Start reading input_elastic ... \n\n";
    } else {
        cout << "ERROR: You need to have an input_elastic file \n\n";
        exit(1);
    }
   
    string line;  //dummy
    getline(in_file_1,line); //read comment line
    getline(in_file_1,line); //read comment line

    for (int i = 0; i < 6; ++i)
         in_file_1 >> Cij[i][0] >> Cij[i][1] >> Cij[i][2] >> Cij[i][3] >> Cij[i][4] >> Cij[i][5];
    
    double eVperAcubic_to_GPa = 160.21766208;

    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            Cij[i][j] /= eVperAcubic_to_GPa;

    in_file_1.close();

}

