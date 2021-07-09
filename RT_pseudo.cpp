#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif
using std::cout;
using std::endl;
#include <string>


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <complex>
#include <cmath>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
#include <map>



// Define "Matrix" as a type. Solve problems using the Eigen package.
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> realMatrix;

using std::cout;
using std::endl;
using std::ifstream;
using namespace std;
using namespace std::complex_literals;

// Parse the one electron integrals files, and return a square matrix
Matrix readOneEFile(int nbasis, const char *filename) {
    FILE *data;
    data = fopen(filename, "r");

    Matrix matrix(nbasis, nbasis);

    for (int i = 0; i < nbasis; i++) {
        for (int j = 0; j < nbasis; j++) {
            fscanf(data, "%lf", &matrix(i, j));
        }
    }
    fclose(data);
    return matrix;
}

double computeElecEnergy(Matrix P, Matrix Hcore, Matrix F, int nbasis){
    double elecE = 0.000;
    for(int mu = 0; mu < nbasis; mu++){
        for(int nu = 0; nu < nbasis; nu++){
            elecE += (P(mu,nu) * (Hcore(mu,nu) + F(mu,nu))).real();
        }
    }
    return elecE;
}

int getCompoundIndex(int i, int j, int k, int l){
    int ij, kl, ijkl;

    if(i>j) {
        ij = i*(i+1)/2 + j;
    }else {
        ij = j*(j+1)/2 + i;
    }

    if(k>l) {
        kl = k*(k+1)/2 + l;
    }else {
        kl = l*(l+1)/2 + k;
    }

    if(ij > kl) {
        ijkl = ij*(ij+1)/2 + kl;
    }else {
        ijkl = kl*(kl+1)/2 + ij;
    }
    return ijkl;
}

// This function returns the field-free hamiltonian F'_0 (t)
// The dependence of F'_0 on t comes from the time dependence of P
Matrix formFock(Matrix Hcore, Matrix P, std::map<int, float> ERI, int nbasis) {
    Matrix F(nbasis, nbasis);
    for (int mu = 0; mu < nbasis; mu++) {
        for (int nu = 0; nu < nbasis; nu++) {
            F(mu, nu) = Hcore(mu, nu);
            for (int sig = 0; sig < nbasis; sig++) {
                for (int lam = 0; lam < nbasis; lam++) {
                    // Need to add one because these loop iterator are 0-indexed
                    // The basis function #s listed in the file are 1-indexed
                    int mu_nu_lam_sig = getCompoundIndex(mu + 1, nu + 1, lam + 1, sig + 1);
                    int mu_lam_nu_sig = getCompoundIndex(mu + 1, lam + 1, nu + 1, sig + 1);
                    F(mu, nu) += P(sig, lam) * std::complex<double>((2 * ERI[mu_nu_lam_sig] - ERI[mu_lam_nu_sig]));
                }
            }
        }
    }
    return  F;
}



int main() {
    int nbasis = 7;
    int nElec = 10;

/***********************************************************************
 * STO-3G RT Calculation for Water:
 ***********************************************************************/


// -----------------------------------------
//
// Step 1: Read in all the files
//
// -----------------------------------------


    // *****************
    // 1.1 Nuclear Repulsion
    // *****************
    FILE *Enuc_data;
    Enuc_data = fopen("enuc.dat", "r");
    double Enuc;
    fscanf(Enuc_data, "%lf", &Enuc);



    // *****************
    // 1.2 One-Electron Integrals
    // *****************

    Matrix S = readOneEFile(nbasis, "cq_s.dat");
    Matrix T = readOneEFile(nbasis, "cq_t.dat");
    Matrix V = readOneEFile(nbasis, "cq_v.dat");

    Matrix H_core = T+V;

    std::cout <<  "Overlap Matrix:"<< "\n";
    cout << S<< endl;
    cout << "" << endl;

    std::cout <<  "Core Hamilton:"<< "\n";
    cout << H_core << endl;
    cout << "" << endl;

    // *****************
    // 1.3 Two-Electron Integrals
    // *****************

    FILE *eri_data;
    int i, j, k, l;
    double val;

    std::map<int, float> ERI {};

    eri_data = fopen("eri.dat", "r");
    while(!feof(eri_data)) {
        fscanf(eri_data, "%d %d %d %d", &i, &j, &k, &l);
        fscanf(eri_data, "%lf",  &val);
        ERI.insert({getCompoundIndex(i,j,k,l), val });
    }
    fclose(eri_data);

    //for (auto& t : ERI)
    //    std::cout << t.first << " " << t.second<< "\n";



    // *****************
    // 1.4 Field-free Dipole Integrals
    // *****************

    Matrix Dx = readOneEFile(nbasis, "dipole_x.dat");
    Matrix Dy = readOneEFile(nbasis, "dipole_y.dat");
    Matrix Dz = readOneEFile(nbasis, "dipole_z.dat");

    Matrix Dipole[3];
    Dipole[0] = Dx;
    Dipole[1] = Dy;
    Dipole[2] = Dz;

    std::cout <<  "Field Free Dipole, x direction:"<< "\n";
    cout << Dx << endl;
    cout << "" << endl;

    std::cout <<  "Field Free Dipole, y direction:"<< "\n";
    cout << Dy << endl;
    cout << "" << endl;

    std::cout <<  "Field Free Dipole, z direction:"<< "\n";
    cout << Dz << endl;
    cout << "" << endl;


    // *****************
    // 1.5 Ground-State Density
    // *****************

    Matrix P0 = readOneEFile(nbasis, "my_density.dat");
    std::cout <<  "Ground-State Density:"<< "\n";
    cout << P0 << endl;
    cout << "" << endl;



// -----------------------------------------
//
// Step #2. Set Parameters:
//
// -----------------------------------------

    // *****************
    // 2.1 Define Electric Field
    // *****************

    // Here I use a delta pulse. I turn it on for the first time step, and turn
    // it off for all future time steps.

    // E is the same for all three directions

    // E is given in atomic units. e/(a0^2) = 5.14x10^9 V/cm

    double E = 0.005;


    // *****************
    // 2.2 Time step & Total time
    // *****************

    double dt = 0.05;
    double totalTime = 15;
    double curTime = 0.00;


// -----------------------------------------
//
// Step #3. First Density matrix propagation:
//
// -----------------------------------------


    // We want to use MMUT to propagate the density matrix. To get P at t_(i+1), we use P at t_(i-1)
    // Obviously this is not doable for first propagation, so we use Trapezoidal rule

    // *****************
    // 3.1 Get field-free Fock matrix
    // *****************

    Matrix F = formFock(H_core, P0, ERI, nbasis);

    std::cout <<  "First Fock matrix from RT program:"<< "\n";
    cout << F << endl;
    cout << "" << endl;

    // *****************
    // 3.2 Get field-dependent Fock matrix
    // *****************

    /*
    for (int direction = 0; direction < 3; direction++){
        for (int i = 0; i < nbasis; i++) {
            for (int j = 0; j < nbasis; j++) {
                F(i,j) += E * Dipole[direction](i,j);
            }
        }
    }
    */

    // *****************
    // 3.3 Transform F into F' (AO basis -> OAO basis)
    // *****************

    // Symmetric Orthogolization:
    // Want to find X+ S X = 1
    // Need to find X = S^(-1/2)
    // Diagonlize S: A+ S A = s, s is diagonal
    // X = S^(-1/2) = A s^(-1/2) A+
    // F' = X+ F X

    // Diagonlize S: We get A+ S A = s, where s is diagonal
    Eigen::SelfAdjointEigenSolver<Matrix> diagS(S);
    Matrix A = diagS.eigenvectors();
    realMatrix s = diagS.eigenvalues();

    // Obtain s^(-1/2)
    // Raise the diagonal matrix s to the power of -1/2
    Matrix s_min_half (nbasis, nbasis);
    Matrix s_half (nbasis, nbasis);
    for(int i = 0; i < nbasis; i++){
        for(int j = 0; j < nbasis; j++){
            if (i == j){
                s_min_half(i,j) = pow(s(i,0), (-1.0/2.0));
            }else{
                s_min_half(i,j) = 0;
            }
        }
    }


    // Build the symmetric orthogonalization matrix S^(-1/2), or X
    // X = S^(-1/2) = A s^(-1/2) A+
    Matrix X = A * s_min_half * A.adjoint();
    Matrix X_inverse = X.inverse();

    // Transform F to F'
    // F' = X+ F X
    Matrix F_prime = X.adjoint() * F * X;

    std::cout <<  "F':"<< "\n";
    cout << F_prime << endl;
    cout << "" << endl;

    // *****************
    // 3.4 Obtain Unitary transformation matrix, U
    // *****************

    // Equation to use: C^dagger F'(t_i) C = epsilon(t_i)
    //                  U = C e^{i dt epsilon(t_i)} C^dagger


    // First we diagonalize F' to get C and epsilon
    Eigen::SelfAdjointEigenSolver<Matrix> diagF(F_prime);
    Matrix C = diagF.eigenvectors();
    realMatrix epsilon = diagF.eigenvalues();

    // Then we get the middle term: e^{i dt epsilon(t_i)}
    Matrix e_i_dt_epsilon (nbasis, nbasis);
    for(int i = 0; i < nbasis; i++){
        for(int j = 0; j < nbasis; j++){
            if (i == j){
                e_i_dt_epsilon(i,j) = exp(1i * dt * epsilon(i,0));
            }else{
                e_i_dt_epsilon(i,j) = 0;
            }
        }
    }

    // U = C e^{i dt epsilon(t_i)} C^dagger
    Matrix U = C * e_i_dt_epsilon * C.adjoint();

    // *****************
    // 3.5 Calculate P for the first time step
    // *****************

    // Equation to use: P(t_(i+1)) = U P_0(t_i) U^dagger

    Matrix curP = U * P0 * U.adjoint();


    std::cout <<  "Density matrix P at the first time step:"<< "\n";
    cout << curP << endl;
    cout << "" << endl;


// -----------------------------------------
//
// Step #4. Propagating the future time steps
//
// -----------------------------------------

    double curE;
    Matrix preP = P0;
    Matrix nextP;
    Matrix AO_P;
    Matrix AO_F = F;
    // To get around the fencepost problem, we need to get P for the next step

    while (curTime < totalTime) {
        AO_P = X * curP * X.adjoint();

        curE = computeElecEnergy(AO_P, H_core, AO_F, nbasis);
        printf("\nAt t= %7.6f, the energy is %17.14f", curTime, curE);

        // Using the current density, we can get the
        // field-free Fock matrix F0 in the AO basis
        AO_F = formFock(H_core, AO_P, ERI, nbasis);

        // We are supposed to do: F = F0 + d * e
        // However, we are using a delta pulse. Electric field is turned off
        // after the first time step. Thus F = F0

        // Orthagonalize F, convert it into OAO basis
        F_prime = X.adjoint() * AO_F * X;

        // Diagonalize F' to get C and epsilon
        Eigen::SelfAdjointEigenSolver<Matrix> diagFprime(F_prime);
        C = diagFprime.eigenvectors();
        epsilon = diagFprime.eigenvalues();

        // Then we get the middle term: e^{i 2dt epsilon(t_i)}
        Matrix e_i_2dt_epsilon (nbasis, nbasis);
        for(int i = 0; i < nbasis; i++){
            for(int j = 0; j < nbasis; j++){
                if (i == j){
                    e_i_2dt_epsilon(i,j) = exp(1i * (2 * dt)  * epsilon(i,0));
                }else{
                    e_i_2dt_epsilon(i,j) = 0.00;
                }
            }
        }

        // U = C e^{i 2 dt epsilon(t_i)} C^dagger
        U = C * e_i_2dt_epsilon * C.adjoint();

        // Equation to use: P(t_(i+1)) = U P_0(t_(i-1)) U^dagger
        nextP = U * preP * U.adjoint();

        // Update:
        preP = curP;
        curP = nextP;
        curTime += dt;
    }
    return 0;

}


