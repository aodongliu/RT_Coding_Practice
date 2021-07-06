#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
#include <map>

// Define "Matrix" as a type. Solve problems using the Eigen package.
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
using std::cout;
using std::endl;
using std::ifstream;
using namespace std;

// Parse the one electron integrals files, and return a square matrix
Matrix readOneEFile(int nbasis, const char *filename) {
    FILE *data;
    data = fopen(filename, "r");

    Matrix matrix(nbasis, nbasis);

    int nline = nbasis * (nbasis + 1) / 2;
    for (int i = 0; i < nline; i++) {
        int row;
        int column;

        // Can not read all three values at the same time!
        fscanf(data, "%d %d", &row, &column);
        fscanf(data, "%lf", &matrix(row - 1, column - 1));
        matrix(column - 1, row - 1) = matrix(row - 1, column - 1);
    }
    fclose(data);
    return matrix;
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

Matrix computeDensity(int nbasis, int nElec, Matrix C){
    Matrix P(nbasis, nbasis);
    for(int sig = 0; sig < nbasis; sig++){
        for(int lam = 0; lam < nbasis; lam++){
            P(sig, lam) = 0.000;
            for (int i = 0; i < nElec/2 ; i++){
                P(sig, lam) += C(sig, i) * C(lam, i);
            }
        }
    }
    return P;
}

double computeElecEnergy(Matrix P, Matrix Hcore, Matrix F, int nbasis){
    double elecE = 0.000;
    for(int mu = 0; mu < nbasis; mu++){
        for(int nu = 0; nu < nbasis; nu++){
            elecE += P(mu,nu) * (Hcore(mu,nu) + F(mu,nu));
        }
    }
    return elecE;
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
                    F(mu, nu) += P(sig, lam) * (2 * ERI[mu_nu_lam_sig] - ERI[mu_lam_nu_sig]);
                }
            }
        }
    }

    return  F;
}

double computeDensityDifference(Matrix curP, Matrix preP, int nbasis){
    double rmsP = 0.000;
    for(int mu = 0; mu < nbasis; mu++){
        for(int nu = 0; nu < nbasis; nu++){
            rmsP += pow((curP(mu, nu) - preP(mu, nu)),2);
        }
    }
    return sqrt(rmsP);
}

int main() {
    int nbasis = 7;
    int nElec = 10;

/***********************************************************************
 * STO-3G RT Calculation for Water:
 ***********************************************************************/

// -----------------------------------------
//
// Step 1: Nuclear Repulsion Energy
//
// -----------------------------------------
    FILE *Enuc_data;
    Enuc_data = fopen("enuc.dat", "r");
    double Enuc;
    fscanf(Enuc_data, "%lf", &Enuc);


// -----------------------------------------
//
// Step 2: One-Electron Integrals
//
// -----------------------------------------

    Matrix S = readOneEFile(nbasis, "s.dat");
    Matrix T = readOneEFile(nbasis, "t.dat");
    Matrix V = readOneEFile(nbasis, "v.dat");

    // *****************
    // TO IMPLEMENT:
    // *****************
    //Matrix dipole = readOneEFile(nbasis, "d.dat");

    // *****************
    // TO IMPLEMENT:
    // *****************
    //Matrix field = ???

    Matrix H_core = T+V;
    cout << H_core;



// -----------------------------------------
//
// Step 3: Two-Electron Integrals
//
// -----------------------------------------
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

// -----------------------------------------
//
// Step 4: Build the Initial Matrices:
//
// -----------------------------------------

    // Equation to use: S L_S = L_S Sigma_S
    // Diagonalize S
    Eigen::SelfAdjointEigenSolver<Matrix> diagS(S);
    Matrix L_S = diagS.eigenvectors();
    Matrix Sig_S = diagS.eigenvalues();

    // Raise the diagonal matrix Sig_S to the power of -1/2
    Matrix Sig_min_half (nbasis, nbasis);
    for(int i = 0; i < nbasis; i++){
        for(int j = 0; j < nbasis; j++){
            if (i == j){
                Sig_min_half(i,j) = pow(Sig_S(i,0), (-1.0/2.0));
            }else{
                Sig_min_half(i,j) = 0;
            }
        }
    }

    // Build the symmetric orthogonalization matrix S^(1/2)
    Matrix S_min_half = L_S * Sig_min_half * L_S.transpose();

    // We form an initial (guess) Fock matrix in the orthonormal AO basis
    // using the core Hamiltonian as a guess:
    Matrix F_prime0 = S_min_half.transpose() * H_core * S_min_half;

    // Diagonalize the fock matrix (F') to get the orthonormal C matrix (C') and
    // the energy matrix (epsilon_0)
    Eigen::SelfAdjointEigenSolver<Matrix> diagFprime(F_prime0);
    Matrix epsilon_0 = diagFprime.eigenvalues();
    Matrix C_prime0 = diagFprime.eigenvectors();

    // Transform the coefficients C' into the original (non-orthogonal) AO basis C:
    Matrix C = S_min_half * C_prime0;

    // Obtain a initial density matrix from the coefficients
    Matrix P_0 = computeDensity(nbasis, nElec, C);

    // E can be computed from P, Hcore and F
    // In our first iteration, our fock matrix F is just Hcore
    double E_0 = computeElecEnergy(P_0, H_core, H_core, nbasis) + Enuc;

// -----------------------------------------
//
// Step #5. First Density matrix propagation:
//
// -----------------------------------------


    //Set the initial parameters

    double time = 0.00;
    double dt = 0.01;
    double totalTime = 12.00;
    double curE = E_0;



    // We want to use MMUT to propagate the density matrix. To get P at t_(i+1), we use P at t_(i-1)
    // Obviously this is not doable for first propagation, so we use Trapezoidal rule

    // First, we get the field-dependent Fock matrix by adding the field contribution our our F'_0
    // *****************
    // TO IMPLEMENT:
    // Equation to use: F'(t) = F'_0(t) + d' e(t)
    // *****************
    // Matrix F_prime = F_prime0 + dipole * field;




    // Next, we want to obtain the unitary transformation matrix, U
    // *****************
    // TO IMPLEMENT:
    // Equation to use: C^dagger F'(t_i) C = epsilon(t_i)
    //                  U = C e^{i dt epsilon(t_i)} C^dagger
    // *****************

    /******** Pseudocode start:
    Eigen::SelfAdjointEigenSolver<Matrix> diagFprime(F_prime);
    Matrix epsilon = diagFprime.eigenvalues();
    Matrix C = diagFprime.eigenvectors();
    Matrix U = formU(Matrix C, Matrix epsilon, double dt, int nbasis);

    // I will probably use a function call to get U
     Matrix formU(Matrix C, Matrix epsilon, double dt, int nbasis) {

     U = Matrix(nbasis, nbasis);
     for(int i = 0; i < nbasis; i++){
         for(int j = 0; j < nbasis; j++){
            if (i == j){
                U(i,j) = e^(i * dt ^epsilon(i,0)) ??? How to multiply by imaginary unit?
            }else{
                U(i,j) = 0;
            }
        }
     }


    return  C * U * C^dagger;
}
    ********* Pseudocode end/



    // Then, we have our P for the next time step:
    // *****************
    // TO IMPLEMENT:
    // Equation to use: P(t_(i+1)) = U P_0(t_i) U^dagger
    // *****************


    /******** Pseudocode start:
    // Matrix nextP = U * P_0 * U.adjoint();
    ********* Pseudocode end/
       */



// -----------------------------------------
//
// Step #6. Propagating the future time steps
//
// -----------------------------------------


    Matrix preP;
    Matrix curP;
    Matrix nextP;
    while (time < totalTime){
        printf("\nAt t= %4.2f, the energy is %17.14f", time, curE);

        //Save the results of the previous run (at time step (i-1) )
        double preE = curE;
        Matrix preP = curP;
        Matrix curP = nextP;
        time += dt;


        // This gets us the field-free Fock matrix in the AO basis
        Matrix F0 = formFock(H_core, preP, ERI, nbasis);

        // Orthagonalize the field-free Fock matrix
        // Now F0 is in the orthonormal basis
        Matrix F_prime0 = S_min_half.transpose() * F0 * S_min_half;

        // Now add the field to the folk matrix
        // *****************
        // TO IMPLEMENT:
        // Equation to use: F'(t) = F'_0(t) + d' e(t)
        // *****************
        //
        // Matrix F_prime = F_prime0 + dipole * field

        // Next, we want to obtain the unitary transformation matrix, U
        // Here we use MMUT  ( 2 dt)
        // *****************
        // TO IMPLEMENT:
        // Equation to use: C^dagger F'(t_i) C = epsilon(t_i)
        //                  U = C e^{i 2dt epsilon(t_i)} C^dagger
        // *****************


        /******** Pseudocode start:
        Eigen::SelfAdjointEigenSolver<Matrix> diagFprime(F_prime);
        Matrix epsilon = diagFprime.eigenvalues();
        Matrix C = diagFprime.eigenvectors();
        Matrix U = formU(Matrix C, Matrix epsilon, double dt, int nbasis);
        ********* Pseudocode end/


        // Then, we can our P at (i+1):
        // *****************
        // TO IMPLEMENT:
        // Equation to use: P(t_(i+1)) = U P_0(t_(i-1)) U^dagger
        // *****************


        /******** Pseudocode start:
        nextP = U * preP * U.adjoint();
        curE = Enuc + computeElecEnergy(curP, H_core, F, nbasis);
        ********* Pseudocode end/
    }
         printf("\nAt t= %4.2f, the energy is %17.14f", time, curE);

    return 0;
}

