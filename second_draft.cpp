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
 * STO-3G Calculation for Water:
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
// Step 4: Build the Orthogonalization Matrix
//
// -----------------------------------------


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
    // Eigen has all kinds of functions that we can utilize
    Matrix S_min_half = L_S * Sig_min_half * L_S.transpose();
    //cout << S_min_half;


// -----------------------------------------
//
// Step 5: Build the Initial Guess Density
//
// -----------------------------------------


    // We form an initial (guess) Fock matrix in the orthonormal AO basis
    // using the core Hamiltonian as a guess:
    Matrix F_prime = S_min_half.transpose() * H_core * S_min_half;

    // Diagonalize the fock matrix (F') to get the orthonormal C matrix (C') and
    // the energy matrix (epsilon_0)
    Eigen::SelfAdjointEigenSolver<Matrix> diagFprime(F_prime);
    Matrix epsilon_0 = diagFprime.eigenvalues();
    Matrix C_prime = diagFprime.eigenvectors();

    // Transform the coefficients C' into the original (non-orthogonal) AO basis C:
    Matrix C = S_min_half * C_prime;

    // Obtain a initial density matrix from the coefficients
    Matrix P = computeDensity(nbasis, nElec, C);

// -----------------------------------------
//
// Step 6: Compute the Inital SCF Energy
//
// -----------------------------------------

    // E can be computed from P, Hcore and F
    // In our first iteration, our fock matrix F is just Hcore
    double elecE = computeElecEnergy(P, H_core, H_core, nbasis);

// -----------------------------------------
//
// Step #7,8,9: SCF procedure
//
// -----------------------------------------

    double dE = 1;
    double rmsP = 1;
    double dEthreshold = 0.000000000001 ;
    double rmsPthreshold = 0.000000000001 ;
    Matrix curP = P;
    double curE = elecE + Enuc;
    int nIter = 1;

    while (dE > dEthreshold || rmsP > rmsPthreshold){
        printf("\nAfter %2u iteration, Total E(SCF) = %17.14f", nIter, curE);

        // Save the previous energy value to test for convergence
        double preE = curE;

        // Save the previous P to test for convergence
        // And also, we need the previous P to form a new Fock matrix
        Matrix preP = curP;

        Matrix curF = formFock(H_core, preP, ERI, nbasis);

        // Orthagonalize the new Fock matrix
        Matrix curFprime = S_min_half.transpose() * curF * S_min_half;

        // Diagonalize F' to get C' and epsilon
        Eigen::SelfAdjointEigenSolver<Matrix> diagcurFprime(curFprime);
        Matrix curEpsilon_0 = diagcurFprime.eigenvalues();
        Matrix curCprime = diagcurFprime.eigenvectors();

        // Back Transform C' to C, we go from orthonormal basis to AO basis
        Matrix curC = S_min_half * curCprime;

        // Compute the new density matrix from C
        curP = computeDensity(nbasis, nElec, curC);

        rmsP = computeDensityDifference(curP, preP, nbasis);

        curE = Enuc + computeElecEnergy(curP, H_core, curF, nbasis);

        dE = curE - preE;

        nIter++;
    }

    printf("\nSCF converged after %2d iteration, Final Total E(SCF) = %17.14f", nIter, curE);

    return 0;
}

