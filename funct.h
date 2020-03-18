#ifndef SUBS_FOR_TRANSPORT
#define SUBS_FOR_TRANSPORT

#include <time.h>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "atom.h"
#include "matmul_cublas.h"

extern "C" { extern
   void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
                       double* w, double* work, int* lwork, int* info );
}

void init_output(ofstream* outfile);

void read_input(unsigned int& natom, unsigned int& n_bias, unsigned int& n_move,
                int& tot_time, int& print_t, double& delta_t, double& beta,
                double& elec_temp, double& phon_temp, double& Vbias,
                double& c_coup, double& sigma, double& d_rate, ifstream& inputf,
                vector<int>& move_at);

void write_output(unsigned int natom, unsigned int n_bias,vector<double>& fock,
                  vector<complex<double> >& rho_AT, vector<complex<double> >& rho_OM,
                  vector<double>& eigen_E, double time, double& currentA,
                  double& currentB, vector<Atom>& at_list, vector<double>& ph_pop,
                  ofstream* outfile);

void atom_creator(vector<Atom>& at_list, unsigned int n_bias,
                  unsigned int natom, unsigned int n_move,
                  vector<int>& move_at);

void hamiltonian_creator(vector<double>& fock, double beta,
                         vector<Atom>& at_list, unsigned int natom,
                         unsigned int n_bias);

void eigenval_elec_calc(vector<double>& mat, vector<double>& eigenval,
                        vector<double>& coef, unsigned int ntotal);

void warm_up_elec(vector<complex<double> >& rho, vector<double>& eigen_E,
                  unsigned int natom, double elec_temp);

void warm_up_ph(vector<double>& ph_pop, vector<Atom>& at_list, double ph_temp,
                unsigned int natom);

double Fop_j(unsigned int jj, unsigned int aa, unsigned int bb,
              unsigned int natom, vector<double>& coef);

double dirac_delta(double Ea, double Eb, double wj, double sigma);

void eliminating_negligible_terms(unsigned int natom, vector<Atom>& at_list,
                                  vector<int>& tri_index,
                                  double sigma, double c_coup,
                                  vector<double>& eigen_E,
                                  vector<double>& coef);

void eta_lambda_calc_phon_evol(unsigned int natom, vector<Atom>& at_list,
                               vector<int>& tri_index,
                               double sigma, double c_coup,
                               vector<double>& eta_term,
                               vector<double>& lambda_term,
                               vector<double>& eigen_E,
                               vector<double>& phon_pop,
                               vector<double>& Dphon_pop,
                               vector<double>& coef,
                               vector<complex<double> >& rho);

void matmul(vector<double>& matA, vector<double>& matB, vector<double>& matC,
            unsigned int dim);
void matmul(vector<complex<double> >& matA, vector<double>& matB,
            vector<complex<double> >& matC, unsigned int dim);
void matmul(vector<double>& matA, vector<complex<double> >& matB,
            vector<complex<double> >& matC, unsigned int dim);
void matmul_sparse(vector<double>& matA, vector<complex<double> >& matB,
                   vector<complex<double> >& matC, unsigned int dim);
void matmul_sparse(vector<complex<double> >& matA, vector<double>& matB,
                   vector<complex<double> >& matC, unsigned int dim);

void LVN_propagation(unsigned int natom, vector<double>& fock,
                     vector<complex<double> >& rho,
                     vector<complex<double> >& Drho);

void apply_potential(vector<double>& fock, vector<Atom>& at_list, double Vpot,
                     unsigned int n_bias, unsigned int natom);

void apply_Driving_term(vector<complex<double> >& rho,
                        vector<complex<double> >& rho_ref,
                        vector<complex<double> >& Drho, double d_rate,
                        unsigned int natom, unsigned int n_bias,
                        double& currentA, double& currentB);

void electron_phonon_correction(unsigned int natom, unsigned int n_bias,
                                vector<Atom>& at_list,
                                vector<int>& tri_index,
                                double sigma, double c_coup,
                                vector<double>& eta_term,
                                vector<double>& lambda_term,
                                vector<double>& eigen_E,
                                vector<double>& phon_pop,
                                vector<double>& Dphon_pop,
                                vector<double>& coef,
                                vector<double>& coefT,
                                vector<complex<double> >& rho,
                                vector<complex<double> >& Drho);

#endif
