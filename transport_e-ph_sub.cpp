#include "funct.h"

//##############################################################################
void init_output(ofstream* outfile){
      outfile[0].open("elec_energy.out");
      outfile[1].open("total_elec_density.out");
      outfile[2].open("elec_eigen_val.out");
      outfile[3].open("elec_pop.out");
      outfile[4].open("time.out");
      outfile[5].open("chargeA.out");
      outfile[6].open("chargeB.out");
      outfile[7].open("chargeM.out");
      outfile[8].open("currentA.out");
      outfile[9].open("currentB.out");
      outfile[10].open("phon_energy.out");
      outfile[11].open("total_phon_density.out");
      outfile[12].open("phon_eigen_val.out");
      outfile[13].open("phon_pop.out");
      outfile[14].open("ref_charge.out");
}
//##############################################################################
void read_input(UNINT& natom, UNINT& n_bias, UNINT& n_move,
                int& tot_time, int& print_t, double& delta_t, double& beta,
                double& elec_temp, double& phon_temp, double& Vbias,
                double& c_coup, double& sigma, double& d_rate, ifstream& inputf,
                vector<int>& move_at){

   inputf.open("input");

   if (!inputf) {
      cout << "Unable to open file";
      exit(1); // terminate with error
    }

    inputf>> natom;
    inputf>> n_bias;
    inputf>> n_move;

    for (int ii=0; ii < n_move; ii++){
      move_at.push_back(0);
      inputf>> move_at[ii];
   }

    inputf>> tot_time;
    inputf>> print_t;
    inputf>> delta_t;
    inputf>> beta;
    inputf>> elec_temp;
    inputf>> phon_temp;
    inputf>> Vbias;
    inputf>> c_coup;
    inputf>> sigma;
    inputf>> d_rate;

    Vbias = Vbias * 0.0367493;

    inputf.close();

}

//##############################################################################
void write_output(UNINT natom, UNINT n_bias,vector<double>& fock,
                  vector<complex<double> >& rho_AT, vector<complex<double> >& rho_OM,
                  vector<double>& eigen_E, double time, double& currentA,
                  double& currentB, vector<Atom>& at_list, vector<double>& ph_pop,
                  ofstream* outfile){

   if(time <= 0.0e0){
      for(int ii=0; ii<natom; ii++){
         outfile[2] << eigen_E[ii] << endl;
         outfile[12] << at_list[ii].GetFreq() << endl;
      }
   }

   UNINT nat2 = natom * natom;
   int    n_mol      = natom - 2*n_bias;
   double tot_elec_E = 0.0;
   double tot_elec   = 0.0;
   double tot_phon_E = 0.0;
   double tot_phon   = 0.0;
   double chargeA    = 0.0;
   double chargeB    = 0.0;
   double chargeM    = 0.0;
   vector<complex<double> > aux_mat1(nat2, 0.0);
   vector<complex<double> > aux_mat2(nat2, 0.0);

   matmul(rho_AT, fock, aux_mat1, natom);

   for(int ii=0; ii<natom; ii++){
      tot_elec_E += real(aux_mat1[ii+ii*natom]);
      tot_elec   += real(rho_OM[ii+ii*natom]);
      tot_phon_E += (ph_pop[ii] + 0.5) * at_list[ii].GetFreq();
      tot_phon   += ph_pop[ii];
   }

   for(int ii=0; ii<n_bias; ii++){
      chargeA += -real(rho_AT[ii+ii*natom]) + 0.5;
      chargeB += -real(rho_AT[(natom-1-ii)+(natom-1-ii)*natom]) + 0.5 ;
   }
   for(int ii=n_bias; ii<n_bias+n_mol; ii++){
      chargeM += real(rho_AT[ii+ii*natom]) - 0.5 ;
   }

   outfile[0]<< tot_elec_E <<endl;
   outfile[1]<< tot_elec <<endl;
   outfile[4]<< time <<endl;
   outfile[5]<< chargeA <<endl;
   outfile[6]<< chargeB <<endl;
   outfile[7]<< chargeM <<endl;
   outfile[8]<< 2*currentA*1.6e-19*1e6/(2.419e-17) <<endl;
   outfile[9]<< 2*currentB*1.6e-19*1e6/(2.419e-17) <<endl;
   outfile[10]<< tot_phon_E <<endl;
   outfile[11]<< tot_phon <<endl;

   for(int ii=0; ii<natom; ii++){
      outfile[3]<<real(rho_OM[ii+ii*natom])<<endl;
      outfile[13]<<ph_pop[ii]<<endl;
   }

}

//##############################################################################

void atom_creator(vector<Atom>& at_list, UNINT n_bias,
                  UNINT natom, UNINT n_move,
                  vector<int>& move_at){

   double mass  = 1822.8885*0.5;//23244.3;
   double Rcoor = -3.0;
   double w_min = 0.00367493 , w_max = 0.0110248;
   double Hii   = 0.0e0;
   double freq;
   bool   move;
   int    jj;


   srand(time(NULL));

   for (int ii = 0; ii<natom; ii++){
      double delta  =  w_max - w_min ;
      double n_rand = ((double) rand() / (RAND_MAX));

      // freq = n_rand * delta + w_min ;

      move = false;

      jj = 0;
      while(!move && jj < n_move){
         move = (ii == move_at[jj]);
         jj += 1;
      }

      freq = 0.0;
      if(move){freq = 0.00734987;}//n_rand * delta + w_min;}//0.004409919;

      Atom new_atom(ii, mass, Rcoor * (ii+1), freq, Hii, move);
      at_list.push_back(new_atom);
   }
}
//##############################################################################

void hamiltonian_creator(vector<double>& fock, double beta,
                         vector<Atom>& at_list, UNINT natom,
                         UNINT n_bias){

   double gamma = 0.012249775;//-0.0183747;

   UNINT nat2 = natom * natom;
   for(int ii=0; ii < nat2; ii++){fock[ii]=0.0e0;}

   for(int ii=0; ii < natom; ii++){
      fock[ii + ii*natom] = at_list[ii].GetHii();
      if(ii != natom-1){
         fock[ii + (ii+1)*natom] = beta;
         fock[(ii+1) + ii*natom] = beta;
      }
   }

   // int indx = n_bias;
   // fock[indx-1 + indx*natom]   = gamma;
   // fock[indx + (indx-1)*natom] = gamma;
   // indx = natom - n_bias;
   // fock[indx-1 + indx*natom]   = gamma;
   // fock[indx + (indx-1)*natom] = gamma;

}
//##############################################################################

void eigenval_elec_calc(vector<double>& mat, vector<double>& eigenval,
                        vector<double>& coef, UNINT ntotal){

   int          info, lwork;
   int          dim = (int) ntotal;
   UNINT n2  = ntotal * ntotal;
   double       wkopt;
   double*      work;
   char         jobz='V';
   char         uplo='U';

   for(int ii=0; ii < n2; ii++){coef[ii]=mat[ii];}

   lwork = -1;
   dsyev_( &jobz,&uplo, &dim, & *coef.begin(), &dim, & *eigenval.begin(),
           &wkopt, &lwork, &info);
   lwork = (int)wkopt;
   work = (double*)malloc( lwork*sizeof(double) );
   dsyev_( &jobz,&uplo, &dim, & *coef.begin(), &dim, & *eigenval.begin(), work,
           &lwork, &info );

}
//##############################################################################
void warm_up_elec(vector<complex<double> >& rho, vector<double>& eigen_E,
                  UNINT natom, double elec_temp){

   for (int ii=0; ii < natom; ii++){
      rho[ii+ii*natom] = 1.0e0/(exp(eigen_E[ii]/elec_temp) + 1.0e0);
   }
}
//##############################################################################
void warm_up_ph(vector<double>& ph_pop, vector<Atom>& at_list, double ph_temp,
                UNINT natom){

   for (int ii=0; ii < natom; ii++){
      if (at_list[ii].GetMove()){
         ph_pop[ii] = 1.0e0/(exp(at_list[ii].GetFreq()/ph_temp) - 1.0e0);
      }
   }
}
//##############################################################################
double Fop_j(UNINT jj, UNINT aa, UNINT bb,
              UNINT natom, vector<double>& coef){

   double l_term;
   double r_term;
   double Fj;

   l_term = coef[jj+aa*natom] * coef[(jj-1)+bb*natom] +
            coef[(jj-1)+aa*natom] * coef[jj+bb*natom];
   r_term = -coef[(jj+1)+aa*natom] * coef[jj+bb*natom] -
            coef[jj+aa*natom] * coef[(jj+1)+bb*natom];

   Fj = pow(l_term + r_term, 2.0e0);

   return Fj;
}
//##############################################################################

double dirac_delta(double Ea, double Eb, double wj, double sigma){

   double norm;
   double arg;
   double dir_del;
   const double pi = 3.141592653589793;

   norm = 1.0 / sqrt(2.0 * pi * pow(sigma, 2.0));
   arg  = -pow((Ea - Eb + wj), 2.0)/(2.0 * pow(sigma, 2.0));

   dir_del = norm * exp(arg);

   return dir_del;
}
//##############################################################################

void eliminating_negligible_terms(UNINT natom, vector<Atom>& at_list,
                                  vector<int>& tri_index,
                                  double sigma, double c_coup,
                                  vector<double>& eigen_E,
                                  vector<double>& coef){
   int nat2 = natom*natom;

   for (int ii=0; ii<natom; ii++){
   for (int jj=0; jj<natom; jj++){
   for (int kk=0; kk<natom; kk++){
      if(at_list[kk].GetMove()){
         double aux1, exp1, exp2;
         double Ea   = eigen_E[ii];
         double Eb   = eigen_E[jj];
         double wj   = at_list[kk].GetFreq();
         aux1 = pow(c_coup, 2.0) * Fop_j(kk, ii, jj, natom, coef);
         exp1 = dirac_delta(Ea, Eb, wj, sigma);
         exp2 = dirac_delta(Ea, Eb, -wj, sigma);
         bool acceptable;
         acceptable = (aux1*exp1 > 1.0e-6) || (aux1*exp2 > 1.0e-6);
         if (acceptable){
            tri_index.push_back(ii+jj*natom+kk*nat2);
         }
      }
   }
   }
   }
}

//##############################################################################
void eta_lambda_calc_phon_evol(UNINT natom, vector<Atom>& at_list,
                               vector<int>& tri_index,
                               double sigma, double c_coup,
                               vector<double>& eta_term,
                               vector<double>& lambda_term,
                               vector<double>& eigen_E,
                               vector<double>& phon_pop,
                               vector<double>& Dphon_pop,
                               vector<double>& coef,
                               vector<complex<double> >& rho){



   int n_iter = tri_index.size();
   int nat2   = natom * natom;
   const double pi = 3.141592653589793;

   for (int ii=0; ii < natom; ii++){
      eta_term[ii]    = 0.0;
      lambda_term[ii] = 0.0;
      Dphon_pop[ii]   = 0.0;
   }

   for (int indx=0; indx < n_iter; indx++){

      int kk = int(tri_index[indx]/nat2);
      int jj = int((tri_index[indx]-kk*nat2)/natom);
      int ii = int(tri_index[indx]-jj*natom-kk*nat2);

      if (ii != jj){
         double aux1, exp1, exp2;
         double Ea   = eigen_E[ii];
         double Eb   = eigen_E[jj];
         double wj   = at_list[kk].GetFreq();
         double fa   = real(rho[ii+ii*natom]);
         double fb   = real(rho[jj+jj*natom]);
         double Nj   = phon_pop[kk];
         double mass = at_list[kk].GetMass();

         aux1 = pow(c_coup, 2.0) * pi *
                Fop_j(kk, ii, jj, natom, coef)/(mass * wj);
         exp1 = dirac_delta(Ea, Eb, wj, sigma);
         exp2 = dirac_delta(Ea, Eb, -wj, sigma);

         eta_term[ii] += aux1 * ((Nj + fb) * exp1 +
                         (Nj + 1.0 - fb) * exp2);
         lambda_term[ii] += aux1 * fb * ((Nj + 1.0) * exp1 + Nj * exp2);
         Dphon_pop[kk] += aux1 * (-fa * (1.0 - fb)* Nj +
                             fb * (1.0 - fa) * (Nj + 1.0)) * (Eb-Ea)/wj
                             * exp1;
      }
   }
}
//##############################################################################
void matmul(vector<double>& matA, vector<double>& matB, vector<double>& matC,
            UNINT dim){

   int            ii, jj, kk;
   int            dim2 = dim * dim;
   vector<double> aux_mat(dim2, 0.0);

   for(ii=0; ii < dim2; ii++){matC[ii] = 0.0;}
   for(ii=0; ii < dim; ii++){
   for(jj=0; jj < dim; jj++){
      aux_mat[ii+jj*dim] = matA[jj+ii*dim];
   }
   }

   for(ii=0; ii < dim; ii++){
   for(jj=0; jj < dim; jj++){
   for(kk=0; kk < dim; kk++){
      matC[ii+jj*dim] += matA[kk+ii*dim] * matB[kk+jj*dim];
   }
   }
   }
}

void matmul(vector<complex<double> >& matA, vector<double>& matB,
            vector<complex<double> >& matC, UNINT dim){

   int            ii;
   int                        dim2 = dim * dim;
   vector< complex<double> >  aux_mat(dim2, 0.0);

   for(ii=0; ii < dim2; ii++){
      matC[ii] = 0.0;
      aux_mat[ii]=matB[ii];
   }

   matcublas(& *matA.begin(), & *aux_mat.begin(), & *matC.begin(), dim);

}

void matmul(vector<double>& matA, vector<complex<double> >& matB,
            vector<complex<double> >& matC, UNINT dim){

   int ii;
   int dim2 = dim * dim;
   vector< complex<double> >  aux_mat(dim2, 0.0);

   for(ii=0; ii < dim2; ii++){
      matC[ii] = 0.0;
      aux_mat[ii]=matA[ii];
   }

   matcublas(& *aux_mat.begin(), & *matB.begin(), & *matC.begin(), dim);

}

void matmul_sparse(vector<double>& matA, vector<complex<double> >& matB,
                   vector<complex<double> >& matC, UNINT dim){

   for(int jj=0; jj < dim; jj++){
   for(int ii=0; ii < dim; ii++){
      matC[ii+jj*dim] = 0.0;
      int kk_min;
      int kk_max;
      if (ii == 0){
         kk_min = 0; kk_max = 2;
      }
      else if(ii == dim-1){
         kk_min = dim - 2; kk_max = dim;
      }
      else{
         kk_min = ii - 1; kk_max = ii + 2;
      }
      for(int kk=kk_min; kk < kk_max; kk++){
         matC[ii+jj*dim] += matA[ii+kk*dim] * matB[kk+jj*dim];
      }
   }
   }

}

void matmul_sparse(vector<complex<double> >& matA, vector<double>& matB,
                   vector<complex<double> >& matC, UNINT dim){

   for(int jj=0; jj < dim; jj++){
   for(int ii=0; ii < dim; ii++){
      matC[ii+jj*dim] = 0.0;
      int kk_min;
      int kk_max;
      if (jj == 0){
         kk_min = 0; kk_max = 2;
      }
      else if(jj == dim-1){
         kk_min = dim - 2; kk_max = dim;
      }
      else{
         kk_min = jj - 1; kk_max = jj + 2;
      }
      for(int kk=kk_min; kk < kk_max; kk++){
         matC[ii+jj*dim] += matA[ii+kk*dim] * matB[kk+jj*dim];
      }
   }
   }

}
//##############################################################################

void LVN_propagation(UNINT natom, vector<double>& fock,
                     vector<complex<double> >& rho,
                     vector<complex<double> >& Drho){

   UNINT              nat2 = natom * natom;
   const complex<double>     i_cmplx (0.0, 1.0);
   vector<complex<double> >  aux1(nat2, 0.0);
   vector<complex<double> >  aux2(nat2, 0.0);
   vector<complex<double> >  fp_pf(nat2, 0.0);

   matmul_sparse(fock, rho, aux1, natom);
   matmul_sparse(rho, fock, aux2, natom);

   for (int ii=0; ii < nat2; ii++){fp_pf[ii]=aux1[ii]-aux2[ii];}
   for (int ii=0; ii < nat2; ii++){Drho[ii] = -i_cmplx * fp_pf[ii];}
}

//##############################################################################
void apply_potential(vector<double>& fock, vector<Atom>& at_list, double Vpot,
                     UNINT n_bias, UNINT natom){

   for (int ii=0; ii < n_bias; ii++){
      double Ha = at_list[ii].GetHii();
      double Hb = at_list[natom-(ii+1)].GetHii();
      fock[ii+ii*natom] = Ha + Vpot/2.0;
      fock[(natom-(ii+1))+(natom-(ii+1))*natom] = Hb - Vpot/2.0;
   }
}

//##############################################################################
void apply_Driving_term(vector<complex<double> >& rho,
                        vector<complex<double> >& rho_ref,
                        vector<complex<double> >& Drho, double d_rate,
                        UNINT natom, UNINT n_bias,
                        double& currentA, double& currentB){

   int n_mol  = natom - 2*n_bias;
   currentA   = 0.0;
   currentB   = 0.0;

   for(int jj=0; jj< n_bias; jj++){
   for(int ii=0; ii< n_bias; ii++){
      int index1 = ii + jj*natom;
      int index2 = natom - 1 - ii + (natom - 1 - jj) * natom;
      Drho[index1] += -d_rate * (rho[index1] - rho_ref[index1]);
      Drho[index2] += -d_rate * (rho[index2] - rho_ref[index2]);
      if (ii==jj){
         currentA += d_rate * real(rho[index1] - rho_ref[index1]);
         currentB += d_rate * real(rho[index2] - rho_ref[index2]);
      }
   }
   }
   for(int jj=n_bias; jj< (n_bias + n_mol); jj++){
   for(int ii=0; ii< n_bias; ii++){
      int index1 = ii + jj*natom;
      int index2 = jj + ii*natom;
      Drho[index1] += -d_rate * 0.5 * (rho[index1] - rho_ref[index1]);
      Drho[index2] += -d_rate * 0.5 * (rho[index2] - rho_ref[index2]);
   }
   }
   for(int jj=(n_bias + n_mol); jj< natom; jj++){
   for(int ii=n_bias; ii< (n_bias + n_mol); ii++){
      int index1 = ii + jj*natom;
      int index2 = jj + ii*natom;
      Drho[index1] += -d_rate * 0.5 * (rho[index1]  - rho_ref[index1]);
      Drho[index2] += -d_rate * 0.5 * (rho[index2]  - rho_ref[index2]);
   }
   }
   for(int jj=(n_bias + n_mol); jj< natom; jj++){
   for(int ii=0; ii < n_bias; ii++){
      int index1 = ii + jj*natom;
      int index2 = jj + ii*natom;
      Drho[index1] += -d_rate * (rho[index1] - rho_ref[index1]);
      Drho[index2] += -d_rate * (rho[index2] - rho_ref[index2]);
   }
   }
}

//##############################################################################
void electron_phonon_correction(UNINT natom, vector<Atom>& at_list,
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
                                vector<complex<double> >& Drho){

   UNINT    nat2 = natom * natom;
   vector<complex<double> > aux_mat1(nat2, 0.0);
   vector<complex<double> > aux_mat2(nat2, 0.0);
   vector<complex<double> > rho_OM(nat2, 0.0);
   vector<complex<double> > Drho_OM(nat2, 0.0);

   matmul(rho, coef, aux_mat1, natom);
   matmul(coefT, aux_mat1, rho_OM, natom);

   eta_lambda_calc_phon_evol(natom, at_list, tri_index, sigma, c_coup, eta_term,
                             lambda_term, eigen_E, phon_pop, Dphon_pop,
                             coef, rho_OM);

   for (int jj=0; jj < natom; jj++){
   for (int ii=0; ii < natom; ii++){
      if (ii == jj){
         Drho_OM[ii + jj*natom] = - eta_term[ii] * rho_OM[ii + jj*natom]
                                + lambda_term[ii];
      }
      else{
         Drho_OM[ii + jj*natom] = -0.5*(eta_term[ii]+eta_term[jj])
                                   * rho_OM[ii + jj*natom];
         // Drho_OM[ii + jj*natom] = -sqrt(eta_term[ii]*eta_term[jj])
                                  // * rho_OM[ii + jj*natom];
      }
   }
   }

   matmul(Drho_OM, coefT, aux_mat1, natom);
   matmul(coef, aux_mat1, aux_mat2, natom);

   // for (int ii=0; ii < natom; ii++){
   //    cout<<real(aux_mat2[ii+ii*natom])<<endl;
   // }
   //
   // cout<<"________________"<<endl;

   for (int ii=0; ii < nat2; ii++){Drho[ii] += aux_mat2[ii];}

}
//##############################################################################
