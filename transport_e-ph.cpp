#include "atom.h"
#include "funct.h"

int main(){

   vector<Atom>          at_list;
   UNINT          natom;
   UNINT          n_bias;
   UNINT          n_move;
   UNINT          output_tot = 16;
   UNINT          nat2;
   vector<int>           move_at;
   vector<int>           tri_index;
   int                   tot_time;
   int                   print_t;
   double                delta_t;
   double                beta;
   double                elec_temp;
   double                phon_temp;
   double                Vbias; // 0.1 eV
   double                V_aux;
   double                c_coup;
   double                sigma;
   double                d_rate;
   double                mass;
   double                mean_freq;
   double                dr_aux ;
   double                aux1;
   double                aux2;
   complex<double>       exp_i;
   complex<double>       exp_j;
   double                currentA = 0.0, currentB = 0.0;
   ifstream              inputf;

   read_input(natom, n_bias, n_move, tot_time, print_t, delta_t, beta,
              elec_temp, phon_temp, Vbias, c_coup, sigma, d_rate,inputf,
              move_at, mass, mean_freq);

   nat2 = natom * natom;

   vector<double>        fock(nat2, 0.0);
   vector<double>        eigen_coef(nat2,0.0);
   vector<double>        eigen_coefT(nat2,0.0);
   vector<double>        aux_coef(nat2,0.0);
   vector<double>        aux_coefT(nat2,0.0);
   vector<double>        eigen_E(natom,0.0);
   vector<double>        eigen_aux(natom,0.0);
   vector<double>        ph_pop(natom,0.0);
   vector<double>        new_phon_pop(natom,0.0);
   vector<double>        Dphon_pop(natom,0.0);
   vector<double>        aux_array1(natom,0.0);
   vector<double>        eta_term(natom,0.0);
   vector<double>        lambda_term(natom,0.0);
   const double          pi = 3.141592653589793;
   vector < complex<double> > rho(nat2,0.0);
   vector < complex<double> > rho_OM(nat2,0.0);
   vector < complex<double> > rho_new(nat2,0.0);
   vector < complex<double> > Drho(nat2,0.0);
   vector < complex<double> > rho_ref(nat2,0.0);
   vector < complex<double> > aux_mat1(nat2,0.0);
   const complex<double> i_cmplx (0.0,1.0);
   ofstream              outfile[output_tot];

   init_output(outfile);
   atom_creator(at_list, n_bias, natom, n_move, move_at, mass, mean_freq);

//Preparing reference state/////////////////////////////////////////////////////
   hamiltonian_creator(fock, beta, at_list, natom, n_bias);
   apply_potential(fock, at_list, Vbias, n_bias, natom);
   eigenval_elec_calc(fock, eigen_E, eigen_coef, natom);
   for(int jj=0; jj < natom; jj++){
   for(int ii=0; ii < natom; ii++){
      eigen_coefT[jj + ii * natom] = eigen_coef[ii + jj * natom];
   }
   }
   warm_up_elec(rho_OM, eigen_E, natom, elec_temp);
   matmul(rho_OM, eigen_coefT, aux_mat1, natom);
   matmul(eigen_coef, aux_mat1, rho_ref, natom);

   double chargeA = 0.0;
   double chargeB = 0.0;

   for(int ii=0; ii<n_bias; ii++){
      chargeA += real(rho_ref[ii+ii*natom]) - 0.5;
      chargeB += real(rho_ref[(natom-1-ii)+(natom-1-ii)*natom]) - 0.5 ;
   }

   outfile[14]<<"A   "<< chargeA <<endl;
   outfile[14]<<"B   "<< chargeB <<endl;

////////////////////////////////////////////////////////////////////////////////

//Preparing initial state///////////////////////////////////////////////////////
   hamiltonian_creator(fock, beta, at_list, natom, n_bias);
   apply_potential(fock, at_list, Vbias, n_bias, natom);
   eigenval_elec_calc(fock, eigen_aux, aux_coef, natom);
   for(int jj=0; jj < natom; jj++){
   for(int ii=0; ii < natom; ii++){
      aux_coefT[jj + ii * natom] = aux_coef[ii + jj * natom];
   }
   }
   warm_up_elec(rho_OM, eigen_aux, natom, elec_temp);

   matmul(rho_OM, aux_coefT, aux_mat1, natom);
   matmul(aux_coef, aux_mat1, rho, natom);

//Kick over the density matrix
    // for (int ii = 1; ii < natom; ii++){
    // for (int jj = 1; jj < natom; jj++){
    //    aux1 =  cos(2*pi*d_rate*ii);
    //    aux2 = -sin(2*pi*d_rate*ii);
    //    exp_i=  aux1 + i_cmplx*aux2;
    //    aux1 =  cos(2*pi*d_rate*jj);
    //    aux2 =  sin(2*pi*d_rate*jj);
    //    exp_j=  aux1 + i_cmplx*aux2;
    //
    //    if( ii != jj ){
    //       rho[ii+jj*natom] = rho[ii+jj*natom] * exp_i * exp_j;
    //    }
    // }
    // }

   warm_up_ph(ph_pop, at_list, phon_temp, natom);

   write_output(natom, n_bias, fock, rho, rho_OM, eigen_E, 0.0, currentA,
                currentB, at_list, ph_pop, outfile);

   eliminating_negligible_terms(natom, at_list, tri_index, sigma, c_coup,
                                eigen_E, eigen_coef);

//##############################################################################
// Time propagation, using 2nd order runge kutta integrator.                   #
//##############################################################################
   for (int tt = 1; tt <= tot_time; tt++){

      dr_aux = d_rate;
      V_aux  = Vbias; //0.0e0;
      if (tt <= 1000){
      //    // dr_aux = d_rate * exp(-pow(((tt-1000)/250),2.0));
         V_aux     = 0.5 * Vbias * (1.0 + cos(pi/1000 * tt));
         apply_potential(fock, at_list, V_aux, n_bias, natom);
      }
      LVN_propagation(natom, fock, rho, Drho);
      electron_phonon_correction(natom, n_bias, at_list, tri_index, sigma, c_coup,
                                    eta_term, lambda_term, eigen_E, ph_pop,
                                    Dphon_pop, eigen_coef, eigen_coefT, rho,
                                    Drho);

      for(int ii=0; ii < nat2; ii++){
         aux_mat1[ii] = rho[ii] + 0.5 * Drho[ii] * delta_t;
      }
      for(int ii=0; ii < natom; ii++){
         aux_array1[ii] = ph_pop[ii] + 0.5 * Dphon_pop[ii] * delta_t;
         // aux_array1[ii] = ph_pop[ii];
      }


      if ((tt+0.5) <= 1000){
         V_aux     = 0.5 * Vbias * (1.0 + cos(pi/1000 * (tt+0.5)));
         apply_potential(fock, at_list, V_aux, n_bias, natom);
      }

      LVN_propagation(natom, fock, aux_mat1, Drho);
      electron_phonon_correction(natom,n_bias,at_list, tri_index, sigma, c_coup,
                                    eta_term, lambda_term, eigen_E, aux_array1,
                                    Dphon_pop, eigen_coef, eigen_coefT,
                                    aux_mat1, Drho);


      apply_Driving_term(rho, rho_ref, Drho, dr_aux, natom, n_bias,
                         currentA, currentB);

      for(int ii=0; ii < nat2; ii++){
         rho_new[ii] = rho[ii] + Drho[ii] * delta_t;
      }
      for(int ii=0; ii < natom; ii++){
         new_phon_pop[ii] = ph_pop[ii] + Dphon_pop[ii] * delta_t;
         // new_phon_pop[ii] = ph_pop[ii];
      }


      for(int ii=0; ii < nat2; ii++){rho[ii] = rho_new[ii];}
      for(int ii=0; ii < natom; ii++){ph_pop[ii] = new_phon_pop[ii];}

      if (tt%print_t == 0){
         double time = tt * delta_t;
         matmul(rho, eigen_coef, aux_mat1, natom);
         matmul(eigen_coefT, aux_mat1, rho_OM, natom);
         write_output(natom, n_bias, fock, rho, rho_OM, eigen_E, time, currentA,
                      currentB, at_list, ph_pop, outfile);
      }

   }

   for (int ii=0; ii<output_tot; ii++){outfile[ii].close();}

   return 0;
}
