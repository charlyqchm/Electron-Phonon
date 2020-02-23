#include "atom.h"

Atom::Atom(int ind, double mass, double Rcoor, double freq, double Hii,
           bool move){

   at_ind   = ind;
   at_mass  = mass;
   at_Rcoor = Rcoor;
   at_freq  = freq;
   at_Hii   = Hii;
   at_move  = move;
}

int Atom::GetInd()const{
   return at_ind;
}

double Atom::GetMass()const{
   return at_mass;
}

double Atom::GetRCoor()const{
   return at_Rcoor;
}

double Atom::GetFreq()const{
   return at_freq;
}

double Atom::GetHii()const{
   return at_Hii;
}

bool Atom::GetMove()const{
   return at_move;
}

void Atom::SetInd(int ind){
   at_ind = ind;
}

void Atom::SetMass(double mass){
   at_mass = mass;
}

void Atom::SetRCoor(double Rcoor){
   at_Rcoor = Rcoor;
}

void Atom::SetFreq(double freq){
   at_freq = freq;
}

void Atom::SetHii(double Hii){
   at_Hii = Hii;
}

void Atom::SetMove(bool move){
   at_move = move;
}
