#ifndef CLASS_Atom
#define CLASS_Atom

#include <iostream>
#include <vector>
#include <string>
#include <complex>
using namespace std;
typedef unsigned int UNINT;

class Atom
{
public:
   Atom(int ind, double mass, double Rcoor, double freq, double Hii,
        bool move);
   int    GetInd()   const;
   double GetMass()  const;
   double GetRCoor() const;
   double GetFreq()  const;
   double GetHii()   const;
   bool   GetMove()  const;
   void SetInd(int ind);
   void SetMass(double mass);
   void SetRCoor(double Rcoor);
   void SetFreq(double freq);
   void SetHii(double Hii);
   void SetMove(bool move);

private:
   int    at_ind;
   double at_mass;
   double at_Rcoor;
   double at_freq;
   double at_Hii;
   bool   at_move;
};

#endif
