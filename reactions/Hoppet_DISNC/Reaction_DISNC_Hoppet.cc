
/*
   @file Reaction_DISNC_Hoppet.cc
   @date 2017-04-08
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-04-08
*/

#include "Reaction_DISNC_Hoppet.h"
#include <iostream>
#include <cstdio>


// the class factories
extern "C" Reaction_DISNC_Hoppet *create()
{
  return new Reaction_DISNC_Hoppet();
}

// Initialize at the start of the computation
void Reaction_DISNC_Hoppet::atStart()
{
}

// Main function to compute results at an iteration
void Reaction_DISNC_Hoppet::compute(TermData *td, valarray<double> &valExternal, map<string, valarray<double>> &errExternal)
{
}

void Reaction_DISNC_Hoppet::atIteration()
{
  // Make sure to call the parent class initialization:
  super::atIteration();
}

//
void Reaction_DISNC_Hoppet::initTerm(TermData *td)
{
  super::initTerm(td);
};

void Reaction_DISNC_Hoppet::F2gamma BASE_PARS
{
  valarray<double> f2u, f2d;
}

void Reaction_DISNC_Hoppet::F2gammaZ BASE_PARS
{
  valarray<double> f2u, f2d;
}

void Reaction_DISNC_Hoppet::F2Z BASE_PARS
{
  valarray<double> f2u, f2d;
}

void Reaction_DISNC_Hoppet::FLgamma BASE_PARS
{
  valarray<double> flu, fld;
}

void Reaction_DISNC_Hoppet::FLgammaZ BASE_PARS
{
  valarray<double> flu, fld;
}

void Reaction_DISNC_Hoppet::FLZ BASE_PARS
{
  valarray<double> flu, fld;
}

void Reaction_DISNC_Hoppet::xF3gammaZ BASE_PARS
{
  valarray<double> xf3u, xf3d;
}

void Reaction_DISNC_Hoppet::xF3Z BASE_PARS
{
  valarray<double> xf3u, xf3d;
}


