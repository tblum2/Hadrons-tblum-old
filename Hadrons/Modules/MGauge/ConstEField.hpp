/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MGauge/ConstEField.hpp

Copyright (C) 2015-2019

Author: Tom Blum

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#ifndef Hadrons_MGauge_ConstEField_hpp_
#define Hadrons_MGauge_ConstEField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                              ConstEField gauge                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

/****************************************************************************
*  ConstEField times a gauge field:
*
*  Ue_mu(x) = U_mu(x)*exp(ieqA_mu(x))
*
*  with
*
*  - gauge: U_mu(x): gauge field
*  - emField: A_mu(x): electromagnetic photon field
*  - e: value for the elementary charge
*  - q: charge in units of e
*
*****************************************************************************/


class ConstEFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ConstEFieldPar,
                                    std::string, gauge,
                                    double, charge);
};

template <typename GImpl>
class TConstEField: public Module<ConstEFieldPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    typedef PhotonR::GaugeField     EmField;
    // constructor
    TConstEField(const std::string name);
    // destructor
    virtual ~TConstEField(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ConstEField, TConstEField<GIMPL>, MGauge);

/******************************************************************************
*                            TConstEField implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TConstEField<GImpl>::TConstEField(const std::string name)
: Module<ConstEFieldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TConstEField<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};

    return in;
}

template <typename GImpl>
std::vector<std::string> TConstEField<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TConstEField<GImpl>::setup(void)
{
    envCreateLat(GaugeField, getName());
    envTmpLat(LatticeComplex, "eiAmu");
    envTmpLat(LatticeComplex, "Amu");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TConstEField<GImpl>::execute(void)
{
    LOG(Message) << "ConstEField times the gauge field " << par().gauge << " with charge q= "  << par().charge << std::endl;
    
    auto &Ue = envGet(GaugeField, getName());
    auto &U = envGet(GaugeField, par().gauge);

    envGetTmp(LatticeComplex, eiAmu);
    //envGetTmp(LatticeComplex, Amu);

    Complex i(0.0,1.0);

    // constant E field in z dir (arXive hep-lat/1001.1131)
    int Nt = env().getDim(Tp);
    int Ns = env().getDim(Xp);
    RealD Elatt = 2*EIGEN_PI/(0.33333333)/Real(Nt*Ns);
    LatticeComplex z(env().getGrid()); LatticeCoordinate(z,2);
    LatticeComplex t(env().getGrid()); LatticeCoordinate(t,3);
    PokeIndex<LorentzIndex>(Ue, PeekIndex<LorentzIndex>(U, 0), 0);
    PokeIndex<LorentzIndex>(Ue, PeekIndex<LorentzIndex>(U, 1), 1);
    // linear in t potential in z dir
    LatticeComplex Amu(env().getGrid());
    Amu = -par().charge * Elatt * t;
    eiAmu = exp(i * Amu);
    PokeIndex<LorentzIndex>(Ue, PeekIndex<LorentzIndex>(U, 2) * eiAmu, 2);
    // boundary terms for t component of potential
    Lattice<iScalar<vInteger>> ti(env().getGrid()); LatticeCoordinate(ti,3);
    LatticeComplex zero(env().getGrid());
    zero = Zero();
    z = where(ti==(Integer)(Nt-1), z, zero);
    eiAmu = exp(i * par().charge * (double)Nt * Elatt * z);
    PokeIndex<LorentzIndex>(Ue, PeekIndex<LorentzIndex>(U, 3) * eiAmu, 3);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_ConstEField_hpp_
