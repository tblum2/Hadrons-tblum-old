/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MAction/NaiveStaggered.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>

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

#ifndef Hadrons_MAction_NaiveStaggered_hpp_
#define Hadrons_MAction_NaiveStaggered_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            TNaiveStaggered quark action                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class NaiveStaggeredPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(NaiveStaggeredPar,
                                    std::string, gauge,
                                    double     , mass,
                                    double     , c1,
                                    double     , tad,
                                    std::string, boundary,
                                    std::string, string,
                                    std::string, twist);
};

template <typename FImpl>
class TNaiveStaggered: public Module<NaiveStaggeredPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TNaiveStaggered(const std::string name);
    // destructor
    virtual ~TNaiveStaggered(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(NaiveStaggered, TNaiveStaggered<STAGIMPL>, MAction);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(NaiveStaggeredF, TNaiveStaggered<STAGIMPLF>, MAction);
#endif

/******************************************************************************
 *                     TNaiveStaggered template implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TNaiveStaggered<FImpl>::TNaiveStaggered(const std::string name)
: Module<NaiveStaggeredPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TNaiveStaggered<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge };
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TNaiveStaggered<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TNaiveStaggered<FImpl>::setup(void)
{
    LOG(Message) << "Setting up NaiveStaggered fermion matrix with m= " << par().mass
                 << " using gauge field '" << par().gauge << "'" << std::endl;
                 
    auto &U      = envGet(GaugeField, par().gauge);
    auto &grid   = *envGetGrid(FermionField);
    auto &gridRb = *envGetRbGrid(FermionField);
    typename NaiveStaggeredFermion<FImpl>::ImplParams implParams;
    
    envCreateDerived(FMat, NaiveStaggeredFermion<FImpl>, getName(), 1,
                     U,
                     grid, gridRb,
                     par().mass, par().c1, par().tad, implParams);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TNaiveStaggered<FImpl>::execute()
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_NaiveStaggered_hpp_
