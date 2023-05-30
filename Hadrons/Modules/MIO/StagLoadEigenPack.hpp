/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MIO/StagLoadEigenPack.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>

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
#ifndef Hadrons_MIO_StagLoadEigenPack_hpp_
#define Hadrons_MIO_StagLoadEigenPack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Solver.hpp> //needed for FMat

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Load eigen vectors/values package                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class StagLoadEigenPackPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagLoadEigenPackPar,
                                    std::string, filestem,
                                    std::string, action,
                                    bool, multiFile,
                                    bool, doubleMemory,
                                    bool, even,
                                    unsigned int, size,
                                    double, mass,
                                    unsigned int, Ls,
                                    std::string, gaugeXform);
};

template <typename Pack, typename GImpl, typename FImpl>
class TStagLoadEigenPack: public Module<StagLoadEigenPackPar>
{
public:
    typedef FermionOperator<FImpl> FMat;
    typedef typename Pack::Field   Field;
    typedef typename Pack::FieldIo FieldIo;
    typedef BaseEigenPack<Field>   BasePack;

public:
    GAUGE_TYPE_ALIASES(GImpl, );
    typedef typename GImpl::GaugeLinkField GaugeMat;
public:
    // constructor
    TStagLoadEigenPack(const std::string name);
    // destructor
    virtual ~TStagLoadEigenPack(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(StagLoadFermionEigenPack, ARG(TStagLoadEigenPack<FermionEigenPack<STAGIMPL>, GIMPL, STAGIMPL>), MIO);

/******************************************************************************
 *                    TStagLoadEigenPack implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl, typename FImpl>
TStagLoadEigenPack<Pack, GImpl, FImpl>::TStagLoadEigenPack(const std::string name)
: Module<StagLoadEigenPackPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack, typename GImpl, typename FImpl>
std::vector<std::string> TStagLoadEigenPack<Pack, GImpl, FImpl>::getInput(void)
{
    std::vector<std::string> in;

    if (!par().gaugeXform.empty())
    {
        in = {par().gaugeXform};
    }
    in.push_back(par().action);
    
    return in;
}

template <typename Pack, typename GImpl, typename FImpl>
std::vector<std::string> TStagLoadEigenPack<Pack, GImpl, FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl, typename FImpl>
void TStagLoadEigenPack<Pack, GImpl, FImpl>::setup(void)
{
    GridBase *gridIo = nullptr, *grid = nullptr;

    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        if(par().Ls > 1){
            gridIo = envGetRbGrid(FieldIo, par().Ls);
        }else{
            gridIo = envGetRbGrid(FieldIo);
        }
    }
    if(par().doubleMemory){
        gridIo = envGetRbGrid(Field);
        grid = envGetGrid(Field);//double storage so we can overwrite with v vector
    }else{
        grid = envGetRbGrid(Field);
    }
    envCreateDerived(BasePack, Pack, getName(),
                     par().Ls, par().size,
                     grid, gridIo);
    
    if (!par().gaugeXform.empty())
    {
        if (par().Ls > 1)
        {
            LOG(Message) << "Setup 5d GaugeMat for Ls = " << par().Ls << std::endl;
            envTmp(GaugeMat,    "tmpXform", par().Ls, envGetGrid5(Field, par().Ls));
            envTmp(GaugeMat, "tmpXformOdd", par().Ls, envGetRbGrid5(Field, par().Ls));
        }
        else
        {
            LOG(Message) << "Setup 4d GaugeMat for Ls = " << par().Ls << std::endl;
            envTmp(GaugeMat,    "tmpXform", par().Ls, envGetGrid(Field));
            envTmp(GaugeMat, "tmpXformOdd", par().Ls, envGetRbGrid(Field));
        }
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl, typename FImpl>
void TStagLoadEigenPack<Pack, GImpl, FImpl>::execute(void)
{
    
    auto &action= envGet(FMat, par().action);
    auto &epack = envGetDerived(BasePack, Pack, getName());

    epack.read(par().filestem, par().multiFile, vm().getTrajectory());
    epack.eval.resize(par().size);
    
    ComplexD minusI(0, -1.0);
    ComplexD cc;
    RealD eval;
    double mass = par().mass;
    Field temp(epack.evec[0].Grid());
    
    // make even evecs from Odd
    if(par().even && !par().doubleMemory){
        for (unsigned int i = 0; i < par().size; i++)
        {
            eval=sqrt(epack.eval[i]-mass*mass);
            epack.evec[i].Checkerboard() = Odd;
            action.Meooe(epack.evec[i], temp);
            cc = minusI/eval;
            epack.evec[i] = cc * temp; // now it's even!
            epack.evec[i].Checkerboard() = Even;
        }
    } else {
        for (unsigned int i = 0; i < par().size; i++)
            epack.evec[i].Checkerboard() = Odd;
    }
    
//    // assign 4d evecs to 5th dim.
//    int Ls = par().Ls;
//    for (unsigned int i = 0; i < par().size; i++){
//        int s=i/Ls;
//        &(epack.evec5d[s+i*Ls]) = &(epack.evec[i]);
//    }
    
    if (!par().gaugeXform.empty())
    {

        LOG(Message) << "Applying gauge transformation to eigenvectors " << getName()
                     << " using " << par().gaugeXform << std::endl;
        auto &xform = envGet(GaugeMat, par().gaugeXform);
        envGetTmp(GaugeMat,    tmpXform);
        envGetTmp(GaugeMat, tmpXformOdd);

        if (par().Ls > 1) 
        {
            LOG(Message) << "Creating 5d GaugeMat from " << par().gaugeXform << std::endl;
            startTimer("5-d gauge transform creation");
            for (unsigned int j = 0; j < par().Ls; j++)
            {
                InsertSlice(xform, tmpXform, j, 0);
            }
            stopTimer("5-d gauge transform creation");
        }

        pickCheckerboard(Odd, tmpXformOdd, tmpXform);
        startTimer("Transform application");
        for (unsigned int i = 0; i < par().size; i++)
        {
            LOG(Message) << "Applying gauge transformation to eigenvector i = " << i << "/" << par().size << std::endl;
            epack.evec[i].Checkerboard() = Odd;
            epack.evec[i] = tmpXformOdd * epack.evec[i];
        }
        stopTimer("Transform application");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_StagLoadEigenPack_hpp_
