/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/StagSparseA2AMesonField.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#ifndef Hadrons_MContraction_StagSparseA2AMesonField_hpp_
#define Hadrons_MContraction_StagSparseA2AMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     All-to-all meson field creation                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class StagSparseA2AMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagSparseA2AMesonFieldPar,
                                    int, cacheBlock,
                                    int, block,
                                    std::string, left,
                                    std::string, right,
                                    std::string, output);
};

class StagSparseA2AMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagSparseA2AMesonFieldMetadata,
                                    std::string,momstr,
                                    std::string,gamstr);
};

template <typename T, typename FImpl>
class StagSparseMesonFieldKernel: public A2AKernel<T, typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    StagSparseMesonFieldKernel(GridBase *grid)
    : grid_(grid)
    {
        vol_ = 1.;
        for (auto &d: grid_->GlobalDimensions())
        {
            vol_ *= d;
        }
    }

    virtual ~StagSparseMesonFieldKernel(void) = default;
    
    void operator()(A2AMatrixSet<T> &m,
                            const FermionField *left,
                            const FermionField *right,
                            const unsigned int orthogDim,
                            double &t)
    {
        A2Autils<FImpl>::StagMesonField(m, left, right, orthogDim, &t);
    }
    
    void operator()(A2AMatrixSet<T> &m,
                    int mu,
                    const LatticeColourMatrix &Umu,
                    const FermionField *left,
                    const FermionField *right,
                    const unsigned int orthogDim,
                    double &t)
    {
        //unimplemented
        assert(0);
    }
    void operator()(A2AMatrixSet<T> &m,
                    int mu,
                    FermionOperator<FImpl> &Dns,
                    //LatticeGaugeField &U,
                    const LatticeColourMatrix &Umu,
                    const FermionField *levec,
                    const FermionField *revec,
                    const Real *eval,
                    //const Real mass,
                    const unsigned int orthogDim, double &time)
    {
        //unimplemented (for Staggered conserved current)
        assert(0);
    }
    
    virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        // needs to be updated for staggered
        return vol_*(2*8.0+6.0+8.0)*blockSizei*blockSizej;
        //return vol_*(22.0+6.0*mom_.size())*blockSizei*blockSizej;
    }

    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        // 3.0 ? for colors
        return vol_*(3.0*sizeof(T))*blockSizei*blockSizej
               +  vol_*(2.0*sizeof(T))*blockSizei*blockSizej;
    }
private:
    GridBase                          *grid_;
    double                            vol_;
};

template <typename FImpl>
class TStagSparseA2AMesonField : public Module<StagSparseA2AMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef A2AMatrixBlockComputation<Complex, 
                                      FermionField, 
                                      StagSparseA2AMesonFieldMetadata,
                                      HADRONS_A2AM_IO_TYPE> Computation;
    typedef StagSparseMesonFieldKernel<Complex, FImpl> Kernel;
public:
    // constructor
    TStagSparseA2AMesonField(const std::string name);
    // destructor
    virtual ~TStagSparseA2AMesonField(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
};

MODULE_REGISTER(StagSparseA2AMesonField, ARG(TStagSparseA2AMesonField<STAGIMPL>), MContraction);

/******************************************************************************
*                  TStagSparseA2AMesonField implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStagSparseA2AMesonField<FImpl>::TStagSparseA2AMesonField(const std::string name)
: Module<StagSparseA2AMesonFieldPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStagSparseA2AMesonField<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().left, par().right};

    return in;
}

template <typename FImpl>
std::vector<std::string> TStagSparseA2AMesonField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagSparseA2AMesonField<FImpl>::setup(void)
{
    envTmp(Computation, "computation", 1, envGetGrid(FermionField),
           env().getNd() - 1, 1, 1, par().block,
           par().cacheBlock, this);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagSparseA2AMesonField<FImpl>::execute(void)
{
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);

    int nt         = env().getDim().back();
    int N_i        = left.size();
    int N_j        = right.size();
    int block      = par().block;
    int cacheBlock = par().cacheBlock;

    LOG(Message) << "Computing all-to-all Sparse meson fields" << std::endl;
    LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right << "'" << std::endl;
    
    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j 
                 << " (filesize " << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE)) 
                 << "/momentum/bilinear)" << std::endl;
 
    auto ionameFn = [this](const unsigned int m, const unsigned int g)
    {
        std::stringstream ss;

        ss << "_";

        return ss.str();
    };
    
    auto filenameFn = [this, &ionameFn](const unsigned int m, const unsigned int g)
    {
        return par().output + "." + std::to_string(vm().getTrajectory())
               + "/" + ionameFn(1,1) + ".h5";
    };

    auto metadataFn = [this](const unsigned int m, const unsigned int g)
    {
        StagSparseA2AMesonFieldMetadata md;
        
        md.momstr='0 0 0';
        md.gamstr = 'n/a';
        
        return md;
    };

    Kernel      kernel(envGetGrid(FermionField));

    envGetTmp(Computation, computation);
    computation.execute(left, right, kernel, ionameFn, filenameFn, metadataFn);
}



END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_StagSparseA2AMesonField_hpp_
