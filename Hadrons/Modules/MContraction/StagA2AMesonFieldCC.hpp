/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/StagA2AMesonFieldCC.hpp

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
#ifndef Hadrons_MContraction_StagA2AMesonFieldCC_hpp_
#define Hadrons_MContraction_StagA2AMesonFieldCC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/EigenPack.hpp>
//#include <Hadrons/utils_memory.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     All-to-all meson field creation                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class StagA2AMesonFieldCCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagA2AMesonFieldCCPar,
                                    int, cacheBlock,
                                    int, block,
                                    //int, Size,
                                    //int, istart,
                                    //int, jstart,
                                    std::string, gauge,
                                    //std::string, left,
                                    //std::string, right,
                                    std::string, eigenPack,
                                    std::string, output,
                                    std::string, gammas,
                                    std::vector<std::string>, mom);
};

class StagA2AMesonFieldCCMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagA2AMesonFieldCCMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra, gamma);
};

template <typename T, typename FImpl>
class StagMesonFieldCCKernel: public A2AKernel<T, typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    StagMesonFieldCCKernel(const std::vector<Gamma::Algebra> &gamma,
                           const std::vector<LatticeComplex> &mom,
                           GridBase *grid)
    : gamma_(gamma), mom_(mom), grid_(grid)
    {
        vol_ = 1.;
        for (auto &d: grid_->GlobalDimensions())
        {
            vol_ *= d;
        }
    }

    virtual ~StagMesonFieldCCKernel(void) = default;

    void operator()(A2AMatrixSet<T> &m,
                            int mu,
                            const LatticeColourMatrix &Umu,
                            const FermionField *left,
                            const FermionField *right,
                            const unsigned int orthogDim,
                            double &t)
    {
        A2Autils<FImpl>::StagMesonFieldCC(m, mu, Umu, left, right, gamma_, mom_, orthogDim, &t);
    }
    void operator()(A2AMatrixSet<T> &m, const FermionField *left,
                    const FermionField *right,
                    const unsigned int orthogDim, double &t)
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
        // updated for staggered
        return vol_*(22.0+6.0*mom_.size())*blockSizei*blockSizej;
    }

    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        // 3.0 ? for colors. updated for staggered
        return vol_*(3.0*sizeof(T))*blockSizei*blockSizej
               +  vol_*(2.0*sizeof(T)*mom_.size())*blockSizei*blockSizej;
    }
private:
    //const LatticeGaugeField &U_;
    const std::vector<Gamma::Algebra> &gamma_;
    const std::vector<LatticeComplex> &mom_;
    GridBase                          *grid_;
    double                            vol_;
};

template <typename FImpl>
class TStagA2AMesonFieldCC : public Module<StagA2AMesonFieldCCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef A2AMatrixBlockComputation<Complex,
                                      FermionField,
                                      StagA2AMesonFieldCCMetadata,
                                      HADRONS_A2AM_IO_TYPE> Computation;
    typedef StagMesonFieldCCKernel<Complex, FImpl> Kernel;
public:
    // constructor
    TStagA2AMesonFieldCC(const std::string name);
    // destructor
    virtual ~TStagA2AMesonFieldCC(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool                               hasPhase_{false};
    std::string                        momphName_;
    std::vector<Gamma::Algebra>        gamma_;
    std::vector<std::vector<Real>>     mom_;
};

MODULE_REGISTER(StagA2AMesonFieldCC, ARG(TStagA2AMesonFieldCC<STAGIMPL>), MContraction);

/******************************************************************************
*                  TStagA2AMesonFieldCC implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStagA2AMesonFieldCC<FImpl>::TStagA2AMesonFieldCC(const std::string name)
: Module<StagA2AMesonFieldCCPar>(name)
, momphName_(name + "_momph")
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStagA2AMesonFieldCC<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge, par().eigenPack};

    return in;
}

template <typename FImpl>
std::vector<std::string> TStagA2AMesonFieldCC<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagA2AMesonFieldCC<FImpl>::setup(void)
{
    //printMem("Begin StagMesonFieldCC setup() ", env().getGrid()->ThisRank());
    gamma_.clear();
    mom_.clear();
    if (par().gammas == "all")
    {
        gamma_ = {
            Gamma::Algebra::Gamma5,
            Gamma::Algebra::Identity,
            Gamma::Algebra::GammaX,
            Gamma::Algebra::GammaY,
            Gamma::Algebra::GammaZ,
            Gamma::Algebra::GammaT,
            Gamma::Algebra::GammaXGamma5,
            Gamma::Algebra::GammaYGamma5,
            Gamma::Algebra::GammaZGamma5,
            Gamma::Algebra::GammaTGamma5,
            Gamma::Algebra::SigmaXY,
            Gamma::Algebra::SigmaXZ,
            Gamma::Algebra::SigmaXT,
            Gamma::Algebra::SigmaYZ,
            Gamma::Algebra::SigmaYT,
            Gamma::Algebra::SigmaZT
        };
    }
    else
    {
        gamma_ = strToVec<Gamma::Algebra>(par().gammas);
    }
    for (auto &pstr: par().mom)
    {
        auto p = strToVec<Real>(pstr);

        if (p.size() != env().getNd() - 1)
        {
            HADRONS_ERROR(Size, "Momentum has " + std::to_string(p.size())
                                + " components instead of "
                                + std::to_string(env().getNd() - 1));
        }
        mom_.push_back(p);
    }
    envCache(std::vector<ComplexField>, momphName_, 1,
             par().mom.size(), envGetGrid(ComplexField));
    //printMem("StagMesonFieldCC setup(): after envCache ", env().getGrid()->ThisRank());
    envTmpLat(ComplexField, "coor");
    //printMem("StagMesonFieldCC setup(): after envTmpLat ", env().getGrid()->ThisRank());
    envTmp(Computation, "computation", 1, envGetGrid(FermionField),
           env().getNd() - 1, mom_.size(), gamma_.size(), par().block,
           par().cacheBlock, this);
    //printMem("StagMesonFieldCC setup() End ", env().getGrid()->ThisRank());
}


// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagA2AMesonFieldCC<FImpl>::execute(void)
{
    auto &epack  = envGet(BaseFermionEigenPack<FImpl>, par().eigenPack);
    //auto &left  = envGet(std::vector<FermionField>, par().left);
    int N_i        = epack.evec.size();
    //auto &right = envGet(std::vector<FermionField>, par().left);// same as left
    int N_j        = epack.evec.size();

    auto &U = envGet(LatticeGaugeField, par().gauge);
    int nt         = env().getDim().back();
    int ngamma     = gamma_.size();
    assert(ngamma==1);// do one at a time
    int nmom       = mom_.size();
    int block      = par().block;
    int cacheBlock = par().cacheBlock;

    LOG(Message) << "Computing all-to-all meson fields" << std::endl;
    LOG(Message) << "Left and right: '" << par().eigenPack << "'" << std::endl;
    LOG(Message) << "Momenta:" << std::endl;
    for (auto &p: mom_)
    {
        LOG(Message) << "  " << p << std::endl;
    }
    LOG(Message) << "Spin bilinears:" << std::endl;
    for (auto &g: gamma_)
    {
        LOG(Message) << "  " << g << std::endl;
    }
    LOG(Message) << "Meson field chunk size: " << nt << "*" << 2*N_i << "*" << 2*N_j
    << " (filesize " << sizeString(nt*2*N_i*2*N_j*sizeof(HADRONS_A2AM_IO_TYPE))
    << "/momentum/bilinear)" << std::endl;

    auto &ph = envGet(std::vector<ComplexField>, momphName_);

    auto ionameFn = [this](const unsigned int m, const unsigned int g)
    {
        std::stringstream ss;

        ss << gamma_[g] << "_";
        for (unsigned int mu = 0; mu < mom_[m].size(); ++mu)
        {
            ss << mom_[m][mu] << ((mu == mom_[m].size() - 1) ? "" : "_");
        }

        return ss.str();
    };

    auto filenameFn = [this, &ionameFn](const unsigned int m, const unsigned int g)
    {
        return par().output + "." + std::to_string(vm().getTrajectory())
        + "/" + ionameFn(m, g) + ".h5";
    };

    auto metadataFn = [this](const unsigned int m, const unsigned int g)
    {
        StagA2AMesonFieldCCMetadata md;

        for (auto pmu: mom_[m])
        {
            md.momentum.push_back(pmu);
        }
        md.gamma = gamma_[g];

        return md;
    };

    // Staggered Phases. Do spatial and temporal gamma only
    Lattice<iScalar<vInteger> > x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(U.Grid()); LatticeCoordinate(z,2);
    Lattice<iScalar<vInteger> > t(U.Grid()); LatticeCoordinate(t,3);
    Lattice<iScalar<vInteger> > lin_z(U.Grid()); lin_z=x+y;
    Lattice<iScalar<vInteger> > lin_t(U.Grid()); lin_t=x+y+z;
    Lattice<iScalar<vInteger> > sum(U.Grid()); sum=lin_z+z+t;
    ph[0] = 1.0;
    ph[0] = where( mod(sum,2)==(Integer)0, ph[0],-ph[0]);

    ComplexField phases(U.Grid());
    phases=1.0;
    int mu;
    if(gamma_[0]==Gamma::Algebra::GammaX)mu=0;
    else if(gamma_[0]==Gamma::Algebra::GammaY){
        mu=1;
        phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
    } else if(gamma_[0]==Gamma::Algebra::GammaZ){
        mu=2;
        phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
    } else if(gamma_[0]==Gamma::Algebra::GammaT){
        mu=3;
        phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);
    } else assert(0);

    LatticeColourMatrix Umu(U.Grid());
    Umu = PeekIndex<LorentzIndex>(U,mu);
    Umu *= phases;

    Kernel      kernel(gamma_, ph, envGetGrid(FermionField));

    envGetTmp(Computation, computation);
    computation.execute(mu, Umu, epack.evec, epack.evec, kernel,
                        ionameFn, filenameFn, metadataFn);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_StagA2AMesonFieldCC_hpp_
