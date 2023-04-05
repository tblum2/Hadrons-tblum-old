/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/StagA2AMesonFieldCCHalfMem.hpp

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
#ifndef Hadrons_MContraction_StagA2AMesonFieldCCHalfMem_hpp_
#define Hadrons_MContraction_StagA2AMesonFieldCCHalfMem_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/EigenPack.hpp>
//#include <Hadrons/utils_memory.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     All-to-all meson field creation                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class StagA2AMesonFieldCCHalfMemPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagA2AMesonFieldCCHalfMemPar,
                                    int, cacheBlock,
                                    int, block,
                                    int, checkerboard,
                                    std::string, gauge,
                                    std::string, action,
                                    std::string, eigenPack,
                                    std::string, output,
                                    std::string, gammas,
                                    std::vector<std::string>, mom,
                                    double, mass);
};

class StagA2AMesonFieldCCHalfMemMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagA2AMesonFieldCCHalfMemMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra, gamma);
};

template <typename T, typename FImpl>
class StagMesonFieldCCHalfMemKernel: public A2AKernel<T, typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    StagMesonFieldCCHalfMemKernel(const std::vector<Gamma::Algebra> &gamma,
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

    virtual ~StagMesonFieldCCHalfMemKernel(void) = default;

    void operator()(A2AMatrixSet<T> &m,
                    int mu,
                    //LatticeGaugeField &U,
                    FermionOperator<FImpl> &Dns,
                    const LatticeColourMatrix &Umu,
                    const FermionField *evec,
                    const Real *eval,
                    //const Real mass,
                    const unsigned int orthogDim,
                    double &t)
    {
        A2Autils<FImpl>::StagMesonFieldCCHalfMem(m, mu, Dns, Umu, evec, eval, orthogDim, &t);
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
    void operator()(A2AMatrixSet<T> &m, const FermionField *left,
                    const FermionField *right,
                    const unsigned int orthogDim, double &t)
    {
       //unimplemented
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
class TStagA2AMesonFieldCCHalfMem : public Module<StagA2AMesonFieldCCHalfMemPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef A2AMatrixBlockComputation<Complex,
                                      FermionField,
                                      StagA2AMesonFieldCCHalfMemMetadata,
                                      HADRONS_A2AM_IO_TYPE> Computation;
    typedef StagMesonFieldCCHalfMemKernel<Complex, FImpl> Kernel;
public:
    // constructor
    TStagA2AMesonFieldCCHalfMem(const std::string name);
    // destructor
    virtual ~TStagA2AMesonFieldCCHalfMem(void){};
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

MODULE_REGISTER(StagA2AMesonFieldCCHalfMem, ARG(TStagA2AMesonFieldCCHalfMem<STAGIMPL>), MContraction);

/******************************************************************************
*                  TStagA2AMesonFieldCCHalfMem implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStagA2AMesonFieldCCHalfMem<FImpl>::TStagA2AMesonFieldCCHalfMem(const std::string name)
: Module<StagA2AMesonFieldCCHalfMemPar>(name)
, momphName_(name + "_momph")
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStagA2AMesonFieldCCHalfMem<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge, par().eigenPack, par().action};

    return in;
}

template <typename FImpl>
std::vector<std::string> TStagA2AMesonFieldCCHalfMem<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagA2AMesonFieldCCHalfMem<FImpl>::setup(void)
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
    envTmpLat(ComplexField, "coor");
    //printMem("StagMesonFieldCC setup(): after envTmpLat ", env().getGrid()->ThisRank());
    envTmp(Computation, "computation", 1, envGetGrid(FermionField),
           env().getNd() - 1, mom_.size(), gamma_.size(), par().block,
           par().cacheBlock, this);
}


// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagA2AMesonFieldCCHalfMem<FImpl>::execute(void)
{
    auto &epack  = envGet(BaseFermionEigenPack<FImpl>, par().eigenPack);
    int N        = epack.evec.size();
    // set checkerboard of evecs
    for(int j=0; j<N; j++)
      epack.evec[j].Checkerboard()=par().checkerboard;

    auto &U = envGet(LatticeGaugeField, par().gauge);
    auto &Dns = envGet(FMat, par().action);
    
    int nt         = env().getDim().back();
    int ngamma     = gamma_.size();
    assert(ngamma==1);// do one at a time
    int nmom       = mom_.size();
    int block      = par().block;
    int cacheBlock = par().cacheBlock;
    double mass = par().mass;

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
    LOG(Message) << "Meson field chunk size: " << nt << "*" << N << "*" << N
    << " (filesize " << sizeString(nt*N*N*sizeof(HADRONS_A2AM_IO_TYPE))
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
        StagA2AMesonFieldCCHalfMemMetadata md;

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
    //ph[0] = 1.0;
    //ph[0] = where( mod(sum,2)==(Integer)0, ph[0],-ph[0]);

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
    std::vector<Real> eval(epack.eval.size());
    
    for (unsigned int i = 0; i < N; i++)
    {
        // imag. part of eval of unpreconditioned Dirac op
        // needed to recontruct even from odd and vs.
        eval[i]=sqrt(epack.eval[i]-mass*mass);
    }

    Kernel      kernel(gamma_, ph, envGetGrid(FermionField));

    envGetTmp(Computation, computation);
    computation.execute(mu, Dns, Umu, epack.evec, eval, kernel,
                        ionameFn, filenameFn, metadataFn);
    
    // save +eval
    std::vector<complex<double>> ev(N);
    for (unsigned int i = 0; i < N; i++)
    {
        ev[i] = complex<double>(mass,eval[i]);
    }
    
    if ( env().getGrid()->IsBoss() ) {
        LOG(Message)<<"Saving evals "<<std::endl;
        std::string eval_filename = par().output + "." + std::to_string(vm().getTrajectory()) + "/evals.h5";
        A2AVectorsIo::initEvalFile(eval_filename,
                               ev.size());// total size
        A2AVectorsIo::saveEvalBlock(eval_filename,
                                    ev.data(),
                                    0,// start of chunk
                                    N);// size of chunk saved
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_StagA2AMesonFieldCCHalfMem_hpp_
