/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2ANucleonField.hpp

Copyright (C) 2015-2018

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
#ifndef Hadrons_MContraction_A2ANucleonField_hpp_
#define Hadrons_MContraction_A2ANucleonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
// MCA - change this to A2AMatrixNucleon.hpp?
#include <Hadrons/A2AMatrixNucleon.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     All-to-all nucleon field creation                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2ANucleonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ANucleonFieldPar,
                                    int, cacheBlock,
                                    int, block,
                                    std::string, left,
                                    std::string, right,
                                    std::string, q3,
                                    std::string, output,
                                    std::string, gammas, // may not need this for 2pt
                                    std::vector<std::string>, mom);
};

class A2ANucleonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ANucleonFieldMetadata,
                                    std::vector<RealF>, momentum
                                    /*Gamma::Algebra, gamma*/);
};

template <typename T, typename FImpl>
class NucleonFieldKernel: public A2AKernelNucleon<T, typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField FermionField;
public:
    NucleonFieldKernel(const std::vector<Gamma::Algebra> &gamma,
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

    virtual ~NucleonFieldKernel(void) = default;
    // MCA - adjust this for nucleon? [A2AMatrixNucleonSet]
    virtual void operator()(A2AMatrixSetNuc<T> &m, const FermionField *left, 
                            const FermionField *right,
                            const FermionField *q3,
                            const unsigned int orthogDim, double &t)
    {
		// MCA - adjust this for nucleon [A2Autils<FImpl>::NucleonField(m, left, right, q3, gamma_, mom_, orthogDim, &t);]
        A2Autils<FImpl>::NucleonField(m, left, right, q3, gamma_, mom_, orthogDim, &t);
    }

    virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej, const unsigned int blockSizek)
    {
        return vol_*(2*8.0+6.0+8.0*mom_.size())*blockSizei*blockSizej*blockSizek*gamma_.size();
    }

    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej, const unsigned int blockSizek)
    {
        return vol_*(12.0*sizeof(T))*blockSizei*blockSizej*blockSizek
               +  vol_*(2.0*sizeof(T)*mom_.size())*blockSizei*blockSizej*blockSizek*gamma_.size();
    }
private:
    const std::vector<Gamma::Algebra> &gamma_;
    const std::vector<LatticeComplex> &mom_;
    GridBase                          *grid_;
    double                            vol_;
};

template <typename FImpl>
class TA2ANucleonField : public Module<A2ANucleonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    // MCA - UPDATE THIS FOR NUCLEON
    typedef A2AMatrixNucleonBlockComputation<Complex, 
                                      FermionField, 
                                      A2ANucleonFieldMetadata, 
                                      HADRONS_A2AN_IO_TYPE> Computation;
    typedef NucleonFieldKernel<Complex, FImpl> Kernel;
public:
    // constructor
    TA2ANucleonField(const std::string name);
    // destructor
    virtual ~TA2ANucleonField(void){};
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

MODULE_REGISTER(A2ANucleonField, ARG(TA2ANucleonField<FIMPL>), MContraction);

/******************************************************************************
*                  TA2ANucleonField implementation                            *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2ANucleonField<FImpl>::TA2ANucleonField(const std::string name)
: Module<A2ANucleonFieldPar>(name)
, momphName_(name + "_momph")
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2ANucleonField<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().left, par().right, par().q3};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2ANucleonField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ANucleonField<FImpl>::setup(void)
{
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
    envTmp(Computation, "computation", 1, envGetGrid(FermionField), 
           env().getNd() - 1, mom_.size(), Ns, par().block, 
           par().cacheBlock, this);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ANucleonField<FImpl>::execute(void)
{
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);
    auto &q3    = envGet(std::vector<FermionField>, par().q3);

    int nt         = env().getDim().back();
    int N_i        = left.size();
    int N_j        = right.size();
    int N_k        = q3.size();
    int ngamma     = Ns;
    int nmom       = mom_.size();
    int block      = par().block;
    int cacheBlock = par().cacheBlock;

    LOG(Message) << "Computing all-to-all nucleon fields" << std::endl;
    LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right << "' Q3: '" << par().q3 << "'" << std::endl;
    LOG(Message) << "Momenta:" << std::endl;
    for (auto &p: mom_)
    {
        LOG(Message) << "  " << p << std::endl;
    }
    // MCA - Disabling this for nucleon
    /*
    LOG(Message) << "Spin bilinears:" << std::endl;
    for (auto &g: gamma_)
    {
        LOG(Message) << "  " << g << std::endl;
    }
    */
    LOG(Message) << "Nucleon field size: " << Ns << "*" << nt << "*" << N_i << "*" << N_j << "*" << N_k
                 << " (filesize " << sizeString(Ns*nt*N_i*N_j*N_k*sizeof(HADRONS_A2AN_IO_TYPE)) 
                 << "/momentum/bilinear)" << std::endl;

    auto &ph = envGet(std::vector<ComplexField>, momphName_);

    if (!hasPhase_)
    {
        startTimer("Momentum phases");
        for (unsigned int j = 0; j < nmom; ++j)
        {
            Complex           i(0.0,1.0);
            std::vector<Real> p;

            envGetTmp(ComplexField, coor);
            ph[j] = Zero;
            for(unsigned int mu = 0; mu < mom_[j].size(); mu++)
            {
                LatticeCoordinate(coor, mu);
                ph[j] = ph[j] + (mom_[j][mu]/env().getDim(mu))*coor;
            }
            ph[j] = exp((Real)(2*M_PI)*i*ph[j]);
        }
        hasPhase_ = true;
        stopTimer("Momentum phases");
    }

    auto ionameFn = [this](const unsigned int m)
    {
        std::stringstream ss;

        ss << "Nucleon_";
        for (unsigned int mu = 0; mu < mom_[m].size(); ++mu)
        {
            ss << mom_[m][mu] << ((mu == mom_[m].size() - 1) ? "" : "_");
        }

        return ss.str();
    };

    auto filenameFn = [this, &ionameFn](const unsigned int m)
    {
        return par().output + "." + std::to_string(vm().getTrajectory()) 
               + "/" + ionameFn(m) + ".h5";
    };

    auto metadataFn = [this](const unsigned int m)
    {
        A2ANucleonFieldMetadata md;

        for (auto pmu: mom_[m])
        {
            md.momentum.push_back(pmu);
        }
        //md.gamma = gamma_[g];
        
        return md;
    };

    Kernel      kernel(gamma_, ph, envGetGrid(FermionField));

    envGetTmp(Computation, computation);
    // Update A2AMatrixBlockComputation for this
    computation.execute(left, right, q3, kernel, ionameFn, filenameFn, metadataFn);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2ANucleonField_hpp_
