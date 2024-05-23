/*
 * A2AVectors.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
 * Author: fionnoh <fionnoh@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MSolver_A2AVectors_hpp_
#define Hadrons_MSolver_A2AVectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Create all-to-all V & W vectors                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2AVectorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsPar,
                                    std::string, noise,
                                    std::string, action,
                                    std::string, gauge,
                                    std::string, eigenPack,
                                    std::string, solver,
                                    std::string, output,
                                    int, inc,
                                    int, tinc,
                                    bool, doubleMemory,
                                    double, mass,
                                    bool,        multiFile);
};

template <typename FImpl, typename Pack>
class TA2AVectors : public Module<A2AVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef HADRONS_DEFAULT_SCHUR_A2A<FImpl> A2A;
public:
    // constructor
    TA2AVectors(const std::string name);
    // destructor
    virtual ~TA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string  solverName_;
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(A2AVectors, 
    ARG(TA2AVectors<FIMPL, BaseFermionEigenPack<FIMPL>>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectors, 
    ARG(TA2AVectors<ZFIMPL, BaseFermionEigenPack<ZFIMPL>>), MSolver);

/******************************************************************************
 *                       TA2AVectors implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TA2AVectors<FImpl, Pack>::TA2AVectors(const std::string name)
: Module<A2AVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectors<FImpl, Pack>::getInput(void)
{
    std::string              sub_string;
    std::vector<std::string> in;

    if (!par().eigenPack.empty())
    {
        in.push_back(par().eigenPack);
        sub_string = (!par().eigenPack.empty()) ? "_subtract" : "";
    }
    in.push_back(par().solver + sub_string);
    in.push_back(par().noise);

    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectors<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "_v", getName() + "_w"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2AVectors<FImpl, Pack>::setup(void)
{
    bool        hasLowModes = (!par().eigenPack.empty());
    std::string sub_string  = (hasLowModes) ? "_subtract" : "";
    auto        &noise      = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto        &action     = envGet(FMat, par().action);
    auto        &solver     = envGet(Solver, par().solver + sub_string);
    int         Ls          = env().getObjectLs(par().action);

    if (hasLowModes)
    {
        auto &epack = envGet(Pack, par().eigenPack);
        Nl_ = epack.evec.size();
    }
    envCreate(std::vector<FermionField>, getName() + "_v", 1, 
              Nl_ + noise.fermSize(), envGetGrid(FermionField));
    envCreate(std::vector<FermionField>, getName() + "_w", 1, 
              Nl_ + noise.fermSize(), envGetGrid(FermionField));
    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action, solver);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2AVectors<FImpl, Pack>::execute(void)
{
    std::string sub_string = (Nl_ > 0) ? "_subtract" : "";
    auto        &action    = envGet(FMat, par().action);
    auto        &solver    = envGet(Solver, par().solver + sub_string);
    auto        &noise     = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto        &v         = envGet(std::vector<FermionField>, getName() + "_v");
    auto        &w         = envGet(std::vector<FermionField>, getName() + "_w");
    int         Ls         = env().getObjectLs(par().action);

    envGetTmp(A2A, a2a);

    if (Nl_ > 0)
    {
        LOG(Message) << "Computing all-to-all vectors "
                     << " using eigenpack '" << par().eigenPack << "' ("
                     << Nl_ << " low modes) and noise '"
                     << par().noise << "' (" << noise.fermSize() 
                     << " noise vectors)" << std::endl;
    }
    else
    {
        LOG(Message) << "Computing all-to-all vectors "
                     << " using noise '" << par().noise << "' (" << noise.fermSize() 
                     << " noise vectors)" << std::endl;
    }
    // Low modes
    for (unsigned int il = 0; il < Nl_; il++)
    {
        auto &epack  = envGet(Pack, par().eigenPack);

        startTimer("V low mode");
        LOG(Message) << "V vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeV(v[il], epack.evec[il], epack.eval[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeV5D(v[il], f5, epack.evec[il], epack.eval[il]);
        }
        stopTimer("V low mode");
        startTimer("W low mode");
        LOG(Message) << "W vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeW(w[il], epack.evec[il], epack.eval[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeW5D(w[il], f5, epack.evec[il], epack.eval[il]);
        }
        stopTimer("W low mode");
    }

    // High modes
    for (unsigned int ih = 0; ih < noise.fermSize(); ih++)
    {
        startTimer("V high mode");
        LOG(Message) << "V vector i = " << Nl_ + ih
                     << " (" << ((Nl_ > 0) ? "high " : "") 
                     << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeV(v[Nl_ + ih], noise.getFerm(ih));
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeV5D(v[Nl_ + ih], f5, noise.getFerm(ih));
        }
        stopTimer("V high mode");
        startTimer("W high mode");
        LOG(Message) << "W vector i = " << Nl_ + ih
                     << " (" << ((Nl_ > 0) ? "high " : "") 
                     << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeW(w[Nl_ + ih], noise.getFerm(ih));
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeW5D(w[Nl_ + ih], f5, noise.getFerm(ih));
        }
        stopTimer("W high mode");
    }

    // Print out v mode norms
    for (unsigned int il = 0; il < Nl_; il++)
        {
          LOG(Message) << "V vector i = " << il << " (low mode)" << " | norm " << norm2(v[il]) << std::endl;
        }
    
    // Print out w mode norms
    for (unsigned int il = 0; il < Nl_; il++)
        {
          LOG(Message) << "W vector i = " << il << " (low mode)" << " | norm " << norm2(w[il]) << std::endl;
    }

    // I/O if necessary
    if (!par().output.empty())
    {
        startTimer("V I/O");
        A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
        stopTimer("V I/O");
        startTimer("W I/O");
        A2AVectorsIo::write(par().output + "_w", w, par().multiFile, vm().getTrajectory());
        stopTimer("W I/O");
    }
}

template <typename FImpl, typename Pack>
class TStagA2AVectors : public Module<A2AVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef A2AVectorsSchurStaggered<FImpl> A2A;
public:
    // constructor
    TStagA2AVectors(const std::string name);
    // destructor
    virtual ~TStagA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string  solverName_;
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(StagA2AVectors,
                    ARG(TStagA2AVectors<STAGIMPL, BaseFermionEigenPack<STAGIMPL>>),
                    MSolver);

/******************************************************************************
 *                       TStagA2AVectors implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TStagA2AVectors<FImpl, Pack>::TStagA2AVectors(const std::string name)
: Module<A2AVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TStagA2AVectors<FImpl, Pack>::getInput(void)
{
    std::string              sub_string;
    std::vector<std::string> in;
    
    if (!par().eigenPack.empty())
    {
        in.push_back(par().eigenPack);
        sub_string = (!par().eigenPack.empty()) ? "_subtract" : "";
    }
    in.push_back(par().solver);
    in.push_back(par().noise);
    
    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TStagA2AVectors<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "_v", getName() + "_w"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagA2AVectors<FImpl, Pack>::setup(void)
{
    bool        hasLowModes = (!par().eigenPack.empty());
    std::string sub_string  = (hasLowModes) ? "_subtract" : "";
    auto        &noise      = envGet(ColorDiagonalNoise<FImpl>, par().noise);
    auto        &action     = envGet(FMat, par().action);
    auto        &solver     = envGet(Solver, par().solver);
    int         Ls          = env().getObjectLs(par().action);
    
    
    if (hasLowModes)
    {
        auto &epack = envGet(Pack, par().eigenPack);
        Nl_ = epack.evec.size();
    }
    envCreate(std::vector<FermionField>, getName() + "_v", 1,
              2*Nl_ + noise.fermSize(), envGetGrid(FermionField));
    envCreate(std::vector<FermionField>, getName() + "_w", 1,
              2*Nl_ + noise.fermSize(), envGetGrid(FermionField));
    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action, solver);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagA2AVectors<FImpl, Pack>::execute(void)
{
    auto        &action    = envGet(FMat, par().action);
    auto        &solver    = envGet(Solver, par().solver);
    auto        &noise     = envGet(ColorDiagonalNoise<FImpl>, par().noise);
    auto        &v         = envGet(std::vector<FermionField>, getName() + "_v");
    auto        &w         = envGet(std::vector<FermionField>, getName() + "_w");
    int         Ls         = env().getObjectLs(par().action);
    double      mass       = par().mass;
    
    envGetTmp(A2A, a2a);
    
    if (Nl_ > 0)
    {
        LOG(Message) << "Computing all-to-all vectors "
        << " using eigenpack '" << par().eigenPack << "' ("
        << 2*Nl_ << " low modes) and noise '"
        << par().noise << "' (" << noise.fermSize()
        << " noise vectors)" << std::endl;
    }
    else
    {
        LOG(Message) << "Computing all-to-all vectors "
        << " using noise '" << par().noise << "' (" << noise.size()
        << " noise vectors)" << std::endl;
    }
    // Low modes
    auto &epack  = envGet(Pack, par().eigenPack);
    for (unsigned int il = 0; il < Nl_; il++)
    {
        // eval of unpreconditioned Dirac op
        std::complex<double> eval(mass,sqrt(epack.eval[il]-mass*mass));
        
        startTimer("V low mode");
        LOG(Message) << "V vector i = " << 2*il << ", " << 2*il+1 << " (low modes)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeV(v[2*il], epack.evec[il], eval);
            // construct -lambda evec
            a2a.makeLowModeV(v[2*il+1], epack.evec[il], eval, 1);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeV5D(v[2*il], f5, epack.evec[il], eval);
            // construct -lambda evec
            a2a.makeLowModeV5D(v[2*il+1], f5, epack.evec[il], eval, 1);
        }
        stopTimer("V low mode");
        startTimer("W low mode");
        LOG(Message) << "W vector i = " << 2*il << ", " << 2*il+1 << " (low modes)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeW(w[2*il], epack.evec[il], eval);
            // construct -lambda evec
            a2a.makeLowModeW(w[2*il+1], epack.evec[il], eval, 1);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeW5D(w[2*il], f5, epack.evec[il], eval);
            // construct -lambda evec
            a2a.makeLowModeW5D(w[2*il+1], f5, epack.evec[il], eval, 1);
        }
        stopTimer("W low mode");
    }
    
    // High modes
#if 1
    
    FermionField sub(env().getGrid());
    FermionField noise_ih(env().getGrid());
    
    for (unsigned int ih = 0; ih < noise.fermSize(); ih++)
    {
        startTimer("V high mode");
        LOG(Message) << "V vector i = " << 2*Nl_ + ih
        << " (" << ((Nl_ > 0) ? "high " : "")
        << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeV(v[2*Nl_ + ih], noise.getFerm(ih));
            // subtract the low modes
            sub = Zero();
            noise_ih = noise.getFerm(ih);
            for (int i=0;i<2*Nl_;i++) {
              const FermionField& tmp = v[i];
              // eval of unpreconditioned Dirac op
              std::complex<double> eval(mass,sqrt(epack.eval[i/2]-mass*mass));
              eval = i%2 ? eval : conjugate(eval);
              // need to subtract |l><l|/lambda_l |src>, so
              // mult, do not / by eval since v already has 1/eval
              axpy(sub,TensorRemove(innerProduct(tmp,noise_ih)) * eval,tmp,sub);
            }
            v[2*Nl_ + ih] -= sub;
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeV5D(v[2*Nl_ + ih], f5, noise.getFerm(ih));
        }
        std::cout << "v high " << ih << std::endl;
        //std::cout << v[2*Nl_+ih] << std::endl;
        
        stopTimer("V high mode");
        startTimer("W high mode");
        LOG(Message) << "W vector i = " << 2*Nl_ + ih
        << " (" << ((Nl_ > 0) ? "high " : "")
        << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeW(w[2*Nl_ + ih], noise.getFerm(ih));
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeW5D(w[2*Nl_ + ih], f5, noise.getFerm(ih));
        }
        stopTimer("W high mode");
    }
#endif
    // I/O if necessary
    if (!par().output.empty())
    {
        startTimer("V I/O");
        A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
        stopTimer("V I/O");
        startTimer("W I/O");
        A2AVectorsIo::write(par().output + "_w", w, par().multiFile, vm().getTrajectory());
        stopTimer("W I/O");
    }
    //printMem("End StagA2AVectors execute() ", env().getGrid()->ThisRank());
}

// Don't divide v vecs by eval. for use with cons. current meson fields.

template <typename FImpl, typename Pack>
class TStagNoEvalA2AVectors : public Module<A2AVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef A2AVectorsSchurStaggeredNoEval<FImpl> A2A;
public:
    // constructor
    TStagNoEvalA2AVectors(const std::string name);
    // destructor
    virtual ~TStagNoEvalA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string  solverName_;
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(StagNoEvalA2AVectors,
                    ARG(TStagNoEvalA2AVectors<STAGIMPL, BaseFermionEigenPack<STAGIMPL>>),
                    MSolver);

/******************************************************************************
 *                       TStagNoEvalA2AVectors implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TStagNoEvalA2AVectors<FImpl, Pack>::TStagNoEvalA2AVectors(const std::string name)
: Module<A2AVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TStagNoEvalA2AVectors<FImpl, Pack>::getInput(void)
{
    std::string              sub_string;
    std::vector<std::string> in;
    
    in.push_back(par().eigenPack);
    sub_string = (!par().eigenPack.empty()) ? "_subtract" : "";
    in.push_back(par().action);
    //in.push_back(par().vname);
    
    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TStagNoEvalA2AVectors<FImpl, Pack>::getOutput(void)
{
    
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagNoEvalA2AVectors<FImpl, Pack>::setup(void)
{
    bool        hasLowModes = (!par().eigenPack.empty());
    std::string sub_string  = (hasLowModes) ? "_subtract" : "";
    auto        &action     = envGet(FMat, par().action);
    int         Ls          = env().getObjectLs(par().action);
    
    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagNoEvalA2AVectors<FImpl, Pack>::execute(void)
{
    //printMem("Begin StagNoEvalA2AVectors execute() ", env().getGrid()->ThisRank());
    std::string sub_string = (Nl_ > 0) ? "_subtract" : "";
    auto        &action    = envGet(FMat, par().action);
    auto &epack  = envGet(Pack, par().eigenPack);
    Nl_ = epack.evec.size();
    std::vector<FermionField> v;
    if(!par().doubleMemory) v.resize(Nl_, envGetGrid(FermionField));
    int         Ls         = env().getObjectLs(par().action);
    double      mass       = par().mass;
    
    envGetTmp(A2A, a2a);
    
    LOG(Message) << "Computing all-to-all vectors "
    << " using eigenpack '" << par().eigenPack << "' ("
    << Nl_ << " low modes) '" << std::endl;

    //save for later
    std::vector<complex<double>> evalM(2*Nl_);
    
    // Low modes only, v(=w) vecs only
    // evecs start at index = start
    // index v vec from 0
    // if doubleMemory construct Even and Odd site (full) evecs
    // and write into evec. This is for meson field calc which computes
    // minus lambda contribution from plus
    
    FermionField temp(env().getGrid());
    FermionField tempRB(env().getRbGrid());
    
    for (unsigned int il = 0; il < Nl_; il++)
    {
        
        // eval of unpreconditioned Dirac op
        std::complex<double> eval(mass,sqrt(epack.eval[il]-mass*mass));
        
        startTimer("V low mode");
        LOG(Message) << "V vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            if(par().doubleMemory){
                // construct full eo evec from Odd
                pickCheckerboard(Odd,tempRB,epack.evec[il]);
                a2a.makeLowModeV(temp, tempRB, eval);
                epack.evec[il]=temp;
            }else{
                a2a.makeLowModeV(v[il], epack.evec[il], eval);
            }
            evalM[2*il] = eval;
            evalM[2*il+1] = conjugate(eval);
        }
        else
        {
            assert(0);
        }
        stopTimer("V low mode");
    }
    
    // I/O if necessary
    if (!par().output.empty())
    {
        std::string dir = dirname(par().output);
        int status = mkdir(dir);
        if (status)
        {
            HADRONS_ERROR(Io, "cannot create directory '" + dir
                          + "' ( " + std::strerror(errno) + ")");
        }
        if ( env().getGrid()->IsBoss() ) {
            
            std::string eval_filename = A2AVectorsIo::evalFilename(par().output,vm().getTrajectory());
            A2AVectorsIo::initEvalFile(eval_filename,
                                       evalM.size());// total size
            A2AVectorsIo::saveEvalBlock(eval_filename,
                                        evalM.data(),
                                        0,// start of chunk
                                        2*Nl_);// size of chunk saved
        }
    }
}
//
// Sparsened A2A vectors for staggered conserved current
//
template <typename FImpl, typename Pack>
class TStagSparseA2AVectors : public Module<A2AVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef A2AVectorsSchurStaggered<FImpl> A2A;
public:
    // constructor
    TStagSparseA2AVectors(const std::string name);
    // destructor
    virtual ~TStagSparseA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string  solverName_;
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(StagSparseA2AVectors,
                    ARG(TStagSparseA2AVectors<STAGIMPL, BaseFermionEigenPack<STAGIMPL>>),MSolver);

/******************************************************************************
 *                       TStagSparseA2AVectors implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TStagSparseA2AVectors<FImpl, Pack>::TStagSparseA2AVectors(const std::string name)
: Module<A2AVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TStagSparseA2AVectors<FImpl, Pack>::getInput(void)
{
    std::vector<std::string> in;
    
    in.push_back(par().eigenPack);
    in.push_back(par().gauge);
    in.push_back(par().action);
    in.push_back(par().solver);
    
    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TStagSparseA2AVectors<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "sparse_w"};
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagSparseA2AVectors<FImpl, Pack>::setup(void)
{
    auto        &action     = envGet(FMat, par().action);
    auto        &solver     = envGet(Solver, par().solver);
    
    auto &epack = envGet(Pack, par().eigenPack);
    Nl_ = epack.evec.size();
    // Sparse Grid
    
    // V, W vecs. 2x Nl_ for +-lambda
    // sparse vectors need sparse grid
    Coordinate sparseLatSize = envGetGrid(FermionField)->FullDimensions();
    sparseLatSize[0] /= par().inc;
    sparseLatSize[1] /= par().inc;
    sparseLatSize[2] /= par().inc;
    sparseLatSize[3] /= par().tinc;
    Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    GridCartesian sparseGrid(sparseLatSize,simd_layout,mpi_layout);
    envCreate(std::vector<FermionField>, "v", 1, 2*Nl_, &sparseGrid);
    envCreate(std::vector<FermionField>, "w", 1, 2*Nl_, &sparseGrid);
    envTmp(A2A, "a2a", 1, action, solver);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagSparseA2AVectors<FImpl, Pack>::execute(void)
{
    auto    &action = envGet(FMat, par().action);
    auto    &U      = envGet(LatticeGaugeField, par().gauge);
    auto    &solver = envGet(Solver, par().solver);
    auto    &epack  = envGet(Pack, par().eigenPack);
    double  mass    = par().mass;
    int     nt      = env().getDim(Tp);
    int     ns      = env().getDim(Xp);
    envGetTmp(A2A, a2a);
    
    auto &v=envGet(std::vector<FermionField>,"v");
    auto &w=envGet(std::vector<FermionField>,"w");
    v[0].Grid()->show_decomposition();
    w[0].Grid()->show_decomposition();
    epack.evec[0].Grid()->show_decomposition();

    LOG(Message) << "Computing all-to-all vectors using eigenpack " << par().eigenPack << " with " << 2*Nl_ << " low modes " << std::endl;

    // Staggered Phases. Do spatial gamma only
    Lattice<iScalar<vInteger> > x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > lin_z(U.Grid()); lin_z=x+y;
        
    ComplexField phases(U.Grid());
    int shift;
    ColourVector vec;
    FermionField temp(U.Grid());
    FermionField temp2(U.Grid());
    Coordinate site;
    Coordinate sparseSite;
    
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<uint32_t> uid(0, nt-1);
    
    for (int mu=0;mu<3;mu++){
        
        phases=1.0;
        if(mu==1){
            phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
        } else if(mu==2){
            phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
        }
        LatticeColourMatrix Umu(U.Grid());
        Umu = PeekIndex<LorentzIndex>(U,mu);
        Umu *= phases;
        
        for (unsigned int il = 0; il < 2*Nl_; il++)
        {
            // eval of unpreconditioned Dirac op
            std::complex<double> eval(mass,sqrt(epack.eval[il/2]-mass*mass));
            
            startTimer("V low mode");
            LOG(Message) << "V vector i = " << il << " (low modes)" << std::endl;
            a2a.makeLowModeV(temp, epack.evec[il/2], eval, il%2);
            stopTimer("V low mode");
            temp2 = Umu*Cshift(temp, mu, 1);
            temp *= conjugate(eval);// w vector, no eval in denom
            
            // Sparsen
            for(int t=0; t<nt;t+=par().tinc){
                int xshift=uid(gen);
                int yshift=uid(gen);
                int zshift=uid(gen);
                for(int z=0; z<ns;z+=par().inc){
                    for(int y=0; y<ns;y+=par().inc){
                        for(int x=0; x<ns;x+=par().inc){

                            site[0]=(x+xshift+ns)%ns;
                            site[1]=(y+yshift+ns)%ns;
                            site[2]=(z+zshift+ns)%ns;
                            site[3]=t;
                            sparseSite[0]=x/par().inc;
                            sparseSite[1]=y/par().inc;
                            sparseSite[2]=z/par().inc;
                            sparseSite[3]=t/par().tinc;
                            if(mu==0){// do v once
                                peekSite(vec,temp,site);
                                pokeSite(vec,v[il],sparseSite);
                            }
                            peekSite(vec,temp2,site);
                            pokeSite(vec,w[il],sparseSite);
                        }
                    }
                }
            }
        }
        // I/O
        startTimer("W I/O");
        A2AVectorsIo::write(par().output + "_w_mu" + std::to_string(mu), w,                         par().multiFile, vm().getTrajectory());
        stopTimer("W I/O");
    }// end mu
    startTimer("W I/O");
    A2AVectorsIo::write(par().output + "_v", v,par().multiFile, vm().getTrajectory());
    stopTimer("W I/O");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectors_hpp_
