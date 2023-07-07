/*
 * RBPrecCG.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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

#ifndef Hadrons_MSolver_RBPrecBCG_hpp_
#define Hadrons_MSolver_RBPrecBCG_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MSolver/Guesser.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     Schur red-black preconditioned BCG                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class RBPrecBCGPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RBPrecBCGPar ,
                                    std::string , action,
                                    unsigned int, maxIteration,
                                    double      , residual,
                                    bool, solInitGuess,
                                    std::string , eigenPack);
};

template <typename FImpl, int nBasis>
class TRBPrecBCG: public Module<RBPrecBCGPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TRBPrecBCG(const std::string name);
    // destructor
    virtual ~TRBPrecBCG(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RBPrecBCG, ARG(TRBPrecBCG<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
MODULE_REGISTER_TMP(ZRBPrecBCG, ARG(TRBPrecBCG<ZFIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
MODULE_REGISTER_TMP(StagRBPrecBCG, ARG(TRBPrecBCG<STAGIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>), MSolver);
/******************************************************************************
 *                      TRBPrecBCG template implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
TRBPrecBCG<FImpl, nBasis>::TRBPrecBCG(const std::string name)
: Module(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
std::vector<std::string> TRBPrecBCG<FImpl, nBasis>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TRBPrecBCG<FImpl, nBasis>::getReference(void)
{
    std::vector<std::string> ref = {par().action};
    
    if (!par().eigenPack.empty())
    {
        ref.push_back(par().eigenPack);
    }

    return ref;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TRBPrecBCG<FImpl, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_subtract"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TRBPrecBCG<FImpl, nBasis>::setup(void)
{
    if (par().maxIteration == 0)
    {
        HADRONS_ERROR(Argument, "zero maximum iteration");
    }

    LOG(Message) << "setting up Schur red-black preconditioned BCG for"
                 << " action '" << par().action << "' with residual "
                 << par().residual << ", maximum iteration " 
                 << par().maxIteration << std::endl;

    auto Ls        = env().getObjectLs(par().action);
    auto &mat      = envGet(FMat, par().action);
    auto guesserPt = makeGuesser<FImpl, nBasis>(par().eigenPack);

    auto makeSolver = [&mat, guesserPt, this](bool subGuess) {
        return [&mat, guesserPt, subGuess, this](FermionField &sol,
                                     const FermionField &source) {
            BlockConjugateGradient<FermionField> bcgrq(BlockCGrQ,
                                                       0,par().residual,
                                                       par().maxIteration);
            
            //HADRONS_DEFAULT_SCHUR_SOLVE<FermionField> schurSolver(cg);
            //schurSolver.subtractGuess(subGuess);
            // use sol as initial guess:
            HADRONS_DEFAULT_SCHUR_SOLVE<FermionField> schurSolver(bcgrq,subGuess,par().solInitGuess);
            schurSolver(mat, source, sol, *guesserPt);
        };
    };
    auto solver = makeSolver(false);
    envCreate(Solver, getName(), Ls, solver, mat);
    auto solver_subtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls, solver_subtract, mat);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TRBPrecBCG<FImpl, nBasis>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_RBPrecBCG_hpp_
