/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/MesonCCLoopHL.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>
Author: Vera Guelpers <vmg1n14@soton.ac.uk>
Author: Tom Blum (conserved currents)

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

#ifndef Hadrons_MContraction_MesonLoopCCHL5D_hpp_
#define Hadrons_MContraction_MesonLoopCCHL5D_hpp_

#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MSource/Point.hpp>
#include <Hadrons/Solver.hpp>
#include <Grid/lattice/Lattice_reduction.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *             TMesonLoopCCHL                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)


class MesonLoopCCHL5DPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonLoopCCHL5DPar,
                                    std::string, gauge,
                                    std::string, output,
                                    std::string, eigenPack,
                                    std::string, solver,
                                    std::string, solver5D,
                                    std::string, action,
                                    std::string, action5D,
                                    double, mass,
                                    int, tinc);
};

template <typename FImpl1, typename FImpl2>
class TStagMesonLoopCCHL5D: public Module<MesonLoopCCHL5DPar>
{
public:
    typedef typename FImpl1::FermionField FermionField;
    typedef A2AVectorsSchurStaggered<FImpl1> A2A;
    typedef FermionOperator<FImpl1>          FMat;
    FERM_TYPE_ALIASES(FImpl1,1);
    FERM_TYPE_ALIASES(FImpl2,2);
    SOLVER_TYPE_ALIASES(FImpl1,);
    //BASIC_TYPE_ALIASES(FImpl1,);
    //BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    //SINK_TYPE_ALIASES(Scalar);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TStagMesonLoopCCHL5D(const std::string name);
    // destructor
    virtual ~TStagMesonLoopCCHL5D(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    protected:
    // execution
    virtual void setup(void);
    // execution
    virtual void execute(void);
    inline bool exists (const std::string& name) {
      struct stat buffer;
      return (stat (name.c_str(), &buffer) == 0);
    }
private:
    //FMat         *action_{nullptr};
    //Solver       *solver_{nullptr};
};

MODULE_REGISTER_TMP(StagMesonLoopCCHL5D, ARG(TStagMesonLoopCCHL5D<STAGIMPL, STAGIMPL>), MContraction);

/******************************************************************************
 *                           TStagMesonLoopCCHL5D implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TStagMesonLoopCCHL5D<FImpl1, FImpl2>::TStagMesonLoopCCHL5D(const std::string name)
: Module<MesonLoopCCHL5DPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonLoopCCHL5D<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().gauge, par().action, par().action5D};
    std::string sub_string;
    
    input.push_back(par().eigenPack);
    //input.push_back(par().solver + "_subtract");
    input.push_back(par().solver);
    input.push_back(par().solver5D);
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonLoopCCHL5D<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {};

    return output;
}

// setup ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonLoopCCHL5D<FImpl1, FImpl2>::setup(void)
{
    
    auto        &action     = envGet(FMat, par().action);
    //auto        &solver     = envGet(Solver, par().solver + "_subtract");
    auto        &solver     = envGet(Solver, par().solver);
    auto        &solver5D     = envGet(Solver, par().solver5D);
    envTmp(A2A, "a2a", 1, action, solver);
    
    auto Ls = env().getObjectLs(par().action5D);
    envTmpLat(FermionField, "source",Ls);
    envTmpLat(FermionField, "sourceshift",Ls);
    envTmpLat(FermionField, "tmp",Ls);
    envTmpLat(FermionField, "tmp2",Ls);
    envTmpLat(FermionField, "sol",Ls);
    envTmpLat(FermionField, "solshift",Ls);
    envTmpLat(FermionField, "w");
    envTmpLat(FermionField, "w2");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonLoopCCHL5D<FImpl1, FImpl2>::execute(void)
{
    LOG(Message) << "Computing Conserved Current Stag meson contractions " << std::endl;

    //printMem("MesonLoopCCHL execute() ", env().getGrid()->ThisRank());

    std::vector<TComplex>  buf;
    std::vector<Result> result(3);
    std::vector<ComplexD>  corr;

    int nt = env().getDim(Tp);
    int ns = env().getDim(Xp);
    
    
    
    for(int mu=0;mu<3;mu++){
        result[mu].corr.resize(nt);
        for(int t=0;t<nt;t++){
            result[mu].corr[t]=0.;
        }
    }

    auto &U       = envGet(LatticeGaugeField, par().gauge);
    auto        &action      = envGet(FMat, par().action);
    auto        &action5D    = envGet(FMat, par().action5D);
    //auto        &solver    = envGet(Solver, par().solver + "_subtract");
    auto        &solver5D    = envGet(Solver, par().solver5D);
    auto &epack   = envGet(BaseFermionEigenPack<FImpl1>, par().eigenPack);
    double mass = par().mass;
    auto Ls = env().getObjectLs(par().action5D);
    LOG(Message) << "Ls=" << Ls << std::endl;
    typedef DeflatedGuesser<FermionField>  Guesser;
    // for solution
    Guesser LMA(epack.evec, epack.eval);
    std::vector<double> mlsq(epack.eval.size());
    for(int i=0;i<epack.eval.size();i++){
        mlsq[i]=(epack.eval[i]-mass*mass) * mass;
    }
    // for subtraction
    Guesser LLsub(epack.evec, mlsq);
    
    envGetTmp(A2A, a2a);
    
    // Do spatial and temporal gammas only
    Lattice<iScalar<vInteger> > x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(U.Grid()); LatticeCoordinate(z,2);
    Lattice<iScalar<vInteger> > t(U.Grid()); LatticeCoordinate(t,3);
    Lattice<iScalar<vInteger> > t5d(env().getGrid(Ls));LatticeCoordinate(t5d,4);

    Lattice<iScalar<vInteger> > lin_z(U.Grid()); lin_z=x+y;
    Lattice<iScalar<vInteger> > lin_t(U.Grid()); lin_t=x+y+z;
    LatticeComplex phases(U.Grid());
    std::vector<LatticeComplex> localphases(3,U.Grid());
    std::vector<LatticeColourMatrix> Umu(3,U.Grid());
    // 5d source, solution
    envGetTmp(FermionField, source);
    envGetTmp(FermionField, sourceshift);
    envGetTmp(FermionField, tmp);
    envGetTmp(FermionField, tmp2);
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, solshift);
    envGetTmp(FermionField, w);
    envGetTmp(FermionField, w2);
    FermionField tmp_e(env().getRbGrid());
    FermionField tmp_o(env().getRbGrid());
    FermionField tmpRB(env().getRbGrid());
    FermionField tmpRB2(env().getRbGrid());
    FermionField sol4d(env().getGrid());
    FermionField sub(env().getGrid());
    
    std::string outFileName;
    
    for(int mu=0;mu<3;mu++){

        //staggered phases go into links
        Umu[mu] = PeekIndex<LorentzIndex>(U,mu);
        localphases[mu]=1.0;
        phases=1.0;
        if(mu==0){
            localphases[0] = where( mod(x    ,2)==(Integer)0, localphases[0],-localphases[0]);
        }else if(mu==1){
            localphases[1] = where( mod(y    ,2)==(Integer)0, localphases[1],-localphases[1]);
            phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
        }else if(mu==2){
            localphases[2] = where( mod(z    ,2)==(Integer)0, localphases[2],-localphases[2]);
            phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
        }else if(mu==3){
            localphases[3] = where( mod(t    ,2)==(Integer)0, localphases[3],-localphases[3]);
            phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);
        }else assert(0);
        Umu[mu] *= phases;
    }

    int Nl_ = epack.evec.size();
    for (unsigned int il = 0; il < Nl_; il+=Ls/2)
    {
        LOG(Message) << "Vector block " << il << std::endl;
        // prepare 5d source of evecs
        int s=0;
        for(int iv=il;iv<il+Ls/2;iv++){
            // abs(eval) of unpreconditioned Dirac op
            std::complex<double> eval(mass,sqrt(epack.eval[iv]-mass*mass));
            // both plus/minus evecs
            for(int pm=0;pm<2;pm++){
                // construct full lattice evec/eval as 4d source
                a2a.makeLowModeW(w, epack.evec[iv], eval, pm);
                InsertSlice(w,source,s,0);
                s++;
            }
        }
        // loop over source time slices
        for(int ts=0; ts<nt;ts+=par().tinc){
            
            LOG(Message) << "StagMesonLoopCCHL5D src_t " << ts << std::endl;

            // note that for 5d vecs, mu=x,y,z are 1,2,3
            for(int mu=0;mu<3;mu++){

                LOG(Message) << "StagMesonLoopCCHL5D src_mu " << mu << std::endl;
                tmp = where((t5d == ts), source, source*0.);
                for(int s=0;s<Ls;s++){
                    ExtractSlice(w,tmp,s,0);
                    w2 = adj(Umu[mu]) * w;
                    InsertSlice(w2,tmp,s,0);
                }
                // shift source at x to x+mu
                tmp2 = Cshift(tmp, mu+1, -1);
                // now low-mode solution on source as initial guess sol
                for(int s=0;s<Ls;s++){
                    ExtractSlice(sol4d,tmp2,s,0);
                    /////////////////////////// from RedBlackSource in Grid: //////////////////////////////////////
                    pickCheckerboard(Even,tmp_e,sol4d);
                    pickCheckerboard(Odd ,tmp_o,sol4d);
                    /////////////////////////////////////////////////////
                    // src_o = (source_o - Moe MeeInv source_e)
                    /////////////////////////////////////////////////////
                    action.MooeeInv(tmp_e,tmpRB);  assert( tmpRB.Checkerboard() ==Even);
                    action.Meooe   (tmpRB,tmpRB2); assert( tmpRB2.Checkerboard() ==Odd);
                    tmpRB2=tmp_o-tmpRB2;           assert( tmpRB2.Checkerboard() ==Odd);
                    action.Mooee(tmpRB2,tmp_o);
                    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // now hit with low-mode prop (guess)
                    LMA(tmp_o, tmpRB);
                    setCheckerboard(sol4d,tmpRB);
                    // zero out even part
                    tmpRB2.Checkerboard()=Even;
                    tmpRB2=Zero();
                    setCheckerboard(sol4d,tmpRB2);
                    InsertSlice(sol4d,sol,s,0);
                }
                solver5D(sol, tmp2);
                // subtract the low modes
                for(int s=0;s<Ls;s++){
                    ExtractSlice(sol4d,sol,s,0);
                    ExtractSlice(w,tmp2,s,0);
                    sub = Zero();
                    pickCheckerboard(Even,tmp_e,w);
                    action.Meooe(tmp_e,tmp_o);
                    LLsub(tmp_o,tmp_e);
                    action.Meooe(tmp_e,tmp_o);// tmp_o is now even
                    setCheckerboard(sub,tmp_o);
                    sol4d += sub;
                    InsertSlice(sol4d,sol,s,0);
                }
                // take inner-product with eigen bra on all time slices
                for(int s=0;s<Ls;s++){
                    ExtractSlice(w,source,s,0);
                    w2=Cshift(w, mu, 1);
                    w = Umu[mu] * w2;
                    std::complex<double> eval(mass,sqrt(epack.eval[s/2+il]-mass*mass));
                    // use conjugate since w conjugated below
                    if(s%2==0)eval=conjugate(eval);
                    Complex cc = 1.0/eval;
                    w2 = cc*w;
                    InsertSlice(w2,tmp,s,0);
                }
                sliceInnerProductVector(corr,tmp,sol,4);
                for(int tsnk=0; tsnk<nt; tsnk++){
                    result[mu].corr[(tsnk-ts+nt)%nt] += corr[tsnk];
                }
                
                solshift=Cshift(sol, mu+1, 1);
                // take inner-product with eigen bra on all time slices
                for(int s=0;s<Ls;s++){
                    ExtractSlice(w,source,s,0);
                    w2 = adj(Umu[mu]) * w;
                    std::complex<double> eval(mass,sqrt(epack.eval[s/2+il]-mass*mass));
                    // use conjugate since w conjugated below
                    if(s%2==0)eval=conjugate(eval);
                    Complex cc = 1.0/eval;
                    w = cc*w2;
                    InsertSlice(w,tmp,s,0);
                }
                sliceInnerProductVector(corr,tmp,solshift,4);
                for(int tsnk=0; tsnk<nt; tsnk++){
                    result[mu].corr[(tsnk-ts+nt)%nt] += corr[tsnk];
                }

                tmp = where((t5d == ts), source, source*0.);
                // shift source
                tmp2 = Cshift(tmp, mu+1, 1);
                for(int s=0;s<Ls;s++){
                    ExtractSlice(w,tmp2,s,0);
                    w2 = Umu[mu] * w;
                    InsertSlice(w2,tmp2,s,0);
                }
                // now low-mode solution on source as initial guess sol
                for(int s=0;s<Ls;s++){
                    ExtractSlice(sol4d,tmp2,s,0);
                    /////////////////////////// from RedBlackSource in Grid: ///////////////////////////////////////////
                    pickCheckerboard(Even,tmp_e,sol4d);
                    pickCheckerboard(Odd ,tmp_o,sol4d);
                    /////////////////////////////////////////////////////
                    // src_o = (source_o - Moe MeeInv source_e)
                    /////////////////////////////////////////////////////
                    action.MooeeInv(tmp_e,tmpRB);  assert( tmpRB.Checkerboard() ==Even);
                    action.Meooe   (tmpRB,tmpRB2); assert( tmpRB2.Checkerboard() ==Odd);
                    tmpRB2=tmp_o-tmpRB2;           assert( tmpRB2.Checkerboard() ==Odd);
                    action.Mooee(tmpRB2,tmp_o);
                    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // now hit with low-mode prop (guess)
                    LMA(tmp_o, tmpRB);
                    setCheckerboard(sol4d,tmpRB);
                    // zero out even part
                    tmpRB2.Checkerboard()=Even;
                    tmpRB2=0.;
                    setCheckerboard(sol4d,tmpRB2);
                    InsertSlice(sol4d,sol,s,0);
                }
                
                solver5D(sol, tmp2);
                // subtract the low modes
                for(int s=0;s<Ls;s++){
                    ExtractSlice(sol4d,sol,s,0);
                    ExtractSlice(w,tmp2,s,0);
                    sub = Zero();
                    pickCheckerboard(Even,tmp_e,w);
                    action.Meooe(tmp_e,tmp_o);
                    LLsub(tmp_o,tmp_e);
                    action.Meooe(tmp_e,tmp_o);// tmp_o is now even
                    setCheckerboard(sub,tmp_o);
                    sol4d += sub;
                    InsertSlice(sol4d,sol,s,0);
                }
                
                // take inner-product with eigen bra on all time slices
                for(int s=0;s<Ls;s++){
                    ExtractSlice(w,source,s,0);
                    w2=Cshift(w, mu, 1);
                    w = Umu[mu] * w2;
                    std::complex<double> eval(mass,sqrt(epack.eval[s/2+il]-mass*mass));
                    // use conjugate since w conjugated below
                    if(s%2==0)eval=conjugate(eval);
                    Complex cc = 1.0/eval;
                    w2 = w*cc;
                    InsertSlice(w2,tmp,s,0);
                }
                sliceInnerProductVector(corr,tmp,sol,4);
                for(int tsnk=0; tsnk<nt; tsnk++){
                    result[mu].corr[(tsnk-ts+nt)%nt] += corr[tsnk];
                }
                
                solshift=Cshift(sol, mu+1, 1);
                // take inner-product with eigen bra on all time slices
                for(int s=0;s<Ls;s++){
                    ExtractSlice(w,source,s,0);
                    w2 = adj(Umu[mu]) * w;
                    std::complex<double> eval(mass,sqrt(epack.eval[s/2+il]-mass*mass));
                    // use conjugate since w conjugated below
                    if(s%2==0)eval=conjugate(eval);
                    Complex cc = 1.0/eval;
                    w = w2*cc;
                    InsertSlice(w,tmp,s,0);
                }
                sliceInnerProductVector(corr,tmp,solshift,4);
                for(int tsnk=0; tsnk<nt; tsnk++){
                    result[mu].corr[(tsnk-ts+nt)%nt] += corr[tsnk];
                }
            }
        }
    }
    
    for (int i = 0; i < 3; ++i){
        if(U.Grid()->IsBoss()){
            makeFileDir(par().output);
            outFileName = par().output+"HLcc_2pt_mu"+std::to_string(i);
            saveResult(outFileName, "HLCC", result[i]);
        }
    }
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif
