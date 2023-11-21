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

#ifndef Hadrons_MContraction_MesonLoopCCHL_hpp_
#define Hadrons_MContraction_MesonLoopCCHL_hpp_

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


class MesonLoopCCHLPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonLoopCCHLPar,
                                    std::string, gauge,
                                    std::string, output,
                                    std::string, eigenPack,
                                    std::string, solver,
                                    std::string, action,
                                    double, mass,
                                    int, numHits,
                                    int, blockSize,
                                    int, tinc);
};

template <typename FImpl1, typename FImpl2>
class TStagMesonLoopCCHL: public Module<MesonLoopCCHLPar>
{
public:
    typedef typename FImpl1::FermionField FermionField;
    typedef A2AVectorsSchurStaggered<FImpl1> A2A;
    typedef FermionOperator<FImpl1>          FMat;
    FERM_TYPE_ALIASES(FImpl1,1);
    FERM_TYPE_ALIASES(FImpl2,2);
    SOLVER_TYPE_ALIASES(FImpl1,);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TStagMesonLoopCCHL(const std::string name);
    // destructor
    virtual ~TStagMesonLoopCCHL(void) {};
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

MODULE_REGISTER_TMP(StagMesonLoopCCHL, ARG(TStagMesonLoopCCHL<STAGIMPL, STAGIMPL>), MContraction);

/******************************************************************************
 *                           TStagMesonLoopCCHL implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TStagMesonLoopCCHL<FImpl1, FImpl2>::TStagMesonLoopCCHL(const std::string name)
: Module<MesonLoopCCHLPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonLoopCCHL<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().gauge, par().action};
    std::string sub_string;
    
    input.push_back(par().eigenPack);
    //input.push_back(par().solver + "_subtract");
    input.push_back(par().solver);
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonLoopCCHL<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {};

    return output;
}

// setup ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonLoopCCHL<FImpl1, FImpl2>::setup(void)
{
    
    auto        &action     = envGet(FMat, par().action);
    auto        &solver     = envGet(Solver, par().solver);
    envTmp(A2A, "a2a", 1, action, solver);
    
    envTmpLat(FermionField, "srcb");
    envTmpLat(FermionField, "snkb");
    envTmpLat(FermionField, "sub");
    envTmpLat(FermionField, "tmp");
    envTmpLat(FermionField, "tmp2");
    envTmpLat(FermionField, "sol");
    envTmpLat(FermionField, "solshift");
    envTmpLat(FermionField, "w");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonLoopCCHL<FImpl1, FImpl2>::execute(void)
{
    LOG(Message) << "Computing High-Low Conserved Current Stag meson contractions " << std::endl;

    std::vector<ComplexD>  corr;
    std::vector<Result> result(3);
    int nt = env().getDim(Tp);
    int ns = env().getDim(Xp);
    // init
    for(int mu=0;mu<3;mu++){
        result[mu].corr.resize(nt);
        for(int t=0;t<nt;t++){
            result[mu].corr[t]=(ComplexD)(0.,0.);
        }
    }

    auto &U       = envGet(LatticeGaugeField, par().gauge);
    auto        &action      = envGet(FMat, par().action);
    //auto        &solver    = envGet(Solver, par().solver + "_subtract");
    auto        &solver    = envGet(Solver, par().solver);
    auto &epack   = envGet(BaseFermionEigenPack<FImpl1>, par().eigenPack);
    double mass = par().mass;
    std::vector<double> mlsq(epack.eval.size());
    for(int i=0;i<epack.eval.size();i++){
        mlsq[i]=(epack.eval[i]-mass*mass) * mass;
    }
    DeflatedGuesser<FermionField> LLsub(epack.evec, mlsq);
    FermionField tmp_e(env().getRbGrid());
    FermionField tmp_o(env().getRbGrid());

    envGetTmp(A2A, a2a);
    
    // Do spatial gammas only
    Lattice<iScalar<vInteger>> x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger>> y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger>> z(U.Grid()); LatticeCoordinate(z,2);
    Lattice<iScalar<vInteger>> t(U.Grid()); LatticeCoordinate(t,3);
    Lattice<iScalar<vInteger>> lin_z(U.Grid()); lin_z=x+y;
    Lattice<iScalar<vInteger>> lin_t(U.Grid()); lin_t=x+y+z;
    LatticeComplex phases(U.Grid());
    std::vector<LatticeColourMatrix> Umu(3,U.Grid());
    // source, solution
    envGetTmp(FermionField, srcb);
    envGetTmp(FermionField, snkb);
    envGetTmp(FermionField, tmp);
    envGetTmp(FermionField, tmp2);
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, solshift);
    envGetTmp(FermionField, sub);
    envGetTmp(FermionField, w);
    
    std::string outFileName;
    int Nl_ = epack.evec.size();

    // save randoms so exact and sub use same
    std::vector<Complex> eta(nt*3*par().numHits*Nl_*2);
    int numb=2*Nl_/par().blockSize;
    for(int ts=0; ts<nt;ts++){
        for(int mu=0;mu<3;mu++){
            for(int ih=0;ih<par().numHits;ih++){
                for (unsigned int ib = 0; ib < numb; ib++){
                    
                    for(int il=0;il<par().blockSize;il++){
                        int idx=il+par().blockSize*(ib+numb*(ih+par().numHits*(mu+3*ts)));
                        bernoulli(rngSerial(), eta[idx]);
                        Complex shift(1., 1.);
                        eta[idx] = (2.*eta[idx] - shift)*(1./::sqrt(2.));
//                        LOG(Message) << "il " << il << std::endl;
//                        LOG(Message) << "ib " << ib << std::endl;
//                        LOG(Message) << "ih " << ih << std::endl;
//                        LOG(Message) << "mu " << mu << std::endl;
//                        LOG(Message) << "ts " << ts << std::endl;
//                        LOG(Message) << "eta idx " << idx << std::endl;
                    }
                }
            }
        }
    }
    for(int mu=0;mu<3;mu++){

        //staggered phases go into links
        Umu[mu] = PeekIndex<LorentzIndex>(U,mu);
        phases=1.0;
        if(mu==0){
        }else if(mu==1){
            phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
        }else if(mu==2){
            phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
        }else if(mu==3){
            phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);
        }else assert(0);
        Umu[mu] *= phases;
    }

    // loop over source time slices
    for(int mu=0;mu<3;mu++){

        LOG(Message) << "StagMesonLoopCCHL src_mu " << mu << std::endl;
        for(int ts=0; ts<nt;ts+=par().tinc){
        
            LOG(Message) << "StagMesonLoopCCHL src_t " << ts << std::endl;

            // loop over hits
            for(int ih=0;ih<par().numHits;ih++){
                
                LOG(Message) << "StagMesonLoopCCHL hit " << ih << std::endl;
                
                // loop over blocks of evecs
                for (unsigned int ib = 0; ib < numb; ib++){
                    
                    // sum evecs over block
                    srcb = 0.;
                    snkb = 0.;
                    for(int il=0;il<par().blockSize;il++){
                        
                        int ivec = il+ib*par().blockSize;
                        int idx=il+par().blockSize*(ib+numb*(ih+par().numHits*(mu+3*ts)));
                        eta[idx] /= pow(epack.eval[ivec/2], 0.25);
                    
                        std::complex<double> eval(mass,sqrt(epack.eval[ivec/2]-mass*mass));
                        double arg=std::arg(eval);
                        std::complex<double> phase(cos(arg),sin(arg));
                        
                        // do plus/minus evecs
                        int pm = ivec%2;
                        
                        // construct full lattice evec as 4d source (no 1/lambda here)
                        a2a.makeLowModeW(w, epack.evec[ivec/2], eval, pm);
                        
                        // -lambda eigenvalue
                        if(pm){
                            phase=conjugate(phase);
                        }
                        // source and sink block vectors
                        w *= eta[idx];
                        snkb += w;
                        w *= 1./phase;
                        srcb += w;
                    }
                    // split grid here?
                    
                    tmp = where(t == ts, srcb, srcb*0.);
                    tmp2 = adj(Umu[mu]) * tmp;
                   
                    // shift source at x to x+mu
                    tmp = Cshift(tmp2, mu, -1);
                    
                    solver(sol, tmp);
                    // subtract the low modes
                    sub = Zero();
                    pickCheckerboard(Even,tmp_e,tmp);
                    action.Meooe(tmp_e,tmp_o);
                    LLsub(tmp_o,tmp_e);
                    action.Meooe(tmp_e,tmp_o);// tmp_o is now even
                    setCheckerboard(sub,tmp_o);
                    sol += sub;
                    // take inner-product with eigenbra on all time slices
                    tmp = Cshift(snkb, mu, 1);
                    tmp2 = Umu[mu] * tmp;
                    sliceInnerProductVector(corr,tmp2,sol,3);
                    for(int tsnk=0; tsnk<nt; tsnk++){
                        ComplexD cc = corr[tsnk];
                        result[mu].corr[(tsnk-ts+nt)%nt] += cc;
                    }
                    
                    solshift=Cshift(sol, mu, 1);
                    // take inner-product with eigenbra on all time slices
                    tmp = adj(Umu[mu]) * snkb;
                    sliceInnerProductVector(corr,tmp,solshift,3);
                    for(int tsnk=0; tsnk<nt; tsnk++){
                        ComplexD cc = corr[tsnk];
                        result[mu].corr[(tsnk-ts+nt)%nt] += cc;
                    }

                    tmp = where(t == ts, srcb, srcb*0.);
                    // shift source
                    tmp2 = Cshift(tmp, mu, 1);
                    tmp = Umu[mu] * tmp2;

                    solver(sol, tmp);
                    sub = Zero();
                    pickCheckerboard(Even,tmp_e,tmp);
                    action.Meooe(tmp_e,tmp_o);
                    LLsub(tmp_o,tmp_e);
                    action.Meooe(tmp_e,tmp_o);// tmp_o is now even
                    setCheckerboard(sub,tmp_o);
                    sol += sub;
                    
                    // take inner-product with eigenmode on all time slices
                    tmp = Cshift(snkb, mu, 1);
                    tmp2 = Umu[mu] * tmp;
                    sliceInnerProductVector(corr,tmp2,sol,3);
                    for(int tsnk=0; tsnk<nt; tsnk++){
                        ComplexD cc = corr[tsnk];
                        result[mu].corr[(tsnk-ts+nt)%nt] += cc;
                    }

                    solshift=Cshift(sol, mu, 1);
                    // take inner-product with eigenmode on all time slices
                    tmp = adj(Umu[mu]) * snkb;
                    sliceInnerProductVector(corr,tmp,solshift,3);
                    for(int tsnk=0; tsnk<nt; tsnk++){
                        ComplexD cc = corr[tsnk];
                        result[mu].corr[(tsnk-ts+nt)%nt] += cc;
                    }
                }
            }
        }
        if(U.Grid()->IsBoss()){
            makeFileDir(par().output);
            outFileName = par().output+"HLcc_2pt_mu"+std::to_string(mu);
            for(int t=0; t<nt; t++)
                result[mu].corr[t] /= par().numHits;
            saveResult(outFileName, "HLCC", result[mu]);
        }
    }
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif
