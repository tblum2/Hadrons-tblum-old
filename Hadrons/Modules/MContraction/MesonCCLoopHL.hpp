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
                                    std::string, solver5D,
                                    std::string, action,
                                    std::string, action5D,
                                    double, mass,
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
    std::vector<std::string> input = {par().gauge, par().action, par().action5D};
    std::string sub_string;
    
    input.push_back(par().eigenPack);
    //input.push_back(par().solver + "_subtract");
    input.push_back(par().solver);
    input.push_back(par().solver5D);
    
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
    // grid can't handle real * prop, so use complex
    envTmpLat(LatticeComplex,  "herm_phase");
    envGetTmp(LatticeComplex, herm_phase);

    // sink
    Lattice<iScalar<vInteger> > x(env().getGrid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(env().getGrid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(env().getGrid()); LatticeCoordinate(z,2);
    Lattice<iScalar<vInteger> > t(env().getGrid()); LatticeCoordinate(t,3);
    Lattice<iScalar<vInteger> > s(env().getGrid());

    //``Hermiticity" phase, (-1)^(x+y)
    // sink only for now
    herm_phase = 1.0;
    s=x+y+z+t;
    herm_phase = where( mod(s,2)==(Integer)0, herm_phase, -herm_phase);
    //printMem("MesonLoopCCHL setup() end", env().getGrid()->ThisRank());
    
    auto        &action     = envGet(FMat, par().action);
    //auto        &solver     = envGet(Solver, par().solver + "_subtract");
    auto        &solver     = envGet(Solver, par().solver);
    auto        &solver5D     = envGet(Solver, par().solver5D);
    envTmp(A2A, "a2a", 1, action, solver);
    
    auto Ls = env().getObjectLs(par().action5D);
    envTmpLat(LatticeComplex, "corr");
    envTmpLat(FermionField, "source",Ls);
    envTmpLat(FermionField, "sourceshift",Ls);
    envTmpLat(FermionField, "tmp",Ls);
    envTmpLat(FermionField, "tmp2",Ls);
    envTmpLat(FermionField, "sol",Ls);
    envTmpLat(FermionField, "solshift",Ls);
    envTmpLat(FermionField, "sol2",Ls);
    envTmpLat(FermionField, "v");
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonLoopCCHL<FImpl1, FImpl2>::execute(void)
{
    LOG(Message) << "Computing Conserved Current Stag meson contractions " << std::endl;

    //printMem("MesonLoopCCHL execute() ", env().getGrid()->ThisRank());

    std::vector<TComplex>  buf;
    Result    result;
    int nt = env().getDim(Tp);
    int ns = env().getDim(Xp);
    
    result.corr.resize(nt);

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
    Guesser LMA(epack.evec, epack.eval);
    
    envGetTmp(LatticeComplex, corr);
    envGetTmp(LatticeComplex, herm_phase);
    envGetTmp(A2A, a2a);
    
    // Do spatial and temporal gammas only
    Lattice<iScalar<vInteger> > x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(U.Grid()); LatticeCoordinate(z,2);
    Lattice<iScalar<vInteger> > t(U.Grid()); LatticeCoordinate(t,3);
    Lattice<iScalar<vInteger> > lin_z(U.Grid()); lin_z=x+y;
    Lattice<iScalar<vInteger> > lin_t(U.Grid()); lin_z=x+y+z;
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
    envGetTmp(FermionField, sol2);
    envGetTmp(FermionField, v);
    Lattice<iScalar<vInteger> > t5d(source.Grid()); LatticeCoordinate(t5d,4);
    FermionField tmp_e(env().getRbGrid());
    FermionField tmp_o(env().getRbGrid());
    FermionField tmpRB(env().getRbGrid());
    FermionField tmpRB2(env().getRbGrid());
    FermionField sol4d(env().getGrid());
    
    Coordinate srcSite;
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

    tmp.Grid()->show_decomposition();
    Umu[0].Grid()->show_decomposition();
    assert(tmp.Grid()==Umu[0].Grid());

    int Nl_ = epack.evec.size();
    corr=Zero();
    for (unsigned int il = 0; il < Nl_; il+=Ls/2)
    {
        LOG(Message) << "Vector block " << il << std::endl;
        // prepare 5d source of evecs
        source = 0.0;
        sol = 0.0;
        int s=0;
        for(int iv=il;iv<il+Ls/2;iv++){
            // abs(eval) of unpreconditioned Dirac op
            std::complex<double> eval(mass,sqrt(epack.eval[iv]-mass*mass));
            // plus/minus evecs
            for(int pm=0;pm<2;pm++){
                // construct full lattice evec/eval as 4d source
                a2a.makeLowModeV(v, epack.evec[iv], eval, pm);
                InsertSlice(v,source,s,0);
                s++;
            }
        }
        // loop over source time slices
        for(int ts=0; ts<nt;ts+=par().tinc){
            
            LOG(Message) << "StagMesonLoopCCHLHL src_t " << ts << std::endl;


//                outFileName = par().output+"/cc_2pt_"+std::to_string(t)+"_mu_";
//                std::string file = resultFilename(outFileName+"0","h5");
//                bool f1 = exists(file);
//                file = resultFilename(outFileName+"1","h5");
//                bool f2 = exists(file);
//                file = resultFilename(outFileName+"2","h5");
//                bool f3 = exists(file);
//                if(f1 and f2 and f3){
//                    std::cout << "Skipping src point " << x << y << z << t << std::endl;
//                            continue;
//                }
            // source only on desired time slice
            // tmp = where((t5d == ts), source, source*0.);
            // now low-mode solution on source as initial guess sol
//            for(int s=0;s<Ls;s++){
//                ExtractSlice(sol4d,tmp,s,0);
///////////////////////////// from RedBlackSource in Grid: /////////////////////////////////////////////////////////////
//                pickCheckerboard(Even,tmp_e,sol4d);
//                pickCheckerboard(Odd ,tmp_o,sol4d);
//                /////////////////////////////////////////////////////
//                // src_o = (source_o - Moe MeeInv source_e)
//                /////////////////////////////////////////////////////
//                action.MooeeInv(tmp_e,tmpRB);     assert( tmpRB.Checkerboard() ==Even);
//                action.Meooe   (tmpRB,tmpRB2);    assert( tmpRB2.Checkerboard() ==Odd);
//                tmpRB2=tmp_o-tmpRB2;              assert( tmpRB2.Checkerboard() ==Odd);
//                action.Mooee(tmpRB2,tmp_o);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                // now hit with low-mode prop (guess)
//                LMA(tmp_o, tmpRB);
//                setCheckerboard(sol4d,tmpRB);
//                // zero out even part
//                tmpRB2.Checkerboard()=Even;
//                tmpRB2=0.;
//                setCheckerboard(sol4d,tmpRB2);
//                InsertSlice(sol4d,sol,s,0);
//            }
//            solver5D(sol, tmp);
            //FermToProp<FImpl1>(q1, sol, c);
        

                for(int mu=0;mu<3;mu++){

                    LOG(Message) << "StagMesonLoopCCHLHL src_mu " << mu << std::endl;
                    tmp = where((t5d == ts), source, source*0.);
                    tmp = adj(Umu[mu]) * tmp;
                    // shift source
                    tmp2 = Cshift(tmp, mu, 1);
                    // now low-mode solution on source as initial guess sol
                    for(int s=0;s<Ls;s++){
                        ExtractSlice(sol4d,tmp2,s,0);
        /////////////////////////// from RedBlackSource in Grid: /////////////////////////////////////////////////////////////
                        pickCheckerboard(Even,tmp_e,sol4d);
                        pickCheckerboard(Odd ,tmp_o,sol4d);
                        /////////////////////////////////////////////////////
                        // src_o = (source_o - Moe MeeInv source_e)
                        /////////////////////////////////////////////////////
                        action.MooeeInv(tmp_e,tmpRB);     assert( tmpRB.Checkerboard() ==Even);
                        action.Meooe   (tmpRB,tmpRB2);    assert( tmpRB2.Checkerboard() ==Odd);
                        tmpRB2=tmp_o-tmpRB2;              assert( tmpRB2.Checkerboard() ==Odd);
                        action.Mooee(tmpRB2,tmp_o);
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        // now hit with low-mode prop (guess)
                        LMA(tmp_o, tmpRB);
                        setCheckerboard(sol4d,tmpRB);
                        // zero out even part
                        tmpRB2.Checkerboard()=Even;
                        tmpRB2=0.;
                        setCheckerboard(sol4d,tmpRB2);
                        InsertSlice(sol4d,sol2,s,0);
                    }
                    solver5D(sol2, tmp2);
                    solshift=Cshift(sol2, mu, 1);
                    solshift = Umu[mu] * solshift;
                    corr += innerProduct(source, solshift);
                    sourceshift=Cshift(source, mu, 1);
                    sourceshift = Umu[mu] * sourceshift;
                    corr += innerProduct(sourceshift, sol2);
                    
                    tmp = where((t5d == ts), source, source*0.);
                    tmp2 = Cshift(tmp,mu,-1);
                    tmp = Umu[mu] * tmp2;
                    // now low-mode solution on source as initial guess sol
                    for(int s=0;s<Ls;s++){
                        ExtractSlice(sol4d,tmp,s,0);
        /////////////////////////// from RedBlackSource in Grid: /////////////////////////////////////////////////////////////
                        pickCheckerboard(Even,tmp_e,sol4d);
                        pickCheckerboard(Odd ,tmp_o,sol4d);
                        /////////////////////////////////////////////////////
                        // src_o = (source_o - Moe MeeInv source_e)
                        /////////////////////////////////////////////////////
                        action.MooeeInv(tmp_e,tmpRB);     assert( tmpRB.Checkerboard() ==Even);
                        action.Meooe   (tmpRB,tmpRB2);    assert( tmpRB2.Checkerboard() ==Odd);
                        tmpRB2=tmp_o-tmpRB2;              assert( tmpRB2.Checkerboard() ==Odd);
                        action.Mooee(tmpRB2,tmp_o);
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        // now hit with low-mode prop (guess)
                        LMA(tmp_o, tmpRB);
                        setCheckerboard(sol4d,tmpRB);
                        // zero out even part
                        tmpRB2.Checkerboard()=Even;
                        tmpRB2=0.;
                        setCheckerboard(sol4d,tmpRB2);
                        InsertSlice(sol4d,sol2,s,0);
                    }
                    solver5D(sol2, tmp);
                    solshift=Cshift(sol2, mu, 1);
                    solshift = Umu[mu] * solshift;
                    corr += innerProduct(source, solshift);
                    //sourceshift=Cshift(source, mu, 1);
                    //sourceshift = Umu[mu] * sourceshift;
                    corr += innerProduct(sourceshift, sol2);
                    
//                    //corr += trace(adj(q1) * Umu[mu] * qshift * adj(UmuSrc));
//                    corr += innerProduct(q1 * UmuSrc, Umu[mu] * qshift);
//
//                    qshift = Cshift(q1, mu, 1);
//                    // minus signs from herm transform of prop
//                    corr -= trace(adj(q2) * Umu[mu] * qshift * UmuSrc); // -1^muhat
//                    corr -= trace(adj(qshift) * adj(Umu[mu]) * q2 * adj(UmuSrc)); //-1^muhat
//
//                    corr *= herm_phase; // sink
//                    if((x+y+z+t)%2)corr *= -1.0; // source
//
//
//                    // do the local current
//                    corr = trace(adj(q1) * q1);
//                    // ks phases (includes hermiticity phase)
//                    corr *= localphases[mu]; // gamma at sink
//                    if(      mu==0 && x % 2 )corr *= -1.0; // gamma phase at source
//                    else if( mu==1 && y % 2 )corr *= -1.0; // gamma phase at source
//                    else if( mu==2 && z % 2 )corr *= -1.0; // gamma phase at source
//                    else if( mu==3 && t % 2 )corr *= -1.0; // gamma phase at source
//
//                    sliceSum(corr, buf, Tp);
//                    for (unsigned int tsnk = 0; tsnk < buf.size(); ++tsnk){
//                        result.corr[tsnk] = TensorRemove(buf[tsnk]);
//                    }
//                    outFileName = par().output+"local_2pt_"+std::to_string(t)+"_mu_"+std::to_string(mu);
//                    saveResult(outFileName, "mesonLL", result);
                }
//                // do the local Goldstone pion
//                corr = trace(adj(q1) * q1);
//                sliceSum(corr, buf, Tp);
//                for (unsigned int tsnk = 0; tsnk < buf.size(); ++tsnk){
//                    result.corr[tsnk] = TensorRemove(buf[tsnk]);
//                }
//                outFileName = par().output+"local_pion_"+
//                          std::to_string(x)+"_"+
//                          std::to_string(y)+"_"+
//                          std::to_string(z)+"_"+
//                          std::to_string(t);
//                saveResult(outFileName, "mesonLL", result);
        }
    }
    sliceSum(corr, buf, Tp);
    for (unsigned int tsnk = 0; tsnk < buf.size(); ++tsnk){
        result.corr[tsnk] = TensorRemove(buf[tsnk]);
    }
    outFileName = par().output+"HLcc_2pt";
    saveResult(outFileName, "HLCC", result);
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_MesonLoopCCHL_hpp_
