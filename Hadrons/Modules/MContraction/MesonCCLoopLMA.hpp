/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/Meson.hpp

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

#ifndef Hadrons_MContraction_MesonCCLoopLMA_hpp_
#define Hadrons_MContraction_MesonCCLoopLMA_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MSource/Point.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 2pt conserved-conserved staggered correlation function
 -----------------------------
 
 * options:
 - q1: input propagator 1 (string) src at y
 - q2: input propagator 2 (string) src at y+hat mu
 
*/

/******************************************************************************
 *                                TMesonCCLoopLMA                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)


class MesonCCLoopLMAPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonCCLoopLMAPar,
                                    std::string, gauge,
                                    std::string, output,
                                    std::string, eigenPack,
                                    std::string, action,
                                    double, mass,
                                    int, inc,
                                    int, tinc);
};

template <typename FImpl1, typename FImpl2>
class TStagMesonCCLoopLMA: public Module<MesonCCLoopLMAPar>
{
public:
    typedef typename FImpl1::FermionField FermionField;
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    //SOLVER_TYPE_ALIASES(FImpl1,);
    BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    SINK_TYPE_ALIASES(Scalar);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TStagMesonCCLoopLMA(const std::string name);
    // destructor
    virtual ~TStagMesonCCLoopLMA(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    inline bool exists (const std::string& name) {
      struct stat buffer;
      return (stat (name.c_str(), &buffer) == 0);
    }
    // low mode approx
    void LMA(const FermionField &src, FermionField &guess)
    // copied original from DeflatedGuesser.h
    {
        
      /*guess = Zero();
      assert(evec.size()==eval.size());
      auto N = evec.size();
      for (int i=0;i<N;i++) {
        const Field& tmp = evec[i];
        axpy(guess,TensorRemove(innerProduct(tmp,src)) / eval[i],tmp,guess);
      }*/
      guess.Checkerboard() = src.Checkerboard();
    };
private:
    //Solver       *solver_{nullptr};
};

MODULE_REGISTER_TMP(StagMesonCCLoopLMA, ARG(TStagMesonCCLoopLMA<STAGIMPL, STAGIMPL>), MContraction);

/******************************************************************************
 *                           TStagMesonCCLoopLMA implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TStagMesonCCLoopLMA<FImpl1, FImpl2>::TStagMesonCCLoopLMA(const std::string name)
: Module<MesonCCLoopLMAPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonCCLoopLMA<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().gauge, par().eigenPack, par().action};
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonCCLoopLMA<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {};
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonCCLoopLMA<FImpl1, FImpl2>::setup(void)
{
    envTmpLat(LatticeComplex, "corr");
    envTmpLat(PropagatorField1, "q1");
    envTmpLat(PropagatorField2, "q2");
    envTmpLat(PropagatorField1, "qshift");
    envTmpLat(FermionField, "source");
    envTmpLat(FermionField, "sol");
    
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
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonCCLoopLMA<FImpl1, FImpl2>::execute(void)
{
    typedef DeflatedGuesser<FermionField>  Guesser;
    
    LOG(Message) << "Computing Conserved Current Stag meson contractions using LMA approx " << std::endl;
    
    std::vector<TComplex>  buf;
    Result    result;
    int nt = env().getDim(Tp);
    int ns = env().getDim(Xp);
    
    result.corr.resize(nt);
    
    auto &U = envGet(LatticeGaugeField, par().gauge);
    auto &epack = envGet(BaseFermionEigenPack<FImpl1>, par().eigenPack);
    auto &action = envGet(FermionOperator<FImpl1>, par().action); // for mult by Meo, Moe
    Guesser LMA(epack.evec, epack.eval);
    
    // need another guesser to mult by m/lambda^2;
    double mass = par().mass;
    std::vector<double> movevalsq(epack.eval.size());
    for(int i=0;i<epack.eval.size();i++){
        movevalsq[i]=(epack.eval[i]-mass*mass) * epack.eval[i] / mass;
    }
    Guesser LMA2(epack.evec, movevalsq);
    
    envGetTmp(LatticeComplex, corr);
    envGetTmp(LatticeComplex, herm_phase);
    envGetTmp(PropagatorField1, q1);
    envGetTmp(PropagatorField2, q2);
    envGetTmp(PropagatorField1, qshift);
    
    // Do spatial gamma's only
    Lattice<iScalar<vInteger> > x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(U.Grid()); LatticeCoordinate(z,2);
    Lattice<iScalar<vInteger> > lin_z(U.Grid()); lin_z=x+y;
    LatticeComplex phases(U.Grid());
    std::vector<LatticeComplex> localphases(3,U.Grid());

    std::vector<LatticeColourMatrix> Umu(3,U.Grid());
    envGetTmp(FermionField, source);
    envGetTmp(FermionField, sol);
    FermionField halfSource(env().getRbGrid());
    FermionField halfSol(env().getRbGrid());
    FermionField tmp(env().getRbGrid());

    Coordinate srcSite;
    ColourMatrix UmuSrc;
    ColourVector Csrc;
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
        }else assert(0);
        Umu[mu] *= phases;
    }
        
    // loop over source position
    // assumed to be Even for now
    for(int t=0; t<nt;t+=par().tinc){
        for(int z=0; z<ns;z+=par().inc){
            for(int y=0; y<ns;y+=par().inc){
                for(int x=0; x<ns;x+=par().inc){
                    
                    srcSite[0]=x;
                    srcSite[1]=y;
                    srcSite[2]=z;
                    srcSite[3]=t;
                    assert((x+y+z+t)%2==0);// must be Even
                    
                    outFileName = par().output+"/cc_2pt_"+
                        std::to_string(x)+"_"+
                        std::to_string(y)+"_"+
                        std::to_string(z)+"_"+
                        std::to_string(t)+"_mu_";
                    std::string file = resultFilename(outFileName+"0","h5");
                    //bool f1 = std::__fs::filesystem::exists(file);
                    bool f1 = exists(file);
                    file = resultFilename(outFileName+"1","h5");
                    //bool f2 = std::__fs::filesystem::exists(file);
                    bool f2 = exists(file);
                    file = resultFilename(outFileName+"2","h5");
                    //bool f3 = std::__fs::filesystem::exists(file);
                    bool f3 = exists(file);
                    if(f1 and f2 and f3){
                        std::cout << "Skipping src point " << x << y << z << t << std::endl;
                        continue;
                    }
                    
                    for (unsigned int c = 0; c < FImpl1::Dimension; ++c){
                        source = Zero();
                        Csrc=Zero();
                        Csrc()()(c)=Complex(1.0,0.0);
                        pokeSite(Csrc,source,srcSite);
                        sol = Zero();
                        // sol on odd sites first
                        pickCheckerboard(Even,halfSource,source);
                        action.Meooe(halfSource,tmp);
                        LMA(tmp, halfSol);
                        setCheckerboard(sol,halfSol);
                        // now even
                        LMA2(tmp, halfSol);
                        action.Meooe(halfSol,tmp);
                        setCheckerboard(sol,tmp);
                        sol *= -1.0;
                        FermToProp<FImpl1>(q1, sol, c);
                    }
                    
                    for(int mu=0;mu<3;mu++){
                        
                        srcSite[0]=x;
                        srcSite[1]=y;
                        srcSite[2]=z;
                        
                        LOG(Message) << "StagMesonCCLoopLMA src_xyzt " << srcSite[0] <<" "<< srcSite[1]<<" "<<srcSite[2] <<" "<< srcSite[3] <<" mu "<< mu << std::endl;
                        
                        peekSite(UmuSrc, Umu[mu], srcSite);
                        
                        srcSite[mu]=(srcSite[mu]+1)%ns;
                        
                        for (unsigned int c = 0; c < FImpl1::Dimension; ++c){
                            source = Zero();
                            Csrc=Zero();
                            Csrc()()(c)=Complex(1.0,0.0);
                            pokeSite(Csrc,source,srcSite);
                            sol = Zero();
                            // sol on odd sites first
                            // assuming the source is Odd from point-splitting
                            pickCheckerboard(Odd,halfSource,source);
                            LMA(halfSource, halfSol);
                            tmp = mass * halfSol;
                            setCheckerboard(sol,tmp);
                            // now even sol
                            action.Meooe(halfSol,tmp);
                            tmp *= -1.0;
                            setCheckerboard(sol,tmp);
                            FermToProp<FImpl1>(q2, sol, c);
                        }
                        
                        qshift = Cshift(q2, mu, 1);
                        corr = trace(adj(qshift) * adj(Umu[mu]) * q1 * UmuSrc);
                        corr += trace(adj(q1) * Umu[mu] * qshift * adj(UmuSrc));

                        qshift = Cshift(q1, mu, 1);
                        // minus signs from herm transform of prop
                        corr -= trace(adj(q2) * Umu[mu] * qshift * UmuSrc); // -1^muhat
                        corr -= trace(adj(qshift) * adj(Umu[mu]) * q2 * adj(UmuSrc)); //-1^muhat

                        corr *= herm_phase; // sink
                        if((x+y+z+t)%2)corr *= -1.0; // source
                        
                        sliceSum(corr, buf, Tp);
                        for (unsigned int tsnk = 0; tsnk < buf.size(); ++tsnk){
                            result.corr[tsnk] = TensorRemove(buf[tsnk]);
                        }
                        outFileName = par().output+"/cc_2pt_"+
                            std::to_string(x)+"_"+
                            std::to_string(y)+"_"+
                            std::to_string(z)+"_"+
                            std::to_string(t)+"_mu_"+
                            std::to_string(mu);
                        saveResult(outFileName, "mesonCC", result);
                        
                        // do the local current
                        corr = trace(adj(q1) * q1);
                        // ks phases (includes hermiticity phase)
                        corr *= localphases[mu]; // gamma at sink
                        if(      mu==0 && x % 2 )corr *= -1.0; // gamma phase at source
                        else if( mu==1 && y % 2 )corr *= -1.0; // gamma phase at source
                        else if( mu==2 && z % 2 )corr *= -1.0; // gamma phase at source
                        
                        sliceSum(corr, buf, Tp);
                        for (unsigned int tsnk = 0; tsnk < buf.size(); ++tsnk){
                            result.corr[tsnk] = TensorRemove(buf[tsnk]);
                        }
                        outFileName = par().output+"/local_2pt_"+
                            std::to_string(x)+"_"+
                            std::to_string(y)+"_"+
                            std::to_string(z)+"_"+
                            std::to_string(t)+"_mu_"+
                            std::to_string(mu);
                        saveResult(outFileName, "mesonLL", result);
                    }
                }
            }
        }
    }
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_MesonCCLoopLMA_hpp_
