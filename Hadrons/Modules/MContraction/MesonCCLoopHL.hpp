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

#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MSource/Point.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/*

 2pt conserved-conserved staggered correlation function
 -----------------------------

 * options:
 - q1: input propagator 1 (string) src at y
 - q2: input propagator 2 (string) src at y+hat mu

*/

/******************************************************************************
 *                                TMesonLoopCCHL                                    *
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
                                    double, mass,
                                    int, inc,
                                    int, tinc);
};

template <typename FImpl1, typename FImpl2>
class TStagMesonLoopCCHL: public Module<MesonLoopCCHLPar>
{
public:
    typedef typename FImpl1::FermionField FermionField;
    typedef A2AVectorsSchurStaggered<FImpl1> A2A;
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    SOLVER_TYPE_ALIASES(FImpl1,);
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
    Solver       *solver_{nullptr};
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
    std::vector<std::string> input = {par().gauge, par().solver, par().eigenPack};

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
    envTmpLat(LatticeComplex, "corr");
    envTmpLat(PropagatorField1, "q1");
    envTmpLat(PropagatorField2, "q2");
    envTmpLat(PropagatorField1, "qshift");
    envTmpLat(FermionField, "source");
    envTmpLat(FermionField, "sol");
    envTmpLat(FermionField, "v");

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
    
    envTmp(A2A, "a2a", 1, action, solver);
    
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
    auto &solver  = envGet(Solver, par().solver);
    auto &mat     = solver.getFMat();
    auto &epack   = envGet(BaseFermionEigenPack<FImpl1>, par().eigenPack);
    double mass = par().mass;
    
    envGetTmp(LatticeComplex, corr);
    envGetTmp(LatticeComplex, herm_phase);
    envGetTmp(PropagatorField1, q1);
    envGetTmp(PropagatorField2, q2);
    envGetTmp(PropagatorField1, qshift);
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
    envGetTmp(FermionField, source);
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, v);
    Coordinate srcSite;
    ColourMatrix UmuSrc;
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
    for (unsigned int il = 0; il < Nl_; il++)
    {
        // eval of unpreconditioned Dirac op
        std::complex<double> eval(mass,sqrt(epack.eval[il]-mass*mass));
        
        LOG(Message) << "V vector i = " << 2*il << ", " << 2*il+1 << " (low modes)" << std::endl;
        
        for (unsigned int eo = 0; eo < 2; eo++)
        {
            // construct full lattice evec
            a2a.makeLowModeV(v, epack.evec[il], eval, eo);
    
            // loop over source time slices
            for(int ts=0; ts<nt;ts+=par().tinc){
                
                LOG(Message) << "StagMesonLoopCCHLHL src_t " << srcSite[3] << std::endl;


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

                //source = 1.;
                source = where(t == ts, v, 0.);
                sol = Zero();
                solver(sol, source);
                //FermToProp<FImpl1>(q1, sol, c);
            

//                for(int mu=0;mu<3;mu++){
//
//                    LOG(Message) << "StagMesonLoopCCHLHL src_mu " << mu << std::endl;
//
//                    peekSite(UmuSrc, Umu[mu], srcSite);
//
//                    srcSite[mu]=(srcSite[mu]+1)%ns;
//
//
//                    source = Zero();
//                    sol = Zero();
//                    solver(sol, source);
//                    //FermToProp<FImpl1>(q2, sol, c);
//
//
//                    qshift = Cshift(q2, mu, 1);
//                    //corr = trace(adj(qshift) * adj(Umu[mu]) * q1 * UmuSrc);
//                    corr = innerProduct(Umu[mu] * qshift, sol * UmuSrc);
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
//                    sliceSum(corr, buf, Tp);
//                    for (unsigned int tsnk = 0; tsnk < buf.size(); ++tsnk){
//                        result.corr[tsnk] = TensorRemove(buf[tsnk]);
//                    }
//                    outFileName = par().output+"/cc_2pt_"+std::to_string(t)+"_mu_"+
//                    std::to_string(mu);
//                    saveResult(outFileName, "mesonCC", result);
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
//                }
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
    }
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_MesonLoopCCHL_hpp_
