#ifndef Hadrons_MUtilities_A2ACovariantSmear_hpp_
#define Hadrons_MUtilities_A2ACovariantSmear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/CovariantSmearing.h>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2ACovariantSmear                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class A2ACovariantSmearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ACovariantSmearPar,
                                    std::string, a2aVectors,
                                    std::string, gauge,
                                    double, alpha,
                                    unsigned int, N,
                                    unsigned int, orthog_axis,
                                    bool, alt_convention,
                                    std::string, output,
                                    bool, multiFile);
};

template <typename GImpl, typename FImpl>
class TA2ACovariantSmear: public Module<A2ACovariantSmearPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2ACovariantSmear(const std::string name);
    // destructor
    virtual ~TA2ACovariantSmear(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2ACovariantSmear, ARG(TA2ACovariantSmear<GIMPL, FIMPL>), MUtilities);

/******************************************************************************
 *                 TA2ACovariantSmear implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
TA2ACovariantSmear<GImpl, FImpl>::TA2ACovariantSmear(const std::string name)
: Module<A2ACovariantSmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
std::vector<std::string> TA2ACovariantSmear<GImpl, FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().a2aVectors, par().gauge};
    
    return in;
}

template <typename GImpl, typename FImpl>
std::vector<std::string> TA2ACovariantSmear<GImpl, FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2ACovariantSmear<GImpl, FImpl>::setup(void)
{
	envCreateLat(FermionField, getName());    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2ACovariantSmear<GImpl, FImpl>::execute(void)
{
	auto &a2aVec = envGet(std::vector<FermionField>, par().a2aVectors);
	unsigned int Ni = a2aVec.size();

	const auto &U = envGet(GaugeField, par().gauge);

	std::vector<LatticeColourMatrix> Umu(Nd,env().getGrid());
	for(int i=0;i<Nd;i++)
	{
        Umu[i] = PeekIndex<LorentzIndex>(U,i);
	}
	
    auto &a2aTmp = envGet(FermionField, getName());
	
	RealD width = par().alpha;
	unsigned int iterations = par().N;
	unsigned int orthog_axis = par().orthog_axis;

    LOG(Message) << "The A2ACovariantSmear module is working." << std::endl;
    LOG(Message) << par().a2aVectors << " has size " << Ni << std::endl;
    
	CovariantSmearing<GImpl> CovSmear;
	
	for (int i = 0; i < Ni; i++)
	{
		a2aTmp = a2aVec[i];

		//LOG(Message) << "Testing a2aVector output before smearing" << std::endl;
		//LOG(Message) << a2aTmp << std::endl;

		if (par().alt_convention == true)
		{
			CovSmear.template GaussianSmearAlt<FermionField>(Umu, a2aTmp, width, iterations, orthog_axis);
		}
		else
		{
			CovSmear.template GaussianSmear<FermionField>(Umu, a2aTmp, width, iterations, orthog_axis);
		}
		
		//LOG(Message) << "Testing a2aVector output after smearing" << std::endl;
		//LOG(Message) << a2aTmp << std::endl;
		
		a2aVec[i] = a2aTmp;
	}

    // I/O if necessary
    if (!par().output.empty())
    {
        startTimer("I/O");
        A2AVectorsIo::write(par().output, a2aVec, par().multiFile, vm().getTrajectory());
        stopTimer("I/O");
    }
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_A2ACovariantSmear_hpp_
