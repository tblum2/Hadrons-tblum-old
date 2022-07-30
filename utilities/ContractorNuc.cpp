/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Utilities/Contractor.cc

Copyright (C) 2015-2018


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
#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/A2AMatrixNucleon.hpp>
#include <Hadrons/DiskVector.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Module.hpp>

using namespace Grid;
//using namespace QCD;
using namespace Hadrons;

#define TIME_MOD(t) (((t) + par.global.nt) % par.global.nt)

namespace Contractor
{
    class TrajRange: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                        unsigned int, start,
                                        unsigned int, end,
                                        unsigned int, step);
    };
    
    class GlobalPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar,
                                        TrajRange, trajCounter,
                                        unsigned int, nt,
                                        std::string, diskVectorDir,
                                        std::string, output);
    };

	// MCA - renamed this for 3-quark fields (i.e. nucleon A2A matrix)
    class A2AMatrixNucPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMatrixNucPar,
                                        std::string, file,
                                        std::string, dataset,
                                        unsigned int, cacheSize,
                                        std::string, name);
    };

// MCA - need to adjust this for nucleon 2pt function

    class ProductPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
                                        std::string, terms,
                                        std::vector<std::string>, times,
                                        std::string, translations,
                                        bool, translationAverage,
                                        int, boundaryT,
                                        std::string, output);
    };


    class CorrelatorResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult,
                                        std::vector<Contractor::A2AMatrixNucPar>,  a2aMatrixNuc,
                                        ProductPar, contraction,
                                        std::vector<unsigned int>, times,
                                        //std::vector<ComplexD>, correlator,
                                        std::vector<SpinMatrix>, corr);
    };
}

struct ContractorPar
{
    Contractor::GlobalPar                       global;
    std::vector<Contractor::A2AMatrixNucPar>    a2aMatrixNuc;
    std::vector<Contractor::ProductPar>    	    product;
};

// MCA - need to go through this
void makeTimeSeq(std::vector<std::vector<unsigned int>> &timeSeq, 
                 const std::vector<std::set<unsigned int>> &times,
                 std::vector<unsigned int> &current,
                 const unsigned int depth)
{
    if (depth > 0)
    {
        for (auto t: times[times.size() - depth])
        {
            current[times.size() - depth] = t;
            makeTimeSeq(timeSeq, times, current, depth - 1);
        }
    }
    else
    {
        timeSeq.push_back(current);
    }
}

// MCA - need to go through this
void makeTimeSeq(std::vector<std::vector<unsigned int>> &timeSeq, 
                 const std::vector<std::set<unsigned int>> &times)
{
    std::vector<unsigned int> current(times.size());

    makeTimeSeq(timeSeq, times, current, times.size());
}

// MCA - need to go through this

void saveCorrelator(const Contractor::CorrelatorResult &result, const std::string dir, 
                    const unsigned int dt, const unsigned int traj)
{
    std::string              fileStem = "", filename;
    std::vector<std::string> terms = strToVec<std::string>(result.contraction.terms);

	fileStem += "nuc2pt";

	// MCA - don't need this for 2pt
    /*for (unsigned int i = 0; i < terms.size() - 1; i++)
    {
        fileStem += terms[i] + "_" + std::to_string(result.times[i]) + "_";
    }
    fileStem += terms.back();*/
    if (!result.contraction.translationAverage)
    {
        fileStem += "_dt_" + std::to_string(dt);
    }
    filename = dir + "/" + ModuleBase::resultFilename(fileStem, traj);
    std::cout << "Saving correlator to '" << filename << "'" << std::endl;
    //std::cout << "Result correlator is " << result.corr << std::endl;
    makeFileDir(dir);
    ResultWriter writer(filename);
    write(writer, fileStem, result);
}


// MCA - need to go through this

std::set<unsigned int> parseTimeRange(const std::string str, const unsigned int nt)
{
    std::regex               rex("([0-9]+)|(([0-9]+)\\.\\.([0-9]+))");
    std::smatch              sm;
    std::vector<std::string> rstr = strToVec<std::string>(str);
    std::set<unsigned int>   tSet;

    for (auto &s: rstr)
    {
        std::regex_match(s, sm, rex);
        if (sm[1].matched)
        {
            unsigned int t;
            
            t = std::stoi(sm[1].str());
            if (t >= nt)
            {
                HADRONS_ERROR(Range, "time out of range (from expression '" + str + "')");
            }
            tSet.insert(t);
        }
        else if (sm[2].matched)
        {
            unsigned int ta, tb;

            ta = std::stoi(sm[3].str());
            tb = std::stoi(sm[4].str());
            if ((ta >= nt) or (tb >= nt))
            {
                HADRONS_ERROR(Range, "time out of range (from expression '" + str + "')");
            }
            for (unsigned int ti = ta; ti <= tb; ++ti)
            {
                tSet.insert(ti);
            }
        }
    }

    return tSet;
}

template <typename C>
void spinMatInit(C &mat)
{
    for (int mu = 0; mu < Ns; mu++)
    for (int nu = 0; nu < Ns; nu++)
    {
        mat(mu, nu) = 0.;
    }
}

struct Sec
{
    Sec(const double usec)
    {
        seconds = usec/1.0e6;
    }
    
    double seconds;
};

inline std::ostream & operator<< (std::ostream& s, const Sec &&sec)
{
    s << std::setw(10) << sec.seconds << " sec";

    return s;
}

struct Flops
{
    Flops(const double flops, const double fusec)
    {
        gFlopsPerSec = flops/fusec/1.0e3;
    }
    
    double gFlopsPerSec;
};

inline std::ostream & operator<< (std::ostream& s, const Flops &&f)
{
    s << std::setw(10) << f.gFlopsPerSec << " GFlop/s";

    return s;
}

struct Bytes
{
    Bytes(const double bytes, const double busec)
    {
        gBytesPerSec = bytes/busec*1.0e6/1024/1024/1024;
    }
    
    double gBytesPerSec;
};

inline std::ostream & operator<< (std::ostream& s, const Bytes &&b)
{
    s << std::setw(10) << b.gBytesPerSec << " GB/s";

    return s;
}

int main(int argc, char* argv[])
{
    
    // parse command line
    std::string   parFilename;

    if (argc != 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file>";
        std::cerr << std::endl;
        
        return EXIT_FAILURE;
    }
    parFilename = argv[1];

    // parse parameter file
    ContractorPar par;
    TimerArray tAr;
    unsigned int  nMat, nCont;
    XmlReader     reader(parFilename);

    read(reader, "global",    par.global);
    read(reader, "a2aMatrixNuc", par.a2aMatrixNuc);
    read(reader, "product",   par.product);
    nMat  = par.a2aMatrixNuc.size();
    nCont = par.product.size();

    // create diskvectors
    std::map<std::string, EigenDiskVectorNuc<HADRONS_A2AN_CALC_TYPE>> a2aMatNuc;
    unsigned int                                        cacheSize;

    for (auto &p: par.a2aMatrixNuc)
    {
        std::string dirName = par.global.diskVectorDir + "/" + p.name;

        a2aMatNuc.emplace(p.name, EigenDiskVectorNuc<HADRONS_A2AN_CALC_TYPE>(dirName, par.global.nt, p.cacheSize));
    }

    // trajectory loop
    for (unsigned int traj = par.global.trajCounter.start; 
         traj < par.global.trajCounter.end; traj += par.global.trajCounter.step)
    {
        std::cout << ":::::::: Trajectory " << traj << std::endl;

		std::cout << "Start loading" << std::endl;
        // load data
        tAr.startTimer("Total");
        tAr.startTimer("Load data from HDF5");
        for (auto &p: par.a2aMatrixNuc)
        {
            std::string filename = p.file;
            double      t, size;

			std::cout << "p.name " << p.name << std::endl;

            tokenReplace(filename, "traj", traj);
            std::cout << "======== Loading '" << filename << "'" << std::endl;

            A2AMatrixNucIo<HADRONS_A2AN_IO_TYPE> a2aNucIo(filename, p.dataset, par.global.nt);

			std::cout << "Still working here before loading" << std::endl;
            a2aNucIo.load(a2aMatNuc.at(p.name), &t);
            std::cout << "Read " << a2aNucIo.getSize() << " bytes in " << t/1.0e6 
                    << " sec, " << a2aNucIo.getSize()/t*1.0e6/1024/1024 << " MB/s" << std::endl;
        }
        tAr.stopTimer("Load data from HDF5");

		std::cout << "Testing a2aMatNuc data" << std::endl;
		
		for (int t = 0; t < par.global.nt; t++)
		{
			const A2AMatrixNuc<HADRONS_A2AN_CALC_TYPE> &ref = a2aMatNuc.at(par.a2aMatrixNuc[0].name)[0];
			std::cout << "t = " << t << " -- Dimensions";
			for (int i = 0; i < 4; i++)
			{
				std::cout << " " << ref.dimension(i);
			}
			std::cout << std::endl;
		}
        
        // contract
        EigenDiskVectorNuc<HADRONS_A2AN_CALC_TYPE>::Tensor buf;

		// Doing contraction testing
		/*
		TimerArray testAr;
		testAr.startTimer("Nucleon Test");
		A2AContractionNucleon::contNucTenTest();
		testAr.stopTimer("Nucleon Test");
		std::cout << Sec(testAr.getDTimer("Nucleon Test")) << std::endl;
		*/

		//A2AContractionNucleon::testProjTPlus();

		// MCA - remove this loop or adjust to do projectors for 2pt nucleon
		
        for (auto &p: par.product)
        {
			std::cout << "In the ProductPar loop" << std::endl;
			
            std::vector<std::string>               term = strToVec<std::string>(p.terms);
            //std::vector<std::set<unsigned int>>    times;
            //std::vector<std::vector<unsigned int>> timeSeq;
            std::set<unsigned int>                 translations;
            std::vector<A2AMatrixNuc<HADRONS_A2AN_CALC_TYPE>>    lastTerm(par.global.nt);
            A2AMatrixNuc<HADRONS_A2AN_CALC_TYPE>                 tenW;
            //A2AMatrixNuc<ComplexD>                 prod;
            //TimerArray                             tAr;
            double                                 fusec, busec, flops, bytes, tusec;
            Contractor::CorrelatorResult           result;
            std::vector<HADRONS_A2AN_CALC_TYPE>                  tmp_corr;
            std::vector<A2AMatrix<HADRONS_A2AN_CALC_TYPE>>       tmp_spinMat(par.global.nt);
            //std::vector<SpinMatrix>                tmp_spinMat(par.global.nt);


			//std::cout << "Dummy line" << std::endl;
			std::cout << "Printing out terms" << std::endl;
			for (int i = 0; i < term.size(); i++)
			{
				std::cout << term[i] << std::endl;
			}
			//std::cout << "Term.front is " << term.front() << std::endl;

            //tAr.startTimer("Total calc");
            // MCA - not needed for 2pt
			/*
            std::cout << "======== Contraction tr(";
            for (unsigned int g = 0; g < term.size(); ++g)
            {
                std::cout << term[g] << ((g == term.size() - 1) ? ')' : '*');
            }
            std::cout << std::endl;
            
            // error checking
            if (term.size() != p.times.size() + 1)
            {
                HADRONS_ERROR(Size, "number of terms (" + std::to_string(term.size()) 
                            + ") different from number of times (" 
                            + std::to_string(p.times.size() + 1) + ")");
            }
            
            // loop over insertion(?) times
            for (auto &s: p.times)
            {
				// add a set of times
                times.push_back(parseTimeRange(s, par.global.nt));
            }
            
            */
            
            for (auto &m: par.a2aMatrixNuc)
            {
                if (std::find(result.a2aMatrixNuc.begin(), result.a2aMatrixNuc.end(), m) == result.a2aMatrixNuc.end())
                {
                    result.a2aMatrixNuc.push_back(m);
                    tokenReplace(result.a2aMatrixNuc.back().file, "traj", traj);
                }
            }
            result.contraction = p;
            //result.correlator.resize(par.global.nt, 0.);
            
            result.corr.resize(par.global.nt);
            for (int t = 0; t < par.global.nt; t++)
            {
                result.corr[t] = Zero();
                //result.corrSpinMat[t].resize(Ns, Ns);
                //spinMatInit(result.corrSpinMat[t]);
            }
            
            tmp_corr.resize(par.global.nt, 0.);
			
			std::vector<unsigned int> debug_times = {0};
			result.times = debug_times; // DEBUG -- to keep saveCorrelator from crashing

            translations = parseTimeRange(p.translations, par.global.nt);

			/*
            // MCA - not needed for 2pt
            makeTimeSeq(timeSeq, times);
            std::cout << timeSeq.size()*translations.size()*(term.size() - 2) << " A*B, "
                    << timeSeq.size()*translations.size()*par.global.nt << " tr(A*B)"
                    << std::endl;
            */
            
			// MCA - may need to remove this due to memory limitations
                    
            std::cout << "* Caching last term sink times" << std::endl;
            for (unsigned int t = 0; t < par.global.nt; ++t)
            {
                tAr.startTimer("Disk vector overhead");
                const A2AMatrixNuc<HADRONS_A2AN_CALC_TYPE> &ref = a2aMatNuc.at(term.front())[t];
                tAr.stopTimer("Disk vector overhead");

                tAr.startTimer("Last term caching");
                lastTerm[t].resize(ref.dimension(0), ref.dimension(1), ref.dimension(2), ref.dimension(3));
                thread_for(mu,ref.dimension(0),{
                    thread_for(k,ref.dimension(3),{
                        thread_for(j,ref.dimension(2),{
                            for (unsigned int i = 0; i < ref.dimension(1); ++i)
                            {
                                lastTerm[t](mu, i, j, k) = ref(mu, i, j, k);
                            }
                        });
                    });
                });
                tAr.stopTimer("Last term caching");
            }
            bytes = par.global.nt*lastTerm[0].dimension(0)*lastTerm[0].dimension(1)*lastTerm[0].dimension(2)
									*lastTerm[0].dimension(3)*sizeof(HADRONS_A2AN_IO_TYPE);
            std::cout << Sec(tAr.getDTimer("Last term caching")) << " " 
                      << Bytes(bytes, tAr.getDTimer("Last term caching")) << std::endl;
             
             
            // MCA - remove this outer loop for nucleon 2pt?
            // it loops over timeSeq which is vector of sets of time insertions(?) defined in the XML file
            /*
            for (unsigned int i = 0; i < timeSeq.size(); ++i)
            {
                unsigned int dti = 0;
                auto         &t = timeSeq[i];

                result.times = t;
                for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                {
                    result.correlator[tLast] = 0.;
                }
             */

				// Initialize the correlator
                for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                {
                    //result.correlator[tLast] = 0.;
                    
                    tmp_spinMat[tLast].resize(Ns,Ns);
                    spinMatInit(tmp_spinMat[tLast]);
                    //tmp_spinMat[tLast] = Zero();
                                        
                }

             
                for (auto &dt: translations)
                {
                    std::cout << "* Step " << dt + 1
                            << "/" << translations.size()
                            << " -- dt= " << dt << std::endl;
                    // MCA - don't need for 2pt
                    /*
                    if (term.size() > 2)
                    {
                        std::cout << std::setw(8) << "products";
                    }
                    */
                    flops  = 0.;
                    bytes  = 0.;
                    fusec  = tAr.getDTimer("A*B algebra"); // MCA - want to change this
                    busec  = tAr.getDTimer("A*B total"); // MCA - this also
                    tAr.startTimer("Linear algebra");
                    tAr.startTimer("Disk vector overhead");
                    const A2AMatrixNuc<HADRONS_A2AN_CALC_TYPE> &ref = a2aMatNuc.at(term.back())[dt];
                    tenW.resize(ref.dimension(0), ref.dimension(1), ref.dimension(2), ref.dimension(3));
 		    thread_for(mu,ref.dimension(0),{
                        thread_for(k,ref.dimension(3),{
                            thread_for(j,ref.dimension(2),{
                                for (unsigned int i = 0; i < ref.dimension(1); ++i)
                                {
                                    tenW(mu, i, j, k) = ref(mu, i, j, k);
                                }
                            });
                        });
                    });
                    tAr.stopTimer("Disk vector overhead");
                    
                    flops  = 0.;
                    bytes  = 0.;
                    fusec  = tAr.getDTimer("tr(A*B)"); // this too
                    busec  = tAr.getDTimer("tr(A*B)"); // same 
                    // MCA - core computation loop -- needs to be cleaned up
                    // that contracts two nucleon LMA/A2A fields leaving a 4x4 spin matrix that will be
                    // projected upon and traced over
                    //for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast){
                    thread_for(tLast, par.global.nt, {
			tmp_corr[TIME_MOD(tLast - dt)] = 0.;
                        spinMatInit(tmp_spinMat[TIME_MOD(tLast - dt)]);
			//tmp_spinMat[TIME_MOD(tLast - dt)] = Zero();
                        //tAr.startTimer("tr(A*B)"); // adjust this
                        
                        // do contractions
                        A2AContractionNucleon::contNucTen(tmp_spinMat[TIME_MOD(tLast - dt)], lastTerm[tLast], tenW);
                        if (tLast < dt)
                        {
			   tmp_spinMat[TIME_MOD(tLast - dt)] *= p.boundaryT;
			}

                        for (int mu = 0; mu < Ns; mu++)
                        for (int nu = 0; nu < Ns; nu++)
                        {
                            result.corr[TIME_MOD(tLast - dt)]()(mu, nu)() += tmp_spinMat[TIME_MOD(tLast - dt)](mu, nu);
                        }
                        //tAr.stopTimer("tr(A*B)");
                    });
                    tAr.stopTimer("Linear algebra");
                    //std::cout << Sec(tAr.getDTimer("tr(A*B)") - busec) << " "
                     //       << Flops(flops, tAr.getDTimer("tr(A*B)") - fusec) << " " 
                      //      << Bytes(bytes, tAr.getDTimer("tr(A*B)") - busec) << std::endl;
                            
                    // if we're not averaging over source times -- save data to file and flush correlator data
                    if (!p.translationAverage)
                    {
                        saveCorrelator(result, p.output, dt, traj);
                        for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                        {
                            //result.correlator[tLast] = 0.;
                            result.corr[tLast] = Zero();
                        }
                    }
                    //dti++;
                }
                
                // if we're averaging over all source times
                if (p.translationAverage)
                {
                    for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                    {
                        //result.correlator[tLast] /= translations.size();
                        for (int mu = 0; mu < Ns; mu++)
                        for (int nu = 0; nu < Ns; nu++)
                        {
                            result.corr[tLast]()(mu,nu)() *= (1. / translations.size());
                        }
                        std::cout << tLast << " -- " << result.corr[tLast] << std::endl;
                    }
                    saveCorrelator(result, p.output, 0, traj);
                    
                    //std::cout << "Boundary condition in T direction is " << p.boundaryT << std::endl;
                    /*
                    for (unsigned int tLast = 0; tLast < par.global.nt; tLast++)
                    {
						std::cout << tLast << " - " << result.corr[tLast] << std::endl;
					}
                    */
                }
            //}  end of times loop (don't need for 2pt)
            //tAr.stopTimer("Total calc");
            
        }
    
    tAr.stopTimer("Total");
    printTimeProfile(tAr.getTimings(), tAr.getTimer("Total"));
    }
    
    return EXIT_SUCCESS;
}
