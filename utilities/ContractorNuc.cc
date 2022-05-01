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

using namespace Grid;
using namespace QCD;
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
                                        std::string, projectors,
                                        int, boundaryT);
    };


    class CorrelatorResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult,
                                        std::vector<Contractor::A2AMatrixNucPar>,  a2aMatrixNuc,
                                        ProductPar, contraction,
                                        std::vector<unsigned int>, times,
                                        std::vector<ComplexD>, correlator);
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
    filename = dir + "/" + RESULT_FILE_NAME(fileStem, traj);
    std::cout << "Saving correlator to '" << filename << "'" << std::endl;
    std::cout << "Result correlator is " << result.correlator << std::endl;
    makeFileDir(dir);
    ResultWriter writer(filename);
    std::cout << "Still working here" << std::endl;
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
    unsigned int  nMat, nCont;
    XmlReader     reader(parFilename);

    read(reader, "global",    par.global);
    read(reader, "a2aMatrixNuc", par.a2aMatrixNuc);
    read(reader, "product",   par.product);
    nMat  = par.a2aMatrixNuc.size();
    nCont = par.product.size();

    // create diskvectors
    std::map<std::string, EigenDiskVectorNuc<ComplexD>> a2aMatNuc;
    unsigned int                                        cacheSize;

    for (auto &p: par.a2aMatrixNuc)
    {
        std::string dirName = par.global.diskVectorDir + "/" + p.name;

        a2aMatNuc.emplace(p.name, EigenDiskVectorNuc<ComplexD>(dirName, par.global.nt, p.cacheSize));
    }

    // trajectory loop
    for (unsigned int traj = par.global.trajCounter.start; 
         traj < par.global.trajCounter.end; traj += par.global.trajCounter.step)
    {
        std::cout << ":::::::: Trajectory " << traj << std::endl;

		std::cout << "Testing nucleon loading" << std::endl;
        // load data
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

		std::cout << "Testing a2aMatNuc data" << std::endl;
		
		for (int t = 0; t < par.global.nt; t++)
		{
			const A2AMatrixNuc<ComplexD> &ref = a2aMatNuc.at("nuc2pt_vvv")[0];
			std::cout << "t = " << t << " -- Dimensions";
			for (int i = 0; i < 4; i++)
			{
				std::cout << " " << ref.dimension(i);
			}
			std::cout << std::endl;
		}
        // contract
        EigenDiskVectorNuc<ComplexD>::Tensor buf;

		// Doing contraction testing
		/*
		TimerArray testAr;
		testAr.startTimer("Nucleon Test");
		A2AContractionNucleon::contNucTenTest();
		testAr.stopTimer("Nucleon Test");
		std::cout << Sec(testAr.getDTimer("Nucleon Test")) << std::endl;
		*/

		std::cout << "Just adding to recompile" << std::endl;
		A2AContractionNucleon::testProjTPlus();

		// MCA - remove this loop or adjust to do projectors for 2pt nucleon
		
        for (auto &p: par.product)
        {
			std::cout << "In the ProductPar loop" << std::endl;
			std::cout << p.projectors << std::endl;
			
            std::vector<std::string>               term = strToVec<std::string>(p.terms);
            //std::vector<std::set<unsigned int>>    times;
            //std::vector<std::vector<unsigned int>> timeSeq;
            std::set<unsigned int>                 translations;
            std::vector<A2AMatrixNuc<ComplexD>>    lastTerm(par.global.nt);
            A2AMatrixNuc<ComplexD>                 tenW;
            //A2AMatrixNuc<ComplexD>                 prod;
            TimerArray                             tAr;
            double                                 fusec, busec, flops, bytes, tusec;
            Contractor::CorrelatorResult           result;
            std::vector<ComplexD>                  tmp_corr;


			//std::cout << "Dummy line" << std::endl;
			std::cout << "Printing out terms" << std::endl;
			for (int i = 0; i < term.size(); i++)
			{
				std::cout << term[i] << std::endl;
			}
			std::cout << "Term.front is " << term.front() << std::endl;

            tAr.startTimer("Total");
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
            result.correlator.resize(par.global.nt, 0.);
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
                const A2AMatrixNuc<ComplexD> &ref = a2aMatNuc.at(term.front())[t];
                tAr.stopTimer("Disk vector overhead");

                tAr.startTimer("Last term caching");
                lastTerm[t].resize(ref.dimension(0), ref.dimension(1), ref.dimension(2), ref.dimension(3));
                parallel_for (unsigned int mu = 0; mu < ref.dimension(0) ; mu++)
                parallel_for (unsigned int k = 0; k < ref.dimension(3); k++)
                parallel_for (unsigned int j = 0; j < ref.dimension(2); ++j)
                for (unsigned int i = 0; i < ref.dimension(1); ++i)
                {
                    lastTerm[t](mu, i, j, k) = ref(mu, i, j, k);
                }
                tAr.stopTimer("Last term caching");
            }
            bytes = par.global.nt*lastTerm[0].dimension(0)*lastTerm[0].dimension(1)*lastTerm[0].dimension(2)
									*lastTerm[0].dimension(3)*sizeof(ComplexD);
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
                    result.correlator[tLast] = 0.;
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
                    const A2AMatrixNuc<ComplexD> &ref = a2aMatNuc.at(term.back())[dt];
                    tenW.resize(ref.dimension(0), ref.dimension(1), ref.dimension(2), ref.dimension(3));
					parallel_for (unsigned int mu = 0; mu < ref.dimension(0) ; mu++)
					parallel_for (unsigned int k = 0; k < ref.dimension(3); k++)
					parallel_for (unsigned int j = 0; j < ref.dimension(2); ++j)
					for (unsigned int i = 0; i < ref.dimension(1); ++i)
					{
						tenW(mu, i, j, k) = ref(mu, i, j, k);
					}
                    tAr.stopTimer("Disk vector overhead");
                    // MCA - don't think I need this, need to double check though
                    /*
                    for (unsigned int j = 1; j < term.size() - 1; ++j)
                    {
                        tAr.startTimer("Disk vector overhead");
                        //const A2AMatrix<ComplexD> &ref = a2aMat.at(term[j])[TIME_MOD(t[j] + dt)];
                        tAr.stopTimer("Disk vector overhead");
                        
                        tAr.startTimer("A*B total");
                        tAr.startTimer("A*B algebra");
                        //A2AContraction::mul(tmp, prod, ref); -- replace this with nuc2pt code
                        tAr.stopTimer("A*B algebra");
                        //flops += A2AContraction::mulFlops(prod, ref); - replace here too
                        //prod   = tmp;
                        tAr.stopTimer("A*B total");
                        //bytes += 3.*tmp.rows()*tmp.cols()*sizeof(ComplexD);
                    }
                    */
                    // MCA - don't need for 2pt
                    /*
                    if (term.size() > 2)
                    {
                        std::cout << Sec(tAr.getDTimer("A*B total") - busec) << " "
                                << Flops(flops, tAr.getDTimer("A*B algebra") - fusec) << " " 
                                << Bytes(bytes, tAr.getDTimer("A*B total") - busec) << std::endl;
                    }
                    */
                    //std::cout << std::setw(8) << "traces"; - may change this
                    flops  = 0.;
                    bytes  = 0.;
                    fusec  = tAr.getDTimer("tr(A*B)"); // this too
                    busec  = tAr.getDTimer("tr(A*B)"); // same 
                    // MCA - core computation loop -- replace this with nucleon code in A2AMatrixNuc
                    // that contracts two nucleon LMA/A2A fields leaving a 4x4 spin matrix that will be
                    // projected upon and traced over
                    for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                    {
						tmp_corr[TIME_MOD(tLast - dt)] = 0.;
                        tAr.startTimer("tr(A*B)"); // adjust this
                        //A2AContraction::accTrMul(result.correlator[TIME_MOD(tLast - dt)], prod, lastTerm[tLast]);
                        A2AContractionNucleon::ContractNucleonTPlus(tmp_corr[TIME_MOD(tLast - dt)], lastTerm[tLast], tenW);
                        //A2AContractionNucleon::ContractNucleonTPlus(tmp_corr[tLast], lastTerm[tLast], tenV);
                        // MCA - antiperiodic boundary conditions -- need to double check this
                        
                        // if antiperiodic (boundaryT = -1) do sign change
                        if (tLast < dt)
                        {
							tmp_corr[TIME_MOD(tLast - dt)] *= p.boundaryT;
						}
						
						std::cout << "tLast " << tLast << " - dt " << dt << " | tsep " << TIME_MOD(tLast - dt) << 
								" == " << tmp_corr[TIME_MOD(tLast - dt)] << std::endl;
						
						//std::cout << "Summing correlator for tsnk - tsrc " << TIME_MOD(tLast - dt) << std::endl;
						//std::cout << "Existing: " << result.correlator[TIME_MOD(tLast - dt)] << " --- "
						//			<< "Adding: " << tmp_corr[TIME_MOD(tLast - dt)] << std::endl;
						/*if (TIME_MOD(tLast - dt) == 6)
						{
							// dummy line
							std::cout << "tLast - tSrc = " << TIME_MOD(tLast - dt) << " : " << tmp_corr[TIME_MOD(tLast - dt)]
										<< std::endl;
						}*/
						result.correlator[TIME_MOD(tLast - dt)] += tmp_corr[TIME_MOD(tLast - dt)];
						//result.correlator[tLast] += tmp_corr[tLast];
                        tAr.stopTimer("tr(A*B)");
                        //flops += A2AContraction::accTrMulFlops(prod, lastTerm[tLast]);
                        //bytes += 2.*prod.rows()*prod.cols()*sizeof(ComplexD);
                    }
                    tAr.stopTimer("Linear algebra");
                    std::cout << Sec(tAr.getDTimer("tr(A*B)") - busec) << " "
                            << Flops(flops, tAr.getDTimer("tr(A*B)") - fusec) << " " 
                            << Bytes(bytes, tAr.getDTimer("tr(A*B)") - busec) << std::endl;
                    if (!p.translationAverage)
                    {
                        saveCorrelator(result, par.global.output, dt, traj);
                        for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                        {
                            result.correlator[tLast] = 0.;
                        }
                    }
                    //dti++;
                }
                if (p.translationAverage)
                {
                    for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                    {
                        result.correlator[tLast] /= translations.size();
                        std::cout << tLast << " -- " << result.correlator[tLast] << std::endl;
                    }
                    saveCorrelator(result, par.global.output, 0, traj);
                    
                    // MCA - DEBUG
                    //std::cout << "Printing correlator" << std::endl;
                    std::cout << "Boundary condition in T direction is " << p.boundaryT << std::endl;
                    //std::cout << "This is a new test" << std::endl;
                    for (unsigned int tLast = 0; tLast < par.global.nt; tLast++)
                    {
						std::cout << tLast << " - " << result.correlator[tLast] << std::endl;
					}
                }
            //}  end of times loop (don't need for 2pt)
            tAr.stopTimer("Total");
            printTimeProfile(tAr.getTimings(), tAr.getTimer("Total"));
            
        } 
        
    }
    
    return EXIT_SUCCESS;
}
