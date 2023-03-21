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

template <typename T>
using A2AMatrix3index = Eigen::Tensor<T,3,Eigen::RowMajor>;

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
                                        std::vector<SpinMatrix>, correlator);
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

void saveCorrelator(const Contractor::CorrelatorResult &result, const std::string dir,
                    const unsigned int dt, const unsigned int traj)
{
    std::string              fileStem = "", filename;
    std::vector<std::string> terms = strToVec<std::string>(result.contraction.terms);

    fileStem += "nuc2pt";

    if (!result.contraction.translationAverage)
    {
        fileStem += "_dt_" + std::to_string(dt);
    }
    filename = dir + "/" + ModuleBase::resultFilename(fileStem, traj);
    std::cout << "Saving correlator to '" << filename << "'" << std::endl;
    //DEBUG: std::cout << "Result correlator is " << result.correlator << std::endl;
    makeFileDir(dir);
    ResultWriter writer(filename);
    write(writer, fileStem, result);
}

// correlation function output file creation
void initCorrFile(const std::string filename,
                  const unsigned int ni)
{
#ifdef HAVE_HDF5

    try{
        H5::H5File file(filename, H5F_ACC_EXCL);

        // Create the data space for the dataset.
        hsize_t dims[1];               // dataset dimension: 1d for time dir
        dims[0] = ni;
        H5::DataSpace dataspace(1, dims);

        // Create the dataset.
        H5::DataSet dataset = file.createDataSet("corrfunc_dset", Hdf5Type<ComplexD>::type(), dataspace);
    }
    // catch failure caused by the H5File operations
    catch(H5::FileIException error)
    {
        error.printErrorStack();
    }

#else
    HADRONS_ERROR(Implementation, "correlation function I/O needs HDF5 library");
#endif
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

    if (argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file> [grid options]";
        std::cerr << std::endl;
        
        return EXIT_FAILURE;
    }
    parFilename = argv[1];

    // initialization
    Grid_init(&argc, &argv);
    
    GridCartesian *Grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(4,1), GridDefaultMpi());
    
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
    unsigned int cacheSize;
    
    // time size per process
    int localNt = Grid->LocalDimensions()[3];
    
    for (auto &p: par.a2aMatrixNuc)
    {
        std::string dirName = par.global.diskVectorDir + std::to_string(Grid->ThisRank()) + "/" + p.name;
        a2aMatNuc.emplace(p.name, EigenDiskVectorNuc<ComplexD>(dirName, localNt, p.cacheSize));
    }

    // trajectory loop
    for (unsigned int traj = par.global.trajCounter.start; 
         traj < par.global.trajCounter.end; traj += par.global.trajCounter.step)
    {
        std::cout << ":::::::: Trajectory " << traj << std::endl;

        // load data
        for (auto &p: par.a2aMatrixNuc)
        {
            std::string filename = p.file;
            double      t, size;

			std::cout << "p.name " << p.name << std::endl;

            tokenReplace(filename, "traj", traj);
            std::cout << "======== Loading '" << filename << "'" << std::endl;

            A2AMatrixNucIo<HADRONS_A2AN_IO_TYPE> a2aNucIo(filename, p.dataset, localNt);

            a2aNucIo.load(a2aMatNuc.at(p.name), Grid, &t);
            std::cout << "Read " << a2aNucIo.getSize() << " bytes in " << t/1.0e6 
                      << " sec, " << a2aNucIo.getSize()/t*1.0e6/1024/1024 << " MB/s" << std::endl;
        }
        
        // contract
        EigenDiskVectorNuc<ComplexD>::Tensor buf;

		A2AContractionNucleon::testProjTPlus();

        for (auto &p: par.product)
        {
            std::vector<std::string>               term = strToVec<std::string>(p.terms);
            std::set<unsigned int>                 translations;
            std::vector<A2AMatrixNuc<ComplexD>>    lastTerm(localNt);
            A2AMatrixNuc<ComplexD>                 tenW;
            TimerArray                             tAr;
            double                                 fusec, busec, flops, bytes, tusec;
            Contractor::CorrelatorResult           result;
 	    std::vector<A2AMatrix<ComplexD>>       tmp_spinMat(par.global.nt);


			std::cout << "Baryon Fields:" << std::endl;
			for (int i = 0; i < term.size(); i++)
			{
				std::cout << term[i] << std::endl;
			}
			std::cout << "Term.front is " << term.front() << std::endl;

            tAr.startTimer("Total");
            
            for (auto &m: par.a2aMatrixNuc)
            {
                if (std::find(result.a2aMatrixNuc.begin(), result.a2aMatrixNuc.end(), m) == result.a2aMatrixNuc.end())
                {
                    result.a2aMatrixNuc.push_back(m);
                    tokenReplace(result.a2aMatrixNuc.back().file, "traj", traj);
                }
            }
            result.contraction = p;
            result.correlator.resize(par.global.nt);
            for (int t = 0; t < par.global.nt; t++)
            {
                result.correlator[t] = Zero();
            }
	    std::vector<unsigned int> debug_times = {0};
	    result.times = debug_times; // DEBUG -- to keep saveCorrelator from crashing

            translations = parseTimeRange(p.translations, par.global.nt);

	    std::cout << "* Caching last term sink times" << std::endl;
            for (unsigned int t = 0; t < localNt; t++) {

                tAr.startTimer("Disk vector overhead");
                const A2AMatrixNuc<ComplexD> &ref = a2aMatNuc.at(term.front())[t];
                tAr.stopTimer("Disk vector overhead");

                tAr.startTimer("Last term caching");
                lastTerm[t].resize(ref.dimension(0), ref.dimension(1), ref.dimension(2), ref.dimension(3));
                thread_for(mu,ref.dimension(0),{
                    thread_for(k,ref.dimension(3),{
                        thread_for(j,ref.dimension(2),{
                            for (int i=0;i<ref.dimension(1);i++)
			    {
                                lastTerm[t]((long)mu, (long)i, (long)j, (long)k) = ref((long)mu, (long)i, (long)j, (long)k);
                            }
                        });
                    });
                });
                tAr.stopTimer("Last term caching");
            }
            bytes = localNt*lastTerm[0].dimension(0)*lastTerm[0].dimension(1)*lastTerm[0].dimension(2)
		   	   *lastTerm[0].dimension(3)*sizeof(ComplexD);
            std::cout << Sec(tAr.getDTimer("Last term caching")) << " " 
                      << Bytes(bytes, tAr.getDTimer("Last term caching")) << std::endl;

            // Initialize the correlator
	    for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
            {
                tmp_spinMat[tLast].resize(Ns,Ns);
                spinMatInit(tmp_spinMat[tLast]);
            }

            for (auto &dt: translations)
            {
                std::cout << "* Step " << dt + 1 << "/" << translations.size() << " -- dt= " << dt << std::endl;
                // MCA - don't need for 2pt
                    
                flops  = 0.;
                bytes  = 0.;
                fusec  = tAr.getDTimer("A*B algebra"); // MCA - want to change this
                busec  = tAr.getDTimer("A*B total"); // MCA - this also
                tAr.startTimer("Linear algebra");
                tAr.startTimer("Disk vector overhead");
                
                // Is src time slice on node?
                int srcNode = dt/localNt;
                std::cout<<" dt= "<<dt<<" src node= "<<srcNode<<std::endl;
                
		uint64_t a2abytes; 
                const A2AMatrixNuc<ComplexD> &ref = a2aMatNuc.at(term.back())[dt%localNt];
                tenW.resize(ref.dimension(0), ref.dimension(1), ref.dimension(2), ref.dimension(3));
		A2AMatrix3index<ComplexD> temp;
                temp.resize(ref.dimension(1), ref.dimension(2), ref.dimension(3));
                a2abytes = sizeof(ComplexD)*ref.dimension(1)*ref.dimension(2)*ref.dimension(3);
		// do 1 spin at a time so MPI_bcast doesn't break
                // thread_for broke something here
                for(int mu=0;mu<ref.dimension(0);mu++)
		{
		    if(srcNode==Grid->ThisRank())
		    {
                       for(int k=0;k<ref.dimension(3);k++){
                           for(int j=0;j<ref.dimension(2);j++){
                               for (unsigned int i = 0; i < ref.dimension(1); i++)
                               {
		         	   temp((long)i,(long)j,(long)k) = ref((long)mu, (long)i, (long)j ,(long)k);
                               }
                            }
                        }
                    }
                    Grid->Broadcast(srcNode, temp.data(), a2abytes);
                    for(int k=0;k<ref.dimension(3);k++){
                        for(int j=0;j<ref.dimension(2);j++){
                           for (unsigned int i = 0; i < ref.dimension(1); i++)
                           {
                               tenW((long)mu, (long)i, (long)j, (long)k) = temp((long)i,(long)j,(long)k);
                           }
                        }
                    }
                }
 
                tAr.stopTimer("Disk vector overhead");
                
                flops  = 0.;
                bytes  = 0.;
                fusec  = tAr.getDTimer("tr(A*B)"); // this too
                busec  = tAr.getDTimer("tr(A*B)"); // same
                for (unsigned int tLast = 0; tLast < localNt; ++tLast)
		{
                
                    int gt = TIME_MOD(tLast+Grid->ThisRank()*localNt - dt);
                    spinMatInit(tmp_spinMat[gt]);
                    tAr.startTimer("tr(A*B)"); // adjust this
		    A2AContractionNucleon::contNucTen(tmp_spinMat[gt], lastTerm[tLast], tenW);
                    // if antiperiodic (boundaryT = -1) do sign change for tsink < tsrc
                    if (tLast + Grid->ThisRank()*localNt < dt)
                    {
                       tmp_spinMat[gt] *= p.boundaryT;
                    }

                    for (int mu = 0; mu < Ns; mu++)
                       for (int nu = 0; nu < Ns; nu++)
                       {
                            result.correlator[gt]()(mu, nu)() += tmp_spinMat[gt](mu, nu);
                       }
				
                    tAr.stopTimer("tr(A*B)");
                }
                tAr.stopTimer("Linear algebra");
                std::cout << Sec(tAr.getDTimer("tr(A*B)") - busec) << " "
                        << Flops(flops, tAr.getDTimer("tr(A*B)") - fusec) << " "
                        << Bytes(bytes, tAr.getDTimer("tr(A*B)") - busec) << std::endl;
                if (!p.translationAverage)
                {
                    for (unsigned int t = 0; t < par.global.nt; ++t)
                    {
                        Grid->GlobalSum(result.correlator[t]);
                    } 
                    if(Grid->IsBoss()){
                        saveCorrelator(result, p.output, dt, traj);
                    }
                    for (unsigned int t = 0; t < par.global.nt; ++t)
                    {
                        result.correlator[t] = Zero();
                    }
                }
            }
            if (p.translationAverage)
            {
                for (unsigned int t = 0; t < par.global.nt; ++t)
                {
                    Grid->GlobalSum(result.correlator[t]);
                } 
                if(Grid->IsBoss()){
                    saveCorrelator(result, p.output, 0, traj);
                }
                /*std::cout << "Boundary condition in T direction is " << p.boundaryT << std::endl;
                for (unsigned int t = 0; t < par.global.nt; t++)
                {
		    for (int mu = 0; mu < Ns; mu++)
                        for (int nu = 0; nu < Ns; nu++)
                        {
                            result.correlator[t]()(mu,nu)() *= (1. / translations.size());
                        }
                    std::cout << t << " - " << result.correlator[t] << std::endl;
                }*/
            }
            tAr.stopTimer("Total");
            printTimeProfile(tAr.getTimings(), tAr.getTimer("Total"));
        }
    }
    return EXIT_SUCCESS;
}
