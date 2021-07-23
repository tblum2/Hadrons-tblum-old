/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Utilities/Contractor.cc

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
        modified for cons. current, T. Blum

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
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DiskVector.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/TimerArray.hpp>

using namespace Grid;
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

    class A2AMatrixPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMatrixPar,
                                        std::string, file,
                                        std::string, dataset,
                                        unsigned int, cacheSize,
                                        std::string, name);
    };

    class EvalPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(EvalPar,
                                        std::string, file,
                                        unsigned int, size);
    };
    
    class ProductPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
                                        std::string, terms,
                                        std::vector<std::string>, times,
                                        std::string, translations,
                                        bool, translationAverage);
    };

    class CorrelatorResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult,
                                        std::vector<Contractor::A2AMatrixPar>,  a2aMatrix,
                                        ProductPar, contraction,
                                        std::vector<unsigned int>, times,
                                        std::vector<ComplexD>, correlator);
    };
}

struct ContractorPar
{
    Contractor::GlobalPar                  global;
    std::vector<Contractor::A2AMatrixPar>  a2aMatrix;
    std::vector<Contractor::ProductPar>    product;
    std::vector<Contractor::EvalPar>       eval;
};

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

    for (unsigned int i = 0; i < terms.size() - 1; i++)
    {
        fileStem += terms[i] + "_" + std::to_string(result.times[i]) + "_";
    }
    fileStem += terms.back();
    if (!result.contraction.translationAverage)
    {
        fileStem += "_dt_" + std::to_string(dt);
    }
    filename = dir + "/" + ModuleBase::resultFilename(fileStem, traj);
    std::cout << "Saving correlator to '" << filename << "'" << std::endl;
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

// save correlation function to existing file
void saveCorrelator(const ComplexD *data,
                    const std::string filename,
                    const unsigned int tsize)
{
#ifdef HAVE_HDF5
    
    // Open existing file and dataset.
    H5::H5File file(filename, H5F_ACC_RDWR);
    
    hsize_t offset[1], count[1], stride[1], block[1], dimsm[1];
    
    H5::DataSet dataset = file.openDataSet("corrfunc_dset");
    
    // Specify location, size, and shape of subset to write.
    offset[0] = 0;
    count[0]  = tsize;
    stride[0] = 1;
    block[0]  = 1;
    
    // Define Memory Dataspace. Get file dataspace and select
    // a subset from the file dataspace.
    dimsm[0] = tsize;
    
    H5::DataSpace memspace(1, dimsm, NULL);
    
    H5::DataSpace dataspace = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    
    // Write a subset of data to the dataset
    dataset.write(data, Hdf5Type<ComplexD>::type(), memspace, dataspace);
    
#else
    HADRONS_ERROR(Implementation, "correlation function I/O needs HDF5 library");
#endif
}

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

    if (argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file> [Grid options]";
        std::cerr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    parFilename = argv[1];
    
    // initialization
    Grid_init(&argc, &argv);
    
    GridCartesian *Grid =
        SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                       GridDefaultSimd(4,1),
                                       GridDefaultMpi());

    // parse parameter file
    ContractorPar par;
    //    unsigned int  nMat,nCont;
    XmlReader     reader(parFilename);

    read(reader, "global",    par.global);
    read(reader, "a2aMatrix", par.a2aMatrix);
    read(reader, "product",   par.product);
    read(reader, "eval",      par.eval);

    // create diskvectors
    std::map<std::string, EigenDiskVector<ComplexD>> a2aMat;
    
    int localNt = Grid->LocalDimensions()[3];
    
    for (auto &p: par.a2aMatrix)
    {
        std::string dirName = par.global.diskVectorDir + std::to_string(Grid->ThisRank()) + "/" + p.name;
        a2aMat.emplace(p.name, EigenDiskVector<ComplexD>(dirName, localNt, p.cacheSize));
    }

    // trajectory loop
    for (unsigned int traj = par.global.trajCounter.start; 
         traj < par.global.trajCounter.end; traj += par.global.trajCounter.step)
    {
        std::cout << ":::::::: Trajectory " << traj << std::endl;

        // load evals
        Eigen::VectorXcd eval;
        
        std::string filename = par.eval[0].file;
        int size = par.eval[0].size;
        eval.resize(size);
        
        if(Grid->IsBoss()){
            std::cout << "======== Loading '" << filename << "'" << std::endl;
            A2AVectorsIo::loadEvalBlock(filename,eval,0,size);
        }
        Grid->Broadcast(0, eval.data(), size*sizeof(ComplexD));
        //Grid->Barrier();
        std::cout<<"past eval broadcast"<<std::endl;
        
        // load data
        for (auto &p: par.a2aMatrix)
        {
            std::string filename = p.file;
            double      t;//reading time
            
            tokenReplace(filename, "traj", traj);
            std::cout << "======== Loading '" << filename << "'" << std::endl;
            
            A2AMatrixIo<HADRONS_A2AM_IO_TYPE> a2aIo(filename, p.dataset, localNt);
            
            a2aIo.load(a2aMat.at(p.name), Grid, &t);
            std::cout << "Read " << a2aIo.getSize() << " bytes in " << t/1.0e6
            << " sec, " << a2aIo.getSize()/t*1.0e6/1024/1024 << " MB/s" << std::endl;
        }

        // contract
        EigenDiskVector<ComplexD>::Matrix buf;

        // loop over correlation functions
        for (auto &p: par.product)
        {
            std::vector<std::string>               term = strToVec<std::string>(p.terms);
            std::vector<std::set<unsigned int>>    times;
            std::vector<std::vector<unsigned int>> timeSeq;
            std::set<unsigned int>                 translations;
            std::vector<A2AMatrixTr<ComplexD>>     lastTerm(localNt);
            A2AMatrix<ComplexD>                    prod, buf, tmp;
            TimerArray                             tAr;
            double                                 fusec, busec, flops, bytes;
	        Contractor::CorrelatorResult           result;

            tAr.startTimer("Total");
            std::cout << "======== Contraction tr(";
            for (unsigned int g = 0; g < term.size(); ++g)
            {
                std::cout << term[g] << ((g == term.size() - 1) ? ')' : '*');
            }
            std::cout << std::endl;
            if (term.size() != p.times.size() + 1)
            {
                HADRONS_ERROR(Size, "number of terms (" + std::to_string(term.size()) 
                            + ") different from number of times (" 
                            + std::to_string(p.times.size() + 1) + ")");
            }
            for (auto &s: p.times)
            {
                times.push_back(parseTimeRange(s, par.global.nt));
            }
            for (auto &m: par.a2aMatrix)
            {
                if (std::find(result.a2aMatrix.begin(), result.a2aMatrix.end(), m) == result.a2aMatrix.end())
                {
                    result.a2aMatrix.push_back(m);
                    tokenReplace(result.a2aMatrix.back().file, "traj", traj);
                }
            }
            result.contraction = p;
            result.correlator.resize(par.global.nt, 0.);

            translations = parseTimeRange(p.translations, par.global.nt);
            makeTimeSeq(timeSeq, times);
            std::cout << timeSeq.size()*translations.size()*(term.size() - 2) << " A*B, "
                    << timeSeq.size()*translations.size()*par.global.nt << " tr(A*B)"
                    << std::endl;

            std::cout << "* Caching transposed last term" << std::endl;
            for (unsigned int t = 0; t < localNt; ++t)
            {
                tAr.startTimer("Disk vector overhead");
                const A2AMatrix<ComplexD> &ref = a2aMat.at(term.back())[t];
                tAr.stopTimer("Disk vector overhead");

                tAr.startTimer("Transpose caching");
                lastTerm[t].resize(ref.rows(), ref.cols());
                thread_for( j,ref.cols(),{
                  for (unsigned int i = 0; i < ref.rows(); ++i)
                  {
                      lastTerm[t](i, j) = ref(i, j);
                  }
                });
                tAr.stopTimer("Transpose caching");
            }
            uint64_t a2abytes = lastTerm[0].rows()*lastTerm[0].cols()*sizeof(HADRONS_A2AM_IO_TYPE);
            std::cout << Sec(tAr.getDTimer("Transpose caching")) << " " 
                      << Bytes(a2abytes*localNt, tAr.getDTimer("Transpose caching")) << std::endl;
            std::cout<<"a2abytes "<<a2abytes<<std::endl;
            for (unsigned int i = 0; i < timeSeq.size(); ++i)
            {
                unsigned int dti = 0;
                auto         &t = timeSeq[i];

                result.times = t;
                for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                {
                    result.correlator[tLast] = 0.;
                }
                for (auto &dt: translations)
                {
                    std::cout << "* Step " << i*translations.size() + dti + 1
                            << "/" << timeSeq.size()*translations.size()
                            << " -- positions= " << t << ", dt= " << dt << std::endl;
                    if (term.size() > 2)
                    {
                        std::cout << std::setw(8) << "products";
                    }
                    flops  = 0.;
                    bytes  = 0.;
                    fusec  = tAr.getDTimer("A*B algebra");
                    busec  = tAr.getDTimer("A*B total");
                    tAr.startTimer("Linear algebra");
                    tAr.startTimer("Disk vector overhead");
                    
                    // Is src time slice on node?
                    int srcNode = TIME_MOD(t[0] + dt)/localNt;
                    prod.resize(lastTerm[0].rows(), lastTerm[0].cols());
                    if(srcNode==Grid->ThisRank()){
                        prod = a2aMat.at(term[0])[(t[0] + dt)%localNt];
                    }
                    std::cout<<"t[0]= "<<t[0]<<" dt= "<<dt<<" src node= "<<srcNode<<std::endl;
                    Grid->Broadcast(srcNode, prod.data(), a2abytes);
                    
                    
                    
                    tAr.stopTimer("Disk vector overhead");
//                    for (unsigned int j = 1; j < term.size() - 1; ++j)
//                    {
//                        tAr.startTimer("Disk vector overhead");
//                        const A2AMatrix<ComplexD> &ref = a2aMat.at(term[j])[TIME_MOD(t[j] + dt)];
//                        tAr.stopTimer("Disk vector overhead");
//
//                        tAr.startTimer("A*B total");
//                        tAr.startTimer("A*B algebra");
//                        A2AContraction::mul(tmp, prod, ref);
//                        tAr.stopTimer("A*B algebra");
//                        flops += A2AContraction::mulFlops(prod, ref);
//                        prod   = tmp;
//                        tAr.stopTimer("A*B total");
//                        bytes += 3.*tmp.rows()*tmp.cols()*sizeof(ComplexD);
//                    }
//                    if (term.size() > 2)
//                    {
//                        std::cout << Sec(tAr.getDTimer("A*B total") - busec) << " "
//                                << Flops(flops, tAr.getDTimer("A*B algebra") - fusec) << " "
//                                << Bytes(bytes, tAr.getDTimer("A*B total") - busec) << std::endl;
//                    }
                    std::cout << std::setw(8) << "traces";
                    flops  = 0.;
                    bytes  = 0.;
                    fusec  = tAr.getDTimer("tr(A*B)");
                    busec  = tAr.getDTimer("tr(A*B)");
                    for (unsigned int tLast = 0; tLast < localNt; ++tLast)
                    {
                        tAr.startTimer("tr(A*B)");
                        A2AContraction::accTrMul(result.correlator[TIME_MOD(tLast+Grid->ThisRank()*localNt - dt)],
                                                 prod,
                                                 lastTerm[tLast],
                                                 eval);
                        tAr.stopTimer("tr(A*B)");
                        flops += A2AContraction::accTrMulCCFlops(prod, lastTerm[tLast]);
                        bytes += 2.*prod.rows()*prod.cols()*sizeof(ComplexD);
                    }
                    tAr.stopTimer("Linear algebra");
                    std::cout << Sec(tAr.getDTimer("tr(A*B)") - busec) << " "
                            << Flops(flops, tAr.getDTimer("tr(A*B)") - fusec) << " " 
                            << Bytes(bytes, tAr.getDTimer("tr(A*B)") - busec) << std::endl;
                
                    
                    if (!p.translationAverage)
                    {
                        std::string fileStem = "", filename;
                        std::vector<std::string> terms = strToVec<std::string>(result.contraction.terms);
                        for (unsigned int i = 0; i < terms.size() - 1; i++)
                        {
                            fileStem += terms[i] + "_" + std::to_string(result.times[i]) + "_";
                        }
                        fileStem += terms.back();
                        fileStem += "_dt_" + std::to_string(dt);
                        filename = par.global.output + "/" + ModuleBase::resultFilename(fileStem, traj);
                        std::cout << "Saving correlator to '" << filename << "'" << std::endl;
                        
                        
                        Grid->GlobalSumVector(result.correlator.data(), par.global.nt);
                        if(Grid->IsBoss()){
                            makeFileDir(par.global.output);
                            initCorrFile(filename, par.global.nt);
                            saveCorrelator(result.correlator.data(), filename, par.global.nt);
                        }
                        for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                        {
                            result.correlator[tLast] = 0.;
                        }
                    }
                    dti++;
                }
                if (p.translationAverage)
                {
                    std::string fileStem = "", filename;
                    std::vector<std::string> terms = strToVec<std::string>(result.contraction.terms);
                    for (unsigned int i = 0; i < terms.size() - 1; i++)
                    {
                        fileStem += terms[i] + "_" + std::to_string(result.times[i]) + "_";
                    }
                    fileStem += terms.back();
                    filename = par.global.output + "/" + ModuleBase::resultFilename(fileStem, traj);
                    std::cout << "Saving correlator to '" << filename << "'" << std::endl;
                    
                    
                    for (unsigned int tLast = 0; tLast < par.global.nt; ++tLast)
                    {
                        result.correlator[tLast] /= translations.size();
                    }
                    Grid->GlobalSumVector(result.correlator.data(), par.global.nt);
                    if(Grid->IsBoss()){
                        makeFileDir(par.global.output);
                        initCorrFile(filename, par.global.nt);
                        saveCorrelator(result.correlator.data(), filename, par.global.nt);
                    }
                }
            }
            tAr.stopTimer("Total");
            printTimeProfile(tAr.getTimings(), tAr.getTimer("Total"));
        }
    }
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
