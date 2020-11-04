
#ifndef PE_PROCESSOR_H
#define PE_PROCESSOR_H

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <functional>
#include <memory.h>
#include <unordered_map>
#include <map>
#include <future>
#include <deque>
#include "fastqreader.h"
#include "util.h"
#include "adaptertrimmer.h"
#include "basecorrector.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "polyx.h"
#include "options.h"
#include "threadconfig.h"
#include "filter.h"
#include "umiprocessor.h"
#include "overlapanalysis.h"
#include "writerthread.h"
#include "duplicate.h"
#include "read.h"
#include "transsearcher.hpp"
#include "bwtfmiDB.h"

using namespace std;

struct ReadPairPack
{
    ReadPair ** data;
    int         count;
};


typedef struct ReadPairPack ReadPairPack;


struct ReadPairRepository
{
    ReadPairPack ** packBuffer;
    atomic_long     readPos;
    atomic_long     writePos;

    // std::mutex mtx;
    // std::mutex readCounterMtx;
    // std::condition_variable repoNotFull;
    // std::condition_variable repoNotEmpty;
};


typedef struct ReadPairRepository ReadPairRepository;


class PairEndProcessor
{
    public:
        PairEndProcessor(Options *  opt,
                         BwtFmiDB * tbwtfmiDB);

        ~PairEndProcessor();

        bool process();

    private:
        bool processPairEnd(ReadPairPack *  pack,
                            ThreadConfig *  config,
                            TransSearcher * transSearcher);

        bool processRead(Read *     r,
                         ReadPair * originalRead,
                         bool       reversed);

        void initPackRepository();

        void destroyPackRepository();

        void producePack(ReadPairPack * pack);

        void consumePack(ThreadConfig *  config,
                         TransSearcher * transSearcher);

        void producerTask();

        void consumerTask(ThreadConfig *  config,
                          TransSearcher * transSearcher);

        void initConfig(ThreadConfig * config);

        void initOutput();

        void closeOutput();

        void statInsertSize(Read *          r1,
                            Read *          r2,
                            OverlapResult & ov,
                            int             frontTrimmed1 = 0,
                            int             frontTrimmed2 = 0);

        int getPeakInsertSize();

        void writeTask(WriterThread * config);

        void S2FReport();

    private:
        ReadPairRepository                                                   mRepo;
        atomic_bool                                                          mProduceFinished;
        atomic_int                                                           mFinishedThreads;
        std::mutex                                                           mOutputMtx;
        std::mutex                                                           mInputMtx;
        std::mutex                                                           mReadsMtx;
        std::mutex                                                           logMtx;
        Options *                                                            mOptions;
        Filter *                                                             mFilter;
        gzFile                                                               mZipFile1;
        gzFile                                                               mZipFile2;
        ofstream *                                                           mOutStream1;
        ofstream *                                                           mOutStream2;
        UmiProcessor *                                                       mUmiProcessor;
        long *                                                               mInsertSizeHist;
        WriterThread *                                                       mLeftWriter;
        WriterThread *                                                       mRightWriter;
        WriterThread *                                                       mUnpairedLeftWriter;
        WriterThread *                                                       mUnpairedRightWriter;
        WriterThread *                                                       mMergedWriter;
        WriterThread *                                                       mFailedWriter;
        Duplicate *                                                          mDuplicate;
        BwtFmiDB *                                                           tbwtfmiDB;
        std::string                                                          fileoutname;
		std::vector<std::pair<std::string, int> >			   sortedSpeciesFreqVector;
};
#endif


//~ Formatted by Jindent --- http://www.jindent.com
