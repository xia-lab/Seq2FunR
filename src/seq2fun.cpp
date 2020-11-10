#include <stdio.h>
#include "fastqreader.h"
#include "unittest.h"
#include <time.h>
#include "cmdline.h"
#include <sstream>
#include "util.h"
#include "options.h"
#include "processor.h"
#include "evaluator.h"
#include <map>
#include <fstream>
#include <cmath>
#include <thread>
#include <vector>
#include <unordered_map>
#include <utility>
#include <valarray> 
#include "bwtfmiDB.h"

std::string command;
std::mutex logmtx;
//' @export
// [[Rcpp::export]]

int seq2fun(std::string sampletable,
             std::string genemap,
             std::string tfmi,
             bool outputMappedCleanReads = false,
             bool profiling = false,
             std::string in1 = "",
             std::string in2 = "",
             std::string prefix = "",
             std::string mode = "tGREEDY",
             int mismatch = 2,
             int minscore = 65,
             int minlength = 11,
             int maxtranslength = 65,
             int nThreads = 8,
             bool verbose = true){
  time_t t_begin = time(NULL);
  
  Options opt;
  
  opt.thread = nThreads;
  int n_t = std::thread::hardware_concurrency();
  opt.thread = std::min(std::min(opt.thread, 16), n_t);
  
  opt.compression = 4;
  opt.readsToProcess = 0;
  opt.phred64 = false;
  opt.verbose = verbose;
  opt.debug = false;
  
  opt.adapter.enabled = true;
  opt.adapter.detectAdapterForPE = false;
  opt.adapter.sequence = "auto";
  opt.adapter.sequenceR2 = "auto";
  opt.adapter.fastaFile = "";
  
  
  opt.adapter.polyA = true;
  opt.trim.front1 = false;
  opt.trim.tail1 = false;
  opt.trim.maxLen1 = false;
  
  opt.trim.front2 = false;
  
  opt.trim.tail2 = false;
  opt.trim.maxLen2 = false;
  
  opt.polyGTrim.enabled = false;
  opt.polyGTrim.minLen = 10;
  opt.polyXTrim.enabled = false;
  opt.polyXTrim.minLen = 10;
  
  opt.qualityCut.enabledFront = false;
  opt.qualityCut.enabledRight = false;
  opt.qualityCut.enabledTail = false;
  
  opt.qualityCut.windowSizeShared = 4;
  opt.qualityCut.qualityShared = 20;
  opt.qualityCut.windowSizeFront = 4;
  opt.qualityCut.qualityFront = 20;
  opt.qualityCut.windowSizeTail = 4;
  opt.qualityCut.qualityTail = 20;
  opt.qualityCut.windowSizeRight = 4;
  opt.qualityCut.qualityRight = 20;
  
  
  opt.qualfilter.enabled = true;
  opt.qualfilter.qualifiedQual = num2qual(15);
  opt.qualfilter.unqualifiedPercentLimit = 40;
  opt.qualfilter.avgQualReq = 0;
  opt.qualfilter.nBaseLimit = 5;
  opt.lengthFilter.enabled  = true;
  opt.lengthFilter.requiredLength = max(60, static_cast<int> (opt.transSearch.minAAFragLength) * 3);
  opt.lengthFilter.maxLength = 0;
  opt.complexityFilter.enabled = true;
  opt.complexityFilter.threshold = (min(100, max(0, 30))) / 100.0;
  opt.correction.enabled = true;
  opt.overlapRequire = 30;
  opt.overlapDiffLimit  = 5;
  opt.overlapDiffPercentLimit = 20;
  
  opt.umi.enabled  = false;
  opt.umi.length  = 0;
  opt.umi.prefix  = "";
  opt.umi.skip  = 0;
  
  opt.overRepAnalysis.enabled  = false;
  opt.overRepAnalysis.sampling  = 20;
  std::string blacklist1  = "";
  std::string blacklist2  = "";
  int indexFilterThreshold  = 0;
  opt.initIndexFiltering(blacklist1, blacklist2, indexFilterThreshold);
  
  opt.mHomoSearchOptions.genemap = genemap;
  opt.mHomoSearchOptions.profiling = profiling;
  
  opt.transSearch.tmode = mode;
  
  
  if(opt.transSearch.tmode == "tGREEDY"){
    opt.transSearch.mode = tGREEDY;
  } else {
    opt.transSearch.mode = tMEM;
  } 
  
  opt.transSearch.misMatches = mismatch;
  opt.transSearch.minScore = minscore;
  opt.transSearch.minAAFragLength = minlength;
  opt.transSearch.maxTransLength = maxtranslength;
  opt.transSearch.maxTransLength = max(opt.transSearch.maxTransLength, opt.transSearch.minAAFragLength);
  opt.transSearch.maxTransLength = min((unsigned)60, opt.transSearch.maxTransLength);   
  
  opt.transSearch.tfmi = tfmi;
  
  //read all database tables, maps;
  opt.readDB();
  
  BwtFmiDB * tbwtfmiDB = new BwtFmiDB(& opt);
  
  // I/O
  opt.mHomoSearchOptions.sampleTable = sampletable;
  opt.mHomoSearchOptions.prefix = prefix;
  make_dir_con(opt.mHomoSearchOptions.prefix);
  
  opt.outputMappedCleanReads = outputMappedCleanReads;
  
  if (opt.mHomoSearchOptions.prefix.empty() && opt.mHomoSearchOptions.sampleTable.empty()) {
    error_exit("You must specify output file using --prefix or using --sampletable, which contains prefix string");
  } else if (!opt.mHomoSearchOptions.prefix.empty() && !opt.mHomoSearchOptions.sampleTable.empty()) {
    error_exit("You must specify output file using either --prefix or --sampletable");
  }

	//for single file;
	if (!opt.mHomoSearchOptions.prefix.empty()) {
		opt.in1 = in1;
		opt.in2 = in2;
		std::string outFName;
		opt.jsonFile = opt.mHomoSearchOptions.prefix + "_report.json";
		if (opt.outputMappedCleanReads) {
			opt.out1 = opt.mHomoSearchOptions.prefix + "_mapped_R1.fastq.gz";
			if (opt.isPaired()) opt.out2 = opt.mHomoSearchOptions.prefix + "_mapped_R2.fastq.gz";
		}

		stringstream ss;
		command = ss.str();

		bool supportEvaluation = !opt.inputFromSTDIN && opt.in1 != "/dev/stdin";
		Evaluator eva(&opt);
		if (supportEvaluation) {
			eva.evaluateSeqLen();
			if (opt.overRepAnalysis.enabled)
				eva.evaluateOverRepSeqs();
		}

		long readNum = 0;

		// using evaluator to guess how many reads in total
		if (opt.shallDetectAdapter(false)) {
			if (!supportEvaluation)
				cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
			else {
				cerr << "Detecting adapter sequence for read1..." << endl;
				string adapt = eva.evalAdapterAndReadNum(readNum, false);
				if (adapt.length() > 60)
					adapt.resize(0, 60);
				if (adapt.length() > 0) {
					opt.adapter.sequence = adapt;
					opt.adapter.detectedAdapter1 = adapt;
				} else {
					cerr << "No adapter detected for read1" << endl;
					opt.adapter.sequence = "";
				}
				cerr << endl;
			}
		}
		if (opt.shallDetectAdapter(true)) {
			if (!supportEvaluation)
				cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
			else {
				cerr << "Detecting adapter sequence for read2..." << endl;
				string adapt = eva.evalAdapterAndReadNum(readNum, true);
				if (adapt.length() > 60)
					adapt.resize(0, 60);
				if (adapt.length() > 0) {
					opt.adapter.sequenceR2 = adapt;
					opt.adapter.detectedAdapter2 = adapt;
				} else {
					cerr << "No adapter detected for read2" << endl;
					opt.adapter.sequenceR2 = "";
				}
				cerr << endl;
			}
		}

		opt.validate();

		// using evaluator to guess how many reads in total
		if (opt.split.needEvaluation && supportEvaluation) {
			// if readNum is not 0, means it is already evaluated by other functions
			if (readNum == 0) {
				eva.evaluateReadNum(readNum);
			}
			opt.split.size = readNum / opt.split.number;
			// one record per file at least
			if (opt.split.size <= 0) {
				opt.split.size = 1;
				cerr << "WARNING: the input file has less reads than the number of files to split" << endl;
			}
		}

		// using evaluator to check if it's two color system
		bool trim_poly_g = false;
		bool disable_trim_poly_g = false;
		if (!trim_poly_g && !disable_trim_poly_g && supportEvaluation) {
		  bool twoColorSystem = eva.isTwoColorSystem();
		  if (twoColorSystem) {
		    opt.polyGTrim.enabled = true;
		  }
		}
		
		opt.transSearch.transReadsTaxonKOPairDeque.clear();
		opt.transSearch.transSearchMappedReads = 0;
		opt.mHomoSearchOptions.totalOrigReads = 0;

		Processor p(&opt);
		p.process(tbwtfmiDB);
		time_t t_finished = time(NULL);
		cerr << endl << command << endl;
		cerr << endl << "Seq2Fun v" << SEQ2FUNR_VER << ", time used: " << (t_finished) - t_begin << " seconds, mapping " << opt.transSearch.transSearchMappedReads << " reads out of " << opt.mHomoSearchOptions.totalOrigReads << " (" << getPercentage(opt.transSearch.transSearchMappedReads, opt.mHomoSearchOptions.totalOrigReads) << " %)" << endl << endl;
		opt.transSearch.transReadsTaxonKOPairDeque.clear();
		opt.transSearch.transSearchMappedReads = 0;
		opt.mHomoSearchOptions.totalOrigReads = 0;
	} else {
		opt.parseSampleTable();
		std::vector< std::pair<std::string, std::unordered_map<std::string, int> > > sampleSpeciesTableVec;
		sampleSpeciesTableVec.reserve(opt.samples.size());
		int count = 0;
		for (auto & it : opt.samples) {
			time_t t_cycleBegin = time(NULL);
			count++;
			std::stringstream msgSS;
			msgSS << "processing sample: " << it.prefix << ", " << count << " out of " << opt.samples.size() << " samples" << "\n";
			loginfo(msgSS.str());
			opt.mHomoSearchOptions.prefix = it.prefix;
			make_dir_con(opt.mHomoSearchOptions.prefix);
			opt.in1 = it.in1;
			opt.in2 = it.in2;
			opt.jsonFile = it.prefix +  "_report.json";
			
			if (opt.outputMappedCleanReads) {
				opt.out1 = it.prefix + "_mapped_R1.fastq.gz";
				if (opt.isPaired()) opt.out2 = it.prefix + "_mapped_R2.fastq.gz";
			}
			std::stringstream ss;
			command = ss.str();

			bool supportEvaluation = !opt.inputFromSTDIN && opt.in1 != "/dev/stdin";
			Evaluator eva(&opt);
			if (supportEvaluation) {
				eva.evaluateSeqLen();
				if (opt.overRepAnalysis.enabled)
					eva.evaluateOverRepSeqs();
			}

			long readNum = 0;

			// using evaluator to guess how many reads in total
			if (opt.shallDetectAdapter(false)) {
				if (!supportEvaluation)
					cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
				else {
					cerr << "Detecting adapter sequence for read1..." << endl;
					string adapt = eva.evalAdapterAndReadNum(readNum, false);
					if (adapt.length() > 60)
						adapt.resize(0, 60);
					if (adapt.length() > 0) {
						opt.adapter.sequence = adapt;
						opt.adapter.detectedAdapter1 = adapt;
					} else {
						cerr << "No adapter detected for read1" << endl;
						opt.adapter.sequence = "";
					}
					cerr << endl;
				}
			}
			if (opt.shallDetectAdapter(true)) {
				if (!supportEvaluation)
					cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
				else {
					cerr << "Detecting adapter sequence for read2..." << endl;
					string adapt = eva.evalAdapterAndReadNum(readNum, true);
					if (adapt.length() > 60)
						adapt.resize(0, 60);
					if (adapt.length() > 0) {
						opt.adapter.sequenceR2 = adapt;
						opt.adapter.detectedAdapter2 = adapt;
					} else {
						cerr << "No adapter detected for read2" << endl;
						opt.adapter.sequenceR2 = "";
					}
					cerr << endl;
				}
			}

			opt.validate();

			// using evaluator to guess how many reads in total
			if (opt.split.needEvaluation && supportEvaluation) {
				// if readNum is not 0, means it is already evaluated by other functions
				if (readNum == 0) {
					eva.evaluateReadNum(readNum);
				}
				opt.split.size = readNum / opt.split.number;
				// one record per file at least
				if (opt.split.size <= 0) {
					opt.split.size = 1;
					cerr << "WARNING: the input file has less reads than the number of files to split" << endl;
				}
			}

			// using evaluator to check if it's two color system
			bool trim_poly_g = false;
			bool disable_trim_poly_g = false;
			if (!trim_poly_g && !disable_trim_poly_g && supportEvaluation) {
			  bool twoColorSystem = eva.isTwoColorSystem();
			  if (twoColorSystem) {
			    opt.polyGTrim.enabled = true;
			  }
			}
			
			opt.transSearch.sampleSpeciesAbunUMap.clear();
			opt.transSearch.transReadsTaxonKOPairDeque.clear();
			opt.transSearch.transSearchMappedReads = 0;
			opt.mHomoSearchOptions.totalOrigReads = 0;

			Processor p(&opt);
			p.process(tbwtfmiDB);
			time_t t_cycleEnd = time(NULL);
			cerr << endl << command << endl;
			cerr << endl << "Seq2Fun v" << SEQ2FUNR_VER << ", time used: " << (t_cycleEnd) - t_cycleBegin << " seconds, mapping " << opt.transSearch.transSearchMappedReads << " reads out of " << opt.mHomoSearchOptions.totalOrigReads << " (" << getPercentage(opt.transSearch.transSearchMappedReads, opt.mHomoSearchOptions.totalOrigReads) << " %)" << endl << endl;
			sampleSpeciesTableVec.push_back(std::make_pair(it.prefix, opt.transSearch.sampleSpeciesAbunUMap));
			opt.transSearch.sampleSpeciesAbunUMap.clear();
			opt.transSearch.transReadsTaxonKOPairDeque.clear();
			opt.transSearch.transSearchMappedReads = 0;
			opt.mHomoSearchOptions.totalOrigReads = 0;
		}

		writeSampleSpeciesTable(sampleSpeciesTableVec, opt);
		sampleSpeciesTableVec.clear();
		time_t t_total = time(NULL);
		cerr << endl << "Seq2Fun v" << SEQ2FUNR_VER << ", time used: " << (t_total) - t_begin << " seconds, processed " << opt.samples.size() << " samples" << endl << endl;
	}
	if (tbwtfmiDB) delete tbwtfmiDB;
	return 0;
}
