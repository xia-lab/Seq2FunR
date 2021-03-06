#include "seprocessor.h"

SingleEndProcessor::SingleEndProcessor(Options* opt, BwtFmiDB * tbwtfmiDB) {
	mOptions = opt;
	mProduceFinished = false;
	mFinishedThreads = 0;
	mFilter = new Filter(opt);
	mOutStream = NULL;
	mZipFile = NULL;
	mUmiProcessor = new UmiProcessor(opt);
	mLeftWriter = NULL;
	mFailedWriter = NULL;

	mDuplicate = NULL;
	if (mOptions->duplicate.enabled) {
		mDuplicate = new Duplicate(mOptions);
	}
	this->tbwtfmiDB = tbwtfmiDB;
	fileoutname.clear();
	sortedSpeciesFreqVector.clear();
}

SingleEndProcessor::~SingleEndProcessor() {
	delete mFilter;
	if (mDuplicate) {
		delete mDuplicate;
		mDuplicate = NULL;
	}
	if (mUmiProcessor) {
		delete mUmiProcessor;
		mUmiProcessor = NULL;
	}
}

void SingleEndProcessor::initOutput() {
	if (!mOptions->failedOut.empty())
		mFailedWriter = new WriterThread(mOptions, mOptions->failedOut);
	if (mOptions->out1.empty())
		return;
	mLeftWriter = new WriterThread(mOptions, mOptions->out1);
}

void SingleEndProcessor::closeOutput() {
	if (mLeftWriter) {
		delete mLeftWriter;
		mLeftWriter = NULL;
	}
	if (mFailedWriter) {
		delete mFailedWriter;
		mFailedWriter = NULL;
	}
}

void SingleEndProcessor::initConfig(ThreadConfig* config) {
	if (mOptions->out1.empty())
		return;

	if (mOptions->split.enabled) {
		config->initWriterForSplit();
	}
}

bool SingleEndProcessor::process() {
	if (!mOptions->split.enabled)
		initOutput();

	initPackRepository();
	std::thread producer(std::bind(&SingleEndProcessor::producerTask, this));

	//TODO: get the correct cycles
	int cycle = 151;
	ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
	TransSearcher** transSearchers = new TransSearcher*[mOptions->thread];
	for (int t = 0; t < mOptions->thread; t++) {
		configs[t] = new ThreadConfig(mOptions, t, false);
		initConfig(configs[t]);
		transSearchers[t] = new TransSearcher(tbwtfmiDB, mOptions);
	}
	std::thread** threads = new thread*[mOptions->thread];
	for (int t = 0; t < mOptions->thread; t++) {
		threads[t] = new std::thread(std::bind(&SingleEndProcessor::consumerTask, this, configs[t], transSearchers[t]));
	}

	std::thread* leftWriterThread = NULL;
	std::thread* failedWriterThread = NULL;
	if (mLeftWriter)
		leftWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mLeftWriter));
	if (mFailedWriter)
		failedWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mFailedWriter));

	producer.join();
	for (int t = 0; t < mOptions->thread; t++) {
		threads[t]->join();
	}

	if (!mOptions->split.enabled) {
		if (leftWriterThread)
			leftWriterThread->join();
		if (failedWriterThread)
			failedWriterThread->join();
	}

	if (mOptions->verbose)
		loginfo("start to generate reports\n");

	// merge stats and read filter results
	vector<Stats*> preStats;
	vector<Stats*> postStats;
	vector<FilterResult*> filterResults;
	for (int t = 0; t < mOptions->thread; t++) {
		preStats.push_back(configs[t]->getPreStats1());
		postStats.push_back(configs[t]->getPostStats1());
		filterResults.push_back(configs[t]->getFilterResult());
	}
	Stats* finalPreStats = Stats::merge(preStats);
	Stats* finalPostStats = Stats::merge(postStats);
	FilterResult* finalFilterResult = FilterResult::merge(filterResults);

	mOptions->mHomoSearchOptions.totalOrigReads = finalPreStats->getReads();

	// read filter results to the first thread's
	for (int t = 1; t < mOptions->thread; t++) {
		preStats.push_back(configs[t]->getPreStats1());
		postStats.push_back(configs[t]->getPostStats1());
	}

	int* dupHist = NULL;
	double* dupMeanTlen = NULL;
	double* dupMeanGC = NULL;
	double dupRate = 0.0;
	if (mOptions->duplicate.enabled) {
		dupHist = new int[mOptions->duplicate.histSize];
		memset(dupHist, 0, sizeof (int) * mOptions->duplicate.histSize);
		dupMeanGC = new double[mOptions->duplicate.histSize];
		memset(dupMeanGC, 0, sizeof (double) * mOptions->duplicate.histSize);
		dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
		cerr << endl;
		cerr << "Duplication rate (may be overestimated since this is SE data): " << dupRate * 100.0 << "%" << endl;
	}
	
	S2FReport();
	sortedSpeciesFreqVector.clear();
	
	// make JSON report
	JsonReporter jr(mOptions);
	jr.setDupHist(dupHist, dupMeanGC, dupRate);
	jr.report(finalFilterResult, finalPreStats, finalPostStats);
	
		// clean up
	for (int t = 0; t < mOptions->thread; t++) {
		delete threads[t];
		threads[t] = NULL;
		delete configs[t];
		configs[t] = NULL;
		delete transSearchers[t];
		transSearchers[t] = NULL;
	}

	delete finalPreStats;
	delete finalPostStats;
	delete finalFilterResult;

	if (mOptions->duplicate.enabled) {
		delete[] dupHist;
		delete[] dupMeanGC;
	}
	delete[] threads;
	delete[] configs;
	delete[] transSearchers;
	if (leftWriterThread)
		delete leftWriterThread;
	if (failedWriterThread)
		delete failedWriterThread;

	if (!mOptions->split.enabled)
		closeOutput();
	return true;
}

bool SingleEndProcessor::processSingleEnd(ReadPack* pack, ThreadConfig* config, TransSearcher * transSearcher) {
	string outstr;
	string failedOut;
	int readPassed = 0;
	std::deque<std::pair<std::string, TaxonKO> > tmpReadsTaxonKODeque;

	for (int p = 0; p < pack->count; p++) {

		// original read1
		Read* or1 = pack->data[p];

		// stats the original read before trimming
		config->getPreStats1()->statRead(or1);

		// handling the duplication profiling
		if (mDuplicate)
			mDuplicate->statRead(or1);

		// filter by index
		if (mOptions->indexFilter.enabled && mFilter->filterByIndex(or1)) {
			delete or1;
			continue;
		}

		// umi processing
		if (mOptions->umi.enabled)
			mUmiProcessor->process(or1);

		int frontTrimmed = 0;
		// trim in head and tail, and apply quality cut in sliding window
		Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1, frontTrimmed);

		if (r1 != NULL) {
			if (mOptions->polyGTrim.enabled)
				PolyX::trimPolyG(r1, config->getFilterResult(), mOptions->polyGTrim.minLen);
		}

		if (r1 != NULL && mOptions->adapter.enabled) {
			bool trimmed = false;
			if (mOptions->adapter.hasSeqR1)
				trimmed = AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.sequence, false);
			bool incTrimmedCounter = !trimmed;
			if (mOptions->adapter.hasFasta) {
				AdapterTrimmer::trimByMultiSequences(r1, config->getFilterResult(), mOptions->adapter.seqsInFasta, false, incTrimmedCounter);
			}

			if (mOptions->adapter.polyA) {
				AdapterTrimmer::trimPolyA(r1, config->getFilterResult(), false, incTrimmedCounter);
			}
		}

		if (r1 != NULL) {
			if (mOptions->polyXTrim.enabled)
				PolyX::trimPolyX(r1, config->getFilterResult(), mOptions->polyXTrim.minLen);
		}

		if (r1 != NULL) {
			if (mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length())
				r1->resize(mOptions->trim.maxLen1);
		}

		int result = mFilter->passFilter(r1);

		config->addFilterResult(result, 1);

		if (r1 != NULL && result == PASS_FILTER) {
			std::pair<std::string, TaxonKO> tmpReadsTaxonKOPair;
			transSearcher->transSearch(r1, tmpReadsTaxonKOPair);

			if (!tmpReadsTaxonKOPair.first.empty()) {
				tmpReadsTaxonKODeque.push_back(tmpReadsTaxonKOPair);
				if (mLeftWriter) {
					outstr += r1->toStringWithTag(tmpReadsTaxonKOPair.second.getToSpecies());
				}
			}

			// stats the read after filtering
			config->getPostStats1()->statRead(r1);
			readPassed++;
		} else if (mFailedWriter) {
			failedOut += or1->toStringWithTag(FAILED_TYPES[result]);
		}


		delete or1;
		// if no trimming applied, r1 should be identical to or1
		if (r1 != or1 && r1 != NULL)
			delete r1;
	}

	mReadsMtx.lock();
	mOptions->transSearch.transReadsTaxonKOPairDeque.insert(mOptions->transSearch.transReadsTaxonKOPairDeque.begin(),
			tmpReadsTaxonKODeque.begin(),
			tmpReadsTaxonKODeque.end());
	mReadsMtx.unlock();
	tmpReadsTaxonKODeque.clear();

	if (mOptions->verbose) {
		logMtx.lock();
		int rCount = 0;
		rCount = mOptions->transSearch.transReadsTaxonKOPairDeque.size();
		logMtx.unlock();
		std::string str = "Mapped " + std::to_string(rCount) + " reads";
		loginfo(str);
	}
	// if splitting output, then no lock is need since different threads write different files
	if (!mOptions->split.enabled)
		mOutputMtx.lock();
	if (mOptions->outputToSTDOUT) {
		fwrite(outstr.c_str(), 1, outstr.length(), stdout);
	} else if (mOptions->split.enabled) {
		// split output by each worker thread
		if (!mOptions->out1.empty())
			config->getWriter1()->writeString(outstr);
	}

	if (mLeftWriter) {
		char* ldata = new char[outstr.size()];
		memcpy(ldata, outstr.c_str(), outstr.size());
		mLeftWriter->input(ldata, outstr.size());
	}
	if (mFailedWriter && !failedOut.empty()) {
		// write failed data
		char* fdata = new char[failedOut.size()];
		memcpy(fdata, failedOut.c_str(), failedOut.size());
		mFailedWriter->input(fdata, failedOut.size());
	}
	if (!mOptions->split.enabled)
		mOutputMtx.unlock();

	if (mOptions->split.byFileLines)
		config->markProcessed(readPassed);
	else
		config->markProcessed(pack->count);

	delete pack->data;
	delete pack;

	return true;
}

void SingleEndProcessor::initPackRepository() {
	mRepo.packBuffer = new ReadPack*[PACK_NUM_LIMIT];
	memset(mRepo.packBuffer, 0, sizeof (ReadPack*) * PACK_NUM_LIMIT);
	mRepo.writePos = 0;
	mRepo.readPos = 0;
	//mRepo.readCounter = 0;

}

void SingleEndProcessor::destroyPackRepository() {
	delete mRepo.packBuffer;
	mRepo.packBuffer = NULL;
}

void SingleEndProcessor::producePack(ReadPack* pack) {
	//std::unique_lock<std::mutex> lock(mRepo.mtx);
	/*while(((mRepo.writePos + 1) % PACK_NUM_LIMIT)
		== mRepo.readPos) {
		//mRepo.repoNotFull.wait(lock);
	}*/

	mRepo.packBuffer[mRepo.writePos] = pack;
	mRepo.writePos++;

	/*if (mRepo.writePos == PACK_NUM_LIMIT)
		mRepo.writePos = 0;*/

	//mRepo.repoNotEmpty.notify_all();
	//lock.unlock();
}

void SingleEndProcessor::consumePack(ThreadConfig* config, TransSearcher * transSearcher) {
	ReadPack* data;
	//std::unique_lock<std::mutex> lock(mRepo.mtx);
	// buffer is empty, just wait here.
	/*while(mRepo.writePos % PACK_NUM_LIMIT == mRepo.readPos % PACK_NUM_LIMIT) {
		if(mProduceFinished){
			//lock.unlock();
			return;
		}
		//mRepo.repoNotEmpty.wait(lock);
	}*/

	mInputMtx.lock();
	while (mRepo.writePos <= mRepo.readPos) {
		usleep(1000);
		if (mProduceFinished) {
			mInputMtx.unlock();
			return;
		}
	}
	data = mRepo.packBuffer[mRepo.readPos];
	mRepo.readPos++;

	/*if (mRepo.readPos >= PACK_NUM_LIMIT)
		mRepo.readPos = 0;*/
	mInputMtx.unlock();

	//lock.unlock();
	//mRepo.repoNotFull.notify_all();

	processSingleEnd(data, config, transSearcher);

}

void SingleEndProcessor::producerTask() {
	if (mOptions->verbose)
		loginfo("start to load data");
	long lastReported = 0;
	int slept = 0;
	long readNum = 0;
	bool splitSizeReEvaluated = false;
	Read** data = new Read*[PACK_SIZE];
	memset(data, 0, sizeof (Read*) * PACK_SIZE);
	FastqReader reader(mOptions->in1, true, mOptions->phred64);
	int count = 0;
	bool needToBreak = false;
	while (true) {
		Read* read = reader.read();
		// TODO: put needToBreak here is just a WAR for resolve some unidentified dead lock issue 
		if (!read || needToBreak) {
			// the last pack
			ReadPack* pack = new ReadPack;
			pack->data = data;
			pack->count = count;
			producePack(pack);
			data = NULL;
			if (read) {
				delete read;
				read = NULL;
			}
			break;
		}
		data[count] = read;
		count++;
		// configured to process only first N reads
		if (mOptions->readsToProcess > 0 && count + readNum >= mOptions->readsToProcess) {
			needToBreak = true;
		}
		if (mOptions->verbose && count + readNum >= lastReported + 1000000) {
			lastReported = count + readNum;
			string msg = "loaded " + to_string((lastReported / 1000000)) + "M reads";
			loginfo(msg);
		}
		// a full pack
		if (count == PACK_SIZE || needToBreak) {
			ReadPack* pack = new ReadPack;
			pack->data = data;
			pack->count = count;
			producePack(pack);
			//re-initialize data for next pack
			data = new Read*[PACK_SIZE];
			memset(data, 0, sizeof (Read*) * PACK_SIZE);
			// if the consumer is far behind this producer, sleep and wait to limit memory usage
			while (mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT) {
				//cerr<<"sleep"<<endl;
				slept++;
				usleep(1000);
			}
			readNum += count;
			// if the writer threads are far behind this producer, sleep and wait
			// check this only when necessary
			if (readNum % (PACK_SIZE * PACK_IN_MEM_LIMIT) == 0 && mLeftWriter) {
				while (mLeftWriter->bufferLength() > PACK_IN_MEM_LIMIT) {
					slept++;
					usleep(1000);
				}
			}
			// reset count to 0
			count = 0;
			// re-evaluate split size
			// TODO: following codes are commented since it may cause threading related conflicts in some systems
			/*if(mOptions->split.needEvaluation && !splitSizeReEvaluated && readNum >= mOptions->split.size) {
				splitSizeReEvaluated = true;
				// greater than the initial evaluation
				if(readNum >= 1024*1024) {
					size_t bytesRead;
					size_t bytesTotal;
					reader.getBytes(bytesRead, bytesTotal);
					mOptions->split.size *=  (double)bytesTotal / ((double)bytesRead * (double) mOptions->split.number);
					if(mOptions->split.size <= 0)
						mOptions->split.size = 1;
				}
			}*/
		}
	}

	//std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
	mProduceFinished = true;
	if (mOptions->verbose)
		loginfo("all reads loaded, start to monitor thread status");
	//lock.unlock();

	// if the last data initialized is not used, free it
	if (data != NULL)
		delete[] data;
}

void SingleEndProcessor::consumerTask(ThreadConfig* config, TransSearcher * transSearcher) {
	while (true) {
		if (config->canBeStopped()) {
			mFinishedThreads++;
			break;
		}
		while (mRepo.writePos <= mRepo.readPos) {
			if (mProduceFinished)
				break;
			usleep(1000);
		}
		//std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
		if (mProduceFinished && mRepo.writePos == mRepo.readPos) {
			mFinishedThreads++;
			if (mOptions->verbose) {
				string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
				loginfo(msg);
			}
			//lock.unlock();
			break;
		}
		if (mProduceFinished) {
			if (mOptions->verbose) {
				string msg = "thread " + to_string(config->getThreadId() + 1) + " is processing the " + to_string(mRepo.readPos) + " / " + to_string(mRepo.writePos) + " pack";
				loginfo(msg);
			}
			consumePack(config, transSearcher);
			//lock.unlock();
		} else {
			//lock.unlock();
			consumePack(config, transSearcher);
		}
	}

	if (mFinishedThreads == mOptions->thread) {
		if (mLeftWriter)
			mLeftWriter->setInputCompleted();
		if (mFailedWriter)
			mFailedWriter->setInputCompleted();
	}

	if (mOptions->verbose) {
		string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
		loginfo(msg);
	}
}

void SingleEndProcessor::writeTask(WriterThread* config) {
	while (true) {
		if (config->isCompleted()) {
			// last check for possible threading related issue
			config->output();
			break;
		}
		config->output();
	}

	if (mOptions->verbose) {
		string msg = config->getFilename() + " writer finished";
		loginfo(msg);
	}
}

void SingleEndProcessor::S2FReport() {
	if (mOptions->mHomoSearchOptions.profiling && mOptions->mHomoSearchOptions.prefix.size() != 0) {
		fileoutname.clear();
		fileoutname = mOptions->mHomoSearchOptions.prefix + "_read_taxon.txt";
		std::ofstream * fout = new std::ofstream();
		fout->open(fileoutname.c_str(), std::ofstream::out);
		if (!fout->is_open()) error_exit("Can not open read taxon gene map file: " + fileoutname);
		if (mOptions->verbose) loginfo("Starting to write read taxon mapping table");

		*fout << "Read_id" << "\t" << "Phylum" << "\t" << "Class" << "\t" << "Order" << "\t" << "Family" << "\t" << "Genus" << "\t" << "Species" << "\n";
		for (auto & it : mOptions->transSearch.transReadsTaxonKOPairDeque) {
			*fout << it.first << "\t" << it.second.getToSpecies() << "\n";
		}
		fout->flush();
		fout->close();
		if (fout) delete fout;
		if (mOptions->verbose) loginfo("Finish to write read taxon mapping table");
	}
	mOptions->transSearch.transSearchMappedReads = 0;
	mOptions->transSearch.transSearchMappedReads = mOptions->transSearch.transReadsTaxonKOPairDeque.size();

	//for microbiome;
	std::unordered_map<std::string, int> tmpPhylumFreqUMap;
	std::unordered_map<std::string, int> tmpClassmFreqUMap;
	std::unordered_map<std::string, int> tmpOrderFreqUMap;
	std::unordered_map<std::string, int> tmpFamilyFreqUMap;
	std::unordered_map<std::string, int> tmpGenusFreqUMap;
	std::unordered_map<std::string, int> tmpSpeciesFreqUMap;

	long totalReadsMappedTaxon = 0;
	for (auto & it : mOptions->transSearch.transReadsTaxonKOPairDeque) {
		if (it.second.phylum != "unclassified") {
			totalReadsMappedTaxon++;
			tmpPhylumFreqUMap[it.second.phylum]++;

			auto classm = it.second.phylum + "\t" + it.second.classm;
			tmpClassmFreqUMap[classm]++;

			auto order = it.second.phylum + "\t" + it.second.classm + "\t" + it.second.order;
			tmpOrderFreqUMap[order]++;

			auto family = it.second.phylum + "\t" + it.second.classm + "\t" + it.second.order + "\t" + it.second.family;
			tmpFamilyFreqUMap[family]++;

			auto genus = it.second.phylum + "\t" + it.second.classm + "\t" + it.second.order + "\t" + it.second.family + "\t" + it.second.genus;
			tmpGenusFreqUMap[genus]++;

			auto species = it.second.phylum + "\t" + it.second.classm + "\t" + it.second.order + "\t" + it.second.family + "\t" + it.second.genus + "\t" + it.second.species;
			tmpSpeciesFreqUMap[species]++;
		}
	}
	mOptions->transSearch.transReadsTaxonKOPairDeque.clear();

	if (mOptions->samples.size() > 0) {
		mOptions->transSearch.sampleSpeciesAbunUMap.insert(tmpSpeciesFreqUMap.begin(), tmpSpeciesFreqUMap.end());
	}

	std::vector<std::string> taxonRankVec{"Phylum", "Classm", "Order", "Family", "Genus", "Species"};
	//std::vector<std::pair< std::string, std::vector<std::pair<std::string, int> > > > rankVec;
	std::vector<std::vector<std::pair<std::string, int> > > rankVec;
	rankVec.reserve(taxonRankVec.size());

	auto sortedPhylumFreqVector = sortUMapToVector(tmpPhylumFreqUMap);
	rankVec.push_back(sortedPhylumFreqVector);
	tmpPhylumFreqUMap.clear();
	sortedPhylumFreqVector.clear();
	//rankVec.push_back(std::make_pair("Phylum", tmpSortedPhylumFreqVec));

	auto sortedClassmFreqVector = sortUMapToVector(tmpClassmFreqUMap);
	rankVec.push_back(sortedClassmFreqVector);
	tmpClassmFreqUMap.clear();
	sortedClassmFreqVector.clear();
	//rankVec.push_back(std::make_pair("Class", tmpSortedClassmFreqVec));

	auto sortedOrderFreqVector = sortUMapToVector(tmpOrderFreqUMap);
	rankVec.push_back(sortedOrderFreqVector);
	tmpOrderFreqUMap.clear();
	sortedOrderFreqVector.clear();
	//rankVec.push_back(std::make_pair("Order", tmpSortedOrderFreqVec));

	auto sortedFamilyFreqVector = sortUMapToVector(tmpFamilyFreqUMap);
	rankVec.push_back(sortedFamilyFreqVector);
	tmpFamilyFreqUMap.clear();
	sortedFamilyFreqVector.clear();
	//rankVec.push_back(std::make_pair("Family", tmpSortedFamilyFreqVec));

	auto sortedGenusFreqVector = sortUMapToVector(tmpGenusFreqUMap);
	rankVec.push_back(sortedGenusFreqVector);
	tmpGenusFreqUMap.clear();
	sortedGenusFreqVector.clear();

	sortedSpeciesFreqVector = sortUMapToVector(tmpSpeciesFreqUMap);
	rankVec.push_back(sortedSpeciesFreqVector);
	tmpSpeciesFreqUMap.clear();
	sortedSpeciesFreqVector.clear();
	//rankVec.push_back(std::make_pair("Species", tmpSortedSpeciesFreqVec));

	for (int i = 0; i < rankVec.size(); i++) {
		auto rank = taxonRankVec.at(i);
		fileoutname.clear();
		fileoutname = mOptions->mHomoSearchOptions.prefix + "_" + rank + "_taxon_abundance.txt";
		std::ofstream * fout = new std::ofstream();
		fout->open(fileoutname.c_str(), std::ofstream::out);
		if (!fout->is_open()) error_exit("Can not open abundance file: " + fileoutname);
		if (mOptions->verbose) loginfo("Starting to write taxon abundance table");
		std::string outputname;
		for (int j = 0; j <= i; j++) {
			outputname += taxonRankVec.at(j) + "\t";
		}
		*fout << outputname << "Reads_count\n";
		auto tempVec = rankVec.at(i);
		for (auto & itr : tempVec) {
			*fout << itr.first << "\t" << itr.second << "\n";
		}
		tempVec.clear();
		fout->flush();
		fout->close();
		if (fout) delete fout;
		if (mOptions->verbose) {
			auto msg = "Finish to write taxon abundance table of " + rank;
			loginfo(msg);
		}
	}
	rankVec.clear();
	taxonRankVec.clear();
}