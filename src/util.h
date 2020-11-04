#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <algorithm>
#include <time.h>
#include <mutex>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <stdio.h>  /* defines FILENAME_MAX */
#include <tuple>
#include <utility>
#include <numeric>
#include <vector>
#include <cmath>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#include <limits.h>
#define GetCurrentDir getcwd
#endif

#ifdef WINDOWS
#include <windows.h>
#define GetExePath _getcwd
#else
#include <unistd.h>
#include <limits.h>
#define GetExePath getcwd
#endif

#include "options.h"


using namespace std;

inline char complement(char base) {
    switch(base){
        case 'A':
        case 'a':
            return 'T';
        case 'T':
        case 't':
            return 'A';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        default:
            return 'N';
    }
}

inline bool starts_with( string const & value,  string const & starting)
{
    if (starting.size() > value.size()) return false;
    return  equal(starting.begin(), starting.end(), value.begin());
}

inline bool ends_with( string const & value,  string const & ending)
{
	if (ending.size() > value.size()) return false;
	return  equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline string trim(const string& str)
{
    string::size_type pos = str.find_first_not_of(' ');
    if (pos == string::npos)
    {
        return string("");
    }
    string::size_type pos2 = str.find_last_not_of(' ');
    if (pos2 != string::npos)
    {
        return str.substr(pos, pos2 - pos + 1);
    }
    return str.substr(pos);
}

inline string trimStr(const string& str)
{
    string::size_type pos = str.find_first_not_of(' ');
    if (pos == string::npos)
    {
        return string("");
    }
    string::size_type pos2 = str.find_last_not_of(' ');
    if (pos2 != string::npos)
    {
        return str.substr(pos, pos2 - pos + 1);
    }
    return str.substr(pos);
}

inline int split(const string& str, vector<string>& ret_, string sep = ",")
{
    if (str.empty())
    {
        return 0;
    }

    string tmp;
    string::size_type pos_begin = str.find_first_not_of(sep);
    string::size_type comma_pos = 0;

    while (pos_begin != string::npos)
    {
        comma_pos = str.find(sep, pos_begin);
        if (comma_pos != string::npos)
        {
            tmp = str.substr(pos_begin, comma_pos - pos_begin);
            pos_begin = comma_pos + sep.length();
        }
        else
        {
            tmp = str.substr(pos_begin);
            pos_begin = comma_pos;
        }

        ret_.push_back(tmp);
        tmp.clear();
    }
    return 0;
}

inline int splitStr(const string& str, vector<string>& ret_, string sep = ",")
{
    if (str.empty())
    {
        return 0;
    }

    string tmp;
    string::size_type pos_begin = str.find_first_not_of(sep);
    string::size_type comma_pos = 0;

    while (pos_begin != string::npos)
    {
        comma_pos = str.find(sep, pos_begin);
        if (comma_pos != string::npos)
        {
            tmp = str.substr(pos_begin, comma_pos - pos_begin);
            pos_begin = comma_pos + sep.length();
        }
        else
        {
            tmp = str.substr(pos_begin);
            pos_begin = comma_pos;
        }

        ret_.push_back(tmp);
        tmp.clear();
    }
    return 0;
}

inline string replace(const string& str, const string& src, const string& dest)
{
    string ret;

    string::size_type pos_begin = 0;
    string::size_type pos       = str.find(src);
    while (pos != string::npos)
    {
        ret.append(str.data() + pos_begin, pos - pos_begin);
        ret += dest;
        pos_begin = pos + 1;
        pos       = str.find(src, pos_begin);
    }
    if (pos_begin < str.length())
    {
        ret.append(str.begin() + pos_begin, str.end());
    }
    return ret;
}

inline string reverse(const string& str) {
    string ret(str.length(), 0);
    for(int pos=0; pos<str.length(); pos++) {
        ret[pos] = str[str.length() - pos - 1];
    }
    return ret;
}

inline string basename(const string& filename){
    string::size_type pos = filename.find_last_of('/');
    if (pos == string::npos)
        return filename;
    else if(pos == filename.length()-1)
        return ""; // a bad filename
    else
        return filename.substr(pos+1, filename.length() - pos - 1);
}

inline string get_current_dir() {
    char buff[FILENAME_MAX]; //create string buffer to hold path
    GetCurrentDir(buff, FILENAME_MAX);
    string current_working_dir(buff);
    return current_working_dir;
}

inline string get_upper_dir() {
    std::string cwd = get_current_dir();
    std::string::size_type bepos = cwd.find_last_of("/");
    cwd.erase(bepos);
    return(cwd);
}

#ifdef WINDOWS
std::string getexepath()
{
  char result[ MAX_PATH ];
  return std::string( result, GetModuleFileName( NULL, result, MAX_PATH ) );
}
#else
inline string GetExePath() {
    char result[ PATH_MAX ];
    ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
    return std::string(result, (count > 0) ? count : 0);
}
inline string get_seq2fun_dir() {
    std::string cwd = GetExePath();
    std::string::size_type bepos = cwd.find("/bin");
    cwd.erase(bepos);
    return(cwd);
}
#endif


inline string dirname(const string& filename){
    string::size_type pos = filename.find_last_of('/');
    if (pos == string::npos) {
        return "./";
    } else {
        return filename.substr(0, pos+1);
    }
}

inline string joinpath(const string& dirname, const string& basename){
    if(dirname[dirname.length()-1] == '/'){
        return dirname + basename;
    } else {
        return dirname + "/" + basename;
    }
}

//Check if a string is a file or directory
inline bool file_exists(const  string& s)
{
    bool exists = false;
    if(s.length() > 0) {
        struct stat status;
        int result = stat( s.c_str(), &status );
        if(result == 0) {
            exists = true;
        }
    }
    return exists;
}


// check if a string is a directory
inline bool is_directory(const  string& path)
{
    bool isdir = false;
    struct stat status;
    // visual studion use _S_IFDIR instead of S_IFDIR
    // http://msdn.microsoft.com/en-us/library/14h5k7ff.aspx
#ifdef _MSC_VER
#define S_IFDIR _S_IFDIR
#endif
    stat( path.c_str(), &status );
    if ( status.st_mode &  S_IFDIR  ) {
        isdir = true;
    }
// #endif
    return isdir;
}

inline void check_file_valid(const  string& s) {
    if(!file_exists(s)){
        cerr << "ERROR: file '" << s << "' doesn't exist, quit now" << endl;
        exit(-1);
    }
    if(is_directory(s)){
        cerr << "ERROR: '" << s << "' is a folder, not a file, quit now" << endl;
        exit(-1);
    }
}

inline void make_dir_con(const string & s){
    string::size_type pos = s.find_last_of('/');
    if (pos != string::npos) {
        auto dirName = s.substr(0, pos);
        std::cout << "dir name is " << dirName << std::endl;
    if(!is_directory(dirName)){
      if(mkdir(dirName.c_str(), 0777) == -1){
        cerr << "ERROR: can not create output dir: " << s << endl;
        exit(-1);
      }
     }
    }
}

inline void check_file_writable(const  string& s) {
    string dir = dirname(s);
    if(!file_exists(dir)) {
        cerr << "ERROR: '" << dir << " doesn't exist. Create this folder and run this command again." << endl;
        exit(-1);
    }
    if(is_directory(s)){
        cerr << "ERROR: '" << s << "' is not a writable file, quit now" << endl;
        exit(-1);
    }
}

// Remove non alphabetic characters from a string
inline  string str_keep_alpha(const  string& s)
{
     string new_str;
    for( size_t it =0; it < s.size(); it++) {
        if(  isalpha(s[it]) ) {
            new_str += s[it];
        }
    }
    return new_str;
}


// Remove invalid sequence characters from a string
inline void str_keep_valid_sequence(  string& s, bool forceUpperCase = false)
{
    size_t total = 0;
    const char case_gap = 'a' - 'A';
    for( size_t it =0; it < s.size(); it++) {
        char c = s[it];
        if(forceUpperCase && c>='a' && c<='z') {
            c -= case_gap;
        }
        if(  isalpha(c) || c == '-' || c == '*' ) {
            s[total] = c;
            total ++;
        }
    }

    s.resize(total);
}

//inline int convert_float2int(int & a, float & b){
//    return((int)round(a * b));
//}

inline int find_with_right_pos(const string& str, const string& pattern, int start=0) {
    int pos = str.find(pattern, start);
    if (pos < 0)
        return -1;
    else
        return pos + pattern.length();
}

inline void str2upper(string& s){
    transform(s.begin(), s.end(), s.begin(), (int (*)(int))toupper);
}

inline void str2lower(string& s){
    transform(s.begin(), s.end(), s.begin(), (int (*)(int))tolower);
}

inline char num2qual(int num) {
    if(num > 127 - 33)
        num = 127 - 33;
    if(num < 0)
        num = 0;

    char c = num + 33;
    return c;
}

inline void error_exit(const string& msg) {
    cerr << "ERROR: " << msg << endl;
    exit(-1);
}

extern mutex logmtx;
inline void loginfo(const string s){
    logmtx.lock();
    time_t tt = time(NULL);
    tm* t= localtime(&tt);
    cerr<<"["<<t->tm_hour<<":"<<t->tm_min<<":"<<t->tm_sec<<"] "<<s<<endl;
    logmtx.unlock();
}


inline string trimName(string& str){
   // string strnew;
    str.erase(str.begin());
    string suffixStartCharacters = " /\t\r";
    size_t n = str.find_first_of(suffixStartCharacters);
    if(n != string::npos){
        return str.erase(n);
    } else {
        return str;
    }
}


inline std::string getMostFreqStrFromVec(std::vector<std::string> & vectorko){
   std::map<std::string, int> freq;
   for(int i = 0; i < vectorko.size(); i ++){
       freq[vectorko[i]]++;
   }
   int max_F = 0;
   std::string res = "";
   for(auto & j : freq){
       if(max_F < j.second){
           res = j.first;
           max_F = j.second;
       }
   }
   freq.clear();
   vectorko.clear();
   return res;
}

inline std::string getUniqStrFromVec(std::vector<std::string> & vectorko) {
    sort(vectorko.begin(), vectorko.end());
    vectorko.erase( unique( vectorko.begin(), vectorko.end() ), vectorko.end() );
    
    std::string uniq = "";
    if (vectorko.size() == 1) {
        uniq = vectorko[1];
    }
    vectorko.clear();
    return uniq;
}

inline string removeStr(const string &str, const string &src)
{
    string ret;
    string::size_type pos_begin = 0;
    string::size_type pos = str.find(src);
    while (pos != string::npos)
    {
        ret.append(str.data() + pos_begin, pos - pos_begin);
        pos_begin = pos + src.length();
        pos = str.find(src, pos_begin);
    }
    if (pos_begin < str.length())
    {
        ret.append(str.begin() + pos_begin, str.end());
    }
    return ret;
}

template<typename T1, typename T2>
inline double getPercentage(T1 v1, T2 v2){
   double ret = double (v1 * 100) / (double) v2;
   return ret;
}
//inline double getPercentage(int v1, long v2){
   //double ret = double (v1 * 100) / (double) v2;
   //return ret;
//}

//inline double getPercentageInt(int v1, int v2){
   //double ret = double (v1 * 100) / (double) v2;
   //return ret;
//}

inline std::vector<std::pair<std::string, int> > sortUMapToVector(std::unordered_map<std::string, int> &unMap)
{
    int size = unMap.size();
    std::vector<std::pair<std::string, int> > sortedVec(size);
    std::partial_sort_copy(unMap.begin(),
                           unMap.end(),
                           sortedVec.begin(),
                           sortedVec.end(),
                           [](std::pair<const std::string, float> const &l,
                              std::pair<const std::string, float> const &r) {
                               return l.second > r.second;
                           });
    return sortedVec;
}

inline std::vector<std::tuple<std::string, double, int, int> > sortTupleVector(std::vector<std::tuple<std::string, double, int, int> > & OriVec)
{
    int size = OriVec.size();
    std::vector<std::tuple<std::string, double, int, int> > sortedVec(size);
    std::partial_sort_copy(OriVec.begin(),
                           OriVec.end(),
                           sortedVec.begin(),
                           sortedVec.end(),
                           [](const std::tuple<std::string, double, int, int > &l,
                              const std::tuple<std::string, double, int, int > &r) {
                               return get<1>(l) > get<1>(r);
                           });
    return sortedVec;
}

inline void getUniqVec(std::vector<std::string> &orgVec)
{
    std::unordered_set<std::string> s;
    std::vector<std::string>::iterator itr = orgVec.begin();
    for (auto curr = orgVec.begin(); curr != orgVec.end(); ++curr)
    {
        if (s.insert(*curr).second)
            *itr++ = *curr;
    }
    orgVec.erase(itr, orgVec.end());
}

inline void stripChar(std::string & s) {
    for (auto it = s.begin(); it != s.end(); ++it) {
        if (!isalpha(*it)) {
            s.erase(it);
            it--;
        }
    }
}

inline void writeSampleSpeciesTable(std::vector< std::pair<std::string, std::unordered_map<std::string, int> > > & sampleSpeciesTableVec, Options & opt){

std::set<std::string> allSpeciesIdsSet;//unique kos
    for(auto & it : sampleSpeciesTableVec){
	std::cout << it.first << std::endl; 
        auto itr = it.second;
        for(auto & itt : itr){
            allSpeciesIdsSet.insert(itt.first);
        }
    }
    auto newDir = dirname(opt.mHomoSearchOptions.prefix);
    std::string fileoutname = newDir + "/All_sample_Species_abundance_table.txt";
    std::ofstream * fout = new std::ofstream();
    fout->open(fileoutname.c_str(), std::ofstream::out);
    if(!fout->is_open()) error_exit("Can not open all_all_sample_species_abundance_table.txt");
    if (opt.verbose) loginfo("Starting to write all samples species abundance table");

for(auto & it : sampleSpeciesTableVec){
}

    *fout << "Phylum" << "\t" << "Classm" << "\t" << "Order" << "\t" << "Family" << "\t" << "Genus" << "\t" << "Species" << "\t";
    for (auto & it : sampleSpeciesTableVec) {
        *fout << it.first << "\t";
    }
    *fout << "\n";

	for (auto & specNm : allSpeciesIdsSet){
		*fout << specNm << "\t";
		for(auto & smp : sampleSpeciesTableVec){
			auto count = smp.second.find(specNm);
			if(count != smp.second.end()){
					*fout << count->second << "\t";
			} else {
				*fout << 0 << "\t";
			}
		}
		*fout << "\n";
    }
	allSpeciesIdsSet.clear();
	fout->flush();
    fout->close();
    if (fout) delete fout;
}

inline std::map<int, int> getRareFaction(std::unordered_map<std::string, int> & preRareUMap, long totalMappedReads, long totalRawReads){
	std::vector<std::string> reshuff_vec;
	reshuff_vec.reserve(totalMappedReads);
	for(auto & it : preRareUMap){
	for (int i = 0; i < it.second; i++){
			reshuff_vec.push_back(it.first);
		}
	}
	std::map<int, int> rarefaction_map_tmp;
	std::random_shuffle(reshuff_vec.begin(), reshuff_vec.end());
	int total = reshuff_vec.size();
	double ratio = totalRawReads / totalMappedReads;
	int step = 100;
	int step_size = floor(total / step);
	auto first = reshuff_vec.begin();
	rarefaction_map_tmp[0] = 0;
	for(int i = 1; i < step; i++){
		auto last = reshuff_vec.begin() + step_size * i;
		std::vector<std::string> rarefaction_vec(first, last);
		std::sort(rarefaction_vec.begin(), rarefaction_vec.end());
		int unic = std::unique(rarefaction_vec.begin(), rarefaction_vec.end()) - rarefaction_vec.begin();
		rarefaction_map_tmp[(int) round((step_size * i) * ratio)] = unic;
		rarefaction_vec.clear();
	}
	std::sort(reshuff_vec.begin(), reshuff_vec.end());
	int unic = std::unique(reshuff_vec.begin(), reshuff_vec.end()) - reshuff_vec.begin();
	rarefaction_map_tmp[totalRawReads] = unic;
	reshuff_vec.clear();
	return(rarefaction_map_tmp);
	rarefaction_map_tmp.clear();
}

#endif /* UTIL_H */
