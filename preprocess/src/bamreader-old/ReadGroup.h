#include <string>
#include <vector>
#include <list>
#include "api/BamAlignment.h"
#include "gzstream.h"

namespace std {

class ReadGroup {
public:
	ReadGroup();
	ReadGroup(BamTools::BamAlignment &al, int max_isize, int isize_samples,
		string prefix, list<string> blacklist);
	~ReadGroup();

	void witness(BamTools::BamAlignment &al);

	/* Also trims the read based on q and min_read_len.  If the read
	 * becomes too short, it is rejected and record() returns false.
	 */
	bool recordSC(BamTools::BamAlignment &a, BamTools::BamAlignment &b,
		int clippedMate, ogzstream &clipfile, int big_s_bps, int frag_size,
		int n_cutoff);

	bool recordUM(BamTools::BamAlignment &a, BamTools::BamAlignment &b,
		int big_s_bps, int n_cutoff);

	bool recordUU(BamTools::BamAlignment &a, BamTools::BamAlignment &b,
		int big_s_bps, int n_cutoff);

	string getName();
	int getNReads();
	vector<int> getReadlens();
	bool is_blacklisted();

private:
	ogzstream f1;
	ogzstream f2;
	string name;
	int nreads;
	vector<int> inserts;
	vector<int> readlens;

	string prefix;
	int max_isize;
	int isize_samples;
	bool extract_scs_and_mates;
	bool blacklisted;
};

ostream & operator<< (ostream &os, ReadGroup &rg);

} // namespace std

BamTools::BamAlignment snip(BamTools::BamAlignment &al, int start, int len);
