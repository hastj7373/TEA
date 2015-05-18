#include <iostream>
#include <istream>
#include <string>
#include <map>
#include <list>
#include <utility>
#include <cstdlib>
#include <cmath>

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "utils/bamtools_utilities.h"

#include "Group.h"
#include "gzstream.h"

using namespace std;
using namespace BamTools;

typedef map< string, map<long, bool> > blacklist_t;

/* Utility function for read_isinfo. */
string
get_rg(string s)
{
	int i = s.find_last_of(" \r\n\t", s.find(":") + 1);
	return s.substr(i+1, s.size());
}

list<string>
read_rg_blacklist(string rg_blacklist_name)
{
	list<string> rg_blacklist;

	if (!rg_blacklist_name.empty()) {
		ifstream f(rg_blacklist_name.c_str());
		if (f.is_open()) {
			string line;
			while(!f.eof()) {
				getline(f, line);
				rg_blacklist.push_back(line);
			}
			f.close();
		} else {
			cout << "warning: cannot open reagroup blacklist '"
			     << rg_blacklist_name << "'" << endl;
		}
	}

	return rg_blacklist;
}

/*
 * Open the insert info file, read and parse it, close it,
 * and return a map of insert size statistics keyed by read
 * group name.
 *
 * FORMAT NOTES:
 *       <field>: <read group>
 *       <value>
 *   We are interested in <field> = "Mean insert size" or
 *   <field> = "Standard deviation of insert size".
 */
map< string, pair<double, double> >
read_isinfo(string file)
{
	ifstream isf(file.c_str());
	if (isf.fail()) {
		cout << "read_isinfo: failed to open insert info file '"
		     << file << "'" << endl;
		exit(1);
	}

	map< string, pair<double, double> > isinfo;
	string line, rg;
	while (true) {
		getline(isf, line);
		if (isf.eof())
			break;

		if (line.find("Median insert size") != string::npos) {
			rg = get_rg(line);
			getline(isf, line);
			double val = strtod(line.c_str(), NULL);
			map< string, pair<double,double> >::iterator iter =
				isinfo.find(rg);
			if (iter != isinfo.end())
				isinfo[rg].first = val;
			else
				isinfo[rg] = make_pair<double,double>(val, 0);
		} else if (line.find("Standard deviation of insert size") !=
		           string::npos) {
			rg = get_rg(line);
			getline(isf, line);
			double val = strtod(line.c_str(), NULL);
			map< string, pair<double,double> >::iterator iter =
				isinfo.find(rg);
			if (iter != isinfo.end())
				isinfo[rg].second = val;
			else
				isinfo[rg] = make_pair<double,double>(0, val);
		}
	}

	isf.close();

	return isinfo;
}

blacklist_t
read_blacklist(string name)
{
	istream *bfile;
	igzstream gzf;
	ifstream f;
	bool is_gz = false;
	blacklist_t blacklist;

	if (name.empty())
		return blacklist;

	if (name.substr(name.size() - 3, 3) == ".gz")
		is_gz = true;

	if (is_gz) {
		gzf.open(name.c_str());
		bfile = (istream *)&gzf;
	} else {
		f.open(name.c_str());
		bfile = (istream *)&f;
	}

	char buf[100];
	long pos, cov;

	while (bfile->good()) {
		*bfile >> buf;
		*bfile >> pos;
		*bfile >> cov;

		string chr(buf);
		blacklist[chr][pos] = true;
	}

	if (is_gz)
		((igzstream *)bfile)->close();
	else
		((ifstream *)bfile)->close();

	return blacklist;
}

/*
 * This function is responsible for determining which reads
 * from the input file are written to the output file.  If
 * this function returns true, the read is saved; else it is
 * discarded.
 */
long not_paired = 0;
long not_both_mapped = 0;
long not_same_chr = 0;
long not_diff_strands = 0;
long weird_insert_size = 0;
long big_insert_size = 0;
long same_pos = 0;
bool
is_discordant(BamAlignment &al, double median, double stdev, int nstdevs)
{
	/* Must be a paired end read */
	if (!al.IsPaired()) {
		++not_paired;
		return false;
	}

	/* Both reads must be mapped */
	if (!al.IsMapped() || !al.IsMateMapped()) {
		++not_both_mapped;
		return false;
	}

	if (al.Position == al.MatePosition && al.RefID == al.MateRefID) {
		++same_pos;
		return false;
	}

	/* If the mates map to different chroms, then discordant */
	if (al.RefID != al.MateRefID) {
		++not_same_chr;
		return true;
	}

	/* If both reads map to the same strand, then discordant */
	if (al.IsReverseStrand() == al.IsMateReverseStrand()) {
		++not_diff_strands;
		return true;
	}

	/* If the insert size is the wrong sign, then discordant */
	if ((al.IsReverseStrand() && al.InsertSize > 0) ||
	    (!al.IsReverseStrand() && al.InsertSize < 0)) {
		++weird_insert_size;
		return true;
	}

	/* If the insert size is too large, then discordant */
	if (abs(al.InsertSize) > ceil(median + nstdevs*stdev)) {
		++big_insert_size;
		return true;
	}

	return false;
}

/* vector refs */
bool
check_blacklist(blacklist_t &blacklist, RefVector &refs, long ref_id, long pos)
{
	blacklist_t::iterator iter = blacklist.find(refs[ref_id].RefName);
	if (iter != blacklist.end()) {
		map<long, bool>::iterator iter2 = iter->second.find(pos);
		if (iter2 != iter->second.end())
			return true;
	}

	return false;
}

bool
is_blacklisted(BamAlignment &al, RefVector &refs, blacklist_t &blist)
{
	/* If either mate is blacklisted, both mates are discarded */
	return check_blacklist(blist, refs, al.RefID, al.Position) ||
	       check_blacklist(blist, refs, al.MateRefID, al.MatePosition);
}

bool
is_rg_blacklisted(string &readgroup, list<string> &blist)
{
	for (list<string>::iterator it = blist.begin(); it != blist.end(); ++it)
		if (*it == readgroup)
			return true;
	return false;
}

/* 
 * This "none" string actually has to match the "none" readgroup
 * from bamreader or else the correct insert size statistics
 * won't be used.
 */
inline string f(string s) { return s.size() == 0 ? "none" : s; }

const string help_msg =
"Usage:       ./dre [options] <BAM input file>\n"
"\n"
"Options:\n"
"  -P         Remove PCR duplicates.  A set of reads are considered\n"
"             PCR duplicates if both mapped positions of the two mates\n"
"             are the same.  Currently, the first read in each duplicate\n"
"             set is kept and the rest are discarded--we could do better\n"
"             here by selecting the highest quality read from each set.\n"
"  -s integer Insert size cutoff for discordant read pairs that map to\n"
"             the same chromosome and have the correct orientation.\n"
"             The value is the number of standard deviations from the\n"
"             median insert size.  [3]\n"
"  -b file    A file of blacklisted genomic coordinates.  Reads mapping\n"
"             to any coordinate in this file will not be considered\n"
"             discordant, even if they are.\n"
"  -i file    Specify the file containing insert size distribution info.\n"
"  -o file    Output file.  Written in BAM format.\n"
"  -d file    PCR duplicate file.  Written in BAM format.\n"
"  -r file    Ignore reads from any read group in this file.  The file\n"
"             should specify one read group per line.\n"
"  -f         Force immediate run.  Do not prompt the user to verify\n"
"             run-time options.\n"
"  -v         Verbose mode.\n" "  -h         Print this help message."; 
int
main(int argc, char **argv)
{
	int option;
	bool no_prompt = false;
	bool verbose = false;
	int nstdevs = 3;
	string insfile_name, outfile_name, bfile_name, dupfile_name;
	string rg_blacklist_name;
	size_t max_chr_size = 300;
	bool filter_dups = false;
	while ((option = getopt(argc, argv, "vhfPs:b:i:d:o:r:")) != -1) {
		switch (option) {
		case 'r':
			rg_blacklist_name = optarg;
			break;
		case 'P':
			filter_dups = true;
			break;
		case 'm':
			max_chr_size = strtol(optarg, NULL, 10);
			break;
		case 'd':
			dupfile_name = string(optarg);
			break;
		case 'b':
			bfile_name = string(optarg);
			break;
		case 'i':
			insfile_name = string(optarg);
			break;
		case 'o':
			outfile_name = string(optarg);
			break;
		case 'v':
			verbose = true;
			break;
		case 's':
			nstdevs = strtol(optarg, NULL, 10);
			break;
		case 'f':
			no_prompt = true;
			break;
		case 'h':
		default:
			cout << help_msg << endl;
			exit(1);
		}
	}

	if (optind >= argc) {
		cout << help_msg << endl;
		exit(1);
	}

	string infile_name(argv[optind]);
	/* XXX: switching this to unordered_map might result in */
	/* significant speedup, depending on the # of RGs in the data */
	map< string, pair<double,double> > isinfo = read_isinfo(insfile_name);
	blacklist_t blacklist = read_blacklist(bfile_name);
	if (!filter_dups)
		dupfile_name = "";  /* discard whatever they might've input */

	list<string> rg_blacklist = read_rg_blacklist(rg_blacklist_name);

	/* Options blurb */
	cout << "Options:" << endl
	     << "   input BAM: " << f(infile_name) << endl
	     << "   insert file: " << f(insfile_name) << endl;
	map< string, pair<double,double> >::iterator iter = isinfo.begin();
	while (iter != isinfo.end()) {
		cout << "      " << iter->first << ": median "
		     << iter->second.first << ", stdev "
		     << iter->second.second << endl;
		++iter;
	}
	cout << "   blacklist file: " << f(bfile_name) << endl;
	if (blacklist.size() > 0) {
		int npos = 0;
		blacklist_t::iterator iter = blacklist.begin();
		while (iter != blacklist.end()) {
			npos += iter->second.size();
			++iter;
		}

		cout << "      " << blacklist.size()
		     << " chromosomes present in blacklist" << endl
		     << "      " << npos << " blacklisted positions" << endl;
	}
	cout << "   ignoring read groups: ";
	for (list<string>::iterator it = rg_blacklist.begin();
	     it != rg_blacklist.end(); ++it)
		cout << *it << " ";
	cout << endl;

	cout << "   output file: " << f(outfile_name) << endl
	     << "   PCR dup filter: " << boolalpha << filter_dups << endl
	     << "   PCR duplicate file: " << f(dupfile_name) << endl
	     << "   n stdevs: " << nstdevs << endl;

	/* User confirmation before run unless -f was specified */
	if (!no_prompt) {
		cout << "Continue? (y/N) ";
		char c;
		cin >> noskipws >> c;
		if (tolower(c) != 'y') {
			cout << "Terminating..." << endl;
			exit(1);
		}
	} else
		cout << "Forcing immediate execution." << endl;

	/* Input BAM */
	BamReader reader;
	if (!reader.Open(infile_name)) {
		cout << "ERROR: could not open input BAM file '"
		     << infile_name << "' for reading" << endl;
		exit(1);
	}

	/* Output BAM */
	BamWriter writer;
	if (!writer.Open(outfile_name,
	                 reader.GetHeaderText(),
	                 reader.GetReferenceData())) {
		cout << "ERROR: could not open output BAM file '"
		     << outfile_name << "' for writing" << endl;
		exit(1);
	}

	BamWriter *dups = NULL;
	if (dupfile_name.size() > 0)  {
		dups = new BamWriter;
		if (!dups->Open(dupfile_name,
		    reader.GetHeaderText(),
		    reader.GetReferenceData())) {
			cout << "ERROR: could not open dup file '"
			     << dupfile_name << "' for writing" << endl;
			exit(1);
		}
	}

	BamAlignment al;
	RefVector refs = reader.GetReferenceData();
	long num_in = 0;
	long num_out = 0;
	long num_blacklisted = 0;
	Group g;
	record prev_dups;

	/* Use the slow method to read tags so we can get the read group */
	if (verbose)
		cout << "Starting file scan..." << endl;
	while (reader.GetNextAlignment(al)) {
		++num_in;
		string rg;
		al.GetReadGroup(rg);
		rg = f(rg);
		pair<double, double> p = isinfo[rg];

		if (is_rg_blacklisted(rg, rg_blacklist))
			continue;

		if (!is_discordant(al, p.first, p.second, nstdevs))
			continue;

		if (is_blacklisted(al, refs, blacklist)) {
			++num_blacklisted;
			continue;
		}


		/* Handle PCR dup filtering if specified, else write */
		if (filter_dups) {
			if (g.should_write(al)) {
				num_out += g.write(writer, prev_dups, dups);
				g.clear();
			}
			g.add(al);
		} else {
			writer.SaveAlignment(al);
			++num_out;
		}
	}
	if (filter_dups)
		num_out += g.write(writer, prev_dups, dups);

	if (verbose)
		cout << "Finished file scan." << endl;

	cout << infile_name << ": " << num_in << " reads read" << endl;
	if (!bfile_name.empty()) {
		cout << bfile_name << ": rejected " << num_blacklisted
		     << " discordant pairs" << endl;
	}
	cout << outfile_name << ": " << num_out << " reads written" << endl;
	cout << "   not paired:             " << not_paired << endl
	     << "   not both mapped:        " << not_both_mapped << endl
	     << "   not same chrom:         " << not_same_chr << endl
	     << "   not diff strand:        " << not_diff_strands << endl
	     << "   weird insert size:      " << weird_insert_size << endl
	     << "   big insert size:        " << big_insert_size << endl
	     << "   mates map to same posn: " << same_pos << endl;
	cout << "discarded " <<
		not_same_chr + not_diff_strands + weird_insert_size + big_insert_size - num_out - num_blacklisted - same_pos
	     << " PCR duplicates" << endl;
	reader.Close();
	writer.Close();
	if (dups != NULL)
		dups->Close();

	return 0;
}
