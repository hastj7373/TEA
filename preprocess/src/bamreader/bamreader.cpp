/*
 * The input BAM must be sorted for only one piece of functionality:
 * the blacklist.  If no blacklist is desired, then the BAM can be
 * unsorted.
 */

#include <ostream>
#include <sstream>
#include <cstdlib>
#include <climits>
#include <cmath>
#include <cstring>
#include <map>
#include <algorithm>
#include <vector>
#include <tr1/unordered_map>

#include <strings.h>

#include <errno.h>
#include <sys/stat.h>

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "utils/bamtools_utilities.h"
#include "gzstream.h"
#include "ReadGroup.h"
#include "Histogram.h"

using namespace std;
using namespace BamTools;


void
print_progress(long cur, long total, time_t start)
{
	static int last_pct = -1;
	int cur_pct = (long)100 * cur/total;  /* Floor of current pct */

	if (last_pct != cur_pct) {
		char bar[100];
		int i;
		time_t d = time(NULL) - start;

		for (i = 0; i < cur_pct; ++i)
			bar[i/5] = '=';
		bar[i/5] = '\0';
		printf("   %3d%% [%-20s] %ldh%ldm%lds\r",
			cur_pct, bar, d / 3600, (d % 3600) / 60, d % 60);
		fflush(stdout);
		if (cur_pct == 100)
			printf("\n");
		last_pct = cur_pct;
	}
}

/*
 * Witness the quality scores for alignment 'al'.  This
 * function tracks some information about the distribution
 * of quality scores in the whole file to try to guess the
 * scoring function.
 */
/* need a 64 bit integer here since there can easily be
 * billions of reads in a large bam, and hence hundreds
 * of billions of quality scores. */
#define MAX_QVAL 255    /* there are 255 possible char values in ASCII */
uint64_t sample_space[MAX_QVAL];
void
witness_qualities(BamAlignment &al)
{
	for (int i = 0; i < al.Length; ++i)
		++sample_space[al.Qualities[i]];
}

typedef struct qencoding {
	const char *name;
	int min;
	int max;
	int offset;
} QENCODING;

QENCODING known_qencodings[] = {
	/* ASCII-33 */
	{ "Sanger",			33,	126 },

	/* All of the Illumina formats have a 64 offset */
	{ "Solexa/Illumina 1.0",	59,	126 },

	/* Illumina 1.3 no longer allows negative values */
	{ "Illumina 1.3",		64,	126 },

	/* The only difference between 1.3 and 1.5 is that 1.5 does
	 * not utilize values 0 or 1, and value 2 is used as a read
	 * quality indicator. */
	{ "Illumina 1.5",		66,	126 },
};

int qenc = -1;
void
report_qualities(bool verbose)
{
	int a, b;
	for (int i = 0; i < MAX_QVAL; ++i) {
		if (sample_space[i] != 0) {
			a = i;
			break;
		}
	}
	for (int i = MAX_QVAL - 1; i >= 0; --i) {
		if (sample_space[i] != 0) {
			b = i;
			break;
		}
	}
	if (verbose) {
		cout << "Quality score distribution:" << endl;
		for (int i = 0; i < MAX_QVAL; ++i)
			cout << i << "\t" << sample_space[i] << endl;
	}

	cout << "   Quality string info:" << endl;
	cout << "      quality score range: [" << a << ", " << b << "]" << endl;

	/* Simple guess at the encoding type */
	int best = 0, best_diff = INT_MAX;
	for (int i = 0; i < sizeof(known_qencodings)/sizeof(QENCODING); ++i) {
		int d = abs(known_qencodings[i].min - a) +
			abs(known_qencodings[i].max - b);
		if (d < best_diff) {
			best = i;
			best_diff = d;
		}
	}
	qenc = best;
	for (int i = 0; i < sizeof(known_qencodings)/sizeof(QENCODING); ++i) {
		if (qenc == i)
			cout << "    * ";
		else
			cout << "      ";

		cout.width(25);
		cout << left << known_qencodings[i].name
		     << "[" << known_qencodings[i].min << ", "
		     << known_qencodings[i].max << "]" << endl;
	}
}

void
trim_read(BamAlignment &al, int q)
{
	int i, best_i, s, best_s;

	int start_idx = al.IsReverseStrand() ? 0 : al.Length - 1;
	int increment = al.IsReverseStrand() ? 1 : -1;

	/* substr always throws errors on 0 length strings */
	if ((al.Qualities[start_idx] - known_qencodings[qenc].min) >= q ||
	    al.Length <= 0)
		return;

	s = 0;
	best_s = INT_MIN;
	best_i = start_idx;
	for (i = start_idx; i >= 0 && i < al.Length; i += increment) {
		s += q - (al.Qualities[i] - known_qencodings[qenc].min);
		if (s > best_s) {
			best_s = s;
			best_i = i;
		}
	}

	int x = al.IsReverseStrand() ? best_i + 1 : al.Length - best_i;

	/* if x >= al.Length, then the new length is 0 */
	int new_start = al.IsReverseStrand() ? min(x, al.Length-1) : 0;
	int new_len = al.Length - x;

	al.Length = new_len;
	try {
		al.QueryBases = al.QueryBases.substr(new_start, new_len);
		al.Qualities = al.Qualities.substr(new_start, new_len);
	} catch (exception &e) {
		cout << "ERROR: substr failed in trim_read" << endl;
		cout << "trim_read(" << al.Name << ", " << q << ")" << endl;
		cout << new_start << " " << new_len << " x: " << x
		     << " best_s: " << best_s << " best_i: " << best_i << endl;
		cout << al.IsReverseStrand() << endl;
		cout << al.QueryBases << endl;
		cout << al.Qualities << endl;
		exit(1);
	}

	/* The trim may extend past one CIGAR op, but we only *
	 * care about S ops in the first or last position.    */
	if (al.CigarData.size() == 0)
		return;

	int idx = al.IsReverseStrand() ? 0 : al.CigarData.size() - 1;
	al.CigarData[idx].Length = max(0, (signed)al.CigarData[idx].Length - x);
}

bool
isBigS(BamAlignment &al, CigarOp &cop, int n)
{
	/* cop *can* be a NULL ref if the alignment has a "*" CIGAR string */
	return (al.CigarData.size() > 0) &&
		(cop.Type == 'S' && min(al.Length, (signed)cop.Length) >= n);
}

bool
isMatePair(BamAlignment &a, BamAlignment &b)
{
	return (a.IsFirstMate() && b.IsSecondMate()) ||
		(a.IsSecondMate() && b.IsFirstMate());
}

int
getMateNumber(BamAlignment &al)
{
	return al.IsFirstMate() ? 1 : 2;
}

/*
 * Write out the substring [start, start+length] (inclusive) of the
 * aligned read 'al' to the file 'f' in FASTQ format.
 *
 * writeFQ() quietly aborts if al contains more than n_cutoff N bases.
 *
 * NOTE: Quite similar to the PrintFastq function in the bamtools
 * utils.
 */
void
writeFQ(ofstream &f, BamAlignment &al)
{
	string quals = al.Qualities;
	string bases = al.QueryBases;
	if (al.IsReverseStrand()) {
		Utilities::Reverse(quals);
		Utilities::ReverseComplement(bases);
	}

	/* Using \n instead of endl to avoid unnecessary stream flushes */
	f << "@" << al.Name << "\n"
	  << bases << "\n"
	  << "+" << "\n"
	  << quals << "\n";
}

void
writeFQpair(ofstream &f1, BamAlignment &a,
            ofstream &f2, BamAlignment &b, int n_cutoff)
{
	int n = 0;
	for (size_t i = 0; i < a.QueryBases.size(); ++i)
		if (a.QueryBases[i] == 'N')
			++n;
	if (n > n_cutoff)
		return;
	for (size_t i = 0; i < b.QueryBases.size(); ++i)
		if (b.QueryBases[i] == 'N')
			++n;
	if (n > n_cutoff)
		return;

	writeFQ(f1, a);
	writeFQ(f2, b);
}

ogzstream split1;
ogzstream split2;
void
split_read(BamAlignment &al, int frag_size, int n_cutoff)
{
	if (al.Length < 2*frag_size) {
		cout << "ERROR: split_read tried to save alignment "
		     << al.Name << " with length " << al.Length
		     << " when -s=" << frag_size << endl;
		exit(1);
	}

	BamAlignment left = snip(al, 0, frag_size);
	ostringstream os;
	os << al.Length;
	left.Name += "_" + os.str();
	BamAlignment right = snip(al, al.Length - frag_size, frag_size);
	right.Name += "_" + os.str();
	writeFQpair((ofstream &)split1, left,
	            (ofstream &)split2, right, n_cutoff);
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

const string help_msg =
"Usage:       ./bamreader [options] <BAM input file>\n"
"\n"
"Options:\n"
"  -D         Generate an additional BAM file named \"mapped_sc_um.bam\"\n"
"             which contains read pairs from non-blacklisted readgroups\n"
"             that are either mapped+softclipped or mapped+unmapped.\n"
"             Reads in this file will be quality-trimmed.\n"
"  -b cutoff  Generate a blacklist of genomic positions with coverage\n"
"             greater than 'cutoff'.  The input BAM must be sorted if\n"
"             -b is specified with cutoff > 0.  Set the cutoff to <= 0\n"
"             to disable blacklist construction.  [0]\n"
"             WARNING: currently, the code does not take CIGAR strings\n"
"             into account when computing coverage.  The 'calDepth'\n"
"             example in samtools does consider the CIGAR string.\n"
"  -u         Do not process UU pairs.\n"
"  -x prefix  All output files are prepended with 'prefix'.  [output]\n"
"  -g isize   Maximum allowed insert size for computing the insert\n"
"             size distribution.  This is useful to filter out large\n"
"             inserts that result from SVs.  [1,000]\n"
"  -c bps     Defines the cutoff for calling soft clips.  A read is\n"
"             soft clipped if it contains an 'S' CIGAR operation of\n"
"             at least 'bps' basepairs.  [20]\n"
"  -s bps     Fragment lengths for the unmapped and softclipped files.\n"
"             Defaults to the size specified by '-c'.\n"
"  -n nreads  Output the insert sizes for the first 'nreads' reads\n"
"             from each read group to determine the insert size\n"
"             distribution.  If nreads < 0, all reads in each group\n"
"             will be used.  (Very large values can cause significant\n"
"             loss of speed!)  [10,000]\n"
"  -q qv      Read trimming parameter.  Equivalent to BWA's -q option.  [0]\n"
"  -N integer Directly before writing a read to any output file (except\n"
"             .unmapped.fq and .softclips.fq), require that the read\n"
"             have < 'integer' Ns in its sequence.  [5]\n"
"  -r file    Ignore reads from any read group in this file.  The file\n"
"             should specify one read group per line.\n"
"  -f         Force immediate run.  Do not prompt the user to verify\n"
"             run-time options.\n"
"  -v         Verbose mode.\n"
"  -h         Print this help message.";

int
main(int argc, char **argv)
{
	int option;
	bool no_prompt = false;
	string prefix = "output";
	string rg_blacklist_fname;
	int max_isize = 1000;
	int big_s_bps = 20;
	int frag_size = -1;
	int isize_samples = 10000;
	int q = 0;
	int min_read_len = 0;
	bool verbose = false;
	bool processUU = true;
	int coverage_cutoff = 0;
	int n_cutoff = 5;
	bool generate_mapped_sc_um_file = false;
	while ((option = getopt(argc, argv, "Dvhufx:g:c:n:q:m:s:b:N:r:")) != -1) {
		switch (option) {
		case 'D':
			generate_mapped_sc_um_file = true;
			break;
		case 'r':
			rg_blacklist_fname = optarg;
			break;
		case 'N':
			n_cutoff = strtol(optarg, NULL, 10);
			break;
		case 'b':
			coverage_cutoff = strtol(optarg, NULL, 10);
			break;
		case 'u':
			processUU = false;
			break;
		case 'v':
			verbose = true;
			break;
		case 'x':
			prefix = optarg;
			break;
		case 'g':
			max_isize = strtol(optarg, NULL, 10);
			break;
		case 'c':
			big_s_bps = strtol(optarg, NULL, 10);
			break;
		case 'n':
			isize_samples = strtol(optarg, NULL, 10);
			break;
		case 'f':
			no_prompt = true;
			break;
		case 'q':
			q = strtol(optarg, NULL, 10);
			break;
		case 's':
			frag_size = strtol(optarg, NULL, 10);
			break;
		case 'h':
		default:
			cout << help_msg << endl;
			exit(1);
		}
	}

	/* No longer an option */
	min_read_len = 2 * big_s_bps;
	if (frag_size == -1)
		frag_size = big_s_bps;  /* -s defaults to -c */

	if (optind >= argc) {
		cout << help_msg << endl;
		exit(1);
	}

	string fname(argv[optind]);
	string sep = prefix != "" ? "." : "";
	const string umfname = prefix + sep + "unmapped.fq.gz";
	const string clipname = prefix + sep + "softclips.fq.gz";
	const string isinfoname = prefix + sep + "isinfo";
	const string umrdistname = prefix + sep + "unmapped.rdist";
	const string scrdistname = prefix + sep + "softclips.rdist";
	const string blistname = prefix + sep + "blacklist.gz";
	const string split1name = prefix + sep + "sr.1.fq.gz";
	const string split2name = prefix + sep + "sr.2.fq.gz";
	const string mapped_um_name = prefix + sep + "mapped_um.bam";
	const string mapped_sc_name = prefix + sep + "mapped_sc.bam";

	list<string> rg_blacklist = read_rg_blacklist(rg_blacklist_fname);

	cout << "Options:" << endl
	     << "   input BAM: " << fname << endl
	     << "   prefix: " << prefix << endl;
	if (coverage_cutoff > 0)
		cout << "   generating blacklist: coverage >= "
		     << coverage_cutoff << " reads" << endl;
	else
		cout << "   no blacklist will be generated" << endl;
	cout << "   max insert size: " << max_isize << " bases" << endl
	     << "   soft clip threshold: " << big_s_bps << " bases" << endl
	     << "   fragment size: " << frag_size << " bases" << endl
	     << "   insert size samples: ";
	if (isize_samples < 0)
		cout << "all";
	else
		cout << isize_samples;
	cout << " reads" << endl;
	cout << "   blacklisted read groups: ";
	for (list<string>::iterator it = rg_blacklist.begin();
	     it != rg_blacklist.end(); ++it)
		cout << *it << " ";
	cout << endl;
	cout << "   trim reads: " << q << endl
	     << "   max N bps per read: " << n_cutoff << endl;
	cout << "   output files: " << endl
	     << "      " << umfname << endl
	     << "      " << clipname << endl
	     << "      " << isinfoname << endl
	     << "      " << umrdistname << endl
	     << "      " << scrdistname << endl
	     << "      " << split1name << endl
	     << "      " << split2name << endl;
	if (coverage_cutoff > 0)
	     cout << "      " << blistname << endl;
	if (generate_mapped_sc_um_file)
		cout << "      " << mapped_sc_name << endl;
		cout << "      " << mapped_um_name << endl;
	cout << "      " << prefix << "/ (directory)" << endl;

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

	if (mkdir(prefix.c_str(), 0777) < 0 && errno != EEXIST) {
		cout << "ERROR: could not create directory '" << prefix << "'" << endl;
		exit(1);
	}

	BamReader reader;
	if (!reader.Open(fname)) {
		cout << "ERROR: could not open BAM file '" << fname << "'" << endl;
		exit(1);
	}

	/*** First run through the BAM: record some basic statistics ***
	 *** and discover all soft clipped reads.  All of the soft   ***
	 *** clipped reads are stored in a map for the second run.   ***
	 *** GetNextAlignmentWithName() is a custom extension to the ***
	 *** BamTools library--it works nearly as fast as *Core, but ***
	 *** also fills in the read name.  We need to store the read ***
	 *** name so that we can identify mated pairs in the second  ***
	 *** pass.                                                   ***/
	BamAlignment al;
	long num_total = 0;
	long num_unmapped = 0;
	long num_clipped = 0;
	long max_read_len = 0;
	int bps = min(big_s_bps, frag_size);
	tr1::unordered_map<string, int> softclips;
	if (verbose)
		cout << "Beginning first pass..." << endl;
	while (reader.GetNextAlignmentWithName(al)) {
		++num_total;
		if (al.Length > max_read_len)
			max_read_len = al.Length;

		witness_qualities(al);

		/* UM+UM pairs and singleton UMs aren't interesting */
		if (!al.IsMapped() && (!al.IsMateMapped() || !al.IsPaired())) {
			++num_unmapped;
			continue;
		}

		/* Save the names of interesting read pairs.  "Interesting"
		 * is defined as: (1) the mate is mapped and (2) either this
		 * read contains a "big" S operation or is unmapped. */
		/* S can only occur in the first or last cigar op */
		bool extract = isBigS(al, al.CigarData.front(), bps) ||
			isBigS(al, al.CigarData.back(), bps) ||
			!al.IsMapped();
if (al.Name == "SRR034962.904741") cout << "extract=" << extract << endl;
		if (extract && al.IsMateMapped()) {
			softclips[al.Name] = getMateNumber(al);
			++num_clipped;
		}
	}
	cout << fname << ": " << endl
	     << "   " << num_total << " reads" << endl
	     << "   " << num_unmapped << " unaligned" << endl
	     << "   " << num_clipped << " clipped (" <<
		softclips.size() << " hashed)" << endl
	     << "   " << max_read_len << " bps in the longest read" << endl;
	report_qualities(verbose);

	/*** Second run through the file.  We need to output soft     ***
	 *** clips and their mates simultaneously to different files. ***
	 *** We use a map to buffer read mates which we know (using   ***
	 *** the previously constructed map) are either soft clipped  ***
	 *** or have soft clipped mates.  Once both mates are read,   ***
	 *** we output both reads and remove them from the buffer.    ***
	 *** All streams are gzipped streams to reduce I/O load.      ***/
	long unmapped_rejected = 0;
	long clipped_rejected = 0;
	long cur_read = 0, cur_pct = 1;
	long cur_pos = -1, cur_chrom = -1;
	long *cov;
	long covsize = sizeof(long)*max_read_len;
	ogzstream blist(blistname.c_str());
	if (coverage_cutoff > 0) {
		cov = (long *)malloc(covsize);
		bzero((char *)cov, covsize);
	}
	time_t start = time(NULL);
	reader.Rewind();
	ogzstream umfile(umfname.c_str());
	ogzstream clipfile(clipname.c_str());
	split1.open(split1name.c_str());
	split2.open(split2name.c_str());
	map<string, BamAlignment> albuf;
	map<string, BamAlignment> umbuf;
	map<string, ReadGroup *> rgs;
	Histogram uhist, chist;
	RefVector refnames = reader.GetReferenceData();
        BamWriter mapped_sc_file;
        BamWriter mapped_um_file;
	if (generate_mapped_sc_um_file) {
        	if (!mapped_sc_file.Open(mapped_sc_name,
                         	reader.GetHeaderText(), refnames)) {
                	cout << "ERROR: could not open output BAM file '"
                     	<< mapped_sc_name << "' for writing" << endl;
                	exit(1);
        	}
        	if (!mapped_um_file.Open(mapped_um_name,
                         	reader.GetHeaderText(), refnames)) {
                	cout << "ERROR: could not open output BAM file '"
                     	<< mapped_um_name << "' for writing" << endl;
                	exit(1);
        	}
	}

	if (verbose)
		cout << "Beginning second pass..." << endl;
	while(reader.GetNextAlignment(al)) {
		++cur_read;
		if (verbose)
			print_progress(cur_read, num_total, start);

		/* XXX: does not account for CIGAR string data like calDepth */
		/* these 4 Is* conditions match calDepth's default mask */
		if (coverage_cutoff > 0 &&
		    al.IsMapped() && al.IsPrimaryAlignment() &&
		    !al.IsFailedQC() && !al.IsDuplicate()) {
			if (al.RefID != cur_chrom && cur_chrom != -1) {
				for (int i = 0; i < max_read_len; ++i) {
					if (cov[i] >= coverage_cutoff) {
						blist << refnames[cur_chrom].RefName
						      << "\t" << cur_pos + i + 1
						      << "\t" << cov[i] << endl;
						cov[i] = 0;
					}
				}
			}
			else if (al.Position != cur_pos && cur_pos != -1) {
				/* compute the number of positions in our cov
				 * buffer that need to be kicked out */
				int x = min(al.Position - cur_pos,
				            max_read_len);
				for (int i = 0; i < x; ++i) {
					if (cov[i] >= coverage_cutoff) {
						blist << refnames[cur_chrom].RefName
						      << "\t" << cur_pos + i + 1
						      << "\t" << cov[i] << endl;
					}
				}

				memmove(cov, cov+x, covsize - sizeof(long)*x);
				bzero((char *)(cov + max_read_len - x),
					sizeof(long)*x);
			}
			cur_pos = al.Position;
			cur_chrom = al.RefID;
			for (int i = 0; i < al.Length; ++i)
				++cov[i];
		}

		/* Determine this alignment's readgroup */
		ReadGroup *rg;
		string rg_name;
		al.GetReadGroup(rg_name);
		map<string, ReadGroup *>::iterator rgit = rgs.find(rg_name);
		if (rgit == rgs.end()) {
			ReadGroup *x = new ReadGroup(al, max_isize,
				isize_samples, prefix, rg_blacklist);
			rgs[rg_name] = x;
			rg = x;
		} else {
			rg = rgit->second;
			rg->witness(al);
		}

if (al.Name == "SRR034962.904741") cout << "got read" << endl;
		if (rg->is_blacklisted())
			continue;

		/* 1. Write all unmapped reads to unmapped.fq */
		if (!al.IsMapped()) {
			trim_read(al, q);
			uhist.add(al);
			if (al.Length >= 2*frag_size) {
				BamAlignment copy(al);
				if (al.IsPaired())
					copy.Name += copy.IsFirstMate() ?
						"_1" : "_2";
				writeFQ((ofstream &)umfile, copy);
				split_read(copy, frag_size, n_cutoff);
			} else
				++unmapped_rejected;
		}

		/* 2. Handle UU pairs */
		if (!al.IsMapped() && !al.IsMateMapped() && processUU) {
			map<string, BamAlignment>::iterator x =
				umbuf.find(al.Name);
			if (x == umbuf.end())
				umbuf[al.Name] = al;

			BamAlignment mate = umbuf[al.Name];
			if (!isMatePair(al, mate))
				continue;

			/* Both reads are trimmed by (1) before we get here */
			rg->recordUU(al, mate, big_s_bps, n_cutoff);
			umbuf.erase(x);
			continue;
		}

		/* 3. Handle mapped+[softclipped|unmapped] read pairs */
		tr1::unordered_map<string, int>::iterator it =
			softclips.find(al.Name);
		if (it != softclips.end()) {
			map<string, BamAlignment>::iterator x =
				albuf.find(al.Name);
			if (x == albuf.end())
				albuf[al.Name] = al;
			else if (isMatePair(al, x->second)) {
				/* If we're here, then we've got both reads for
				 * some pair that was stored in softclips */
				BamAlignment mate = x->second;

				trim_read(al, q);
				trim_read(mate, q);
				chist.add(getMateNumber(al) == it->second ?
					al : mate);

				bool did_output = false;
				if (al.IsMapped() && mate.IsMapped()) {
					did_output = rg->recordSC(al, mate, it->second,
				    	    clipfile, big_s_bps, frag_size, n_cutoff);
					if (!did_output)
						++clipped_rejected;
					if (did_output && generate_mapped_sc_um_file) {
						mapped_sc_file.SaveAlignment(al);
						mapped_sc_file.SaveAlignment(mate);
					}
				} else {
					did_output = rg->recordUM(al, mate, big_s_bps, n_cutoff);
					if (!did_output)
						++clipped_rejected;
					if (generate_mapped_sc_um_file &&
					    al.Length >= 2*frag_size &&
					    mate.Length >= 2*frag_size) {
						mapped_um_file.SaveAlignment(al);
						mapped_um_file.SaveAlignment(mate);
					}
				}

				albuf.erase(x);
			}
		}
	}
	/* XXX: need to kick out all of the remaining values in cov */
	umfile.close();
	reader.Close();
	blist.close();
	split1.close();
	split2.close();

	cout << "Final read counts after rejection due to quality (-q "
	     << q << " -m " << min_read_len << "):" << endl
	     << "   unmapped reads: "
		<< num_unmapped - unmapped_rejected
		<< " (" << unmapped_rejected << " rejected)" << endl
	     << "   clipped reads:  "
		<< num_clipped - clipped_rejected
		<< " (" << clipped_rejected << " rejected)" << endl;

	ofstream umrdist(umrdistname.c_str());
	uhist.print(umrdist);
	umrdist.close();
	ofstream scrdist(scrdistname.c_str());
	chist.print(scrdist);
	scrdist.close();

	if (generate_mapped_sc_um_file) {
		mapped_um_file.Close();
		mapped_sc_file.Close();
	}


	ofstream isinfo(isinfoname.c_str());
	map<string, ReadGroup *>::iterator rgit;
	for (rgit = rgs.begin(); rgit != rgs.end(); ++rgit) {
		ReadGroup *x = rgit->second;
		cout << *x << endl;
		isinfo << "Read length:\t" << x->getName() << endl;
		isinfo << x->getReadlens()[0] << endl;
		delete x;
		rgs.erase(rgit);
	}
	isinfo.close();

	return 0;
}
