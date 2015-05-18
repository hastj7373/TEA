#include <iostream>
#include <algorithm>
#include <vector>

#include "ReadGroup.h"
#include "bamreader.h"

using namespace std;
using namespace BamTools;

ReadGroup::ReadGroup() { }

ReadGroup::~ReadGroup()
{
	if (!blacklisted) {
		f1.close();
		f2.close();
	}

	/* Write the insert size distribution to file */
	ofstream f((prefix + "/" + name + ".is").c_str());
	for (int i = 0; i < inserts.size(); ++i)
		f << inserts[i] << "\n";
	f.close();
}

ReadGroup::ReadGroup(BamAlignment &al, int max_isize, int isize_samples,
	string prefix, list<string> blacklist) :
	max_isize(max_isize),
	isize_samples(isize_samples),
	prefix(prefix),
	blacklisted(false)
{
	if (!al.GetReadGroup(name))
		name = "none";

	nreads = 0;

	/* Determine if this read group is in the blacklist */
	for (list<string>::iterator it = blacklist.begin();
	     it != blacklist.end(); ++it) {
		if (*it == name) {
			blacklisted = true;
			break;
		}
	}

	if (!blacklisted) {
		f1.open((prefix + "/" + name + "_1.fq.gz").c_str());
		f2.open((prefix + "/" + name + "_2.fq.gz").c_str());
	}

	witness(al);
}

void
ReadGroup::witness(BamAlignment &al)
{
	++nreads;

	/* Selects first 'isize_samples' samples instead of random sampling */
	/* Don't use negative insert sizes.  If it's negative, then the mate  *
	 * has a positive insert of the same size.  If both are recorded, the *
	 * average will always be 0.                                          */
	if (al.InsertSize > 0 && al.InsertSize <= max_isize &&
	    inserts.size() < isize_samples)
		inserts.push_back(al.InsertSize);

	if (find(readlens.begin(), readlens.end(), al.Length) == readlens.end()) {
		if (readlens.size() > 0) {
			cerr << "ERROR: differing read lengths in read group '"
			     << name << "'" << endl;
		}
		readlens.push_back(al.Length);
	}
}

/* Assumes at least one of 'a' and 'b' are 'S' ops */
bool
copcomp(CigarOp &a, CigarOp &b)
{
	if (b.Type != 'S')
		return false;
	if (a.Type != 'S')
		return true;

	return a.Length < b.Length;
}

void
clipAlignment(BamAlignment &al)
{
	int offset, length;
	CigarOp cop1 = al.CigarData[0];
	CigarOp cop2 = al.CigarData[al.CigarData.size() - 1];

	if (copcomp(cop2, cop1)) {
		offset = 0;
		length = min(al.Length, (signed)cop1.Length);
	} else {
		offset = al.Length - min(al.Length, (signed)cop2.Length);
		length = min(al.Length, (signed)cop2.Length);
	}

	try {
		al.Qualities = al.Qualities.substr(offset, length);
		al.QueryBases = al.QueryBases.substr(offset, length);
	} catch (exception &e) {
		cout << "ERROR: substr failed in clipAlignment()" << endl;
		cout << al.Name << " " << (al.IsReverseStrand() ? "(-)" : "(+)");
		cout << " offset: " << offset << " length: " << length
		     << " taglen: " << al.Length << endl;
		cout << "cop1: " << cop1.Length << cop1.Type << endl;
		cout << "cop2: " << cop2.Length << cop2.Type << endl;
		exit(1);
	}
}

/*
 * snip() doesn't leave a valid BamAlignment; it contains
 * correct FASTQ data.  Handles negative strand alignments:
 * 'start=0' will always correspond to the 5'-most basepair
 * in the read.
 */
BamAlignment
snip(BamAlignment &a, int start, int len)
{
	BamAlignment copy(a);

	/* Handle reverse strand mappings */
	int converted_start = copy.IsReverseStrand() ?
		copy.Length - start - len : start;
	copy.Length = len;
	try {
		copy.QueryBases = copy.QueryBases.substr(converted_start, len);
		copy.Qualities = copy.Qualities.substr(converted_start, len);
	} catch (exception &e) {
		cout << "ERROR: substr failed in snip(" << a.Name << ", "
		     << start << ", " << len << ")" << endl;
		cout << (a.IsReverseStrand() ? "(-)" : "(+)")
		     << ", converted_start: " << converted_start << endl;
		cout << a.QueryBases << endl;
		cout << a.Qualities << endl;
		exit(1);
	}
	return copy;
}

/*
 * The alignments a and b are already trimmed by quality score.
 * Create read pairs out of the front and back 'big_s_bps' bps
 * of each read, provided the reads are big enough to supply
 * 2*big_s_bps nonoverlapping bases.  This can produce up to 4
 * pairs: FF, FB, BF, BB (F - front fragment, B - back fragment).
 * Lixing wants these named uu1-uu4.
 */
bool
ReadGroup::recordUU(BamAlignment &a, BamAlignment &b, int big_s_bps, int n_cutoff)
{
	/* Not even an FF pair can be produced */
	if (a.Length < big_s_bps || b.Length < big_s_bps)
		return false;

	/* Produce the FF pair */
	BamAlignment af = snip(a, 0, big_s_bps);
	af.Name += "uu1";
	BamAlignment bf = snip(b, 0, big_s_bps);
	bf.Name += "uu1";
	writeFQpair((ofstream &)f1, af, (ofstream &)f2, bf, n_cutoff);

	/* If the read is long enough, try back fragments */
	BamAlignment ab;
	BamAlignment bb;
	bool aback = false, bback = false;
	if (a.Length >= 2*big_s_bps) {
		/* BF pair */
		aback = true;
		ab = snip(a, a.Length - big_s_bps, big_s_bps);
		ab.Name += "uu2";
		bf.Name = b.Name + "uu2";
		writeFQpair((ofstream &)f1, ab, (ofstream &)f2, bf, n_cutoff);
	}
	if (b.Length >= 2*big_s_bps) {
		/* FB pair */
		bback = true;
		bb = snip(b, b.Length - big_s_bps, big_s_bps);
		bb.Name += "uu3";
		af.Name = a.Name + "uu3";
		writeFQpair((ofstream &)f1, af, (ofstream &)f2, bb, n_cutoff);
	}
	if (aback && bback) {
		/* BB pair */
		ab.Name = a.Name + "uu4";
		bb.Name = b.Name + "uu4";
		writeFQpair((ofstream &)f1, ab, (ofstream &)f2, bb, n_cutoff);
	}

	return true;
}

bool
ReadGroup::recordUM(BamAlignment &a, BamAlignment &b, int big_s_bps, int n_cutoff)
{
	BamAlignment mapped = a.IsMapped() ? a : b;
	BamAlignment unmapped = a.IsMapped() ? b : a;

	if (unmapped.Length < big_s_bps || mapped.Length < 2*big_s_bps)
		return false;

	mapped.Name += "mu";
	unmapped.Name += "mu";

	/* unmapped >= big_s_bps.  Make a fragment from the 5' end */
	BamAlignment cp1 = snip(unmapped, 0, big_s_bps);
	cp1.Name += "1";
	mapped.Name += "1";
	writeFQpair((ofstream &)f1, mapped, (ofstream &)f2, cp1, n_cutoff);

	if (unmapped.Length < 2*big_s_bps)
		return true;

	/* unmapped >= 2*big_s_bps.  Make another fragment from the 3' end */
	BamAlignment cp2 =
		snip(unmapped, unmapped.Length - big_s_bps, big_s_bps);
	cp2.Name += "2";
	mapped.Name[mapped.Name.size() - 1] = '2';
	writeFQpair((ofstream &)f1, mapped, (ofstream &)f2, cp2, n_cutoff);

	return true;
}

bool
ReadGroup::recordSC(BamAlignment &a, BamAlignment &b, int clippedMate,
	ogzstream &clipfile, int big_s_bps, int frag_size, int n_cutoff)
{
	bool did_write = false;

	/* Make a copy of each alignment so we don't modify the base data */
	BamAlignment clipped, mate;
	if (clippedMate == getMateNumber(a)) {
		clipped = a;
		mate = b;
	} else {
		clipped = b;
		mate = a;
	}

	/* The trim could've ruined the big S op */
	if (clipped.CigarData.empty())
		return false;

	/* Handle two completely different filesets and parameters. */
	/* The softclips.fq file uses 'frag_size' to determine if a
	 * read is a valid softclip. */
	if ((isBigS(clipped, clipped.CigarData.front(), frag_size) ||
	     isBigS(clipped, clipped.CigarData.back(), frag_size)) &&
	    clipped.Length >= 2*frag_size && mate.Length >= 2*frag_size) {
		did_write = true;
		writeFQ((ofstream &)clipfile, clipped);
		split_read(clipped, frag_size, n_cutoff);
	}

	/* The _1 and _2 files use 'big_s_bps' to determine if a read
	 * is a valid softclip. */
	if ((isBigS(clipped, clipped.CigarData.front(), big_s_bps) ||
	     isBigS(clipped, clipped.CigarData.back(), big_s_bps)) &&
	    mate.Length >= 2*big_s_bps) {
		// For now, the -D option code only wants the previous case
		//did_write = true;
		BamAlignment copy(clipped);
		clipAlignment(copy);
		copy.Name += "sc";
		mate.Name += "sc";
		writeFQpair((ofstream &)f1, copy, (ofstream &)f2, mate, n_cutoff);
	}

	return did_write;
}

ostream &
std::operator<< (ostream &out, ReadGroup &rg)
{
	return out << "Readgroup [" + rg.getName() + "]: "
	    << rg.getNReads() << " reads";
}

string ReadGroup::getName() {
	return name;
}

int ReadGroup::getNReads() {
	return nreads;
}

vector<int> ReadGroup::getReadlens() {
	return readlens;
}

bool ReadGroup::is_blacklisted() {
	return blacklisted;
}
