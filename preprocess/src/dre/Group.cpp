/*
 * A "group" is similar to a single position in a pile up.
 * The group contains all reads starting at one position;
 * but these reads may be on either strand and the mates
 * could be anywhere.
 */
#include <algorithm>

#include "Group.h"

using namespace BamTools;

Group::Group() : n(0)
{
}

bool
is_dup(BamAlignment &a, BamAlignment &b)
{
	long a_mpos = a.IsReverseStrand() ?
		-a.MatePosition : a.MatePosition;
	long b_mpos = b.IsReverseStrand() ?
		-b.MatePosition : b.MatePosition;
	return a.MateRefID == b.MateRefID && a_mpos == b_mpos;
}

void
Group::add(BamAlignment &al)
{
	if (n == 0)
		pos = al.Position;
	else if (al.Position != pos) {
		cout << "ERROR: added read with pos " << al.Position
		     << " to group with pos=" << pos << endl;
		exit(1);
	}

	/* Is this read a dup of any existing read in the group? */
	for (size_t i = 0; i < reads.size(); ++i) {
		if (is_dup(reads[i], al)) {
			/* Yep, make a new PCRDup group */
			dup_groups.push_back(PCRDups(reads[i], al));
			/* Remove read i from the main list */
			reads.erase(reads.begin() + i);
			return;
		}
	}

	/* Is this read a dup of any existing dup group? */
	for (size_t i = 0; i < dup_groups.size(); ++i) {
		if (is_dup(dup_groups[i].rep(), al)) {
			dup_groups[i].add(al);
			return;
		}
	}

	/* Ok, it isn't a dup of any kind.  Store it. */
	reads.push_back(al);
	++n;
}

/*
 * Call this function when all reads at this position have
 * been cycled through.  This will write out all non-dup
 * reads.
 * Reads in dup groups will be handled differently: if the
 * dup group has been seen before (e.g., we're looking at
 * mates of a previous dup group), consult the record to
 * determine which read should be written; otherwise, we
 * choose a 'best' read, write it out and save it to the
 * record for when we see this dup block's mate block.
 */
int
Group::write(BamWriter &writer, record &prev_dups, BamWriter *dups)
{
	int n = 0;
	long pos;

	pos = reads[0].Position;

	/* Easy part: write out non dups */
	for (size_t i = 0; i < reads.size(); ++i) {
		if (reads[i].Position != pos) {
			cout << "ERROR: group with pos=" << pos
			     << "has read with pos=" << reads[i].Position
			     << endl;
			exit(1);
		}
		writer.SaveAlignment(reads[i]);
		++n;
	}

	/* Write the best read from each dup group */
	for (size_t i = 0; i < dup_groups.size(); ++i) {
		BamAlignment best = dup_groups[i].best(prev_dups);
		if (best.Position != pos) {
			cout << "ERROR: group with pos=" << pos
			     << "has dup group with pos=" << best.Position
			     << endl;
			exit(1);
		}
		/* if this is the mate block of a PCR dup block, then we */
		/* -should- find the same-named mate.                    */
		record::iterator it = prev_dups.find(make_frag_coords(best));
		if (it != prev_dups.end() && best.Name != it->second) {
			cout << "adding dup (" << best.Name
			     << ") to record where one already exists! ("
			     << it->second << ")" << endl;
			exit(1);
		}
		prev_dups[make_frag_coords(best)] = best.Name;
		writer.SaveAlignment(best);
		if (dups != NULL) 
			dup_groups[i].write_except(dups, best);
		++n;
	}

	return n;
}

void
Group::clear()
{
	n = 0;
	reads.clear();
	dup_groups.clear();
}

bool
Group::should_write(BamAlignment &al)
{
	if (reads.empty() && dup_groups.empty())
		return false;

	BamAlignment x = reads.empty() ? dup_groups[0].rep() : reads[0];
	if (al.RefID != x.RefID || al.Position != x.Position)
		return true;

	return false;
}

frag_coords
make_frag_coords(BamTools::BamAlignment &al)
{
	/* The positions need to be tagged with + and - for discordant */
	/* mate pairs that are on the same strand (i.e., ++ pairs and --pairs */
	coords a(al.RefID, (al.IsReverseStrand() ? -1:1) * al.Position),
	       b(al.MateRefID, (al.IsMateReverseStrand() ? -1:1) * al.MatePosition);
	
	/* If the two reads are on different strands */
	if (al.IsReverseStrand() != al.IsMateReverseStrand())
		return al.IsReverseStrand() ? frag_coords(b, a) : frag_coords(a, b);
	/* If the two reads are on the same strand... */
	else
		return a < b ? frag_coords(a, b) : frag_coords(b, a);
}

