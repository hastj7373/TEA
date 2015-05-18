#include <vector>
#include <utility>
#include <map>
#include "api/BamAlignment.h"
#include "api/BamWriter.h"

using namespace std;

typedef pair<long, long> coords;
typedef pair<coords, coords> frag_coords;
typedef map<frag_coords, string> record;

frag_coords make_frag_coords(BamTools::BamAlignment &al);

#include "PCRDups.h"

class Group {
public:
	Group();
	void add(BamTools::BamAlignment &al);

	/* Filter out PCR duplicates and write out non-dup reads */
	int write(BamTools::BamWriter &writer, record &prev_dups,
		BamTools::BamWriter *dups);

	/* Empty out this group */
	void clear();

	/* Should we write out the current group and make a new one? */
	bool should_write(BamTools::BamAlignment &al);

private:
	vector<BamTools::BamAlignment> reads;
	vector<PCRDups> dup_groups;
	int n;
	long pos;
};
