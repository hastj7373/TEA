#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include "api/BamAlignment.h"

class PCRDups {
public:
	PCRDups(BamTools::BamAlignment &a, BamTools::BamAlignment &b) {
		dups.push_back(a);
		dups.push_back(b);
	}

	/* Does 'al' belong to this PCR duplicate group? */
	void add(BamTools::BamAlignment &al) {
		dups.push_back(al);
	}

	/* always exists, because PCRDups must have >= 2 elements */
	BamTools::BamAlignment &rep() {
		return dups[0];
	}

	/* For now, the "best" alignment is just the first one in the list */
	BamTools::BamAlignment &best(record &prev_dups) {
		frag_coords x = make_frag_coords(rep());
		record::iterator it = prev_dups.find(x);

		/* No previous matching duplicate block, just return rep */
		if (it == prev_dups.end())
			return rep();

		/* We already saw the matching dup block, pick the mate of
		 * the same read we picked before. */
		string prev_name = it->second;
		for (size_t i = 0; i < dups.size(); ++i) {
			if (dups[i].Name == prev_name)
				return dups[i];
		}

		cout << "ERROR: could not find a best read in dup group at <"
		     << rep().RefID << ", " << rep().Position << ">" << endl;
		cout << "Previous dup block contained mate: " << it->second << endl;
		cout << "This block contains dups:" << endl;
		for (size_t i = 0; i < dups.size(); ++i)
			cout << dups[i].Name << "\t" << dups[i].Position
			     << "\tmate:" << dups[i].MatePosition << endl;
		exit(1);
	}

	void
	write_except(BamTools::BamWriter *out, BamTools::BamAlignment &best) {
		for (size_t i = 0; i < dups.size(); ++i)
			if (dups[i].Name != best.Name)
				out->SaveAlignment(dups[i]);
	}

private:
	vector<BamTools::BamAlignment> dups;
};
