#include <iostream>
#include <map>
#include "api/BamAlignment.h"

class Histogram {
public:
	Histogram();

	void add(BamTools::BamAlignment &al);
	void print(std::ofstream &f);

private:
	std::map<int, int> len_freqs;
	int total;
};
