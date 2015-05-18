#include "Histogram.h"

Histogram::Histogram() :
	total(0)
{ }

void Histogram::add(BamTools::BamAlignment &al)
{
	++len_freqs[al.Length];
	++total;
}

void Histogram::print(std::ofstream &f)
{
	std::map<int, int>::iterator it = len_freqs.begin();

	while (it != len_freqs.end()) {
		f << it->first << "\t" << it->second << "\t"
		  << (double)it->second / (double)total << std::endl;
		++it;
	}
}
