#include <iostream>
#include "qgrs.h"
using namespace std;

int main() {
	string input = "GGGGAGGGGGAGGGGGGGAGGGGGG";
  	vector<G4> results = find(input, true, 2, 17);
  	for (std::vector<G4>::const_iterator i = results.begin(); i != results.end(); ++i) {
    	std::cout << i->sequence << endl;
    	for (auto &overlap : i->overlaps) {
    		cout << " ~ " << overlap.sequence << endl;
    	}
    }
}