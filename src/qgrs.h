
#include <string>
#include <queue>
#include <cmath>
using namespace std;

typedef unsigned long nt;
class G4Candidate;
class G4;

string findJSON(string, bool overlaps = false, short min_tetrads=2, short min_score=17);
vector<G4> find(string sequence, bool overlaps = false, short min_tetrads = 2, short min_score = 17);

class G4Candidate {
public:
    G4Candidate(const string & sequence, short tetrads, nt start_pos);
    short score() {
        double gavg = (abs(y1-y2) + abs(y2-y3) + abs(y1-y3))/3.0;
        return floor(gmax() - gavg + gmax() * (numTetrads-2));
    }
    short gmax(){
        return maxLength - (numTetrads * 4 + 1);
    }

    short length() {
       return 4 * numTetrads + y1 + y2 + y3;
    }

    nt t1(){
        return start;
    }

    nt t2() {
        return t1() + numTetrads + y1;
    }

    nt t3() {
        return t2() + numTetrads + y2;
    }

    nt t4() {
        return t3() + numTetrads + y3;
    }

    nt cursor() {
        if (y1 < 0 ) return t1() + numTetrads;
        else if (y2 < 0) return t2() + numTetrads;
        else if (y3 < 0) return t3() + numTetrads;
        else return -1;
    }

    short partialLength();
    short minAcceptableLoopLength(){
        if (y1 == 0 || y2 == 0 || y3 == 0) return 1;
        else return 0;
    }

    bool complete() {
        if (y1 < 0 || y2 < 0 || y3 < 0 ) return false;
        return true;
    }

    short y1 ;
    short y2 ;
    short y3 ;
    string sequence;
    short numTetrads;
    nt start ;
    string tstring;
    short maxLength;

    string toString();
    bool viable(int min_score);
    void findLoopLengthsFrom(queue<int> & ys, int i) ;
    void expand(queue<G4Candidate> &cands) ;
};

class G4 {
public:
    G4();
    G4(G4Candidate &candidate);
    bool isequal(const G4 & other);
    void makejson(string name, nt value, stringstream &out); 
    void makejson(string name, short value, stringstream &out);
    void makejson(string name, string value, stringstream &out);
    string toJSON(bool print_overlaps=true);
    nt start;
    nt tetrad1;
    nt tetrad2;
    nt tetrad3;
    nt tetrad4;
    short y1;
    short y2;
    short y3;
    short tetrads;
    short length;
    short gscore;
    string sequence;
    vector<G4> overlaps;
};

bool operator< (const G4 &left, const G4 &right);