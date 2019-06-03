//
// Borja :)
//
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;
/*
 * Class to build the tetranucleotide frequency, aka genomic composition
 */
class Composition
{
public:
    Composition(){}
    int getVal(char c)
    {
        return _dnaTable[(size_t) c];
    }
protected:
    static vector<int> _dnaTable;
    struct TableFiller
    {
        TableFiller()
        {
            static bool tableFilled = false;
            if (!tableFilled)
            {
                tableFilled = true;
                _dnaTable.assign(256, -1);
                _dnaTable[(size_t)'A'] = 0;
                _dnaTable[(size_t)'a'] = 0;
                _dnaTable[(size_t)'C'] = 1;
                _dnaTable[(size_t)'c'] = 1;
                _dnaTable[(size_t)'G'] = 2;
                _dnaTable[(size_t)'g'] = 2;
                _dnaTable[(size_t)'T'] = 3;
                _dnaTable[(size_t)'t'] = 3;
            }
        }
    };
    static TableFiller _filler;
};

class NaiveComposition: public Composition
{
public:
    NaiveComposition(string contig)
    {
        vector<int> curr_state(4,-1);
        for (auto c:contig)
        {
            if (getVal(c) != -1)
            {
                for (size_t i = 0; i < (curr_state.size()-1); ++i)
                {
                    curr_state[i+1] = curr_state[i]+getVal(c)*pow(4,i);
                }
                if (curr_state[3] > 0)
                    _composition[curr_state[3]]++;
                curr_state[0] = getVal(c);
            }else
                curr_state[0] = -10000;
        }
    }

    size_t getFreqPos(uint8_t pos)
    {
        return _composition[pos];
    }
private:
    vector<size_t> _composition = vector<size_t>(256,0);
};