//
// Created by borja on 9/06/19.
//

#ifndef ESPRELA_PROBTABLE_H
#define ESPRELA_PROBTABLE_H
    #ifndef  UTILS_H
    #include "../utils/utils.h"
    #endif
    #define KMER_SIZE 50
    #define RANGE_VALUES 9 /*Always equal or higher the number of environments*/
    #define MAX_REDUNCANCY 999999
    #define MAX_MI 999999
    #define MINIMAL -99999

    #include <math.h>
    #include <vector>
    #include <unordered_map>
    #include "../sequence/kmer.h"

    using namespace std;
    typedef uint8_t bit;
    typedef double probability;

    class FSContext{
    public:
        FSContext(){}
        FSContext(size_t columns, size_t rows);
        ~FSContext()
        {
            for (size_t i = 0; i < _num_columns;++i)
                free(_matrix_features[i]);
            free(_matrix_features);
            free(_environment);
            free(_classes);
        }
        bit update(size_t, size_t, probability, bit, bit);
        pair<bit *, size_t> getFeatureValues(size_t) const;
        pair<bit *, size_t> getClasses() const
        {
            return pair<bit*, size_t>(_classes, _num_rows);
        }
        vector<size_t> getFeature();
        void increaseMatrix(size_t, size_t);
    private:
        bit ** _matrix_features, * _classes, * _environment;
        size_t _num_columns, _num_rows;
    };

    class ProbTable{
    public:
        ProbTable(){}
        ProbTable(size_t);
        ~ProbTable()
        {
            free(_prob_matrix);
            free(_prob_classes);
        }
        probability * getProbabilities(size_t);
        probability * getClassesProbabilities();
        size_t getNumColumns()
        {
            return _num_columns;
        }
        void calculateProbTable(const FSContext &);
    private:
        probability ** _prob_matrix;
        probability * _prob_classes;
        size_t _num_columns;
    };

    class JointProbability{
    public:
        static probability * getJointProbability(bit *, bit *, size_t
                , size_t = (RANGE_VALUES+1), size_t = (RANGE_VALUES+1));
    };

    class MutualInfo{
    public:
        static probability getMutualInformation(bit *, bit *, probability *, probability *
                , size_t, size_t = (RANGE_VALUES+1), size_t = (RANGE_VALUES+1));
    };

    class fastmRMR
    {
    public:
        fastmRMR(ProbTable *, FSContext *, size_t);
        ~fastmRMR()
        {
            free(_relevance);
            free(_redundance);
        }
        vector<size_t> getFeatures(probability);
    private:
        probability _updateRedundances(vector<size_t>&);
        void _getRelevances();
        size_t _num_features, _best_feature;
        probability * _relevance, * _redundance, _max_MI;
        ProbTable * _pt;
        FSContext * _fsc;
    };
#endif