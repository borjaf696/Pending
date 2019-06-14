//
// Created by borja on 9/06/19.
//
#include "probTable.h"

/*
 * FSContext
 */
FSContext::FSContext(size_t columns, size_t rows):_num_columns(columns),_num_rows(rows)
{
    /*
     * Create matrix
     */
    _matrix_features = (bit**) malloc(_num_columns*sizeof(bit*));
    for (size_t i = 0; i < _num_columns; ++i) {
        _matrix_features[i] = (bit *) malloc (_num_rows*sizeof(bit));
        for (size_t j = 0; j < _num_rows; ++j) {
            _matrix_features[i][j] = 0;
        }
    }
    _classes = (bit*)malloc(_num_rows);
    _environment = (bit*)malloc(_num_rows);
    for (size_t i = 0; i < _num_rows; ++i)
    {
        _classes[i] = 0;_environment[i] = 0;
    }
}

bit FSContext::update(size_t column, size_t row, probability value, bit type, bit environment)
{
    _matrix_features[column][row] = ceil(value*10);
    _classes[row] = type;
    _environment[row] = environment;
    return ceil(value*10);
}

pair<bit*, size_t>  FSContext::getFeatureValues(size_t column) const
{
    return pair<bit*, size_t>(_matrix_features[column], _num_rows);
}

/*
 * ProbTable
 */
ProbTable::ProbTable(size_t num_columns):_num_columns(num_columns)
{
    _prob_matrix = (probability**) malloc(_num_columns*sizeof(probability*));
    for (size_t i = 0; i < _num_columns; ++i)
    {
        _prob_matrix[i] = (probability*)malloc(RANGE_VALUES*sizeof(probability));
        for (size_t j = 0; j < RANGE_VALUES; ++j)
            _prob_matrix[i][j] = 0.0;
    }
    _prob_classes = (probability*)malloc(RANGE_VALUES*sizeof(probability));
}

void ProbTable::calculateProbTable(const FSContext & fs)
{
    for (size_t i = 0; i < _num_columns; ++i)
    {
        pair<bit*, size_t> feature_info = fs.getFeatureValues(i);
        probability jump = 1 /(double)feature_info.second;
        for (size_t j = 0; j < feature_info.second; ++j)
            _prob_matrix[i][j] += jump;
    }
    /*
     * Classes -.-"
     */
    {
        pair < bit * , size_t > classes_info = fs.getClasses();
        probability jump = 1 / (double) classes_info.second;
        for (size_t j = 0; j < classes_info.second; ++j)
            _prob_classes[j] += jump;
    }
}

probability * ProbTable::getProbabilities(size_t column)
{
    return _prob_matrix[column];
}

probability * ProbTable::getClassesProbabilities()
{
    return _prob_classes;
}

/*
 * Joint Probability
 */
probability * JointProbability::getJointProbability(bit * feature_one, bit * feature_two
        , size_t num_rows, size_t range1, size_t range2)
{
    probability * jointProb = new probability[range1*range2], jump = 1/(double)num_rows;
    for (size_t i = 0; i < num_rows; ++i)
        jointProb[feature_one[i]*range2+feature_two[i]]+=jump;
    return jointProb;
}

/*
 * Mutual Information
 */
probability MutualInfo::getMutualInformation(bit * feature1, bit * feature2
        , probability * probs1 , probability * probs2, size_t num_rows, size_t range1, size_t range2)
{
    probability * jointProb = JointProbability::getJointProbability(feature1, feature2, num_rows), mutualInfo = 0;
    for (size_t i = 0; i < range1; ++i)
        for (size_t j = 0; j < range2; ++j) {
            probability jProb = jointProb[i * range2 + j], marg1 = probs1[i], marg2 = probs2[j];
            if (jProb != 0 && marg1 != 0 && marg2 != 0)
                mutualInfo +=  jProb* log2(jProb / (marg1 * marg2));
        }
    /*
     * Clear jointProb
     */
    delete jointProb;
    return mutualInfo;
}
/*
 * Fast-mRMR
 */

fastmRMR::fastmRMR(ProbTable * pt, FSContext * fsc):_num_features(pt->getNumColumns()),_pt(pt), _fsc(fsc)
{
    _relevance = (bit*) malloc(_num_features*sizeof(bit));
    _redundance = (bit*) malloc(_num_features*sizeof(bit));
    _getRelevances();
}
void fastmRMR::_getRelevances()
{
    /*
     * Progress
     */
    {
        Progress::get().size_total = _num_features;
        Progress::get().show = true;
    }
    for (size_t i = 0; i < _num_features; ++i)
    {
        Progress::update(i);
        pair<bit *, size_t> feature = _fsc->getFeatureValues(i), classes = _fsc->getClasses();
        probability * prob1 = _pt->getProbabilities(i), * probl_class = _pt->getClassesProbabilities();
        _relevance[i] = MutualInfo::getMutualInformation(feature.first, classes.first, prob1, probl_class,feature.second);
    }
    Progress::update(_num_features);
}
