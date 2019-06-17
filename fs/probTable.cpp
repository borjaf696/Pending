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
    _matrix_features[column][row] = ceil(value*RANGE_VALUES);
    _classes[row] = type;
    _environment[row] = environment;
    return ceil(value*RANGE_VALUES);
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
        _prob_matrix[i] = (probability*)malloc((RANGE_VALUES+1)*sizeof(probability));
        for (size_t j = 0; j < RANGE_VALUES; ++j)
            _prob_matrix[i][j] = 0.0;
    }
    _prob_classes = (probability*)malloc((RANGE_VALUES+1)*sizeof(probability));
}

void ProbTable::calculateProbTable(const FSContext & fs)
{
    for (size_t i = 0; i < _num_columns; ++i)
    {
        pair<bit*, size_t> feature_info = fs.getFeatureValues(i);
        probability jump = 1 /(double)feature_info.second;
        for (size_t j = 0; j < feature_info.second; ++j)
            _prob_matrix[i][feature_info.first[j]] += jump;
    }
    /*
     * Classes -.-"
     */
    {
        pair < bit * , size_t > classes_info = fs.getClasses();
        probability jump = 1 / (double) classes_info.second;
        for (size_t j = 0; j < classes_info.second; ++j)
            _prob_classes[classes_info.first[j]] += jump;
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
    for (size_t i = 0; i < range1*range2; i++)
        jointProb[i] = 0;
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
    bool same_feature = true;
    probability * jointProb = JointProbability::getJointProbability(feature1, feature2, num_rows), mutualInfo = 0;
    for (size_t i = 0; i < range1; ++i)
        for (size_t j = 0; j < range2; ++j) {
            probability jProb = jointProb[i * range2 + j], marg1 = probs1[i], marg2 = probs2[j];
            if (same_feature)
                same_feature = (jProb == (double)1);
            if (jProb != 0 && marg1 != 0 && marg2 != 0)
                mutualInfo +=  jProb* log2(jProb / (marg1 * marg2));
        }
    /*
     * Clear jointProb
     */
    delete jointProb;
    return (same_feature)?MAX_REDUNCANCY:mutualInfo;
}
/*
 * Fast-mRMR
 */

fastmRMR::fastmRMR(ProbTable * pt, FSContext * fsc, size_t num_env):_num_features(pt->getNumColumns()),_best_feature(0),_pt(pt), _fsc(fsc)
{
    _relevance = (probability *) malloc(_num_features*sizeof(probability));
    _redundance = (probability *) malloc(_num_features*sizeof(probability));
    for (size_t i = 0; i < _num_features; ++i)
    {
        _relevance[i] = 0;
        _redundance[i] = 0;
    }
    _getRelevances();
    /*
     * Maximal expected MI num_env*log(num_env*num_env)
     */
    if (num_env < 2)
        _max_MI = MAX_MI;
    _max_MI = num_env*log2(num_env*num_env);
    cout << "Maximal expected relevance: "<<_max_MI<<endl;
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
    probability best_relevance = MINIMAL;
    for (size_t i = 0; i < _num_features; ++i)
    {
        Progress::update(i);
        pair<bit *, size_t> feature = _fsc->getFeatureValues(i), classes = _fsc->getClasses();
        probability * prob1 = _pt->getProbabilities(i), * probl_class = _pt->getClassesProbabilities();
        _relevance[i] = MutualInfo::getMutualInformation(feature.first, classes.first, prob1, probl_class,feature.second);
        if (best_relevance < _relevance[i])
        {
            _best_feature = i;
            best_relevance = _relevance[i];
        }
    }
    Progress::update(_num_features);
}

vector<size_t> fastmRMR::getFeatures(probability rate)
{
    vector<size_t> features;
    probability first_relevance = _relevance[_best_feature], last_relevance = _relevance[_best_feature];
    while (last_relevance > (first_relevance*rate))
    {
        features.push_back(_best_feature);
        _redundance[_best_feature] = MAX_REDUNCANCY;
        last_relevance = _updateRedundances(features);
    }
    return features;
}

probability fastmRMR::_updateRedundances(vector <size_t> & selected_features)
{
    /*
     * Progress
     */
    {
        Progress::get().size_total = _num_features;
        Progress::get().show = true;
    }
    probability best_relevance = 0, best_relation = MINIMAL;
    pair<bit *, size_t> feature1 = _fsc->getFeatureValues(_best_feature);
    probability  * prob1 = _pt->getProbabilities(_best_feature);
    for (size_t i = 0; i < _num_features; i++)
    {
        Progress::update(i);
        if (_redundance[i] <= MAX_REDUNCANCY)
        {
            pair < bit * , size_t > feature2 = _fsc->getFeatureValues(i);
            probability * prob2 = _pt->getProbabilities(i);
            probability redundance = MutualInfo::getMutualInformation(feature1.first, feature2.first, prob1, prob2, feature1.second);
            if (redundance == MAX_REDUNCANCY)
                selected_features.push_back(i);
            _redundance[i] += redundance;
            if ((_relevance[i] - _redundance[i]) > best_relation)
            {
                _best_feature = i;
                best_relation = (_relevance[i] -_redundance[i]);
                best_relevance = _relevance[i];
            }
        }
    }
    Progress::update(_num_features);
    return best_relevance;
}
