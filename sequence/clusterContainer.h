//
// Created by borja on 7/06/19.
//

#ifndef ESPRELA_CLUSTERCONTAINER_H
    #define ESPRELA_CLUSTERCONTAINER_H
#endif //ESPRELA_CLUSTERCONTAINER_H

#include "../utils/utils.h"

using namespace sdsl;
using namespace std;

class ClusterContainer{
public:
    ClusterContainer(vector<string> reads, vector<string> f_environments, string directory, string binner, bool clustering, rank_support_v<1> rs_ds):
        _reads(reads),_f_environments(f_environments),_binner(binner),_directory(directory),_clustering(clustering),_rs_ds(rs_ds)
    {
        _doClustering();
        _loadClusters();
    }
    vector<vector<int>> getClusters();
    void classifyClusters();
    void selectKmers()
    {
        _processHybridClusters();
    }
private:
    void _doClustering()
    {
        if (_clustering) {
            cout << "Clustering with "<<_binner<<endl;
            if (_binner == "maxbin")
                _binnerExecution(_directory+"output");
            else if (_binner == "metabat")
                _binnerExecution(_directory+"output");
            else
            {
                cout << "Binner not accepted"<<endl;
                exit(1);
            }
        }else
            cout << "Skipping clustering"<<endl;
    }

    void _loadClusters()
    {
        set<string> cluster_files = System::getAllFaFqFiles(_directory);
        for (auto f: cluster_files)
        {
            cout << "File: "<<f<<endl;
            vector<int> current_cluster;
            ifstream  contigsFile (f, ios::out | ios::app | ios::binary);
            if (contigsFile.is_open())
            {
                string line;
                while ( getline (contigsFile,line) )
                {
                    if ( line[0] == '>'){
                        current_cluster.push_back(atoi(line.substr(1, line.size()).c_str()));
                    }
                }
            }
            _clusters.push_back(current_cluster);
        }
    }

    void _binnerExecution(string output)
    {
        string contig_file = "contigs.fasta", instruction;
        if (_binner == "maxbin")
        {
            instruction = "bash -c \"perl third-party/MaxBin/run_MaxBin.pl -contig " + contig_file + " ";
            instruction += "-reads " + _reads[0] + " ";
            for (size_t i = 1; i < _reads.size(); ++i)
                instruction += "-reads" + to_string(i + 1) + " " + _reads[i] + " ";
            instruction += "-thread 32 -out " + output + "\"";
            cout << "Instruction: " << instruction << endl;
        }
        if (_binner == "metabat")
        {
            string contig_file = "contigs.fasta";
            string instruction = "export PATH=$PATH:~/metabat/bin/;bash -c \"metabat2 -i "+contig_file+" -o "+output+"\"";
            cout << "Instruction: "<<instruction<<endl;
        }
        System::execute(instruction);
    }

    void _processHybridClusters()
    {
        set<string> cluster_files = System::getAllFaFqFiles(_directory);
        vector<string> cluster_files_vect;
        cluster_files_vect.insert(cluster_files_vect.end(), cluster_files.begin(), cluster_files.end());
        for (size_t i = 0; i < _hybridClusters.size(); i++)
            if (_hybridClusters[i])
                _selectSignificantKmersForClusters(cluster_files_vect[i]);
    }

    void _selectSignificantKmersForClusters(string clusterFile)
    {
        cout << "Process: "<<clusterFile<<endl;
    }
    vector<string> _reads,_f_environments;
    const string _binner, _directory;
    bool _clustering;
    rank_support_v<1> _rs_ds;
    vector<vector<int>> _clusters;
    vector<bool> _hybridClusters;
};
