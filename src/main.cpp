#include "../operations/GenomicComposition.hpp"
#include "../sequence/sequence_container.h"
#include "../sequence/clusterContainer.h"

namespace po = boost::program_options;
namespace
{
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
}
using namespace sdsl;
using namespace std;

void parse_args(int argc, char ** argv, vector<string> * directories, vector<string> * environment,
        vector<string> * output, vector<string> * preffixes, bool * clustering, string * binner)
{
    po::options_description des("Options");
    des.add_options()
            ("help,h","Help message")
            ("dir,d",po::value<vector<string>>(directories)->multitoken(),"Reads directories one per environment")
            ("environment,e",po::value<vector<string>>(environment)->multitoken(),"Environment which belong each sample")
            ("output,o",po::value<vector<string>>(output)->multitoken(),"Preffixes Output directories")
            ("clusters,c",po::bool_switch(clustering),"Internal parameter")
            ("preffix,p",po::value<vector<string>>(preffixes)->multitoken(),"Preffixes where metaspades has been run - environment 1")
            ("binner,b",po::value<string>(binner), "Binner to produce clusters (metabat/maxbin)");
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc,argv).options(des).run(),vm);
        if ( vm.count("help")  )
        {
            std::cout << "EGTs" << endl << endl;
            rad::OptionPrinter::printStandardAppDesc("EGTs",
                                                     cout,
                                                     des);
            exit(SUCCESS);
        }
        po::notify(vm); /*Required commands*/
    }catch(boost::program_options::required_option& e)
    {
        rad::OptionPrinter::formatRequiredOptionError(e);
        std::cerr << "ERROR: " << e.what() << endl << endl;
        rad::OptionPrinter::printStandardAppDesc("EGTs",
                                                 cout,
                                                 des);
        exit(ERROR_IN_COMMAND_LINE);
    }
    catch(boost::program_options::error& e)
    {
        std::cerr << "ERROR: " << e.what() << endl << endl;
        rad::OptionPrinter::printStandardAppDesc("EGTs",
                                                 cout,
                                                 des);
        exit(ERROR_IN_COMMAND_LINE);
    }
}

vector<string> spades_execution(vector<string> files, string output)
{
    vector<string> output_dirs;
    size_t jump = 1;
    for (size_t i = 0; i < files.size(); i+=jump)
    {
        string output_file = output+to_string(i);
        output_dirs.push_back(output_file);
        string instruction = "";
        if (Bio::pairedEnd(files[i], files[i+1])) {
            instruction =
                    "bash -c \"python third-party/spades/bin/metaspades.py -1 " + files[i] + " -2 " + files[i + 1] +
                    " -o " + output_file + "\"";
            jump = 2;
        } else {
            if (Bio::isFasta(files[i]))
                instruction = "bash -c \"python third-party/spades/bin/metaspades.py -s " + files[i] + " -o " + output_file +
                              " --only-assembler\"";
            else
                instruction = "bash -c \"python third-party/spades/bin/metaspades.py -s " + files[i] + " -o " + output_file +
                        "\"";
            jump = 1;
        }
        cout << "Instruction: "<<instruction<<endl;
        System::execute(instruction);
    }
    return output_dirs;
}

vector<Composition> getCompositionsContigs(vector<string> paths)
{
    vector<Composition> output;
    for (auto f:paths)
    {
        string contig = "";
        cout << "File: "<<f<<endl;
        ifstream contigsFile(f,ios::out | ios::app | ios::binary);
        if (contigsFile.is_open())
        {
            string line;
            while(getline(contigsFile, line))
            {
                if (line[0] != '>')
                    contig += line;
                else
                {
                    cout << "Full contig: "<<contig<<endl;
                    cout << "Next Line: "<<line<<endl;
                    exit(1);
                    output.push_back(NaiveComposition(contig));
                    contig = "";
                }
            }
        }
    }
    return output;
}

void buildStructures(string file, sdsl::lcp_bitcompressed<> & lcp,
        sdsl::csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> fm_index)
{
    sdsl::construct(fm_index, file, 1);
    cout << "CSA completed"<<endl;
    sdsl::construct(lcp, file, 1);
    cout << "LCP completed" << endl;
}

size_t buildFullContigFile(vector<string> path, bit_vector & ds)
{
    /*
     * Open the created files and concat contigs
     */
    size_t separator_pos = 0, num_contig = 0;
    bool skip = false;
    string name_file = "contigs.fasta";
    if (remove(name_file) != 0)
        cout << "The file does not exist yet"<<endl;
    else
        cout << "Replacing "<<name_file<<" content"<<endl;
    ofstream singleFile(name_file, ios::app);
    for (auto f: path)
    {
        cout << "File: "<<f<<endl;
        if (separator_pos > 0)
            ds[num_contig] = 1;
        ifstream  contigsFile (f, ios::out | ios::app | ios::binary);
        if (contigsFile.is_open() and singleFile.is_open())
        {
            string line;
            while ( getline (contigsFile,line) )
            {
                if (line[0] != '>' && line[0] != '@')
                {
                    if (!skip)
                    {
                        singleFile << line <<endl;
                        separator_pos += line.length();
                    }
                }else{
                    if (line[0] == '@')
                        skip = false;
                    if (line[0] == '+')
                        skip = true;
                    /*if ((separator_pos > 0) && ((line[0] == '@') || (line[0] == '>')))
                    {
                        singleFile << '>'+to_string(num_contig)<<endl;
                        cs[separator_pos] = 1;
                        separator_pos++;
                        num_contig++;
                    }*/
                    if ((line[0] == '@') || (line[0] == '>'))
                    {
                        singleFile << '>'+to_string(num_contig++)<<endl;
                    }
                }
            }
        }
        contigsFile.close();
    }
    ds[num_contig] = 1;
    singleFile.close();
    return num_contig;
}

void processHybridClusters(string && directory, vector<bool> hybridClusters)
{
    set<string> cluster_files = System::getAllFaFqFiles(directory);
    vector<string> cluster_files_vect;
    cluster_files_vect.insert(cluster_files_vect.end(), cluster_files.begin(), cluster_files.end());
    for (size_t i = 0; i < cluster_files_vect.size();i++)
    {
        if (hybridClusters[i])
        {
            cout << "File: " << cluster_files_vect[i] << endl;
            SequenceContainer sc;
            sc.loadFromFile(cluster_files_vect[i], false);
        }
    }
}

int main(int argc, char** argv)
{
    /*
     * Separator = Space
     * Each environment should be in different directories
     * For each especific sample just two PE files
     */
    cout << "Starting:"<<endl;
    string binner = "metabat";
    vector<string> preffixes, directories,output, environment, f_environments, reads;
    bool clustering = false;
    parse_args(argc, argv,&directories,&environment,&output,&preffixes,&clustering,&binner);
    vector<set<string>> ordered_files;
    for (auto d: directories)
    {
        cout << "Directory: "<<d<<endl;
        ordered_files.push_back(System::getAllFaFqFiles(d,"",true));
    }
    vector<vector<string>> environmental_samples;
    for (auto s: ordered_files)
        environmental_samples.push_back(vector<string>(s.begin(),s.end()));
    /*
     * Envirmental Directories ordered
     */
    for (auto s: environmental_samples)
        for (auto d: s) {
            cout << "File: " << d << endl;
            reads.push_back(d);
        }
    /*
     * MetaSpades Execution + build full chain
     */
    vector<string> t_assembly, t_assembly2, assembly_files;
    if (!output.empty())
    {
        for (size_t i = 0; i < environmental_samples.size(); ++i)
        {
            t_assembly2 = spades_execution(environmental_samples[i], output[i]);
            for (auto s: t_assembly2)
                f_environments.push_back(environment[i]);
            t_assembly.insert(t_assembly.end(), t_assembly2.begin(), t_assembly2.end());
        }
    }else
    {
        for (size_t i = 0; i < preffixes.size(); ++i)
        {
            bool valid = true;
            size_t j = 0;
            while ( valid )
            {
                if (System::exist(preffixes[i]+to_string(j))) {
                    t_assembly.push_back(preffixes[i] + to_string(j));
                    f_environments.push_back(environment[i]);
                }else
                    valid = false;
                j += 2;
            }
        }
    }
    for (size_t i = 0; i < t_assembly.size(); ++i)
    {
        cout << "File (meta): " << t_assembly[i]+"/contigs.fasta" << endl;
        cout << "Environment: "<<f_environments[i]<<endl;
        assembly_files.push_back(t_assembly[i]+"/contigs.fasta");
    }
    /*
     * Genomic Composition
     */
    //vector<Composition> compositions = getCompositionsContigs(assembly_files);
    bit_vector dataset_separator(pow(2, 30), 0); /*BitVectors*/
    buildFullContigFile(assembly_files, dataset_separator);
    /*
     * Rank/Select supporting
     */
    rank_support_v<1> rs_ds(&dataset_separator);
    select_support_mcl<1> ss_ds(&dataset_separator);
    cout << "Number of datasets: "<<rs_ds(pow(2,30))<<endl;
    for (size_t i = 1; i < rs_ds(pow(2,30))+1; ++i)
    {
        (i > 1 )?cout << "Number of contigs of sample "<<i<<": "<<ss_ds(i)-ss_ds(i-1)<<endl
            :cout << "Number of contigs of sample "<<i<<": "<<ss_ds(i)<<endl;
    }
    /*
     * Clustering
     */
    ClusterContainer cc(reads,f_environments,"output/",binner, clustering, rs_ds);
    vector<vector<int>> clusters = cc.getClusters();
    int num_cluster = 0;
    for (auto v:clusters)
    {
        cout << "Number of cluster: "<<num_cluster++<<endl;
        for (auto c:v)
            cout << c<<",";
        cout << endl;
    }
    /*
     * Classify clusters
     */
    cc.classifyClusters();
    //CompactDataStructures
    /*csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> fm_index;
    sdsl::lcp_bitcompressed<> lcp;
    buildStructures("toLCPfile.txt",lcp,fm_index);*/
}
