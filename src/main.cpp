#include "../operations/GenomicComposition.hpp"
#include "../utils/utils.h"

namespace po = boost::program_options;
namespace
{
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
}
using namespace sdsl;
using namespace std;

void parse_args(int argc, char ** argv, vector<string> * directories, vector<string> * environment,
        vector<string> * output, vector<string> * preffixes)
{
    po::options_description des("Options");
    des.add_options()
            ("help,h","Help message")
            ("dir,d",po::value<vector<string>>(directories)->multitoken(),"Reads directories one per environment")
            ("environment,e",po::value<vector<string>>(environment)->multitoken(),"Environment which belong each sample")
            ("output,o",po::value<vector<string>>(output)->multitoken(),"Preffixes Output directories")
            ("preffix,p",po::value<vector<string>>(preffixes)->multitoken(),"Preffixes where metaspades has been run - environment 1");
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
    for (size_t i = 0; i < files.size(); i+=2)
    {
        string output_file = output+to_string(i);
        output_dirs.push_back(output_file);
        string instruction = "bash -c \"python third-party/spades/bin/metaspades.py -1 "+files[i]+" -2 "+files[i+1]+" -o "+output_file+"\"";
        cout << "Instruction: "<<instruction<<endl;
        System::execute(instruction);
    }
    return output_dirs;
}

void maxbin_execution(vector<string> files, string output_preffix)
{
    string contig_file = "contigs.fasta";
    string instruction = "bash -c \"perl third-party/MaxBin/run_MaxBin.pl -contig "+contig_file+" ";
    int num_read = 1;
    for (auto s:files)
        instruction+= "-read"+to_string(num_read++)+" "+s+" ";
    instruction += "-thread 32 -output "+output_preffix+"\"";
    cout << "Instruction: "<<instruction<<endl;
    System::execute(instruction);
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

void buildFullContigFile(vector<string> path, bit_vector & ds)
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
                        singleFile << '>'+to_string(num_contig)<<endl;
                        num_contig++;
                    }
                }
            }
        }
        contigsFile.close();
    }
    ds[num_contig] = 1;
    singleFile.close();
}

int main(int argc, char** argv)
{
    /*
     * Separator = Space
     * Each environment should be in different directories
     * For each especific sample just two PE files
     */
    cout << "Starting:"<<endl;
    vector<string> preffixes, directories,output, environment, f_environments;
    parse_args(argc, argv,&directories,&environment,&output,&preffixes);
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
        for (auto d: s)
            cout << "Directories: "<<d<<endl;
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
                if (System::exist(preffixes[i]+to_string(j)))
                    t_assembly.push_back(preffixes[i]+to_string(j));
                else
                    valid = false;
                j += 2;
            }
        }
    }
    for (auto t: t_assembly)
    {
        cout << "File (meta): " << t+"/contigs.fasta" << endl;
        assembly_files.push_back(t+"/contigs.fasta");
    }
    /*
     * Genomic Composition
     */
    //vector<Composition> compositions = getCompositionsContigs(assembly_files);
    bit_vector dataset_separator(pow(2, 30), 0); /*BitVectors*/
    buildFullContigFile(assembly_files, dataset_separator);
    /*
     * MaxBin Execution
     */
    maxbin_execution(assembly_files, "output");
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

    //CompactDataStructures
    /*csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> fm_index;
    sdsl::lcp_bitcompressed<> lcp;
    buildStructures("toLCPfile.txt",lcp,fm_index);*/
}
