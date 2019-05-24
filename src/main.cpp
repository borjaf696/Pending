#include "../utils/utils.h"

namespace po = boost::program_options;
namespace
{
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
}
using namespace sdsl;
using namespace std;

void parse_args(int argc, char ** argv, string * dir1, string * dir2)
{
    po::options_description des("Options");
    des.add_options()
            ("help,h","Help message")
            ("dir1,d1",po::value<std::string>(dir1)->required(),"Thermophile reads directory")
            ("dir2,d2",po::value<std::string>(dir2)->required(),"Non-thermophile reads directory");
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

void spades_execution(string path)
{
	string instruction = "bash -c python ../thrid-party/spades.py --meta ";	
}

void buildStructures(string file)
{
    sdsl::csa_bitcompressed<> csa;
    sdsl::construct(csa, file, 1);
    cout << "CSA completed"<<endl;
    sdsl::lcp_bitcompressed<> lcp;
    sdsl::construct(lcp, file, 1);
    cout << "LCP completed" << endl;
    print_array(lcp, "lcp");
}

void buildFullContigFile(vector<string> path, bit_vector & cs, bit_vector & ds)
{
    /*
     * Open the created files and concat contigs
     */
    size_t separator_pos = 0;
    string name_file = "toLCPfile.txt";
    if (remove(name_file) != 0)
        cout << "The file does not exist yet"<<endl;
    else
        cout << "Replacing "<<name_file<<" content"<<endl;
    ofstream singleFile(name_file, ios::app);
    for (auto f: path)
    {
        cout << "File: "<<f<<endl;
        if (separator_pos > 0)
            ds[separator_pos] = 1;
        ifstream  contigsFile (f, ios::out | ios::app | ios::binary);
        if (contigsFile.is_open() and singleFile.is_open())
        {
            string line;
             while ( getline (contigsFile,line) )
            {
                if (line[0] != '>')
                {
                    singleFile << line;
                    separator_pos += line.length();
                }else{
                    if (separator_pos > 0)
                    {
                        singleFile << '$';
                        cs[separator_pos] = 1;
                        separator_pos++;
                    }
                }
            }
        }
        contigsFile.close();
    }
    singleFile.close();
}

int main(int argc, char** argv)
{
    /*
     * Pairs of paired_end reads
     * -1 1,2 -2 1,2
     */
    string dir1 = "", dir2 = "";
    parse_args(argc, argv,&dir1, &dir2);
    bit_vector contig_separator(pow(2, 30), 0), dataset_separator(pow(2, 30), 0);
    string path = "";
    vector <string> contigsPath;
    //Ashoc tests
    contigsPath.push_back("datasets/MetagenomicTrial/contigs.fasta");
    contigsPath.push_back("datasets/MetagenomicTrial2/contigs.fasta");
    //BuildFullFile
    buildFullContigFile(contigsPath, contig_separator, dataset_separator);
    rank_support_v<1> rs_cs(&contig_separator), rs_ds(&dataset_separator);
    select_support_mcl<1> ss_cs(&contig_separator), ss_ds(&dataset_separator);
    cout << "Number of datasets: "<<rs_ds(pow(2,30))<<endl;
    for (size_t i = 1; i < 10; ++i)
        cout << "Sequence: "<<i<<" Length: "<<ss_cs(i+1)-ss_cs(i)<<endl;
    //Building LCP
    buildStructures("toLCPfile.txt");
}
