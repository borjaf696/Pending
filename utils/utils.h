#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <boost/config.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/detail/config_file.hpp>
#include <boost/program_options/parsers.hpp>
#include "sdsl/construct.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/csa_bitcompressed.hpp"
#include "sdsl/lcp_bitcompressed.hpp"
#include "OptionPrinter.hpp"

using namespace std;

template <class t_array>
void print_array(const t_array &array, string name = ""){
    for(uint64_t i = 0; i < array.size(); ++i){
        cout << name << "[" << i << "]=" << array[i] << endl;
    }
}

template <class t_array>
void check_lcp (const t_array &array, size_t threshold)
{
    for (size_t i = 0; i < array.size(); ++i)
    {
        if (array[i] < threshold)
            cout << "[" << i << "]=" <<array[i] <<endl;
    }
}

/*
 * System operations
 */
struct System{
    /*
     * System utils
     */
    static bool exist(string && s)
    {
        struct stat buffer;
        return (stat(s.c_str(), &buffer) == 0);
    }
    static set<string> getAllFaFqFiles(string path, string chain = "", bool recursive = false)
    {
        namespace fs = boost::filesystem;
        set<string> output;
        fs::path fs_path(path);
        if (fs::is_regular_file(fs_path))
            output.emplace(path);
        else
        {
            for (auto & p : fs::directory_iterator(path))
            {
                if (fs::is_regular_file(p))
                {
                    bool add = true;
                    if (chain != "" && (chain != p.path().string()))
                        add = false;
                    if (add)
                    {
                        ostringstream oss;
                        oss << p;
                        string converted_path = oss.str().substr(1, oss.str().size() - 2);
                        string extension = converted_path.substr(converted_path.rfind('.'));
                        if (extension == ".fastq" || extension == ".fq" || extension == ".fasta" || extension == ".fa")
                            output.emplace(converted_path);
                    }
                }else if (recursive)
                {
                    set<string> files = getAllFaFqFiles(p.path().string(),chain,recursive);
                    output.insert(files.begin(),files.end());
                }
            }
        }
        return output;
    }

    static void execute(std::string instruction)
    {
        if (system(instruction.c_str()))
        {
            cout << "Fail on: "<<instruction<<"\n";
            exit(1);
        }
    }

    static std::string appendFiles(std::vector<string> files, std::string newFile)
    {
        std::string instruction = "cat ";
        for (auto s: files)
            instruction += s+" ";
        instruction += ">"+newFile;
        if (system(instruction.c_str())){
            cout << "Problem executing: "<<instruction<<"\n";
            cout << "Exiting\n";
            exit(1);
        }
        return newFile;
    }
    template<typename T>
    static void write_histogram(std::unordered_map<T,std::vector<size_t>> histogram, std::string file_name)
    {
        ofstream file;
        file.open(file_name);
        for (auto k: histogram)
        {
            file << (k.first.str()+"\n");
            for (auto t: k.second)
                file << (std::to_string(t)+" ");
            file << endl;
        }
        file.close();
    }
};
/*
 * Basic Operations
 */
struct Basics
{
    float mean(vector<int> contigsLengthVector)
    {
        int count = 0;
        for (auto v:contigsLengthVector)
            count += v;
        return count / contigsLengthVector.size();
    }

    float standardDeviation(vector<int> contigsLengthVector, float mean)
    {
        float var = 0;
        for (auto v: contigsLengthVector)
            var += (v-mean)*(v-mean);
        var /= contigsLengthVector.size();
        return sqrt(var);
    }
};