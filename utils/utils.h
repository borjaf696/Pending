#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
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
    static std::vector<string> getAllFaFqFiles(std::string path)
    {
        namespace fs = boost::filesystem;
        std::vector<string> output;
        fs::path fs_path(path);
        if (fs::is_regular_file(fs_path))
            output.push_back(path);
        else
        {
            for (auto & p : fs::directory_iterator(path))
            {
                std::ostringstream oss;
                oss << p;
                std::string converted_path = oss.str().substr(1, oss.str().size() - 2);
                output.push_back(converted_path);
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
