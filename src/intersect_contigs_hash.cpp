#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include "contigs.h"

void progress(int step) {
    std::cout << "[";
    for (int i = 0; i < step; i++) {
        std::cout << "#";
    }
    for (int i = step; i < 20; i++) {
        std::cout << ".";
    }
    std::cout << '\n';
}

int64_t PRIME = 1000000007;
constexpr int64_t base = 31;

int64_t calculateHashString(const std::string &str) {
    int64_t hashValue = 0;
    for (char c : str) {
        int data;
        switch (c) {
            case 'A':
                data = 0;
                break;
            case 'C':
                data = 1;
                break;
            case 'G':
                data = 2;
                break;
            default:
                data = 3;
        }

        hashValue = (hashValue * base + data) % PRIME;
    }
    return hashValue;
}

int main(int argc, char* argv[]) {
    // Check if the number of arguments is correct
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " paternal_dump maternal_dump child_rep paternal_origin_rep maternal_origin_rep" << std::endl;
        return 1;
    }

    const std::string paternal_dump = argv[1];
    const std::string maternal_dump = argv[2];
    const std::string child_rep = argv[3];
    const std::string paternal_origin_rep = argv[4];
    const std::string maternal_origin_rep = argv[5];

    // Open files for reading
    std::ifstream paternalFile(paternal_dump);
    std::ifstream maternalFile(maternal_dump);

    // Check if files opened successfully
    if (!paternalFile.is_open() || !maternalFile.is_open()) {
        std::cerr << "Failed to open one or more files." << std::endl;
        return 1;
    }

    std::unordered_map<unsigned long long, int> maternalKmerMap;
    std::cout << "Reading data with maternal origin...\n";

    std::string kmer;
    int number = 1;
    while (maternalFile >> kmer >> number) {
        unsigned long long hashkmer = calculateHashString(kmer);
        maternalKmerMap[hashkmer] = number;
    }

    std::unordered_map<unsigned long long, int> paternalKmerMap;
    std::cout << "Reading data with paternal origin...\n";

    while (paternalFile >> kmer >> number) {
        unsigned long long hashkmer = calculateHashString(kmer);
        paternalKmerMap[hashkmer] = number;
    }

    paternalFile.close();
    maternalFile.close();

    std::cout << "============= counting =============\n";

    for (int i = 0; i < contigs_HG06808.size(); i++) {
        std::ifstream child_contig(child_rep + "HG06808_" + contigs_HG06808[i] + "_sorted.txt");
        if (!child_contig.is_open()) {
            std::cerr << "Failed to open child contig number " << i << ' ' << contigs[i] << '\n';
            return 1;
        }

        std::ofstream paternal_origin(paternal_origin_rep + "HG06808_" + contigs_HG06808[i] + "_pan010.txt");
        if (!paternal_origin.is_open()) {
            std::cerr << "Failed to open paternal origin number " << i << ' ' << contigs[i] << '\n';
            return 1;
        }

        std::ofstream maternal_origin(maternal_origin_rep + "HG06808_" + contigs_HG06808[i] + "_pan011.txt");
        if (!maternal_origin.is_open()) {
            std::cerr << "Failed to open maternal origin number " << i << ' ' << contigs[i] << '\n';
            return 1;
        }

        // Read child until we have kmers in it
        while (!child_contig.eof()) {
            // Go to next
            child_contig >> kmer >> number;

            unsigned long long hashkmer = calculateHashString(kmer);

            if (maternalKmerMap.find(hashkmer) != maternalKmerMap.end() && paternalKmerMap.find(hashkmer) != paternalKmerMap.end()) {
                continue;
            }

            // Mother catches the child up
            if (paternalKmerMap.find(hashkmer) != paternalKmerMap.end()) {
                paternal_origin << kmer << '\n';
            }
            if (maternalKmerMap.find(hashkmer) != maternalKmerMap.end()) {
                maternal_origin << kmer << '\n';
            }
        }

        child_contig.close();
        maternal_origin.close();
        paternal_origin.close();
    }

    return 0;
}
