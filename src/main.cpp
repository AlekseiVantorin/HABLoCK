#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

int main() {
    std::string data = "ACAGTCGACTGCTAGCTCGTACGTCGATCGACCAAAAAACATATCATCTACGGCAGACTACTGCTGACTGCATGCATCGATCGATCGATCGCTAGCATCGATCA";

    std::vector<std::string> kmers;
    for (int i = 0; i < data.size() - 27; i++) {
        std::string kmer;
        for (int j = 0; j < 27; j++) {
            kmer += data[i + j];
        }
        kmers.push_back(kmer);
    }
    std::sort(kmers.begin(), kmers.end());

    for (auto kmer : kmers) {
        std::cout << kmer << '\n';
    }
    return 0;
}
