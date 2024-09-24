#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {
    // Check if the number of arguments is correct
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " counted_file" << std::endl;
        return 1;
    }

    const std::string filefile = argv[1];
    std::ifstream file(filefile);

    std::string kmer;
    std::string prev_kmer;
    int number;

    file >> kmer >> number;
    prev_kmer = kmer;

    int answer = 0;
    while (!file.eof()) {
        file >> kmer >> number;
        if (kmer < prev_kmer) {
            answer ++;
        }
        prev_kmer = kmer;
    }

    std::cout << answer;

    return 0;
}
