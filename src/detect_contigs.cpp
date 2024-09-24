#include <iostream>
#include <fstream>
#include <string>

void findPotentialCrossingOvers(const std::string &input, const std::string &output) {

    std::cout << "============= counting =============\n";

    // Open files
    std::ifstream childFile(input);
    std::ofstream outFile(output);

    std::string contig;
    char c;
    bool skip_enter = false;
    while (childFile.get(c)) {
//        if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N' || c == '>') {
        if (c == 'u' || c == 'm' || c == 'p') {
            continue;
        }
        if (c == '\n') {
            if (skip_enter) {
                skip_enter = false;
                continue;
            } else {
                skip_enter = true;
            }
        }
        outFile << c;
    }
}

int main(int argc, char* argv[]) {

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input_file.txt output_file.txt" << std::endl;
        return 1;
    }

    std::string input = argv[1];
    std::string output = argv[2];

    findPotentialCrossingOvers(input, output);
    return 0;
}
