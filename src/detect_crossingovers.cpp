#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#define AREA 1000000

// Function to find neighboring numbers and store the first and last number of each sequence
void findNeighbors(const std::string& inputFile, const std::string& outputFile) {
    // Open the input file
    std::ifstream inFile(inputFile);
    if (!inFile.is_open()) {
        std::cerr << "Error: Unable to open " << inputFile << std::endl;
        return;
    }

    std::ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open " << outputFile << std::endl;
        return;
    }

    std::cout << "======================= counting ======================\n";

    // Read numbers from the file and process
    unsigned long long num;
    inFile >> num; // Read the first number
    unsigned long long start = num; // Start of the current sequence
    unsigned long long prev = num; // Previous number
    //outFile << num << '\n';
    while (inFile >> num) {
        // std::cout << num << ' ' << prev << ' ' << start << '\n';
        if (num - prev > AREA) { // Check if current number is outside the area
            if (prev - start <= AREA && prev != start) { // Check if a sequence of neighbors was found
                outFile << start << '\n';
                outFile << prev << '\n';
            }
            start = num; // Start a new sequence
        }
        prev = num; // Update previous number
    }
    inFile.close();

    // Check if the last sequence is a sequence of neighbors
    if (prev - start >= AREA) {
        outFile << start << '\n';
        outFile << prev << '\n';
    }

    std::cout << "====================== finished ========================\n";
}

int main(int argc, char* argv[]) {
    // Check if the number of arguments is correct
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input_file output_file" << std::endl;
        return 1;
    }

    const std::string inputFile = argv[1];
    const std::string outputFile = argv[2];

    // Find neighboring numbers and store the first and last number of each sequence
    findNeighbors(inputFile, outputFile);

    return 0;
}
