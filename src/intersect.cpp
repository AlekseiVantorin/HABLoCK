#include <iostream>
#include <fstream>
#include <string>

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

int main(int argc, char* argv[]) {
    // Check if the number of arguments is correct
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " paternal_dump maternal_dump child_dump paternal_origin maternal_origin" << std::endl;
        return 1;
    }

    const std::string paternal_dump = argv[1];
    const std::string maternal_dump = argv[2];
    const std::string child_dump = argv[3];
    const std::string paternal_origin = argv[4];
    const std::string maternal_origin = argv[5];

    // Open files for reading
    std::ifstream paternalFile(paternal_dump);
    std::ifstream maternalFile(maternal_dump);
    std::ifstream childFile(child_dump);

    // Open files for writing
    std::ofstream paternalOriginFile(paternal_origin);
    std::ofstream maternalOriginFile(maternal_origin);

    // Check if files opened successfully
    if (!paternalFile.is_open() || !maternalFile.is_open() || !childFile.is_open() ||
        !paternalOriginFile.is_open() || !maternalOriginFile.is_open()) {
        std::cerr << "Failed to open one or more files." << std::endl;
        return 1;
    }

    std::string kmer_pat;
    std::string kmer_mat;
    std::string kmer_child;
    int number_pat = 0;
    int number_mat = 0;
    int number_child = 0;

    // Read the initial kmer and number from each file
    paternalFile >> kmer_pat >> number_pat;
    maternalFile >> kmer_mat >> number_mat;
    childFile >> kmer_child >> number_child;

    std::cout << "============= counting =============\n";

    // Read child until we have kmers in it
    while (!childFile.eof()) {
        // Mother catches the child up
        while (kmer_child > kmer_mat && !maternalFile.eof()) {
            maternalFile >> kmer_mat >> number_mat;
        }

        // Father catches the child up
        while (kmer_child > kmer_pat && !paternalFile.eof()) {
            paternalFile >> kmer_pat >> number_pat;
        }

        // Skip if we have kmer in both parents
        if (kmer_child == kmer_mat && kmer_child == kmer_pat) {
            if (!paternalFile.eof()) {
                paternalFile >> kmer_pat >> number_pat;
            }
            if (!maternalFile.eof()) {
                maternalFile >> kmer_mat >> number_mat;
            }
            if (!childFile.eof()) {
                childFile >> kmer_child >> number_child;
            }
            continue;
        }

        // Kmer has maternal origin
        if (kmer_child == kmer_mat) {
            maternalOriginFile << kmer_mat << '\n';
            if (!maternalFile.eof()) {
                maternalFile >> kmer_mat >> number_mat;
            }
        }

        // Kmer has paternal origin
        if (kmer_child == kmer_pat) {
            paternalOriginFile << kmer_pat << '\n';
            if (!paternalFile.eof()) {
                paternalFile >> kmer_pat >> number_pat;
            }
        }

        // Go to next
        childFile >> kmer_child >> number_child;
    }

    // Close all files
    paternalFile.close();
    maternalFile.close();
    childFile.close();
    paternalOriginFile.close();
    maternalOriginFile.close();

    return 0;
}
