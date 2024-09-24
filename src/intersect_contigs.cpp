#include <iostream>
#include <fstream>
#include <string>
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

    std::vector<std::ifstream> child_contigs(333);
    for (int i = 0; i < 333; i++) {
        child_contigs[i] = std::ifstream (child_rep + "HG06807_" + contigs[i] + "_sorted.txt");
        if (!child_contigs[i].is_open()) {
            std::cerr << "Failed to open child contig number " << i << ' ' << contigs[i] << '\n';
            return 1;
        }
    }

    // Open files for writing
    std::vector<std::ofstream> paternal_origins(333);
    for (int i = 0; i < 333; i++) {
        paternal_origins[i] = std::ofstream (paternal_origin_rep + "HG06807_" + contigs[i] + "_pan010.txt");
        if (!paternal_origins[i].is_open()) {
            std::cerr << "Failed to open paternal origin number " << i << ' ' << contigs[i] << '\n';
            return 1;
        }
    }
    std::vector<std::ofstream> maternal_origins(333);
    for (int i = 0; i < 333; i++) {
        maternal_origins[i] = std::ofstream (maternal_origin_rep + "HG06807_" + contigs[i] + "_pan011.txt");
        if (!maternal_origins[i].is_open()) {
            std::cerr << "Failed to open maternal origin number " << i << ' ' << contigs[i] << '\n';
            return 1;
        }
    }

    // Check if files opened successfully
    if (!paternalFile.is_open() || !maternalFile.is_open()) {
        std::cerr << "Failed to open one or more files." << std::endl;
        return 1;
    }

    std::string kmer_pat;
    std::string kmer_mat;
    std::vector<std::string> kmer_contig(333);

    int number_pat = 0;
    int number_mat = 0;
    int number_child = 0;

    // Read the initial kmer and number from each file
    paternalFile >> kmer_pat >> number_pat;
    maternalFile >> kmer_mat >> number_mat;
    for (int i = 0; i < kmer_contig.size(); i++) {
        child_contigs[i] >> kmer_contig[i] >> number_child;
    }

    std::cout << "============= counting =============\n";

    // Read child until we have kmers in it
    while (!maternalFile.eof()) {
        // Mother catches the child up
        for (int i = 0; i < kmer_contig.size(); i++) {
            while (kmer_contig[i] < kmer_mat && !child_contigs[i].eof()) {
                child_contigs[i] >> kmer_contig[i] >> number_child;
            }
        }

        // Kmer has maternal origin
        for (int i = 0; i < kmer_contig.size(); i++) {
            if (kmer_contig[i] == kmer_mat) {
                maternal_origins[i] << kmer_mat << '\n';
                if (!child_contigs[i].eof()) {
                    child_contigs[i] >> kmer_contig[i] >> number_child;
                }
            }
        }

        // Go to next
        maternalFile >> kmer_mat >> kmer_mat;
    }

    while (!paternalFile.eof()) {
        // Mother catches the child up
        for (int i = 0; i < kmer_contig.size(); i++) {
            while (kmer_contig[i] < kmer_pat && !child_contigs[i].eof()) {
                child_contigs[i] >> kmer_contig[i] >> number_child;
            }
        }

        // Kmer has maternal origin
        for (int i = 0; i < kmer_contig.size(); i++) {
            if (kmer_contig[i] == kmer_pat) {
                paternal_origins[i] << kmer_pat << '\n';
                if (!child_contigs[i].eof()) {
                    child_contigs[i] >> kmer_contig[i] >> number_child;
                }
            }
        }

        // Go to next
        paternalFile >> kmer_pat >> kmer_pat;
    }

    // Close all files
    paternalFile.close();
    maternalFile.close();
    for (int i = 0; i < 333; i++) {
        child_contigs[i].close();
        maternal_origins[i].close();
        paternal_origins[i].close();
    }





// That shouldn't work like that and there must be better solution, but this is life

    std::ifstream secondPaternalFile(paternal_dump);
    std::ifstream secondMaternalFile(maternal_dump);

    for (int i = 0; i < 333; i++) {
        child_contigs[i] = std::ifstream (child_rep + "HG06807_" + contigs[i + 333] + "_sorted.txt");
        if (!child_contigs[i].is_open()) {
            std::cerr << "Failed to open child contig number " << i << ' ' << contigs[i + 333] << '\n';
            return 1;
        }
    }

    // Open files for writing
    for (int i = 0; i < 333; i++) {
        paternal_origins[i] = std::ofstream (paternal_origin_rep + "HG06807_" + contigs[i + 333] + "_pan010.txt");
        if (!paternal_origins[i].is_open()) {
            std::cerr << "Failed to open paternal origin number " << i << ' ' << contigs[i + 333] << '\n';
            return 1;
        }
    }
    for (int i = 0; i < 333; i++) {
        maternal_origins[i] = std::ofstream (maternal_origin_rep + "HG06807_" + contigs[i + 333] + "_pan011.txt");
        if (!maternal_origins[i].is_open()) {
            std::cerr << "Failed to open maternal origin number " << i << ' ' << contigs[i + 333] << '\n';
            return 1;
        }
    }

    // Check if files opened successfully
    if (!secondPaternalFile.is_open() || !secondMaternalFile.is_open()) {
        std::cerr << "Failed to open one or more files." << std::endl;
        return 1;
    }

    // Read the initial kmer and number from each file
    secondPaternalFile >> kmer_pat >> number_pat;
    secondMaternalFile >> kmer_mat >> number_mat;
    for (int i = 0; i < kmer_contig.size(); i++) {
        child_contigs[i] >> kmer_contig[i] >> number_child;
    }

    std::cout << "============= counting =============\n";

    // Read child until we have kmers in it
    while (!secondMaternalFile.eof()) {
        // Mother catches the child up
        for (int i = 0; i < kmer_contig.size(); i++) {
            while (kmer_contig[i] < kmer_mat && !child_contigs[i].eof()) {
                child_contigs[i] >> kmer_contig[i] >> number_child;
            }
        }

        // Kmer has maternal origin
        for (int i = 0; i < kmer_contig.size(); i++) {
            if (kmer_contig[i] == kmer_mat) {
                maternal_origins[i] << kmer_mat << '\n';
                if (!child_contigs[i].eof()) {
                    child_contigs[i] >> kmer_contig[i] >> number_child;
                }
            }
        }

        // Go to next
        secondMaternalFile >> kmer_mat >> kmer_mat;
    }

    while (!secondPaternalFile.eof()) {
        // Mother catches the child up
        for (int i = 0; i < kmer_contig.size(); i++) {
            while (kmer_contig[i] < kmer_pat && !child_contigs[i].eof()) {
                child_contigs[i] >> kmer_contig[i] >> number_child;
            }
        }

        // Kmer has maternal origin
        for (int i = 0; i < kmer_contig.size(); i++) {
            if (kmer_contig[i] == kmer_pat) {
                paternal_origins[i] << kmer_pat << '\n';
                if (!child_contigs[i].eof()) {
                    child_contigs[i] >> kmer_contig[i] >> number_child;
                }
            }
        }

        // Go to next
        secondPaternalFile >> kmer_pat >> kmer_pat;
    }

    // Close all files
    secondMaternalFile.close();
    secondPaternalFile.close();
    for (int i = 0; i < 333; i++) {
        child_contigs[i].close();
        maternal_origins[i].close();
        paternal_origins[i].close();
    }





    // And the last time


    std::ifstream thirdPaternalFile(paternal_dump);
    std::ifstream thirdMaternalFile(maternal_dump);

    child_contigs = std::vector<std::ifstream> (276);
    for (int i = 0; i < 276; i++) {
        child_contigs[i] = std::ifstream (child_rep + "HG06807_" + contigs[i + 666] + "_sorted.txt");
        if (!child_contigs[i].is_open()) {
            std::cerr << "Failed to open child contig number " << i << ' ' << contigs[i + 666] << '\n';
            return 1;
        }
    }

    // Open files for writing
    paternal_origins = std::vector<std::ofstream> (276);
    for (int i = 0; i < 276; i++) {
        paternal_origins[i] = std::ofstream (paternal_origin_rep + "HG06807_" + contigs[i + 666] + "_pan010.txt");
        if (!paternal_origins[i].is_open()) {
            std::cerr << "Failed to open paternal origin number " << i << ' ' << contigs[i + 666] << '\n';
            return 1;
        }
    }
    maternal_origins = std::vector<std::ofstream> (276);
    for (int i = 0; i < 276; i++) {
        maternal_origins[i] = std::ofstream (maternal_origin_rep + "HG06807_" + contigs[i + 666] + "_pan011.txt");
        if (!maternal_origins[i].is_open()) {
            std::cerr << "Failed to open maternal origin number " << i << ' ' << contigs[i + 666] << '\n';
            return 1;
        }
    }

    // Check if files opened successfully
    if (!thirdPaternalFile.is_open() || !thirdMaternalFile.is_open()) {
        std::cerr << "Failed to open one or more files." << std::endl;
        return 1;
    }

    kmer_contig = std::vector<std::string> (276);

    // Read the initial kmer and number from each file
    thirdPaternalFile >> kmer_pat >> number_pat;
    maternalFile >> kmer_mat >> number_mat;
    for (int i = 0; i < kmer_contig.size(); i++) {
        child_contigs[i] >> kmer_contig[i] >> number_child;
    }

    std::cout << "============= counting =============\n";

    // Read child until we have kmers in it
    while (!thirdMaternalFile.eof()) {
        // Mother catches the child up
        for (int i = 0; i < kmer_contig.size(); i++) {
            while (kmer_contig[i] < kmer_mat && !child_contigs[i].eof()) {
                child_contigs[i] >> kmer_contig[i] >> number_child;
            }
        }

        // Kmer has maternal origin
        for (int i = 0; i < kmer_contig.size(); i++) {
            if (kmer_contig[i] == kmer_mat) {
                maternal_origins[i] << kmer_mat << '\n';
                if (!child_contigs[i].eof()) {
                    child_contigs[i] >> kmer_contig[i] >> number_child;
                }
            }
        }

        // Go to next
        thirdMaternalFile >> kmer_mat >> kmer_mat;
    }

    while (!thirdPaternalFile.eof()) {
        // Mother catches the child up
        for (int i = 0; i < kmer_contig.size(); i++) {
            while (kmer_contig[i] < kmer_pat && !child_contigs[i].eof()) {
                child_contigs[i] >> kmer_contig[i] >> number_child;
            }
        }

        // Kmer has maternal origin
        for (int i = 0; i < kmer_contig.size(); i++) {
            if (kmer_contig[i] == kmer_pat) {
                paternal_origins[i] << kmer_pat << '\n';
                if (!child_contigs[i].eof()) {
                    child_contigs[i] >> kmer_contig[i] >> number_child;
                }
            }
        }

        // Go to next
        thirdPaternalFile >> kmer_pat >> kmer_pat;
    }

    // Close all files
    thirdPaternalFile.close();
    thirdMaternalFile.close();
    for (int i = 0; i < 333; i++) {
        child_contigs[i].close();
        maternal_origins[i].close();
        paternal_origins[i].close();
    }


    return 0;
}
