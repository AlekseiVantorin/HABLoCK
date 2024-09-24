#include <iostream>
#include <fstream>
#include <string>

void findPotentialCrossingOvers(const std::string &input, const std::string &output) {

    std::cout << "============= counting =============\n";

    // Open files
    std::ifstream childFile(input);
    std::ofstream outFile(output);

    char c;
    int64_t position = 0;
    int cur_pos = 0;
    int border = 50000;

    int maternal = 0;
    int paternal = 0;
    int unknown = 0;

    int number_paternal = 0;
    int number_maternal = 0;

    std::string prev_contig;
    childFile >> prev_contig;

    // Read file with origins
    while (childFile.get(c)) {
        // Count origins
        switch (c) {
            case 'm':
                maternal++;
                break;
            case 'p':
                paternal++;
                break;
            case 'u':
                unknown++;
                break;
            default:
                std::string new_contig;
                childFile.get(c);
                while (c != '\n') {
                    new_contig += c;
                    childFile.get(c);
                }

                outFile << prev_contig << " - ";
                prev_contig = new_contig;

                if (maternal > paternal) {
                    outFile << "m - " << maternal << " " << paternal << " " << unknown << " - ";

                    outFile << maternal * 100.0 / border << " " << paternal * 100.0 / border << " - ";

                    if (paternal * 99 > maternal) {
                        outFile << "mistake in counting\n";
                    } else {
                        outFile << "ok\n";
                    }
                    number_maternal++;
                } else if (maternal < paternal) {
                    outFile << "p - " << maternal << " " << paternal << " " << unknown << " - ";

                    outFile << maternal * 100.0 / border << " " << paternal * 100.0 / border << " - ";

                    if (maternal * 99 > paternal) {
                        outFile << "mistake in counting\n";
                    } else {
                        outFile << "ok\n";
                    }
                    number_paternal++;
                } else {
                    outFile << "a - " << maternal << " " << paternal << " " << unknown << " - ";

                    outFile << maternal * 100.0 / border << " " << paternal * 100.0 / border << " - ";

                    if ((maternal * 99 > paternal) || (paternal * 99 > maternal)) {
                        outFile << "mistake in counting\n";
                    } else {
                        outFile << "ok\n";
                    }
                }

                maternal = 0;
                paternal = 0;
                unknown = 0;
        }

        cur_pos++;
        position++;
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
