#include <iostream>
#include <fstream>
#include <string>
#include "contigs.h"

void catch_next_contig(std::ifstream &origins) {
    char c;
    while (origins.get(c) && c != '\n');
}

struct Block {
    char origin;
    int start;
    int end;

    Block(char orig, int st, int en) : origin(orig), start(st), end(en) {}
};

void detectHaploblocks(const std::string &origins, const std::string &answer, const std::string &output) {

    std::cout << "============= counting =============\n";

    // Open files
    std::ifstream originsFile(origins);
    std::ifstream answerFile(answer);
    std::ofstream outFile(output);

    if (!originsFile.is_open()) {
        std::cout << "Problem with origin file\n";
    }
    if (!answerFile.is_open()) {
        std::cout << "Problem with answer file\n";
    }

    char c;
    for (int j = 0; j < contigs_HG06808.size(); j++) {
        //std::cout << "taking a look on contig " << contigs[i] << '\n';
        if (originsFile.eof()) {
            return;
        }
        //std::string name = contigs_with_size[i].first;
        std::string name = contigs_HG06808[j];
        //long long size = contigs_with_size[i].second;
        long long size = 0;

        // look at answer according to this contig
        // if we have big amount of parents -> we want to look for haploblocks
        // otherwise we are not interested

        std::string contig_name;
        std::string contig_name_from_origins;
        int maternal = 0;
        int paternal = 0;
        int unknown = 0;

        if (j != 0) {
            answerFile >> contig_name;
            //originsFile >> contig_name_from_origins;
            //while (contig_name == "") {
            //    answerFile >> contig_name;
           // }
            //while (contig_name_from_origins == "") {
            //    originsFile >> contig_name_from_origins;
            //}
            //contig_name.erase(0, 1);
            //contig_name_from_origins.erase(0, 1);
            //if (contig_name != contig_name_from_origins) {
            //    std::cout << "Problem with reading contig " << contig_name << ' ' << contig_name_from_origins << '\n';
            //    catch_next_contig(originsFile);
            //    answerFile >> maternal;
            //    answerFile >> paternal;
            //    answerFile >> unknown;
            //    contig++;
            //}
        }
        answerFile >> maternal;
        answerFile >> paternal;
        answerFile >> unknown;

        outFile << name << ":\n";
        std::cout << name << ":\n";

        long long iterator = 0;

//        size = unknown + maternal + paternal;
        size = maternal + paternal;

        //std::cout << name << " - " << contig_name << " - " << " - " << maternal << " - " << paternal << " - " << unknown << '\n';

        if (maternal < 20 && paternal < 20) {
            while (c != '>') {
                originsFile.get(c);
            }
            while (c != '\n') {
                originsFile.get(c);
            }
            outFile << "    low amount kmer with origin in contig\n";
            std::cout << "    low amount kmer with origin in contig\n";
            continue;
        }
        // too low maternal and paternal k-mers to have any suggestions

        if (maternal * 100 < paternal || paternal * 100 < maternal) {
            while (c != '>') {
                originsFile.get(c);
            }
            while (c != '\n') {
                originsFile.get(c);
            }
            outFile << "    full contig is ";
            std::cout << "    full contig is ";
            if (maternal > paternal) {
                outFile << "m\n";
                std::cout << "m\n";
            } else {
                outFile << "p\n";
                std::cout << "p\n";
            }
            continue;
        }
        // that means that the whole contig is fully maternal or fully paternal

        int maternal_counted = 0;
        int paternal_counted = 0;
        int block_size = 20;
        int blocks_num = (size - 1) / block_size + 1;
        //int blocks_num = 1;
        std::vector<std::pair<int, int>> origins_in_block(blocks_num);
        // first is paternal, second is maternal
        int iterator_blocks = 0;
        int iterator_chars = 0;

        while (originsFile.get(c) && c != '>') {
            if (c == 'p') {
                origins_in_block[iterator_blocks].first++;
                paternal_counted++;
            } else if (c == 'm') {
                origins_in_block[iterator_blocks].second++;
                maternal_counted++;
            }

            if (iterator_chars == block_size) {
                iterator_chars = 0;
                iterator_blocks++;
                if (iterator_blocks > origins_in_block.size()) {
                    origins_in_block.resize(iterator_blocks);
                }
            }

            if (c != 'u')
            iterator_chars++;
        }

        blocks_num = origins_in_block.size();
        std::vector<char> block_originality(blocks_num);
        for (int i = 0; i < blocks_num; i++) {
            if (origins_in_block[i].first > 5 && origins_in_block[i].second > 5) {
                block_originality[i] = 'c';
            } else if (origins_in_block[i].first > origins_in_block[i].second) {
                block_originality[i] = 'p';
            } else if (origins_in_block[i].first < origins_in_block[i].second) {
                block_originality[i] = 'm';
            } else {
                block_originality[i] = 'a';
            }
        }

        for (int i = 1; i < blocks_num - 1; i++) {
            if (block_originality[i] == 'a' || block_originality[i] == 'c') {
                if (block_originality[i - 1] == block_originality[i + 1] && block_originality[i - 1] != 'a' && block_originality[i - 1] != 'c') {
                    block_originality[i] = block_originality[i - 1];
                } else if (block_originality[i - 1] == block_originality[i + 2] && block_originality[i - 1] != 'a' && block_originality[i - 1] != 'c' && (block_originality[i + 1] == 'a'
                    || block_originality[i + 1] == 'c')) {
                    block_originality[i] = block_originality[i - 1];
                    block_originality[i + 1] = block_originality[i - 1];
                }
            }
        }

        int mat = 0;
        int pat = 0;
        int cro = 0;
        int unk = 0;
        for (char i : block_originality) {
            switch (i) {
                case 'm':
                    mat++;
                    break;
                case 'p':
                    pat++;
                    break;
                case 'c':
                    cro++;
                    break;
                default:
                    unk++;
            }
        }
        std::cout << mat << ' ' << pat << ' ' << cro << ' ' << unk << '\n';
        //std::cout << maternal_counted << ' ' << paternal_counted << '\n';

        int start_block = 0;
        int end_block = 0;
        char prev_block = block_originality[0];
        std::vector<Block> result;
        for (int i = 1; i < blocks_num; i++) {
            if (block_originality[i] == 'c') {
                bool is_crossingover = false;
                if (block_originality[i - 1] == 'm' && block_originality[i + 1] == 'p') {
                    is_crossingover = true;
                } else if (block_originality[i - 1] == 'p' && block_originality[i + 1] == 'm') {
                    is_crossingover = true;
                }

                //if (i == 1 || block_originality[i - 1] != block_originality[i - 2]) {
                //    is_crossingover = false;
                //}

                //if (i == blocks_num - 1 || block_originality[i + 1] != block_originality[i + 2]) {
                //    is_crossingover = false;
                //}

                if (is_crossingover) {
                    if (block_originality[i] == 'p') {
                        outFile << "    Crossingover between positions [" << i * block_size << " and " << (i + 1) * block_size << "] - p -> c -> m\n";
                        std::cout << "    Crossingover between positions [" << i * block_size << " and " << (i + 1) * block_size << "] - p -> c -> m\n";
                    } else {
                        outFile << "    Crossingover between positions [" << i * block_size << " and " << (i + 1) * block_size << "] - m -> c -> p\n";
                        std::cout << "    Crossingover between positions [" << i * block_size << " and " << (i + 1) * block_size << "] - m -> c -> p\n";
                    }
                }
            }

            if (block_originality[i] == prev_block) {
                end_block = i;
            } else {
                if (block_originality[i] == 'm') {
                    if (block_originality[i - 1] == 'p') {
                        outFile << "    Crossingover between positions [" << i * block_size << " and " << (i + 1) * block_size << "] - p -> m\n";
                        std::cout << "    Crossingover between positions [" << i * block_size << " and " << (i + 1) * block_size << "] - p -> m\n";
                    }
                }
                if (block_originality[i] == 'p') {
                    if (block_originality[i - 1] == 'm') {
                        outFile << "    Crossingover between positions [" << i * block_size << " and " << (i + 1) * block_size << "] - m -> p\n";
                        std::cout << "    Crossingover between positions [" << i * block_size << " and " << (i + 1) * block_size << "] - m -> p\n";
                    }
                }
                if (start_block == end_block) {
                    start_block = i;
                    end_block = i;
                    prev_block = block_originality[i];
                } else {
                    if (block_originality[i] != 'a') {
                        result.push_back({prev_block, start_block, end_block});
                    }
                    start_block = i;
                    end_block = i;
                    prev_block = block_originality[i];
                }
            }
        }

//        if (mat == 0) {
//            outFile << "    Origin of contig is paternal, number of maternal blocks is 0\n";
//            std::cout << "    Origin of contig is paternal, number of maternal blocks is 0\n";
//            continue;
//        }
//
//        if (pat == 0) {
//            outFile << "    Origin of contig is maternal, number of paternal blocks is 0\n";
//            std::cout << "    Origin of contig is maternal, number of paternal blocks is 0\n";
//            continue;
//        }

        for (auto res : result) {
            if (res.end - res.start > 1 && res.origin != 'a' && res.origin != 'c') {
                outFile << "    Origin of parent " << res.origin << " between positions [" << res.start * block_size << " and " << (res.end + 1) * block_size << "]\n";
                //std::cout << "    Origin of parent " << res.origin << " between positions [" << res.start * block_size << " and " << (res.end + 1) * block_size << "]\n";
            }
        }

        //if (originsFile.eof()) {
        //    return;
        //}

        while (c != '\n') {
            originsFile.get(c);
        }

        if (originsFile.eof()) {
            return;
        }
    }
}

int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " origins_file answer_file output_file" << std::endl;
        return 1;
    }

    std::string origins = argv[1];
    std::string answer = argv[2];
    std::string output = argv[3];

    detectHaploblocks(origins, answer, output);
    return 0;
}
