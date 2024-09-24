#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include "contigs.h"

int64_t PRIME = 1000000007;
constexpr int64_t base = 31;

void progress(int step) {
    std::cout << "[";
    for (int i = 0; i < step; i++) {
        std::cout << "#";
    }
    for (int i = step; i < 20; i++) {
        std::cout << ".";
    }
    std::cout << "]\n";
}

struct Node {
    int data;
    Node* next;

    Node(int data) : data(data), next(nullptr) {}
};

// Calculate hash on given Nodes
int64_t calculateHash(Node* head) {
    int64_t hashValue = 0;
    while (head != nullptr) {
        hashValue = (hashValue * base + head->data) % PRIME;
        head = head->next;
    }
    return hashValue;
}

// Some mathematical numbers according to hash recounting
// 31 ^ 26 = 559050727
// for A : 0
// for C : 559050727
// for G : 559050727 * 2 % 1000000007 = 118101447
// for T : 559050727 * 3 % 1000000007 = 677152174
// default: 559050727 * 4 % 1000000007 = 236202894

// -559050727 % 1000000007 = 440949280
// -118101447 % 1000000007 = 881898560
// -677152174 % 1000000007 = 322847833
// -236202894 % 1000000007 = 763797113

// Recalculate given hash with deleting given head and adding given tail
int64_t recalculateHash(unsigned long long cur_hash, Node* head, Node* tail) {
    int64_t hashValue = cur_hash;
    switch (head->data) {
        case 0:
            break;
        case 1:
            hashValue += 440949280;
            hashValue %= PRIME;
            break;
        case 2:
            hashValue += 881898560;
            hashValue %= PRIME;
            break;
        case 4:
            hashValue += 322847833;
            hashValue %= PRIME;
            break;
        default:
            hashValue += 763797113;
            hashValue %= PRIME;
    }
    hashValue = (hashValue * base + tail->data) % PRIME;

    return hashValue;
}

// Calculate hash of given string
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
            case 'T':
                data = 3;
                break;
            default:
                data = 4;
        }

        hashValue = (hashValue * base + data) % PRIME;
    }
    return hashValue;
}

std::pair<Node*, Node*> starting_reading(std::ifstream &childFile) {
    Node* head = nullptr;
    Node* tail = nullptr;

    char c;
    for (int i = 0; i < 26; ++i) {
        childFile.get(c);
        int data = 0;
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
            case 'T':
                data = 3;
                break;
            default:
                data = 4;
        }

        if (head == nullptr) {
            head = new Node(data);
            tail = head;
        } else {
            tail->next = new Node(data);
            tail = tail->next;
        }
    }

    return {head, tail};
}

//std::string nodeToString(Node *head) {
//    std::string ans;
//    while (head != nullptr) {
//        switch (head->data) {
//            case 0:
//                ans += 'A';
//                break;
//            case 1:
//                ans += 'C';
//                break;
//            case 2:
//                ans += 'G';
//                break;
//            default:
//                ans += 'T';
//        }
//        head = head->next;
//    }
//    return ans;
//}

unsigned long long mother = 0;
unsigned long long father = 0;
unsigned long long unknown = 0;

void findOriginInsideContig(std::ifstream &childFile, std::ofstream &outFile, std::ofstream &answerFile, int contig) {
    std::string kmer;
    int number = 1;

    std::string repo = "../../mnt/projects/assembly/project2/crossingover/counted_intersections/HG06808/pan01";

    std::unordered_map<unsigned long long, int> maternalKmerMap;
    std::unordered_map<unsigned long long, int> paternalKmerMap;

    std::unordered_map<std::string, int> maternalKmerMapString;
    std::unordered_map<std::string, int> paternalKmerMapString;

    Node* head = nullptr;
    Node* tail = nullptr;
    char c;

    std::ifstream maternalFile(repo + "1/HG06808_" + contigs_HG06808[contig] + "_pan011.txt");
    std::ifstream paternalFile(repo + "0/HG06808_" + contigs_HG06808[contig] + "_pan010.txt");

    maternalKmerMap.clear();
    paternalKmerMap.clear();

    while (maternalFile >> kmer) {
        unsigned long long hashkmer = calculateHashString(kmer);
        maternalKmerMap[hashkmer] = number;
        //maternalKmerMapString[kmer] = number;
    }

    while (paternalFile >> kmer) {
        unsigned long long hashkmer = calculateHashString(kmer);
        paternalKmerMap[hashkmer] = number;
        //paternalKmerMapString[kmer] = number;
    }

    unsigned long long position = 0;
    unsigned long long mother_contig = 0;
    unsigned long long father_contig = 0;
    unsigned long long unknown_contig = 0;

    // Reading first 26 symbols in child fasta
    auto start = starting_reading(childFile);
    head = start.first;
    tail = start.second;

    int64_t value = 0;

    // Reading child fasta
    while (childFile.peek() != EOF && childFile.get(c)) {
        // Add new symbol
        int data = 0;
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
            case 'T':
                data = 3;
                break;
            case 'N':
                data = 4;
                break;
            default:
                answerFile << mother_contig << ' ' << father_contig << ' ' << unknown_contig << '\n';
                mother_contig = 0;
                father_contig = 0;
                unknown_contig = 0;

                outFile << c;
                answerFile << c;
                while (childFile.peek() != '\n') {
                    childFile.get(c);
                    outFile << c;
                    answerFile << c;
                }
                childFile.get(c);
                outFile << c << '\n';
                answerFile << c << '\n';
                position = 0;

                while (head != nullptr) {
                    Node* temp = head;
                    head = head->next;
                    delete temp;
                }
                return;
        }
        // Update tail of list
        tail->next = new Node(data);
        tail = tail->next;

        // (re)calculating hash of list
        if (position == 0) {
            value = calculateHash(head);
        } else {
            value = recalculateHash(value, head, tail);
        }

        // Update head of list
        if (position != 0) {
            Node* temp = head;
            head = head->next;
            delete temp;
        }

        // Find origin of current kmer
        if (paternalKmerMap.find(value) != paternalKmerMap.end() && maternalKmerMap.find(value) != maternalKmerMap.end()) {
            outFile << 'u';
            unknown++;
            unknown_contig++;
        } else if (paternalKmerMap.find(value) != paternalKmerMap.end()) {
            outFile << 'p';
            father++;
            father_contig++;
        } else if (maternalKmerMap.find(value) != maternalKmerMap.end()) {
            outFile << 'm';
            mother++;
            mother_contig++;
        } else {
            outFile << 'u';
            unknown++;
            unknown_contig++;
        }

        position += 1;
    }

    // Clean up memory
    while (head != nullptr) {
        Node* temp = head;
        head = head->next;
        delete temp;
    }

    maternalFile.close();
    paternalFile.close();

    answerFile << mother_contig << ' ' << father_contig << ' ' << unknown_contig << '\n';
}

// Function to find originality of kmers in given sample
void findOrigin(std::string &maternalPath, std::string &paternalPath, const std::string &childFastaFile, const std::string &outputFile, const std::string &answer) {
    // Open files
    std::ifstream childFile(childFastaFile);
    std::ofstream outFile(outputFile);
    std::ofstream answerFile(answer);

    std::cout << "============= counting =============\n";

    for (int contig = 0; contig < contigs_HG06808.size(); contig++) {
        findOriginInsideContig(childFile, outFile, answerFile, contig);
    }

    answerFile << "\n\n ========= summary ========\n";
    answerFile << mother << ' ' << father << ' ' << unknown << '\n';
}

int main(int argc, char* argv[]) {
    // Check if the number of arguments is correct
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " maternal_intersection paternal_intersection child.fasta potential_crossingovers.txt answer.txt" << std::endl;
        return 1;
    }

    std::string maternalOriginFile = argv[1];
    std::string paternalOriginFile = argv[2];
    std::string childFastaFile = argv[3];
    std::string outputFile = argv[4];
    std::string answer = argv[5];

    findOrigin(maternalOriginFile, paternalOriginFile, childFastaFile, outputFile, answer);
    return 0;
}
