#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cmath>

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

// -559050727 % 1000000007 = 440949280
// -118101447 % 1000000007 = 881898560
// -677152174 % 1000000007 = 322847833

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
        default:
            hashValue += 322847833;
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
            default:
                data = 3;
        }

        hashValue = (hashValue * base + data) % PRIME;
    }
    return hashValue;
}

//std::unordered_map<unsigned long long, int> maternalKmerMap;
//std::unordered_map<std::string, int> maternalKmerMapStrings;
//while (maternalFile >> kmer >> number) {
//unsigned long long hashkmer = calculateHashString(kmer);
//maternalKmerMap[hashkmer] = number;
//maternalKmerMapStrings[kmer] = number;
////std::cout << hashkmer << '\n';
//}
//
//std::unordered_map<unsigned long long, int> paternalKmerMap;
//std::unordered_map<std::string, int> paternalKmerMapStrings;
//while (paternalFile >> kmer >> number) {
//unsigned long long hashkmer = calculateHashString(kmer);
//paternalKmerMap[hashkmer] = number;
//paternalKmerMapStrings[kmer] = number;
////std::cout << hashkmer << '\n';
//}

// Function to find originality of kmers in given sample
void findPotentialCrossingOvers(const std::string &maternalOriginFile, const std::string &paternalOriginFile, const std::string &childFastaFile, const std::string &outputFile) {
    // Open files
    std::ifstream maternalFile(maternalOriginFile);
    std::ifstream paternalFile(paternalOriginFile);

    std::string kmer;
    int number = 1;

    // Reading kmers with maternal origin
    // Counting their hash
    std::unordered_map<unsigned long long, int> maternalKmerMap;
    std::cout << "Reading data with maternal origin...\n";
    while (maternalFile >> kmer) {
        unsigned long long hashkmer = calculateHashString(kmer);
        maternalKmerMap[hashkmer] = number;

        //std::cout << kmer << " - " << hashkmer << '\n';
    }

    // Reading kmers with paternal origin
    // Counting their hash
    std::cout << "Reading data with paternal origin...\n";
    std::unordered_map<unsigned long long, int> paternalKmerMap;
    while (paternalFile >> kmer) {
        unsigned long long hashkmer = calculateHashString(kmer);
        paternalKmerMap[hashkmer] = number;

       // std::cout << kmer << " - " << hashkmer << '\n';
    }

    std::cout << maternalKmerMap.size() << " " << paternalKmerMap.size() << '\n';

    std::cout << "============= counting =============\n";

    // Open files
    std::ifstream childFile(childFastaFile);
    std::ofstream outFile(outputFile);

    Node* head = nullptr;
    Node* tail = nullptr;
    char c;

    unsigned long long position = 0;
    unsigned long long mother = 0;
    unsigned long long father = 0;
    unsigned long long unknown = 0;

    // Reading first 26 symbols in child fasta
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
            default:
                data = 3;
        }

        if (head == nullptr) {
            head = new Node(data);
            tail = head;
        } else {
            tail->next = new Node(data);
            tail = tail->next;
        }
    }

    int64_t value = 0;

    // Reading child fasta
    while (childFile.get(c)) {
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
            default:
                data = 3;
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
        if (paternalKmerMap.find(value) != paternalKmerMap.end()) {
            outFile << 'p';
            father++;
        } else if (maternalKmerMap.find(value) != maternalKmerMap.end()) {
            outFile << 'm';
            mother++;
        } else {
            outFile << 'u';
            unknown++;
        }

        position += 1;
    }

    // Clean up memory
    while (head != nullptr) {
        Node* temp = head;
        head = head->next;
        delete temp;
    }

    std::ofstream answer("answer.txt");
    answer << "maternal - " << mother << '\n';
    answer << "paternal - " << father << '\n';
    answer << "unknown - " << unknown << '\n';
}

int main(int argc, char* argv[]) {
    // Check if the number of arguments is correct
    if (argc != 5 && argc != 6) {
        std::cerr << "Usage: " << argv[0] << " maternal_intersection paternal_intersection child.fasta potential_crossingovers.txt" << std::endl;
        return 1;
    }

    std::string maternalOriginFile = argv[1];
    std::string paternalOriginFile = argv[2];
    std::string childFastaFile = argv[3];
    std::string outputFile = argv[4];

    findPotentialCrossingOvers(maternalOriginFile, paternalOriginFile, childFastaFile, outputFile);
    return 0;
}
