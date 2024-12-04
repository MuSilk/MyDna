#include <iostream>
#include <fstream>
#include <unordered_map>
#include <deque>
#include <string>

#include "Fasta.h"
#include "utils.h"

using namespace std;
typedef long long ll;

int main() {
	ifstream in("sample.txt");
	string sample; in >> sample;
    string rc_sample = rc(sample);

    Fasta reference_database("reference/reference.fasta");

    unordered_map<ll, ll> sample_hash;
    auto st = clock();
    Fasta::get_hash(sample, 50, sample_hash);
    auto ed = clock();
    std::cout << "hash time:" << ed - st << "ms" << std::endl;
       
    st = clock();
    auto [pos, score] = reference_database.calc_score(5, sample_hash, 50, sample.length() * 1.5);
    ed = clock();
    std::cout << "deal time:" << ed - st << "ms" << std::endl;

    /*vector<int64_t> lookup;
    auto st = clock();
    Fasta::get_hash(sample, 50, lookup);
    auto ed = clock();
    std::cout << "hash time:" << ed - st << "ms" << std::endl;

    st = clock();
    auto [pos, score] = reference_database.calc_score(5, lookup, 50, sample.length() * 1.5);
    ed = clock();
    std::cout << "deal time:" << ed - st << "ms" << std::endl;*/

    //auto [id, pos, score] = reference_database.calc_score(sample_hash, 50, sample.length() * 1.5);
       
    //reference_database.write_seq(id, pos, sample.length() * 2, "reference.txt");

   /* in = ifstream("reference.txt");
    string reference; in >> reference;
    cout << sample.length() << " " << reference.length() << "\n";

    auto seg_match = DNA_match(reference, sample, 50, 30, 7);
    MatchGroup matches;
    SegMatchToMatchGroup(seg_match, matches);
    std::cout << calculate_value(matches, reference, sample, 30) << " ";
    printMatchGroup(matches, pos);
    MatchGroup_optimize(reference, sample, matches, 30, pos);*/
    return 0;
}