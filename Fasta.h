#pragma once

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>



class Fasta {
public:
	std::string file_path;
	std::ofstream logging;

	typedef	std::unordered_map<std::string, std::string> SeqInfo;

	std::vector<SeqInfo> seqinfo;
	std::vector<int64_t> seqdata_offset;


	Fasta(const std::string& file_path);

	std::array<int64_t, 3> calc_score(std::unordered_map<int64_t, int64_t>& scorelist, int64_t seglen, int64_t max_len);
	
	std::pair<int64_t, int64_t> calc_score(int x, std::unordered_map<int64_t, int64_t>& scorelist, int64_t seg_len, int64_t max_len);
	std::pair<int64_t, int64_t> calc_score(int x, std::vector<int64_t>& scorelist, int64_t seglen, int64_t max_len);

	std::string get_seq(int x, int64_t st, int64_t len);
	void write_seq(int x, int64_t st, int64_t len,std::string file_path);

	static void get_hash(const std::string& s, int seg_len, std::unordered_map<int64_t, int64_t>& h);
	static void get_hash(const std::string& s, int seg_len, std::vector<int64_t>& lookup);
	static bool find(const std::vector<int64_t>& lookup, int64_t val);
private:
	void string_to_info(const std::string& s, SeqInfo& res);
};