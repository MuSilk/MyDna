#include "Fasta.h"

#include <deque>
#include <array>
#include <algorithm>
#include <functional>

const int64_t M = 998244353;
const int64_t INF = 7;

std::unordered_map<char, int> char_map = {
	{'A' , 1},
	{'T' , 2},
	{'C' , 3},
	{'G' , 4}
};

int64_t qpow(int64_t a, int64_t b) {
	int64_t ans = 1;
	while (b) {
		if (b & 1)ans = ans * a % M;
		a = a * a % M;
		b >>= 1;
	}
	return ans;
}

void Fasta::get_hash(const std::string& s, int seg_len, std::unordered_map<int64_t, int64_t>& h) {
	h.clear();
	const int64_t LEFT_HASH = qpow(INF, seg_len - 1);
	std::deque<int64_t> buffer;
	int64_t last = 0;
	for (int i = 0; i < seg_len; i++) {
		last = (last * INF + char_map[s[i]]) % M;
		buffer.push_back(char_map[s[i]]);
	}
	h[last] += 1;

	for (int i = seg_len; i < s.length(); i++) {
		last = (last - buffer[0] * LEFT_HASH % M + M) % M;
		buffer.pop_front();
		last = (last * INF + char_map[s[i]]) % M;
		buffer.push_back(char_map[s[i]]);
		h[last] += 1;
	}
}

void Fasta::get_hash(const std::string& s, int seg_len, std::vector<int64_t>& lookup) {
	std::vector<int64_t> h;
	const int64_t LEFT_HASH = qpow(INF, seg_len - 1);
	std::deque<int64_t> buffer;
	int64_t last = 0;
	for (int i = 0; i < seg_len; i++) {
		last = (last * INF + char_map[s[i]]) % M;
		buffer.push_back(char_map[s[i]]);
	}
	h.push_back(last);

	for (int i = seg_len; i < s.length(); i++) {
		last = (last - buffer[0] * LEFT_HASH % M + M) % M;
		buffer.pop_front();
		last = (last * INF + char_map[s[i]]) % M;
		buffer.push_back(char_map[s[i]]);
		h.push_back(last);
	}
	std::sort(h.begin(), h.end());
	h.erase(std::unique(h.begin(), h.end()), h.end());

	lookup = std::vector<int64_t>(h.size() << 2, -1);
	std::function<void(int,int,int)> init = [&]( int l, int r, int cur)->void {
		if (r < l)return;
		if (l == r) {
			if (cur > lookup.size())std::cout << l << " " << cur << " " << lookup.size() << std::endl;
			lookup[cur] = h[l];
			return;
		}
		int m = (l + r) >> 1;
		lookup[cur] = h[m];
		init(l, m - 1, cur << 1);
		init(m + 1, r, cur << 1 | 1);
	};
	init(0, h.size() - 1, 1);
}

bool Fasta::find(const std::vector<int64_t>& lookup, int64_t val) {
	for (int i = 1; i < lookup.size() && lookup[i] != -1;) {
		if (val == lookup[i])return 1;
		else if (val < lookup[i])i <<= 1;
		else i = i << 1 | 1;
	}
	return 0;
}

Fasta::Fasta(const std::string& file_path) {
	logging = std::ofstream("out.log");
	this->file_path = file_path;

	FILE* fp;
	fopen_s(&fp, file_path.c_str(), "r");
	
	if (fp == NULL) {
		std::cout << "Fail to open the file" << std::endl;
		return;
	}

	char buf[512];
	int64_t st = -1;
	while (fgets(buf, 512, fp)) {
		if (buf[0] == '>') {

			SeqInfo info;
			string_to_info(buf, info);
			seqinfo.emplace_back(info);

			if (st != -1) seqdata_offset.emplace_back(st);
			st = _ftelli64(fp);
		}
	}
	seqdata_offset.emplace_back(st);
	fclose(fp);
}

//info start with '>'
void Fasta::string_to_info(const std::string& s, std::unordered_map<std::string, std::string>& res) {
	res.clear();
	std::string key = "name", value;
	bool cur = 1;
	std::cout << s << " " << std::endl;
	for (int i = 1; i < s.size(); i++) {
		if (cur == 0) {
			if (s[i] == ':')cur = 1;
			else if (s[i] == ' ')continue;
			else key.push_back(s[i]);
		}
		else {
			if (s[i] == ' ' || s[i] == '\n') {
				cur = 0, res[key] = value;
				key = value = "";
			}
			else value.push_back(s[i]);
		}
	}
}

std::array<int64_t,3> Fasta::calc_score(std::unordered_map<int64_t, int64_t>& scorelist, int64_t seglen, int64_t max_len) {
	int64_t best_pos = 0;
	int64_t best_score = 0;
	int64_t best_id = 0;
	for (int i = 0; i < seqinfo.size(); i++) {
		auto p = calc_score(i, scorelist, seglen, max_len);
		if (p.second > best_score) {
			best_pos = p.first;
			best_score = p.second;
			best_id = i;
		}
		logging << "best id:" << best_id << " best pos:" << best_pos << " best score:" << best_score << std::endl;
	}
	return { best_id,best_pos,best_score };
}

std::pair<int64_t,int64_t> Fasta::calc_score(int x, std::unordered_map<int64_t, int64_t>& scorelist, int64_t seg_len, int64_t max_len) {
	
	int64_t pos = 0;
	std::deque<int64_t> hash_buf;
	std::deque<std::pair<int64_t,int64_t>> score_buf;
	const int64_t LEFT_HASH = qpow(INF, seg_len - 1);
	int64_t last = 0;
	int64_t score = 0;

	int64_t best_pos = 0;
	int64_t best_score = 0;

	FILE* fp;
	fopen_s(&fp, file_path.c_str(), "r");
	_fseeki64(fp, seqdata_offset[x], SEEK_SET);
	char buf[512];
	while (fgets(buf, 512, fp)) {
		if (buf[0] == '>')break;

		for (int i = 0; i < 512 && buf[i] != '\n'; i++) {
			while (score_buf.size() > 0 && score_buf[0].first < pos - max_len + 1) {
				score -= score_buf[0].second;
				score_buf.pop_front();
			}

			if (hash_buf.size() >= seg_len) {
				last = (last - hash_buf[0] * LEFT_HASH % M + M) % M;
				hash_buf.pop_front();
			}
			last = (last * INF + char_map[buf[i]]) % M;
			hash_buf.push_back(char_map[buf[i]]);


			if (hash_buf.size() == seg_len) {
				auto p = scorelist.find(last);
				if (p != scorelist.end()) {
					score += p->second;
					score_buf.push_back({ pos - seg_len + 1 ,p->second });
				}
			}
			
			if (score > best_score) {
				best_score = score;
				best_pos = pos - max_len + 1;
			}

			if (pos % 1000000 == 0) {
				logging << x << " " << pos << " " << best_pos << " " << best_score << std::endl;
			}
			pos++;
		}
	}
	fclose(fp);
	return { best_pos,best_score };
}

std::pair<int64_t, int64_t> Fasta::calc_score(int x, std::vector<int64_t>& scorelist, int64_t seg_len, int64_t max_len) {

	int64_t pos = 0;
	std::deque<int64_t> hash_buf;
	std::deque<std::pair<int64_t, int64_t>> score_buf;
	const int64_t LEFT_HASH = qpow(INF, seg_len - 1);
	int64_t last = 0;
	int64_t score = 0;

	int64_t best_pos = 0;
	int64_t best_score = 0;

	FILE* fp;
	fopen_s(&fp, file_path.c_str(), "r");
	_fseeki64(fp, seqdata_offset[x], SEEK_SET);
	char buf[512];
	while (fgets(buf, 512, fp)) {
		if (buf[0] == '>')break;

		for (int i = 0; i < 512 && buf[i] != '\n'; i++) {
			while (score_buf.size() > 0 && score_buf[0].first < pos - max_len + 1) {
				score -= score_buf[0].second;
				score_buf.pop_front();
			}

			if (hash_buf.size() >= seg_len) {
				last = (last - hash_buf[0] * LEFT_HASH % M + M) % M;
				hash_buf.pop_front();
			}
			last = (last * INF + char_map[buf[i]]) % M;
			hash_buf.push_back(char_map[buf[i]]);


			if (hash_buf.size() == seg_len) {
				if (find(scorelist, last)) {
					score += 1;
					score_buf.push_back({ pos - seg_len + 1 ,1 });
				}
			}

			if (score > best_score) {
				best_score = score;
				best_pos = pos - max_len + 1;
			}

			if (pos % 1000000 == 0) {
				logging << x << " " << pos << " " << best_pos << " " << best_score << std::endl;
			}
			pos++;
		}
	}
	fclose(fp);
	return { best_pos,best_score };
}

std::string Fasta::get_seq(int x, int64_t st, int64_t len) {
	FILE* fp;
	fopen_s(&fp, file_path.c_str(), "r");
	_fseeki64(fp, seqdata_offset[x], SEEK_SET);
	std::string res;

	int64_t pos = 0;
	char buf[512];
	while (fgets(buf, 512, fp)) {
		if (buf[0] == '>')break;

		for (int i = 0; i < 512 && buf[i] != '\n'; i++) {
			if (pos >= st + len)break;
			if (pos >= st)res.push_back(buf[i]);
			pos++;
		}
		if (pos >= st + len)break;
	}
	fclose(fp);
	return res;
}

void Fasta::write_seq(int x, int64_t st, int64_t len, std::string file_path) {
	std::string res = get_seq(x, st, len);
	std::ofstream out(file_path);
	out << res;
}



