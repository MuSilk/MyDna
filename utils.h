#pragma once
#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include <array>
#include <random>

#include "3rdparty/edlib.h"

std::mt19937_64 gen(time(NULL));

template<typename T,typename updateF,typename scoringF>
std::pair<T, int64_t> simulated_annealing(
    T origin,
    updateF update,
    scoringF scoring,
    double t0, double kt = 0.97, int round = 20) {

    static std::uniform_real_distribution<double> rnd(0,1);

    T data = origin;
    int64_t score = scoring(origin);
    T final_data = data;
    int64_t final_score = score;
    while (t0 >= 1) {
        for (int i = 0; i < round; i++) {
            T nxt_data = update(data, t0);
            int64_t nxt_score = scoring(nxt_data);

            int64_t delta = nxt_score - score;
            
            if (exp(delta / t0) > rnd(gen)) {
                data = nxt_data;
                score = nxt_score;
                if (score > final_score) {
                    final_data = data;
                    final_score = score;
                }
            }
        }
        data = final_data;
        score = final_score;
        t0 *= kt;
    }
    return { final_data,final_score };
}

std::string rc(std::string seq) {
    reverse(seq.begin(), seq.end());
    for (int i = 0; i < seq.size(); i++) {
        switch (seq[i]){
            case 'A':seq[i] = 'T'; break;
            case 'T':seq[i] = 'A'; break;
            case 'C':seq[i] = 'G'; break;
            case 'G':seq[i] = 'C'; break;
        default:
            break;
        }
    }
    return seq;
}

int64_t calculate_distance(
    const std::string& ref,
    const std::string& query,
    int64_t ref_st, int64_t ref_en,
    int64_t query_st, int64_t query_en) {

    EdlibAlignResult result = edlibAlign(ref.c_str() + ref_st, ref_en - ref_st, query.c_str() + query_st, query_en - query_st, edlibDefaultAlignConfig());
    int res = result.editDistance;
    edlibFreeAlignResult(result);
    return res;
}    

typedef std::array<int64_t,4> Match;
typedef std::vector<std::array<int64_t, 4>> MatchGroup;
int64_t calculate_value(
    const MatchGroup& points,
    const std::string& ref, 
    const std::string& query,
    int64_t penalty) {

    int64_t editdistance = 0;
    int64_t aligned = 0;
    int64_t preend = 0;
    for (auto [query_st, query_en, ref_st, ref_en] : points) {
        if (preend > query_st) return 0;
        if (query_en < query_st)return 0;
        if (ref_en < ref_st)return 0;
        preend = query_en;
        editdistance += calculate_distance(ref, query, ref_st, ref_en, query_st, query_en);
        aligned += query_en - query_st;
    }
    return std::max(aligned - editdistance - (int64_t)points.size() * penalty, 0ll);
}
    

typedef std::array<int64_t,3> PhaseScore;

class SegMatch {
public:
    int64_t l, r;
    std::vector<PhaseScore> scorelist;
    SegMatch(int64_t l, int64_t r, const std::vector<PhaseScore>& ml) {
        this->l = l;
        this->r = r;
        this->scorelist = ml;
    }
};

std::vector<SegMatch> DNA_match(
    const std::string& reference,
    const std::string& sample,
    int64_t seglen, int64_t penalty, int64_t topk) {

    std::vector<SegMatch> match_list;
    for (int64_t i = 0; i < sample.length(); i += seglen) {

        if (i % 1000 == 0)std::cout << "\r" << i << ":" << sample.length();

        int64_t ri = std::min(i + seglen, (int64_t)sample.length());
        std::vector<PhaseScore> scorelist;
        for (int j = 0; j < reference.length(); j += seglen) {
            int64_t rj = std::min(j + seglen, (int64_t)reference.length());
            scorelist.push_back({ ri - i - calculate_distance(reference, sample, j, rj, i, ri), j, rj });
        }
        nth_element(scorelist.begin(), scorelist.begin() + topk, scorelist.end(), [](const PhaseScore& a, const PhaseScore& b) {
            return a[0] > b[0];
            });
        scorelist.resize(std::min((size_t)topk, scorelist.size()));

        for (int j = 0; j < scorelist.size(); j++) {
            MatchGroup origin = { {
                i,ri,scorelist[j][1],scorelist[j][2]
            } };
            auto update = [&reference](MatchGroup data, double t0) {
                int64_t dv = std::max((int64_t)t0, 5ll);
                std::uniform_int_distribution<int64_t> rnd(-dv, dv);
                data[0][2] += rnd(gen);
                data[0][2] = std::clamp(data[0][2], 0ll, (int64_t)reference.length());
                data[0][3] += rnd(gen);
                data[0][3] = std::clamp(data[0][3], 0ll, (int64_t)reference.length());
                return data;
            };

            auto scoring = [&reference, &sample, &penalty](MatchGroup data) {
                return calculate_value(data, reference, sample, penalty);
            };

            auto [data, score] = simulated_annealing(origin, update, scoring, 3.0, 0.97, 5);
            scorelist[j] = { score,data[0][2],data[0][3] };

        }
        sort(scorelist.begin(), scorelist.end(), [](const PhaseScore& a, const PhaseScore& b) {
            return a[0] > b[0];
            });
        while (scorelist.size() > 0 && scorelist.back()[0] <= 0)scorelist.pop_back();
        if (scorelist.size() > 0) {
            match_list.push_back({ i,ri,scorelist });
        }
    }
    std::cout << "\ninit finished" << std::endl;

    std::vector<bool> tag(match_list.size(), true);
    while (true) {
        int64_t mergel = 0;
        std::vector<SegMatch> new_match_list;
        std::vector<bool> new_tag;
        while (mergel < match_list.size()) {
            std::cout << "\r" << mergel << ":" << match_list.size();
            if (mergel + 1 == match_list.size()) {
                new_match_list.push_back(match_list[mergel]);
                new_tag.push_back(0);
                break;
            }
            if (tag[mergel] == 0 && tag[mergel + 1] == 0) {
                new_match_list.push_back(match_list[mergel]);
                new_tag.push_back(0);
                mergel += 1;
                continue;
            }

            int64_t l = match_list[mergel].l;
            int64_t r = match_list[mergel + 1].r;
            int64_t lenl = match_list[mergel].r - match_list[mergel].l;
            int64_t lenr = match_list[mergel + 1].r - match_list[mergel + 1].l;

            std::vector<PhaseScore> scorelist;
            int64_t scorebound = match_list[mergel].scorelist[0][0] + match_list[mergel + 1].scorelist[0][0];

            auto merge = [&](int64_t refl, int64_t refr) {
                MatchGroup origin = { {l,r,refl,refr} };
                int64_t ref_length = reference.length();

                auto update = [&ref_length](MatchGroup data, double t0) {
                    int64_t dv = std::max((int64_t)t0, 10ll);
                    std::uniform_int_distribution<int64_t> rnd(-dv, dv);
                    data[0][2] += rnd(gen);
                    data[0][2] = std::clamp(data[0][2], 0ll, ref_length);
                    data[0][3] += rnd(gen);
                    data[0][3] = std::clamp(data[0][3], 0ll, ref_length);
                    return data;
                };

                auto scoring = [&reference, &sample, &penalty](MatchGroup data) {
                    return calculate_value(data, reference, sample, penalty);
                };

                auto [data, score] = simulated_annealing(origin, update, scoring, std::max((r - l) / 50, 5ll), 0.96, 10);
                if (score > scorebound) {

                    auto IoU = [](int64_t l1, int64_t r1, int64_t l2, int64_t r2) {
                        int64_t d1 = std::max(r1 - l1, 0ll);
                        int64_t d2 = std::max(r2 - l2, 0ll);
                        int64_t d3 = std::max(std::min(r1, r2) - std::max(l1, l2), 0ll);
                        if (d3 > 0.9 * (d1 + d2 - d3))return false;
                        return true;
                    };

                    bool flag = true;
                    for (auto& i : scorelist) {
                        if (!IoU(i[1], i[2], data[0][2], data[0][3])) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag)scorelist.push_back({ score,data[0][2],data[0][3] });
                }
            };

            for (auto& j : match_list[mergel].scorelist)
                merge(j[1], std::min(j[2] + lenr, (int64_t)reference.length()));
            for (auto& j : match_list[mergel + 1].scorelist)
                merge(std::max(0ll, j[1] - lenl), j[2]);

            if (scorelist.size() > 0) {
                sort(scorelist.begin(), scorelist.end(), [&](PhaseScore& a, PhaseScore& b) {return a[0] > b[0]; });
                scorelist.resize(std::min(scorelist.size(), (size_t)topk));
                new_match_list.push_back({ l,r,scorelist });
                new_tag.push_back(1);
                mergel += 2;
            }
            else {
                new_match_list.push_back(match_list[mergel]);
                new_tag.push_back(0);
                mergel += 1;
            }
        }
        if (new_match_list.size() < match_list.size()) {
            match_list = new_match_list;
            tag = new_tag;
            std::cout << "segment count: " << match_list.size() << "\n";
        }
        else {
            std::cout << "segment count: " << match_list.size() << "\n";
            break;
        }
    }
    return match_list;
}

void SegMatchToMatchGroup(
    const std::vector<SegMatch>& match_list,
    MatchGroup& res) {

    res.clear();
    for (auto& match : match_list) {
        res.push_back({ match.l, match.r, match.scorelist[0][1], match.scorelist[0][2] });
    }
}

void printMatchGroup(const MatchGroup& match_list, int64_t base = 0) {
    std::cout << "[";
    for (auto [i,ri,j,rj] : match_list) {
        std::cout << "(" << i << " " << ri << " " << j + base << " " << rj + base << ")";
    }
    std::cout << "]" << std::endl;
}

void MatchGroup_optimize(const std::string& reference, const std::string& sample, MatchGroup& match_list, int64_t penalty, int64_t base = 0) {

    int64_t ref_length = reference.length();
    int64_t sam_length = sample.length();
    
    for (size_t i = 0; i < match_list.size() * 5; i++) {
        size_t l = i % (match_list.size() - 1);

        int64_t saml = (l == 0) ? 0 : match_list[l - 1][1];
        int64_t samr = (l + 2 == match_list.size()) ? sample.length() : match_list[l + 2][0];

        auto scoring = [&](MatchGroup data) {
            return calculate_value(data, reference, sample, penalty);
        };
        auto update = [&](MatchGroup data, double t0) {
            int64_t dv = std::max((int64_t)t0, 10ll);
            std::uniform_int_distribution<int64_t> rnd(-dv, dv);
            data[0][0] += rnd(gen);
            data[0][0] = std::clamp(data[0][0], saml, samr);
            data[0][1] += rnd(gen);
            data[0][1] = std::clamp(data[0][1], saml, samr);

            data[1][0] = data[0][1];
            data[1][1] += rnd(gen);
            data[1][1] = std::clamp(data[1][1], saml, samr);

            auto update1 = [&](MatchGroup data, double t0) {
                int64_t dv = std::max((int64_t)t0, 10ll);
                for (int i = 0; i < data.size(); i++) {
                    data[i][2] += rnd(gen);
                    data[i][2] = std::clamp(data[i][2], 0ll, ref_length);
                    data[i][3] += rnd(gen);
                    data[i][3] = std::clamp(data[i][3], 0ll, ref_length);
                }
                return data;
            };

            auto [newdata, score] = simulated_annealing(data, update1, scoring, t0, 0.96, 10);

            return newdata;
        };

        MatchGroup data = {
            match_list[l],match_list[l + 1]
        };

        auto [newdata, score] = simulated_annealing(data, update, scoring, 50, 0.96, 10);
        match_list[l] = newdata[0];
        match_list[l + 1] = newdata[1];
        std::cout << scoring(match_list) << " ";
        printMatchGroup(match_list);
    }
               
}