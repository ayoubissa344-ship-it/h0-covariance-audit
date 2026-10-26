// h0audit/io.hpp -- CSV reader for the compilation (columns: id,group,h0,sigma[,label]).
#pragma once
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "linalg.hpp"

namespace h0 {

struct Dataset {
    std::vector<std::string> id;
    Vec yL, sL, yR, sR;
};

inline std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> f;
    std::string cur;
    for (char c : line) {
        if (c == ',') { f.push_back(cur); cur.clear(); }
        else if (c != '\r') cur.push_back(c);
    }
    f.push_back(cur);
    return f;
}

inline Dataset read_dataset(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string line;
    do {
        if (!std::getline(in, line)) throw std::runtime_error("empty file");
    } while (line.empty() || line[0] == '#');  // skip leading comment lines
    const auto hdr = split_csv_line(line);
    int ig = -1, ih = -1, is = -1, ii = -1;
    for (std::size_t k = 0; k < hdr.size(); ++k) {
        if (hdr[k] == "group") ig = static_cast<int>(k);
        else if (hdr[k] == "h0") ih = static_cast<int>(k);
        else if (hdr[k] == "sigma") is = static_cast<int>(k);
        else if (hdr[k] == "id") ii = static_cast<int>(k);
    }
    if (ig < 0 || ih < 0 || is < 0) throw std::runtime_error("header must contain group,h0,sigma");
    Dataset d;
    std::size_t ln = 1;
    while (std::getline(in, line)) {
        ++ln;
        if (line.empty() || line[0] == '#') continue;
        const auto f = split_csv_line(line);
        if (f.size() <= static_cast<std::size_t>(std::max(ig, std::max(ih, is))))
            throw std::runtime_error("line " + std::to_string(ln) + ": too few fields");
        const double y = std::stod(f[ih]), s = std::stod(f[is]);
        if (f[ig] == "L") { d.yL.push_back(y); d.sL.push_back(s); }
        else if (f[ig] == "R") { d.yR.push_back(y); d.sR.push_back(s); }
        else throw std::runtime_error("line " + std::to_string(ln) + ": group must be L or R");
        d.id.push_back(ii >= 0 ? f[ii] : std::to_string(ln));
    }
    if (d.yL.size() < 2 || d.yR.size() < 2) throw std::runtime_error("need >= 2 rows per group");
    return d;
}

}  // namespace h0
