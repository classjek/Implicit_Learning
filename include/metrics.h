#pragma once

#include <chrono>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

namespace metrics {

// Current resident set size (bytes)
std::size_t current_rss_bytes();

// Peak resident set size since process start (bytes)
std::size_t peak_rss_bytes();

// human readable byte string
std::string human_bytes(std::size_t b);

class Checkpoint {
public:
    explicit Checkpoint(std::string label = "start");
    void tick(const std::string& name);
    void print(std::ostream& os = std::cout) const;

private:
    using Clock = std::chrono::steady_clock;
    Clock::time_point start_;
    Clock::time_point last_;
    std::string label_;
    std::vector<std::string> lines_; // buffered output
};

}
