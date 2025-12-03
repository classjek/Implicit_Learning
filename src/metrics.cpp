#include "metrics.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <sys/resource.h>
#include <unistd.h>

namespace metrics {

std::size_t current_rss_bytes() {
    std::ifstream f("/proc/self/statm");
    long pages_res = 0, dummy = 0;
    if (f) { f >> dummy >> pages_res; }
    const long page = sysconf(_SC_PAGESIZE);
    return pages_res > 0 ? static_cast<std::size_t>(pages_res) * static_cast<std::size_t>(page) : 0;
}

std::size_t peak_rss_bytes() {
    rusage ru{};
    if (getrusage(RUSAGE_SELF, &ru) == 0)
        return static_cast<std::size_t>(ru.ru_maxrss) * 1024ULL; // KiB -> bytes
    return 0;
}

std::string human_bytes(std::size_t b) {
    const char* u[] = {"B","KiB","MiB","GiB","TiB"};
    int i = 0;
    double x = static_cast<double>(b);
    while (x >= 1024.0 && i < 4) { x /= 1024.0; ++i; }
    char buf[48];
    std::snprintf(buf, sizeof(buf), "%.2f %s", x, u[i]);
    return std::string(buf);
}


Checkpoint::Checkpoint(std::string label)
    : start_(Clock::now()), last_(start_), label_(std::move(label)) {
    std::ostringstream ss;
    ss << "[Start] " << label_ << " | RSS="  << human_bytes(current_rss_bytes()) << " | Peak=" << human_bytes(peak_rss_bytes());
    lines_.push_back(ss.str());
}

void Checkpoint::tick(const std::string& name) {
    auto now = Clock::now();
    double dt_last  = std::chrono::duration<double>(now - last_).count();
    double dt_total = std::chrono::duration<double>(now - start_).count();
    last_ = now;

    std::ostringstream ss;
    ss << "[Tick] " << name << " | +dt="   << dt_last  << "s" << " | total=" << dt_total << "s"
       << " | RSS="   << human_bytes(current_rss_bytes()) << " | Peak="  << human_bytes(peak_rss_bytes());
    lines_.push_back(ss.str());
}

void Checkpoint::print(std::ostream& os) const {
    for (const auto& line : lines_) os << line << '\n';
}

} 
