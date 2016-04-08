#include <algorithm>
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <random>
#include <utility>
#include <set>
#include <sys/time.h>

#include "pretty_printing.h"

using namespace std;


#define debug(x) \
    cerr << #x " = " << (x) << endl
#define debug2(x, y) \
    cerr << #x " = " << (x) \
    << ", " #y " = " << (y) << endl
#define debug3(x, y, z) \
    cerr << #x " = " << (x) \
    << ", " #y " = " << (y) \
    << ", " #z " = " << (z) << endl


#ifdef LOCAL
const double TIME_LIMIT = 7.0;
#else
const double TIME_LIMIT = 9.0;
#endif

int get_time_cnt = 0;
double get_time() {
    get_time_cnt++;
    timeval tv;
    gettimeofday(&tv, nullptr);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}


vector<pair<double*, double> > undo_log;
void record(double &d) {
    undo_log.emplace_back(&d, d);
}

class Checkpoint {
    int n;
public:
    Checkpoint() {
        n = undo_log.size();
    }
    ~Checkpoint() {
        for (int i = undo_log.size() - 1; i >= n; i--) {
            *undo_log[i].first = undo_log[i].second;
        }
        undo_log.erase(undo_log.begin() + n, undo_log.end());
    }
};


// What they call 'root' in the problem statement.
struct Segment {

    Segment *parent;
    int x1, y1, x2, y2;
    vector<Segment*> children;
    double fraction;

    // Stuff below describes the whole subtree.
    double min_x, min_y, max_x, max_y;
    double total_length;

    void update() {
        record(total_length);
        record(min_x);
        record(min_y);
        record(max_x);
        record(max_y);

        int dx = x2 - x1;
        int dy = y2 - y1;
        total_length = fraction * sqrt(dx*dx + dy*dy);

        min_x = min(x1, x2);
        min_y = min(y1, y2);
        max_x = max(x1, x2);
        max_y = max(y1, y2);

        if (fraction == 1.0) {
            for (const Segment *child : children) {
                total_length += child->total_length;
                min_x = min(min_x, child->min_x);
                min_y = min(min_y, child->min_y);
                max_x = max(max_x, child->max_x);
                max_y = max(max_y, child->max_y);
            }
        }
    }
};

default_random_engine rnd_gen(42);

int NP;
vector<Segment> segments;

struct Line {
    // defined by the equation
    //   nx*x + ny*y + nn = 0

    Line() = default;
    Line(int x1, int y1, int x2, int y2)
            : x1(x1), y1(y1), x2(x2), y2(y2) {
        nx = y2 - y1;
        ny = x1 - x2;
        nn = - nx * x1 - ny * y1;
    }

    int eval(int x, int y) const {
        return nx * x + ny * y + nn;
    }

    int x1, y1, x2, y2, nx, ny, nn;
};

ostream& operator<<(ostream &out, const Line &line) {
    out << "Line("
        << line.x1 << ", " << line.y1 << ", "
        << line.x2 << ", " << line.y2 << ")";
    return out;
}


Line get_random_separator(int x1, int y1, int x2, int y2) {
    int min_x = min(x1, x2) - 2;
    if (min_x < 0) min_x = 0;
    int min_y = min(y1, y2) - 2;
    if (min_y < 0) min_y = 0;
    int max_x = max(x1, x2) + 2;
    if (max_x > 1024) max_x = 1024;
    int max_y = max(y1, y2) + 2;
    if (max_y > 1024) max_y = 1024;

    while (true) {
        int ax = uniform_int_distribution<>(min_x, max_x)(rnd_gen);
        int ay = uniform_int_distribution<>(min_y, max_y)(rnd_gen);
        int bx = uniform_int_distribution<>(0, 1024)(rnd_gen);
        int by = uniform_int_distribution<>(0, 1024)(rnd_gen);
        Line line(ax, ay, bx, by);
        int e1 = line.eval(x1, y1);
        int e2 = line.eval(x2, y2);
        if (e1 < 0 && e2 > 0 || e1 > 0 && e2 < 0)
            return line;
    }
}


double cut(Segment *seg, const Line &line) {
    double min_n = line.nn +
        min(line.nx * seg->min_x, line.nx * seg->max_x) +
        min(line.ny * seg->min_y, line.ny * seg->max_y);
    double max_n = line.nn +
        max(line.nx * seg->min_x, line.nx * seg->max_x) +
        max(line.ny * seg->min_y, line.ny * seg->max_y);
    assert(min_n <= max_n);

    if (max_n < 0 || min_n > 0) {
        // no cut, no traversal
        return 0.0;
    }

    double t = seg->total_length;

    int n1 = line.eval(seg->x1, seg->y1);
    if (n1 == 0) {
        record(seg->fraction);
        seg->fraction = 0.0;
        seg->update();
        return t - seg->total_length;
    }

    int n2 = line.eval(seg->x2, seg->y2);

    if (n1 > 0 && n2 <= 0 || n1 < 0 && n2 >= 0)  {
        double new_fraction = 1.0 * n1 / (n1 - n2) - 1e-8;
        if (new_fraction < seg->fraction) {
            record(seg->fraction);
            seg->fraction = new_fraction;
            seg->update();
        }
        return t - seg->total_length;
    }

    // no cut
    if (seg->fraction == 1.0) {
        bool need_update = false;
        for (Segment *child : seg->children) {
            if (cut(child, line) > 0.0)
                need_update = true;
        }
        if (need_update)
            seg->update();
    }
    return t - seg->total_length;
}

double cut(const Line &line) {
    double result = 0.0;
    for (int i = 0; i < NP; i++)
        result += cut(&segments[i], line);
    return result;
}


struct Separation {
    int segment1;
    int segment2;
    double best_loss;
    Line best_line;

    Separation() = default;

    bool operator<(const Separation &other) const {
        return best_loss > other.best_loss ||
               best_loss == other.best_loss &&
               segment1 * 200 + segment2 < other.segment1 * 200 + other.segment2;
    }

    bool is_split(const Line &line) const {
        int e1 = line.eval(segments[segment1].x1, segments[segment1].y1);
        int e2 = line.eval(segments[segment2].x1, segments[segment2].y1);
        return e1 <= 0 && e2 >= 0 || e1 >= 0 && e2 <= 0;
    }
};

ostream &operator<<(ostream &out, const Separation &sep) {
    out << "Sep("
        << sep.segment1 << ", "
        << sep.segment2 << ", loss=" << sep.best_loss << ")";
    return out;
}


class CutTheRoots {
public:
    vector<int> makeCuts(int NP, vector<int> points, vector<int> roots) {
        (void)roots;  // unused

#ifndef LOCAL
        assert(false && "asserts should be disabled");
#endif

        double start_time = get_time();
        ::NP = NP;

        assert(roots.size() % 2 == 0);
        segments.resize(roots.size() / 2 + NP);

        for (int i = 0; i < NP; i++) {
            Segment &seg = segments[i];
            seg.x1 = seg.x2 = points[i * 2];
            seg.y1 = seg.y2 = points[i * 2 + 1];
            seg.parent = nullptr;
        }
        for (int i = 0; i < roots.size() / 2; i++) {
            assert(roots[2 * i + 1] == i + NP);

            Segment &parent = segments[roots[2 * i]];
            Segment &seg = segments[i + NP];

            parent.children.push_back(&seg);
            seg.parent = &parent;

            seg.x1 = parent.x2;
            seg.y1 = parent.y2;
            seg.x2 = points[(i + NP) * 2];
            seg.y2 = points[(i + NP) * 2 + 1];
        }

        for (int i = segments.size() - 1; i >= 0; i--) {
            Segment &seg = segments[i];
            seg.fraction = 1.0;
            seg.update();
        }
        double total_length = 0.0;
        for (int i = 0; i < NP; i++)
            total_length += segments[i].total_length;

        vector<Line> solution;

        set<Separation> separations;

        for (int i = 0; i < NP; i++) {
            for (int j = 0; j < i; j++) {
                Separation sep;
                sep.segment1 = j;
                sep.segment2 = i;
                sep.best_loss = 1e10;
                separations.insert(sep);
            }
        }

        while (!separations.empty()) {
            for (int i = 0; i < 2000; i++) {
                Separation sep = *separations.begin();

                Line line = get_random_separator(
                    segments[sep.segment1].x1, segments[sep.segment1].y1,
                    segments[sep.segment2].x1, segments[sep.segment2].y1);
                double loss;
                {
                    Checkpoint cp;
                    loss = cut(line);
                }

                set<Separation>::iterator it = separations.begin();
                while (it != separations.end()) {
                    if (loss >= it->best_loss)
                        break;

                    if (it->is_split(line)) {
                        Separation sep2 = *it;
                        sep2.best_loss = loss;
                        sep2.best_line = line;
                        it = separations.erase(it);
                        separations.insert(sep2);
                    } else {
                        ++it;
                    }
                }
            }

            double loss = separations.begin()->best_loss;
            assert(loss < 1e10);
            Line line = separations.begin()->best_line;
            // debug2(loss, separations.size());

            set<Separation>::iterator it = separations.begin();
            while (it != separations.end()) {
                if (it->is_split(line)) {
                    it = separations.erase(it);
                } else {
                    ++it;
                }
            }

            cut(line);
            solution.push_back(line);
        }

        debug(separations.size());

        debug(solution.size());
        vector<int> result;
        for (Line line : solution) {
            result.push_back(line.x1);
            result.push_back(line.y1);
            result.push_back(line.x2);
            result.push_back(line.y2);
        }

        double remaining_length = 0.0;
        for (int i = 0; i < NP; i++)
            remaining_length += segments[i].total_length;
        debug(remaining_length);

        double predicted_score = 1e6 * remaining_length / total_length;
        debug(predicted_score);

        cerr << "it took " << get_time() - start_time << endl;
        return result;
    }
};
