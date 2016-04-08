#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <chrono>
#include <thread>

using namespace std;

#define LOCAL

// 'cause topcoder requires single file submission
#include "solution.cpp"


template<class T> void getVector(vector<T>& v) {
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}


int main() {
    int NP;
    cin >> NP;

    int Npoints;
    cin >> Npoints;
    vector<int> points(Npoints);
    getVector(points);

    int Nroots;
    cin >> Nroots;
    vector<int> roots(Nroots);
    getVector(roots);

    CutTheRoots cr;
    vector<int> ret = cr.makeCuts(NP, points, roots);

    // To sort of ensure that ErrorReader thread in tester will get a chance
    // to pick all stderr up.
    cerr.flush();
    std::this_thread::sleep_for(std::chrono::milliseconds(200));


    cout << ret.size() << endl;
    for (int i = 0; i < ret.size(); ++i) {
        cout << ret[i] << endl;
    }
    cout.flush();

    return 0;
}
