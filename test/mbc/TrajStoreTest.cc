//
// Created by Chuang on 2019/6/11.
//
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <tuple>
#include <map>
#include<stdlib.h>
#include<time.h>
#include <cmath>
#include <spatialindex/SpatialIndex.h>
#include "storagemanager/TrajStore.h"

using namespace std;
using namespace SpatialIndex;


#define sourceFile "D://t1000.txt"
#define maxLinesToRead 1e10
#define dimension 2

template<class Type>
Type stringToNum(const std::string &str) {
    std::istringstream iss(str);
    Type num;
    iss >> num;
    return num;
}

struct xyt {
    double x;
    double y;
    double t;
};

xyt makemid(xyt p1, xyt p2, double t) {
    if (t > p2.t)
        cout << p1.x << " " << p2.x << endl <<
             p1.y << " " << p2.y << endl <<
             p1.t << " " << p2.t << " " << t << endl;
    assert(p1.t <= t);
    assert(t <= p2.t);
    double h1 = (t - p1.t) / (p2.t - p1.t);
    double h2 = (p2.t - t) / (p2.t - p1.t);
    double x = h2 * p1.x + h1 * p2.x;
    double y = h2 * p1.y + h1 * p2.y;
    xyt ret = {x, y, t};
    return ret;
}

vector<pair<id_type, Trajectory> > loadGTToTrajs() {
    //first level: vector of time period
    //second level: vector of segments in the time period
    cerr << "loading generated trajectories from txt to trajectories" << endl;

    ifstream inFile(sourceFile, ios::in);
    string lineStr;
    set<id_type> ids;
    multimap<id_type, xyt> trajs;
    vector<pair<id_type, Trajectory> > res;
    int curLine = 0;
    double minx = 40000, maxx = 0, miny = 40000, maxy = 0;
    while (getline(inFile, lineStr) && curLine < maxLinesToRead) {
        string str;
        stringstream ss(lineStr);
        getline(ss, str, '\t');
        getline(ss, str, '\t');
        int id = stringToNum<int>(str);
        getline(ss, str, '\t');
        getline(ss, str, '\t');
        getline(ss, str, '\t');
        double t = stringToNum<double>(str);
        getline(ss, str, '\t');
        double x = stringToNum<double>(str);
        getline(ss, str, '\t');
        double y = stringToNum<double>(str);
        getline(ss, str, '\t');
        double speed = stringToNum<double>(str);
        xyt p = {x, y, t};
        if (x > maxx) maxx = x;
        if (x < minx) minx = x;
        if (y > maxy) maxy = y;
        if (y < miny) miny = y;
        ids.insert(id);
        trajs.insert(make_pair(id, p));
        curLine++;
    }
    cout << minx << " " << maxx << " " << miny << " " << maxy << endl;
    for (auto id:ids) {
        multimap<id_type, xyt>::iterator beg, end, iter;
        vector<xyt> traj;
        beg = trajs.lower_bound(id);
        end = trajs.upper_bound(id);
        for (iter = beg; iter != end; iter++) {
            traj.emplace_back(iter->second);
        }
        trajs.erase(id);
        if (traj.size() >= 10) {
            vector<STPoint> tps;
            for (auto p:traj) {
                double xy[] = {p.x, p.y};
                tps.emplace_back(STPoint(xy, p.t, dimension));
            }
            if (!tps.empty()) {
                res.emplace_back(make_pair(id, Trajectory(tps)));
            }
        }
    }
    return res;
}

int main() {
    try {
        srand((int) time(NULL));
        vector<pair<id_type, Trajectory> > trajs = loadGTToTrajs();
        vector<pair<id_type, Trajectory> > empty1;
        vector<pair<id_type, vector<Trajectory>>> segs;
        for (auto traj:trajs) {
            segs.emplace_back(make_pair(traj.first, traj.second.getSegments(2000)));
//            vector<Trajectory> tmp(1);
//            tmp[0]=traj.second;
//            segs.push_back(make_pair(traj.first,tmp));
        }
        cout << trajs[15].second << endl;
//        trajs.swap(empty1);
        string name0 = "name0";
        IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096);
        StorageManager::IBuffer *file0 = StorageManager::createNewRandomEvictionsBuffer(*diskfile0, 10, false);
        auto ts = new TrajStore(name0, diskfile0, 4096);
        ts->loadSegments(segs);
//        for(auto e:ts.m_entries){
//            cout<<e.first<<endl;
//            cout<<e.second->m_page<<" "<<e.second->m_start<<" "<<e.second->m_len<<endl;
//        }
        id_type id = 1500;
        Trajectory tj;
        tj = ts->getTrajByTime(id, 0, 1000);
        cout << trajs[15].second << endl;
        cout << tj;
    }
    catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
}