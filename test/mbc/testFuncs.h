//
// Created by Chuang on 2019/9/1.
//

#ifndef SPATIALINDEX_TESTFUNCS_H
#define SPATIALINDEX_TESTFUNCS_H

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
#include <dirent.h>
#include <stdlib.h>
#include <stdio.h>
#include <strings.h>

#if !WIN32

#include <sys/mman.h>
#include <unistd.h>

#endif

#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

#define random(x, y) (((double)rand()/RAND_MAX)*(y-x)+x)

#include <spatialindex/SpatialIndex.h>


#include "storagemanager/TrajStore.h"
#include "../../src/storagemanager/DiskStorageManager.h"
#include "../../src/mbcrtree/MBCRTree.h"
#include "../../src/rtree/RTree.h"
//#define sourceFile "D://t1000.txt"
#define genFile "D://00.txt"
#define GLFile "/root/GLSC.csv"
#define fileFolder "/root/out/"
#define maxLinesToRead 1e10
#define testtime 100
#define dimension 2
#define indexcap 10
#define leafcap 10000
#define QueryType 2
//1 for range, 2 for 5-NN

using namespace std;
using namespace SpatialIndex;

int drop_cache(int drop) {
    int ret = 0;
#if !WIN32
    int fd = 0;
    sync();
    fd = open("/proc/sys/vm/drop_caches", O_RDWR);
    if (fd < 0) {
        return -1;
    }
    char dropData[32] = {0};
    int dropSize = snprintf(dropData, sizeof(dropData), "%d", drop);

    ret = write(fd, dropData, dropSize);
    close(fd);
#endif
    return ret;

}

class MyVisitor : public IVisitor {
public:
    size_t m_indexIO;
    size_t m_leafIO;
    size_t m_indexvisited;
    size_t m_leafvisited;
    size_t m_resultGet;
    id_type m_lastResult;
    double m_lastDist = 0;
    IShape *m_query;
    TrajStore *ts = nullptr;

public:
    MyVisitor() : m_indexIO(0), m_leafIO(0), m_resultGet(0), m_indexvisited(0), m_leafvisited(0) {}

    void visitNode(const INode &n) {
//        if (n.isLeaf()) m_leafIO++;
//        else m_indexIO++;
        uint32_t size = n.getIndexByteArraySize();

        if (n.isLeaf()) {
            m_leafvisited++;
            m_leafIO += size;
        }
        else {
            m_indexvisited++;
            m_indexIO += size;
        }
    }

    void visitData(const IData &d) {
        m_resultGet++;
        IShape *pS;
        d.getShape(&pS);
        // do something.
        delete pS;
//        cout<<"data"<<endl;

        // data should be an array of characters representing a Region as a string.
        uint8_t *pData = 0;
        uint32_t cLen = 0;
        d.getData(cLen, &pData);
        m_leafIO += cLen;
        // do something.
//        double *s = reinterpret_cast<double*>(pData);
//        cout << *s << endl;

        m_lastResult = d.getIdentifier();
//        cerr << d.getIdentifier()<<std::endl;
//        auto mou=dynamic_cast<const MBCRTree::MBCRTree::simpleData*>(&d);
//        if(mou!=nullptr){
//            m_lastDist=mou->m_dist;
//            cerr << d.getIdentifier()<<"\t"<<mou->m_dist << endl;
//        }
//        auto mou2=dynamic_cast<const RTree::RTree::simpleData*>(&d);
//        if(mou2!=nullptr){
//            m_lastDist=mou2->m_dist;
//            cerr << d.getIdentifier()<<"\t"<<mou2->m_dist << endl;
//        }
//        //id of the data
//        if(ts== nullptr)
//            cerr << d.getIdentifier() << endl;
//        else
//            cerr << d.getIdentifier() /100<< endl;
//        //the traj
//        Trajectory traj;
//        if(ts== nullptr) {
//            traj.loadFromByteArray(pData);
//            double mindist=m_query->getMinimumDistance(traj);
//            cerr<<"traj dist is"<<mindist<<"\n\n";
//        }
//        else{
//            id_type tid=d.getIdentifier();
//            traj=ts->getTrajByTime(tid,0,1000);
//            auto brs=ts->getMBRsByTime(tid,0,1000);
//            auto bcs=ts->getMBCsByTime(tid,0,1000);
//            cerr<<"traj dist is"<<m_query->getMinimumDistance(traj)<<"\n"
//                    <<m_query->getMinimumDistance(brs)<<"\n"
//                    <<m_query->getMinimumDistance(bcs)<<"\n\n";
//        }

        delete[] pData;
    }

    void visitData(std::vector<const IData *> &v) {
        cout << v[0]->getIdentifier() << " " << v[1]->getIdentifier() << endl;
    }
};

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

template<class Type>
Type stringToNum(const std::string &str) {
    std::istringstream iss(str);
    Type num;
    iss >> num;
    return num;
}

vector<pair<id_type, Trajectory> > loadGTToTrajs(string filename = genFile) {
    //first level: vector of time period
    //second level: vector of segments in the time period
#ifndef NDEBUG
    cerr << "loading generated trajectories from txt to trajectories" << endl;
#endif
    auto stat = trajStat::instance();
    ifstream inFile(filename, ios::in);
    string lineStr;
    set<id_type> ids;
    multimap<id_type, xyt> trajs;
    vector<pair<id_type, Trajectory> > res;
    int curLine = 0;
    while (getline(inFile, lineStr) && curLine < maxLinesToRead) {
        try {
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
            if (x > stat->maxx) stat->maxx = x;
            if (x < stat->minx) stat->minx = x;
            if (y > stat->maxy) stat->maxy = y;
            if (y < stat->miny) stat->miny = y;
            if (t > stat->maxt) stat->maxt = t;
            if (t < stat->mint) stat->mint = t;
            ids.insert(id);
            trajs.insert(make_pair(id, p));
            curLine++;
        }
        catch (...) {
            break;
        }
    }
    stat->Dx = stat->maxx - stat->minx;
    stat->Dy = stat->maxy - stat->miny;
    stat->Dt = stat->maxt - stat->mint;
    for (auto id:ids) {
        multimap<id_type, xyt>::iterator beg, end, iter;
        vector<xyt> traj;
        beg = trajs.lower_bound(id);
        end = trajs.upper_bound(id);
        for (iter = beg; iter != end; iter++) {
            traj.emplace_back(iter->second);
        }
        trajs.erase(id);
        if (traj.size() >= 2) {
            vector<STPoint> tps;
            xyt lastpoint;
            for (int l = 0; l < traj.size(); l++) {
                tps.emplace_back(STPoint(traj[l].x, traj[l].y, traj[l].t));
                if (l != 0) {
                    stat->dist += std::sqrt(sq(traj[l].x - lastpoint.x) + sq(traj[l].y - lastpoint.y));
                }
                lastpoint = traj[l];
            }
            if (!tps.empty()) {
                Trajectory tmp(tps);
                stat->lineCount += tps.size() - 1;
                stat->trajCount += 1;
                stat->M += tmp.m_endTime() - tmp.m_startTime();
                res.emplace_back(make_pair(id, tmp));
            }
        }
    }
#ifndef NDEBUG
    std::cerr << "load data finished\n";
#endif
    stat->tl = stat->M / stat->lineCount;
    stat->jt = stat->M / stat->trajCount;
    stat->v = stat->dist / stat->M;
    std::cerr << *stat;
//    drop_cache(3);
    return res;
}

vector<pair<id_type, Trajectory> > loadDumpedFiledToTrajs(string filename = genFile) {
    auto stat = trajStat::instance();
    ifstream inFile(filename, ios::in);
    string lineStr;
    set<id_type> ids;
    vector<pair<id_type, Trajectory> > res;
    Trajectory tj;
    Region r;
    getline(inFile, lineStr);
    stat->fromString(lineStr);
    int curLine = 0;
    while (getline(inFile, lineStr) && curLine < maxLinesToRead) {
        try {
            string str;
            stringstream ss(lineStr);
            getline(ss, str);
            id_type id = stringToNum<id_type>(str);
            getline(inFile, str);
            tj.loadFromString(str);
            if (tj.m_points.size() >= 2) {
                ids.insert(id);
                res.emplace_back(make_pair(id, tj));
                curLine++;
//                tj.getMBR(r);
//                if (r.m_pHigh[0] > stat->maxx) stat->maxx = r.m_pHigh[0];
//                if (r.m_pLow[0] < stat->minx) stat->minx = r.m_pLow[0];
//                if (r.m_pHigh[1] > stat->maxy) stat->maxy = r.m_pHigh[1];
//                if (r.m_pLow[1] < stat->miny) stat->miny = r.m_pLow[1];
//                if (r.m_pHigh[2] > stat->maxt) stat->maxt = r.m_pHigh[2];
//                if (r.m_pLow[2] < stat->mint) stat->mint = r.m_pLow[2];
//                stat->dist += tj.m_dist();
//                stat->lineCount += tj.m_points.size() - 1;
//                stat->trajCount += 1;
//                stat->M += tj.m_endTime() - tj.m_startTime();
            }
        }
        catch (...) {
            break;
        }
    }
//    stat->Dx = stat->maxx - stat->minx;
//    stat->Dy = stat->maxy - stat->miny;
//    stat->Dt = stat->maxt - stat->mint;
//    stat->tl = stat->M / stat->lineCount;
//    stat->jt = stat->M / stat->trajCount;
//    stat->v = stat->dist / stat->M;
    std::cerr << *stat;
//    drop_cache(3);
    return res;
}


void dumpToFile(vector<pair<id_type, Trajectory> > &trajs, string filename = "dumpedtraj.txt") {
    ofstream outFile(filename, ios::out);
    auto stat = trajStat::instance();
    outFile<<stat->toString()<<"\n";
    for (auto traj : trajs) {
        outFile << traj.first << "\n";
        outFile << traj.second.toString() << "\n";
    }
    outFile.close();
}


vector<pair<id_type, Trajectory> > loadGLToTrajs(string filename = GLFile) {
    //first level: vector of time period
    //second level: vector of segments in the time period
    cerr << "loading geolife trajectories from txt to trajectories" << endl;
    auto stat = trajStat::instance();
    ifstream inFile(filename, ios::in);
    string lineStr;
    set<id_type> ids;
    multimap<id_type, xyt> trajs;
    vector<pair<id_type, Trajectory> > res;
    int curLine = 0;
    getline(inFile, lineStr);

    while (getline(inFile, lineStr) && curLine < maxLinesToRead) {
        try {
            string str;
            stringstream ss(lineStr);
            getline(ss, str, ',');
            int id = stringToNum<int>(str);
            getline(ss, str, ',');
            double x = stringToNum<double>(str);
            getline(ss, str, ',');
            double y = stringToNum<double>(str);
            getline(ss, str, ',');
            double t = stringToNum<double>(str);
            xyt p = {x, y, t};
            if (x > stat->maxx) stat->maxx = x;
            if (x < stat->minx) stat->minx = x;
            if (y > stat->maxy) stat->maxy = y;
            if (y < stat->miny) stat->miny = y;
            if (t > stat->maxt) stat->maxt = t;
            if (t < stat->mint) stat->mint = t;
            ids.insert(id);
            trajs.insert(make_pair(id, p));
            curLine++;
        }
        catch (...) {
            break;
        }
    }

    stat->Dx = stat->maxx - stat->minx;
    stat->Dy = stat->maxy - stat->miny;
    stat->Dt = stat->maxt - stat->mint;
    double gdist = 0, gtime = 0;
    for (auto id:ids) {
        multimap<id_type, xyt>::iterator beg, end, iter;
        vector<xyt> traj;
        beg = trajs.lower_bound(id);
        end = trajs.upper_bound(id);
        for (iter = beg; iter != end; iter++) {
            traj.emplace_back(iter->second);
        }
        trajs.erase(id);
        if (traj.size() >= 2) {
            vector<STPoint> tps;
            xyt lastpoint;
            for (int l = 0; l < traj.size(); l++) {
                tps.emplace_back(STPoint(traj[l].x, traj[l].y, traj[l].t));
                if (l != 0) {
                    stat->dist += std::sqrt(sq(traj[l].x - lastpoint.x) + sq(traj[l].y - lastpoint.y));
                    if (std::sqrt(sq(traj[l].x - lastpoint.x) + sq(traj[l].y - lastpoint.y)) /
                        (traj[l].t - lastpoint.t) > 1e-5) {
                        gdist += std::sqrt(sq(traj[l].x - lastpoint.x) + sq(traj[l].y - lastpoint.y));
                        gtime += (traj[l].t - lastpoint.t);
                    }
                }
                lastpoint = traj[l];
            }
            if (!tps.empty()) {
                Trajectory tmp(tps);
                stat->lineCount += tps.size() - 1;
                stat->trajCount += 1;
                stat->M += tmp.m_endTime() - tmp.m_startTime();
                res.emplace_back(make_pair(id, tmp));
            }
        }
    }
    std::cerr << "vh is" << gdist / gtime << "\n";
#ifndef NDEBUG
    std::cerr << "load data finished\n";
#endif
    stat->tl = stat->M / stat->lineCount;
    stat->jt = stat->M / stat->trajCount;
    stat->v = stat->dist / stat->M;
    std::cerr << *stat;
    stat->Dx /= 6.5;
    stat->Dy /= 5;
    stat->v *= 2.5;
//    drop_cache(3);
    return res;
}

vector<pair<id_type, Trajectory> > loadGTFolder(int num = 10, string folder = fileFolder) {
    vector<pair<id_type, Trajectory> > res;
    vector<string> files;
    struct dirent *ptr;
    DIR *dir;
    string PATH = fileFolder;
    dir = opendir(PATH.c_str());
    while ((ptr = readdir(dir)) != NULL && files.size() < num) {
        if (ptr->d_name[0] == '.')
            continue;
        //cout << ptr->d_name << endl;
        files.emplace_back(PATH + ptr->d_name);
    }
    for (auto file:files) {
        vector<pair<id_type, Trajectory> > tmptrajs = loadGTToTrajs(file);
        res.insert(res.begin(), tmptrajs.begin(), tmptrajs.end());
    }
//    drop_cache(3);
    return res;
}

double rn(double n, double S, double Nq, double qt, double v2) {
    double a = std::sqrt(n * S / M_PI / Nq) * qt / 2;
    double b = qt / 4 * std::sqrt(2 * v2 * qt * v2 * qt + 4 * n * S / M_PI / Nq);
    return a + b;
}


double knncost(double bt, int k, double qt, int f, bool useMBR, double _rk) {
    auto stat = trajStat::instance();
    double v2 = stat->v * 2;
    double Nt = stat->trajCount;
    double Nq = Nt * (qt + 2 * stat->jt) / stat->Dt;
    Nq = min(Nq, Nt);
    double Ltc = pow(stat->Dx * stat->Dy * stat->Dt / stat->M * bt * f / v2 / v2, 1.0 / 3);
    double Lxc = Ltc * v2;
    double Lt = Ltc + bt;
//    double vv=min(-7e-9*Lt+0.0001,1e-4)+5e-5;
    double vv = v2;
    double Lx = vv * Lt;

    double r = sqrt(k * stat->Dx * stat->Dy / M_PI / Nq);
    double Rk = (r * qt / 2 + qt * sqrt(v2 * v2 * qt * qt * 2 + 4 * r * r) / 4);
    double rk = Rk / qt;

    double Rcost = (pow((2 * rk + Lx), 2) * (qt + Lt) + Lx * qt * qt * v2) / (pow(Ltc, 3) * v2 * v2);

    double rnq = rn(5, stat->Dx * stat->Dy, Nq, qt, v2) + stat->v / 2 * bt * qt / 2;
    double nmin = 5.0;
    double nmax = 10000.0;
    while (nmax - nmin > 0.1) {
        double mid = (nmax + nmin) / 2;
        if (rn(mid, stat->Dx * stat->Dy, Nq, qt, v2) > rnq) {
            nmax = mid;
        } else {
            nmin = mid;
        }
    }
    double Tcost = (nmax + nmin) / 4;
    return Rcost + Tcost;
}

struct d4 {
    d4(double _low, double _high, double _vlow, double _vhigh) {
        low = _low, high = _high, vlow = _vlow, vhigh = _vhigh;
    }

    double low, high, vlow, vhigh;
};

double biSearchMax(int k, double qt, int f, bool useMBR, double rk = -1, double low = 1, double high = 10000) {
    std::stack<d4> st;
    double vlow, vhigh;
    vlow = knncost(low, k, qt, f, useMBR, rk);
    vhigh = knncost(high, k, qt, f, useMBR, rk);
    st.push(d4(low, high, vlow, vhigh));
    d4 first(0, 0, 0, 0);
    auto stat = trajStat::instance();
    while (!st.empty()) {
        first = st.top();
        st.pop();
        double mid = (first.low + first.high) / 2;
        if (first.high - first.low < 1) {
            return std::max(mid, stat->tl);
        }
        double vmid = knncost(mid, k, qt, f, useMBR, rk);
        bool rising = knncost(mid + 0.1, k, qt, f, useMBR, rk) > vmid;
//        std::cerr<<mid<<"\t"<<vmid<<"\n";
        if (rising) {
            st.push(d4(first.low, mid, first.vlow, vmid));
        } else {
            st.push(d4(mid, first.high, vmid, first.vhigh));
        }
    }
    throw Tools::IllegalStateException("biSearch:failed searching");
}

void TreeQueryBatch(ISpatialIndex *tree, const vector<IShape *> &queries, TrajStore *ts = nullptr, int thennk = 5) {
    MyVisitor vis;
    vis.ts = ts;
    auto start = std::chrono::system_clock::now();
//    drop_cache(3);
    for (int i = 0; i < queries.size(); i++) {
//        drop_cache(3);
//        cerr<<"Query is "<<queries.at(i)->toString();
        if (QueryType == 1) {
            tree->intersectsWithQuery(*queries[i], vis);
        } else if (QueryType == 2) {
            vis.m_query = queries[i];
            tree->nearestNeighborQuery(thennk, *queries[i], vis);
//            cerr<<"finished "<<i<<"already\n";
        }
    }
    double time;
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    time = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
    cerr << "Querying time: " << time << endl;
    cerr << "VISIT NODE " << vis.m_indexvisited << "\t" << vis.m_leafvisited << endl;
    cerr << "TrajStore Statistic" << ts->m_indexIO << "\t" << ts->m_trajIO << endl;
}

double kNNQueryBatch(ISpatialIndex *tree, const vector<IShape *> &queries, TrajStore *ts = nullptr, int thennk = 5,
                     bool reportEnd = false) {
    ts->cleanStatistic();
    ts->flush();
    drop_cache(3);
    int num = queries.size();
    MyVisitor vis;
    vis.ts = ts;
    auto start = std::chrono::system_clock::now();
    double rad = 0;
    int indio = 0;
    std::vector<int> indios;
    for (int i = 0; i < queries.size(); i++) {
        vis.m_query = queries[i];
        tree->nearestNeighborQuery(thennk, *queries[i], vis);
        rad += vis.m_lastDist;
        if (reportEnd) std::cerr << "end\n";
        indios.emplace_back(ts->m_indexIO - indio);
        indio = ts->m_indexIO;
    }
    rad /= queries.size();
    double time;
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    time = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
//    cerr <<"Average Querying time: "<< time/num<<endl;
//    cerr <<"Averaged VISIT NODE "<<1.0*vis.m_indexvisited/num<<"\t"<<1.0*vis.m_leafvisited/num<<endl;
//    cerr <<"TrajStore Statistic"<< 1.0*ts->m_indexIO/num<<"\t"<<1.0*ts->m_trajIO/num<<endl;
    sort(indios.begin(), indios.end());
    int mid1 = queries.size() * 0.1, mid2 = queries.size() * 0.9;
    double sum = 0;
    for (int i = mid1; i <= mid2; i++) {
        sum += indios[i];
    }
    sum = sum / (mid2 - mid1);
    cerr << "time\tindexVisit\tLeafVisit\t leaf1\tleaf2\tindexIO\ttrajIO\tresult\n";
    cerr << time / num << "\t" << 1.0 * vis.m_indexvisited / num << "\t" << 1.0 * vis.m_leafvisited / num << "\t"
         << 1.0 * ts->m_leaf1 / num << "\t" << 1.0 * ts->m_leaf2 / num << "\t" << 1.0 * ts->m_indexIO / num << "\t"
         << 1.0 * ts->m_trajIO / num << "\t" << 1.0 * ts->m_loadedTraj / num << "\t" << sum << endl;
//    cerr <<time/num<<"\n";
    return rad;
}

void rangeQueryBatch(ISpatialIndex *tree, const vector<IShape *> &queries, TrajStore *ts = nullptr) {
    ts->cleanStatistic();
    int num = queries.size();
    MyVisitor vis;
    vis.ts = ts;
    sb = 0, sbb = 0;
    auto start = std::chrono::system_clock::now();
    for (int i = 0; i < queries.size(); i++) {
        vis.m_query = queries[i];
        tree->intersectsWithQuery(*queries[i], vis);
    }
    double time;
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    time = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
//    cerr <<"Average Querying time: "<< time/num<<endl;
//    cerr <<"Averaged VISIT NODE "<<1.0*vis.m_indexvisited/num<<"\t"<<1.0*vis.m_leafvisited/num<<endl;
//    cerr <<"TrajStore Statistic"<< 1.0*ts->m_indexIO/num<<"\t"<<1.0*ts->m_trajIO/num<<endl;
    cerr << time / num << "\t" << 1.0 * vis.m_indexvisited / num << "\t" << 1.0 * vis.m_leafvisited / num << "\t"
         << 1.0 * ts->m_indexIO / num << "\t" << 1.0 * ts->m_trajIO / num << "\t" << double(sbb) / sb << "\t" << sb
         << "\t" << sbb << endl;
//    cerr <<time/num<<"\n";
}

int TreeQuery(ISpatialIndex *tree, IShape *query, TrajStore *ts = nullptr) {
    clock_t start, end;
    MyVisitor vis;
    if (ts != nullptr)
        vis.ts = ts;
    vis.m_query = query;
    start = clock();
    if (QueryType == 1) {
        tree->intersectsWithQuery(*query, vis);
    } else if (QueryType == 2) {
        vis.m_query = query;
        tree->nearestNeighborQuery(5, *query, vis);
    }
    end = clock();
    if (QueryType == 1) {
        return vis.m_resultGet;
    } else if (ts != nullptr) {
        return vis.m_lastResult;
    } else return vis.m_lastResult;
}


#endif //SPATIALINDEX_TESTFUNCS_H
