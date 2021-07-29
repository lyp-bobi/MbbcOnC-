//
// Created by Chuang on 2020/11/24.
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
#include "storagemanager/xStore.h"
#include "../../src/xrtree/xRTree.h"
#include <storagemanager/tjsql.h>
#if !WIN32

#include <sys/mman.h>
#include <unistd.h>

#endif

#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

#define xRP shared_ptr<xRTree>


#include <spatialindex/SpatialIndex.h>


#include "storagemanager/xStore.h"
#include "../../src/storagemanager/DiskStorageManager.h"
//#define sourceFile "D://t1000.txt"
#define genFile "D://00.txt"
#define GLFile "/root/GLSC.csv"
#ifndef WIN32
    #define fileFolder "/root/out/"
#else
    #define fileFolder "D://out/"
#endif
#define maxLinesToRead 1e10
#define indexcap 10
#define leafcap 10000
//extern int QueryType;
//1 for range, 2 for 5-NN
#include <thread>
using namespace std;
using namespace SpatialIndex;
using namespace xRTreeNsp;
#if defined(TJDEBUG) || defined(WIN32) || !defined(NDEBUG)
#define NUMCORE 1

extern int NUMTHREAD=NUMCORE;
extern double testtime = 100;
#else
#define NUMCORE 4
extern double testtime = 4000;
extern int NUMTHREAD=NUMCORE;
#endif
extern bool testxfirstOutput = true;


using namespace std;

static int drop_cache(int drop) {
#if !WIN32
    if(drop == 3) {
        sync();
    }
    {
        ofstream ofs("/proc/sys/vm/drop_caches");
        ofs<< drop<<endl;
        ofs.close();
    }
    {
        ofstream ofs("/sys/block/vda/queue/nr_requests");
        ofs<< 32<<endl;
        ofs.close();
    }
    {
        ofstream ofs("/sys/block/vda/queue/max_sectors_kb");
        ofs<< 16<<endl;
        ofs.close();
    }

#endif
    return 0;
}


class MyVisitor : public IVisitor {
public:
    size_t m_indexvisited=0;
    size_t m_leafvisited=0;
    size_t m_resultGet;
    id_type m_lastResult;
    double m_lastDist = 0;
    IShape *m_query;
    xStore *ts = nullptr;

public:
    MyVisitor() : m_resultGet(0), m_indexvisited(0), m_leafvisited(0) {}

    void visitNode(const INode &n) {
        if (n.isLeaf()) {
            m_leafvisited++;
        }
        else {
            m_indexvisited++;
        }
    }
    void visitData(std::vector<const IData*>& v){}
    void visitData(const IData &d) {
        m_resultGet++;
        m_lastResult = d.getIdentifier();
        auto mou=dynamic_cast<const xRTreeNsp::xRTree::simpleData*>(&d);
        if(mou!=nullptr){
            m_lastDist=mou->m_dist;
#if !defined(NDEBUG) || defined(TJDEBUG)
            cerr <<"result" << d.getIdentifier()<<"\t"<<mou->m_dist << endl;
#endif
        }
    }

    void clear(){
        m_indexvisited = m_leafvisited = m_resultGet=0;
    }
};

struct xyt {
    double x;
    double y;
    double t;
};

static xyt makemid(xyt p1, xyt p2, double t) {
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

static vector<pair<id_type, xTrajectory> > loadGTToTrajs(string filename = genFile) {
    //first level: vector of time period
    //second level: vector of segments in the time period
    cerr << "loading generated trajectories from"<<filename<<" to trajectories" << endl;

    ifstream inFile(filename, ios::in);
    string lineStr;
    set<id_type> ids;
    multimap<id_type, xyt> trajs;
    vector<pair<id_type, xTrajectory> > res;
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
            if (x > tjstat->maxx) tjstat->maxx = x;
            if (x < tjstat->minx) tjstat->minx = x;
            if (y > tjstat->maxy) tjstat->maxy = y;
            if (y < tjstat->miny) tjstat->miny = y;
            if (t > tjstat->maxt) tjstat->maxt = t;
            if (t < tjstat->mint) tjstat->mint = t;
            ids.insert(id);
            trajs.insert(make_pair(id, p));
            curLine++;
        }
        catch (...) {
            break;
        }
    }
    tjstat->Dx = tjstat->maxx - tjstat->minx;
    tjstat->Dy = tjstat->maxy - tjstat->miny;
    tjstat->Dt = tjstat->maxt - tjstat->mint;
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
            vector<xPoint> tps;
            xyt lastpoint;
            for (int l = 0; l < traj.size(); l++) {
                tps.emplace_back(xPoint(traj[l].x, traj[l].y, traj[l].t));
                if (l != 0) {
                    tjstat->dist += std::sqrt(sq(traj[l].x - lastpoint.x) + sq(traj[l].y - lastpoint.y));
                }
                lastpoint = traj[l];
            }
            if (!tps.empty()) {
                xTrajectory tmp(tps);
                tjstat->lineCount += tps.size() - 1;
                tjstat->trajCount += 1;
                tjstat->M += tmp.m_endTime() - tmp.m_startTime();
                res.emplace_back(make_pair(id, tmp));
            }
        }
    }
#ifndef NDEBUG
    std::cerr << "load data finished\n";
#endif
    tjstat->tl = tjstat->M / tjstat->lineCount;
    tjstat->jt = tjstat->M / tjstat->trajCount;
    tjstat->v = tjstat->dist / tjstat->M;
    std::cerr << *tjstat;
//    tjstat->usedata("od");
    inFile.close();
//    drop_cache(3);
    return res;
}

static vector<pair<id_type, xTrajectory> > loadDumpedFiledToTrajs(string filename = genFile) {

    tjstat->init();
    ifstream inFile(filename, ios::in);
    string lineStr;
    set<id_type> ids;
    vector<pair<id_type, xTrajectory> > res;
    xTrajectory tj;
    xMBR r;
//    tjstat->fromString(lineStr);
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
                tj.getxMBR(r);
                if (r.m_xmax > tjstat->maxx) tjstat->maxx = r.m_xmax;
                if (r.m_xmin < tjstat->minx) tjstat->minx = r.m_xmin;
                if (r.m_ymax > tjstat->maxy) tjstat->maxy = r.m_ymax;
                if (r.m_ymin < tjstat->miny) tjstat->miny = r.m_ymin;
                if (r.m_tmax > tjstat->maxt) tjstat->maxt = r.m_tmax;
                if (r.m_tmin < tjstat->mint) tjstat->mint = r.m_tmin;
                tjstat->dist += tj.m_dist();
                tjstat->lineCount += tj.m_points.size() - 1;
                tjstat->trajCount += 1;
                if(tjstat->vmax < tj.maxSpeed()){
                    tjstat->vmax=tj.maxSpeed();
                }
                tjstat->M += tj.m_endTime() - tj.m_startTime();
            }
        }
        catch (...) {
            break;
        }
    }
    tjstat->Dx = tjstat->maxx - tjstat->minx;
    tjstat->Dy = tjstat->maxy - tjstat->miny;
    tjstat->Dt = tjstat->maxt - tjstat->mint;
    tjstat->tl = tjstat->M / tjstat->lineCount;
    tjstat->jt = tjstat->M / tjstat->trajCount;
    tjstat->v = tjstat->dist / tjstat->M;
    tjstat->Sr = (tjstat->Dx + tjstat->Dy) / 2;
    tjstat->P = tjstat->Dt;
    inFile.close();
    std::cerr << *tjstat;
    std::cerr << tjstat->toString();
//    drop_cache(3);
    std::cerr<<filename<<endl;
    if (filename.find("tdexpand")!=filename.npos) tjstat->usedata("tdexpand");
    else if (filename.find("td")!=filename.npos) tjstat->usedata("td");
    else if (filename.find("gl")!=filename.npos) tjstat->usedata("glexpand");
    else if (filename.find("gl")!=filename.npos) tjstat->usedata("gl");
    return res;
}


static id_type dumpToFile(vector<pair<id_type, xTrajectory> > &trajs, string filename = "dumpedtraj.txt", int num =-1,id_type id=0) {
    ofstream outFile(filename, ios::out);

    //outFile << tjstat->toString() << "\n";
    if( num ==-1){
        num = trajs.size();
    }
    else{
        num = min(num,int(trajs.size()));
    }
    for(int i=0;i<num;i++){
        auto traj = trajs[i];
        outFile << id << "\n";
        outFile << traj.second.toString() << "\n";
        id++;
    }
    cerr<<"dumping to "<<filename<<endl;
    outFile.close();
    return id;
}


static int getLastId(string s)
{
    string filename = s;
    ifstream fin;

    fin.open(filename);
    if(fin.is_open()) {
        fin.seekg(0,ios_base::end);                // go to one spot before the EOF
        for(int i=0;i<2;i++) {
            fin.seekg(-2,ios_base::cur);
            bool keepLooping = true;
            while (keepLooping) {
                char ch;
                ch = fin.peek();                            // Get current byte's data
//                std::cerr<<ch;
                if ((long long) fin.tellg() <= 1) {             // If the data was at or before the 0th byte
                    fin.seekg(0);                       // The first line is the last line
                    keepLooping = false;                // So stop there
                } else if (ch == '\n') {                   // If the data was a newline
                    keepLooping = false;                // Stop at the current position.
                } else {                                  // If the data was neither a newline nor at the 0 byte
                    fin.seekg(-1,ios_base::cur);
                    // Move to the front of that data, then to the front of the data before it
                }
            }
        }
        string lastLine;
        getline(fin,lastLine);                      // Read the current line
        getline(fin,lastLine);                      // Read the current line
        fin.close();
        cerr<<"last id is "<<lastLine<<endl;
        return stoll(lastLine);
    }

    return 0;
}

static id_type dumpToFile_append(vector<pair<id_type, xTrajectory> > &trajs, string filename = "dumpedtraj.txt", int num =-1, id_type id = -1) {
    if(id==-1) id = getLastId(filename) + 1;
    ofstream outFile(filename, ios::app);
//    outFile << tjstat->toString() << "\n";
    if( num == -1){
        num = trajs.size();
    }
    else{
        num = min(num,int(trajs.size()));
    }
    for(int i=0;i<num;i++){
        auto traj = trajs[i];
        outFile << id << endl;
        outFile << traj.second.toString() << endl;
        id++;
    }
    cerr<<"dumping to "<<filename<<endl;
    outFile.flush();
    outFile.close();
    return id;
}


static vector<pair<id_type, xTrajectory> > loadGLToTrajs(string filename = GLFile) {
    //first level: vector of time period
    //second level: vector of segments in the time period
    cerr << "loading geolife trajectories from txt to trajectories" << endl;

    ifstream inFile(filename, ios::in);
    string lineStr;
    set<id_type> ids;
    multimap<id_type, xyt> trajs;
    vector<pair<id_type, xTrajectory> > res;
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
            if (x > tjstat->maxx) tjstat->maxx = x;
            if (x < tjstat->minx) tjstat->minx = x;
            if (y > tjstat->maxy) tjstat->maxy = y;
            if (y < tjstat->miny) tjstat->miny = y;
            if (t > tjstat->maxt) tjstat->maxt = t;
            if (t < tjstat->mint) tjstat->mint = t;
            ids.insert(id);
            trajs.insert(make_pair(id, p));
            curLine++;
        }
        catch (...) {
            break;
        }
    }

    tjstat->Dx = tjstat->maxx - tjstat->minx;
    tjstat->Dy = tjstat->maxy - tjstat->miny;
    tjstat->Dt = tjstat->maxt - tjstat->mint;
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
            vector<xPoint> tps;
            xyt lastpoint;
            for (int l = 0; l < traj.size(); l++) {
                tps.emplace_back(xPoint(traj[l].x, traj[l].y, traj[l].t));
                if (l != 0) {
                    tjstat->dist += std::sqrt(sq(traj[l].x - lastpoint.x) + sq(traj[l].y - lastpoint.y));
                    if (std::sqrt(sq(traj[l].x - lastpoint.x) + sq(traj[l].y - lastpoint.y)) /
                        (traj[l].t - lastpoint.t) > 1e-5) {
                        gdist += std::sqrt(sq(traj[l].x - lastpoint.x) + sq(traj[l].y - lastpoint.y));
                        gtime += (traj[l].t - lastpoint.t);
                    }
                }
                lastpoint = traj[l];
            }
            if (!tps.empty()) {
                xTrajectory tmp(tps);
                tjstat->lineCount += tps.size() - 1;
                tjstat->trajCount += 1;
                tjstat->M += tmp.m_endTime() - tmp.m_startTime();
                res.emplace_back(make_pair(id, tmp));
            }
        }
    }
    std::cerr << "vh is" << gdist / gtime << "\n";
#ifndef NDEBUG
    std::cerr << "load data finished\n";
#endif
    tjstat->tl = tjstat->M / tjstat->lineCount;
    tjstat->jt = tjstat->M / tjstat->trajCount;
    tjstat->v = tjstat->dist / tjstat->M;
    std::cerr << *tjstat;
    tjstat->Dx /= 6.5;
    tjstat->Dy /= 5;
    tjstat->v *= 2.5;
    inFile.close();
//    drop_cache(3);
    return res;
}

static vector<pair<id_type, xTrajectory> > loadGTFolder(int num = 10, string folder = fileFolder) {
    vector<pair<id_type, xTrajectory> > res;
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
        vector<pair<id_type, xTrajectory> > tmptrajs = loadGTToTrajs(file);
        res.insert(res.begin(), tmptrajs.begin(), tmptrajs.end());
    }
//    drop_cache(3);
    return res;
}

static double rn(double n, double S, double Nq, double qt, double v2) {
    double a = std::sqrt(n * S / M_PI / Nq) * qt / 2;
    double b = qt / 4 * std::sqrt(2 * v2 * qt * v2 * qt + 4 * n * S / M_PI / Nq);
    return a + b;
}


static double knncost(double bt, int k, double qt, int f, bool useMBR, double _rk) {

    double v2 = tjstat->v * 2;
    double Nt = tjstat->trajCount;
    double Nq = Nt * (qt + 2 * tjstat->jt) / tjstat->Dt;
    Nq = min(Nq, Nt);
    double Ltc = pow(tjstat->Dx * tjstat->Dy * tjstat->Dt / tjstat->M * bt * f / v2 / v2, 1.0 / 3);
    double Lxc = Ltc * v2;
    double Lt = Ltc + bt;
//    double vv=min(-7e-9*Lt+0.0001,1e-4)+5e-5;
    double vv = v2;
    double Lx = vv * Lt;

    double r = sqrt(k * tjstat->Dx * tjstat->Dy / M_PI / Nq);
    double Rk = (r * qt / 2 + qt * sqrt(v2 * v2 * qt * qt * 2 + 4 * r * r) / 4);
    double rk = Rk / qt;

    double Rcost = (pow((2 * rk + Lx), 2) * (qt + Lt) + Lx * qt * qt * v2) / (pow(Ltc, 3) * v2 * v2);

    double rnq = rn(5, tjstat->Dx * tjstat->Dy, Nq, qt, v2) + tjstat->v / 2 * bt * qt / 2;
    double nmin = 5.0;
    double nmax = 10000.0;
    while (nmax - nmin > 0.1) {
        double mid = (nmax + nmin) / 2;
        if (rn(mid, tjstat->Dx * tjstat->Dy, Nq, qt, v2) > rnq) {
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

static double biSearchMax(int k, double qt, int f, bool useMBR, double rk = -1, double low = 50, double high = 10000) {

    double mincost=1e300, bestbt;
    double gap=low;
    if(tjstat->dataset == "od") gap = 5;
    for (double i =low;i<high;i+=gap){
        double c = tjstat->knncost(i, k, qt);
//        std::cerr<<i<<"\t"<<c<<"\n";
        if(c<mincost){
            mincost = c;
            bestbt = i;
        }
    }
    tjstat->knncost(bestbt, k, qt, true);
    if(tjstat->tl > bestbt){
        bestbt = int(tjstat->tl / gap) * gap;
    }
    return bestbt;
}
//
//void TreeQueryBatch(xRTree *tree, const vector<IShape *> &queries, xStore *ts = nullptr, int thennk = 5) {
//    MyVisitor vis;
//    vis.ts = ts;
//    auto start = std::chrono::system_clock::now();
////    drop_cache(3);
//    for (int i = 0; i < queries.size(); i++) {
////        drop_cache(3);
////        cerr<<"Query is "<<queries.at(i)->toString();
//        if (QueryType == 1) {
//            tree->intersectsWithQuery(*queries[i], vis);
//        } else if (QueryType == 2) {
//            vis.m_query = queries[i];
//            tree->nearestNeighborQuery(thennk, *queries[i], vis);
////            cerr<<"finished "<<i<<"already\n";
//        }
//    }
//    double time;
//    auto end = std::chrono::system_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//    time = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
//    cerr << "Querying time: " << time << endl;
//    cerr << "VISIT NODE " << vis.m_indexvisited << "\t" << vis.m_leafvisited << endl;
//    cerr << "xStore Statistic" << ts->m_indexIO << "\t" << ts->m_trajIO << endl;
//}

static double kNNQueryBatch(xRTree *tree, const vector<xTrajectory> &queries, xStore *ts = nullptr, int thennk = 5,
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
        vis.m_query = (IShape *) &(queries[i]);
        tree->nearestNeighborQuery(thennk, queries[i], vis);
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
//    cerr <<"xStore Statistic"<< 1.0*ts->m_indexIO/num<<"\t"<<1.0*ts->m_trajIO/num<<endl;
//    sort(indios.begin(), indios.end());
//    int mid1 = queries.size() * 0.1, mid2 = queries.size() * 0.9;
//    double sum = 0;
//    for (int i = mid1; i <= mid2; i++) {
//        sum += indios[i];
//    }
//    sum = sum / (mid2 - mid1);
    cerr << "time\tindexVisit\tLeafVisit\t leaf1\tleaf2\tindexIO\ttrajIO\n";
    cerr << time / num << "\t" << 1.0 * vis.m_indexvisited / num << "\t" << 1.0 * vis.m_leafvisited / num << "\t"
         << 1.0 * ts->m_leaf1 / num << "\t" << 1.0 * ts->m_leaf2 / num << "\t" << 1.0 * ts->m_indexIO / num << "\t"
         << 1.0 * ts->m_trajIO / num << endl;
//    cerr <<time/num<<"\n";
    return rad;
}

static double kNNQueryBatch(xRP tree, const vector<xTrajectory> &queries, xStore *ts = nullptr, int thennk = 5,
                            bool reportEnd = false){
    return kNNQueryBatch(tree.get(), queries,ts,thennk,reportEnd);
}

static void rangeQueryBatch(xRTree *tree, const vector<xCylinder *> &queries, xStore *ts = nullptr, MyVisitor* vis = nullptr) {
    ts->cleanStatistic();
    int num = queries.size();
    if(vis == nullptr){
        vis = new MyVisitor();
    }
    vis->ts = ts;
    auto start = std::chrono::system_clock::now();
    for (int i = 0; i < queries.size(); i++) {
        vis->m_query = queries[i];
        tree->intersectsWithQuery(*queries[i], *vis);
    }
    double time;
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    time = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
//    cerr <<"Average Querying time: "<< time/num<<endl;
//    cerr <<"Averaged VISIT NODE "<<1.0*vis.m_indexvisited/num<<"\t"<<1.0*vis.m_leafvisited/num<<endl;
//    cerr <<"xStore Statistic"<< 1.0*ts->m_indexIO/num<<"\t"<<1.0*ts->m_trajIO/num<<endl;
    cerr << "average time\tIndexVisit\tLeafVisit\tIndexIO\ttrajIO\tprevalidateRate\tinternum\tcontainNum\n";
    cerr << time / num << "\t" << 1.0 * vis->m_indexvisited / num << "\t" << 1.0 * vis->m_leafvisited / num << "\t"
         << 1.0 * ts->m_indexIO / num << "\t" << 1.0 * ts->m_trajIO / num << "\t" << endl;
//    cerr <<vis->m_resultGet<<"\n";
//    cerr <<time/num<<"\n";
}

static void affine_transform(vector<pair<id_type, xTrajectory>> &ts, xPoint center, double angle, xPoint trans){
    for(auto &t:ts){
        for(auto &p:t.second.m_points){
            p.rotate(center,angle);
            p=p+trans;
        }
    }
}

// multithreaded




struct queryRet{
    double time=0;
    double indexVisit=0;
    double leafVisit=0;
    double leaf1 = 0;
    double leaf2 = 0;
    double indexIO=0;
    double trajIO=0;
    double nresult=0;
    double qps=0;
    queryRet operator+(const queryRet& r) const{
        queryRet res;
        res.time= time+r.time;
        res.indexVisit= indexVisit+r.indexVisit;
        res.leafVisit= leafVisit+r.leafVisit;
        res.leaf1= leaf1+r.leaf1;
        res.leaf2= leaf2+r.leaf2;
        res.indexIO= indexIO+r.indexIO;
        res.trajIO= trajIO+r.trajIO;
        res.nresult = nresult+r.nresult;
        res.qps = qps+r.qps;
        return res;
    }
    string toString() const{
        stringstream s;
        if(testxfirstOutput){
            testxfirstOutput =false;
            s<<"qps\ttime\tindexVisit\tleafVisit\tindexIO\ttrajIO\tleaf1\tleaf2\tnresult\ttotalIO\n";
        }
        s<<qps<<"\t"<<time<<"\t"<<indexVisit<<"\t"<<leafVisit<<"\t"<<indexIO<<"\t"<<trajIO<<"\t"<<leaf1<<"\t"<<leaf2<<"\t"<<nresult<<"\t"<<indexIO+trajIO<<"\n";
        return s.str();
    }
};


static queryRet average(vector<queryRet> &sum){
    queryRet res;
    for(auto &s:sum) res=res+s;
    res.time/=sum.size();
    res.indexVisit/=sum.size();
    res.leafVisit/=sum.size();
    res.leaf1/=sum.size();
    res.leaf2/=sum.size();
    res.indexIO/=sum.size();
    res.trajIO/=sum.size();
    res.nresult/=sum.size();
    return res;
}

enum queryType{
    qt_range,
    qt_knn
};

struct queryInput{
    xRTreeQueryObject * tree = nullptr;
    xStore * store;
    queryType type = qt_knn;
    vector<xTrajectory> knn_queries;
    vector<xCylinder> range_queries;
    int nnk =5;
};


static void QueryBatchThread(queryInput inp, queryRet *res) {
    xStore * ts = inp.store;
    ts->cleanStatistic();
    int num;
    MyVisitor vis;
    vis.ts = ts;
    auto start = std::chrono::system_clock::now();
    double rad = 0;
    if(inp.type == qt_knn) {
        num = inp.knn_queries.size();
        for (int i = 0; i < inp.knn_queries.size(); i++) {
            try {
//                drop_cache(1);
                vis.m_query = (IShape *) &(inp.knn_queries.at(i));
                inp.tree->nearestNeighborQuery(inp.nnk, inp.knn_queries.at(i), vis);
                rad += vis.m_lastDist;
            } catch (Tools::Exception &e) {
                cerr<<"error occurs at query \n "<<inp.knn_queries.at(i)<<endl;
                cerr << "******ERROR******" << endl;
                std::string s = e.what();
                cerr << s << endl;
            }
        }
    }else if(inp.type == qt_range){
        num = inp.range_queries.size();
        for (int i = 0; i < inp.range_queries.size(); i++) {
//            drop_cache(1);
            try {
                vis.m_query = (IShape *) &(inp.range_queries.at(i));
                inp.tree->intersectsWithQuery(inp.range_queries.at(i), vis);
                rad += vis.m_lastDist;
            } catch (Tools::Exception &e) {
                cerr<<"error occurs at query \n "<<inp.knn_queries.at(i)<<endl;
                cerr << "******ERROR******" << endl;
                std::string s = e.what();
                cerr << s << endl;
            }
        }
    }
    double time;
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    time = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
    res->time = time/num;
    res->qps = num/time;
    res->indexVisit = 1.0*vis.m_indexvisited/num;
    res->leafVisit = 1.0 * vis.m_leafvisited/ num;
    res->leaf1 = 1.0 * ts->m_leaf1/num;
    res->leaf2 = 1.0 * ts->m_leaf2 / num;
    res->indexIO = 1.0 * ts->m_indexIO / num;
    res->trajIO = 1.0 * ts->m_trajIO / num;
    res->nresult = 1.0 * vis.m_resultGet / num;
    dbt.release();
    return;
}

/**
 * multi-threaded queries
 */
class MTQ{
public:
    vector<xStore*> m_stores;
    vector<xRTreeQueryObject*> m_trees;
    vector<queryInput> m_queries;
    vector<queryRet> m_res;
    vector<thread> m_ths;
    int nthread=NUMTHREAD;
    ~MTQ(){
        for(auto &a:m_trees){
            delete a;
        }
        for(auto &a:m_stores){
            delete a;
        }
    }
    void prepareTrees(xStore* x,
                      const function<xRTree*(IStorageManager*)> &treeBuilder){
        delete treeBuilder(x);
        for(int i=0;i<nthread;i++) {
            m_stores.emplace_back(x->clone());
            m_trees.emplace_back(treeBuilder(m_stores.back()));
            queryInput q;
            q.tree = m_trees.back();
            q.store=m_stores.back();
            m_queries.emplace_back(q);
            m_res.emplace_back(queryRet());
        }
    }
    void prepareForest(xStore* x,map<pair<double, double>, double> &lens, double slab=1e300){
        delete buildSBBForest(x,xTrajectory::OPTS,lens,slab);
        for(int i=0;i<nthread;i++) {
            m_stores.emplace_back(x->clone());
            m_trees.emplace_back(buildSBBForest(m_stores.back(),xTrajectory::OPTS,lens,slab));
            queryInput q;
            q.tree = m_trees.back();
            q.store=m_stores.back();
            m_queries.emplace_back(q);
            m_res.emplace_back(queryRet());
        }
    }
    void appendQueries(vector<xTrajectory> &knnq, int nnk=6){
        int i=0;
        for(auto &q:knnq)
        {
            m_queries[i % nthread].knn_queries.emplace_back(q);
            i++;
        }
        for(auto &qs:m_queries){
            qs.nnk = nnk;
        }
    }

    void appendQueries(vector<xTrajectory> &knnq, vector<int> &nnk){
        int i=0;
        for(auto &q:knnq)
        {
            m_queries[i % nthread].knn_queries.emplace_back(q);
            i++;
        }
        for(int i=0;i<m_queries.size();i++){
            m_queries[i].nnk=nnk[i];
        }
    }


    void appendQueries(vector<xCylinder> &rq){
        int i=0;
        for(auto &q:rq)
        {
            m_queries[i % nthread].range_queries.emplace_back(q);
            i++;
        }
        for(auto &qs:m_queries){
            qs.type=qt_range;
        }
    }
    queryRet runQueries(){
        drop_cache(3);
        usleep(5000);
        for(int i=0;i<nthread;i++) {
            m_ths.emplace_back(thread(QueryBatchThread, m_queries[i], &m_res[i]));
        }
        for(int i=0;i<nthread;i++) {
            m_ths[i].join();
        }
        return average(m_res);
    }
};

string testFileName(string s){
#if (defined _WIN32 || defined _WIN64 || defined WIN32 || defined WIN64)
    return "D://TRI-framework/dumpedtraj.txt";
#else
    return "/root/"+s;
#endif
}

#include "random"
void fillQuerySet(vector<xTrajectory>& list, xStore& x, double len, double var=0, int num=testtime){
    if(var!=0){
        default_random_engine e;
        auto queryLen =normal_distribution<double>(len,var);
        for (int i = 0; i < num; i++) {
            list.emplace_back(x.randomSubtraj(queryLen(e)));
        }
    }
    else{
        for (int i = 0; i < num; i++) {
            list.emplace_back(x.randomSubtraj(len));
        }
    }
    cerr<<"sampled queries from "<<db_last_trajid()<<"trajs\n";
}

void fillQuerySetRand(vector<xTrajectory>& list, xStore& x){
    default_random_engine e;
    auto queryLen =uniform_real_distribution<double>(300,10800);
    for (int i = 0; i < testtime; i++) {
        list.emplace_back(x.randomSubtraj(queryLen(e)));
    }
}

void fillQuerySet(vector<xCylinder>& list, xStore& x, double rd,double qt, double var=0, int num=testtime){
    if(var!=0){
        default_random_engine e;
        auto rdd =normal_distribution<double>(rd,var);
        for (int i = 0; i < num; i++) {
            xPoint p =x.randomPoint();
            list.emplace_back(xCylinder(p,rdd(e),p.m_t-qt/2,p.m_t+qt/2));
        }
    }
    else{
        for (int i = 0; i < num; i++) {
            xPoint p =x.randomPoint();
            list.emplace_back(xCylinder(p,rd,p.m_t-qt/2,p.m_t+qt/2));
        }
    }
}

#endif //SPATIALINDEX_TESTFUNCS_H
