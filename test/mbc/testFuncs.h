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
#define random(x,y) (((double)rand()/RAND_MAX)*(y-x)+x)
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
//1 for time-slice range, 2 for 5-NN

using namespace std;
using namespace SpatialIndex;

class MyVisitor : public IVisitor
{
public:
    size_t m_indexIO;
    size_t m_leafIO;
    size_t m_indexvisited;
    size_t m_leafvisited;
    size_t m_resultGet;
    id_type m_lastResult;
    IShape *m_query;
    TrajStore *ts=nullptr;

public:
    MyVisitor() : m_indexIO(0), m_leafIO(0),m_resultGet(0),m_indexvisited(0),m_leafvisited(0) {}

    void visitNode(const INode& n)
    {
//        if (n.isLeaf()) m_leafIO++;
//        else m_indexIO++;
        uint32_t size=n.getIndexByteArraySize();

        if (n.isLeaf()) {m_leafvisited++;m_leafIO+=size;}
        else {m_indexvisited++;m_indexIO+=size;}
    }

    void visitData(const IData& d)
    {
        m_resultGet++;
        IShape* pS;
        d.getShape(&pS);
        // do something.
        delete pS;
//        cout<<"data"<<endl;

        // data should be an array of characters representing a Region as a string.
        uint8_t* pData = 0;
        uint32_t cLen = 0;
        d.getData(cLen, &pData);
        m_leafIO+=cLen;
        // do something.
//        double *s = reinterpret_cast<double*>(pData);
//        cout << *s << endl;

        m_lastResult=d.getIdentifier();
//        auto mou=dynamic_cast<const MBCRTree::MBCRTree::simpleData*>(&d);
//        cerr << d.getIdentifier()<<"\t"<<mou->m_dist << endl;
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

    void visitData(std::vector<const IData*>& v)
    {
        cout << v[0]->getIdentifier() << " " << v[1]->getIdentifier() << endl;
    }
};


class TrajMbrStream: public IDataStream{
public:
    vector<pair<id_type,Region> > mbrs;
    vector<pair<id_type,Trajectory> > *trajs;

    int i=0;
    void feedTraj(vector<pair<id_type ,Trajectory> > *period){
        cerr<<"feeding traj to TrajMbrStream\n";
        mbrs.clear();
        for(auto &idt:*period){
            Region mbr;
            idt.second.getMBRfull(mbr);
            mbrs.emplace_back(make_pair(idt.first,mbr));
        }
        trajs=period;
        rewind();
    }
    virtual bool hasNext() override
    {
        return i<mbrs.size();
    }
    virtual IData* getNext() override{
        uint8_t *data;
        uint32_t len;
        trajs->at(i).second.storeToByteArray(&data,len);
        RTree::Data* d=new RTree::Data(len, data, mbrs[i].second, mbrs[i].first);
        i++;
        return d;
    }
    virtual uint32_t size()
    {
        return mbrs.size();
    }

    virtual void rewind(){i=0;}
};

struct xyt{
    double x;
    double y;
    double t;
};
xyt makemid(xyt p1, xyt p2, double t){
    if(t>p2.t)
        cout<<p1.x<<" "<<p2.x<<endl<<
            p1.y<<" "<<p2.y<<endl<<
            p1.t<<" "<<p2.t<<" "<<t<<endl;
    assert(p1.t<=t);
    assert(t<=p2.t);
    double h1= (t-p1.t)/(p2.t-p1.t);
    double h2= (p2.t-t)/(p2.t-p1.t);
    double x=h2*p1.x+h1*p2.x;
    double y=h2*p1.y+h1*p2.y;
    xyt ret ={x,y,t};
    return ret;
}

template <class Type>
Type stringToNum(const std::string& str)
{
    std::istringstream iss(str);
    Type num;
    iss >> num;
    return num;
}
vector<pair<id_type ,Trajectory> >  loadGTToTrajs(string filename=genFile){
    //first level: vector of time period
    //second level: vector of segments in the time period
#ifndef NDEBUG
    cerr<<"loading generated trajectories from txt to trajectories"<<endl;
#endif
    auto stat = trajStat::instance();
    ifstream inFile(filename, ios::in);
    string lineStr;
    set<id_type> ids;
    multimap<id_type,xyt> trajs;
    vector<pair<id_type ,Trajectory> > res;
    int curLine=0;
    while (getline(inFile, lineStr)&&curLine<maxLinesToRead){
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
        catch(...) {
            break;
        }
    }
    stat->Dx=stat->maxx-stat->minx;
    stat->Dy=stat->maxy-stat->miny;
    stat->Dt=stat->maxt-stat->mint;
    for(auto id:ids){
        multimap<id_type ,xyt>::iterator beg,end,iter;
        vector<xyt> traj;
        beg = trajs.lower_bound(id);
        end = trajs.upper_bound(id);
        for(iter=beg;iter!=end;iter++){
            traj.emplace_back(iter->second);
        }
        trajs.erase(id);
        if(traj.size()>=2){
            vector<STPoint> tps;
            xyt lastpoint;
            for(int l=0;l<traj.size();l++){
                tps.emplace_back(STPoint(traj[l].x, traj[l].y,traj[l].t));
                if(l!=0){
                    stat->dist+=std::sqrt(sq(traj[l].x-lastpoint.x)+sq(traj[l].y-lastpoint.y));
                }
                lastpoint=traj[l];
            }
            if(!tps.empty()){
                Trajectory tmp(tps);
                stat->lineCount+=tps.size()-1;
                stat->trajCount+=1;
                stat->M+=tmp.m_endTime()-tmp.m_startTime();
                res.emplace_back(make_pair(id,tmp));
            }
        }
    }
#ifndef NDEBUG
    std::cerr<<"load data finished\n";
#endif
    stat->tl=stat->M/stat->lineCount;
    stat->jt=stat->M/stat->trajCount;
    stat->v=stat->dist/stat->M;
    std::cerr<<*stat;
    return res;
}

vector<pair<id_type ,Trajectory> >  loadGLToTrajs(string filename=GLFile){
    //first level: vector of time period
    //second level: vector of segments in the time period
    cerr<<"loading geolife trajectories from txt to trajectories"<<endl;
    auto stat = trajStat::instance();
    ifstream inFile(filename, ios::in);
    string lineStr;
    set<id_type> ids;
    multimap<id_type,xyt> trajs;
    vector<pair<id_type ,Trajectory> > res;
    int curLine=0;
    getline(inFile, lineStr);
    while (getline(inFile, lineStr)&&curLine<maxLinesToRead){
        try {
            string str;
            stringstream ss(lineStr);
            getline(ss, str, ',');
            int id= stringToNum<int>(str);
            getline(ss, str, ',');
            double x= stringToNum<double>(str);
            getline(ss, str, ',');
            double y= stringToNum<double>(str);
            getline(ss, str, ',');
            double t= stringToNum<double>(str);
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
        catch(...) {
            break;
        }
    }
    stat->Dx=stat->maxx-stat->minx;
    stat->Dy=stat->maxy-stat->miny;
    stat->Dt=stat->maxt-stat->mint;
    for(auto id:ids){
        multimap<id_type ,xyt>::iterator beg,end,iter;
        vector<xyt> traj;
        beg = trajs.lower_bound(id);
        end = trajs.upper_bound(id);
        for(iter=beg;iter!=end;iter++){
            traj.emplace_back(iter->second);
        }
        trajs.erase(id);
        if(traj.size()>=2){
            vector<STPoint> tps;
            xyt lastpoint;
            for(int l=0;l<traj.size();l++){
                tps.emplace_back(STPoint(traj[l].x, traj[l].y,traj[l].t));
                if(l!=0){
                    stat->dist+=std::sqrt(sq(traj[l].x-lastpoint.x)+sq(traj[l].y-lastpoint.y));
                }
                lastpoint=traj[l];
            }
            if(!tps.empty()){
                Trajectory tmp(tps);
                stat->lineCount+=tps.size()-1;
                stat->trajCount+=1;
                stat->M+=tmp.m_endTime()-tmp.m_startTime();
                res.emplace_back(make_pair(id,tmp));
            }
        }
    }
#ifndef NDEBUG
    std::cerr<<"load data finished\n";
#endif
    stat->tl=stat->M/stat->lineCount;
    stat->jt=stat->M/stat->trajCount;
    stat->v=stat->dist/stat->M;
    std::cerr<<*stat;
    return res;
}

double knncost(double bt,int k,double qt,int f,bool useMBR,double _rk){
    auto stat=trajStat::instance();
    double nq=std::min(double(stat->trajCount),stat->M/stat->jt*(qt+2*stat->jt)/stat->Dt);
    double rk=sqrt(k*stat->Dx*stat->Dy/M_PI/nq);
    if(_rk>0) rk=_rk;
    double d=bt*stat->v;
    double lmd;
    if(useMBR) lmd=1+4*d/M_PI/M_PI/rk+d*d/rk/rk/M_PI/M_PI;
    else lmd=sq(1+0.15*d/rk);
    double ltc=pow(stat->Dx*stat->Dy*stat->Dt/stat->M*bt*f/stat->v/stat->v,1.0/3);
    double lt=ltc+bt;
    double lx=lt*stat->v;
    double leafnum=stat->M/bt/f;
    double interleafnum=leafnum*lt/stat->Dt;
    double hc=1.1*(1+qt/lt)*sq(2*rk+lx)*interleafnum/stat->Dx/stat->Dy,
            lc=k*(1+qt/bt)*lmd*std::ceil(bt/stat->tl/160);
    //std::cerr<<bt<<"\t"<<hc<<"\t"<<lc<<"\n";
    return hc+lc;
}
struct d4{
    d4(double _low,double _high,double _vlow,double _vhigh){
        low=_low,high=_high,vlow=_vlow,vhigh=_vhigh;
    }
    double low,high,vlow,vhigh;
};
double biSearchMax(int k,double qt,int f,bool useMBR, double rk=-1){
    double low=1,high=10000;
    std::stack<d4> st;
    double vlow,vhigh;
    vlow=knncost(low,k,qt,f,useMBR,rk);
    vhigh=knncost(high,k,qt,f,useMBR,rk);
    st.push(d4(low,high,vlow,vhigh));
    d4 first(0,0,0,0);
    auto stat=trajStat::instance();
    while(!st.empty()){
        first=st.top();
        st.pop();
        double mid=(first.low+first.high)/2;
        if(first.high-first.low<1){
            return std::max(mid,stat->tl);
        }
        double vmid=knncost(mid,k,qt,f,useMBR,rk);
        bool rising=knncost(mid+0.1,k,qt,f,useMBR,rk)>vmid;
//        std::cerr<<mid<<"\t"<<vmid<<"\n";
        if(rising){
            st.push(d4(first.low,mid,first.vlow,vmid));
        }else{
            st.push(d4(mid,first.high,vmid,first.vhigh));
        }
    }
    throw Tools::IllegalStateException("biSearch:failed searching");
}

void TreeQueryBatch(ISpatialIndex* tree,const vector<IShape*> &queries,TrajStore *ts= nullptr,int thennk=5){
    MyVisitor vis;
    vis.ts=ts;
    auto start = std::chrono::system_clock::now();
    for(int i=0;i<queries.size();i++){
//        cerr<<"Query is "<<queries.at(i)->toString();
        if(QueryType==1){
            tree->intersectsWithQuery(*queries[i],vis);
        }else if(QueryType==2){
            vis.m_query=queries[i];
            tree->nearestNeighborQuery(thennk,*queries[i],vis);
//            cerr<<"finished "<<i<<"already\n";
        }
    }
    double time;
    auto end = std::chrono::system_clock::now();auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    time=double(duration.count()) * std::chrono::microseconds::period::num/ std::chrono::microseconds::period::den;
    cerr <<"Querying time: "<< time<<endl;
    cerr <<"VISIT NODE "<<vis.m_indexvisited<<"\t"<<vis.m_leafvisited<<endl;
    cerr <<"TrajStore Statistic"<< ts->m_indexIO<<"\t"<<ts->m_trajIO<<endl;
}

void kNNQueryBatch(ISpatialIndex* tree,const vector<IShape*> &queries,TrajStore *ts= nullptr,int thennk=5){
    ts->cleanStatistic();
    int num=queries.size();
    MyVisitor vis;
    vis.ts=ts;
    auto start = std::chrono::system_clock::now();
    for(int i=0;i<queries.size();i++){
        vis.m_query=queries[i];
        tree->nearestNeighborQuery(thennk,*queries[i],vis);
    }
    double time;
    auto end = std::chrono::system_clock::now();auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    time=double(duration.count()) * std::chrono::microseconds::period::num/ std::chrono::microseconds::period::den;
//    cerr <<"Average Querying time: "<< time/num<<endl;
//    cerr <<"Averaged VISIT NODE "<<1.0*vis.m_indexvisited/num<<"\t"<<1.0*vis.m_leafvisited/num<<endl;
//    cerr <<"TrajStore Statistic"<< 1.0*ts->m_indexIO/num<<"\t"<<1.0*ts->m_trajIO/num<<endl;
    cerr <<time/num<<"\t"<<1.0*vis.m_indexvisited/num<<"\t"<<1.0*vis.m_leafvisited/num<<"\t"<<1.0*ts->m_leaf1/num<<"\t"<<1.0*ts->m_leaf2/num<<"\t"<< 1.0*ts->m_indexIO/num<<"\t"<<1.0*ts->m_trajIO/num<<"\t"<<1.0*ts->m_loadedTraj/num<<endl;
//    cerr <<time/num<<"\n";
}
void rangeQueryBatch(ISpatialIndex* tree,const vector<IShape*> &queries,TrajStore *ts= nullptr,int thennk=5){
    ts->cleanStatistic();
    int num=queries.size();
    MyVisitor vis;
    vis.ts=ts;
    auto start = std::chrono::system_clock::now();
    for(int i=0;i<queries.size();i++){
        vis.m_query=queries[i];
        tree->intersectsWithQuery(*queries[i],vis);
    }
    double time;
    auto end = std::chrono::system_clock::now();auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    time=double(duration.count()) * std::chrono::microseconds::period::num/ std::chrono::microseconds::period::den;
//    cerr <<"Average Querying time: "<< time/num<<endl;
//    cerr <<"Averaged VISIT NODE "<<1.0*vis.m_indexvisited/num<<"\t"<<1.0*vis.m_leafvisited/num<<endl;
//    cerr <<"TrajStore Statistic"<< 1.0*ts->m_indexIO/num<<"\t"<<1.0*ts->m_trajIO/num<<endl;
    cerr <<time/num<<"\t"<<1.0*vis.m_indexvisited/num<<"\t"<<1.0*vis.m_leafvisited/num<<"\t"<< 1.0*ts->m_indexIO/num<<"\t"<<1.0*ts->m_trajIO/num<<endl;
//    cerr <<time/num<<"\n";
}

int TreeQuery(ISpatialIndex* tree,IShape* query,TrajStore *ts= nullptr){
    clock_t start,end;
    MyVisitor vis;
    if(ts!= nullptr)
        vis.ts=ts;
    vis.m_query=query;
    start=clock();
    if(QueryType==1){
        tree->intersectsWithQuery(*query,vis);
    }else if(QueryType==2){
        vis.m_query=query;
        tree->nearestNeighborQuery(5,*query,vis);
    }
    end=clock();
    if(QueryType==1){
        return vis.m_resultGet;
    }else if(ts!= nullptr){
        return vis.m_lastResult;
    }
    else return vis.m_lastResult;
}




#endif //SPATIALINDEX_TESTFUNCS_H
