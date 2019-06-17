//
// Created by Chuang on 2019/5/21.
//

//
// Created by chuang on 4/9/19.
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
#define random(x,y) (((double)rand()/RAND_MAX)*(y-x)+x)
#include <spatialindex/SpatialIndex.h>
#include "storagemanager/TrajStore.h"
//#define sourceFile "D://t200n100s.txt"
#define sourceFile "D://t1000.txt"
#define maxLinesToRead 1e10
#define testtime 1000
#define dimension 2
#define indexcap 10
#define leafcap 10000
#define QueryType 1
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

public:
    MyVisitor() : m_indexIO(0), m_leafIO(0),m_resultGet(0),m_indexvisited(0),m_leafvisited(0) {}

    void visitNode(const INode& n)
    {
//        if (n.isLeaf()) m_leafIO++;
//        else m_indexIO++;
        uint32_t size=n.getIndexByteArraySize();

        if (n.isLeaf()) {m_indexvisited++;m_indexIO+=size;}
        else {m_leafvisited++;m_indexIO+=size;}
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
//        cout << d.getIdentifier() << endl;
////         the ID of this data entry is an answer to the query. I will just print it to stdout.
//        Trajectory traj;
//        traj.loadFromByteArray(pData);
//        Region br;
//        Mbbc bc;
//        traj.getMBR(br);
//        traj.getMbbc(bc);
//        double mindist=m_query->getMinimumDistance(traj);
//        cout<<"traj dist is"<<mindist<<"\t";
//        cout<<"br loose"<<mindist-m_query->getMinimumDistance(br)<<"\t";
//        cout<<"bc loose is"<<mindist-m_query->getMinimumDistance(bc)<<"\n";
//        if(mindist-m_query->getMinimumDistance(bc)<0){
//            std::cerr<<m_query->toString()<<traj.toString();
//            system("pause");
//        }
//        cout<<traj.toString()<<endl;
//        double pLow[2]={39.993017,116.320135};
//        double pHigh[2]={39.994017,116.321135};
//        traj.intersectsShape(TimeRegion(pLow,pHigh,855,855,2));
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
class TrajMbcStream: public IDataStream{
public:
    vector<pair<id_type,MBC> > mbcs;
    vector<pair<id_type,Trajectory> > *trajs;

    int i=0;
    void feedTraj(vector<pair<id_type ,Trajectory> > *period){
        cerr<<"feeding traj to TrajMbcStream\n";
        mbcs.clear();
        MBC mbc;
        for(auto idt:*period){
            idt.second.getMBC(mbc);
            mbcs.emplace_back(make_pair(idt.first,mbc));
        }
        trajs=period;
        rewind();
    }
    virtual bool hasNext() override
    {
        return i<mbcs.size();
    }
    virtual IData* getNext() override{
        uint8_t *data;
        uint32_t len;
        trajs->at(i).second.storeToByteArray(&data,len);
        MBCRTree::Data* d=new MBCRTree::Data(len, data, mbcs[i].second, mbcs[i].first);
        i++;
        return d;
    }
    virtual uint32_t size()
    {
        return mbcs.size();
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

vector<pair<id_type ,Trajectory> >  loadGTToTrajs(){
    //first level: vector of time period
    //second level: vector of segments in the time period
    cerr<<"loading generated trajectories from txt to trajectories"<<endl;

    ifstream inFile(sourceFile, ios::in);
    string lineStr;
    set<id_type> ids;
    multimap<id_type,xyt> trajs;
    vector<pair<id_type ,Trajectory> > res;
    int curLine=0;
    double minx=40000,maxx=0,miny=40000,maxy=0;
    while (getline(inFile, lineStr)&&curLine<maxLinesToRead){
        string str;
        stringstream ss(lineStr);
        getline(ss, str, '\t');
        getline(ss, str, '\t');
        int id= stringToNum<int>(str);
        getline(ss, str, '\t');
        getline(ss, str, '\t');
        getline(ss, str, '\t');
        double t= stringToNum<double>(str);
        getline(ss, str, '\t');
        double x= stringToNum<double>(str);
        getline(ss, str, '\t');
        double y= stringToNum<double>(str);
        getline(ss, str, '\t');
        double speed=stringToNum<double>(str);
        xyt p={x,y,t};
        if(x>maxx) maxx=x;
        if(x<minx) minx=x;
        if(y>maxy) maxy=y;
        if(y<miny) miny=y;
        ids.insert(id);
        trajs.insert(make_pair(id,p));
        curLine++;
    }
    cout<<minx<<" "<<maxx<<" "<<miny<<" "<<maxy<<endl;
    for(auto id:ids){
        multimap<id_type ,xyt>::iterator beg,end,iter;
        vector<xyt> traj;
        beg = trajs.lower_bound(id);
        end = trajs.upper_bound(id);
        for(iter=beg;iter!=end;iter++){
            traj.emplace_back(iter->second);
        }
        trajs.erase(id);
        if(traj.size()>=10){
            vector<TimePoint> tps;
            for(auto p:traj){
                double xy[]={p.x,p.y};
                tps.emplace_back(TimePoint(xy, p.t, p.t, dimension));
            }
            if(!tps.empty()){
                res.emplace_back(make_pair(id,Trajectory(tps)));
            }
        }
    }
    return res;
}

void TreeQueryBatch(ISpatialIndex* tree,const vector<IShape*> &queries){
    clock_t start,end;

    MyVisitor vis;
    start=clock();
    for(int i=0;i<queries.size();i++){
//        cerr<<"Query is "<<queries.at(i)->toString();
        if(QueryType==1){
            tree->intersectsWithQuery(*queries[i],vis);
        }else if(QueryType==2){
            tree->nearestNeighborQuery(5,*queries[i],vis);
        }
//        cout<<"finished "<<i<<"already\n";
    }
    end=clock();
    cerr<< "=================\n\n";
    cerr<<"Querying time: "<< end-start<<endl;
    cerr<<"VISIT NODE "<<vis.m_indexvisited<<"\t"<<vis.m_leafvisited<<endl;
    cerr<<"Byte loaded "<<vis.m_indexIO<<"\t"<<vis.m_leafIO<<endl;
    cerr << *tree;
    cerr << "================\n\n";
}

int TreeQuery(ISpatialIndex* tree,IShape* query){
    clock_t start,end;
    MyVisitor vis;
    vis.m_query=query;
    start=clock();
    if(QueryType==1){
        tree->intersectsWithQuery(*query,vis);
    }else if(QueryType==2){
        tree->nearestNeighborQuery(5,*query,vis);
    }
    end=clock();
    return vis.m_indexIO;
}


int main(){
    try {
        srand((int) time(NULL));
        vector<pair<id_type, Trajectory> > trajs = loadGTToTrajs();
//        TrajMbrStream ds1;
//        TrajMbcStream ds2;
//        ds1.feedTraj(&trajs);
//        ds2.feedTraj(&trajs);

        vector<pair<id_type, vector<Trajectory>>> segs;
        vector<pair<id_type, Trajectory> > empty1;
        for(auto traj:trajs){
            segs.push_back(make_pair(traj.first,traj.second.getSegments(3000)));
        }
        trajs.swap(empty1);
        vector<IShape *> queries;
//        double pLow[]={19476.912748,26935.636464,913.000000};
//        double pHigh[]={20982.272042,27436.950774,913.000000};
//        Region *rg = new Region(pLow, pHigh, 3);
//        queries.emplace_back(rg);
        for (int i = 0; i < testtime; i++) {
            if (QueryType == 1) {
                double t = int(random(0, PeriodLen));
                double pLow[3] = {random(0, 25000), random(0, 30000), t};
                double pHigh[3] = {pLow[0] + random(1, 2000), pLow[1] + random(1, 2000), t};
                Region *rg = new Region(pLow, pHigh, 3);
                queries.emplace_back(rg);
            }
        }

        string name0 = "name0", name1 = "name1", name2 = "name2", name3 = "name3", name4 = "name4",
                name5 = "name5", name6 = "name6", name7 = "name7";
        id_type indexIdentifier0, indexIdentifier1, indexIdentifier2, indexIdentifier3, indexIdentifier4,
                indexIdentifier5, indexIdentifier6, indexIdentifier7;
        IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096),
                *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
                *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096),
                *diskfile3 = StorageManager::createNewDiskStorageManager(name3, 4096),
                *diskfile4 = StorageManager::createNewDiskStorageManager(name4, 4096),
                *diskfile5 = StorageManager::createNewDiskStorageManager(name5, 4096),
                *diskfile6 = StorageManager::createNewDiskStorageManager(name6, 4096),
                *diskfile7 = StorageManager::createNewDiskStorageManager(name7, 4096);
        // Create a new storage manager with the provided base name and a 4K page size.
        StorageManager::IBuffer *file0 = StorageManager::createNewRandomEvictionsBuffer(*diskfile0, 10, false),
                *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false),
                *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false),
                *file3 = StorageManager::createNewRandomEvictionsBuffer(*diskfile3, 10, false),
                *file4 = StorageManager::createNewRandomEvictionsBuffer(*diskfile4, 10, false),
                *file5 = StorageManager::createNewRandomEvictionsBuffer(*diskfile5, 10, false),
                *file6 = StorageManager::createNewRandomEvictionsBuffer(*diskfile6, 10, false),
                *file7 = StorageManager::createNewRandomEvictionsBuffer(*diskfile7, 10, false);
        TrajStore ts1(file1,4096);
        ts1.loadSegments(segs);
        TrajStore ts2(file2,4096);
        ts2.loadSegments(segs);
        ISpatialIndex *r=RTree::createAndBulkLoadNewRTreeWithTrajStore(&ts1,5,3,indexIdentifier1);
        ISpatialIndex *rc=MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(&ts2,5,3,indexIdentifier2);

        cerr << "start query!" << endl << endl << endl;
        TreeQueryBatch(r, queries);
        TreeQueryBatch(rc, queries);
    }
    catch (Tools::Exception& e)
    {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
}