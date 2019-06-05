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
//#define sourceFile "D://geolifedatasimplify.csv"
//#define sourceFile "D://geolifedata.csv"
#define sourceFile "D://t200n10ks.txt"
#define linesToRead 1e10
#define testtime 1000
#define dimension 2
#define indexcap 5
#define leafcap 5
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
    MyVisitor() : m_indexIO(0), m_leafIO(0),m_resultGet(0) {}

    void visitNode(const INode& n)
    {
//        if (n.isLeaf()) m_leafIO++;
//        else m_indexIO++;
        uint32_t size=n.getByteArraySize();

        if (n.isLeaf()) {m_indexvisited++;m_leafIO+=size;}
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
        for(auto idt:*period){
            Region mbr;
            idt.second.getMBR(mbr);
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

class TrajMBRkStream: public IDataStream{
public:
    vector<pair<id_type,MBRk> > mbrks;
    vector<pair<id_type,Trajectory> > *trajs;

    int i=0;
    void feedTraj(vector<pair<id_type ,Trajectory> > *period,int k){
        cerr<<"feeding traj to TrajMBR"<<k<<" Stream\n";
        mbrks.clear();
        for(auto idt:*period){
            MBRk brk(k);
//            cout<<idt.first<<endl<<idt.second.toString()<<endl;
            idt.second.getMBRk(k,brk);

            mbrks.emplace_back(make_pair(idt.first,brk));
        }
        trajs=period;
        rewind();
        cerr<<"feeding traj to TrajMBR"<<k<<" Stream finished\n";
    }
    virtual bool hasNext() override
    {
        return i<mbrks.size();
    }
    virtual IData* getNext() override{
        uint8_t* data;
        uint32_t len;
        trajs->at(i).second.storeToByteArray(&data,len);
        PAARTree::Data* d=new PAARTree::Data(len, data, mbrks[i].second, mbrks[i].first);
        i++;
        return d;
    }
    virtual uint32_t size()
    {
        return mbrks.size();
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
vector< vector<xyt> > cuttraj(vector<xyt> traj){
    vector< vector<xyt> > segments(getMaxPeriod());
    int oldpd=getPeriod(traj.at(0).t);
    for(int i=0;i<traj.size();i++){
        if(i>1&&traj[i].t<traj[i-1].t) break;//stop reading when coming to a new day
        int newpd=getPeriod(traj[i].t);
        if(newpd-1==oldpd){
            xyt mid1=makemid(traj[i-1],traj[i],getPeriodEnd(traj[i-1].t));
            segments.at(oldpd).emplace_back(mid1);
            if(traj[i].t-getPeriodStart(traj[i].t)>=0.1){
                xyt mid2=makemid(traj[i-1],traj[i],getPeriodStart(traj[i].t));
                segments.at(newpd).emplace_back(mid2);
            }
            oldpd=newpd;
        }
        segments.at(newpd).emplace_back(traj[i]);
//        cout<<traj[i].t<<"\n";
    }
    return segments;
}



list<vector<pair<id_type ,Trajectory> > > loadGLToTrajs(){
    //first level: vector of time period
    //second level: vector of segments in the time period
    cerr<<"loading csv to trajectories"<<endl;
    ifstream inFile(sourceFile, ios::in);
    string lineStr;
    getline(inFile, lineStr);
    set<id_type> ids;
    multimap<id_type,xyt> trajs;
    list<vector<pair<id_type ,Trajectory> > > res(getMaxPeriod());
    int curLine=0;
    while (getline(inFile, lineStr)&&curLine<linesToRead){
        string str;
        stringstream ss(lineStr);
        getline(ss, str, ',');
        int id= stringToNum<int>(str);
        getline(ss, str, ',');
        double x= stringToNum<double>(str);
        getline(ss, str, ',');
        double y= stringToNum<double>(str);
        getline(ss, str, ',');
        getline(ss, str, ',');
        double t=naivetime(str);
        xyt p={x,y,t};
        ids.insert(id);
        trajs.insert(make_pair(id,p));
        curLine++;
    }
    for(auto id:ids){
        multimap<id_type ,xyt>::iterator beg,end,iter;
        vector<xyt> traj;
        beg = trajs.lower_bound(id);
        end = trajs.upper_bound(id);
        for(iter=beg;iter!=end;iter++){
            traj.emplace_back(iter->second);
        }
        trajs.erase(id);
        if(!traj.empty()){
//            cout<<id<<endl;
            vector< vector<xyt> > segs = cuttraj(traj);
            auto iperiod=res.begin();
            for(int j =0;j<getMaxPeriod();j++){
                vector<TimePoint> tps;
                for(auto p:segs[j]){
                    double xy[]={p.x,p.y};
                    double faket=p.t-getPeriodStart(p.t);
                    tps.emplace_back(TimePoint(xy, faket, faket, dimension));
                }
                if(!tps.empty()){
                    iperiod->emplace_back(make_pair(id*1000+j,Trajectory(tps)));
//                    iperiod++;
                }
            }
        }
    }
    return res;
}
list<vector<pair<id_type ,Trajectory> > > loadGTToTrajs(){
    //first level: vector of time period
    //second level: vector of segments in the time period
    cerr<<"loading generated trajectories from txt to trajectories"<<endl;

    ifstream inFile(sourceFile, ios::in);
    string lineStr;
    set<id_type> ids;
    multimap<id_type,xyt> trajs;
    list<vector<pair<id_type ,Trajectory> > > res(getMaxPeriod());
    int curLine=0;
    double minx=40000,maxx=0,miny=40000,maxy=0;
    while (getline(inFile, lineStr)&&curLine<linesToRead){
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
//        cout<<id<<","<<x<<","<<y<<","<<t<<endl;
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
        if(!traj.empty()){
            vector< vector<xyt> > segs = cuttraj(traj);
            auto iperiod=res.begin();
            for(int j =0;j<getMaxPeriod();j++){
                vector<TimePoint> tps;
                for(auto p:segs[j]){
                    double xy[]={p.x,p.y};
                    double faket=p.t-getPeriodStart(p.t);
                    tps.emplace_back(TimePoint(xy, faket, faket, dimension));
                }
                if(!tps.empty()){
                    iperiod->emplace_back(make_pair(id*1000+j,Trajectory(tps)));
//                    iperiod++;
                }
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
    cerr<<"Querying time: "<< end-start<<endl;
    cerr<<"VISIT NODE "<<vis.m_indexvisited<<"\t"<<vis.m_leafvisited<<endl;
    cerr<<"Byte loaded "<<vis.m_indexIO<<"\t"<<vis.m_leafIO<<endl;
    cerr << *tree;
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
//    cerr<<"Querying time: "<< end-start<<endl;
//    cerr<<"VISIT NODE "<<vis.m_indexIO<<","<<vis.m_leafIO<<endl;
//    cerr << *tree;
//    return vis.m_resultGet;
    return vis.m_indexIO;
}


int main(){
    try {
        srand((int) time(NULL));
        list<vector<pair<id_type, Trajectory> > > trajs = loadGTToTrajs();
        auto traj1=*trajs.begin();
        cout<<traj1[0].second.toString();
//        TrajMbrStream ds1;
//        TrajMbbcStream ds2;
        int a=8,b=20,c=25,d=50,e=100;
        TrajMBRkStream ds3;
        TrajMBRkStream ds34;
        TrajMBRkStream ds38;
        TrajMBRkStream ds316;
        TrajMBRkStream ds332;
        cout<<"total trajectories:"<<trajs.front().size()<<endl;
//        ds1.feedTraj(&traj1);
//        ds2.feedTraj(&trajs.front());
        ds3.feedTraj(&traj1, a);
        ds34.feedTraj(&traj1, b);
        ds38.feedTraj(&traj1, c);
        ds316.feedTraj(&traj1, d);
        ds332.feedTraj(&traj1, e);
        vector<IShape *> queries;
        for (int i = 0; i < testtime; i++) {
            if (QueryType == 1) {
//          double pLow[2] = {random(39.992560,39.996714), random(116.319681,116.321562)};
                double pLow[2] = {random(0, 25000), random(0, 30000)};
                double pHigh[2] = {pLow[0] + random(1, 500), pLow[1] + random(1, 500)};
                Region r(pLow, pHigh, 2);
//          cout<<pLow[0]<<","<<pLow[1]<<endl;
                int t = int(random(0, PeriodLen));
                TimeRegion *tr = new TimeRegion(pLow, pHigh, t, t, 2);
                queries.emplace_back(tr);
            } else if (QueryType == 2) {
                queries.emplace_back(&traj1[int(random(0,10*traj1.size()))%traj1.size()].second);
            }
        }
        cout<<"here";

        string name0 = "name0", name1 = "name1", name2 = "name2", name3 = "name3", name4 = "name4",
                name5 = "name5", name6 = "name6", name7 = "name7";
        id_type indexIdentifier0, indexIdentifier1, indexIdentifier2, indexIdentifier3, indexIdentifier4,
                indexIdentifier5, indexIdentifier6, indexIdentifier7;
        IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name1, 4096),
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
        // applies a main memory random buffer on top of the persistent storage manager
        // (LRU buffer, etc can be created the same way).
//    ISpatialIndex* real = R2Tree::createAndBulkLoadNewR2Tree(
//            R2Tree::BulkLoadMethod::BLM_STR, ds2, *file0, 0.9, indexcap,100000, 2, indexIdentifier0);
//    ds1.rewind();
//        ISpatialIndex *r = RTree::createAndBulkLoadNewRTree(
//                RTree::BulkLoadMethod::BLM_STR, ds1, *file1, 0.9, indexcap, leafcap, 2, SpatialIndex::RTree::RV_RSTAR,
//                indexIdentifier1);
        ISpatialIndex *paar1 = PAARTree::createAndBulkLoadNewPAARTree(
                PAARTree::BulkLoadMethod::BLM_STR2, ds3, *file2, 0.9, indexcap, leafcap, 2, a, indexIdentifier2);
        ISpatialIndex *paar2 = PAARTree::createAndBulkLoadNewPAARTree(
                PAARTree::BulkLoadMethod::BLM_STR2, ds34, *file3, 0.9, indexcap, leafcap, 2, b, indexIdentifier3);
        ISpatialIndex *paar3 = PAARTree::createAndBulkLoadNewPAARTree(
                PAARTree::BulkLoadMethod::BLM_STR2, ds38, *file4, 0.9, indexcap, leafcap, 2, c, indexIdentifier4);
        ISpatialIndex *paar4 = PAARTree::createAndBulkLoadNewPAARTree(
                PAARTree::BulkLoadMethod::BLM_STR2, ds316, *file5, 0.9, indexcap, leafcap, 2, d, indexIdentifier5);
        ISpatialIndex *paar5 = PAARTree::createAndBulkLoadNewPAARTree(
                PAARTree::BulkLoadMethod::BLM_STR2, ds332, *file5, 0.9, indexcap, leafcap, 2, e, indexIdentifier6);

//    real->m_DataType=TrajectoryType;
//        r->m_DataType = TrajectoryType;
        paar1->m_DataType = TrajectoryType;
        paar2->m_DataType = TrajectoryType;
        paar3->m_DataType = TrajectoryType;
        paar4->m_DataType = TrajectoryType;
        paar5->m_DataType = TrajectoryType;
        cerr << "start query!" << endl << endl << endl;
//        TreeQueryBatch(r, queries);
        cerr << "\n\n\npaa-"<<a<<"\n";
        TreeQueryBatch(paar1, queries);
        cerr << "\n\n\npaa-"<<b<<"\n";
        TreeQueryBatch(paar2, queries);
        cerr << "\n\n\npaa-"<<c<<"\n";
        TreeQueryBatch(paar3, queries);
        cerr << "\n\n\npaa-"<<d<<"\n";
        TreeQueryBatch(paar4, queries);
        cerr << "\n\n\npaa-"<<e<<"\n";
        TreeQueryBatch(paar5, queries);
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