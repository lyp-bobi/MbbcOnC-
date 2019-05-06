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
#define sourceFile "D://geolifedata.csv"
#define testtime 100000
#define dimension 2
#define indexcap 5
#define leafcap 5

using namespace std;
using namespace SpatialIndex;



class RangeVisitor : public IVisitor
{
public:
    size_t m_indexIO;
    size_t m_leafIO;

public:
    RangeVisitor() : m_indexIO(0), m_leafIO(0) {}

    void visitNode(const INode& n)
    {
        if (n.isLeaf()) m_leafIO++;
        else m_indexIO++;

    }

    void visitData(const IData& d)
    {
        IShape* pS;
        d.getShape(&pS);
        // do something.
        delete pS;
//        cout<<"data"<<endl;

        // data should be an array of characters representing a Region as a string.
        byte* pData = 0;
        uint32_t cLen = 0;
        d.getData(cLen, &pData);
        // do something.
        //string s = reinterpret_cast<char*>(pData);
        //cout << s << endl;
        delete[] pData;

//        cout << d.getIdentifier() << endl;
        // the ID of this data entry is an answer to the query. I will just print it to stdout.
    }

    void visitData(std::vector<const IData*>& v)
    {
        cout << v[0]->getIdentifier() << " " << v[1]->getIdentifier() << endl;
    }
};

class MbrStream: public IDataStream{
public:
    vector<pair<int,Region>> mbrs;
    int i=0;
    void feedTraj(const vector<pair<id_type ,Trajectory> > *period){
        mbrs.clear();
        for(auto idt:*period){
            Region br;
            idt.second.getMBR(br);
            mbrs.emplace_back(make_pair(idt.first,br));
        }
        rewind();
    }
    virtual bool hasNext() override
    {
        return i<mbrs.size();
    }
    virtual IData* getNext() override{
        RTree::Data* d=new RTree::Data(sizeof(double), reinterpret_cast<byte*>(mbrs[i].second.m_pLow), mbrs[i].second, mbrs[i].first);
        i++;
        return d;
    }
    virtual uint32_t size()
    {
        return mbrs.size();
    }

    virtual void rewind(){i=0;}
};

class MbbcStream: public IDataStream{
public:
    vector<pair<int,Mbbc>> mbbcs;
    int i=0;
    void feedTraj(const vector<pair<id_type ,Trajectory> > *period){
        mbbcs.clear();
        for(auto idt:*period){
            Mbbc bc;
            idt.second.getMbbc(bc);
            mbbcs.emplace_back(make_pair(idt.first,bc));
        }
        rewind();
    }
    virtual bool hasNext() override
    {
        return i<mbbcs.size();
    }
    virtual IData* getNext() override{
        R2Tree::Data* d=new R2Tree::Data(sizeof(double), reinterpret_cast<byte*>(mbbcs[i].second.m_smbr.m_pLow), mbbcs[i].second, mbbcs[i].first);
        i++;
        return d;
    }
    virtual uint32_t size()
    {
        return mbbcs.size();
    }

    virtual void rewind(){i=0;}
};

class TrajMbrStream: public IDataStream{
public:
    vector<pair<id_type,Region> > mbrs;
    vector<pair<id_type,Trajectory> > trajs;

    int i=0;
    TrajMbrStream(vector<pair<id_type,Region> > otherr,vector<pair<id_type,Trajectory> > othert){
        assert(mbrs.size()==trajs.size());
        mbrs=otherr;
        trajs=othert;
    }
    TrajMbrStream(vector<pair<id_type ,Trajectory> > period){
        mbrs.clear();
        for(auto idt:period){
            Region mbr;
            idt.second.getMBR(mbr);
            mbrs.emplace_back(make_pair(idt.first,mbr));
        }
        trajs=period;
    }
    virtual bool hasNext() override
    {
        return i<mbrs.size();
    }
    virtual IData* getNext() override{
        byte* data;
        uint32_t len;
        trajs[i].second.storeToByteArray(&data,len);
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

class TrajMbbcStream: public IDataStream{
public:
    vector<pair<id_type,Mbbc> > mbbcs;
    vector<pair<id_type,Trajectory> > trajs;

    int i=0;
    TrajMbbcStream(vector<pair<id_type,Mbbc> > otherbc,vector<pair<id_type,Trajectory> > othert){
        assert(mbbcs.size()==trajs.size());
        mbbcs=otherbc;
        trajs=othert;
    }
    TrajMbbcStream(vector<pair<id_type ,Trajectory> > period){
        mbbcs.clear();
        for(auto idt:period){
            Mbbc bc;
            idt.second.getMbbc(bc);
            mbbcs.emplace_back(make_pair(idt.first,bc));
        }
        trajs=period;
    }
    virtual bool hasNext() override
    {
        return i<mbbcs.size();
    }
    virtual IData* getNext() override{
        byte* data;
        uint32_t len;
        trajs[i].second.storeToByteArray(&data,len);
        R2Tree::Data* d=new R2Tree::Data(len, data, mbbcs[i].second, mbbcs[i].first);
        i++;
        return d;
    }
    virtual uint32_t size()
    {
        return mbbcs.size();
    }

    virtual void rewind(){i=0;}
};

template <class Type>
Type stringToNum(const string& str)
{
    istringstream iss(str);
    Type num;
    iss >> num;
    return num;
}
//time division related
double naivetime(string l){
    int h;
    int m;
    int s;
    if(l.size()==8){
        h = stringToNum<int>(l.substr(0,2));
        m = stringToNum<int>(l.substr(3,5));
        s = stringToNum<int>(l.substr(6,8));
    } else{
        h = stringToNum<int>(l.substr(0,1));
        m = stringToNum<int>(l.substr(2,4));
        s = stringToNum<int>(l.substr(5,7));
    }
    return 10000*h+160*m+1.6*s;
}
int getPeriod(double time){
    int pd= int(floor(time))/PeriodLen;
    return pd;
}
int getMaxPeriod(){
    return 240000/PeriodLen;
}
double getPeriodStart(double time){
    int pd=getPeriod(time);
    return pd*PeriodLen;
}
double getPeriodEnd(double time){
    int pd=getPeriod(time);
    return (pd+1)*PeriodLen-0.00001;
}

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
    double x=h1*p1.x+h2*p2.x;
    double y=h1*p1.y+h2*p2.y;
    xyt ret ={x,y,t};
    return ret;
}
vector< vector<xyt> > cuttraj(vector<xyt> traj){
    vector< vector<xyt> > segments(getMaxPeriod());
    int oldpd=getPeriod(traj.at(0).t);
    for(int i=0;i<traj.size();i++){
        if(traj[i].t<traj[i-1].t) break;//stop reading when coming to a new day
        int newpd=getPeriod(traj[i].t);
        if(newpd-1==oldpd){
            xyt mid1=makemid(traj[i-1],traj[i],getPeriodEnd(traj[i-1].t));
            segments.at(oldpd).emplace_back(mid1);
            if(traj[i].t-getPeriodStart(traj[i].t!=0)){
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



list<vector<pair<id_type ,Trajectory> > > loadCsvToTrajs(){
    //first level: vector of time period
    //second level: vector of segments in the time period
    ifstream inFile(sourceFile, ios::in);
    string lineStr;

    getline(inFile, lineStr);
    set<id_type> ids;
    multimap<id_type,xyt> trajs;
    list<vector<pair<id_type ,Trajectory> > > res(getMaxPeriod());
    while (getline(inFile, lineStr)){
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
//        if(id==14718) cout<<x<<" "<<y<<" "<<t<<" "<<endl;
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
                    tps.emplace_back(TimePoint(xy,faket,faket,dimension));
//                    tps.emplace_back(TimePoint(xy,p.t,p.t,dimension));
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
void TreeQuery(ISpatialIndex* tree,const vector<IShape*> &queries){
    clock_t start,end;

    RangeVisitor vis;
    start=clock();
    for(int i=0;i<queries.size();i++){
        tree->intersectsWithQuery(*queries[i],vis);
//        tree->nearestNeighborQuery(5,queries[i],vis);
    }
    end=clock();
    cerr<<"Querying time: "<< end-start<<endl;
    cerr<<"VISIT NODE "<<vis.m_indexIO<<","<<vis.m_leafIO<<endl;
    cerr << *tree;
}

int main(){
    srand((int)time(NULL));
    list<vector<pair<id_type ,Trajectory> > > trajs=loadCsvToTrajs();
    vector<pair<id_type ,Trajectory> >   empty;
    MbrStream ds1;
    MbbcStream ds2;
    ds1.feedTraj(&trajs.front());
    ds2.feedTraj(&trajs.front());
//    for(auto period:trajs){
//        period.swap(empty);
//    }
    vector<IShape*> queries;
    for (int i = 0; i < testtime; i++){
        double pLow[2] = {random(31,40.5), random(110,122)};
        double pHigh[2] = {pLow[0]+random(0.01,0.1), pLow[1]+random(0.01,0.1)};
        Region r(pLow, pHigh, 2);
//        cout<<pLow[0]<<","<<pLow[1]<<endl;
        int t =int(random(0,PeriodLen));
        TimeRegion *tr=new TimeRegion(pLow, pHigh, t, t, 2);
        queries.emplace_back(tr);
//        queries.emplace_back(&trajs[0][i].second);
    }


    string name1 = "name1",name2 = "name2",name3 = "name3",name4 = "name4";
    id_type indexIdentifier1,indexIdentifier2,indexIdentifier3,indexIdentifier4;
    IStorageManager *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
        *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096),
        *diskfile3 = StorageManager::createNewDiskStorageManager(name3, 4096),
        *diskfile4 = StorageManager::createNewDiskStorageManager(name4, 4096);
    // Create a new storage manager with the provided base name and a 4K page size.
    StorageManager::IBuffer *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false),
        *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false),
        *file3 = StorageManager::createNewRandomEvictionsBuffer(*diskfile3, 10, false),
        *file4 = StorageManager::createNewRandomEvictionsBuffer(*diskfile4, 10, false);
    // applies a main memory random buffer on top of the persistent storage manager
    // (LRU buffer, etc can be created the same way).
    ISpatialIndex* r = RTree::createAndBulkLoadNewRTree(
            RTree::BulkLoadMethod::BLM_STR, ds1, *file1, 0.9, indexcap,leafcap, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier1);
    ISpatialIndex* r21 = R2Tree::createAndBulkLoadNewR2Tree(
            R2Tree::BulkLoadMethod::BLM_STR, ds2, *file2, 0.9, indexcap,leafcap,2, indexIdentifier2);
    ISpatialIndex* r22 = R2Tree::createAndBulkLoadNewR2Tree(
            R2Tree::BulkLoadMethod::BLM_STR2, ds2, *file3, 0.9, indexcap,leafcap,2, indexIdentifier3);
    ISpatialIndex* r23 = R2Tree::createAndBulkLoadNewR2Tree(
            R2Tree::BulkLoadMethod::BLM_STR3, ds2, *file4, 0.9, indexcap,leafcap,2, indexIdentifier4);

    cerr<<"start query!"<<endl<<endl<<endl;
    TreeQuery(r,queries);
    cout<<"\n\n\n\n";
    TreeQuery(r21,queries);
    cout<<"\n\n\n\n";
    TreeQuery(r22,queries);
    cout<<"\n\n\n\n";
    TreeQuery(r23,queries);
    return 0;
}