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
#define random(x,y) (((double)rand()/RAND_MAX)*(y-x+1)+x)
#include <spatialindex/SpatialIndex.h>
//#define sourceFile "/home/chuang/geolifedatasimplify.csv"
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
    MbrStream(vector<pair<int,Region>> other){mbrs=other;}
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
        throw Tools::NotSupportedException("Operation not supported.");
    }

    virtual void rewind(){i=0;}
};
class MbbcStream: public IDataStream{
public:
    vector<pair<int,Mbbc>> mbbcs;
    int i=0;
    MbbcStream(vector<pair<int,Mbbc>> other){mbbcs=other;}
    MbbcStream(vector<pair<id_type ,Trajectory> > period){
        mbbcs.clear();
        for(auto idt:period){
            Mbbc bc;
            idt.second.getMbbc(bc);
            mbbcs.emplace_back(make_pair(idt.first,bc));
        }
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
        throw Tools::NotSupportedException("Operation not supported.");
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
double naivetime(string l){
    int h;
    int m;
    int s;
    if(l.size()==9){
        h = stringToNum<int>(l.substr(0,2));
        m = stringToNum<int>(l.substr(3,5));
        s = stringToNum<int>(l.substr(6,8));
    } else{
        h = stringToNum<int>(l.substr(0,1));
        m = stringToNum<int>(l.substr(2,4));
        s = stringToNum<int>(l.substr(5,7));
    }
    return 10000*h+10000*m/60+100*s/60;

}
struct xyt{
        double x;
        double y;
        double t;
        };
void swap(double &x,double &y){
    double z;
    z=x;x=y;y=z;
}
xyt makemid(xyt p1, xyt p2, double t){
    if(p1.t>p2.t){
        swap(p1.x,p2.x);
        swap(p1.y,p2.y);
        swap(p1.t,p2.t);
    }
    double h1= (t-p1.t)/(p2.t-p1.t);
    double h2= (p2.t-t)/(p2.t-p1.t);
    double x=h1*p1.x+h2*p2.x;
    double y=h1*p1.y+h2*p2.y;
    xyt ret ={x,y,t};
    return ret;
}
vector< vector<xyt> > cuttraj(vector<xyt> traj){
    vector< vector<xyt> > segments(24);
    int oldpd=traj.at(0).t/10000;
    for(int i=0;i<traj.size();i++){
        int newpd=traj[i].t/10000;
        if(newpd!=oldpd){
            xyt mid1=makemid(traj.at(i-1),traj[i],(oldpd+1)*10000);
            segments.at(oldpd).emplace_back(mid1);
            if(int(traj[i].t)%10000!=0){
                xyt mid2=makemid(traj.at(i-1),traj[i],newpd*10000);
                segments.at(newpd).emplace_back(mid2);
            }
            oldpd=newpd;
        }
        segments.at(newpd).emplace_back(traj[i]);
    }
    return segments;
}

Mbbc toMbbc(vector<xyt> seg){
    if(seg.empty()) cout<<"no! a empty MBBC!"<<endl;
    double startx=seg.begin()->x,starty=seg.begin()->y,startt=seg.begin()->t;
    double endx=seg.back().x,endy=seg.back().y,endt=seg.back().t;
    double maxvxP=0,maxvxN=0,maxvyP=0,maxvyN=0;
    double minx=startx,maxx=startx,miny=starty,maxy=starty;
    for(int i=0;i<seg.size();i++){
        if(seg[i].t!=startt){
            double vx=(seg[i].x-startx)/(seg[i].t-startt);
            if(vx>maxvxP) maxvxP=vx;
            if(vx<maxvxN) maxvxN=vx;

            double vy=(seg[i].y-starty)/(seg[i].t-startt);
            if(vy>maxvyP) maxvyP=vy;
            if(vy<maxvyN) maxvyN=vy;
        }
        if(seg[i].t!=endt){
            double vx=(endx-seg[i].x)/(endt-seg[i].t);
            if(vx>maxvxP) maxvxP=vx;
            if(vx<maxvxN) maxvxN=vx;
            double vy=(endy-seg[i].y)/(endt-seg[i].t);
            if(vy>maxvyP) maxvyP=vy;
            if(vy<maxvyN) maxvyN=vy;
        }

        if(seg[i].x<minx) minx=seg[i].x;
        if(seg[i].x>maxx) maxx=seg[i].x;
        if(seg[i].y<miny) miny=seg[i].y;
        if(seg[i].y>maxy) maxy=seg[i].y;
    }
    double sLow[2]={startx,starty};
    double sHigh[2]={startx,starty};
    double eLow[2]={endx,endy};
    double eHigh[2]={endx,endy};
    double vLow[2]={maxvxN,maxvyN};
    double vHigh[2]={maxvxP,maxvyP};
    double wLow[2]={minx,miny};
    double wHigh[2]={maxx,maxy};
    double stime=int(startt/10000)*10000;
    return Mbbc(Region(sLow,sHigh,2),Region(eLow,eHigh,2),
            Region(vLow,vHigh,2),Region(wLow,wHigh,2),stime,stime+10000);
}
Region toMbr(vector<xyt> seg){
    if(seg.empty()) cout<<"no! a empty MBR!"<<endl;
    double startx=seg.begin()->x,starty=seg.begin()->y,startt=seg.begin()->t;
    double minx=startx,maxx=startx,miny=starty,maxy=starty;
    for(int i=0;i<seg.size();i++){
        if(seg[i].x<minx) minx=seg[i].x;
        if(seg[i].x>maxx) maxx=seg[i].x;
        if(seg[i].y<miny) miny=seg[i].y;
        if(seg[i].y>maxy) maxy=seg[i].y;
    }
    double pLow[2]={minx,miny};
    double pHigh[2]={maxx,maxy};
    return Region(pLow,pHigh,2);
}

void loadCsvToMbr(IDataStream &ds1,const vector<IShape*> &queries){
    set<int> ids;
    multimap<int,xyt> trajs;
    vector< vector< pair<int,Region> > > liar(24);
    string name = "name";
    id_type indexIdentifier;
    IStorageManager* diskfile = StorageManager::createNewDiskStorageManager(name, 4096);
    // Create a new storage manager with the provided base name and a 4K page size.

    StorageManager::IBuffer* file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
    // applies a main memory random buffer on top of the persistent storage manager
    // (LRU buffer, etc can be created the same way).
    clock_t start,end;
    start=clock();

    ISpatialIndex* tree = RTree::createAndBulkLoadNewRTree(
            RTree::BLM_STR, ds1, *file, 0.9, indexcap,leafcap, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
    end=clock();
    cout<<"Tree Building Time: "<< end-start<<endl;
    RangeVisitor vis;
    start=clock();
    for(int i=0;i<queries.size();i++){
        tree->intersectsWithQuery(*queries[i],vis);
//        tree->nearestNeighborQuery(5,queries[i],vis);
    }
    end=clock();
    cout<<"Querying time: "<< end-start<<endl;
    cout<<"vis"<<vis.m_indexIO<<","<<vis.m_leafIO<<endl;
    std::cerr << *tree;
    std::cerr << "Buffer hits: " << file->getHits() << std::endl;
    std::cerr << "Index ID: " << indexIdentifier << std::endl;

    bool ret = tree->isIndexValid();
    if (ret == false) std::cerr << "ERROR: Structure is invalid!" << std::endl;
    else std::cerr << "The stucture seems O.K." << std::endl;

    delete tree;
    delete file;
    delete diskfile;
}


void loadCsvToMbbc(R2Tree::BulkLoadMethod blm,IDataStream &ds2,const vector<IShape*> &queries){
    set<int> ids;
    multimap<int,xyt> trajs;
    vector< vector< pair<int,Mbbc> > > liar(24);
    string name = "name";
    id_type indexIdentifier;
    IStorageManager* diskfile = StorageManager::createNewDiskStorageManager(name, 4096);
    // Create a new storage manager with the provided base name and a 4K page size.

    StorageManager::IBuffer* file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
    // applies a main memory random buffer on top of the persistent storage manager
    // (LRU buffer, etc can be created the same way).

    clock_t start,end;
    start=clock();
    ISpatialIndex* tree = R2Tree::createAndBulkLoadNewR2Tree(
            blm, ds2, *file, 0.9, indexcap,leafcap,2, indexIdentifier);
    end=clock();
    cout<<"Tree Building time: "<< end-start<<endl;
    bool ret = tree->isIndexValid();
    if (ret == false) std::cerr << "ERROR: Structure is invalid!" << std::endl;
    else std::cerr << "The stucture seems O.K." << std::endl;

    RangeVisitor vis;
    start=clock();
    for(int i=0;i<queries.size();i++){
        tree->intersectsWithQuery(*queries[i],vis);
//        tree->nearestNeighborQuery(5,queries[i],vis);
    }
    end=clock();
    cout<<"Querying time: "<< end-start<<endl;
    cout<<"vis"<<vis.m_indexIO<<","<<vis.m_leafIO<<endl;
    std::cerr << *tree;
    std::cerr << "Buffer hits: " << file->getHits() << std::endl;
    std::cerr << "Index ID: " << indexIdentifier << std::endl;


    delete tree;
    delete file;
    delete diskfile;

}


list<vector<pair<id_type ,Trajectory> > > loadCsvToTrajs(){
    //first level: vector of time period
    //second level: vector of segments in the time period
    ifstream inFile(sourceFile, ios::in);
    string lineStr;

    getline(inFile, lineStr);
    set<id_type> ids;
    multimap<id_type,xyt> trajs;
    list<vector<pair<id_type ,Trajectory> > > res(24);
    auto iperiod=res.begin();
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
            vector< vector<xyt> > segs = cuttraj(traj);//size 24
            for(int j =0;j<24;j++){
                vector<xyt> seg;
                vector<TimePoint> tps;
                seg=segs[j];
                for(auto p:seg){
                    double xy[]={p.x,p.y};
                    double faket=int(p.t)%10000;
                    tps.emplace_back(TimePoint(xy,faket,faket,dimension));
//                    tps.emplace_back(TimePoint(xy,p.t,p.t,dimension));
                }
                if(!tps.empty()){
                    iperiod->emplace_back(make_pair(id,Trajectory(tps)));
//                    iperiod++;
                }
            }
        }
    }
    return res;
}


int main(){
    srand((int)time(NULL));
    list<vector<pair<id_type ,Trajectory> > > trajs=loadCsvToTrajs();
    vector<pair<id_type ,Trajectory> >   empty;
//    vector<pair<id_type ,Trajectory> > tjtjtj;
//    tjtjtj.insert(tjtjtj.end(),   trajs[0].begin(),   trajs[0].end());
//    tjtjtj.insert(tjtjtj.end(),   trajs[0].begin(),   trajs[0].end());
//    tjtjtj.insert(tjtjtj.end(),   trajs[0].begin(),   trajs[0].end());
//    tjtjtj.insert(tjtjtj.end(),   trajs[0].begin(),   trajs[0].end());
//    TrajMbrStream ds1(trajs[0]);
//    TrajMbbcStream ds2(trajs[0]);
    MbbcStream ds2(trajs.front());
    for(auto period:trajs){
        period.swap(empty);
    }
    vector<IShape*> queries;
    for (int i = 0; i < testtime; i++){
        double pLow[2] = {random(31,40.5), random(110,122)};
        double pHigh[2] = {pLow[0]+random(0,0.01), pLow[1]+random(0,0.01)};
        Region r(pLow, pHigh, 2);
        int t =int(random(0,10000));
        TimeRegion *tr=new TimeRegion(pLow, pHigh, t, t, 2);
        queries.emplace_back(tr);
//        queries.emplace_back(&trajs[0][i].second);
    }
//    loadCsvToMbr(ds1,queries);
//    cout<<"\n\n\n\n";
    loadCsvToMbbc(R2Tree::BulkLoadMethod::BLM_STR,ds2,queries);
    cout<<"\n\n\n\n";
    loadCsvToMbbc(R2Tree::BulkLoadMethod::BLM_STR2,ds2,queries);
    cout<<"\n\n\n\n";
    loadCsvToMbbc(R2Tree::BulkLoadMethod::BLM_STR3,ds2,queries);
    return 0;
}