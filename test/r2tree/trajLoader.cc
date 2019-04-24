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
#define sourceFile "/home/chuang/geolifedata.csv"
#define testtime 10000
#define dimension 2

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

class MyDataStream1: public IDataStream{
public:
    vector<pair<int,Region>> mbrs;
    int i=0;
    MyDataStream1(vector<pair<int,Region>> other){mbrs=other;}
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
class MyDataStream2: public IDataStream{
public:
    vector<pair<int,Mbbc>> mbbcs;
    int i=0;
    MyDataStream2(vector<pair<int,Mbbc>> other){mbbcs=other;}
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

template <class Type>
Type stringToNum(const string& str)
{
    istringstream iss(str);
    Type num;
    iss >> num;
    return num;
}
double naivetime(string l){
    if(l.size()==9){
        int h = stringToNum<int>(l.substr(0,2));
        int m = stringToNum<int>(l.substr(3,5));
        int s = stringToNum<int>(l.substr(6,8));
        return 10000*h+100*m+s;
    } else{
        int h = stringToNum<int>(l.substr(0,1));
        int m = stringToNum<int>(l.substr(2,4));
        int s = stringToNum<int>(l.substr(5,7));
        return 10000*h+100*m+s;
    }

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
            segments.at(oldpd).push_back(mid1);
            if(int(traj[i].t)%10000!=0){
                xyt mid2=makemid(traj.at(i-1),traj[i],newpd*10000);
                segments.at(newpd).push_back(mid2);
            }
            oldpd=newpd;
        }
        segments.at(newpd).push_back(traj[i]);
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

void loadCsvToMbr(const vector<IShape*> &queries){
    ifstream inFile(sourceFile, ios::in);
    string lineStr;
    //cout<<"hi"<<endl;
    getline(inFile, lineStr);
    set<int> ids;
    multimap<int,xyt> trajs;
    vector< vector< pair<int,Region> > > liar(24);
    while (getline(inFile, lineStr)){
        //cout<<lineStr<<endl;
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
    //cout<<ids.size();
    for(auto id:ids){
        multimap<int,xyt>::iterator beg,end,iter;
        vector<xyt> traj;
        beg = trajs.lower_bound(id);
        end = trajs.upper_bound(id);
        for(iter=beg;iter!=end;iter++){
            traj.push_back(iter->second);
        }
        trajs.erase(id);
        if(!traj.empty()){
            vector< vector<xyt> > seg = cuttraj(traj);//size 24

            for(int j =0;j<24;j++){
                if(!seg[j].empty()){
                    liar.at(0).push_back(make_pair(id,toMbr(seg[j])));
                }
            }
        }
    }
    MyDataStream1 ds1(liar.at(0));
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
            RTree::BLM_STR, ds1, *file, 0.9, 4,4, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
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


void loadCsvToMbbc(R2Tree::BulkLoadMethod blm,const vector<IShape*> &queries){
    ifstream inFile(sourceFile, ios::in);
    string lineStr;

    getline(inFile, lineStr);
    set<int> ids;
    multimap<int,xyt> trajs;
    vector< vector< pair<int,Mbbc> > > liar(24);
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
        multimap<int,xyt>::iterator beg,end,iter;
        vector<xyt> traj;
        beg = trajs.lower_bound(id);
        end = trajs.upper_bound(id);
        for(iter=beg;iter!=end;iter++){
            traj.push_back(iter->second);
        }
        trajs.erase(id);
        if(!traj.empty()){
            vector< vector<xyt> > seg = cuttraj(traj);//size 24

            for(int j =0;j<24;j++){
                if(!seg[j].empty()){
                    Mbbc bc=toMbbc(seg[j]);
                    bc.m_startTime=0;
                    bc.m_endTime=10000;
                    liar.at(0).push_back(make_pair(id,bc));
//                    cout<<id<<"\n"<<bc.toString()<<"\n";
                }
            }
        }

    }
//    for (int i = 0; i < 24; ++i) {
//        cout<<liar[i].size()<<endl;
//    }
    MyDataStream2 ds2(liar.at(0));
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
            blm, ds2, *file, 0.9, 4,4,2, indexIdentifier);
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


vector<vector<pair<id_type ,Trajectory> > > loadCsvToTrajs(){
    //first level: vector of time period
    //second level: vector of segments in the time period
    ifstream inFile(sourceFile, ios::in);
    string lineStr;

    getline(inFile, lineStr);
    set<id_type> ids;
    multimap<id_type,xyt> trajs;
    vector<vector<pair<id_type ,Trajectory> > > res(24);
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
            traj.push_back(iter->second);
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
                    tps.push_back(TimePoint(xy,p.t,p.t,dimension));
                }
                if(!tps.empty()){
                    res[j].push_back(make_pair(id,Trajectory(tps)));
                }
            }
        }
    }
    return res;
}

int main(){
    srand((int)time(NULL));
    auto trajs=loadCsvToTrajs();
    vector<IShape*> queries;
    for (int i = 0; i < testtime; i++){
        double pLow[2] = {random(31,40.5), random(110,122)};
        double pHigh[2] = {pLow[0]+random(0,0.01), pLow[1]+random(0,0.01)};
        Region r(pLow, pHigh, 2);
        int t =int(random(0,10000));
        TimeRegion *tr=new TimeRegion(pLow, pHigh, t, t, 2);
        queries.push_back(tr);
//        queries.push_back(&trajs[0][i].second);
    }

//    loadCsvToMbr(queries);
//    cout<<"\n\n\n\n";
//    loadCsvToMbbc(R2Tree::BulkLoadMethod::BLM_STR,queries);
//    cout<<"\n\n\n\n";
    loadCsvToMbbc(R2Tree::BulkLoadMethod::BLM_STR2,queries);
    cout<<"\n\n\n\n";
    loadCsvToMbbc(R2Tree::BulkLoadMethod::BLM_STR3,queries);
    return 0;
}