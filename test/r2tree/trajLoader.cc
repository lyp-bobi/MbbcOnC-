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
#include <spatialindex/SpatialIndex.h>

using namespace std;
using namespace SpatialIndex;

template <class Type>
Type stringToNum(const string& str)
{
    istringstream iss(str);
    Type num;
    iss >> num;
    return num;
}
double naivetime(string l){
    int h = stringToNum<int>(l.substr(0,2));
    int m = stringToNum<int>(l.substr(3,5));
    int s = stringToNum<int>(l.substr(6,8));
    return 10000*h+100*m+s;
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
        int newpd=traj.at(i).t/10000;
        //cout<<newpd<<endl;
        if(newpd!=oldpd){
            xyt mid1=makemid(traj.at(i-1),traj.at(i),(oldpd+1)*10000);
            xyt mid2=makemid(traj.at(i-1),traj.at(i),newpd*10000);
            segments.at(oldpd).push_back(mid1);
            segments.at(newpd).push_back(mid2);
        }
        segments.at(newpd).push_back(traj.at(i));
    }
    return segments;
}

Mbbc toMbbc(vector<xyt> seg){
    if(seg.empty()) cout<<"no! a empty MBBC!"<<endl;
    double startx=seg.begin()->x,starty=seg.begin()->y,startt=seg.begin()->t;
    double endx=seg.end()->x,endy=seg.end()->y,endt=seg.end()->t;
    double maxvxP=0,maxvxN=0,maxvyP=0,maxvyN=0;
    double minx=startx,maxx=startx,miny=starty,maxy=starty;
    for(int i=0;i<seg.size();i++){
        if(seg.at(i).t!=startt){
            double vx=(seg.at(i).x-startx)/(seg.at(i).t-startt);
            if(vx>maxvxP) maxvxP=vx;
            if(vx<maxvxN) maxvxN=vx;
            double vy=(seg.at(i).y-starty)/(seg.at(i).t-startt);
            if(vy>maxvyP) maxvyP=vy;
            if(vy<maxvyN) maxvyN=vy;
        }
        if(seg.at(i).t!=endt){
            double vx=(endx-seg.at(i).x)/(endt-seg.at(i).t);
            if(vx>maxvxP) maxvxP=vx;
            if(vx<maxvxN) maxvxN=vx;
            double vy=(endy-seg.at(i).y)/(endt-seg.at(i).t);
            if(vy>maxvyP) maxvyP=vy;
            if(vy<maxvyN) maxvyN=vy;
        }
        if(seg.at(i).x<minx) minx=seg.at(i).x;
        if(seg.at(i).x>maxx) maxx=seg.at(i).x;
        if(seg.at(i).y<miny) miny=seg.at(i).y;
        if(seg.at(i).y>maxy) maxy=seg.at(i).y;
    }
    double sLow[2]={startx,starty};
    double sHigh[2]={startx,starty};
    double eLow[2]={endx,endy};
    double eHigh[2]={endx,endy};
    double vLow[2]={maxvxN,maxvyN};
    double vHigh[2]={maxvxP,maxvyP};
    double pLow[2]={minx,miny};
    double pHigh[2]={maxx,maxy};
    double stime=int(startt/10000)*10000;
    return *new Mbbc(*new Region(sLow,sHigh,2),*new Region(eLow,eHigh,2),
            *new Region(vLow,vHigh,2),*new Region(pLow,pHigh,2),stime,stime+10000);
}

void loadCsvToMbbc(){
    ifstream inFile("/home/chuang/geolifedata.csv", ios::in);
    string lineStr;
    //cout<<"hi"<<endl;
    getline(inFile, lineStr);
    vector<int> ids;
    multimap<int,xyt> trajs;
    vector< vector<Mbbc> > liar(24);
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
        ids.push_back(id);
        trajs.insert(make_pair(id,p));
    }
    //cout<<ids.size();
    for(int i=0;i<ids.size();i++){
        int id= ids.at(i);
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
                if(!seg.at(j).empty()){
                    liar.at(j).push_back(toMbbc(seg.at(j)));
                }
            }
        }

    }
    for (int i = 0; i < 24; ++i) {
        cout<<liar.at(i).size()<<endl;
    }


}

int main(){
    loadCsvToMbbc();
    return 0;
}