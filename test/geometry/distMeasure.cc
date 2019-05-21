//
// Created by Chuang on 2019/5/20.
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
#include "spatialindex/SpatialIndex.h"


#define linesToRead 1e7
#define dimension 2


using namespace std;
using namespace SpatialIndex;

//Given an array of points, "findMaximumDistance" calculates the GPS point which have largest distance from the line formed by first and last points in RDP algorithm. Returns the index of the point in the array and the distance.

const std::pair<int, double> findMaximumDistance(const vector<SpatialIndex::TimePoint>& points) {
    SpatialIndex::TimePoint firstpoint=points[0];
    SpatialIndex::TimePoint lastpoint=points[points.size()-1];
    int index=0;  //index to be returned
    double Mdist=-1; //the Maximum distance to be returned

    //distance calculation
    for(int i=1;i<points.size()-1;i++){ //traverse through second point to second last point
        double Dist=SpatialIndex::TimePoint::makemid(firstpoint,lastpoint,points[i].m_startTime).getMinimumDistance(points[i]);
        if (Dist>Mdist){
            Mdist=Dist;
            index=i;
        }
    }
    return std::make_pair(index, Mdist);
}



vector<SpatialIndex::TimePoint> simplifyWithRDP(vector<SpatialIndex::TimePoint>& Points, double epsilon) {
    if(Points.size()<3){  //base case 1
        return Points;
    }
    std::pair<int, double> maxDistance=findMaximumDistance(Points);
    if(maxDistance.second>=epsilon){
        int index=maxDistance.first;
        vector<SpatialIndex::TimePoint>::iterator it=Points.begin();
        vector<SpatialIndex::TimePoint> path1(Points.begin(),it+index+1); //new path l1 from 0 to index
        vector<SpatialIndex::TimePoint> path2(it+index,Points.end()); // new path l2 from index to last

        vector<SpatialIndex::TimePoint> r1 =simplifyWithRDP(path1,epsilon);
        vector<SpatialIndex::TimePoint> r2=simplifyWithRDP(path2,epsilon);

        //Concat simplified path1 and path2 together
        vector<SpatialIndex::TimePoint> rs(r1);
        rs.pop_back();
        rs.insert(rs.end(),r2.begin(),r2.end());
        return rs;
    }
    else { //base case 2, all points between are to be removed.
        vector<SpatialIndex::TimePoint> r(1,Points[0]);
        r.push_back(Points[Points.size()-1]);
        return r;
    }
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
    double x=h2*p1.x+h1*p2.x;
    double y=h2*p1.y+h1*p2.y;
    xyt ret ={x,y,t};
    return ret;
}

vector< vector<xyt> > cuttraj(vector<xyt> traj){
    vector< vector<xyt> > segments(getMaxPeriod());
    int oldpd=getPeriod(traj.at(0).t);
    for(int i=0;i<traj.size();i++){
        if(i>0&&traj[i].t<traj[i-1].t) break;//stop reading when coming to a new day
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

list<vector<pair<id_type ,Trajectory> > > loadGTToTrajs(){
    //first level: vector of time period
    //second level: vector of segments in the time period
    cerr<<"loading generated trajectories from txt to trajectories"<<endl;
    ifstream inFile("D://the99trajs.txt", ios::in);
    string lineStr;
    set<id_type> ids;
    multimap<id_type,xyt> trajs;
    list<vector<pair<id_type ,Trajectory> > > res(getMaxPeriod());
    int curLine=0;
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
//        cout<<id<<x<<y<<t<<endl;
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

int main(){
    auto trajs=loadGTToTrajs();
    vector<pair<id_type ,Trajectory> > traj1=trajs.front();
    auto q=traj1[0].second;
    double rate1=0,rate2=0,rate3=0;
    for(int i=1;i<traj1.size();i++){
        auto t=traj1[i];
        double real=q.getMinimumDistance(t.second);
//        cout<<"origin "<<real<<endl;
        MBRk brk;
        t.second.getMBRk(4,brk);
        rate1+=q.getMinimumDistance(brk)/real;
//        cout<<"mbr-4 "<<q.getMinimumDistance(brk)/real<<endl;
        t.second.getMBRk(8,brk);
        rate2+=q.getMinimumDistance(brk)/real;
//        cout<<"mbr-8 "<<q.getMinimumDistance(brk)/real<<endl;
        MBBCk bbck;
        t.second.getMBBCk(4,bbck,5000);
//        cout<<bbck;
        rate3+=q.getMinimumDistance(bbck)/real;
//        cout<<"mbbc-4 "<<q.getMinimumDistance(bbck)/real<<endl;

    }
    rate1/=traj1.size();
    rate2/=traj1.size();
    rate3/=traj1.size();
    cout<<rate1<<rate2<<rate3;
}