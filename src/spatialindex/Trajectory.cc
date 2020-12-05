//
// Created by chuang on 4/23/19.
//

#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>
#include <fstream>


#include <spatialindex/SpatialIndex.h>
#include <storagemanager/TrajStore.h>

double calcuTime[10]={0,0,0,0,0,0,0,0};
int testPhase=0;

int disttype = 0;
bool bUsingVerboseLB = false;
//0 for IED, 1 for MaxSED

using namespace SpatialIndex;
using std::vector;
using std::cout;
using std::endl;
using std::sqrt;

Trajectory::Trajectory() {

}
Trajectory::Trajectory(std::vector<SpatialIndex::STPoint>& in)
{
    m_points=in;
}

Trajectory::Trajectory(bool fakehead, bool fakeback,std::vector<SpatialIndex::STPoint> &in) {
#ifndef NDEBUG
    if(in.front().m_time>=in.back().m_time){
        std::cerr<<"wrong trajectory!\n"<<*this;
    }
#endif
    m_points=in;
    m_fakehead=fakehead;
    m_fakeback=fakeback;
}

Trajectory::Trajectory(const SpatialIndex::Trajectory &in) {
    m_points=in.m_points;
    m_fakehead=in.m_fakehead;
    m_fakeback=in.m_fakeback;
}

Trajectory& Trajectory::operator=(const Trajectory& r)
{
    if(this != &r)
    {
        m_points=r.m_points;
    }

    return *this;
}

bool Trajectory::operator==(const SpatialIndex::Trajectory &r) const {
    if (m_points==r.m_points)
        return true;
    return false;
}
//
// IObject interface
//
Trajectory* Trajectory::clone() {
    return new Trajectory(*this);
}
//
// ISerializable interface
//
uint32_t Trajectory::getByteArraySize() const {
    return 2* sizeof(bool)+sizeof(unsigned long)+(m_dimension+1)* sizeof(double)*m_points.size();
}

void Trajectory::loadFromByteArray(const uint8_t* ptr) {
    memcpy(&m_fakehead, ptr, sizeof(bool));
    ptr += sizeof(bool);
    memcpy(&m_fakeback, ptr, sizeof(bool));
    ptr += sizeof(bool);
    unsigned long size;
    memcpy(&size, ptr, sizeof(unsigned long));
    ptr += sizeof(unsigned long);
    m_points.clear();
    STPoint p;
    p.makeDimension(m_dimension);
    for(int i=0;i<size;i++){
        memcpy(&p.m_time, ptr, sizeof(double));
        ptr += sizeof(double);
        memcpy(p.m_pCoords, ptr, m_dimension * sizeof(double));
        m_points.emplace_back(p);
        if(i!=size-1){
            ptr += m_dimension * sizeof(double);
        }
    }
}

void Trajectory::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    memcpy(ptr, &m_fakehead, sizeof(bool));
    ptr += sizeof(bool);
    memcpy(ptr, &m_fakeback, sizeof(bool));
    ptr += sizeof(bool);
    unsigned long size=m_points.size();
    memcpy(ptr, &size, sizeof(unsigned long));
    ptr += sizeof(unsigned long);
    for(int i=0;i<size;i++){
        memcpy(ptr, &m_points[i].m_time, sizeof(double));
        ptr += sizeof(double);
        memcpy(ptr, m_points[i].m_pCoords, m_dimension * sizeof(double));
        if(i!=size-1){
            ptr += m_dimension * sizeof(double);
        }
    }
}


//
// IShape interface
//
bool Trajectory::intersectsShape(const SpatialIndex::IShape& s) const {
    const Trajectory* pTrajectory = dynamic_cast<const Trajectory*>(&s);
    if (pTrajectory != 0) return intersectsTrajectory(*pTrajectory);

    const TimeRegion* ptr = dynamic_cast<const TimeRegion*>(&s);
    if (ptr != 0) return intersectsTimeRegion(*ptr);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return intersectsRegion(*pr);

    const Cylinder* pcy = dynamic_cast<const Cylinder*>(&s);
    if (pcy != 0) return intersectsCylinder(*pcy);

    throw Tools::IllegalStateException(
            "Trajectory::intersectsShape: Not implemented yet!"
    );
}

STPoint Trajectory::getPointAtTime(const double time) const {
    if(time<m_points.front().m_time||time>m_points.back().m_time){
        throw Tools::IllegalArgumentException(
                "Trajectory::getPointAtTime: time"+std::to_string(time)+"is illegal."
                );}
    if(m_points.size()==1){
        return m_points[0];
    }
    auto pre =m_points.begin(),next=m_points.begin();
    while(next->m_time-pre->m_time<1e-7&&next!=m_points.end()) next++;
    while(next->m_time<time&&next!=m_points.end()){
        pre=next;
        while(next->m_time-pre->m_time<1e-7&&next!=m_points.end()) next++;
    }
    double h1= (time-pre->m_time)/(next->m_time-pre->m_time);
    double h2= (next->m_time-time)/(next->m_time-pre->m_time);
    STPoint ret;
    ret.makeDimension(m_dimension);
    for (int i = 0; i < m_dimension; ++i) {
        ret.m_pCoords[i] = h2 * pre->m_pCoords[i] + h1 * next->m_pCoords[i];
    }
    ret.m_time=time;
    return ret;
}

bool Trajectory::intersectsTimeRegion(const SpatialIndex::TimeRegion &in) const {
    if(in.m_startTime==in.m_endTime){//time slice
        if(in.m_startTime<m_points.front().m_time||in.m_startTime>m_points.back().m_time){
            return false;
        }
        STPoint tp=getPointAtTime(in.m_startTime);
        return tp.intersectsShape(in);
    }else{
        throw Tools::NotSupportedException("time interval range not supported");
    }
}
bool Trajectory::intersectsRegion(const Region& in) const{
    if(m_dimension==in.m_dimension-1){
        double ts=std::max(m_startTime(),in.m_pLow[in.m_dimension-1]),
                te=std::min(m_endTime(),in.m_pHigh[in.m_dimension-1]);
        if(ts>te) return false;
        if(ts==te){
            Region mbr2d(in.m_pLow,in.m_pHigh,in.m_dimension-1);
            return mbr2d.containsPoint(getPointAtTime(ts));
        }
        fakeTpVector timedTraj(&m_points,ts,te);
        for(int i=0;i<timedTraj.m_size-1;i++){
            if(line2MBRMinSED(timedTraj[i],timedTraj[i+1],in)<1e-7){
                return true;
            }
        }
        return false;
//        Region spatial(in.m_pLow,in.m_pHigh,m_dimension);
//        if(in.m_pHigh[m_dimension]<m_points.front().m_time||in.m_pLow[m_dimension]>m_points.back().m_time){
//            return false;
//        }
//        STPoint tp=getPointAtTime(in.m_pLow[m_dimension]);
//        return tp.intersectsShape(spatial);
    }
    else if(m_dimension==in.m_dimension) {
        std::cerr<<"this part should not be used. must something wrong\n";
        for (int i = 0; i < m_points.size(); i++) {
            if (m_points[i].intersectsShape(in)) {
                return true;
            }
        }
    }
    return false;
}
bool Trajectory::intersectsCylinder(const SpatialIndex::Cylinder &in) const {
    double ts=std::max(m_startTime(),in.m_startTime),
            te=std::min(m_endTime(),in.m_endTime);
    if(ts>te) return false;
    if(ts==te){
        return Point(in.m_p,2).getMinimumDistance(getPointAtTime(ts))<in.m_r;
    }
    fakeTpVector timedTraj(&m_points,ts,te);
    for(int i=0;i<timedTraj.m_size-1;i++){
        if(line2lineMinSED(timedTraj[i],timedTraj[i+1],STPoint(in.m_p,timedTraj[i].m_time,2),STPoint(in.m_p,timedTraj[i+1].m_time,2))<in.m_r){
            return true;
        }
    }
    return false;
}
bool Trajectory::intersectsTrajectory(const Trajectory& in) const{
    throw Tools::NotSupportedException(
            "Trajectory:::getMinimumDistance Not implemented yet!"
    );
}

bool Trajectory::containsShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException(
            "Trajectory::...: Not implemented yet!"
    );
}
bool Trajectory::touchesShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException(
            "Trajectory::...: Not implemented yet!"
    );
}
void Trajectory::getCenter(Point& out) const{
    throw Tools::NotSupportedException("not supported now");
}
uint32_t Trajectory::getDimension() const{return 2;}
void Trajectory::getMBR(Region& out) const{
    out.makeInfinite(m_dimension+1);
    Region mid;
    mid.makeInfinite(m_dimension);
    for(int i=0;i<m_points.size();i++){
        mid.combinePoint(m_points[i]);
    }
    out.m_pLow[0] = mid.m_pLow[0];
    out.m_pLow[1] = mid.m_pLow[1];
    out.m_pLow[2] = m_points.front().m_time;
    out.m_pHigh[0] = mid.m_pHigh[0];
    out.m_pHigh[1] = mid.m_pHigh[1];
    out.m_pHigh[2] = m_points.back().m_time;
}


void Trajectory::getMBC(SpatialIndex::MBC &out) const {
    if(m_points.size()<=1){
        std::cerr<<"WARNING: getting MBC at a Trajectory with 0 or 1 points\n";
        out.makeInfinite(m_dimension+1);
        return;
    }
    if(std::fabs(m_endTime()-m_startTime())<1e-7){
        std::cerr<<"probably wrong segment at"<<*this<<"\n";
    }
    double startx=m_points[0].m_pCoords[0],starty=m_points[0].m_pCoords[1],startt=m_points[0].m_time;
    double endx=m_points.back().m_pCoords[0],endy=m_points.back().m_pCoords[1],endt=m_points.back().m_time;
    double avx=(endx-startx)/(endt-startt),avy=(endy-starty)/(endt-startt);

    STPoint p1=*m_points.begin(),p2=m_points.back();
    double rd=0,rv=0;
    double vx,vy;
    for(int i=0;i<m_points.size();i++){
        STPoint p;
        double ptime=m_points[i].m_time;
        double prd,prv;
        if(ptime-startt>0){
            vx=(m_points[i].m_pCoords[0]-startx)/(ptime-startt);
            vy=(m_points[i].m_pCoords[1]-starty)/(ptime-startt);
            prv=std::sqrt((vx-avx)*(vx-avx)+(vy-avy)*(vy-avy));
            if(prv>rv) rv=prv;
        }
        if(endt-m_points[i].m_time>0) {
            vx = (endx - m_points[i].m_pCoords[0]) / (endt - ptime);
            vy = (endy - m_points[i].m_pCoords[1]) / (endt - ptime);
            prv=std::sqrt((vx-avx)*(vx-avx)+(vy-avy)*(vy-avy));
            if(prv>rv) rv=prv;
        }
        prd=m_points[i].getMinimumDistance(STPoint::makemid(p1,p2,ptime));
        if(prd>rd) rd=prd;
    }
    out=MBC(p1.m_pCoords,p2.m_pCoords,startt,endt,m_dimension,rd,rv);
}

void Trajectory::getMBRfull(SpatialIndex::Region &out) const {
    out.makeInfinite(m_dimension+1);
    double *pc=new double[m_dimension+1];
    for(int i=0;i<m_points.size();i++){
        for(int d=0;d<m_dimension;d++) pc[d]=m_points[i].m_pCoords[d];
        pc[m_dimension]=m_points[i].m_time;
        out.combinePoint(Point(pc,m_dimension+1));
//        std::cout<<out<<std::endl;
    }
    delete[](pc);
}

void Trajectory::getTimeMBR(SpatialIndex::TimeRegion &out) const {
    out.makeInfinite(m_dimension);
    for(int i=0;i<m_points.size();i++){
        out.combinePoint(m_points[i]);
        if(m_points[i].m_time<out.m_startTime) out.m_startTime=m_points[i].m_time;
        if(m_points[i].m_time>out.m_endTime) out.m_endTime=m_points[i].m_time;
    }
}

const std::pair<int, double> findMaximumDistance(const vector<SpatialIndex::STPoint>& points) {
    SpatialIndex::STPoint firstpoint=points[0];
    SpatialIndex::STPoint lastpoint=points[points.size()-1];
    int index=0;  //index to be returned
    double Mdist=-1; //the Maximum distance to be returned

    //distance calculation
    for(int i=1;i<points.size()-1;i++){ //traverse through second point to second last point
        double Dist=SpatialIndex::STPoint::makemid(firstpoint,lastpoint,points[i].m_time).getMinimumDistance(points[i]);
        if (Dist>Mdist){
            Mdist=Dist;
            index=i;
        }
    }
    return std::make_pair(index, Mdist);
}
std::vector<SpatialIndex::STPoint> Trajectory::simplifyWithRDP(const std::vector<SpatialIndex::STPoint> &Points,
                                                                 double epsilon){
    if(Points.size()<3){  //base case 1
        return Points;
    }
    std::pair<int, double> maxDistance=findMaximumDistance(Points);
    if(maxDistance.second>=epsilon){
        int index=maxDistance.first;
        vector<SpatialIndex::STPoint>::const_iterator it=Points.begin();
        vector<SpatialIndex::STPoint> path1(Points.begin(),it+index+1); //new path l1 from 0 to index
        vector<SpatialIndex::STPoint> path2(it+index,Points.end()); // new path l2 from index to last

        vector<SpatialIndex::STPoint> r1 =simplifyWithRDP(path1,epsilon);
        vector<SpatialIndex::STPoint> r2=simplifyWithRDP(path2,epsilon);

        //Concat simplified path1 and path2 together
        vector<SpatialIndex::STPoint> rs(r1);
        rs.pop_back();
        rs.insert(rs.end(),r2.begin(),r2.end());
        return rs;
    }
    else { //base case 2, all points between are to be removed.
        vector<SpatialIndex::STPoint> r(1,Points[0]);
        r.emplace_back(Points[Points.size()-1]);
        return r;
    }
}

std::vector<std::vector<SpatialIndex::STPoint>> Trajectory::simplifyWithRDPN(const std::vector<SpatialIndex::STPoint> &Points,
                                                               int numPart){
    if(Points.size()<numPart+1){  //base case 1
        throw Tools::IllegalStateException("simplifyWithRDPN:Error");
    }
    std::pair<int, double> maxDistance=findMaximumDistance(Points);
    vector<vector<SpatialIndex::STPoint>> paths;
    paths.emplace_back(Points);
    while(paths.size()<numPart){
        int pathIndex=0;
        std::pair<int, double> md={0,0};
        for(int i=0;i<paths.size();i++){
            maxDistance=findMaximumDistance(paths[i]);
            if(maxDistance.second>md.second){
                pathIndex=i;
                md=maxDistance;
            }
        }
        if(md.second==0) break;
        int placeIndex=md.first;
        vector<SpatialIndex::STPoint> *p=&paths[pathIndex];
        vector<SpatialIndex::STPoint>::const_iterator it1=p->begin();
        vector<SpatialIndex::STPoint>::const_iterator it2=p->end();
        auto indexIt=paths.begin()+pathIndex;
        vector<SpatialIndex::STPoint> path1(it1,it1+placeIndex+1); //new path l1 from 0 to index
        vector<SpatialIndex::STPoint> path2(it1+placeIndex,it2); // new path l2 from index to last
        paths.insert(indexIt+1,path2);
        indexIt=paths.begin()+pathIndex;
        paths.insert(indexIt+1,path1);
        indexIt=paths.begin()+pathIndex;
        paths.erase(indexIt);
    }
    return paths;
}

std::vector<Trajectory> Trajectory::cuttraj(std::vector<SpatialIndex::STPoint> mask){
    vector<STPoint> seg;
    vector<Trajectory> res;
    if(mask.size()==2){
        res.emplace_back(*this);
        return res;
    }
    auto iter1=m_points.begin();
    auto iter2=mask.begin();
    assert(m_points[0]==mask[0]);
    iter2++;
    for(;iter2!=mask.end();iter2++){
        seg.clear();
        while(iter1!=m_points.end()&&iter1->m_time<=iter2->m_time){
//            std::cerr<<"placed"<<*iter1<<"\n";
            seg.emplace_back(*iter1);
            iter1++;
        }
        if(seg.size()>=2) {
            if(!res.empty()){
                if(seg.front().m_time!=res.back().m_endTime()){
                    std::cerr<<"heelo\n";
                    for(const auto &t:mask) std::cerr<<t<<"\n";
                    std::cerr<<*this<<"\n";
                }
            }
            res.emplace_back(Trajectory(seg));
            iter1--;
//            std::cerr<<iter1->m_startTime<<" "<<iter2->m_startTime;
        }
    }


//    for(int i=0;i<res.size();i++){
//        std::cerr<<"part "<<i<<"is\n";
//        std::cerr<<res[i];
//    }
    return res;
}

std::vector<Trajectory> Trajectory::getRDPSegments(double len) const {
    int segNum=std::ceil((m_endTime()-m_startTime())/len);
    std::vector<Trajectory> res;
    auto m=simplifyWithRDPN(m_points,std::min(int(m_points.size()-1),segNum));
    for(auto &pts:m){
        if(pts.size()<2){
            std::cerr<<"error on getting segments with "<<len<<"and \n"<<*this;
            throw Tools::IllegalStateException("bad RDPN");
        }
        res.emplace_back(Trajectory(pts));
    }
    return res;
}

std::vector<Trajectory> Trajectory::getSegments(double len) const {

    tjstat->bt=len;
    return getGlobalSegmentsCut(len);

//    return getHybridSegments(len);
//    auto m=getStaticSegments(len*sqrt(double(segNum)));
//    for(auto &traj:m){
//        auto seg=traj.getRDPSegments(len);
//        res.insert(res.end(),seg.begin(),seg.end());
//    }
//    return res;

//    return getStaticSegmentsCut(len);
}

std::vector<Trajectory> Trajectory::getHybridSegments(double len) const {
    int segNum=std::ceil((m_endTime()-m_startTime())/len);
    std::vector<Trajectory> res;
    if(segNum==1) {res.emplace_back(*this);return res;}
//    std::cerr<<int(std::floor(std::pow(segNum,1.0/4)))<<"\t";
    auto m=simplifyWithRDPN(m_points,std::min(int(std::sqrt(m_points.size()-1)),int(segNum/5)));
    for(auto &pts:m){
        if(pts.size()<2){
            std::cerr<<"error on getting segments with "<<len<<"and \n"<<*this;
            throw Tools::IllegalStateException("bad RDPN");
        }
        auto seg=Trajectory(pts).getGlobalSegmentsCut(len);
        res.insert(res.end(),seg.begin(),seg.end());
    }
    return res;
}

std::vector<Trajectory> Trajectory::getStaticSegmentsCut(double len) const{
    vector<STPoint> seg;
    vector<Trajectory> res;
    bool fakehead=false,fakeback=false;
    if(m_points.size()<2) {
        throw Tools::IllegalStateException("getStatic:seg with 0 or 1 point");
    }
    double segStart=m_startTime();
    seg.emplace_back(m_points[0]);
    for(int i=1;i<m_points.size();i++){
        if(m_points[i].m_time<segStart+len){
            seg.emplace_back(m_points[i]);
        }
        else if(m_points[i].m_time==segStart+len){
            fakeback=false;
            seg.emplace_back(m_points[i]);
            res.emplace_back(Trajectory(fakehead,fakeback,seg));
            fakehead=false;
            seg.clear();
            seg.emplace_back(m_points[i]);
            segStart=m_points[i].m_time;
        }
        else{
            STPoint mid=STPoint::makemid(m_points[i-1],m_points[i],segStart+len);
            fakeback=true;
            seg.emplace_back(mid);
            res.emplace_back(Trajectory(fakehead,fakeback,seg));
            fakehead=true;
            seg.clear();
            seg.emplace_back(mid);
            segStart=mid.m_time;
            i--;
        }
    }
    if(seg.size()>1){
        res.emplace_back(Trajectory(seg));
        seg.clear();
    }
    return res;
}

std::vector<Trajectory> Trajectory::getGlobalSegmentsCut(double len) const {
    vector<STPoint> seg;
    vector<Trajectory> res;
    bool fakehead=false,fakeback=false;
    if(m_points.size()<2) {
        throw Tools::IllegalStateException("getStatic:seg with 0 or 1 point");
    }

    double segStart=double(int(tjstat->mint));
    seg.emplace_back(m_points[0]);
    while(segStart+len<m_points[0].m_time+1e-7) segStart+=len;
    for(int i=1;i<m_points.size();i++){
        if(m_points[i].m_time<segStart+len){
            seg.emplace_back(m_points[i]);
        }
        else if(fabs(m_points[i].m_time-segStart-len)<1e-7){
            fakeback=false;
            seg.emplace_back(m_points[i]);
            res.emplace_back(Trajectory(fakehead,fakeback,seg));
            fakehead=false;
            seg.clear();
            seg.emplace_back(m_points[i]);
            segStart+=len;
        }
        else{
            STPoint mid=STPoint::makemid(m_points[i-1],m_points[i],segStart+len);
            fakeback=true;
            seg.emplace_back(mid);
            res.emplace_back(Trajectory(fakehead,fakeback,seg));
            fakehead=true;
            seg.clear();
            seg.emplace_back(mid);
            segStart+=len;
            i--;
        }
    }
    if(seg.size()>1){
        res.emplace_back(Trajectory(fakehead,true,seg));
        seg.clear();
    }
    return res;
}


std::vector<Trajectory> Trajectory::getStaticSegments(double len) const{
    vector<STPoint> seg;
    vector<Trajectory> res;
    if(m_points.size()<2) {
        throw Tools::IllegalStateException("getStatic:seg with 0 or 1 point");
    }
    double segStart=m_startTime();
    seg.emplace_back(m_points[0]);
    for(int i=1;i<m_points.size();i++){
        seg.emplace_back(m_points[i]);
        if(i==m_points.size()-1){
            res.emplace_back(Trajectory(seg));
            seg.clear();
            break;
        }
        if(m_points[i+1].m_time-segStart>=len){
            res.emplace_back(Trajectory(seg));
            seg.clear();
            seg.emplace_back(m_points[i]);
            segStart=m_points[i].m_time;
        }
    }
    return res;
}



std::vector<Trajectory> Trajectory::getFixedSegments(int len) const {
    if(len<2) len=2;
    vector<STPoint> seg;
    vector<Trajectory> res;
    if(m_points.size()<2) {
        throw Tools::IllegalStateException("getFixed:seg with 0 or 1 point");
    }
    for(int i=0;i<m_points.size();i++){
        seg.emplace_back(m_points[i]);
        if(i==m_points.size()-1){
            res.emplace_back(Trajectory(seg));
            seg.clear();
            break;
        }
        if(seg.size()>=len){
            res.emplace_back(Trajectory(seg));
            seg.clear();
            seg.emplace_back(m_points[i]);
        }
    }
    return res;
}


std::vector<Trajectory> Trajectory::getItself() const {
    vector<Trajectory> res;
    res.emplace_back(*this);
    return res;
}



double Trajectory::getArea() const{ return 0;}

inline double theF(double c1,double c2,double c3,double c4,double t){
    //the c4 should be the length of that time period
    double delta=4*c1*c3-c2*c2;
    if(delta<=1e-7){
        return (2*c1*t+c2)*sqrt(std::max(0.0,c1*t*t+c2*t+c3))/4/c1/c4;
    }
    else {
        return asinh((2 * c1 * t + c2) / sqrt(delta)) * delta / 8 / c1 / sqrt(c1) / c4
               + (2 * c1 * t + c2) * sqrt(std::max(0.0,c1 * t * t + c2 * t + c3)) / 4 / c1 / c4;
    }
}
inline double theD(double c1,double c2,double c3,double c4,double t){
    //the c4 should be the length of that time period
    return sqrt(c1*t*t+c2*t+c3)/c4;
}
inline double theDdd(double c1,double c2,double c3,double c4,double t){
    return 0;
}

inline double Trajectory::line2lineIED(const SpatialIndex::STPoint &p1s, const SpatialIndex::STPoint &p1e,
                                const SpatialIndex::STPoint &p2s, const SpatialIndex::STPoint &p2e) {
    if(p1s.m_time!=p2s.m_time|p1e.m_time!=p2e.m_time)
        throw Tools::IllegalStateException("line2lineIED: time period not the same");
    double ts = p1s.m_time, te = p1e.m_time;
    double dxs=p1s.m_pCoords[0]-p2s.m_pCoords[0];
    double dys=p1s.m_pCoords[1]-p2s.m_pCoords[1];
    double dxe=p1e.m_pCoords[0]-p2e.m_pCoords[0];
    double dye=p1e.m_pCoords[1]-p2e.m_pCoords[1];
    double c1=sq(dxs-dxe)+sq(dys-dye),
            c2=2*((dxe*ts-dxs*te)*(dxs-dxe)+(dye*ts-dys*te)*(dys-dye)),
            c3=sq(dxe*ts-dxs*te)+sq(dye*ts-dys*te),
            c4=te-ts;
    if(c1<1e-7){
        return std::sqrt(sq(dxs)+sq(dys))*c4;
    }else{
            return (theF(c1,c2,c3,c4,te)-theF(c1,c2,c3,c4,ts));
    }
}

inline double Trajectory::line2lineMinSED(const SpatialIndex::STPoint &p1s, const SpatialIndex::STPoint &p1e,
                                const SpatialIndex::STPoint &p2s, const SpatialIndex::STPoint &p2e) {
    if(p1s.m_time!=p2s.m_time|p1e.m_time!=p2e.m_time)
        throw Tools::IllegalStateException("line2lineMinSED: time period not the same");
    double ts = p1s.m_time, te = p1e.m_time;
    double dxs=p1s.m_pCoords[0]-p2s.m_pCoords[0];
    double dys=p1s.m_pCoords[1]-p2s.m_pCoords[1];
    double dxe=p1e.m_pCoords[0]-p2e.m_pCoords[0];
    double dye=p1e.m_pCoords[1]-p2e.m_pCoords[1];
    double c1=sq(dxs-dxe)+sq(dys-dye),
            c2=2*((dxe*ts-dxs*te)*(dxs-dxe)+(dye*ts-dys*te)*(dys-dye)),
            c3=sq(dxe*ts-dxs*te)+sq(dye*ts-dys*te),
            c4=te-ts;
    if(c1<1e-7){
        return std::sqrt(sq(dxs)+sq(dys));
    }else{
        double middle=-c2/c1/2;
        if(middle>ts&&middle<te){
            return sqrt((4*c1*c3-c2*c2)/4/c1)/c4;
        }
        else{
            return std::min(std::sqrt(sq(dxs)+sq(dys)),std::sqrt(sq(dxe)+sq(dye)));
        }
    }
}

inline int Trajectory::getPhase(const SpatialIndex::Region &r,const Point &p1,const Point &p2){
    // 7 8 9
    // 4 5 6
    // 1 2 3
    if(p1.m_dimension!=2)
        throw Tools::NotSupportedException("Trajectory::getPhase:only for dimension 2");
    double xd1=r.m_pLow[0],xd2=r.m_pHigh[0],yd1=r.m_pLow[1],yd2=r.m_pHigh[1];
    int res=0;
    if(p1.m_pCoords[0]<=xd2+1e-7&&p2.m_pCoords[0]<=xd2+1e-7&&p1.m_pCoords[0]>=xd1-1e-7&&p2.m_pCoords[0]>=xd1-1e-7)
        res+=2;
    else if(p1.m_pCoords[0]<=xd1+1e-7&&p2.m_pCoords[0]<=xd1+1e-7)
        res+=1;
    else if(p1.m_pCoords[0]>=xd2-1e-7&&p2.m_pCoords[0]>=xd2-1e-7)
        res+=3;
    else return -1;
    if(p1.m_pCoords[1]<=yd2+1e-7&&p2.m_pCoords[1]<=yd2+1e-7&&p1.m_pCoords[1]>=yd1-1e-7&&p2.m_pCoords[1]>=yd1-1e-7)
        res+=3;
    else if(p1.m_pCoords[1]<=yd1+1e-7&&p2.m_pCoords[1]<=yd1+1e-7)
        res+=0;
    else if(p1.m_pCoords[1]>=yd2-1e-7&&p2.m_pCoords[1]>=yd2-1e-7)
        res+=6;
    else return -1;
    return res;
}

STPoint* cutByLine(const SpatialIndex::STPoint &ps, const SpatialIndex::STPoint &pe,double value,int axis){
    int otheraxis=1-axis;
    double axisv1=ps.m_pCoords[axis],axisv2=pe.m_pCoords[axis];
    if((axisv1<value)==(axisv2<value)||std::fabs(axisv1-value)<1e-7||std::fabs(axisv2-value)<1e-7) //need no cut
        return nullptr;
    else {
        double d = std::fabs(ps.m_pCoords[axis] - pe.m_pCoords[axis]);
        double d1 = std::fabs(ps.m_pCoords[axis] - value) / d;
        double d2 = std::fabs(pe.m_pCoords[axis] - value) / d;
        //get p=d2*ps+d1*pe
        double xyt[3];
        if(axis==0) {
            xyt[0]=value;
            xyt[1]=d2 * ps.m_pCoords[otheraxis] + d1 * pe.m_pCoords[otheraxis];
            xyt[2]=d2 * ps.m_time + d1 * pe.m_time;
        }else{
            xyt[0]=d2 * ps.m_pCoords[otheraxis] + d1 * pe.m_pCoords[otheraxis];
            xyt[1]=value;
            xyt[2]=d2 * ps.m_time + d1 * pe.m_time;
        }
//        Tools::SmartPointer<STPoint> sp(new STPoint(xyt, xyt[2], 2));
        auto res=new STPoint(xyt, xyt[2], 2);
        return res;
//        return sp.get();
    }
}
std::vector<std::pair<STPoint,STPoint>> Trajectory::cutByPhase(const SpatialIndex::STPoint &ps, const SpatialIndex::STPoint &pe,
                                                     const SpatialIndex::Region &r){
    double xd1=r.m_pLow[0],xd2=r.m_pHigh[0],yd1=r.m_pLow[1],yd2=r.m_pHigh[1];
    std::vector<std::pair<STPoint,STPoint>> res;
    std::vector<std::pair<STPoint,STPoint>> tmp;
    res.emplace_back(std::make_pair(ps,pe));
    for(const auto &line:res){
        STPoint* stp=cutByLine(line.first,line.second,xd1,0);
        if(stp!= nullptr){
            tmp.emplace_back(std::make_pair(line.first,*stp));
            tmp.emplace_back(std::make_pair(*stp,line.second));
            delete stp;
        }
        else{
            tmp.emplace_back(line);
        }
    }
    res=tmp;tmp.clear();
    for(const auto &line:res){
        STPoint *stp=cutByLine(line.first,line.second,xd2,0);
        if(stp!= nullptr){
            tmp.emplace_back(std::make_pair(line.first,*stp));
            tmp.emplace_back(std::make_pair(*stp,line.second));
            delete stp;
        }
        else{
            tmp.emplace_back(line);
        }
    }
    res=tmp;tmp.clear();
    for(const auto &line:res){
        STPoint *stp=cutByLine(line.first,line.second,yd1,1);
        if(stp!= nullptr){
            tmp.emplace_back(std::make_pair(line.first,*stp));
            tmp.emplace_back(std::make_pair(*stp,line.second));
            delete stp;
        }
        else{
            tmp.emplace_back(line);
        }
    }
    res=tmp;tmp.clear();
    for(const auto &line:res){
        STPoint *stp=cutByLine(line.first,line.second,yd2,1);
        if(stp!= nullptr){
            tmp.emplace_back(std::make_pair(line.first,*stp));
            tmp.emplace_back(std::make_pair(*stp,line.second));
            delete stp;
        }
        else{
            tmp.emplace_back(line);
        }
    }
    res=tmp;tmp.clear();
//    for(auto seg:res){
//        cout<<seg.first<<"\n"<<seg.second<<"\n";
//    }
    return res;
}

inline double Trajectory::line2MBRMinSED_impl(const SpatialIndex::STPoint &ps, const SpatialIndex::STPoint &pe,
                                       const SpatialIndex::Region &r, int sr) {
    double ts = ps.m_time, te = pe.m_time;
//    if(std::fabs(te-ts)<1e-7) return 0;
    Region r_2d=Region(r.m_pLow,r.m_pHigh,2);
    double res;
    if (sr == 5) return 0;
    else if (sr % 2 == 0){
        return std::min(ps.getMinimumDistance(r_2d),pe.getMinimumDistance(r_2d));
    }
    else {
        double px, py;
        if (sr == 1 || sr == 7) px = r.m_pLow[0];
        else px = r.m_pHigh[0];
        if (sr == 1 || sr == 3) py = r.m_pLow[1];
        else py = r.m_pHigh[1];
        double coord[2]={px,py};
        return line2lineMinSED(ps, pe, STPoint(coord, ts, 2), STPoint(coord, te, 2));
    }
}

inline double Trajectory::line2MBRMinSED(const SpatialIndex::STPoint &ps, const SpatialIndex::STPoint &pe,
                                  const SpatialIndex::Region &r) {
    assert(r.m_pLow[m_dimension]<=ps.m_time);
    assert(r.m_pHigh[m_dimension]>=pe.m_time);
    int sr = getPhase(r, ps, pe);
    if (sr > 0) {
        return line2MBRMinSED_impl(ps, pe, r, sr);
    } else {
        double min = 1e300;
        auto part = cutByPhase(ps, pe, r);
        for (const auto &p:part) {
            int tmpsr = getPhase(r, p.first, p.second);
            double tmpres = line2MBRMinSED_impl(p.first, p.second, r, tmpsr);
//            cout<<tmpsr<<" "<<tmpres<<"\n";
            min=std::min(min,tmpres);
        }
        return min;
    }
}

inline double Trajectory::line2MBRIED_impl(const SpatialIndex::STPoint &ps, const SpatialIndex::STPoint &pe,
                                    const SpatialIndex::Region &r, int sr) {
    double ts = ps.m_time, te = pe.m_time;
    if(std::fabs(te-ts)<1e-7) return 0;
    Region r_2d=Region(r.m_pLow,r.m_pHigh,2);
    double res;
    if (sr == 5) return 0;
    else if (sr % 2 == 0){
        return 0.5 * (ps.getMinimumDistance(r_2d) + pe.getMinimumDistance(r_2d)) * (te - ts);
    }
    else {
        double px, py;
        if (sr == 1 || sr == 7) px = r.m_pLow[0];
        else px = r.m_pHigh[0];
        if (sr == 1 || sr == 3) py = r.m_pLow[1];
        else py = r.m_pHigh[1];
        double coord[2]={px,py};
        return line2lineIED(ps, pe, STPoint(coord, ts, 2), STPoint(coord, te, 2));
    }
}

inline double Trajectory::line2MBRMaxSED(const SpatialIndex::STPoint &ps, const SpatialIndex::STPoint &pe,
                                  const SpatialIndex::Region &r) {
    assert(r.m_pLow[m_dimension]<=ps.m_time);
    assert(r.m_pHigh[m_dimension]>=pe.m_time);
    Region mbr2d(r.m_pLow,r.m_pHigh,r.m_dimension-1);
    return std::max(ps.getMinimumDistance(mbr2d),pe.getMinimumDistance(mbr2d));
}

inline double Trajectory::line2MBRDistance(const SpatialIndex::STPoint &ps, const SpatialIndex::STPoint &pe,
                                    const SpatialIndex::Region &r) {
    //the line's time period should be in the MBR's time period
    assert(r.m_pLow[m_dimension]<=ps.m_time);
    assert(r.m_pHigh[m_dimension]>=pe.m_time);
    //check if need cutting
    if(disttype==0) {
        int sr = getPhase(r, ps, pe);
        if (sr > 0) {
            return line2MBRIED_impl(ps, pe, r, sr);
        } else {
            double sum = 0;
            auto part = cutByPhase(ps, pe, r);
            for (const auto &p:part) {
                int tmpsr = getPhase(r, p.first, p.second);
                double tmpres = line2MBRIED_impl(p.first, p.second, r, tmpsr);
//            cout<<tmpsr<<" "<<tmpres<<"\n";
                sum += tmpres;
            }
            return sum;
        }
    }else if(disttype==1){
        return line2MBRMinSED(ps,pe,r);
    }
    else{throw Tools::NotSupportedException("Wrong distance");}
}

inline double mbcArea(double t0,double t3,double rd,double rv,double ts,double te){
    double t1=t0+rd/rv,t2=t3-rd/rv;
    double sum=0;
    double tlow,thigh,rlow,rhigh;
    if(ts<t1){
        tlow=ts;
        thigh=std::min(te,t1);
        rlow=(tlow-t0)*rv;
        rhigh=(thigh-t0)*rv;
        sum+=(thigh-tlow)*0.5*(rlow+rhigh);
    }
    if(ts<t2&&te>t1) sum+=rd*(std::min(t2,te)-std::max(t1,ts));
    if(te>t2){
        tlow=std::max(ts,t2);
        thigh=te;
        rlow=(t3-tlow)*rv;
        rhigh=(t3-thigh)*rv;
        sum+=(thigh-tlow)*0.5*(rlow+rhigh);
    }
    return sum;
}

inline double Trajectory::line2MBCDistance(const SpatialIndex::STPoint &ps, const SpatialIndex::STPoint &pe,
                                    const SpatialIndex::MBC &r) {
    if(disttype==0) {
        double x1 = r.m_pLow[0], y1 = r.m_pLow[1], t1 = r.m_startTime;
        double x2 = r.m_pHigh[0], y2 = r.m_pHigh[1], t2 = r.m_endTime;
        double xts = makemidmacro(x1, t1, x2, t2, ps.m_time), xte = makemidmacro(x1, t1, x2, t2, pe.m_time);
        double yts = makemidmacro(y1, t1, y2, t2, ps.m_time), yte = makemidmacro(y1, t1, y2, t2, pe.m_time);
        double p2Low[2] = {xts, yts}, p2High[2] = {xte, yte};
        double mbcr1 = std::min(std::min(r.m_rd, (ps.m_time - r.m_startTime) * r.m_rv),
                                (r.m_endTime - ps.m_time) * r.m_rv);
        double mbcr2 = std::min(std::min(r.m_rd, (pe.m_time - r.m_startTime) * r.m_rv),
                                (r.m_endTime - pe.m_time) * r.m_rv);
        STPoint p2s(p2Low, ps.m_time, 2), p2e(p2High, pe.m_time, 2);
        double s = line2lineIED(ps, pe, p2s, p2e);
        double sm = mbcArea(r.m_startTime, r.m_endTime, r.m_rd, r.m_rv, ps.m_time, pe.m_time);
        if (ps.getMinimumDistance(p2s) > mbcr1 && pe.getMinimumDistance(p2e) > mbcr2) return s - sm;
        else return 0;
    }else if(disttype==1){
        return std::max(r.getMinimumDistance(ps),r.getMinimumDistance(pe));
    }
    else{throw Tools::NotSupportedException("Wrong distance");}
}

double Trajectory::getMinimumDistance(const IShape& s) const{
    const Trajectory* pTrajectory = dynamic_cast<const Trajectory*>(&s);
    if (pTrajectory != nullptr) return getMinimumDistance(*pTrajectory);

//    const STPoint* ptp = dynamic_cast<const STPoint*>(&s);
//    if (ptp != 0) return getMinimumDistance(*ptp);

    const MBC* pmbc= dynamic_cast<const MBC*>(&s);
    if(pmbc!=nullptr) return getMinimumDistance(*pmbc);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != nullptr) return getMinimumDistance(*pr);

    throw Tools::NotSupportedException(
            "Trajectory::...: Not implemented yet!"
    );
}


double Trajectory::getMinimumDistance(const SpatialIndex::Region &in) const {
        //note here we only calculate the partial distance
        //so just calculate intersection of time period
        //this function is used for inner nodes, not leaf nodes
        double tstart, tend;
        tstart = std::max(m_startTime(), in.m_pLow[in.m_dimension - 1]);
        tend = std::min(m_endTime(), in.m_pHigh[in.m_dimension - 1]);
        if(tstart>=tend) return 1e300;
        fakeTpVector timedTraj(&m_points,tstart,tend);
        double sum = 0,min=0;
        for (int i = 0; i < timedTraj.m_size-1; i++) {
            if(disttype==0){
                double pd = line2MBRDistance(timedTraj[i],timedTraj[i+1],in);
                sum+=pd;
            }
            else if(disttype==1){
                double pd = line2MBRMinSED(timedTraj[i],timedTraj[i+1],in);
                min=std::min(min,pd);
            }
        }
        if(disttype==0)
            return sum;
        else if(disttype==1)
            return min;
        else{throw Tools::NotSupportedException("Wrong distance");}
}

double Trajectory::getMaxSED(const SpatialIndex::Region &in) const {
    //used for single trajectory
    Region mbr2d(in.m_pLow,in.m_pHigh,in.m_dimension-1);
    double tstart, tend;
    tstart = std::max(m_startTime(), in.m_pLow[in.m_dimension - 1]);
    tend = std::min(m_endTime(), in.m_pHigh[in.m_dimension - 1]);
    if (tstart >= tend) return 1e300;
    fakeTpVector timedTraj(&m_points, tstart, tend);
    double max = 0,pd;
    for (int i = 0; i < timedTraj.m_size; i++) {
        pd = timedTraj[i].getMinimumDistance(mbr2d);
        max = std::max(max, pd);
    }
    return max;
}

double Trajectory::getMaxSED(const SpatialIndex::MBC &in) const {
    //used for single trajectory
    double tstart, tend;
    tstart = std::max(m_startTime(), in.m_startTime);
    tend = std::min(m_endTime(), in.m_endTime);
    if (tstart >= tend) return 1e300;
    fakeTpVector timedTraj(&m_points, tstart, tend);
    double max = 0,pd;
    for (int i = 0; i < timedTraj.m_size; i++) {
        pd = in.getMinimumDistance(timedTraj[i]);
        max = std::max(max, pd);
    }
    STPoint tmpp;
    double tmpt=in.m_startTime+in.m_rd/in.m_rv;
    if(tmpt>m_startTime()&&tmpt<m_endTime()){
        tmpp=getPointAtTime(tmpt);
        pd = in.getMinimumDistance(tmpp);
        max = std::max(max, pd);
    }
    tmpt=in.m_endTime-in.m_rd/in.m_rv;
    if(tmpt>m_startTime()&&tmpt<m_endTime()){
        tmpp=getPointAtTime(tmpt);
        pd = in.getMinimumDistance(tmpp);
        max = std::max(max, pd);
    }
    return max;
}

double Trajectory::getMinimumDistance(const SpatialIndex::MBC &in) const {
    assert(m_dimension==in.m_dimension);
    double tstart, tend;
    tstart = std::max(m_startTime(), in.m_startTime);
    tend = std::min(m_endTime(), in.m_endTime);
    if(tstart>=tend) return 1e300;
    fakeTpVector timedTraj(&m_points,tstart,tend);
    if(disttype==0) {
        double sum = 0;
        for (int i = 0; i < timedTraj.m_size - 1; i++) {
            double pd = line2MBCDistance(timedTraj[i], timedTraj[i + 1], in);
            sum += pd;
        }
        return sum;
    }
    else if(disttype==1){
        return getMaxSED(in);
    }
    else{throw Tools::NotSupportedException("Wrong distance");}
}

double Trajectory::getMinimumDistance(const SpatialIndex::STPoint &in) const {
    int last;
    for(last=1;m_points[last].m_time<in.m_time;last++);
    return STPoint::makemid(m_points[last-1],m_points[last],in.m_time).getMinimumDistance(in);
}

double Trajectory::getMinimumDistance(const SpatialIndex::Trajectory &in) const {
    if(m_startTime()>=in.m_endTime()||m_endTime()<=in.m_startTime())
        return 1e300;
    fakeTpVector timedTraj2(&in.m_points,m_startTime(),m_endTime());
    double cut1=timedTraj2[0].m_time,cut2=timedTraj2[timedTraj2.m_size-1].m_time;
    double sum=0;
    double max=0;
    fakeTpVector midTraj(&m_points,cut1,cut2);
    if(m_startTime()<cut1){
        double pd;
        if(disttype==0) {
            pd= getStaticIED(timedTraj2[0].m_pCoords[0], timedTraj2[0].m_pCoords[1], m_startTime(), cut1);
            sum += pd;
        }else if(disttype==1){
            pd= getStaticMaxSED(timedTraj2[0].m_pCoords[0], timedTraj2[0].m_pCoords[1], m_startTime(), cut1);
            max=std::max(max,pd);
        }
    }
    if(m_endTime()>cut2){
        double pd;
        if(disttype==0) {
            pd= getStaticIED(timedTraj2[timedTraj2.m_size - 1].m_pCoords[0],
                              timedTraj2[timedTraj2.m_size - 1].m_pCoords[1], cut2, m_endTime());
            sum += pd;
        }else if(disttype==1){
            pd= getStaticMaxSED(timedTraj2[timedTraj2.m_size - 1].m_pCoords[0],
                              timedTraj2[timedTraj2.m_size - 1].m_pCoords[1], cut2, m_endTime());
            max=std::max(max,pd);
        }
    }
    if(midTraj.m_size!=0) {
        double newtime = midTraj[0].m_time, lasttime = midTraj[0].m_time;
        auto iter1 = midTraj.m_vectorPointer->begin()+midTraj.m_is;
        auto iter2 = timedTraj2.m_vectorPointer->begin()+timedTraj2.m_is;
        STPoint lastp1 = midTraj[0], lastp2 = timedTraj2[0], newp1, newp2;
        newp1.makeInfinite(2);newp2.makeInfinite(2);
        if(disttype==1) {max=std::max(max,iter1->getMinimumDistance(*iter2));}
        while (lasttime != timedTraj2[timedTraj2.m_size-1].m_time) {
            if ((iter1 + 1)->m_time == (iter2 + 1)->m_time) {
                newtime = (iter1 + 1)->m_time;
                newp1 = *(iter1 + 1);
                newp2 = *(iter2 + 1);
                iter1++;
                iter2++;
            } else if ((iter1 + 1)->m_time < (iter2 + 1)->m_time) {
                newtime = (iter1 + 1)->m_time;
                newp1 = *(iter1 + 1);
                double x=makemidmacro(iter2->m_pCoords[0],iter2->m_time,(iter2 + 1)->m_pCoords[0],(iter2 + 1)->m_time,newtime);
                double y=makemidmacro(iter2->m_pCoords[1],iter2->m_time,(iter2 + 1)->m_pCoords[1],(iter2 + 1)->m_time,newtime);
                newp2.m_pCoords[0]=x;
                newp2.m_pCoords[1]=y;
                newp2.m_time=newtime;
                iter1++;
            } else {
                newtime = (iter2 + 1)->m_time;
                double x=makemidmacro(iter1->m_pCoords[0],iter1->m_time,(iter1 + 1)->m_pCoords[0],(iter1 + 1)->m_time,newtime);
                double y=makemidmacro(iter1->m_pCoords[1],iter1->m_time,(iter1 + 1)->m_pCoords[1],(iter1 + 1)->m_time,newtime);
                newp1.m_pCoords[0]=x;
                newp1.m_pCoords[1]=y;
                newp1.m_time=newtime;
                newp2 = *(iter2 + 1);
                iter2++;
            }
            lasttime = newtime;
            if(disttype==0) {
                double pd = line2lineIED(lastp1, newp1, lastp2, newp2);
//                std::cerr<<"distance\n"<<lastp1<<"\t"<<newp1<<"\n"<<lastp2<<"\t"<<newp2<<"\n"<<pd<<"\n";
                sum += pd;
            }
            else if(disttype==1){
                max=std::max(newp1.getMinimumDistance(newp2),max);
            }
            else{throw Tools::NotSupportedException("Wrong distance");}
            lastp1 = newp1;
            lastp2 = newp2;
        }
    }
    if(disttype==0) {
        return sum;
    }
    else if(disttype==1){
        return max;
    }
    else{throw Tools::NotSupportedException("Wrong distance");}
}

double Trajectory::getMinimumDistanceInner(const SpatialIndex::Trajectory &in) const {
    if(m_startTime()>=in.m_endTime()||m_endTime()<=in.m_startTime())
        return 1e300;
    fakeTpVector timedTraj2(&in.m_points,m_startTime(),m_endTime());
    double cut1=timedTraj2[0].m_time,cut2=timedTraj2[timedTraj2.m_size-1].m_time;
    double sum=0;
    double max=0;
    fakeTpVector midTraj(&m_points,cut1,cut2);
    if(midTraj.m_size!=0) {
        double newtime = midTraj[0].m_time, lasttime = midTraj[0].m_time;
        auto iter1 = midTraj.m_vectorPointer->begin()+midTraj.m_is;
        auto iter2 = timedTraj2.m_vectorPointer->begin()+timedTraj2.m_is;
        STPoint lastp1 = midTraj[0], lastp2 = timedTraj2[0], newp1, newp2;
        newp1.makeInfinite(2);newp2.makeInfinite(2);
        if(disttype==1) {max=std::max(max,iter1->getMinimumDistance(*iter2));}
        while (lasttime != timedTraj2[timedTraj2.m_size-1].m_time) {
            if ((iter1 + 1)->m_time == (iter2 + 1)->m_time) {
                newtime = (iter1 + 1)->m_time;
                newp1 = *(iter1 + 1);
                newp2 = *(iter2 + 1);
                iter1++;
                iter2++;
            } else if ((iter1 + 1)->m_time < (iter2 + 1)->m_time) {
                newtime = (iter1 + 1)->m_time;
                newp1 = *(iter1 + 1);
                double x=makemidmacro(iter2->m_pCoords[0],iter2->m_time,(iter2 + 1)->m_pCoords[0],(iter2 + 1)->m_time,newtime);
                double y=makemidmacro(iter2->m_pCoords[1],iter2->m_time,(iter2 + 1)->m_pCoords[1],(iter2 + 1)->m_time,newtime);
                newp2.m_pCoords[0]=x;
                newp2.m_pCoords[1]=y;
                newp2.m_time=newtime;
                iter1++;
            } else {
                newtime = (iter2 + 1)->m_time;
                double x=makemidmacro(iter1->m_pCoords[0],iter1->m_time,(iter1 + 1)->m_pCoords[0],(iter1 + 1)->m_time,newtime);
                double y=makemidmacro(iter1->m_pCoords[1],iter1->m_time,(iter1 + 1)->m_pCoords[1],(iter1 + 1)->m_time,newtime);
                newp1.m_pCoords[0]=x;
                newp1.m_pCoords[1]=y;
                newp1.m_time=newtime;
                newp2 = *(iter2 + 1);
                iter2++;
            }
            lasttime = newtime;
            if(disttype==0) {
                double pd = line2lineIED(lastp1, newp1, lastp2, newp2);
//                std::cerr<<"distance\n"<<lastp1<<"\t"<<newp1<<"\n"<<lastp2<<"\t"<<newp2<<"\n"<<pd<<"\n";
                sum += pd;
            }
            else if(disttype==1){
                max=std::max(newp1.getMinimumDistance(newp2),max);
            }
            else{throw Tools::NotSupportedException("Wrong distance");}
            lastp1 = newp1;
            lastp2 = newp2;
        }
    }
    if(disttype==0) {
        return sum;
    }
    else if(disttype==1){
        return max;
    }
    else{throw Tools::NotSupportedException("Wrong distance");}
}

double ldd(double d,double v,double dt){
    if(d+v*dt>0) return dt*(d+v*dt/2);
    else return d*d/2*std::fabs(v);
}
DISTE Trajectory::getFrontIED(const SpatialIndex::Region sbr, double MaxVelocity) const {
    double opti=0,pessi=0;
    double ints=sbr.m_pLow[m_dimension];
    if(bUsingVerboseLB){
        fakeTpVector frontTraj(&m_points,m_startTime(),ints);
        for(int i=0;i<frontTraj.m_size-1;i++) {
            sbr.m_pLow[m_dimension] = frontTraj[i].m_time;
            sbr.m_pHigh[m_dimension] = frontTraj[i + 1].m_time;
            double pd = line2MBRDistance(frontTraj[i], frontTraj[i + 1], sbr);
            double r1 = (ints - frontTraj[i].m_time) * MaxVelocity;
            double r2 = (ints - frontTraj[i + 1].m_time) * MaxVelocity;
            double minus = 0.5 * (frontTraj[i + 1].m_time - frontTraj[i].m_time) * (r1 + r2);
            opti += std::max(0.0, pd - minus);
            pessi += pd+minus;
        }
    }else{
        double ds = pos(ints).getMinimumDistance(Region(sbr.m_pLow,sbr.m_pHigh,sbr.m_dimension-1));
        opti= ldd(ds,-2*MaxVelocity,ints-m_startTime());
        pessi= ldd(ds,2*MaxVelocity,ints-m_startTime());
    }
    return DISTE(opti,pessi,0,true);
}
DISTE Trajectory::getFrontIED(double x, double y, double ints, double MaxVelocity) const {
    double opti=0,pessi=0;
    if(bUsingVerboseLB){
        fakeTpVector frontTraj(&m_points, m_startTime(), ints);
        STPoint sPoint(x,y,ints),sPoint2(x,y,ints);
        for (int i = 0; i < frontTraj.m_size - 1; i++) {
            sPoint.m_time = frontTraj[i].m_time;
            sPoint2.m_time = frontTraj[i + 1].m_time;
            double pd = line2lineIED(sPoint, sPoint2, frontTraj[i], frontTraj[i + 1]);
            double r1 = (ints-frontTraj[i].m_time) * MaxVelocity;
            double r2 = (ints-frontTraj[i + 1].m_time) * MaxVelocity;
            double minus = 0.5 * (frontTraj[i + 1].m_time - frontTraj[i].m_time) * (r1 + r2);
            opti += std::max(0.0, pd - minus);
            pessi += pd+minus;
        }
    }else{
        double xy[]={x,y};
        double ds = pos(ints).getMinimumDistance(Point(xy,2));
        opti= ldd(ds,-2*MaxVelocity,ints-m_startTime());
        pessi= ldd(ds,2*MaxVelocity,ints-m_startTime());
    }
    return DISTE(opti,pessi,0,true);
}
DISTE Trajectory::getBackIED(const SpatialIndex::Region ebr, double MaxVelocity) const {
    double opti=0,pessi=0;
    double inte=ebr.m_pHigh[m_dimension];
    if(bUsingVerboseLB){
        fakeTpVector backTraj(&m_points,inte,m_endTime());
        for(int i=0;i<backTraj.m_size-1;i++){
            ebr.m_pLow[m_dimension]=backTraj[i].m_time;
            ebr.m_pHigh[m_dimension]=backTraj[i+1].m_time;
            double pd=line2MBRDistance(backTraj[i],backTraj[i+1],ebr);
            double r1=(backTraj[i].m_time-inte)*MaxVelocity;
            double r2=(backTraj[i+1].m_time-inte)*MaxVelocity;
            double minus = 0.5*(backTraj[i+1].m_time-backTraj[i].m_time)*(r1+r2);
            opti+=std::max(0.0,pd-minus);
            pessi += pd+minus;
        }
    }else{
        double de = pos(inte).getMinimumDistance(Region(ebr.m_pLow,ebr.m_pHigh,ebr.m_dimension-1));
        opti= ldd(de,-2*MaxVelocity,m_endTime()- inte);
        pessi= ldd(de,2*MaxVelocity,m_endTime()- inte);
    }
    return DISTE(opti,pessi,0,true);
}
DISTE Trajectory::getBackIED(double x, double y,  double inte, double MaxVelocity) const {
    double opti=0,pessi=0;
    if(bUsingVerboseLB){
        fakeTpVector backTraj(&m_points, inte, m_endTime());
        STPoint ePoint(x,y,inte),ePoint2(x,y,inte);
        for (int i = 0; i < backTraj.m_size - 1; i++) {
            ePoint.m_time = backTraj[i].m_time;
            ePoint2.m_time = backTraj[i + 1].m_time;
            double pd = line2lineIED(ePoint, ePoint2, backTraj[i], backTraj[i + 1]);
            double r1 = (  backTraj[i].m_time-inte) * MaxVelocity;
            double r2 = ( backTraj[i + 1].m_time-inte) * MaxVelocity;
            double minus = 0.5 * (backTraj[i + 1].m_time - backTraj[i].m_time) * (r1 + r2);
            opti += std::max(0.0, pd - minus);
            pessi += pd+minus;
        }
    }else{
        double xy[]={x,y};
        double de = pos(inte).getMinimumDistance(Point(xy,2));
        opti= ldd(de,-2*MaxVelocity,m_endTime()- inte);
        pessi= ldd(de,2*MaxVelocity,m_endTime()- inte);
    }
    return DISTE(opti,pessi,0,true);
}

DISTE Trajectory::getMidIED(const SpatialIndex::Region &sbr, const SpatialIndex::Region &ebr,
                            double MaxVelocity,double queryVelocity) {
    double opti=0,pessi=0;
    double ints=sbr.m_pHigh[m_dimension],inte=ebr.m_pLow[m_dimension];
    if(bUsingVerboseLB){
        Region sbr2d(sbr.m_pLow,sbr.m_pHigh,sbr.m_dimension-1),
                ebr2d(ebr.m_pLow,ebr.m_pHigh,ebr.m_dimension-1);
        fakeTpVector timedTraj(&m_points, ints, inte);
        for(int i=0;i<timedTraj.m_size-1;i++){
            double t1=timedTraj[i].m_time,t2=timedTraj[i+1].m_time;
            double l1 = (t1 - ints) * MaxVelocity,
                    l2 = (t2 - ints) * MaxVelocity,
                    h1 = (inte - t1) * MaxVelocity,
                    h2 = (inte - t2) * MaxVelocity;
            double ds1=timedTraj[i].getMinimumDistance(sbr2d),
                    ds2=timedTraj[i+1].getMinimumDistance(sbr2d),
                    de1=timedTraj[i].getMinimumDistance(ebr2d),
                    de2=timedTraj[i+1].getMinimumDistance(ebr2d);
            double a1=l1-ds1,a2=l2-ds2,b1=h1-de1,b2=h2-de2;
            double d1=std::min(a1,b1),d2=std::min(a2,b2);
            double minus,pd;
            if(a2>=b2){
                minus=0.5*(l1+l2)*(t2-t1);
                pd= getStaticIED(sbr, t1, t2);
                opti+=std::max(0.0,pd-minus);
            }else if(b1>=a1){
                minus=0.5*(h1+h2)*(t2-t1);
                pd= getStaticIED(ebr, t1, t2);
                opti+=std::max(0.0,pd-minus);
            }else{
                double vmax;
                if(queryVelocity==-1) vmax=2*MaxVelocity;
                else vmax=MaxVelocity+queryVelocity;
                double t0=(t1+t2+(d2-d1)/(vmax))/2;
                opti+=ldd(d1,-vmax,t0-t1)+ldd(d2,vmax,inte-t2);
            }
        }
    }else{
        double ds = pos(ints).getMinimumDistance(Region(sbr.m_pLow,sbr.m_pHigh,sbr.m_dimension-1));
        double de = pos(inte).getMinimumDistance(Region(ebr.m_pLow,ebr.m_pHigh,ebr.m_dimension-1));
        double to=(ints+inte+(de-ds)/(2*MaxVelocity))/2;
        double tp = (ints+inte+(ds-de)/(2*MaxVelocity))/2;
        opti= ldd(ds,-2*MaxVelocity,to-ints)+ldd(de,2*MaxVelocity,inte-to);
        pessi = ldd(ds,2*MaxVelocity,tp-ints)+ldd(de,-2*MaxVelocity,inte-tp);
    }
    return DISTE(opti,pessi,0,true);
}
DISTE Trajectory::getMidIED(const MBC &sbc, const MBC &ebc,
                            double MaxVelocity,double queryVelocity) {
    double opti=0,pessi=0;
    double ints=sbc.m_endTime,inte=ebc.m_startTime;
    if(bUsingVerboseLB){
        STPoint sPoint(sbc.m_pHigh,sbc.m_endTime,m_dimension),ePoint(ebc.m_pLow,ebc.m_startTime,m_dimension);
        fakeTpVector timedTraj(&m_points, ints, inte);
        for(int i=0;i<timedTraj.m_size-1;i++){
            double t1=timedTraj[i].m_time,t2=timedTraj[i+1].m_time;
            double l1 = (t1 - ints) * MaxVelocity,
                    l2 = (t2 - ints) * MaxVelocity,
                    h1 = (inte - t1) * MaxVelocity,
                    h2 = (inte - t2) * MaxVelocity;
            double ds1=timedTraj[i].getMinimumDistance(sPoint),
                    ds2=timedTraj[i+1].getMinimumDistance(sPoint),
                    de1=timedTraj[i].getMinimumDistance(ePoint),
                    de2=timedTraj[i+1].getMinimumDistance(ePoint);
            double a1=l1-ds1,a2=l2-ds2,b1=h1-de1,b2=h2-de2;
            double d1=std::min(a1,b1),d2=std::min(a2,b2);
            double minus,pd;
            if(a2>=b2){
                minus=0.5*(l1+l2)*(timedTraj[i + 1].m_time-timedTraj[i].m_time);
                pd= getStaticIED(sPoint.m_pCoords[0],sPoint.m_pCoords[1], t1, t2);
                opti+=std::max(0.0,pd-minus);
            }else if(b1>=a1){
                minus=0.5*(h1+h2)*(timedTraj[i + 1].m_time-timedTraj[i].m_time);
                pd= getStaticIED(ePoint.m_pCoords[0],ePoint.m_pCoords[1], t1, t2);
                opti+=std::max(0.0,pd-minus);
            }else{
                double vmax;
                if(queryVelocity==-1) vmax=2*MaxVelocity;
                else vmax=MaxVelocity+queryVelocity;
                double t0=(t1+t2+(d2-d1)/(vmax))/2;
                opti+=ldd(d1,-vmax,t0-t1)+ldd(d2,vmax,inte-t2);
            }
        }
    }else{
        double ds = pos(ints).getMinimumDistance(Point(sbc.m_pLow,2));
        double de = pos(inte).getMinimumDistance(Point(ebc.m_pLow,2));
        double to=(ints+inte+(de-ds)/(2*MaxVelocity))/2;
        double tp = (ints+inte+(ds-de)/(2*MaxVelocity))/2;
        opti=ldd(ds,-2*MaxVelocity,to-ints)+ldd(de,2*MaxVelocity,inte-to);
        pessi = ldd(ds,2*MaxVelocity,tp-ints)+ldd(de,-2*MaxVelocity,inte-tp);
    }
    return DISTE(opti,pessi,0,true);
}

DISTE Trajectory::getMidIED(const STPoint &sPoint, const STPoint &ePoint,
                            double MaxVelocity,double queryVelocity) {
    double opti=0,pessi=0;
    double ints=sPoint.m_time,inte=ePoint.m_time;
    if(bUsingVerboseLB){
        fakeTpVector timedTraj(&m_points, ints, inte);
        for(int i=0;i<timedTraj.m_size-1;i++){
            double t1=timedTraj[i].m_time,t2=timedTraj[i+1].m_time;
            double l1 = (t1 - ints) * MaxVelocity,
                    l2 = (t2 - ints) * MaxVelocity,
                    h1 = (inte - t1) * MaxVelocity,
                    h2 = (inte - t2) * MaxVelocity;
            double ds1=timedTraj[i].getMinimumDistance(sPoint),
                    ds2=timedTraj[i+1].getMinimumDistance(sPoint),
                    de1=timedTraj[i].getMinimumDistance(ePoint),
                    de2=timedTraj[i+1].getMinimumDistance(ePoint);
            double a1=l1-ds1,a2=l2-ds2,b1=h1-de1,b2=h2-de2;
            double d1=std::min(a1,b1),d2=std::min(a2,b2);
            double minus,pd;
            if(a2>=b2){
                minus=0.5*(l1+l2)*(timedTraj[i + 1].m_time-timedTraj[i].m_time);
                pd= getStaticIED(sPoint.m_pCoords[0],sPoint.m_pCoords[1], t1, t2);
                opti+=std::max(0.0,pd-minus);
            }else if(b1>=a1){
                minus=0.5*(h1+h2)*(timedTraj[i + 1].m_time-timedTraj[i].m_time);
                pd= getStaticIED(ePoint.m_pCoords[0],ePoint.m_pCoords[1], t1, t2);
                opti+=std::max(0.0,pd-minus);
            }else{
                double vmax;
                if(queryVelocity==-1) vmax=2*MaxVelocity;
                else vmax=MaxVelocity+queryVelocity;
                double t0=(t1+t2+(d2-d1)/(vmax))/2;
                opti+=ldd(d1,-vmax,t0-t1)+ldd(d2,vmax,inte-t2);
                //            double x1=timedTraj[i].m_pCoords[0],x2=timedTraj[i+1].m_pCoords[0],
                //                    y1=timedTraj[i].m_pCoords[1],y2=timedTraj[i+1].m_pCoords[1];
                //            double psx=sPoint.m_pCoords[0],psy=sPoint.m_pCoords[1],
                //                pex=ePoint.m_pCoords[0],pey=ePoint.m_pCoords[1];
                //            double t0=(2*(pex-psx)*(inte))
            }
        }
    }else{
        double ds = pos(ints).getMinimumDistance(sPoint);
        double de = pos(inte).getMinimumDistance(ePoint);
        double to=(ints+inte+(de-ds)/(2*MaxVelocity))/2;
        double tp = (ints+inte+(ds-de)/(2*MaxVelocity))/2;
        opti=ldd(ds,-2*MaxVelocity,to-ints)+ldd(de,2*MaxVelocity,inte-to);
        pessi = ldd(ds,2*MaxVelocity,tp-ints)+ldd(de,-2*MaxVelocity,inte-tp);
    }
    return DISTE(opti,pessi,0,true);
}

double Trajectory::getStaticIED(double x, double y, double t1, double t2) const {
    double tstart, tend;
    tstart = std::max(m_startTime(), t1);
    tend = std::min(m_endTime(), t2);
    if(tstart>=tend) return 1e300;
    fakeTpVector timedTraj(&m_points,tstart,tend);
    double sum = 0;
    double xy[2]={x,y};
    STPoint ps(xy,0,2),pe(xy,0,2);
    for (int i = 0; i < timedTraj.m_size-1; i++) {
            ps.m_time=timedTraj[i].m_time;
            pe.m_time=timedTraj[i+1].m_time;
            double pd = line2lineIED(timedTraj[i], timedTraj[i + 1], ps, pe);
            sum += pd;
    }
    return sum;
}


double Trajectory::getStaticIED(const SpatialIndex::Region in,double ints, double inte) const {
    assert(ints>=m_startTime()&&inte<=m_endTime());
    in.m_pLow[m_dimension]=ints;
    in.m_pHigh[m_dimension]=inte;
    fakeTpVector timedTraj(&m_points,ints,inte);
    double sum = 0;
    for (int i = 0; i < timedTraj.m_size-1; i++) {
        double pd = line2MBRDistance(timedTraj[i],timedTraj[i+1],in);
        sum+=pd;
    }
    return sum;
}
double Trajectory::getStaticMaxSED(double x, double y, double t1, double t2) const {
    double tstart, tend;
    tstart = std::max(m_startTime(), t1);
    tend = std::min(m_endTime(), t2);
    if(tstart>=tend) return 1e300;
    double max=0;
    fakeTpVector timedTraj(&m_points,tstart,tend);
    STPoint p(x,y,0);
    for(int i=0;i<timedTraj.m_size;i++){
        double pd=timedTraj[i].getMinimumDistance(p);
        if(pd>max) max=pd;
    }
    return max;
}


double Trajectory::getStaticMaxSED(const SpatialIndex::Region in,double ints, double inte) const {
    assert(ints>=m_startTime()&&inte<=m_endTime());
    in.m_pLow[m_dimension]=ints;
    in.m_pHigh[m_dimension]=inte;
    fakeTpVector timedTraj(&m_points,ints,inte);
    double max = 0;
    for (int i = 0; i < timedTraj.m_size; i++) {
        double pd = timedTraj[i].getMinimumDistance(in);
        if(pd>max) max=pd;
    }
    return max;
}

#define sign(x) (x>0?1.0:-1.0)
double Trajectory::getInterTime(const STPoint &ps, const STPoint &pe,Region &r,double t_o, double vmax) const {
    //vmax could be positive or negative
    int sr = getPhase(r, ps, pe);
    if (sr > 0) {
        if(sr==5){
            std::cout<<",...";
        }
        assert(sr!=5);//should be handled previously
        if(sr%2==0){
            double d1=ps.getMinimumDistance(r)-(ps.m_time-t_o)*vmax;
            double d2=pe.getMinimumDistance(r)-(pe.m_time-t_o)*vmax;
            assert(sign(d1)*sign(d2)<=0);
            return (d1*pe.m_time+d2*ps.m_time)/(pe.m_time-ps.m_time);
        }
        else{
            double px, py;
            if (sr == 1 || sr == 7) px = r.m_pLow[0];
            else px = r.m_pHigh[0];
            if (sr == 1 || sr == 3) py = r.m_pLow[1];
            else py = r.m_pHigh[1];
            double ts = ps.m_time, te = pe.m_time;
            double dxs=ps.m_pCoords[0]-px;
            double dys=ps.m_pCoords[1]-py;
            double dxe=pe.m_pCoords[0]-px;
            double dye=pe.m_pCoords[1]-py;
            double c1=sq(dxs-dxe)+sq(dys-dye),
                    c2=2*((dxe*ts-dxs*te)*(dxs-dxe)+(dye*ts-dys*te)*(dys-dye)),
                    c3=sq(dxe*ts-dxs*te)+sq(dye*ts-dys*te),
                    c4=te-ts;
            double a=c1-c4*c4*vmax*vmax,
                b=c2+2*c4*c4*vmax*vmax*t_o,
                c=c3-c4*c4*vmax*vmax*t_o*t_o;
            double t=(-b+sign(a)*sign(vmax)*std::sqrt(b*b-4*a*c))/2/a;
            return t;
        }
    }else {
        auto part = cutByPhase(ps, pe, r);
        for (const auto &p:part) {
            double d1=p.first.getMinimumDistance(r)-(p.first.m_time-t_o)*vmax;
            double d2=p.second.getMinimumDistance(r)-(p.second.m_time-t_o)*vmax;
            if(sign(d1)*sign(d2)<0) return getInterTime(p.first,p.second,r,t_o,vmax);
        }
        throw Tools::IllegalStateException("getInterTime");
    }
}

double Trajectory:: getNodeMinimumDistance(const SpatialIndex::Region &in,double MaxVelocity) const {
    if(m_startTime()>=in.m_pHigh[in.m_dimension-1]||m_endTime()<=in.m_pLow[in.m_dimension-1]) return 1e300;
    if(disttype==0){
        double ints=std::max(m_startTime(),in.m_pLow[in.m_dimension-1]),
                inte=std::min(m_endTime(),in.m_pHigh[in.m_dimension-1]);
        Region copy(in);
        copy.m_pLow[copy.m_dimension-1]=ints;
        copy.m_pHigh[copy.m_dimension-1]=inte;
        double min=1e300;
        fakeTpVector timedTraj(&m_points,ints,inte);
        for (int i = 0; i < timedTraj.m_size-1; i++) {
            double pd = line2MBRMinSED(timedTraj[i],timedTraj[i+1],copy);
            min=std::min(min,pd);
        }
        if(min>=1e300) return 0;
        return min*(m_endTime()-m_startTime());
    }
    else
        return getMinimumDistance(in);
}

double Trajectory::getLeafMinimumDistance(const SpatialIndex::Region &in, double MaxVelocity) const {
    if(disttype==0) {
        double sum=0;
        sum = getMinimumDistance(in);//midTraj
        //frontTraj
        Region sembr=in;
        if(in.m_pLow[m_dimension]>m_startTime()){
            fakeTpVector frontTraj(&m_points,m_startTime(),in.m_pLow[m_dimension]);
            for(int i=0;i<frontTraj.m_size-1;i++){
                sembr.m_pLow[m_dimension]=frontTraj[i].m_time;
                sembr.m_pHigh[m_dimension]=frontTraj[i+1].m_time;
                double pd=line2MBRDistance(frontTraj[i],frontTraj[i+1],sembr);
                double r1=(in.m_pLow[m_dimension]-frontTraj[i].m_time)*MaxVelocity;
                double r2=(in.m_pLow[m_dimension]-frontTraj[i+1].m_time)*MaxVelocity;
                double minus = 0.5*(frontTraj[i+1].m_time-frontTraj[i].m_time)*(r1+r2);
                sum+=std::max(0.0,pd-minus);
            }
        }
        //backTraj
        if(in.m_pHigh[m_dimension]<m_endTime()){
            fakeTpVector backTraj(&m_points,in.m_pHigh[m_dimension],m_endTime());
            for(int i=0;i<backTraj.m_size-1;i++){
                sembr.m_pLow[m_dimension]=backTraj[i].m_time;
                sembr.m_pHigh[m_dimension]=backTraj[i+1].m_time;
                double pd=line2MBRDistance(backTraj[i],backTraj[i+1],sembr);
                double r1=(backTraj[i].m_time-in.m_pHigh[m_dimension])*MaxVelocity;
                double r2=(backTraj[i+1].m_time-in.m_pHigh[m_dimension])*MaxVelocity;
                double minus = 0.5*(backTraj[i+1].m_time-backTraj[i].m_time)*(r1+r2);
                sum+=std::max(0.0,pd-minus);
            }
        }
        return sum;
    }
    else if(disttype==1){
        double min=getMaxSED(in);
        return std::max(0.0,min);
    }else
        throw Tools::NotSupportedException("Wrong distance");
}

double Trajectory::getLeafMinimumDistance(const SpatialIndex::MBC &in, double MaxVelocity) const {
    if(disttype==0) {
        double sum = getMinimumDistance(in);//midTraj
        //frontTraj
        STPoint sPoint(in.m_pLow, 0, 2), sPoint2(in.m_pLow, 0, 2), ePoint(in.m_pHigh, 0, 2), ePoint2(in.m_pHigh, 0, 2);
        if (in.m_startTime > m_startTime()) {
            fakeTpVector frontTraj(&m_points, m_startTime(), in.m_startTime);
            for (int i = 0; i < frontTraj.m_size - 1; i++) {
                sPoint.m_time = frontTraj[i].m_time;
                sPoint2.m_time = frontTraj[i + 1].m_time;
                double pd = line2lineIED(sPoint, sPoint2, frontTraj[i], frontTraj[i + 1]);
                double r1 = (in.m_startTime-frontTraj[i].m_time) * MaxVelocity;
                double r2 = (in.m_startTime-frontTraj[i + 1].m_time) * MaxVelocity;
                double minus = 0.5 * (frontTraj[i + 1].m_time - frontTraj[i].m_time) * (r1 + r2);
                sum += std::max(0.0, pd - minus);
            }
        }
        if (in.m_endTime < m_endTime()) {
            fakeTpVector backTraj(&m_points, in.m_endTime, m_endTime());
            for (int i = 0; i < backTraj.m_size - 1; i++) {
                ePoint.m_time = backTraj[i].m_time;
                ePoint2.m_time = backTraj[i + 1].m_time;
                double pd = line2lineIED(ePoint, ePoint2, backTraj[i], backTraj[i + 1]);
                double r1 = (  backTraj[i].m_time-in.m_endTime) * MaxVelocity;
                double r2 = ( backTraj[i + 1].m_time-in.m_endTime) * MaxVelocity;
                double minus = 0.5 * (backTraj[i + 1].m_time - backTraj[i].m_time) * (r1 + r2);
                sum += std::max(0.0, pd - minus);
            }
        }
        return sum;
    }
    else if(disttype==1) {
        double min=getMaxSED(in);
        return min;
    }
    else
        throw Tools::NotSupportedException("Wrong distance");
}
void Trajectory::makeInfinite(uint32_t dimension)
{
    m_points.clear();
}


void Trajectory::combineTrajectory(const Trajectory& r)
{
    throw Tools::NotSupportedException(
            "Trajectory::...: Not implemented yet!"
    );
}

bool Trajectory::containsTrajectory(const SpatialIndex::Trajectory &r) {
    throw Tools::NotSupportedException(
            "Trajectory::...: Not implemented yet!"
    );
}

void Trajectory::getCombinedTrajectory(Trajectory& out, const Trajectory& in) const
{
    out = *this;
    out.combineTrajectory(in);
}
std::ostream& SpatialIndex::operator<<(std::ostream& os, const Trajectory& r) {
    std::string s;
    s = "Trajectory length:" + std::to_string(r.m_points.size()) + "\n" +
        "m_points are" + "\n";
    for (const auto &p:r.m_points) {
        s += std::to_string(p.m_pCoords[0]) + "," + std::to_string(p.m_pCoords[1]) +
             "," + std::to_string(p.m_time) + " ";
    }
    s += "\n";
    os<<s;

//    os<<"Trajectory length:"<<r.m_points.size()<<"\n"<<"m_points are"<< "\n"<<r.m_points.front()<<"\t"<<r.m_points.back()<<endl;
    return os;
}

std::string Trajectory::toString() const{
    std::string s;
    for (const auto &p:m_points) {
        s += std::to_string(p.m_pCoords[0]) + "," + std::to_string(p.m_pCoords[1]) +
             "," + std::to_string(p.m_time) + " ";
    }
    s.pop_back();
    return s;
}

void Trajectory::loadFromString(std::string str) {
    m_points.clear();
    std::vector<std::string> points=split(str,' ');
    for(const auto &p: points){
        std::vector<std::string> xyt=split(p,',');
        double xy[2];
        xy[0]=std::stod(xyt[0]);xy[1]=std::stod(xyt[1]);
        m_points.emplace_back(STPoint(xy,std::stod(xyt[2]),2));
    }
    if(m_points.size()>1 && m_points.front().m_time>=m_points.back().m_time){
        m_points.clear();
    }
}

void Trajectory::linkTrajectory(SpatialIndex::Trajectory &other) {
    if(m_points.back().m_time==other.m_points.front().m_time){
        if(m_fakeback) {
            m_points.pop_back();
            m_points.insert(m_points.end(),++other.m_points.begin(),other.m_points.end());
            m_fakeback=other.m_fakeback;
        }
        else {
            m_points.insert(m_points.end(), ++other.m_points.begin(), other.m_points.end());
        }
    }
    else if (other.m_points.back().m_time==m_points.front().m_time){
        if(m_fakehead) {
            m_points.erase(m_points.begin());
            m_points.insert(m_points.begin(), other.m_points.begin(), other.m_points.end() - 1);
            m_fakehead=other.m_fakehead;
        } else {
            m_points.insert(m_points.begin(), other.m_points.begin(), other.m_points.end() - 1);
        }
    }
    else{
//        std::cerr<<*this<<other<<"\n";
        std::cerr<<m_points.back()<<" "<<other.m_points.front()<<"\n"
            <<other.m_points.back()<<" "<<m_points.front();
        throw Tools::IllegalStateException("Trajectory::linkTrajectory: the two trajectories to be linked should have a common point.");
    }
}

Point Trajectory::pos(double t) const{
    if(t<=m_startTime()) return m_points[0];
    if(t>=m_endTime()) return m_points[m_points.size()-1];
    int is=0;
    while(m_points[is].m_time<t) is++;
    double x,y;
    if(m_points[is-1].m_time==t) return m_points[is=1];
    x=makemidmacro(m_points[is-1].m_pCoords[0],m_points[is-1].m_time,
                   m_points[is].m_pCoords[0],m_points[is].m_time,t);
    y=makemidmacro(m_points[is-1].m_pCoords[1],m_points[is-1].m_time,
                   m_points[is].m_pCoords[1],m_points[is].m_time,t);
    return STPoint(x,y,t);
}

void Trajectory::getPartialTrajectory(double tstart, double tend, SpatialIndex::Trajectory &out) const {
    //may produce non-exist points through makemid
    //get the inner or equal points
    out.makeInfinite(2);
    if(tstart==tend) return;
    if(tstart>m_endTime()||tend<m_startTime()) return;
    int is=0,ie=m_points.size()-1;
    while(m_points[is].m_time<tstart) is++;
    while(m_points[ie].m_time>tend) ie--;
    double x,y;
    if(is!=0&&m_points[is].m_time!=tstart){
        x=makemidmacro(m_points[is-1].m_pCoords[0],m_points[is-1].m_time,
                       m_points[is].m_pCoords[0],m_points[is].m_time,tstart);
        y=makemidmacro(m_points[is-1].m_pCoords[1],m_points[is-1].m_time,
                       m_points[is].m_pCoords[1],m_points[is].m_time,tstart);
        double xy[2]={x,y};
        out.m_points.push_back(STPoint(xy,tstart,2));
    }
    for(int i=is;i<=ie;i++){
        out.m_points.push_back(STPoint(m_points[i]));
    }
    if(ie!=m_points.size()-1&&m_points[ie].m_time!=tend){
        x=makemidmacro(m_points[ie].m_pCoords[0],m_points[ie].m_time,
                       m_points[ie+1].m_pCoords[0],m_points[ie+1].m_time,tstart);
        y=makemidmacro(m_points[ie].m_pCoords[1],m_points[ie].m_time,
                       m_points[ie+1].m_pCoords[1],m_points[ie+1].m_time,tstart);
        double xy[2]={x,y};
        out.m_points.push_back(STPoint(xy,tend,2));
    }
}

int Trajectory::cutTrajsIntoFile(std::vector<std::pair<SpatialIndex::id_type, SpatialIndex::Trajectory>> &trajs,
                                 double segLen, int strat, std::string filename) {

    tjstat->bt = segLen;
    double totallen=0;
    int maxseg=0;
    int totalseg=0;
    std::ofstream file(filename,std::ios::out);
    for (const auto &traj:trajs) {
        totallen += traj.second.m_points.size();
        std::vector<Trajectory> seg;
        switch(strat){
            case 0:
                seg = traj.second.getHybridSegments(segLen);
                break;
            case 1:
                seg = traj.second.getFixedSegments(int(segLen)+1);
                break;
            default:
                break;
        }
        totalseg += seg.size();
        maxseg = std::max(int(seg.size()), maxseg);
        file<<"SubTraj\n"<<traj.first<<"\n";
        for(auto s:seg){
            file<<s.toString()<<"\n";
        }
        file<<"ESubTraj"<<endl;
    }
    file<<"END"<<endl;
    double avgSegLen=double(totallen)/totalseg;
    std::cerr<<"total sub-trajs num:"<<totalseg<<"\n";
    std::cerr<<"segments' average length is "<<totallen*1.0/totalseg<<"\n";
    file.flush();
    file.close();
    return maxseg;
}

double Trajectory::maxSpeed() const {
    double res = 0;
    for(int i=0;i<m_points.size()-1;i++){
        double v = m_points[i].getMinimumDistance(m_points[i+1]) /
                   (m_points[i+1].m_time - m_points[i].m_time);
        res= max(v,res);
    }
    return res;
}