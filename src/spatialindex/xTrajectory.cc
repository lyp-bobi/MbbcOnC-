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
//0 for IED, 1 for MaxSED

using namespace SpatialIndex;
using std::vector;
using std::cout;
using std::endl;
using std::sqrt;

xTrajectory::xTrajectory() {

}
xTrajectory::xTrajectory(std::vector<SpatialIndex::xPoint>& in)
{
    m_points=in;
}

xTrajectory::xTrajectory(bool fakehead, bool fakeback,std::vector<SpatialIndex::xPoint> &in) {
#ifndef NDEBUG
    if(in.front().m_t>=in.back().m_t){
        std::cerr<<"wrong xTrajectory!\n"<<*this;
    }
#endif
    m_points=in;
    m_fakehead=fakehead;
    m_fakeback=fakeback;
}

xTrajectory::xTrajectory(const SpatialIndex::xTrajectory &in) {
    m_points=in.m_points;
    m_fakehead=in.m_fakehead;
    m_fakeback=in.m_fakeback;
}

xTrajectory& xTrajectory::operator=(const xTrajectory& r)
{
    if(this != &r)
    {
        m_points=r.m_points;
    }

    return *this;
}

bool xTrajectory::operator==(const SpatialIndex::xTrajectory &r) const {
    if (m_points==r.m_points)
        return true;
    return false;
}
//
// IObject interface
//
xTrajectory* xTrajectory::clone() {
    return new xTrajectory(*this);
}
//
// ISerializable interface
//
uint32_t xTrajectory::getByteArraySize() const {
    return 2* sizeof(bool)+sizeof(unsigned long)+(m_dimension+1)* sizeof(double)*m_points.size();
}

void xTrajectory::loadFromByteArray(const uint8_t* ptr) {
    memcpy(&m_fakehead, ptr, sizeof(bool));
    ptr += sizeof(bool);
    memcpy(&m_fakeback, ptr, sizeof(bool));
    ptr += sizeof(bool);
    unsigned long size;
    memcpy(&size, ptr, sizeof(unsigned long));
    ptr += sizeof(unsigned long);
    m_points.clear();
    xPoint p;
    p.makeDimension(m_dimension);
    for(int i=0;i<size;i++){
        p.loadFromByteArray(ptr);
        m_points.emplace_back(p);
        if(i!=size-1){
            ptr += p.getByteArraySize();
        }
    }
}

void xTrajectory::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    uint32_t tmp;
    memcpy(ptr, &m_fakehead, sizeof(bool));
    ptr += sizeof(bool);
    memcpy(ptr, &m_fakeback, sizeof(bool));
    ptr += sizeof(bool);
    unsigned long size=m_points.size();
    memcpy(ptr, &size, sizeof(unsigned long));
    ptr += sizeof(unsigned long);
    for(int i=0;i<size;i++){
        m_points[i].storeToByteArray(&ptr,tmp);
        if(i!=size-1){
            ptr += tmp;
        }
    }
}


//
// IShape interface
//
bool xTrajectory::intersectsShape(const SpatialIndex::IShape& s) const {
    const xTrajectory* pxTrajectory = dynamic_cast<const xTrajectory*>(&s);
    if (pxTrajectory != 0) return intersectsxTrajectory(*pxTrajectory);

    const xMBR* pr = dynamic_cast<const xMBR*>(&s);
    if (pr != 0) return intersectsxMBR(*pr);

    const xCylinder* pcy = dynamic_cast<const xCylinder*>(&s);
    if (pcy != 0) return intersectsxCylinder(*pcy);

    throw Tools::IllegalStateException(
            "xTrajectory::intersectsShape: Not implemented yet!"
    );
}

xPoint xTrajectory::getPointAtTime(const double time) const {
    if(time<m_points.front().m_t||time>m_points.back().m_t){
        throw Tools::IllegalArgumentException(
                "xTrajectory::getPointAtTime: time"+std::to_string(time)+"is illegal."
                );}
    if(m_points.size()==1){
        return m_points[0];
    }
    auto pre =m_points.begin(),next=m_points.begin();
    while(next->m_t-pre->m_t<1e-7&&next!=m_points.end()) next++;
    while(next->m_t<time&&next!=m_points.end()){
        pre=next;
        while(next->m_t-pre->m_t<1e-7&&next!=m_points.end()) next++;
    }
    double h1= (time-pre->m_t)/(next->m_t-pre->m_t);
    double h2= (next->m_t-time)/(next->m_t-pre->m_t);
    xPoint ret;
    ret.m_x = h2 * pre->m_x + h1 * next->m_x;
    ret.m_y = h2 * pre->m_y + h1 * next->m_y;
    ret.m_t=time;
    return ret;
}

bool xTrajectory::intersectsxMBR(const xMBR& in) const{
        double ts=std::max(m_startTime(),in.m_tmin),
                te=std::min(m_endTime(),in.m_tmax);
        if(ts>te) return false;
        if(ts==te){
            return in.containsxPoint(getPointAtTime(ts));
        }
        fakeTpVector timedTraj(&m_points,ts,te);
        for(int i=0;i<timedTraj.m_size-1;i++){
            if(line2MBRMinSED(timedTraj[i],timedTraj[i+1],in)<1e-7){
                return true;
            }
        }
        return false;
}
bool xTrajectory::intersectsxCylinder(const SpatialIndex::xCylinder &in) const {
    double ts=std::max(m_startTime(),in.m_startTime),
            te=std::min(m_endTime(),in.m_endTime);
    if(ts>te) return false;
    if(ts==te){
        return xPoint(in.m_px,in.m_py,ts).getMinimumDistance(getPointAtTime(ts))<in.m_r;
    }
    fakeTpVector timedTraj(&m_points,ts,te);
    for(int i=0;i<timedTraj.m_size-1;i++){
        if(line2lineMinSED(timedTraj[i],timedTraj[i+1],
                           xPoint(in.m_px,in.m_py,timedTraj[i].m_t),
                           xPoint(in.m_px,in.m_py,timedTraj[i+1].m_t))<in.m_r){
            return true;
        }
    }
    return false;
}
bool xTrajectory::intersectsxTrajectory(const xTrajectory& in) const{
    throw Tools::NotSupportedException(
            "xTrajectory:::getMinimumDistance Not implemented yet!"
    );
}

bool xTrajectory::containsShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException(
            "xTrajectory::...: Not implemented yet!"
    );
}
bool xTrajectory::touchesShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException(
            "xTrajectory::...: Not implemented yet!"
    );
}
void xTrajectory::getCenter(Point& out) const{
    throw Tools::NotSupportedException("not supported now");
}
uint32_t xTrajectory::getDimension() const{return 2;}
void xTrajectory::getMBR(Region& out) const{
}


void xTrajectory::getxMBC(SpatialIndex::xMBC &out) const {
    if(m_points.size()<=1){
        std::cerr<<"WARNING: getting MBC at a xTrajectory with 0 or 1 points\n";
        out.makeInfinite(m_dimension+1);
        return;
    }
    if(std::fabs(m_endTime()-m_startTime())<1e-7){
        std::cerr<<"probably wrong segment at"<<*this<<"\n";
    }
    double startx=m_points[0].m_x,starty=m_points[0].m_y,startt=m_points[0].m_t;
    double endx=m_points.back().m_x,endy=m_points.back().m_y,endt=m_points.back().m_t;
    double avx=(endx-startx)/(endt-startt),avy=(endy-starty)/(endt-startt);

    xPoint p1=*m_points.begin(),p2=m_points.back();
    double rd=0,rv=0;
    double vx,vy;
    for(int i=0;i<m_points.size();i++){
        xPoint p;
        double ptime=m_points[i].m_t;
        double prd,prv;
        if(ptime-startt>0){
            vx=(m_points[i].m_x-startx)/(ptime-startt);
            vy=(m_points[i].m_y-starty)/(ptime-startt);
            prv=std::sqrt((vx-avx)*(vx-avx)+(vy-avy)*(vy-avy));
            if(prv>rv) rv=prv;
        }
        if(endt-m_points[i].m_t>0) {
            vx = (endx - m_points[i].m_x) / (endt - ptime);
            vy = (endy - m_points[i].m_y) / (endt - ptime);
            prv=std::sqrt((vx-avx)*(vx-avx)+(vy-avy)*(vy-avy));
            if(prv>rv) rv=prv;
        }
        prd=m_points[i].getMinimumDistance(xPoint::makemid(p1,p2,ptime));
        if(prd>rd) rd=prd;
    }
    out=xMBC(p1,p2,rd,rv);
}


const std::pair<int, double> findMaximumDistance(const vector<SpatialIndex::xPoint>& points) {
    SpatialIndex::xPoint firxPoint=points[0];
    SpatialIndex::xPoint laxPoint=points[points.size()-1];
    int index=0;  //index to be returned
    double Mdist=-1; //the Maximum distance to be returned

    //distance calculation
    for(int i=1;i<points.size()-1;i++){ //traverse through second point to second last point
        double Dist=SpatialIndex::xPoint::makemid(firxPoint,laxPoint,points[i].m_t).getMinimumDistance(points[i]);
        if (Dist>Mdist){
            Mdist=Dist;
            index=i;
        }
    }
    return std::make_pair(index, Mdist);
}
std::vector<SpatialIndex::xPoint> xTrajectory::simplifyWithRDP(const std::vector<SpatialIndex::xPoint> &Points,
                                                                 double epsilon){
    if(Points.size()<3){  //base case 1
        return Points;
    }
    std::pair<int, double> maxDistance=findMaximumDistance(Points);
    if(maxDistance.second>=epsilon){
        int index=maxDistance.first;
        vector<SpatialIndex::xPoint>::const_iterator it=Points.begin();
        vector<SpatialIndex::xPoint> path1(Points.begin(),it+index+1); //new path l1 from 0 to index
        vector<SpatialIndex::xPoint> path2(it+index,Points.end()); // new path l2 from index to last

        vector<SpatialIndex::xPoint> r1 =simplifyWithRDP(path1,epsilon);
        vector<SpatialIndex::xPoint> r2=simplifyWithRDP(path2,epsilon);

        //Concat simplified path1 and path2 together
        vector<SpatialIndex::xPoint> rs(r1);
        rs.pop_back();
        rs.insert(rs.end(),r2.begin(),r2.end());
        return rs;
    }
    else { //base case 2, all points between are to be removed.
        vector<SpatialIndex::xPoint> r(1,Points[0]);
        r.emplace_back(Points[Points.size()-1]);
        return r;
    }
}

std::vector<std::vector<SpatialIndex::xPoint>> xTrajectory::simplifyWithRDPN(const std::vector<SpatialIndex::xPoint> &Points,
                                                               int numPart){
    if(Points.size()<numPart+1){  //base case 1
        throw Tools::IllegalStateException("simplifyWithRDPN:Error");
    }
    std::pair<int, double> maxDistance=findMaximumDistance(Points);
    vector<vector<SpatialIndex::xPoint>> paths;
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
        vector<SpatialIndex::xPoint> *p=&paths[pathIndex];
        vector<SpatialIndex::xPoint>::const_iterator it1=p->begin();
        vector<SpatialIndex::xPoint>::const_iterator it2=p->end();
        auto indexIt=paths.begin()+pathIndex;
        vector<SpatialIndex::xPoint> path1(it1,it1+placeIndex+1); //new path l1 from 0 to index
        vector<SpatialIndex::xPoint> path2(it1+placeIndex,it2); // new path l2 from index to last
        paths.insert(indexIt+1,path2);
        indexIt=paths.begin()+pathIndex;
        paths.insert(indexIt+1,path1);
        indexIt=paths.begin()+pathIndex;
        paths.erase(indexIt);
    }
    return paths;
}

std::vector<xTrajectory> xTrajectory::cuttraj(std::vector<SpatialIndex::xPoint> mask){
    vector<xPoint> seg;
    vector<xTrajectory> res;
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
        while(iter1!=m_points.end()&&iter1->m_t<=iter2->m_t){
//            std::cerr<<"placed"<<*iter1<<"\n";
            seg.emplace_back(*iter1);
            iter1++;
        }
        if(seg.size()>=2) {
            if(!res.empty()){
                if(seg.front().m_t!=res.back().m_endTime()){
                    std::cerr<<"heelo\n";
                    for(const auto &t:mask) std::cerr<<t<<"\n";
                    std::cerr<<*this<<"\n";
                }
            }
            res.emplace_back(xTrajectory(seg));
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

std::vector<xTrajectory> xTrajectory::getRDPSegments(double len) const {
    int segNum=std::ceil((m_endTime()-m_startTime())/len);
    std::vector<xTrajectory> res;
    auto m=simplifyWithRDPN(m_points,std::min(int(m_points.size()-1),segNum));
    for(auto &pts:m){
        if(pts.size()<2){
            std::cerr<<"error on getting segments with "<<len<<"and \n"<<*this;
            throw Tools::IllegalStateException("bad RDPN");
        }
        res.emplace_back(xTrajectory(pts));
    }
    return res;
}

std::vector<xTrajectory> xTrajectory::getSegments(double len) const {
    auto stat=trajStat::instance();
    stat->bt=len;
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

std::vector<xTrajectory> xTrajectory::getHybridSegments(double len) const {
    int segNum=std::ceil((m_endTime()-m_startTime())/len);
    std::vector<xTrajectory> res;
    if(segNum==1) {res.emplace_back(*this);return res;}
//    std::cerr<<int(std::floor(std::pow(segNum,1.0/4)))<<"\t";
    auto m=simplifyWithRDPN(m_points,std::min(int(std::sqrt(m_points.size()-1)),int(segNum/5)));
    for(auto &pts:m){
        if(pts.size()<2){
            std::cerr<<"error on getting segments with "<<len<<"and \n"<<*this;
            throw Tools::IllegalStateException("bad RDPN");
        }
        auto seg=xTrajectory(pts).getGlobalSegmentsCut(len);
        res.insert(res.end(),seg.begin(),seg.end());
    }
    return res;
}

std::vector<xTrajectory> xTrajectory::getStaticSegmentsCut(double len) const{
    vector<xPoint> seg;
    vector<xTrajectory> res;
    bool fakehead=false,fakeback=false;
    if(m_points.size()<2) {
        throw Tools::IllegalStateException("getStatic:seg with 0 or 1 point");
    }
    double segStart=m_startTime();
    seg.emplace_back(m_points[0]);
    for(int i=1;i<m_points.size();i++){
        if(m_points[i].m_t<segStart+len){
            seg.emplace_back(m_points[i]);
        }
        else if(m_points[i].m_t==segStart+len){
            fakeback=false;
            seg.emplace_back(m_points[i]);
            res.emplace_back(xTrajectory(fakehead,fakeback,seg));
            fakehead=false;
            seg.clear();
            seg.emplace_back(m_points[i]);
            segStart=m_points[i].m_t;
        }
        else{
            xPoint mid=xPoint::makemid(m_points[i-1],m_points[i],segStart+len);
            fakeback=true;
            seg.emplace_back(mid);
            res.emplace_back(xTrajectory(fakehead,fakeback,seg));
            fakehead=true;
            seg.clear();
            seg.emplace_back(mid);
            segStart=mid.m_t;
            i--;
        }
    }
    if(seg.size()>1){
        res.emplace_back(xTrajectory(seg));
        seg.clear();
    }
    return res;
}

std::vector<xTrajectory> xTrajectory::getGlobalSegmentsCut(double len) const {
    vector<xPoint> seg;
    vector<xTrajectory> res;
    bool fakehead=false,fakeback=false;
    if(m_points.size()<2) {
        throw Tools::IllegalStateException("getStatic:seg with 0 or 1 point");
    }
    auto stat=trajStat::instance();
    double segStart=double(int(stat->mint));
    seg.emplace_back(m_points[0]);
    while(segStart+len<m_points[0].m_t+1e-7) segStart+=len;
    for(int i=1;i<m_points.size();i++){
        if(m_points[i].m_t<segStart+len){
            seg.emplace_back(m_points[i]);
        }
        else if(fabs(m_points[i].m_t-segStart-len)<1e-7){
            fakeback=false;
            seg.emplace_back(m_points[i]);
            res.emplace_back(xTrajectory(fakehead,fakeback,seg));
            fakehead=false;
            seg.clear();
            seg.emplace_back(m_points[i]);
            segStart+=len;
        }
        else{
            xPoint mid=xPoint::makemid(m_points[i-1],m_points[i],segStart+len);
            fakeback=true;
            seg.emplace_back(mid);
            res.emplace_back(xTrajectory(fakehead,fakeback,seg));
            fakehead=true;
            seg.clear();
            seg.emplace_back(mid);
            segStart+=len;
            i--;
        }
    }
    if(seg.size()>1){
        res.emplace_back(xTrajectory(fakehead,true,seg));
        seg.clear();
    }
    return res;
}


std::vector<xTrajectory> xTrajectory::getStaticSegments(double len) const{
    vector<xPoint> seg;
    vector<xTrajectory> res;
    if(m_points.size()<2) {
        throw Tools::IllegalStateException("getStatic:seg with 0 or 1 point");
    }
    double segStart=m_startTime();
    seg.emplace_back(m_points[0]);
    for(int i=1;i<m_points.size();i++){
        seg.emplace_back(m_points[i]);
        if(i==m_points.size()-1){
            res.emplace_back(xTrajectory(seg));
            seg.clear();
            break;
        }
        if(m_points[i+1].m_t-segStart>=len){
            res.emplace_back(xTrajectory(seg));
            seg.clear();
            seg.emplace_back(m_points[i]);
            segStart=m_points[i].m_t;
        }
    }
    return res;
}



std::vector<xTrajectory> xTrajectory::getFixedSegments(int len) const {
    if(len<2) len=2;
    vector<xPoint> seg;
    vector<xTrajectory> res;
    if(m_points.size()<2) {
        throw Tools::IllegalStateException("getFixed:seg with 0 or 1 point");
    }
    for(int i=0;i<m_points.size();i++){
        seg.emplace_back(m_points[i]);
        if(i==m_points.size()-1){
            res.emplace_back(xTrajectory(seg));
            seg.clear();
            break;
        }
        if(seg.size()>=len){
            res.emplace_back(xTrajectory(seg));
            seg.clear();
            seg.emplace_back(m_points[i]);
        }
    }
    return res;
}


std::vector<xTrajectory> xTrajectory::getItself() const {
    vector<xTrajectory> res;
    res.emplace_back(*this);
    return res;
}



double xTrajectory::getArea() const{ return 0;}

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

inline double xTrajectory::line2lineIED(const SpatialIndex::xPoint &p1s, const SpatialIndex::xPoint &p1e,
                                const SpatialIndex::xPoint &p2s, const SpatialIndex::xPoint &p2e) {
    if(p1s.m_t!=p2s.m_t|p1e.m_t!=p2e.m_t)
        throw Tools::IllegalStateException("line2lineIED: time period not the same");
    double ts = p1s.m_t, te = p1e.m_t;
    double dxs=p1s.m_x-p2s.m_x;
    double dys=p1s.m_y-p2s.m_y;
    double dxe=p1e.m_x-p2e.m_x;
    double dye=p1e.m_y-p2e.m_y;
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

inline double xTrajectory::line2lineMinSED(const SpatialIndex::xPoint &p1s, const SpatialIndex::xPoint &p1e,
                                const SpatialIndex::xPoint &p2s, const SpatialIndex::xPoint &p2e) {
    if(p1s.m_t!=p2s.m_t|p1e.m_t!=p2e.m_t)
        throw Tools::IllegalStateException("line2lineMinSED: time period not the same");
    double ts = p1s.m_t, te = p1e.m_t;
    double dxs=p1s.m_x-p2s.m_x;
    double dys=p1s.m_y-p2s.m_y;
    double dxe=p1e.m_x-p2e.m_x;
    double dye=p1e.m_y-p2e.m_y;
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

inline int xTrajectory::getPhase(const SpatialIndex::xMBR &r,const xPoint &p1,const xPoint &p2){
    // 7 8 9
    // 4 5 6
    // 1 2 3
    double xd1=r.m_xmin,xd2=r.m_xmax,yd1=r.m_ymin,yd2=r.m_ymax;
    int res=0;
    if(p1.m_x<=xd2+1e-7&&p2.m_x<=xd2+1e-7&&p1.m_x>=xd1-1e-7&&p2.m_x>=xd1-1e-7)
        res+=2;
    else if(p1.m_x<=xd1+1e-7&&p2.m_x<=xd1+1e-7)
        res+=1;
    else if(p1.m_x>=xd2-1e-7&&p2.m_x>=xd2-1e-7)
        res+=3;
    else return -1;
    if(p1.m_y<=yd2+1e-7&&p2.m_y<=yd2+1e-7&&p1.m_y>=yd1-1e-7&&p2.m_y>=yd1-1e-7)
        res+=3;
    else if(p1.m_y<=yd1+1e-7&&p2.m_y<=yd1+1e-7)
        res+=0;
    else if(p1.m_y>=yd2-1e-7&&p2.m_y>=yd2-1e-7)
        res+=6;
    else return -1;
    return res;
}

xPoint* cutByLine(const SpatialIndex::xPoint &ps, const SpatialIndex::xPoint &pe,double value,int axis){
    int otheraxis=1-axis;
    double axisv1=ps.m_pCoords(axis),axisv2=pe.m_pCoords(axis);
    if((axisv1<value)==(axisv2<value)||std::fabs(axisv1-value)<1e-7||std::fabs(axisv2-value)<1e-7) //need no cut
        return nullptr;
    else {
        double d = std::fabs(ps.m_pCoords(axis) - pe.m_pCoords(axis));
        double d1 = std::fabs(ps.m_pCoords(axis) - value) / d;
        double d2 = std::fabs(pe.m_pCoords(axis) - value) / d;
        //get p=d2*ps+d1*pe
        double xyt[3];
        if(axis==0) {
            xyt[0]=value;
            xyt[1]=d2 * ps.m_pCoords(otheraxis) + d1 * pe.m_pCoords(otheraxis);
            xyt[2]=d2 * ps.m_t + d1 * pe.m_t;
        }else{
            xyt[0]=d2 * ps.m_pCoords(otheraxis) + d1 * pe.m_pCoords(otheraxis);
            xyt[1]=value;
            xyt[2]=d2 * ps.m_t + d1 * pe.m_t;
        }
//        Tools::SmartPointer<xPoint> sp(new xPoint(xyt, xyt[2], 2));
        auto res=new xPoint(xyt[0],xyt[1],xyt[2]);
        return res;
//        return sp.get();
    }
}
std::vector<std::pair<xPoint,xPoint>> xTrajectory::cutByPhase(const SpatialIndex::xPoint &ps, const SpatialIndex::xPoint &pe,
                                                     const SpatialIndex::xMBR &r){
    double xd1=r.m_xmin,xd2=r.m_xmax,yd1=r.m_ymin,yd2=r.m_ymax;
    std::vector<std::pair<xPoint,xPoint>> res;
    std::vector<std::pair<xPoint,xPoint>> tmp;
    res.emplace_back(std::make_pair(ps,pe));
    for(const auto &line:res){
        xPoint* stp=cutByLine(line.first,line.second,xd1,0);
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
        xPoint *stp=cutByLine(line.first,line.second,xd2,0);
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
        xPoint *stp=cutByLine(line.first,line.second,yd1,1);
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
        xPoint *stp=cutByLine(line.first,line.second,yd2,1);
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

inline double xTrajectory::line2MBRMinSED_impl(const SpatialIndex::xPoint &ps, const SpatialIndex::xPoint &pe,
                                       const SpatialIndex::xMBR &r, int sr) {
    double ts = ps.m_t, te = pe.m_t;
//    if(std::fabs(te-ts)<1e-7) return 0;
    double res;
    if (sr == 5) return 0;
    else if (sr % 2 == 0){
        return std::min(ps.getMinimumDistance(r),pe.getMinimumDistance(r));
    }
    else {
        double px, py;
        if (sr == 1 || sr == 7) px = r.m_xmin;
        else px = r.m_xmax;
        if (sr == 1 || sr == 3) py = r.m_ymin;
        else py = r.m_ymax;
        return line2lineMinSED(ps, pe, xPoint(px,py,ts), xPoint(px,py,te));
    }
}

inline double xTrajectory::line2MBRMinSED(const SpatialIndex::xPoint &ps, const SpatialIndex::xPoint &pe,
                                  const SpatialIndex::xMBR &r) {
    assert(r.m_tmin<=ps.m_t);
    assert(r.m_tmax>=pe.m_t);
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

inline double xTrajectory::line2MBRIED_impl(const SpatialIndex::xPoint &ps, const SpatialIndex::xPoint &pe,
                                    const SpatialIndex::xMBR &r, int sr) {
    double ts = ps.m_t, te = pe.m_t;
    if(std::fabs(te-ts)<1e-7) return 0;
    double res;
    if (sr == 5) return 0;
    else if (sr % 2 == 0){
        return 0.5 * (ps.getMinimumDistance(r) + pe.getMinimumDistance(r)) * (te - ts);
    }
    else {
        double px, py;
        if (sr == 1 || sr == 7) px = r.m_xmin;
        else px = r.m_xmax;
        if (sr == 1 || sr == 3) py = r.m_ymin;
        else py = r.m_ymax;
        return line2lineIED(ps, pe, xPoint(px,py,ts), xPoint(px,py,te));
    }
}

inline double xTrajectory::line2MBRMaxSED(const SpatialIndex::xPoint &ps, const SpatialIndex::xPoint &pe,
                                  const SpatialIndex::xMBR &r) {
    assert(r.m_tmin<=ps.m_t);
    assert(r.m_tmax>=pe.m_t);
    return std::max(ps.getMinimumDistance(r),pe.getMinimumDistance(r));
}

inline double xTrajectory::line2MBRDistance(const SpatialIndex::xPoint &ps, const SpatialIndex::xPoint &pe,
                                    const SpatialIndex::xMBR &r) {
    //the line's time period should be in the MBR's time period
    assert(r.m_tmin<=ps.m_t);
    assert(r.m_tmax>=pe.m_t);
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

inline double xTrajectory::line2MBCDistance(const SpatialIndex::xPoint &ps, const SpatialIndex::xPoint &pe,
                                    const SpatialIndex::xMBC &r) {
    if(disttype==0) {
        double x1 = r.m_ps.m_x, y1 = r.m_ps.m_y, t1 = r.m_ps.m_t;
        double x2 = r.m_pe.m_x, y2 = r.m_pe.m_y, t2 = r.m_pe.m_t;
        double xts = makemidmacro(x1, t1, x2, t2, ps.m_t), xte = makemidmacro(x1, t1, x2, t2, pe.m_t);
        double yts = makemidmacro(y1, t1, y2, t2, ps.m_t), yte = makemidmacro(y1, t1, y2, t2, pe.m_t);
        double mbcr1 = std::min(std::min(r.m_rd, (ps.m_t - r.m_ps.m_t) * r.m_rv),
                                (r.m_pe.m_t - ps.m_t) * r.m_rv);
        double mbcr2 = std::min(std::min(r.m_rd, (pe.m_t - r.m_ps.m_t) * r.m_rv),
                                (r.m_pe.m_t - pe.m_t) * r.m_rv);
        xPoint p2s(xts,yts, ps.m_t), p2e(xte,yte, pe.m_t);
        double s = line2lineIED(ps, pe, p2s, p2e);
        double sm = mbcArea(r.m_ps.m_t, r.m_pe.m_t, r.m_rd, r.m_rv, ps.m_t, pe.m_t);
        if (ps.getMinimumDistance(p2s) > mbcr1 && pe.getMinimumDistance(p2e) > mbcr2) return s - sm;
        else return 0;
    }else if(disttype==1){
        return std::max(r.getMinimumDistance(ps),r.getMinimumDistance(pe));
    }
    else{throw Tools::NotSupportedException("Wrong distance");}
}

double xTrajectory::getMinimumDistance(const IShape& s) const{
    //using Integrative Squared Euclidean Distance
    const xTrajectory* pxTrajectory = dynamic_cast<const xTrajectory*>(&s);
    if (pxTrajectory != nullptr) return getMinimumDistance(*pxTrajectory);

//    const xPoint* ptp = dynamic_cast<const xPoint*>(&s);
//    if (ptp != 0) return getMinimumDistance(*ptp);

    const MBC* pmbc= dynamic_cast<const MBC*>(&s);
    if(pmbc!=nullptr) return getMinimumDistance(*pmbc);

    const xMBR* pr = dynamic_cast<const xMBR*>(&s);
    if (pr != nullptr) return getMinimumDistance(*pr);

    const ShapeList *psl = dynamic_cast<const ShapeList*>(&s);
    if (psl != nullptr) return getMinimumDistance(*psl);

    throw Tools::NotSupportedException(
            "xTrajectory::...: Not implemented yet!"
    );
}


double xTrajectory::getMinimumDistance(const SpatialIndex::xMBR &in) const {
        //note here we only calculate the partial distance
        //so just calculate intersection of time period
        //this function is used for inner nodes, not leaf nodes
        double tstart, tend;
        tstart = std::max(m_startTime(), in.m_tmin);
        tend = std::min(m_endTime(), in.m_tmax);
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

double xTrajectory::getMaxSED(const SpatialIndex::xMBR &in) const {
    //used for single xTrajectory
    double tstart, tend;
    tstart = std::max(m_startTime(), in.m_tmin);
    tend = std::min(m_endTime(), in.m_tmax);
    if (tstart >= tend) return 1e300;
    fakeTpVector timedTraj(&m_points, tstart, tend);
    double max = 0,pd;
    for (int i = 0; i < timedTraj.m_size; i++) {
        pd = timedTraj[i].getMinimumDistance(in);
        max = std::max(max, pd);
    }
    return max;
}

double xTrajectory::getMaxSED(const SpatialIndex::xMBC &in) const {
    //used for single xTrajectory
    double tstart, tend;
    tstart = std::max(m_startTime(), in.m_ps.m_t);
    tend = std::min(m_endTime(), in.m_pe.m_t);
    if (tstart >= tend) return 1e300;
    fakeTpVector timedTraj(&m_points, tstart, tend);
    double max = 0,pd;
    for (int i = 0; i < timedTraj.m_size; i++) {
        pd = in.getMinimumDistance(timedTraj[i]);
        max = std::max(max, pd);
    }
    xPoint tmpp;
    double tmpt=in.m_ps.m_t+in.m_rd/in.m_rv;
    if(tmpt>m_startTime()&&tmpt<m_endTime()){
        tmpp=getPointAtTime(tmpt);
        pd = in.getMinimumDistance(tmpp);
        max = std::max(max, pd);
    }
    tmpt=in.m_pe.m_t-in.m_rd/in.m_rv;
    if(tmpt>m_startTime()&&tmpt<m_endTime()){
        tmpp=getPointAtTime(tmpt);
        pd = in.getMinimumDistance(tmpp);
        max = std::max(max, pd);
    }
    return max;
}

double xTrajectory::getMinimumDistance(const SpatialIndex::xMBC &in) const {
    double tstart, tend;
    tstart = std::max(m_startTime(), in.m_ps.m_t);
    tend = std::min(m_endTime(), in.m_pe.m_t);
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

double xTrajectory::getMinimumDistance(const SpatialIndex::xPoint &in) const {
    int last;
    for(last=1;m_points[last].m_t<in.m_t;last++);
    return xPoint::makemid(m_points[last-1],m_points[last],in.m_t).getMinimumDistance(in);
}

double xTrajectory::getMinimumDistance(const SpatialIndex::xTrajectory &in) const {
    if(m_startTime()>=in.m_endTime()||m_endTime()<=in.m_startTime())
        return 1e300;
    fakeTpVector timedTraj2(&in.m_points,m_startTime(),m_endTime());
    double cut1=timedTraj2[0].m_t,cut2=timedTraj2[timedTraj2.m_size-1].m_t;
    double sum=0;
    double max=0;
    fakeTpVector midTraj(&m_points,cut1,cut2);
    if(m_startTime()<cut1){
        double pd;
        if(disttype==0) {
            pd= getStaticIED(timedTraj2[0].m_x, timedTraj2[0].m_y, m_startTime(), cut1);
            sum += pd;
        }else if(disttype==1){
            pd= getStaticMaxSED(timedTraj2[0].m_x, timedTraj2[0].m_y, m_startTime(), cut1);
            max=std::max(max,pd);
        }
    }
    if(m_endTime()>cut2){
        double pd;
        if(disttype==0) {
            pd= getStaticIED(timedTraj2[timedTraj2.m_size - 1].m_x,
                              timedTraj2[timedTraj2.m_size - 1].m_y, cut2, m_endTime());
            sum += pd;
        }else if(disttype==1){
            pd= getStaticMaxSED(timedTraj2[timedTraj2.m_size - 1].m_x,
                              timedTraj2[timedTraj2.m_size - 1].m_y, cut2, m_endTime());
            max=std::max(max,pd);
        }
    }
    if(midTraj.m_size!=0) {
        double newtime = midTraj[0].m_t, lasttime = midTraj[0].m_t;
        auto iter1 = midTraj.m_vectorPointer->begin()+midTraj.m_is;
        auto iter2 = timedTraj2.m_vectorPointer->begin()+timedTraj2.m_is;
        xPoint lastp1 = midTraj[0], lastp2 = timedTraj2[0], newp1, newp2;
        newp1.makeInfinite(2);newp2.makeInfinite(2);
        if(disttype==1) {max=std::max(max,iter1->getMinimumDistance(*iter2));}
        while (lasttime != timedTraj2[timedTraj2.m_size-1].m_t) {
            if ((iter1 + 1)->m_t == (iter2 + 1)->m_t) {
                newtime = (iter1 + 1)->m_t;
                newp1 = *(iter1 + 1);
                newp2 = *(iter2 + 1);
                iter1++;
                iter2++;
            } else if ((iter1 + 1)->m_t < (iter2 + 1)->m_t) {
                newtime = (iter1 + 1)->m_t;
                newp1 = *(iter1 + 1);
                double x=makemidmacro(iter2->m_x,iter2->m_t,(iter2 + 1)->m_x,(iter2 + 1)->m_t,newtime);
                double y=makemidmacro(iter2->m_y,iter2->m_t,(iter2 + 1)->m_y,(iter2 + 1)->m_t,newtime);
                newp2.m_x=x;
                newp2.m_y=y;
                newp2.m_t=newtime;
                iter1++;
            } else {
                newtime = (iter2 + 1)->m_t;
                double x=makemidmacro(iter1->m_x,iter1->m_t,(iter1 + 1)->m_x,(iter1 + 1)->m_t,newtime);
                double y=makemidmacro(iter1->m_y,iter1->m_t,(iter1 + 1)->m_y,(iter1 + 1)->m_t,newtime);
                newp1.m_x=x;
                newp1.m_y=y;
                newp1.m_t=newtime;
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

double xTrajectory::getMinimumDistanceInner(const SpatialIndex::xTrajectory &in) const {
    if(m_startTime()>=in.m_endTime()||m_endTime()<=in.m_startTime())
        return 1e300;
    fakeTpVector timedTraj2(&in.m_points,m_startTime(),m_endTime());
    double cut1=timedTraj2[0].m_t,cut2=timedTraj2[timedTraj2.m_size-1].m_t;
    double sum=0;
    double max=0;
    fakeTpVector midTraj(&m_points,cut1,cut2);
    if(midTraj.m_size!=0) {
        double newtime = midTraj[0].m_t, lasttime = midTraj[0].m_t;
        auto iter1 = midTraj.m_vectorPointer->begin()+midTraj.m_is;
        auto iter2 = timedTraj2.m_vectorPointer->begin()+timedTraj2.m_is;
        xPoint lastp1 = midTraj[0], lastp2 = timedTraj2[0], newp1, newp2;
        newp1.makeInfinite(2);newp2.makeInfinite(2);
        if(disttype==1) {max=std::max(max,iter1->getMinimumDistance(*iter2));}
        while (lasttime != timedTraj2[timedTraj2.m_size-1].m_t) {
            if ((iter1 + 1)->m_t == (iter2 + 1)->m_t) {
                newtime = (iter1 + 1)->m_t;
                newp1 = *(iter1 + 1);
                newp2 = *(iter2 + 1);
                iter1++;
                iter2++;
            } else if ((iter1 + 1)->m_t < (iter2 + 1)->m_t) {
                newtime = (iter1 + 1)->m_t;
                newp1 = *(iter1 + 1);
                double x=makemidmacro(iter2->m_x,iter2->m_t,(iter2 + 1)->m_x,(iter2 + 1)->m_t,newtime);
                double y=makemidmacro(iter2->m_y,iter2->m_t,(iter2 + 1)->m_y,(iter2 + 1)->m_t,newtime);
                newp2.m_x=x;
                newp2.m_y=y;
                newp2.m_t=newtime;
                iter1++;
            } else {
                newtime = (iter2 + 1)->m_t;
                double x=makemidmacro(iter1->m_x,iter1->m_t,(iter1 + 1)->m_x,(iter1 + 1)->m_t,newtime);
                double y=makemidmacro(iter1->m_y,iter1->m_t,(iter1 + 1)->m_y,(iter1 + 1)->m_t,newtime);
                newp1.m_x=x;
                newp1.m_y=y;
                newp1.m_t=newtime;
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

double xTrajectory::getFrontIED(SpatialIndex::xMBR smbr, double MaxVelocity) const {
    double sum=0;
    double ints=smbr.m_tmin;
    fakeTpVector frontTraj(&m_points,m_startTime(),ints);
    for(int i=0;i<frontTraj.m_size-1;i++) {
        smbr.m_tmin = frontTraj[i].m_t;
        smbr.m_tmax = frontTraj[i + 1].m_t;
        double pd = line2MBRDistance(frontTraj[i], frontTraj[i + 1], smbr);
        double r1 = (ints - frontTraj[i].m_t) * MaxVelocity;
        double r2 = (ints - frontTraj[i + 1].m_t) * MaxVelocity;
        double minus = 0.5 * (frontTraj[i + 1].m_t - frontTraj[i].m_t) * (r1 + r2);
        sum += std::max(0.0, pd - minus);
    }
    return sum;
}
double xTrajectory::getFrontIED(double x, double y, double ints, double MaxVelocity) const {
    fakeTpVector frontTraj(&m_points, m_startTime(), ints);
    xPoint sPoint(x,y,ints),sPoint2(x,y,ints);
    double sum=0;
    for (int i = 0; i < frontTraj.m_size - 1; i++) {
        sPoint.m_t = frontTraj[i].m_t;
        sPoint2.m_t = frontTraj[i + 1].m_t;
        double pd = line2lineIED(sPoint, sPoint2, frontTraj[i], frontTraj[i + 1]);
        double r1 = (ints-frontTraj[i].m_t) * MaxVelocity;
        double r2 = (ints-frontTraj[i + 1].m_t) * MaxVelocity;
        double minus = 0.5 * (frontTraj[i + 1].m_t - frontTraj[i].m_t) * (r1 + r2);
        sum += std::max(0.0, pd - minus);
    }
    return sum;
}
double xTrajectory::getBackIED(SpatialIndex::xMBR embr, double MaxVelocity) const {
    double sum=0;
    double inte=embr.m_tmax;
    fakeTpVector backTraj(&m_points,inte,m_endTime());
    for(int i=0;i<backTraj.m_size-1;i++){
        embr.m_tmin=backTraj[i].m_t;
        embr.m_tmax=backTraj[i+1].m_t;
        double pd=line2MBRDistance(backTraj[i],backTraj[i+1],embr);
        double r1=(backTraj[i].m_t-inte)*MaxVelocity;
        double r2=(backTraj[i+1].m_t-inte)*MaxVelocity;
        double minus = 0.5*(backTraj[i+1].m_t-backTraj[i].m_t)*(r1+r2);
        sum+=std::max(0.0,pd-minus);
    }
    return sum;
}
double xTrajectory::getBackIED(double x, double y,  double inte, double MaxVelocity) const {
    fakeTpVector backTraj(&m_points, inte, m_endTime());
    xPoint ePoint(x,y,inte),ePoint2(x,y,inte);
    double sum=0;
    for (int i = 0; i < backTraj.m_size - 1; i++) {
        ePoint.m_t = backTraj[i].m_t;
        ePoint2.m_t = backTraj[i + 1].m_t;
        double pd = line2lineIED(ePoint, ePoint2, backTraj[i], backTraj[i + 1]);
        double r1 = (  backTraj[i].m_t-inte) * MaxVelocity;
        double r2 = ( backTraj[i + 1].m_t-inte) * MaxVelocity;
        double minus = 0.5 * (backTraj[i + 1].m_t - backTraj[i].m_t) * (r1 + r2);
        sum += std::max(0.0, pd - minus);
    }
    return sum;
}
double ldd(double d,double v,double dt){
    if(d+v*dt>0) return dt*(d+v*dt/2);
    else return d*d/2*std::fabs(v);
}
double xTrajectory::getMidIED(const SpatialIndex::xMBR &sbr, const SpatialIndex::xMBR &ebr,
                             double MaxVelocity,double queryVelocity) {
    double sum=0;
    double ints=sbr.m_tmax,inte=ebr.m_tmin;
    fakeTpVector timedTraj(&m_points, ints, inte);
    for(int i=0;i<timedTraj.m_size-1;i++){
        double t1=timedTraj[i].m_t,t2=timedTraj[i+1].m_t;
        double l1 = (t1 - ints) * MaxVelocity,
                l2 = (t2 - ints) * MaxVelocity,
                h1 = (inte - t1) * MaxVelocity,
                h2 = (inte - t2) * MaxVelocity;
        double ds1=timedTraj[i].getMinimumDistance(sbr),
                ds2=timedTraj[i+1].getMinimumDistance(sbr),
                de1=timedTraj[i].getMinimumDistance(ebr),
                de2=timedTraj[i+1].getMinimumDistance(ebr);
        double a1=l1-ds1,a2=l2-ds2,b1=h1-de1,b2=h2-de2;
        double d1=std::min(a1,b1),d2=std::min(a2,b2);
        double minus,pd;
        if(a2>=b2){
            minus=0.5*(l1+l2)*(t2-t1);
            pd= getStaticIED(sbr, t1, t2);
            sum+=std::max(0.0,pd-minus);
        }else if(b1>=a1){
            minus=0.5*(h1+h2)*(t2-t1);
            pd= getStaticIED(ebr, t1, t2);
            sum+=std::max(0.0,pd-minus);
        }else{
            double vmax;
            if(queryVelocity==-1) vmax=2*MaxVelocity;
            else vmax=MaxVelocity+queryVelocity;
            double t0=(t1+t2+(d2-d1)/(vmax))/2;
            sum+=ldd(d1,-vmax,t0-t1)+ldd(d2,vmax,inte-t2);
        }
    }
    return sum;
}
double xTrajectory::getMidIED(const xMBC &sbc, const xMBC &ebc,
                             double MaxVelocity,double queryVelocity) {
    double sum=0;
    xPoint sPoint(sbc.m_pe),ePoint(ebc.m_ps);
    double ints=sPoint.m_t,inte=ePoint.m_t;
    fakeTpVector timedTraj(&m_points, ints, inte);
    for(int i=0;i<timedTraj.m_size-1;i++){
        double t1=timedTraj[i].m_t,t2=timedTraj[i+1].m_t;
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
            minus=0.5*(l1+l2)*(timedTraj[i + 1].m_t-timedTraj[i].m_t);
            pd= getStaticIED(sPoint.m_x,sPoint.m_y, t1, t2);
            sum+=std::max(0.0,pd-minus);
        }else if(b1>=a1){
            minus=0.5*(h1+h2)*(timedTraj[i + 1].m_t-timedTraj[i].m_t);
            pd= getStaticIED(ePoint.m_x,ePoint.m_y, t1, t2);
            sum+=std::max(0.0,pd-minus);
        }else{
            double vmax;
            if(queryVelocity==-1) vmax=2*MaxVelocity;
            else vmax=MaxVelocity+queryVelocity;
            double t0=(t1+t2+(d2-d1)/(vmax))/2;
            sum+=ldd(d1,-vmax,t0-t1)+ldd(d2,vmax,inte-t2);
//            double x1=timedTraj[i].m_x,x2=timedTraj[i+1].m_x,
//                    y1=timedTraj[i].m_y,y2=timedTraj[i+1].m_y;
//            double psx=sPoint.m_x,psy=sPoint.m_y,
//                pex=ePoint.m_x,pey=ePoint.m_y;
//            double t0=(2*(pex-psx)*(inte))
        }
    }
    return sum;
}

double xTrajectory::getMidIED(const xPoint &sPoint, const xPoint &ePoint,
                             double MaxVelocity,double queryVelocity) {
    double sum=0;
    double ints=sPoint.m_t,inte=ePoint.m_t;
    fakeTpVector timedTraj(&m_points, ints, inte);
    for(int i=0;i<timedTraj.m_size-1;i++){
        double t1=timedTraj[i].m_t,t2=timedTraj[i+1].m_t;
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
            minus=0.5*(l1+l2)*(timedTraj[i + 1].m_t-timedTraj[i].m_t);
            pd= getStaticIED(sPoint.m_x,sPoint.m_y, t1, t2);
            sum+=std::max(0.0,pd-minus);
        }else if(b1>=a1){
            minus=0.5*(h1+h2)*(timedTraj[i + 1].m_t-timedTraj[i].m_t);
            pd= getStaticIED(ePoint.m_x,ePoint.m_y, t1, t2);
            sum+=std::max(0.0,pd-minus);
        }else{
            double vmax;
            if(queryVelocity==-1) vmax=2*MaxVelocity;
            else vmax=MaxVelocity+queryVelocity;
            double t0=(t1+t2+(d2-d1)/(vmax))/2;
            sum+=ldd(d1,-vmax,t0-t1)+ldd(d2,vmax,inte-t2);
//            double x1=timedTraj[i].m_x,x2=timedTraj[i+1].m_x,
//                    y1=timedTraj[i].m_y,y2=timedTraj[i+1].m_y;
//            double psx=sPoint.m_x,psy=sPoint.m_y,
//                pex=ePoint.m_x,pey=ePoint.m_y;
//            double t0=(2*(pex-psx)*(inte))
        }
    }
    return sum;
}

double xTrajectory::getStaticIED(double x, double y, double t1, double t2) const {
    double tstart, tend;
    tstart = std::max(m_startTime(), t1);
    tend = std::min(m_endTime(), t2);
    if(tstart>=tend) return 1e300;
    fakeTpVector timedTraj(&m_points,tstart,tend);
    double sum = 0;
    xPoint ps(x,y,0),pe(x,y,0);
    for (int i = 0; i < timedTraj.m_size-1; i++) {
            ps.m_t=timedTraj[i].m_t;
            pe.m_t=timedTraj[i+1].m_t;
            double pd = line2lineIED(timedTraj[i], timedTraj[i + 1], ps, pe);
            sum += pd;
    }
    return sum;
}


double xTrajectory::getStaticIED(SpatialIndex::xMBR in,double ints, double inte) const {
    assert(ints>=m_startTime()&&inte<=m_endTime());
    in.m_tmin=ints;
    in.m_tmax=inte;
    fakeTpVector timedTraj(&m_points,ints,inte);
    double sum = 0;
    for (int i = 0; i < timedTraj.m_size-1; i++) {
        double pd = line2MBRDistance(timedTraj[i],timedTraj[i+1],in);
        sum+=pd;
    }
    return sum;
}
double xTrajectory::getStaticMaxSED(double x, double y, double t1, double t2) const {
    double tstart, tend;
    tstart = std::max(m_startTime(), t1);
    tend = std::min(m_endTime(), t2);
    if(tstart>=tend) return 1e300;
    double max=0;
    fakeTpVector timedTraj(&m_points,tstart,tend);
    xPoint p(x,y,0);
    for(int i=0;i<timedTraj.m_size;i++){
        double pd=timedTraj[i].getMinimumDistance(p);
        if(pd>max) max=pd;
    }
    return max;
}


double xTrajectory::getStaticMaxSED(SpatialIndex::xMBR in,double ints, double inte) const {
    assert(ints>=m_startTime()&&inte<=m_endTime());
    in.m_tmin=ints;
    in.m_tmax=inte;
    fakeTpVector timedTraj(&m_points,ints,inte);
    double max = 0;
    for (int i = 0; i < timedTraj.m_size; i++) {
        double pd = timedTraj[i].getMinimumDistance(in);
        if(pd>max) max=pd;
    }
    return max;
}

#define sign(x) (x>0?1.0:-1.0)
double xTrajectory::getInterTime(const xPoint &ps, const xPoint &pe,xMBR &r,double t_o, double vmax) const {
    //vmax could be positive or negative
    int sr = getPhase(r, ps, pe);
    if (sr > 0) {
        if(sr==5){
            std::cout<<",...";
        }
        assert(sr!=5);//should be handled previously
        if(sr%2==0){
            double d1=ps.getMinimumDistance(r)-(ps.m_t-t_o)*vmax;
            double d2=pe.getMinimumDistance(r)-(pe.m_t-t_o)*vmax;
            assert(sign(d1)*sign(d2)<=0);
            return (d1*pe.m_t+d2*ps.m_t)/(pe.m_t-ps.m_t);
        }
        else{
            double px, py;
            if (sr == 1 || sr == 7) px = r.m_xmin;
            else px = r.m_xmax;
            if (sr == 1 || sr == 3) py = r.m_ymin;
            else py = r.m_ymax;
            double ts = ps.m_t, te = pe.m_t;
            double dxs=ps.m_x-px;
            double dys=ps.m_y-py;
            double dxe=pe.m_x-px;
            double dye=pe.m_y-py;
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
            double d1=p.first.getMinimumDistance(r)-(p.first.m_t-t_o)*vmax;
            double d2=p.second.getMinimumDistance(r)-(p.second.m_t-t_o)*vmax;
            if(sign(d1)*sign(d2)<0) return getInterTime(p.first,p.second,r,t_o,vmax);
        }
        throw Tools::IllegalStateException("getInterTime");
    }
}

double xTrajectory::getInferredNodeMinimumIED(const SpatialIndex::xMBR &in, double MaxVelocity,double queryMaxVelocity) const {


    xMBR sembr=in;
    double vmax=std::max(MaxVelocity,queryMaxVelocity);
    double t0=in.m_tmin,t3=in.m_tmax;
    double tmid=(t0+t3)/2;
    double d1=m_points.front().getMinimumDistance(in)-vmax*(m_points.front().m_t-t0);
    double d2=m_points.back().getMinimumDistance(in)-vmax*(t3-m_points.back().m_t);
    if(d1>0&&d2>0){//if it has no chance to intersect
        double pd=getStaticIED(in,m_startTime(),m_endTime());
        double minus;
        if(m_endTime()<tmid){//only lower side
            minus=(m_endTime()-m_startTime())*((2*t3-m_points.front().m_t-m_points.back().m_t)*vmax)/2;
        }else if(m_startTime()>tmid){//only higher side
            minus=(m_endTime()-m_startTime())*((m_points.front().m_t+m_points.back().m_t-2*t0))*vmax/2;
        }else{//two sides
            minus=(m_endTime()-tmid)*(vmax*(t3-t0)/2+vmax*(m_points.back().m_t-t0)*vmax)/2+
                    (tmid-m_startTime())*(vmax*t3-t0+vmax*(t3-m_points.front().m_t)*vmax)/2;
        }
        return pd-minus;
    }
    if(m_startTime()>tmid&&d1<0) return 0;
    if(m_endTime()<tmid&&d2<0) return 0;
    double ts=std::max(t0,m_startTime()),te=std::min(t3,m_endTime());
    if(ts>=te) return 1e300; //the distance could be measured, but this means "some other part is better related"
    double inter1=m_startTime(),inter2=m_endTime();//the period that have positive distance
    for(int i=0;i<m_points.size()-1;i++) {
        double t1=m_points[i].m_t,t2=m_points[i+1].m_t;
        if(t1>inter2||t2<inter1) break;
        double l1 = (t1 - t0) * vmax,//smaller
                l2 = (t2 - t0) * vmax, //larger
                h1 = (t3 - t1) * vmax,//larger
                h2 = (t3 - t2) * vmax;//smaller
        double ds1=m_points[i].getMinimumDistance(in),
                ds2=m_points[i+1].getMinimumDistance(in);
        if(t1<tmid&&tmid<t2){
            double x,y;
            x=makemidmacro(m_points[i].m_x,m_points[i].m_t,
                           m_points[i+1].m_x,m_points[i+1].m_t,tmid);
            y=makemidmacro(m_points[i].m_y,m_points[i].m_t,
                           m_points[i+1].m_y,m_points[i+1].m_t,tmid);
            xPoint pt(x,y,tmid);
            double dm=pt.getMinimumDistance(in)-vmax*(te-ts)/2;
            if(dm<=0) return 0;
            else{
                if(ds1-h1<=0) inter1=getInterTime(m_points[i],pt,sembr,t3,-vmax);
                if(ds2-l2<=0) inter2=getInterTime(pt,m_points[i+1],sembr,t0,vmax);
            }
        }
        if(t2>tmid){// inter2 might be updated
            if(ds2-l2<=0&&ds1-l1>=0){
                inter2=getInterTime(m_points[i],m_points[i+1],sembr,t0,vmax);
            }
        }
        if(t1<tmid){// inter1 might be updated
            if(ds2-h2>=0&&ds1-h1<=0){
                inter1=getInterTime(m_points[i],m_points[i+1],sembr,t3,-vmax);
            }
        }
    }
    fakeTpVector timedTraj(&m_points,inter1,inter2);
    sembr.m_tmin=inter1;
    sembr.m_tmax=inter2;
    double sum = 0;
    for (int i = 0; i < timedTraj.m_size-1; i++) {
        double pd = line2MBRDistance(timedTraj[i],timedTraj[i+1],sembr);
        sum+=pd;
    }
    double minus;
    if(inter2<tmid){//only lower side
        minus=(inter2-inter1)*((2*t3-timedTraj.front().m_t-timedTraj.back().m_t)*vmax)/2;
    }else if(inter1>tmid){//only higher side
        minus=(inter2-inter1)*((timedTraj.front().m_t+timedTraj.back().m_t-2*t0))*vmax/2;
    }else{//two sides
        minus=(inter2-tmid)*(vmax*(t3-t0)/2+vmax*(timedTraj.back().m_t-t0)*vmax)/2+
              (tmid-inter1)*(vmax*t3-t0+vmax*(t3-timedTraj.front().m_t)*vmax)/2;
    }
    return sum-minus;
}

double xTrajectory::getNodeMinimumDistance(const SpatialIndex::xMBR &in,double MaxVelocity) const {
    if(m_startTime()>=in.m_tmax||m_endTime()<=in.m_tmin) return 1e300;
    if(disttype==0){
        double ints=std::max(m_startTime(),in.m_tmin),
                inte=std::min(m_endTime(),in.m_tmax);
        xMBR copy(in);
        copy.m_tmin=ints;
        copy.m_tmax=inte;
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

double xTrajectory::getLeafMinimumDistance(const SpatialIndex::xMBR &in, double MaxVelocity) const {
    if(disttype==0) {
        double sum=0;
        sum = getMinimumDistance(in);//midTraj
        //frontTraj
        xMBR sembr=in;
        if(in.m_tmin>m_startTime()){
            fakeTpVector frontTraj(&m_points,m_startTime(),in.m_tmin);
            for(int i=0;i<frontTraj.m_size-1;i++){
                sembr.m_tmin=frontTraj[i].m_t;
                sembr.m_tmax=frontTraj[i+1].m_t;
                double pd=line2MBRDistance(frontTraj[i],frontTraj[i+1],sembr);
                double r1=(in.m_tmin-frontTraj[i].m_t)*MaxVelocity;
                double r2=(in.m_tmin-frontTraj[i+1].m_t)*MaxVelocity;
                double minus = 0.5*(frontTraj[i+1].m_t-frontTraj[i].m_t)*(r1+r2);
                sum+=std::max(0.0,pd-minus);
            }
        }
        //backTraj
        if(in.m_tmax<m_endTime()){
            fakeTpVector backTraj(&m_points,in.m_tmax,m_endTime());
            for(int i=0;i<backTraj.m_size-1;i++){
                sembr.m_tmin=backTraj[i].m_t;
                sembr.m_tmax=backTraj[i+1].m_t;
                double pd=line2MBRDistance(backTraj[i],backTraj[i+1],sembr);
                double r1=(backTraj[i].m_t-in.m_tmax)*MaxVelocity;
                double r2=(backTraj[i+1].m_t-in.m_tmax)*MaxVelocity;
                double minus = 0.5*(backTraj[i+1].m_t-backTraj[i].m_t)*(r1+r2);
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

double xTrajectory::getLeafMinimumDistance(const SpatialIndex::xMBC &in, double MaxVelocity) const {
    if(disttype==0) {
        double sum = getMinimumDistance(in);//midTraj
        //frontTraj
        xPoint sPoint(in.m_ps), sPoint2(in.m_ps), ePoint(in.m_pe), ePoint2(in.m_pe);
        if (in.m_ps.m_t > m_startTime()) {
            fakeTpVector frontTraj(&m_points, m_startTime(), in.m_ps.m_t);
            for (int i = 0; i < frontTraj.m_size - 1; i++) {
                sPoint.m_t = frontTraj[i].m_t;
                sPoint2.m_t = frontTraj[i + 1].m_t;
                double pd = line2lineIED(sPoint, sPoint2, frontTraj[i], frontTraj[i + 1]);
                double r1 = (in.m_ps.m_t-frontTraj[i].m_t) * MaxVelocity;
                double r2 = (in.m_ps.m_t-frontTraj[i + 1].m_t) * MaxVelocity;
                double minus = 0.5 * (frontTraj[i + 1].m_t - frontTraj[i].m_t) * (r1 + r2);
                sum += std::max(0.0, pd - minus);
            }
        }
        if (in.m_pe.m_t < m_endTime()) {
            fakeTpVector backTraj(&m_points, in.m_pe.m_t, m_endTime());
            for (int i = 0; i < backTraj.m_size - 1; i++) {
                ePoint.m_t = backTraj[i].m_t;
                ePoint2.m_t = backTraj[i + 1].m_t;
                double pd = line2lineIED(ePoint, ePoint2, backTraj[i], backTraj[i + 1]);
                double r1 = (  backTraj[i].m_t-in.m_pe.m_t) * MaxVelocity;
                double r2 = ( backTraj[i + 1].m_t-in.m_pe.m_t) * MaxVelocity;
                double minus = 0.5 * (backTraj[i + 1].m_t - backTraj[i].m_t) * (r1 + r2);
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
void xTrajectory::makeInfinite(uint32_t dimension)
{
    m_points.clear();
}


void xTrajectory::combinexTrajectory(const xTrajectory& r)
{
    throw Tools::NotSupportedException(
            "xTrajectory::...: Not implemented yet!"
    );
}

bool xTrajectory::containsxTrajectory(const SpatialIndex::xTrajectory &r) {
    throw Tools::NotSupportedException(
            "xTrajectory::...: Not implemented yet!"
    );
}

void xTrajectory::getCombinedxTrajectory(xTrajectory& out, const xTrajectory& in) const
{
    out = *this;
    out.combinexTrajectory(in);
}
std::ostream& SpatialIndex::operator<<(std::ostream& os, const xTrajectory& r) {
    std::string s;
    s = "xTrajectory length:" + std::to_string(r.m_points.size()) + "\n" +
        "m_points are" + "\n";
    for (const auto &p:r.m_points) {
        s += std::to_string(p.m_x) + "," + std::to_string(p.m_y) +
             "," + std::to_string(p.m_t) + " ";
    }
    s += "\n";
    os<<s;

//    os<<"xTrajectory length:"<<r.m_points.size()<<"\n"<<"m_points are"<< "\n"<<r.m_points.front()<<"\t"<<r.m_points.back()<<endl;
    return os;
}



std::vector<std::string> split(const std::string &strtem,char a)
{
    std::vector<std::string> strvec;

    std::string::size_type pos1, pos2;
    pos2 = strtem.find(a);
    pos1 = 0;
    while (std::string::npos != pos2)
    {
        strvec.emplace_back(strtem.substr(pos1, pos2 - pos1));

        pos1 = pos2 + 1;
        pos2 = strtem.find(a, pos1);
    }
    strvec.emplace_back(strtem.substr(pos1));
    return strvec;
}

std::string xTrajectory::toString() const{
    std::string s;
    for (const auto &p:m_points) {
        s += std::to_string(p.m_x) + "," + std::to_string(p.m_y) +
             "," + std::to_string(p.m_t) + " ";
    }
    s.pop_back();
    return s;
}

void xTrajectory::loadFromString(std::string str) {
    m_points.clear();
    std::vector<std::string> points=split(str,' ');
    for(const auto &p: points){
        std::vector<std::string> xyt=split(p,',');
        m_points.emplace_back(xPoint(std::stod(xyt[0]),std::stod(xyt[1]),std::stod(xyt[2])));
    }
    if(m_points.size()>1 && m_points.front().m_t>=m_points.back().m_t){
        m_points.clear();
    }
}

void xTrajectory::linkxTrajectory(SpatialIndex::xTrajectory &other) {
    if(m_points.back().m_t==other.m_points.front().m_t){
        if(m_fakeback) {
            m_points.pop_back();
            m_points.insert(m_points.end(),++other.m_points.begin(),other.m_points.end());
            m_fakeback=other.m_fakeback;
        }
        else {
            m_points.insert(m_points.end(), ++other.m_points.begin(), other.m_points.end());
        }
    }
    else if (other.m_points.back().m_t==m_points.front().m_t){
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
        throw Tools::IllegalStateException("xTrajectory::linkxTrajectory: the two trajectories to be linked should have a common point.");
    }
}


void xTrajectory::getPartialxTrajectory(double tstart, double tend, SpatialIndex::xTrajectory &out) const {
    //may produce non-exist points through makemid
    //get the inner or equal points
    out.makeInfinite(2);
    if(tstart==tend) return;
    if(tstart>m_endTime()||tend<m_startTime()) return;
    int is=0,ie=m_points.size()-1;
    while(m_points[is].m_t<tstart) is++;
    while(m_points[ie].m_t>tend) ie--;
    double x,y;
    if(is!=0&&m_points[is].m_t!=tstart){
        x=makemidmacro(m_points[is-1].m_x,m_points[is-1].m_t,
                       m_points[is].m_x,m_points[is].m_t,tstart);
        y=makemidmacro(m_points[is-1].m_y,m_points[is-1].m_t,
                       m_points[is].m_y,m_points[is].m_t,tstart);
        out.m_points.push_back(xPoint(x,y,tstart));
    }
    for(int i=is;i<=ie;i++){
        out.m_points.push_back(xPoint(m_points[i]));
    }
    if(ie!=m_points.size()-1&&m_points[ie].m_t!=tend){
        x=makemidmacro(m_points[ie].m_x,m_points[ie].m_t,
                       m_points[ie+1].m_x,m_points[ie+1].m_t,tstart);
        y=makemidmacro(m_points[ie].m_y,m_points[ie].m_t,
                       m_points[ie+1].m_y,m_points[ie+1].m_t,tstart);
        out.m_points.push_back(xPoint(x,y,tend));
    }
}

int xTrajectory::cutTrajsIntoFile(std::vector<std::pair<SpatialIndex::id_type, SpatialIndex::xTrajectory>> &trajs,
                                 double segLen, int strat, std::string filename) {
    auto stat = trajStat::instance();
    stat->bt = segLen;
    double totallen=0;
    int maxseg=0;
    int totalseg=0;
    std::ofstream file(filename,std::ios::out);
    for (const auto &traj:trajs) {
        totallen += traj.second.m_points.size();
        std::vector<xTrajectory> seg;
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
