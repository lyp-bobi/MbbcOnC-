//
// Created by chuang on 4/23/19.
//

#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

Trajectory::Trajectory() {
}
Trajectory::Trajectory(std::vector<SpatialIndex::TimePoint>& in) {
    m_points=in;
}

Trajectory::Trajectory(const SpatialIndex::Trajectory &in) {
    m_points=in.m_points;
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
uint32_t Trajectory::getByteArraySize() {
    return sizeof(unsigned long)+m_points[0].getByteArraySize()*m_points.size();
}

void Trajectory::loadFromByteArray(const uint8_t* ptr) {
    unsigned long size;
    memcpy(&size, ptr, sizeof(unsigned long));
    ptr += sizeof(unsigned long);
    std::vector<TimePoint> p(size);
    for(int i=0;i<size;i++){
        p[i].loadFromByteArray(ptr);
        if(i!=size-1){
            ptr+=p[i].getByteArraySize();
        }
    }
    m_points.clear();
    m_points=p;
}

void Trajectory::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    uint8_t* tmpb;
    uint32_t tmplen;
    unsigned long size=m_points.size();
    memcpy(ptr, &size, sizeof(unsigned long));
    ptr += sizeof(unsigned long);
    for(int i=0;i<size;i++){
        m_points[i].storeToByteArray(&tmpb,tmplen);
        memcpy(ptr, tmpb, tmplen);
        if(i!=size-1){
            ptr += tmplen;
        }
    }

    assert(len==(ptr - *data)+tmplen);
}


//
// IShape interface
//
bool Trajectory::intersectsShape(const SpatialIndex::IShape& s) const {
    const Trajectory* pTrajectory = dynamic_cast<const Trajectory*>(&s);
    if (pTrajectory != 0) return intersectsTrajectory(*pTrajectory);

    const Mbbc* pbc = dynamic_cast<const Mbbc*>(&s);
    if (pbc != 0) return intersectsMbbc(*pbc);

    const TimeRegion* ptr = dynamic_cast<const TimeRegion*>(&s);
    if (ptr != 0) return intersectsTimeRegion(*ptr);

    const MBRk* pmbrk= dynamic_cast<const MBRk*>(&s);
    if(pmbrk!=0) return intersectsMBRk(*pmbrk);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return intersectsRegion(*pr);

    throw Tools::IllegalStateException(
            "Trajectory::intersectsShape: Not implemented yet!"
    );
}

TimePoint Trajectory::getPointAtTime(const double time) const {
    if(time<m_points.front().m_startTime||time>m_points.back().m_endTime){
        throw Tools::IllegalArgumentException(
                "Trajectory::getPointAtTime: time"+std::to_string(time)+"is illegal."
                );}
    if(m_points.size()==1){
        return m_points[0];
    }
    auto pre =m_points.begin(),next=m_points.begin();
    while(next->m_startTime-pre->m_startTime<0.01) next++;
    while(next->m_startTime<time&&next!=m_points.end()){
        pre=next;
        while(next->m_startTime-pre->m_startTime<0.01) next++;
    }
    double h1= (time-pre->m_startTime)/(next->m_startTime-pre->m_startTime);
    double h2= (next->m_startTime-time)/(next->m_startTime-pre->m_startTime);
    double *coords= new double(m_dimension);
    for (int i = 0; i < m_dimension; ++i) {
        coords[i]=h2*pre->m_pCoords[i]+h1*next->m_pCoords[i];
    }
    return TimePoint(coords,time,time,m_dimension);
}

bool Trajectory::intersectsMbbc(const SpatialIndex::Mbbc &in) const {
    for(int i=0;i<m_points.size()-1;i++){
        if(in.intersectsTimePoint(m_points[i])){
            return true;
        }
    }
    return false;
}

bool Trajectory::intersectsMBRk(const SpatialIndex::MBRk &in) const {
    for(int i=0;i<m_points.size()-1;i++){
        if(in.intersectsTimePoint(m_points[i])){
            return true;
        }
    }
    return false;
}

bool Trajectory::intersectsTimeRegion(const SpatialIndex::TimeRegion &in) const {
    if(in.m_startTime==in.m_endTime){//time slice
        if(in.m_startTime<m_points.front().m_startTime||in.m_startTime>m_points.back().m_endTime){
            return false;
        }
        TimePoint tp=getPointAtTime(in.m_startTime);
        return tp.intersectsShape(in);
    }else{
        throw Tools::NotSupportedException("time interval range not supported");
    }
}
bool Trajectory::intersectsRegion(const Region& in) const{
    for(int i=0;i<m_points.size();i++){
        if(m_points[i].intersectsShape(in)){
            return true;
        }
    }
    return false;
}

bool Trajectory::intersectsTrajectory(const Trajectory& in) const{
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}

bool Trajectory::containsShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}
bool Trajectory::touchesShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}
void Trajectory::getCenter(Point& out) const{
    throw Tools::NotSupportedException("not supported now");
}
uint32_t Trajectory::getDimension() const{return 2;}
void Trajectory::getMBR(Region& out) const{
    out.makeInfinite(m_dimension);
    for(int i=0;i<m_points.size();i++){
        out.combinePoint(m_points[i]);
    }
}
//todo: should have implement a time divivsion class!
void Trajectory::getMbbc(Mbbc& out) const{
    out.makeInfinite(m_dimension);

    double startx=m_points.begin()->m_pCoords[0],starty=m_points.begin()->m_pCoords[1],startt=m_points.begin()->m_startTime;
    double endx=m_points.back().m_pCoords[0],endy=m_points.back().m_pCoords[1],endt=m_points.back().m_startTime;
    startt=getPeriodStart(startt);
    endt=getPeriodEnd(endt);
    double maxvxP=-std::numeric_limits<double>::max(),
        maxvxN=std::numeric_limits<double>::max(),
        maxvyP=-std::numeric_limits<double>::max(),
        maxvyN=std::numeric_limits<double>::max();
    double minx=startx,maxx=startx,miny=starty,maxy=starty;
    for(int i=0;i<m_points.size();i++){
        if(m_points[i].m_startTime-startt>0){
            double vx=(m_points[i].m_pCoords[0]-startx)/(m_points[i].m_startTime-startt);
            if(vx>maxvxP) maxvxP=vx;
            if(vx<maxvxN) maxvxN=vx;

            double vy=(m_points[i].m_pCoords[1]-starty)/(m_points[i].m_startTime-startt);
            if(vy>maxvyP) maxvyP=vy;
            if(vy<maxvyN) maxvyN=vy;
        }
        if(endt-m_points[i].m_startTime>0){
            double vx=(endx-m_points[i].m_pCoords[0])/(endt-m_points[i].m_startTime);
            if(vx>maxvxP) maxvxP=vx;
            if(vx<maxvxN) maxvxN=vx;
            double vy=(endy-m_points[i].m_pCoords[1])/(endt-m_points[i].m_startTime);
            if(vy>maxvyP) maxvyP=vy;
            if(vy<maxvyN) maxvyN=vy;
        }

        if(m_points[i].m_pCoords[0]<minx) minx=m_points[i].m_pCoords[0];
        if(m_points[i].m_pCoords[0]>maxx) maxx=m_points[i].m_pCoords[0];
        if(m_points[i].m_pCoords[1]<miny) miny=m_points[i].m_pCoords[1];
        if(m_points[i].m_pCoords[1]>maxy) maxy=m_points[i].m_pCoords[1];
    }
    double sLow[2]={startx,starty};
    double sHigh[2]={startx,starty};
    double eLow[2]={endx,endy};
    double eHigh[2]={endx,endy};
    double vLow[2]={maxvxN,maxvyN};
    double vHigh[2]={maxvxP,maxvyP};
    double wLow[2]={minx,miny};
    double wHigh[2]={maxx,maxy};
//    double nstartx,nendx,nstarty,nendy;
//    nstartx=startx-(endx-startx)/(endt-startt)*startt;
//    nstarty=starty-(endy-starty)/(endt-startt)*startt;
//    nendx=endx+(endx-startx)/(endt-startt)*(PeriodLen-endt);
//    nendy=endy+(endy-starty)/(endt-startt)*(PeriodLen-endt);
//    double sLow[2]={nstartx,nstarty};
//    double sHigh[2]={nstartx,nstarty};
//    double eLow[2]={nendx,nendy};
//    double eHigh[2]={nendx,nendy};
//    double vLow[2]={maxvxN,maxvyN};
//    double vHigh[2]={maxvxP,maxvyP};
//    double wLow[2]={minx,miny};
//    double wHigh[2]={maxx,maxy};
//    double stime=int(startt/PeriodLen)*PeriodLen;
    out= Mbbc(Region(sLow,sHigh,2),Region(eLow,eHigh,2),
                Region(vLow,vHigh,2),Region(wLow,wHigh,2),startt,endt);

}
void Trajectory::getMBRk(int k, SpatialIndex::MBRk &out) const {
    out.m_k=k;
    out.makeInfinite(m_dimension,k);
    int oldPhase=0;
    for(int i=0;i<m_points.size();i++){
        int newPhase=out.getPhase(m_points[i].m_startTime);
        if(i==m_points.size()-1) newPhase=k-1;
        out.m_mbrs[newPhase].combinePoint(m_points[i]);
        if(oldPhase!=newPhase){
            for(int j=oldPhase;j<=newPhase;j++){
                out.m_mbrs[j].combinePoint(m_points[i]);
                if(i!=0) out.m_mbrs[j].combinePoint(m_points[i-1]);
            }
            oldPhase=newPhase;
        }
    }
}

double Trajectory::getArea() const{ return 0;}
double Trajectory::getMinimumDistance(const IShape& s) const{
    const Trajectory* pTrajectory = dynamic_cast<const Trajectory*>(&s);
    if (pTrajectory != 0) return getMinimumDistance(*pTrajectory);

    const Mbbc* pbc = dynamic_cast<const Mbbc*>(&s);
    if (pbc != 0) return getMinimumDistance(*pbc);

    const MBRk* pmbrk= dynamic_cast<const MBRk*>(&s);
    if(pmbrk!=0) return getMinimumDistance(*pmbrk);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return getMinimumDistance(*pr);

    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}

double Trajectory::getMinimumDistance(const SpatialIndex::Region &in) const {
    double sum=0;
    sum+=m_points[0].getMinimumDistance(in)*(m_points[1].m_startTime-m_points[0].m_startTime);
    for(int i=1;i<m_points.size()-1;i++){
        sum+=m_points[i].getMinimumDistance(in)*(m_points[i+1].m_startTime-m_points[i-1].m_startTime);
    }
    sum+=m_points[m_points.size()-1].getMinimumDistance(in)*(m_points[m_points.size()-1].m_startTime-m_points[m_points.size()-2].m_startTime);
    sum/=2;
    return sum;
}

double Trajectory::getMinimumDistance(const SpatialIndex::MBRk &in) const {
    double sum=0;
    double dist=in.getMinimumDistance(m_points[0]);
    sum+=dist*(m_points[1].m_startTime-m_points[0].m_startTime);
    for(int i=1;i<m_points.size()-1;i++){
        dist=in.getMinimumDistance(m_points[i]);
        sum+=dist*(m_points[i+1].m_startTime-m_points[i-1].m_startTime);
    }
    dist=in.getMinimumDistance(m_points[m_points.size()-1]);
    sum+=dist*(m_points[m_points.size()-1].m_startTime-m_points[m_points.size()-2].m_startTime);
    sum=sum/2;
    return sum;
}

double Trajectory::getMinimumDistance(const SpatialIndex::Mbbc &in) const {
    double sum=0;
    double dist=in.getMinimumDistance(m_points[0]);
    sum+=dist*(m_points[1].m_startTime-m_points[0].m_startTime);
    for(int i=1;i<m_points.size()-1;i++){
        dist=in.getMinimumDistance(m_points[i]);
        sum+=dist*(m_points[i+1].m_startTime-m_points[i-1].m_startTime);
    }
    dist=in.getMinimumDistance(m_points[m_points.size()-1]);
    sum+=dist*(m_points[m_points.size()-1].m_startTime-m_points[m_points.size()-2].m_startTime);
    sum=sum/2;
    return sum;
}

double Trajectory::getMinimumDistance(const SpatialIndex::Trajectory &in) const {
    //NOTICE: this function is not symmetry!
    int cursor=0;//cursor for the other trajectory
    double time1,time2;
    time1=m_points[0].m_startTime;
    time2=in.m_points[0].m_startTime;
    double dist;
    double sum=0;
    while(cursor<in.m_points.size()-1&&time1>in.m_points[cursor+1].m_startTime){
        cursor++;
        time2=in.m_points[cursor].m_startTime;
    }
    if(m_points.size()==1){
        return m_points[0].getMinimumDistance(in.m_points[cursor]);
    }
    TimePoint tmp;
    if(time1<time2||cursor==in.m_points.size()-1){//if no corresponding data
        tmp=in.m_points[cursor];
    }else{
        tmp=TimePoint::makemid(in.m_points[cursor],in.m_points[cursor+1],time1);
    }
    dist=m_points[0].getMinimumDistance(tmp);
    sum+=dist*(m_points[1].m_startTime-m_points[0].m_startTime);
    for(int i=1;i<m_points.size()-1;i++){
        time1=m_points[i].m_startTime;
        while(cursor<in.m_points.size()-1&&time1>in.m_points[cursor+1].m_startTime){
            cursor++;
            time2=in.m_points[cursor].m_startTime;
        }
        if(time1<time2||cursor==in.m_points.size()-1){//if no corresponding data
            tmp=in.m_points[cursor];
        }else{
            tmp=TimePoint::makemid(in.m_points[cursor],in.m_points[cursor+1],time1);
        }
        dist=m_points[i].getMinimumDistance(tmp);
        sum+=dist*(m_points[i+1].m_startTime-m_points[i-1].m_startTime);
    }
    time1=m_points[m_points.size()-1].m_startTime;
    while(cursor<in.m_points.size()-1&&time1>in.m_points[cursor+1].m_startTime){
        cursor++;
        time2=in.m_points[cursor].m_startTime;
    }
    if(time1<time2||cursor==in.m_points.size()-1){//if no corresponding data
        dist=m_points[m_points.size()-1].getMinimumDistance(in.m_points[cursor]);
    }else{
        dist=m_points[m_points.size()-1].getMinimumDistance(TimePoint::makemid(in.m_points[cursor],in.m_points[cursor+1],time1));
    }
    sum+=dist*(m_points[m_points.size()-1].m_startTime-m_points[m_points.size()-2].m_startTime);
    sum=sum/2;
    if(_isnan(sum))
        std::cerr<<"returning nan in MinDist between Trajectories\n"<<this->toString()<<in.toString();
    return sum;
}


void Trajectory::makeInfinite(uint32_t dimension)
{
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}


void Trajectory::combineTrajectory(const Trajectory& r)
{
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}

bool Trajectory::containsTrajectory(const SpatialIndex::Trajectory &r) {
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}

void Trajectory::getCombinedTrajectory(Trajectory& out, const Trajectory& in) const
{
    out = *this;
    out.combineTrajectory(in);
}
const std::string Trajectory::toString() const{
    std::string s;
    s="Trajectory length:"+std::to_string(m_points.size())+"\n"+
            "m_points are"+"\n";
    for(auto p:m_points){
        s+=std::to_string(p.m_pCoords[0])+","+std::to_string(p.m_pCoords[1])+
                ","+std::to_string(p.m_startTime)+" ";
    }
    s+="\n";
    return s;
}
std::vector<std::string> split(std::string strtem,char a)
{
    std::vector<std::string> strvec;

    std::string::size_type pos1, pos2;
    pos2 = strtem.find(a);
    pos1 = 0;
    while (std::string::npos != pos2)
    {
        strvec.push_back(strtem.substr(pos1, pos2 - pos1));

        pos1 = pos2 + 1;
        pos2 = strtem.find(a, pos1);
    }
    strvec.push_back(strtem.substr(pos1));
    return strvec;
}

void Trajectory::loadFromString(std::string str) {
    std::vector<std::string> points=split(str,' ');
    for(auto p: points){
        std::vector<std::string> xyt=split(p,',');
        double xy[2];
        xy[0]=std::stod(xyt[0]);xy[1]=std::stod(xyt[1]);
        m_points.push_back(TimePoint(xy,std::stod(xyt[2]),std::stod(xyt[2]),2));
    }
}

