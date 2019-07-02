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
using std::vector;
using std::cout;
using std::endl;
using std::sqrt;

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
uint32_t Trajectory::getByteArraySize() const {
    if(m_points.size()<=0) throw Tools::IllegalStateException("traj with length 0!");
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

    const TimeRegion* ptr = dynamic_cast<const TimeRegion*>(&s);
    if (ptr != 0) return intersectsTimeRegion(*ptr);

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
    double *coords= new double[m_dimension];
    for (int i = 0; i < m_dimension; ++i) {
        coords[i]=h2*pre->m_pCoords[i]+h1*next->m_pCoords[i];
    }
    return TimePoint(coords,time,time,m_dimension);
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
    if(m_dimension==in.m_dimension-1){
        Region spatial(in.m_pLow,in.m_pHigh,m_dimension);
        if(in.m_pHigh[m_dimension]<m_points.front().m_startTime||in.m_pLow[m_dimension]>m_points.back().m_endTime){
            return false;
        }
        TimePoint tp=getPointAtTime(in.m_pLow[m_dimension]);
        return tp.intersectsShape(spatial);
    }
    else if(m_dimension==in.m_dimension) {
        for (int i = 0; i < m_points.size(); i++) {
            if (m_points[i].intersectsShape(in)) {
                return true;
            }
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
    out.makeInfinite(m_dimension);
    for(int i=0;i<m_points.size();i++){
        out.combinePoint(m_points[i]);
    }
}


void Trajectory::getMBC(SpatialIndex::MBC &out) const {
    if(m_points.size()<=1){
        std::cerr<<"WARNING: getting MBC at a Trajectory with 0 or 1 points\n";
        out.makeInfinite(m_dimension+1);
        return;
    }
    double startx=m_points[0].m_pCoords[0],starty=m_points[0].m_pCoords[1],startt=m_points[0].m_startTime;
    double endx=m_points.back().m_pCoords[0],endy=m_points.back().m_pCoords[1],endt=m_points.back().m_startTime;
    double avx=(endx-startx)/(endt-startt),avy=(endy-starty)/(endt-startt);
    TimePoint p1=*m_points.begin(),p2=m_points.back();
    double rd=0,rv=0;
    double vx,vy;
    for(int i=0;i<m_points.size();i++){
        TimePoint p;
        double ptime=m_points[i].m_startTime;
        double prd,prv;
        if(ptime-startt>0){
            vx=(m_points[i].m_pCoords[0]-startx)/(ptime-startt);
            vy=(m_points[i].m_pCoords[1]-starty)/(ptime-startt);
            prv=std::sqrt((vx-avx)*(vx-avx)+(vy-avy)*(vy-avy));
            if(prv>rv) rv=prv;
        }
        if(endt-m_points[i].m_startTime>0) {
            vx = (endx - m_points[i].m_pCoords[0]) / (endt - ptime);
            vy = (endy - m_points[i].m_pCoords[1]) / (endt - ptime);
            prv=std::sqrt((vx-avx)*(vx-avx)+(vy-avy)*(vy-avy));
            if(prv>rv) rv=prv;
        }
        prd=m_points[i].getMinimumDistance(*TimePoint::makemid(p1,p2,ptime));
        if(prd>rd) rd=prd;
    }
    out=MBC(p1.m_pCoords,p2.m_pCoords,startt,endt,m_dimension,rd,rv);
}

void Trajectory::getMBRfull(SpatialIndex::Region &out) const {
    out.makeInfinite(m_dimension+1);
    double *pc=new double[m_dimension+1];
    for(int i=0;i<m_points.size();i++){
        for(int d=0;d<m_dimension;d++) pc[d]=m_points[i].m_pCoords[d];
        pc[m_dimension]=m_points[i].m_startTime;
        out.combinePoint(Point(pc,m_dimension+1));
//        std::cout<<out<<std::endl;
    }
    delete[](pc);
}

void Trajectory::getTimeMBR(SpatialIndex::TimeRegion &out) const {
    out.makeInfinite(m_dimension);
    for(int i=0;i<m_points.size();i++){
        out.combinePoint(m_points[i]);
        if(m_points[i].m_startTime<out.m_startTime) out.m_startTime=m_points[i].m_startTime;
        if(m_points[i].m_startTime>out.m_endTime) out.m_endTime=m_points[i].m_startTime;
    }
}

const std::pair<int, double> findMaximumDistance(const vector<SpatialIndex::TimePoint>& points) {
    SpatialIndex::TimePoint firstpoint=points[0];
    SpatialIndex::TimePoint lastpoint=points[points.size()-1];
    int index=0;  //index to be returned
    double Mdist=-1; //the Maximum distance to be returned

    //distance calculation
    for(int i=1;i<points.size()-1;i++){ //traverse through second point to second last point
        double Dist=SpatialIndex::TimePoint::makemid(firstpoint,lastpoint,points[i].m_startTime)->getMinimumDistance(points[i]);
        if (Dist>Mdist){
            Mdist=Dist;
            index=i;
        }
    }
    return std::make_pair(index, Mdist);
}
std::vector<SpatialIndex::TimePoint> Trajectory::simplifyWithRDP(std::vector<SpatialIndex::TimePoint> &Points,
                                                                 double epsilon) {
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
        r.emplace_back(Points[Points.size()-1]);
        return r;
    }
}

std::vector<Trajectory> Trajectory::cuttraj(std::vector<SpatialIndex::TimePoint> mask){
    vector<TimePoint> seg;
    vector<Trajectory> res;
    auto iter1=m_points.begin();
    auto iter2=mask.begin();
    assert(m_points[0]==mask[0]);
    iter2++;
    for(;iter2!=mask.end();iter2++){
        seg.clear();
        while(iter1!=m_points.end()&&iter1->m_startTime<=iter2->m_startTime){
//            std::cerr<<"placed"<<*iter1<<"\n";
            seg.emplace_back(*iter1);
            iter1++;
        }
        if(seg.size()>0) {
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

std::vector<Trajectory> Trajectory::getSegments(double threshold) {
    auto mask=simplifyWithRDP(m_points,threshold);
    return cuttraj(mask);
}

double Trajectory::getArea() const{ return 0;}

inline double theF(double c1,double c2,double c3,double c4,double t){
    //the c4 should be the length of that time period
    double delta=4*c1*c3-c2*c2;
    if(delta<=1e-3){
        return (2*c1*t+c2)*sqrt(std::max(0.0,c1*t*t+c2*t+c3))/4/c1/c4;
    }
    else {
        return asinh((2 * c1 * t + c2) / sqrt(delta)) * delta / 8 / c1 / sqrt(c1) / c4
               + (2 * c1 * t + c2) * sqrt(std::max(0.0,c1 * t * t + c2 * t + c3)) / 4 / c1 / c4;
    }
}
double theD(double c1,double c2,double c3,double c4,double t){
    //the c4 should be the length of that time period
    return sqrt(c1*t*t+c2*t+c3)/c4;
}


double Trajectory::line2lineDistance(const SpatialIndex::TimePoint &p1s, const SpatialIndex::TimePoint &p1e,
                                const SpatialIndex::TimePoint &p2s, const SpatialIndex::TimePoint &p2e) {
    if(p1s.m_startTime!=p2s.m_startTime|p1e.m_startTime!=p2e.m_startTime)
        throw Tools::IllegalStateException("line2lineDistance: time period not the same");
    double ts = p1s.m_startTime, te = p1e.m_endTime;
    double dxs=p1s.m_pCoords[0]-p2s.m_pCoords[0];
    double dys=p1s.m_pCoords[1]-p2s.m_pCoords[1];
    double dxe=p1e.m_pCoords[0]-p2e.m_pCoords[0];
    double dye=p1e.m_pCoords[1]-p2e.m_pCoords[1];
    double c1=sq(dxs-dxe)+sq(dys-dye),
            c2=2*((dxe*ts-dxs*te)*(dxs-dxe)+(dye*ts-dys*te)*(dys-dye)),
            c3=sq(dxe*ts-dxs*te)+sq(dye*ts-dys*te),
            c4=te-ts;
    if(c1==0){
        return std::sqrt(sq(dxs)+sq(dys));
    }else{
            return (theF(c1,c2,c3,c4,te)-theF(c1,c2,c3,c4,ts));
    }
}

int getRealm(const SpatialIndex::Region &r,const Point &p1,const Point &p2){
    // 7 8 9
    // 4 5 6
    // 1 2 3
    if(p1.m_dimension!=2)
        throw Tools::NotSupportedException("Trajectory::getRealm:only for dimension 2");
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

TimePoint* cutByLine(const SpatialIndex::TimePoint &ps, const SpatialIndex::TimePoint &pe,double value,int axis){
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
            xyt[2]=d2 * ps.m_startTime + d1 * pe.m_endTime;
        }else{
            xyt[0]=d2 * ps.m_pCoords[otheraxis] + d1 * pe.m_pCoords[otheraxis];
            xyt[1]=value;
            xyt[2]=d2 * ps.m_startTime + d1 * pe.m_endTime;
        }
        return new TimePoint(xyt, xyt[2], xyt[2], 2);
    }
}
std::vector<std::pair<TimePoint,TimePoint>> cutByRealm(const SpatialIndex::TimePoint &ps, const SpatialIndex::TimePoint &pe,
                                                     const SpatialIndex::Region &r){
    double xd1=r.m_pLow[0],xd2=r.m_pHigh[0],yd1=r.m_pLow[1],yd2=r.m_pHigh[1];
    std::vector<std::pair<TimePoint,TimePoint>> res;
    std::vector<std::pair<TimePoint,TimePoint>> tmp;
    res.emplace_back(std::make_pair(ps,pe));
    for(auto line:res){
        TimePoint *tp=cutByLine(line.first,line.second,xd1,0);
        if(tp!= nullptr){
            tmp.emplace_back(std::make_pair(line.first,*tp));
            tmp.emplace_back(std::make_pair(*tp,line.second));
        }
        else{
            tmp.push_back(line);
        }
    }
    res=tmp;tmp.clear();
    for(auto line:res){
        TimePoint *tp=cutByLine(line.first,line.second,xd2,0);
        if(tp!= nullptr){
            tmp.emplace_back(std::make_pair(line.first,*tp));
            tmp.emplace_back(std::make_pair(*tp,line.second));
        }
        else{
            tmp.push_back(line);
        }
    }
    res=tmp;tmp.clear();
    for(auto line:res){
        TimePoint *tp=cutByLine(line.first,line.second,yd1,1);
        if(tp!= nullptr){
            tmp.emplace_back(std::make_pair(line.first,*tp));
            tmp.emplace_back(std::make_pair(*tp,line.second));
        }
        else{
            tmp.push_back(line);
        }
    }
    res=tmp;tmp.clear();
    for(auto line:res){
        TimePoint *tp=cutByLine(line.first,line.second,yd2,1);
        if(tp!= nullptr){
            tmp.emplace_back(std::make_pair(line.first,*tp));
            tmp.emplace_back(std::make_pair(*tp,line.second));
        }
        else{
            tmp.push_back(line);
        }
    }
    res=tmp;tmp.clear();
//    for(auto seg:res){
//        cout<<seg.first<<"\n"<<seg.second<<"\n";
//    }
    return res;
}
double Trajectory::line2MBRDistance_impl(const SpatialIndex::TimePoint &ps, const SpatialIndex::TimePoint &pe,
                                    const SpatialIndex::Region &r,int sr) {
    double ts = ps.m_startTime, te = pe.m_startTime;
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
        return line2lineDistance(ps,pe,TimePoint(coord,ts,ts,2),TimePoint(coord,te,te,2));
    }
}

double Trajectory::line2MBRDistance(const SpatialIndex::TimePoint &ps, const SpatialIndex::TimePoint &pe,
                                    const SpatialIndex::Region &r) {
    //the line's time period should be in the MBR's time period
    assert(r.m_pLow[m_dimension]<=ps.m_startTime);
    assert(r.m_pHigh[m_dimension]>=pe.m_startTime);
    //check if need cutting
    int sr=getRealm(r,ps,pe);
    if(sr>0){
        return line2MBRDistance_impl(ps,pe,r,sr);
    }else{
        double sum =0;
        auto part=cutByRealm(ps,pe,r);
        for(const auto &p:part) {
            int tmpsr=getRealm(r,p.first,p.second);
            double tmpres=line2MBRDistance_impl(p.first,p.second,r,tmpsr);
//            cout<<tmpsr<<" "<<tmpres<<"\n";
            sum+=tmpres;
        }
        return sum;
    }
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

double Trajectory::line2MBCDistance(const SpatialIndex::TimePoint &ps, const SpatialIndex::TimePoint &pe,
                                    const SpatialIndex::MBC &r) {
////    auto pr1=r.getCenterRdAtTime(ps.m_startTime),
////        pr2=r.getCenterRdAtTime(pe.m_startTime);
//  simple version
    double x1=r.m_pLow[0],y1=r.m_pLow[1],t1=r.m_startTime;
    double x2=r.m_pHigh[0],y2=r.m_pHigh[1],t2=r.m_endTime;
    double xts=makemidmacro(x1,t1,x2,t2,ps.m_startTime),xte=makemidmacro(x1,t1,x2,t2,pe.m_startTime);
    double yts=makemidmacro(y1,t1,y2,t2,ps.m_startTime),yte=makemidmacro(y1,t1,y2,t2,pe.m_startTime);
    double p2Low[2]={xts,yts},p2High[2]={xte,yte};
    double mbcr1=std::min(std::min(r.m_rd,(ps.m_startTime-r.m_startTime)*r.m_rv),(r.m_endTime-ps.m_startTime)*r.m_rv);
    double mbcr2=std::min(std::min(r.m_rd,(pe.m_startTime-r.m_startTime)*r.m_rv),(r.m_endTime-pe.m_startTime)*r.m_rv);
    TimePoint p2s(p2Low,ps.m_startTime,ps.m_startTime,2),p2e(p2High,pe.m_startTime,pe.m_startTime,2);
    double s=line2lineDistance(ps,pe,p2s,p2e);
    double sm=mbcArea(r.m_startTime,r.m_endTime,r.m_rd,r.m_rv,ps.m_startTime,pe.m_startTime);
    if(ps.getMinimumDistance(p2s)>mbcr1&&pe.getMinimumDistance(p2e)>mbcr2) return s-sm;
    else return 0;
//    else {
//        assert(std::isfinite(r.m_rv));
//        double t0=r.m_startTime,t1=t0+r.m_rd/r.m_rv,t3=r.m_endTime,t2=t3-r.m_rd/r.m_rv;
//        double ts=ps.m_startTime,te=pe.m_startTime;
//        TimePoint p1s=ps,p1e=pe,p2s=r.getCenterRdAtTime(ts).first,
//                p2e=r.getCenterRdAtTime(te).first;
//        double dxs=p1s.m_pCoords[0]-p2s.m_pCoords[0];
//        double dys=p1s.m_pCoords[1]-p2s.m_pCoords[1];
//        double dxe=p1e.m_pCoords[0]-p2e.m_pCoords[0];
//        double dye=p1e.m_pCoords[1]-p2e.m_pCoords[1];
//        double c1=sq(dxs-dxe)+sq(dys-dye),
//                c2=2*((dxe*ts-dxs*te)*(dxs-dxe)+(dye*ts-dys*te)*(dys-dye)),
//                c3=sq(dxe*ts-dxs*te)+sq(dye*ts-dys*te),
//                c4=te-ts;
//        double rootdelta=c2*c2-4*c1*(c3-sq(r.m_rd*c4));
//        double interdelta;
//        if(rootdelta<0){
//            //root not exist
//            return theF(c1,c2,c3,c4,te)-theF(c1,c2,c3,c4,ts)-mbcArea(r,ts,te);
//        }else{
//    //        double root1=(-c2-sqrt(c2*c2-4*c1*(c3-sq(r.m_rd*c4))))/2/c1;
//    //        double root2=(-c2+sqrt(c2*c2-4*c1*(c3-sq(r.m_rd*c4))))/2/c1;
//
//            double inter1,inter2;
//            if((2*c1*t2+c2<0&&rootdelta<sq(2*c1*t2+c2))||(2*c1*t1+c2<0&&rootdelta<sq(2*c1*t1+c2))) {
//                //root1<t2
//                if (2*c1*t2+c2<0&&rootdelta<sq(2*c1*t2+c2)) {
//                    //check if inter1 and inter2 in range
//                    double a = c1 - sq(c4 * r.m_rv),
//                            b = c2 + sq(c4 * r.m_rv) * 2 * t3,
//                            c = c3 - sq(c4 * r.m_rv) * sq(t3);
////                    assert(a>0);
//                    interdelta=b*b-4*a*c;
//                    inter1 = (-b -(std::fabs(a)/a)*sqrt(b * b - 4 * a * c)) / 2 / a;
//                    inter2 = (-b + (std::fabs(a)/a)*sqrt(b * b - 4 * a * c)) / 2 / a;
//                    if(inter1<t3&&inter1>t2){}
//                    else {
//                        //inter1 and inter2 don't exist
//                        return theF(c1, c2, c3, c4, te) - theF(c1, c2, c3, c4, ts) - mbcArea(r, ts, te);
//                    }
//                }
//                if (2*c1*t1+c2<0&&rootdelta<sq(2*c1*t1+c2)) {
//                    double a = c1 - sq(c4 * r.m_rv),
//                            b = c2 + sq(c4 * r.m_rv) * 2 * t0,
//                            c = c3 - sq(c4 * r.m_rv) * sq(t0);
//                    inter1 = (-b -(std::fabs(a)/a)*sqrt(b * b - 4 * a * c)) / 2 / a;
//                    inter2 = (-b + (std::fabs(a)/a)*sqrt(b * b - 4 * a * c)) / 2 / a;
//                    if(inter2<t1&&inter2>t0){}
//                    else{
//                        //inter1 and inter2 don't exist
//                        return theF(c1, c2, c3, c4, te) - theF(c1, c2, c3, c4, ts) - mbcArea(r, ts, te);
//                    }
//                }
//                else{
//                    throw Tools::IllegalStateException("???");
//                }
//            }
//            else {
//                double root1=(-c2-sqrt(c2*c2-4*c1*(c3-sq(r.m_rd*c4))))/2/c1;
//                double root2=(-c2+sqrt(c2*c2-4*c1*(c3-sq(r.m_rd*c4))))/2/c1;
//                if (root1 < t1) {
//                    //calculte inter1
//                    double a = c1 - sq(c4 * r.m_rv),
//                            b = c2 + sq(c4 * r.m_rv) * 2 * t0,
//                            c = c3 - sq(c4 * r.m_rv) * sq(t0);
//                    if (a > 0) inter1 = (-b - sqrt(b * b - 4 * a * c)) / 2 / a;
//                    else inter1 = (-b + sqrt(b * b - 4 * a * c)) / 2 / a;
//                } else {//if(root1>=t1&&root1<t2)
//                    inter1 = root1;
//                }
//                if (root2 > t2) {
//                    //calculate inter2
//                    double a = c1 - sq(c4 * r.m_rv),
//                            b = c2 + sq(c4 * r.m_rv) * 2 * t3,
//                            c = c3 - sq(c4 * r.m_rv) * sq(t3);
//                    if (a > 0) inter2 = (-b + sqrt(b * b - 4 * a * c)) / 2 / a;
//                    else inter2 = (-b - sqrt(b * b - 4 * a * c)) / 2 / a;
//                } else {//if(root2>t1&&root2<=t2)
//                    inter2 = root2;
//                }
//            }
//            double sum=0;
//            if(ts<inter1){
//                double tlow=ts;
//                double thigh=std::min(inter1,te);
//                sum+= theF(c1,c2,c3,c4,thigh)-theF(c1,c2,c3,c4,tlow)-mbcArea(r,tlow,thigh);
//            }
//            if(te>inter2){
//                double tlow=std::max(inter2,ts);
//                double thigh=te;
//                sum+= theF(c1,c2,c3,c4,thigh)-theF(c1,c2,c3,c4,tlow)-mbcArea(r,tlow,thigh);
//            }
//            return sum;
//        }
//    }
}

double Trajectory::getMinimumDistance(const IShape& s) const{
    //using Integrative Squared Euclidean Distance
    const Trajectory* pTrajectory = dynamic_cast<const Trajectory*>(&s);
    if (pTrajectory != 0) return getMinimumDistance(*pTrajectory);

//    const TimePoint* ptp = dynamic_cast<const TimePoint*>(&s);
//    if (ptp != 0) return getMinimumDistance(*ptp);

    const MBC* pmbc= dynamic_cast<const MBC*>(&s);
    if(pmbc!=0) return getMinimumDistance(*pmbc);

    const SBR* psbr= dynamic_cast<const SBR*>(&s);
    if(psbr!=0) return getMinimumDistance(*psbr);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return getMinimumDistance(*pr);

    const ShapeList *psl = dynamic_cast<const ShapeList*>(&s);
    if (psl != 0) return getMinimumDistance(*psl);

    throw Tools::NotSupportedException(
            "Trajectory::...: Not implemented yet!"
    );
}


double Trajectory::getMinimumDistance(const SpatialIndex::Region &in) const {
        //note here we only calculate the partial distance
        //so just calculate intersection of time period
        double tstart, tend;
        tstart = std::max(m_startTime(), in.m_pLow[in.m_dimension - 1]);
        tend = std::min(m_endTime(), in.m_pHigh[in.m_dimension - 1]);
        if(tstart>=tend) return 1e300;
        fakeTpVector timedtraj(&m_points,tstart,tend);
        double sum = 0;
        for (int i = 0; i < timedtraj.m_size-1; i++) {
            double pd = line2MBRDistance(timedtraj[i],timedtraj[i+1],in);
//            std::cerr<<"MBR dist for"<<timedtraj.m_points[i].m_startTime<<"is"<<pd<<endl;
            sum+=pd;
        }
        if(std::isnan(sum)){
            std::cerr<<"producing nan here!\n";
        }
        return sum;
}

double Trajectory::getMinimumDistance(const SpatialIndex::MBC &in) const {
    assert(m_dimension==in.m_dimension);
    double tstart, tend;
    tstart = std::max(m_startTime(), in.m_startTime);
    tend = std::min(m_endTime(), in.m_endTime);
    if(tstart>=tend) return 1e300;
    fakeTpVector timedtraj(&m_points,tstart,tend);
    double sum = 0;
    for (int i = 0; i < timedtraj.m_size-1; i++) {
        double pd = line2MBCDistance(timedtraj[i],timedtraj[i+1],in);
        sum+=pd;
    }
    return sum;
}
double Trajectory::getMinimumDistance(const SpatialIndex::SBR &in) const {
    double tstart, tend;
    tstart = std::max(m_startTime(), in.m_startTime);
    tend = std::min(m_endTime(), in.m_endTime);
    if(tstart>=tend) return 1e300;
    fakeTpVector timedtraj(&m_points,tstart,tend);
    double sum = 0;
    Region rg1,rg2;
    for (int i = 0; i < timedtraj.m_size-1; i++) {
        in.getMBRAtTime(timedtraj[i].m_startTime,rg1);
        in.getMBRAtTime(timedtraj[i+1].m_startTime,rg2);
        double pd = 0.5*(timedtraj[i].getMinimumDistance(rg1)+timedtraj[i+1].getMinimumDistance(rg2));
        sum+=pd;
    }
    return sum;
}

double Trajectory::getMinimumDistance(const SpatialIndex::TimePoint &in) const {
    int last;
    for(last=1;m_points[last].m_startTime<in.m_startTime;last++);
    return in.getMinimumDistance(*TimePoint::makemid(m_points[last-1],m_points[last],in.m_startTime));
}

double Trajectory::getMinimumDistance(const SpatialIndex::Trajectory &in) const {
    if(m_startTime()>=in.m_endTime()||m_endTime()<=in.m_startTime())
        return 1e300;
    fakeTpVector timedTraj2(&in.m_points,m_startTime(),m_endTime());
    double cut1=timedTraj2[0].m_startTime,cut2=timedTraj2[timedTraj2.m_size-1].m_startTime;
    double sum=0;
    fakeTpVector frontTraj(&m_points,m_startTime(),cut1);
    fakeTpVector midTraj(&m_points,cut1,cut2);
    fakeTpVector backTraj(&m_points,cut2,m_endTime());
    if(frontTraj.m_size!=0){
        for(int i=0;i<frontTraj.m_size-1;i++){
            double pd=0.5*(frontTraj[i].getMinimumDistance(timedTraj2[0])
                    +frontTraj[i+1].getMinimumDistance(timedTraj2[0]))
                            *(frontTraj[i+1].m_startTime-frontTraj[i].m_startTime);
//            std::cerr<<"point dist for"<<frontTraj.m_points[i].m_startTime<<"is"<<pd<<endl;
            sum+=pd;
        }
    }
    if(backTraj.m_size!=0){
        for(int i=0;i<backTraj.m_size-1;i++){
            double pd=0.5*(backTraj[i].getMinimumDistance(timedTraj2[timedTraj2.m_size-1])
                      +backTraj[i+1].getMinimumDistance(timedTraj2[timedTraj2.m_size-1]))
                 *(backTraj[i+1].m_startTime-backTraj[i].m_startTime);
//            std::cerr<<"point dist for"<<backTraj.m_points[i].m_startTime<<"is"<<pd<<endl;
            sum+=pd;
        }
    }
    if(midTraj.m_size!=0) {
        double newtime = midTraj[0].m_startTime, lasttime = midTraj[0].m_startTime;
        auto iter1 = midTraj.m_vectorPointer->begin()+midTraj.m_is;
        auto iter2 = timedTraj2.m_vectorPointer->begin()+timedTraj2.m_is;
        TimePoint lastp1 = *iter1, lastp2 = *iter2, newp1, newp2;
        while (lasttime != timedTraj2[timedTraj2.m_size-1].m_startTime) {
            if ((iter1 + 1)->m_startTime == (iter2 + 1)->m_startTime) {
                newtime = (iter1 + 1)->m_startTime;
                newp1 = *(iter1 + 1);
                newp2 = *(iter2 + 1);
                iter1++;
                iter2++;
            } else if ((iter1 + 1)->m_startTime < (iter2 + 1)->m_startTime) {
                newtime = (iter1 + 1)->m_startTime;
                newp1 = *(iter1 + 1);
                newp2 = *TimePoint::makemid(*iter2, *(iter2 + 1), newtime);
                iter1++;
            } else {
                newtime = (iter2 + 1)->m_startTime;
                newp1 = *TimePoint::makemid(*iter1, *(iter1 + 1), newtime);
                newp2 = *(iter2 + 1);
                iter2++;
            }
            lasttime = newtime;
            double pd= line2lineDistance(lastp1, newp1, lastp2, newp2);
            sum+=pd;
            lastp1 = newp1;
            lastp2 = newp2;
        }
    }
    return sum;
}

double Trajectory::getMinimumDistance(const ShapeList &in) const {
    if(in.m_MBRList.size()==0&&in.m_MBCList.size()==0){
        std::cerr<<"???\n";
    }
    double sum=0;
    double ints,inte;
    double pd;
    if(in.m_shapeType==SpatialIndex::LeafBoundByMBR||!in.m_MBRList.empty()) {
        ints=in.m_MBRList.front()->m_pLow[m_dimension];
        inte=in.m_MBRList.back()->m_pHigh[m_dimension];
        Region smbr=*in.m_MBRList.front(),embr=*in.m_MBRList.back();
        for(auto &br:in.m_MBRList){
            pd=getMinimumDistance(*br);
            if(pd!=1e300) sum+=pd;
        }
        if(m_startTime()<ints){
            smbr.m_pLow[m_dimension]=m_startTime();
            smbr.m_pHigh[m_dimension]=ints;
            pd=getMinimumDistance(smbr);
            if(pd!=1e300) sum+=pd;
        }
        if(m_endTime()>inte){
            embr.m_pLow[m_dimension]=inte;
            embr.m_pHigh[m_dimension]=m_endTime();
            pd=getMinimumDistance(embr);
            if(pd!=1e300) sum+=pd;
        }
    }else if(in.m_shapeType==SpatialIndex::LeafBoundByMBC||!in.m_MBCList.empty()){
        Point sPoint,ePoint;
        ints=in.m_MBCList.front()->m_startTime;
        inte=in.m_MBCList.back()->m_endTime;
        sPoint=Point(in.m_MBCList.front()->m_pLow,m_dimension);
        ePoint=Point(in.m_MBCList.back()->m_pHigh,m_dimension);
        for(auto &bc:in.m_MBCList){
            pd=getMinimumDistance(*bc);
            if(pd!=1e300) sum+=pd;
        }
        if(m_startTime()<ints){
            double pLow[m_dimension+1],pHigh[m_dimension+1];
            for(int i=0;i<m_dimension;i++){
                pLow[i]=sPoint.m_pCoords[i];
                pHigh[i]=sPoint.m_pCoords[i];
            }
            pLow[m_dimension]=m_startTime();
            pHigh[m_dimension]=ints;
            Region sbr(pLow,pHigh,m_dimension+1);
            pd=getMinimumDistance(sbr);
            if(pd!=1e300) sum+=pd;
        }
        if(m_endTime()>inte){
            double pLow[m_dimension+1],pHigh[m_dimension+1];
            for(int i=0;i<m_dimension;i++){
                pLow[i]=ePoint.m_pCoords[i];
                pHigh[i]=ePoint.m_pCoords[i];
            }
            pLow[m_dimension]=inte;
            pHigh[m_dimension]=m_endTime();
            Region ebr(pLow,pHigh,m_dimension+1);
            pd=getMinimumDistance(ebr);
            if(pd!=1e300) sum+=pd;
        }
    } else{
        std::cerr<<"ShapeList without DataType?\n"<<in<<"\n";
    }
    return sum;
}

double Trajectory::getPeriodMinimumDistance(const SpatialIndex::Region &in,double MaxVelocity) const {
    double sum=getMinimumDistance(in);//midTraj
    //frontTraj
    Region mbr2d(in.m_pLow,in.m_pHigh,2);
    if(in.m_pLow[m_dimension]>m_startTime()){
        fakeTpVector frontTraj(&m_points,m_startTime(),in.m_pLow[m_dimension]);
        for(int i=0;i<frontTraj.m_size-1;i++){
            double r1=frontTraj[i].getMinimumDistance(mbr2d)-
                    frontTraj[i].m_startTime*MaxVelocity;
            double r2=frontTraj[i+1].getMinimumDistance(mbr2d)-
                    frontTraj[i+1].m_startTime*MaxVelocity;
            if(r1>=0&&r2>=0) sum+=0.5*(frontTraj[i+1].m_startTime-frontTraj[i].m_startTime)
                    *(r1+r2);
        }
    }
    if(in.m_pHigh[m_dimension]<m_endTime()){
        fakeTpVector backTraj(&m_points,in.m_pHigh[m_dimension],m_endTime());
        for(int i=0;i<backTraj.m_size-1;i++){
            double r1=backTraj[i].getMinimumDistance(mbr2d)-
                    backTraj[i].m_startTime*MaxVelocity;
            double r2=backTraj[i+1].getMinimumDistance(mbr2d)-
                    backTraj[i+1].m_startTime*MaxVelocity;
            if(r1>=0&&r2>=0) sum+=0.5*(backTraj[i+1].m_startTime-backTraj[i].m_startTime)
                                  *(r1+r2);
        }
    }
    return sum;
}

double Trajectory::getPeriodMinimumDistance(const SpatialIndex::MBC &in,double MaxVelocity) const {
    double sum=getMinimumDistance(in);//midTraj
    //frontTraj
    Point sPoint(in.m_pLow,2),ePoint(in.m_pHigh,2);
    if(in.m_startTime>m_startTime()){
        fakeTpVector frontTraj(&m_points,m_startTime(),in.m_startTime);
        for(int i=0;i<frontTraj.m_size-1;i++){
            double r1=frontTraj[i].getMinimumDistance(sPoint)
                    -frontTraj[i].m_startTime*MaxVelocity;
            double r2=frontTraj[i+1].getMinimumDistance(sPoint)
                    -frontTraj[i+1].m_startTime*MaxVelocity;
            if(r1>=0&&r2>=0) sum+=0.5*(frontTraj[i+1].m_startTime-frontTraj[i].m_startTime)
                                  *(r1+r2);
        }
    }
    if(in.m_endTime<m_endTime()){
        fakeTpVector backTraj(&m_points,in.m_endTime,m_endTime());
        for(int i=0;i<backTraj.m_size-1;i++){
            double r1=backTraj[i].getMinimumDistance(ePoint)
                    -backTraj[i].m_startTime*MaxVelocity;
            double r2=backTraj[i+1].getMinimumDistance(ePoint)
                    -backTraj[i+1].m_startTime*MaxVelocity;
            if(r1>=0&&r2>=0) sum+=0.5*(backTraj[i+1].m_startTime-backTraj[i].m_startTime)
                                  *(r1+r2);
        }
    }
    return sum;
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
             "," + std::to_string(p.m_startTime) + " ";
    }
    s += "\n";
    os<<s;
    return os;
}
std::vector<std::string> split(const std::string strtem,char a)
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

void Trajectory::loadFromString(std::string str) {
    std::vector<std::string> points=split(str,' ');
    for(auto p: points){
        std::vector<std::string> xyt=split(p,',');
        double xy[2];
        xy[0]=std::stod(xyt[0]);xy[1]=std::stod(xyt[1]);
        m_points.emplace_back(TimePoint(xy,std::stod(xyt[2]),std::stod(xyt[2]),2));
    }
}

void Trajectory::linkTrajectory(SpatialIndex::Trajectory other) {
    if(m_points.back().m_startTime==other.m_points.front().m_startTime){
        m_points.insert(m_points.end(),++other.m_points.begin(),other.m_points.end());
    }
    else if (other.m_points.back().m_startTime==m_points.front().m_startTime){
        m_points.insert(m_points.begin(),other.m_points.begin(),other.m_points.end()-1);
    }
    else{
        std::cerr<<m_points.back()<<" "<<other.m_points.front()<<"\n"
            <<other.m_points.back()<<" "<<m_points.front();
        throw Tools::IllegalStateException("Trajectory::linkTrajectory: the two trajectories to be linked should have a common point.");
    }
}