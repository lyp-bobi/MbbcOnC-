//
// Created by chuang on 4/23/19.
//

#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>

#include <spatialindex/SpatialIndex.h>
double calcuTime[2]={0,0};
int testPhase=0;

int disttype = 0;
//0 for IED, 1 for MaxSED

using namespace SpatialIndex;
using std::vector;
using std::cout;
using std::endl;
using std::sqrt;

Trajectory::Trajectory() {
}
Trajectory::Trajectory(std::vector<SpatialIndex::STPoint>& in) {
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
    std::vector<STPoint> p(size);
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

STPoint Trajectory::getPointAtTime(const double time) const {
    if(time<m_points.front().m_time||time>m_points.back().m_time){
        throw Tools::IllegalArgumentException(
                "Trajectory::getPointAtTime: time"+std::to_string(time)+"is illegal."
                );}
    if(m_points.size()==1){
        return m_points[0];
    }
    auto pre =m_points.begin(),next=m_points.begin();
    while(next->m_time-pre->m_time<0.01) next++;
    while(next->m_time<time&&next!=m_points.end()){
        pre=next;
        while(next->m_time-pre->m_time<0.01) next++;
    }
    double h1= (time-pre->m_time)/(next->m_time-pre->m_time);
    double h2= (next->m_time-time)/(next->m_time-pre->m_time);
    double *coords= new double[m_dimension];
    for (int i = 0; i < m_dimension; ++i) {
        coords[i]=h2*pre->m_pCoords[i]+h1*next->m_pCoords[i];
    }
    return STPoint(coords,time,m_dimension);
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
        Region spatial(in.m_pLow,in.m_pHigh,m_dimension);
        if(in.m_pHigh[m_dimension]<m_points.front().m_time||in.m_pLow[m_dimension]>m_points.back().m_time){
            return false;
        }
        STPoint tp=getPointAtTime(in.m_pLow[m_dimension]);
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
        prd=m_points[i].getMinimumDistance(*STPoint::makemid(p1,p2,ptime));
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
        double Dist=SpatialIndex::STPoint::makemid(firstpoint,lastpoint,points[i].m_time)->getMinimumDistance(points[i]);
        if (Dist>Mdist){
            Mdist=Dist;
            index=i;
        }
    }
    return std::make_pair(index, Mdist);
}
std::vector<SpatialIndex::STPoint> Trajectory::simplifyWithRDP(std::vector<SpatialIndex::STPoint> &Points,
                                                                 double epsilon) {
    if(Points.size()<3){  //base case 1
        return Points;
    }
    std::pair<int, double> maxDistance=findMaximumDistance(Points);
    if(maxDistance.second>=epsilon){
        int index=maxDistance.first;
        vector<SpatialIndex::STPoint>::iterator it=Points.begin();
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

std::vector<Trajectory> Trajectory::cuttraj(std::vector<SpatialIndex::STPoint> mask){
    vector<STPoint> seg;
    vector<Trajectory> res;
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
        if(!seg.empty()) {
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


double Trajectory::line2lineIED(const SpatialIndex::STPoint &p1s, const SpatialIndex::STPoint &p1e,
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
        return new STPoint(xyt, xyt[2], 2);
    }
}
std::vector<std::pair<STPoint,STPoint>> cutByRealm(const SpatialIndex::STPoint &ps, const SpatialIndex::STPoint &pe,
                                                     const SpatialIndex::Region &r){
    double xd1=r.m_pLow[0],xd2=r.m_pHigh[0],yd1=r.m_pLow[1],yd2=r.m_pHigh[1];
    std::vector<std::pair<STPoint,STPoint>> res;
    std::vector<std::pair<STPoint,STPoint>> tmp;
    res.emplace_back(std::make_pair(ps,pe));
    for(const auto &line:res){
        STPoint *tp=cutByLine(line.first,line.second,xd1,0);
        if(tp!= nullptr){
            tmp.emplace_back(std::make_pair(line.first,*tp));
            tmp.emplace_back(std::make_pair(*tp,line.second));
        }
        else{
            tmp.push_back(line);
        }
    }
    res=tmp;tmp.clear();
    for(const auto &line:res){
        STPoint *tp=cutByLine(line.first,line.second,xd2,0);
        if(tp!= nullptr){
            tmp.emplace_back(std::make_pair(line.first,*tp));
            tmp.emplace_back(std::make_pair(*tp,line.second));
        }
        else{
            tmp.push_back(line);
        }
    }
    res=tmp;tmp.clear();
    for(const auto &line:res){
        STPoint *tp=cutByLine(line.first,line.second,yd1,1);
        if(tp!= nullptr){
            tmp.emplace_back(std::make_pair(line.first,*tp));
            tmp.emplace_back(std::make_pair(*tp,line.second));
        }
        else{
            tmp.push_back(line);
        }
    }
    res=tmp;tmp.clear();
    for(const auto &line:res){
        STPoint *tp=cutByLine(line.first,line.second,yd2,1);
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
double Trajectory::line2MBRIED_impl(const SpatialIndex::STPoint &ps, const SpatialIndex::STPoint &pe,
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

double Trajectory::line2MBRDistance(const SpatialIndex::STPoint &ps, const SpatialIndex::STPoint &pe,
                                    const SpatialIndex::Region &r) {
    //the line's time period should be in the MBR's time period
    assert(r.m_pLow[m_dimension]<=ps.m_time);
    assert(r.m_pHigh[m_dimension]>=pe.m_time);
    //check if need cutting
    if(disttype==0) {
        int sr = getRealm(r, ps, pe);
        if (sr > 0) {
            return line2MBRIED_impl(ps, pe, r, sr);
        } else {
            double sum = 0;
            auto part = cutByRealm(ps, pe, r);
            for (const auto &p:part) {
                int tmpsr = getRealm(r, p.first, p.second);
                double tmpres = line2MBRIED_impl(p.first, p.second, r, tmpsr);
//            cout<<tmpsr<<" "<<tmpres<<"\n";
                sum += tmpres;
            }
            return sum;
        }
    }else if(disttype==1){
        Region mbr2d(r.m_pLow,r.m_pHigh,r.m_dimension-1);
        return std::max(ps.getMinimumDistance(mbr2d),pe.getMinimumDistance(mbr2d));
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

double Trajectory::line2MBCDistance(const SpatialIndex::STPoint &ps, const SpatialIndex::STPoint &pe,
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
        double t1=r.m_startTime+r.m_rd/r.m_rv,t2=r.m_endTime-r.m_rd/r.m_rv;
        double ma=0;
        if(t1>ps.m_time&&t1<pe.m_time){
            double x=makemidmacro(ps.m_pCoords[0],ps.m_time,pe.m_pCoords[0],pe.m_time,t1);
            double y=makemidmacro(ps.m_pCoords[1],ps.m_time,pe.m_pCoords[1],pe.m_time,t1);
            double coord[2]={x,y};
            ma=r.getMinimumDistance(STPoint(coord,t1,2));
        }
        if(t2>ps.m_time&&t2<pe.m_time){
            double x=makemidmacro(ps.m_pCoords[0],ps.m_time,pe.m_pCoords[0],pe.m_time,t2);
            double y=makemidmacro(ps.m_pCoords[1],ps.m_time,pe.m_pCoords[1],pe.m_time,t2);
            double coord[2]={x,y};
            ma=r.getMinimumDistance(STPoint(coord,t2,2));
        }
        return std::max(ma,std::max(r.getMinimumDistance(ps),r.getMinimumDistance(pe)));
    }
    else{throw Tools::NotSupportedException("Wrong distance");}

}

double Trajectory::getMinimumDistance(const IShape& s) const{
    //using Integrative Squared Euclidean Distance
    const Trajectory* pTrajectory = dynamic_cast<const Trajectory*>(&s);
    if (pTrajectory != nullptr) return getMinimumDistance(*pTrajectory);

//    const STPoint* ptp = dynamic_cast<const STPoint*>(&s);
//    if (ptp != 0) return getMinimumDistance(*ptp);

    const MBC* pmbc= dynamic_cast<const MBC*>(&s);
    if(pmbc!=nullptr) return getMinimumDistance(*pmbc);

    const SBR* psbr= dynamic_cast<const SBR*>(&s);
    if(psbr!=nullptr) return getMinimumDistance(*psbr);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != nullptr) return getMinimumDistance(*pr);

    const ShapeList *psl = dynamic_cast<const ShapeList*>(&s);
    if (psl != nullptr) return getMinimumDistance(*psl);

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
        fakeTpVector timedtraj(&m_points,tstart,tend);
        double sum = 0,max=0;
        for (int i = 0; i < timedtraj.m_size-1; i++) {
            double pd = line2MBRDistance(timedtraj[i],timedtraj[i+1],in);
            if(disttype==0){
                sum+=pd;
            }
            else if(disttype==1){
                max=std::max(max,pd);
            }
        }
        if(disttype==0)
            return sum;
        else if(disttype==1)
            return max;
        else{throw Tools::NotSupportedException("Wrong distance");}
}

double Trajectory::getMinimumDistance(const SpatialIndex::MBC &in) const {
    assert(m_dimension==in.m_dimension);
    double tstart, tend;
    tstart = std::max(m_startTime(), in.m_startTime);
    tend = std::min(m_endTime(), in.m_endTime);
    if(tstart>=tend) return 1e300;
    fakeTpVector timedtraj(&m_points,tstart,tend);
    if(disttype==0) {
        double sum = 0;
        for (int i = 0; i < timedtraj.m_size - 1; i++) {
            double pd = line2MBCDistance(timedtraj[i], timedtraj[i + 1], in);
            sum += pd;
        }
        return sum;
    }
    else if(disttype==1){
        double max=0;
        for (int i = 0; i < timedtraj.m_size - 1; i++) {
            double pd = line2MBCDistance(timedtraj[i], timedtraj[i + 1], in);
            max=std::max(max,pd);
        }
        return max;
    }
    else{throw Tools::NotSupportedException("Wrong distance");}
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
        in.getMBRAtTime(timedtraj[i].m_time,rg1);
        in.getMBRAtTime(timedtraj[i+1].m_time,rg2);
        double pd = 0.5*(timedtraj[i].getMinimumDistance(rg1)+timedtraj[i+1].getMinimumDistance(rg2));
        sum+=pd;
    }
    return sum;
}

double Trajectory::getMinimumDistance(const SpatialIndex::STPoint &in) const {
    int last;
    for(last=1;m_points[last].m_time<in.m_time;last++);
    return STPoint::makemid(m_points[last-1],m_points[last],in.m_time)->getMinimumDistance(in);
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
        double pd=getMinimumDistance(timedTraj2[0].m_pCoords[0],timedTraj2[0].m_pCoords[1],m_startTime(),cut1);
        if(disttype==0) {
            sum += pd;
        }else if(disttype==1){
            max=std::max(max,pd);
        }
    }
    if(m_endTime()>cut2){
        double pd=getMinimumDistance(timedTraj2[timedTraj2.m_size-1].m_pCoords[0],timedTraj2[timedTraj2.m_size-1].m_pCoords[1],cut2,m_endTime());
        if(disttype==0) {
            sum += pd;
        }else if(disttype==1){
            max=std::max(max,pd);
        }
    }
    if(midTraj.m_size!=0) {
        double newtime = midTraj[0].m_time, lasttime = midTraj[0].m_time;
        auto iter1 = midTraj.m_vectorPointer->begin()+midTraj.m_is;
        auto iter2 = timedTraj2.m_vectorPointer->begin()+timedTraj2.m_is;
        STPoint lastp1 = *iter1, lastp2 = *iter2, newp1, newp2;
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

double Trajectory::getMinimumDistance(const ShapeList &in) const {
    double sum=0,max=0;
    double ints,inte;
    double pd;
    if(in.m_shapeType==SpatialIndex::LeafBoundByMBR||!in.m_MBRList.empty()) {
        ints=in.m_MBRList.front()->m_pLow[m_dimension];
        inte=in.m_MBRList.back()->m_pHigh[m_dimension];
        Region smbr=*in.m_MBRList.front(),embr=*in.m_MBRList.back();
        if(disttype==0) {
            for (const auto &br:in.m_MBRList) {
                pd = getMinimumDistance(*br);
                if (pd != 1e300) sum += pd;
            }
            if (m_startTime() < ints) {
                smbr.m_pLow[m_dimension] = m_startTime();
                smbr.m_pHigh[m_dimension] = ints;
                pd = getMinimumDistance(smbr);
                if (pd != 1e300) sum += pd;
            }
            if (m_endTime() > inte) {
                embr.m_pLow[m_dimension] = inte;
                embr.m_pHigh[m_dimension] = m_endTime();
                pd = getMinimumDistance(embr);
                if (pd != 1e300) sum += pd;
            }
        }
        else if(disttype==1){
            for (const auto &br:in.m_MBRList) {
                pd = getMinimumDistance(*br);
                if (pd != 1e300) max=std::max(max,pd);
            }
            if (m_startTime() < ints) {
                smbr.m_pLow[m_dimension] = m_startTime();
                smbr.m_pHigh[m_dimension] = ints;
                pd = getMinimumDistance(smbr);
                if (pd != 1e300) max=std::max(max,pd);
            }
            if (m_endTime() > inte) {
                embr.m_pLow[m_dimension] = inte;
                embr.m_pHigh[m_dimension] = m_endTime();
                pd = getMinimumDistance(embr);
                if (pd != 1e300) max=std::max(max,pd);
            }
        }
        else{throw Tools::NotSupportedException("Wrong distance");}
    }else if(in.m_shapeType==SpatialIndex::LeafBoundByMBC||!in.m_MBCList.empty()){
        if(disttype==0) {
            Point sPoint, ePoint;
            ints = in.m_MBCList.front()->m_startTime;
            inte = in.m_MBCList.back()->m_endTime;
            sPoint = Point(in.m_MBCList.front()->m_pLow, m_dimension);
            ePoint = Point(in.m_MBCList.back()->m_pHigh, m_dimension);
            for (const auto &bc:in.m_MBCList) {
                pd = getMinimumDistance(*bc);
                if (pd != 1e300) sum += pd;
            }
            if (m_startTime() < ints) {
                pd = getMinimumDistance(sPoint.m_pCoords[0], sPoint.m_pCoords[1], m_startTime(), ints);
                if (pd != 1e300) sum += pd;
            }
            if (m_endTime() > inte) {
                pd = getMinimumDistance(ePoint.m_pCoords[0], ePoint.m_pCoords[1], inte, m_endTime());
                if (pd != 1e300) sum += pd;
            }
        }
        else if(disttype==1){
            Point sPoint, ePoint;
            ints = in.m_MBCList.front()->m_startTime;
            inte = in.m_MBCList.back()->m_endTime;
            sPoint = Point(in.m_MBCList.front()->m_pLow, m_dimension);
            ePoint = Point(in.m_MBCList.back()->m_pHigh, m_dimension);
            for (const auto &bc:in.m_MBCList) {
                pd = getMinimumDistance(*bc);
                if (pd != 1e300) max=std::max(max,pd);
            }
            if (m_startTime() < ints) {
                pd = getMinimumDistance(sPoint.m_pCoords[0], sPoint.m_pCoords[1], m_startTime(), ints);
                if (pd != 1e300) max=std::max(max,pd);
            }
            if (m_endTime() > inte) {
                pd = getMinimumDistance(ePoint.m_pCoords[0], ePoint.m_pCoords[1], inte, m_endTime());
                if (pd != 1e300) max=std::max(max,pd);
            }
        }
        else{throw Tools::NotSupportedException("Wrong distance");}
    }
    if(disttype==0){
        return sum;
    }
    else if(disttype==1){
        return max;
    }else{
        std::cerr<<"ShapeList without DataType?\n"<<in<<"\n";
    }
}

double Trajectory::getMinimumDistance(double x,double y, double t1,double t2) const {
    double tstart, tend;
    tstart = std::max(m_startTime(), t1);
    tend = std::min(m_endTime(), t2);
    if(tstart>=tend) return 1e300;
    fakeTpVector timedtraj(&m_points,tstart,tend);
    double sum = 0,max=0;
    double xy[2]={x,y};
    STPoint ps(xy,0,2),pe(xy,0,2);
    max=timedtraj[0].getMinimumDistance(ps);
    for (int i = 0; i < timedtraj.m_size-1; i++) {
        if(disttype==0) {
            ps.m_time=timedtraj[i].m_time;
            pe.m_time=timedtraj[i+1].m_time;
            double pd = line2lineIED(timedtraj[i], timedtraj[i + 1], ps, pe);
            sum += pd;
        }else if (disttype==1){
            max=std::max(max,timedtraj[i+1].getMinimumDistance(ps));
        }
    }
    if(disttype==0)
        return sum;
    else if(disttype==1)
        return max;
    else
        throw Tools::NotSupportedException("Wrong distance");
}

double Trajectory::minSED(const SpatialIndex::Region &in) const {
    assert(in.m_dimension==3);
    Region rg2d(in.m_pLow,in.m_pHigh,in.m_dimension-1);
    fakeTpVector part(&m_points,in.m_pLow[in.m_dimension-1],in.m_pHigh[in.m_dimension-1]);
    double min = 1e300;
    for(int i=0;i<part.m_size;i++){
        double dist=part[i].getMinimumDistance(rg2d);
        if(dist<min) min=dist;
    }
    // WARNING: This is not accurate, only a higher bound!
    if(min==1e300)
        return -1;
    return min;
}

double Trajectory::minSED(const SpatialIndex::MBC &in) const {
    assert(in.m_dimension==2);
    fakeTpVector part(&m_points,in.m_startTime,in.m_endTime);
    double min = 1e300;
    for(int i=0;i<part.m_size;i++){
        double dist=part[i].getMinimumDistance(in);
        if(dist<min) min=dist;
    }
    // WARNING: This is not accurate, only a higher bound!
    if(min==1e300)
        return -1;
    return min;
}

double Trajectory::getLeafMinimumDistance(const SpatialIndex::Region &in, double MaxVelocity) const {
    double sum=getMinimumDistance(in);//midTraj
    //frontTraj
    Region sembr=in;
    if(in.m_pLow[m_dimension]>m_startTime()){
        fakeTpVector frontTraj(&m_points,m_startTime(),in.m_pLow[m_dimension]);
        for(int i=0;i<frontTraj.m_size-1;i++){
            sembr.m_pLow[m_dimension]=frontTraj[i].m_time;
            sembr.m_pHigh[m_dimension]=frontTraj[i+1].m_time;
            double pd=line2MBRDistance(frontTraj[i],frontTraj[i+1],sembr);
            double r1=(frontTraj[i].m_time-m_startTime())*MaxVelocity;
            double r2=(frontTraj[i+1].m_time-m_startTime())*MaxVelocity;
            double minus = 0.5*(frontTraj[i+1].m_time-frontTraj[i].m_time)*(r1+r2);
            sum+=std::max(0.0,pd-minus);
        }
    }
    if(in.m_pHigh[m_dimension]<m_endTime()){
        fakeTpVector backTraj(&m_points,in.m_pHigh[m_dimension],m_endTime());
        for(int i=0;i<backTraj.m_size-1;i++){
            sembr.m_pLow[m_dimension]=backTraj[i].m_time;
            sembr.m_pHigh[m_dimension]=backTraj[i+1].m_time;
            double pd=line2MBRDistance(backTraj[i],backTraj[i+1],sembr);
            double r1=(m_endTime()-backTraj[i].m_time)*MaxVelocity;
            double r2=(m_endTime()-backTraj[i+1].m_time)*MaxVelocity;
            double minus = 0.5*(backTraj[i+1].m_time-backTraj[i].m_time)*(r1+r2);
            sum+=std::max(0.0,pd-minus);
        }
    }
    return sum;
}

double Trajectory::getLeafMinimumDistance(const SpatialIndex::MBC &in, double MaxVelocity) const {
    double sum=getMinimumDistance(in);//midTraj
    //frontTraj
    STPoint sPoint(in.m_pLow,0,2),sPoint2(in.m_pLow,0,2),ePoint(in.m_pHigh,0,2),ePoint2(in.m_pHigh,0,2);
    if(in.m_startTime>m_startTime()){
        fakeTpVector frontTraj(&m_points,m_startTime(),in.m_startTime);
        for(int i=0;i<frontTraj.m_size-1;i++){
            sPoint.m_time=frontTraj[i].m_time;
            sPoint2.m_time=frontTraj[i+1].m_time;
            double pd=line2lineIED(sPoint,sPoint2,frontTraj[i],frontTraj[i+1]);
            double r1=(frontTraj[i].m_time-m_startTime())*MaxVelocity;
            double r2=(frontTraj[i+1].m_time-m_startTime())*MaxVelocity;
            double minus = 0.5*(frontTraj[i+1].m_time-frontTraj[i].m_time)*(r1+r2);
            sum+=std::max(0.0,pd-minus);
        }
    }
    if(in.m_endTime<m_endTime()){
        fakeTpVector backTraj(&m_points,in.m_endTime,m_endTime());
        for(int i=0;i<backTraj.m_size-1;i++){
            ePoint.m_time=backTraj[i].m_time;
            ePoint2.m_time=backTraj[i+1].m_time;
            double pd=line2lineIED(ePoint,ePoint2,backTraj[i],backTraj[i+1]);
            double r1=(m_endTime()-backTraj[i].m_time)*MaxVelocity;
            double r2=(m_endTime()-backTraj[i+1].m_time)*MaxVelocity;
            double minus = 0.5*(backTraj[i+1].m_time-backTraj[i].m_time)*(r1+r2);
            sum+=std::max(0.0,pd-minus);
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
             "," + std::to_string(p.m_time) + " ";
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
    for(const auto &p: points){
        std::vector<std::string> xyt=split(p,',');
        double xy[2];
        xy[0]=std::stod(xyt[0]);xy[1]=std::stod(xyt[1]);
        m_points.emplace_back(STPoint(xy,std::stod(xyt[2]),2));
    }
}

void Trajectory::linkTrajectory(SpatialIndex::Trajectory other) {
    if(m_points.back().m_time==other.m_points.front().m_time){
        m_points.insert(m_points.end(),++other.m_points.begin(),other.m_points.end());
    }
    else if (other.m_points.back().m_time==m_points.front().m_time){
        m_points.insert(m_points.begin(),other.m_points.begin(),other.m_points.end()-1);
    }
    else{
        std::cerr<<m_points.back()<<" "<<other.m_points.front()<<"\n"
            <<other.m_points.back()<<" "<<m_points.front();
        throw Tools::IllegalStateException("Trajectory::linkTrajectory: the two trajectories to be linked should have a common point.");
    }
}