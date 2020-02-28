//
// Created by chuang on 4/23/19.
//


#pragma once

#include "ShapeList.h"
#define subTrajFile "./subTrajFile.stj"
extern double calcuTime[10];
extern int testPhase;
extern int disttype;

#define bip auto start = std::chrono::system_clock::now();
#define bbip auto end = std::chrono::system_clock::now();auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);calcuTime[testPhase]+=double(duration.count()) * std::chrono::microseconds::period::num/ std::chrono::microseconds::period::den;
using std::vector;
using std::pair;
using std::string;
using namespace SpatialIndex;

namespace SpatialIndex
{
    class SIDX_DLL Trajectory: public Tools::IObject, public virtual IShape{

        public:
        class fakeTpVector{
        public:
            int m_is=0;
            int m_size=0;
            STPoint *m_front= nullptr;
            STPoint *m_back = nullptr;
            std::vector<STPoint> *m_vectorPointer;
            ~fakeTpVector(){
                delete(m_front);
                delete(m_back);
            }
            fakeTpVector(const std::vector<STPoint> *vector,double ts,double te)
            {
                m_vectorPointer= const_cast<std::vector<STPoint>*>(vector);
                if(ts==te) return;
                if(ts>m_vectorPointer->back().m_time||te<m_vectorPointer->front().m_time) return;
                m_is=0;
                int m_ie=m_vectorPointer->size()-1;
                while(m_vectorPointer->at(m_is).m_time<ts) m_is++;
                while(m_vectorPointer->at(m_ie).m_time>te) m_ie--;
                if(m_is!=0&&m_vectorPointer->at(m_is).m_time!=ts){
                    double x=makemidmacro(m_vectorPointer->at(m_is-1).m_pCoords[0],m_vectorPointer->at(m_is-1).m_time,
                                          m_vectorPointer->at(m_is).m_pCoords[0],m_vectorPointer->at(m_is).m_time,ts);
                    double y=makemidmacro(m_vectorPointer->at(m_is-1).m_pCoords[1],m_vectorPointer->at(m_is-1).m_time,
                                          m_vectorPointer->at(m_is).m_pCoords[1],m_vectorPointer->at(m_is).m_time,ts);
                    double xy[2]={x,y};
                    m_front=new STPoint(xy,ts,2);
                    m_is--;
                }
                if(m_ie!=m_vectorPointer->size()-1&&m_vectorPointer->at(m_ie).m_time!=te){
                    double x=makemidmacro(m_vectorPointer->at(m_ie).m_pCoords[0],m_vectorPointer->at(m_ie).m_time,
                                          m_vectorPointer->at(m_ie+1).m_pCoords[0],m_vectorPointer->at(m_ie+1).m_time,te);
                    double y=makemidmacro(m_vectorPointer->at(m_ie).m_pCoords[1],m_vectorPointer->at(m_ie).m_time,
                                          m_vectorPointer->at(m_ie+1).m_pCoords[1],m_vectorPointer->at(m_ie+1).m_time,te);
                    double xy[2]={x,y};
                    m_back=new STPoint(xy,te,2);
                    m_ie++;
                }
                m_size=m_ie-m_is+1;
            }
            STPoint& operator[](int n){
                if(n==0&&m_front!= nullptr) return *m_front;
                if(n==m_size-1&&m_back!= nullptr) return *m_back;
                return m_vectorPointer->at(m_is+n);
            }
            STPoint& front(){
                if(m_front!= nullptr) return *m_front;
                else return m_vectorPointer->at(m_is);
            }
            STPoint& back(){
                if(m_back!= nullptr) return *m_back;
                return m_vectorPointer->at(m_is+m_size-1);
            }
        };
        Trajectory();
        explicit Trajectory(std::vector<STPoint> &in);
        explicit Trajectory(bool fakehead,bool fakeback,std::vector<STPoint> &in);
        Trajectory(const Trajectory& in);
        Trajectory &operator=(const Trajectory &r);

        virtual bool operator==(const Trajectory &r) const;


        //
        // IObject interface
        //
        virtual Trajectory *clone();

        //
        // ISerializable interface
        //
        virtual uint32_t getByteArraySize() const;

        virtual void loadFromByteArray(const uint8_t *data);

        virtual void storeToByteArray(uint8_t **data, uint32_t &len);



        //
        // IShape interface
        //
        virtual bool intersectsShape(const IShape& in) const;
        virtual bool containsShape(const IShape& in) const;
        virtual bool touchesShape(const IShape& in) const;
        virtual void getCenter(Point& out) const;
        virtual uint32_t getDimension() const;
        virtual void getMBR(Region& out) const;

        virtual void getPartialTrajectory(double tstart, double tend, SpatialIndex::Trajectory &out) const;

        virtual double getArea() const;
        virtual double getMinimumDistance(const IShape& in) const;



        virtual double getMinimumDistance(const STPoint& in) const;
        virtual double getMinimumDistance(const Region& in) const;
        virtual double getMaxSED(const Region& in) const;
        virtual double getMaxSED(const MBC& in) const;
        virtual double getMinimumDistance(const MBC& in) const;
        virtual double getMinimumDistance(const Trajectory& in) const;
        virtual double getMinimumDistance(const ShapeList& in,bool hasPrev=false,bool hasNext=false,double MaxVelocity=1e300) const;

        static double line2lineIED(const STPoint &p1s, const STPoint &p1e, const STPoint &p2s, const STPoint &p2e);
        static double line2lineMinSED(const STPoint &p1s, const STPoint &p1e, const STPoint &p2s, const STPoint &p2e);
        static double line2MBRDistance(const STPoint &ps,const STPoint &pe,const Region &r);
        static double line2MBRIED_impl(const STPoint &ps, const STPoint &pe, const Region &r, int sr);
        static double line2MBRMinSED(const STPoint &ps, const STPoint &pe, const Region &r);
        static double line2MBRMaxSED(const STPoint &ps, const STPoint &pe, const Region &r);
        static double line2MBRMinSED_impl(const STPoint &ps, const STPoint &pe, const Region &r, int sr);
        static double line2MBCDistance(const STPoint &ps,const STPoint &pe,const MBC &r);

        virtual double getFrontIED(const Region smbr, double MaxVelocity) const;
        virtual double getBackIED(const Region embr, double MaxVelocity) const;
        virtual double getFrontIED(double x,double y,double ts, double MaxVelocity) const;
        virtual double getBackIED(double x,double y,double te, double MaxVelocity) const;

        virtual double getMidIED(const Region &sbr,const Region &ebr,double MaxVelocity,double queryVelocity=-1);
        virtual double getMidIED(const MBC &sbc, const MBC &ebc,double MaxVelocity,double queryVelocity=-1);
        virtual double getMidIED(const STPoint &sp, const STPoint &ep,double MaxVelocity,double queryVelocity=-1);

        virtual double getStaticIED(double x, double y, double t1, double t2) const;
        virtual double getStaticIED(const Region in,double ints, double inte) const;
        virtual double getStaticMaxSED(double x, double y, double t1, double t2) const;
        virtual double getStaticMaxSED(const Region in,double ints, double inte) const;

        virtual double getInferredNodeMinimumIED(const Region &in, double MaxVelocity,double queryMaxVelocity=0) const;

        virtual double getNodeMinimumDistance(const Region &in,double MaxVelocity) const;
        virtual double getLeafMinimumDistance(const Region &in, double MaxVelocity) const;
        virtual double getLeafMinimumDistance(const MBC &in, double MaxVelocity) const;


        virtual bool intersectsTimeRegion(const TimeRegion& in) const;
        virtual bool intersectsRegion(const Region& in) const;
        virtual bool intersectsCylinder(const Cylinder& in) const;
        virtual bool intersectsTrajectory(const Trajectory& in) const;

        virtual void combineTrajectory(const Trajectory& in);
        virtual bool containsTrajectory(const Trajectory& in);
        virtual void getCombinedTrajectory(Trajectory& out, const Trajectory& in) const;
        double getInterTime(const STPoint &ps, const STPoint &pe,Region &r,double t_o, double vmax) const;

        virtual void getMBC(MBC& out) const;
        virtual void getMBRfull(Region& out) const;
        virtual void getTimeMBR(TimeRegion& out) const;
        STPoint getPointAtTime(double time) const;
        static std::vector<std::vector<SpatialIndex::STPoint>> simplifyWithRDPN(const std::vector<SpatialIndex::STPoint>& Points, int numPart);
        static std::vector<SpatialIndex::STPoint> simplifyWithRDP(const std::vector<SpatialIndex::STPoint>& Points, double threshold);
        std::vector<Trajectory> cuttraj(std::vector<SpatialIndex::STPoint>);
        std::vector<Trajectory> getSegments(double len) const;
        std::vector<Trajectory> getHybridSegments(double len) const;
        std::vector<Trajectory> getRDPSegments(double len) const;
        std::vector<Trajectory> getStaticSegments(double len) const;
        std::vector<Trajectory> getStaticSegmentsCut(double len) const;
        std::vector<Trajectory> getFixedSegments(int len=140) const;
        std::vector<Trajectory> getGlobalSegmentsCut(double len) const;
        void linkTrajectory(Trajectory &other);

        static int cutTrajsIntoFile(std::vector<std::pair<id_type, Trajectory> > &trajs,double segLen,std::string filename=subTrajFile);
        class subTrajStream{
        public:
            std::ifstream inFile;
            string nextline="";
            subTrajStream(std::string str=subTrajFile){
                inFile.open(str,std::ios::in);
            }
            bool hasNext(){
                return nextline!="END";
            }
            pair<id_type,vector<Trajectory>> getNext(){
                id_type id;
                vector<Trajectory> v;
                if(!hasNext()){
                    return std::make_pair(id,v);
                }
                std::string s;
                s=nextline;
                Trajectory tmp;
                bool isId=false;
                while(true){
                    if(s=="SubTraj"){
                        isId=true;
                    }
                    else if(s=="ESubTraj"){
                        std::getline(inFile,s);
                        break;
                    }
                    else if(isId){
                        id=std::stoll(s);
                        isId=false;
                    }
                    else if(!s.empty()){
                        tmp.loadFromString(s);
                        v.emplace_back(tmp);
                    } else{
                        if(!s.empty()) std::cout<<"some err";
                    }
                    std::getline(inFile,s);
                }
                nextline=s;
                return std::make_pair(id,v);
            }
            ~subTrajStream(){
                inFile.close();
            }
        };
        
        inline double m_startTime() const{return m_points.front().m_time;}
        inline double m_endTime() const { return m_points.back().m_time;}

        std::string toString() const ;
        void loadFromString(std::string s);

        static int getPhase(const SpatialIndex::Region &r,const Point &p1,const Point &p2);
        static std::vector<std::pair<STPoint,STPoint>> cutByPhase(const SpatialIndex::STPoint &ps, const SpatialIndex::STPoint &pe,
                                                                  const SpatialIndex::Region &r);

        std::vector<STPoint> m_points;
        bool m_fakehead=false,m_fakeback=false;

        friend SIDX_DLL std::ostream& operator<<(std::ostream& os,const Trajectory &r);

        static const uint32_t m_dimension=2;
        virtual void makeInfinite(uint32_t dimension);
    private:
};
typedef Tools::PoolPointer<Trajectory> TrajectoryPtr;
SIDX_DLL std::ostream& operator<<(std::ostream& os, const Trajectory& r);
}
std::vector<std::string> split(const std::string &strtem,char a);