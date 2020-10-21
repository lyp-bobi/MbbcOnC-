//
// Created by chuang on 4/23/19.
//


#pragma once
#include <cmath>
#define subTrajFile "./subTrajFile.stj"

#define random(x, y) (((double)rand()/RAND_MAX)*(y-x)+x)
extern double calcuTime[10];
extern int testPhase;
extern int disttype;

#define Tristat

#define bip auto start = std::chrono::system_clock::now();
#define bbip auto end = std::chrono::system_clock::now();auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);calcuTime[testPhase]+=double(duration.count()) * std::chrono::microseconds::period::num/ std::chrono::microseconds::period::den;
using std::vector;
using std::pair;
using std::string;
using namespace SpatialIndex;

namespace SpatialIndex
{
    class SIDX_DLL xTrajectory: public Tools::IObject, public virtual IShape{

        public:
        class fakeTpVector{
        public:
            int m_is=0;
            int m_size=0;
            xPoint *m_front= nullptr;
            xPoint *m_back = nullptr;
            std::vector<xPoint> *m_vectorPointer;
            ~fakeTpVector(){
                delete(m_front);
                delete(m_back);
            }
            fakeTpVector(const std::vector<xPoint> *vector,double ts,double te)
            {
                m_vectorPointer= const_cast<std::vector<xPoint>*>(vector);
                if(ts==te) return;
                if(ts>m_vectorPointer->back().m_t || te < m_vectorPointer->front().m_t) return;
                m_is=0;
                int m_ie=m_vectorPointer->size()-1;
                while(m_vectorPointer->at(m_is).m_t < ts) m_is++;
                while(m_vectorPointer->at(m_ie).m_t > te) m_ie--;
                if(m_is!=0&& m_vectorPointer->at(m_is).m_t != ts){
                    double x=makemidmacro(m_vectorPointer->at(m_is-1).m_x, m_vectorPointer->at(m_is-1).m_t,
                                          m_vectorPointer->at(m_is).m_x, m_vectorPointer->at(m_is).m_t, ts);
                    double y=makemidmacro(m_vectorPointer->at(m_is-1).m_y, m_vectorPointer->at(m_is-1).m_t,
                                          m_vectorPointer->at(m_is).m_y, m_vectorPointer->at(m_is).m_t, ts);
                    m_front=new xPoint(x,y,ts);
                    m_is--;
                }
                if(m_ie!=m_vectorPointer->size()-1&& m_vectorPointer->at(m_ie).m_t != te){
                    double x=makemidmacro(m_vectorPointer->at(m_ie).m_x, m_vectorPointer->at(m_ie).m_t,
                                          m_vectorPointer->at(m_ie+1).m_x, m_vectorPointer->at(m_ie+1).m_t, te);
                    double y=makemidmacro(m_vectorPointer->at(m_ie).m_y, m_vectorPointer->at(m_ie).m_t,
                                          m_vectorPointer->at(m_ie+1).m_y, m_vectorPointer->at(m_ie+1).m_t, te);
                    m_back=new xPoint(x,y,te);
                    m_ie++;
                }
                m_size=m_ie-m_is+1;
            }
            xPoint& operator[](int n){
                if(n==0&&m_front!= nullptr) return *m_front;
                if(n==m_size-1&&m_back!= nullptr) return *m_back;
                return m_vectorPointer->at(m_is+n);
            }
            xPoint& front(){
                if(m_front!= nullptr) return *m_front;
                else return m_vectorPointer->at(m_is);
            }
            xPoint& back(){
                if(m_back!= nullptr) return *m_back;
                return m_vectorPointer->at(m_is+m_size-1);
            }
        };
        xTrajectory();
        explicit xTrajectory(std::vector<xPoint> &in);
        explicit xTrajectory(bool fakehead,bool fakeback,std::vector<xPoint> &in);
        xTrajectory(const xTrajectory& in);
        xTrajectory &operator=(const xTrajectory &r);

        virtual bool operator==(const xTrajectory &r) const;


        //
        // IObject interface
        //
        virtual xTrajectory *clone();

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

        virtual void getPartialxTrajectory(double tstart, double tend, SpatialIndex::xTrajectory &out) const;

        virtual double getArea() const;
        virtual double getMinimumDistance(const IShape& in) const;



        virtual double getMinimumDistance(const xPoint& in) const;
        virtual double getMinimumDistance(const xMBR& in) const;
        virtual double getMaxSED(const xMBR& in) const;
        virtual double getMaxSED(const xMBC& in) const;
        virtual double getMinimumDistance(const xMBC& in) const;
        virtual double getMinimumDistance(const xTrajectory& in) const;
        virtual double getMinimumDistanceInner(const xTrajectory& in) const;

        static double line2lineIED(const xPoint &p1s, const xPoint &p1e, const xPoint &p2s, const xPoint &p2e);
        static double line2lineMinSED(const xPoint &p1s, const xPoint &p1e, const xPoint &p2s, const xPoint &p2e);
        static double line2MBRDistance(const xPoint &ps,const xPoint &pe,const xMBR &r);
        static double line2MBRIED_impl(const xPoint &ps, const xPoint &pe, const xMBR &r, int sr);
        static double line2MBRMinSED(const xPoint &ps, const xPoint &pe, const xMBR &r);
        static double line2MBRMaxSED(const xPoint &ps, const xPoint &pe, const xMBR &r);
        static double line2MBRMinSED_impl(const xPoint &ps, const xPoint &pe, const xMBR &r, int sr);
        static double line2MBCDistance(const xPoint &ps,const xPoint &pe,const xMBC &r);

        virtual double getFrontIED(xMBR smbr, double MaxVelocity) const;
        virtual double getBackIED(xMBR embr, double MaxVelocity) const;
        virtual double getFrontIED(double x,double y,double ts, double MaxVelocity) const;
        virtual double getBackIED(double x,double y,double te, double MaxVelocity) const;

        virtual double getMidIED(const xMBR &sbr,const xMBR &ebr,double MaxVelocity,double queryVelocity=-1);
        virtual double getMidIED(const xMBC &sbc, const xMBC &ebc,double MaxVelocity,double queryVelocity=-1);
        virtual double getMidIED(const xPoint &sp, const xPoint &ep,double MaxVelocity,double queryVelocity=-1);

        virtual double getStaticIED(double x, double y, double t1, double t2) const;
        virtual double getStaticIED(xMBR in,double ints, double inte) const;
        virtual double getStaticMaxSED(double x, double y, double t1, double t2) const;
        virtual double getStaticMaxSED(xMBR in,double ints, double inte) const;

        virtual double getInferredNodeMinimumIED(const xMBR &in, double MaxVelocity,double queryMaxVelocity=0) const;

        virtual double getNodeMinimumDistance(const xMBR &in,double MaxVelocity) const;
        virtual double getLeafMinimumDistance(const xMBR &in, double MaxVelocity) const;
        virtual double getLeafMinimumDistance(const xMBC &in, double MaxVelocity) const;


        virtual bool intersectsxMBR(const xMBR& in) const;
        virtual bool intersectsxCylinder(const xCylinder& in) const;
        virtual bool intersectsxTrajectory(const xTrajectory& in) const;

        virtual void combinexTrajectory(const xTrajectory& in);
        virtual bool containsxTrajectory(const xTrajectory& in);
        virtual void getCombinedxTrajectory(xTrajectory& out, const xTrajectory& in) const;
        double getInterTime(const xPoint &ps, const xPoint &pe,xMBR &r,double t_o, double vmax) const;

        virtual void getxMBC(xMBC& out) const;
        xPoint getPointAtTime(double time) const;
        static std::vector<std::vector<SpatialIndex::xPoint>> simplifyWithRDPN(const std::vector<SpatialIndex::xPoint>& Points, int numPart);
        static std::vector<SpatialIndex::xPoint> simplifyWithRDP(const std::vector<SpatialIndex::xPoint>& Points, double threshold);
        std::vector<xTrajectory> cuttraj(std::vector<SpatialIndex::xPoint>);
        std::vector<xTrajectory> getSegments(double len) const;
        std::vector<xTrajectory> getHybridSegments(double len) const;
        std::vector<xTrajectory> getRDPSegments(double len) const;
        std::vector<xTrajectory> getStaticSegments(double len) const;
        std::vector<xTrajectory> getStaticSegmentsCut(double len) const;
        std::vector<xTrajectory> getFixedSegments(int len=169) const;
        std::vector<xTrajectory> getGlobalSegmentsCut(double len) const;
        std::vector<xTrajectory> getItself() const;

        void linkxTrajectory(xTrajectory &other);

        static int cutTrajsIntoFile(std::vector<std::pair<id_type, xTrajectory> > &trajs,double segLen, int strat=0,std::string filename=subTrajFile);
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
            pair<id_type,vector<xTrajectory>> getNext(){
                id_type id;
                vector<xTrajectory> v;
                if(!hasNext()){
                    return std::make_pair(id,v);
                }
                std::string s;
                s=nextline;
                xTrajectory tmp;
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
        
        inline double m_startTime() const{return m_points.front().m_t;}
        inline double m_endTime() const { return m_points.back().m_t;}
        inline double m_dist() const{
            double res=0;
            for(int i=1;i<m_points.size();i++){
                res += std::sqrt(sq(m_points[i].m_x-m_points[i-1].m_x)
                        +sq(m_points[i].m_y-m_points[i-1].m_y));
            }
            return res;
        }

        xPoint randomPoint(){
            int i = int(random(0,m_points.size()-1));
            return m_points[i];
        }

        std::string toString() const ;
        void loadFromString(std::string s);

        static int getPhase(const SpatialIndex::xMBR &r,const xPoint &p1,const xPoint &p2);
        static std::vector<std::pair<xPoint,xPoint>> cutByPhase(const SpatialIndex::xPoint &ps, const SpatialIndex::xPoint &pe,
                                                                  const SpatialIndex::xMBR &r);

        std::vector<xPoint> m_points;
        bool m_fakehead=false,m_fakeback=false;

        friend SIDX_DLL std::ostream& operator<<(std::ostream& os,const xTrajectory &r);

        static const uint32_t m_dimension=2;
        virtual void makeInfinite(uint32_t dimension);
    private:
};
typedef Tools::PoolPointer<xTrajectory> xTrajectoryPtr;
SIDX_DLL std::ostream& operator<<(std::ostream& os, const xTrajectory& r);
}
std::vector<std::string> split(const std::string &strtem,char a);