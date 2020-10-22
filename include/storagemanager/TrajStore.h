//
// Created by Chuang on 2019/6/10.
//

#pragma once


#include <spatialindex/SpatialIndex.h>
#include <cstring>

using namespace SpatialIndex;
using namespace SpatialIndex::StorageManager;
using std::vector;

extern bool bSecondaryIndex;

namespace SpatialIndex
{
    class trajStat{
    private:
        trajStat(){};
        ~trajStat(){delete singleton;}
        static trajStat* singleton;
    public:
        double bt=0;
        double M=0;//total time
        long lineCount=0;
        long trajCount=0;
        double tl=0;//time len of a segment
        double jt=0;
        double v=0;
        double vmax = 0;
        double minx=1e300,maxx=-1e300,miny=1e300,maxy=-1e300,mint=1e300,maxt=-1e300;
        double Dx=0,Dy=0,Dt=0;
        double Sr=0;
        double P=0;
        double dist=0;
        double Df =2;
        double f = 50,fp=170;
        string dataset = "";
        static trajStat* instance();
        void init(){
            bt=0;
            M=0;
            long lineCount=0;
            long trajCount=0;
            tl=0;//time len of a segment
            jt=0;
            v=0;
            minx=1e300,maxx=-1e300,miny=1e300,maxy=-1e300,mint=1e300,maxt=-1e300;
            Dx=0,Dy=0,Dt=0;
            Sr=0;
            P=0;
            dist=0;
            Df =2;
            f = 50;
        }
        void usedata(string str){
            dataset = str;
            if (dataset == "gl"){
                M = 1.31808e+08;
                lineCount = 20069299;
                trajCount = 20619;
                tl = 6.56764;
                jt = 6392.55;
                Sr = 0.5;
                Df = 1.25;
                P = 60000;
                std::cerr<<"use gl sta\n";
            }else if(dataset == "td"){
                M = 5.20574e+09;
                lineCount = 16205956;
                trajCount = 10267;
                tl = 321.224;
                jt = 507037;
                Sr = 0.5;
                Df = 1.56;
                P = 533315;
                std::cerr<<"use td sta\n";
            }
            else if (dataset == "od"){
                Sr = 10000;
                Df = 1.4;
                P = 5000;
                std::cerr<<"use od sta\n";
            }
        }

        double rd(double x){
            if (dataset == "gl"){
                return 2e-05 * x + 0.0023;
            }else if(dataset == "td"){
                return 1E-05 * x + 0.0003;
            }
            else if (dataset == "od"){
                return -0.0326 * x*x + 20.289*x - 78.702;
            }
            return vv()*x;
        }

        double dd(double x){
            if (dataset == "gl"){
                return 4.5e-5*x - 1e-5/4000*x*x;
            }else if(dataset == "td"){
                return 4e-5*x;
            }
            else if (dataset == "od"){
                return 55*x;
            }
            return vv()*x;
        }

        double vv(){
            if (dataset == "gl"){
                return dd(bt)/bt;
            }else if(dataset == "td"){
                return dd(bt)/bt;
            }
            else if (dataset == "od"){
                return dd(bt)/bt;
            }
            return v;
        }

        double Lx(){
            if (dataset == "gl"){

            }else if(dataset == "td"){

            }
            else if (dataset == "od"){

            }else{
                return pow(M_PI*Sr*Sr*P*M*dd(bt)/bt/bt/f,0.33) + dd(bt);
            }
        }

        double Lt(){
            if (dataset == "gl"){

            }else if(dataset == "td"){

            }
            else if (dataset == "od"){

            }else{
                return (pow(M_PI*Sr*Sr*P*M*dd(bt)/bt/bt/f,0.33) + dd(bt))*bt/dd(bt);
            }
        }

        double knncost(double _bt,int k,double qt, bool out = false){
            bt = _bt;
            double leafs = M/bt/f;
            double nt = pow(leafs * sq(vv())*sq(P) / sq((M_PI * sq(Sr))),1.0/3);
            double ltc = P/nt;
            double lt = ltc+bt;
            double lxc = ltc * vv();
            double lx = lt * vv();
            double Nq = min(1.0,(qt+jt) /P)  * trajCount;
            double rk =  pow(k/Nq, Df/2)*Sr;
            double Rk = rk*qt/2 + qt/4* sqrt(2*sq(dd(qt)+4*sq(rk)));
            double Rnq = Rk +rd(bt)*qt;

            double nmin = k;
            double nmax = 100*k;
            while (nmax - nmin > 0.1) {
                double mid = (nmax + nmin) / 2;
                double rn =  pow(mid/Nq, Df/2)*Sr;
                double Rn = rn*qt/2 + qt/4* sqrt(2*sq(dd(qt)+4*sq(rn)));
                if (Rn > Rnq) {
                    nmax = mid;
                } else {
                    nmin = mid;
                }
            }
            double nq = (nmax+nmin)/2;

            double pruneCost = sq(lx + 2*Rk/qt + dd(qt))*(lt+qt)   *(M/bt/f/P/pow(M_PI*sq(Sr),Df/2));
            double filterCost = nq*(1+qt/fp/lt);
            if (out){
                std::cerr<<"nt "<<nt <<" lt "<<lt<<" lx "<<lx<<"\n";
                std::cerr<<"expect cost is "<<pruneCost<<"\t"<<filterCost<<"\n";
            }
            return pruneCost+filterCost/8;
        }

        void set(double _bt, double _M, long _lineCount, long _trajCount,
                 double _tl, double _jt, double _v, double _minx, double _maxx,
                 double _miny, double _maxy, double _mint, double _maxt,
                 double _Dx, double _Dy, double _Dt, double _dist)
        {
            bt=_bt; M=_M; lineCount=_lineCount; trajCount=_trajCount;
            tl=_tl; jt=_jt; v=_v; minx=_minx; maxx=_maxx;
            miny=_miny; maxy=_maxy; mint=_mint; maxt=_maxt;
            Dx=_Dx; Dy=_Dy; Dt=_Dt; dist=_dist;
        }
        void output(){
            std::cerr<<bt<<","<< M<<","<< lineCount<<","<< trajCount<<","<<
                     tl<<","<< jt<<","<< v<<","<< minx<<","<< maxx<<","<<
                     miny<<","<< maxy<<","<< mint<<","<< maxt<<","<<
                     Dx<<","<< Dy<<","<< Dt<<","<< dist<<"\n";
        }
        string toString(){
            ostringstream ostream;
            ostream<<bt<<" "<< M<<" "<< lineCount<<" "<< trajCount<<" "<<
                   tl<<" "<< jt<<" "<< v<<" "<< minx<<" "<< maxx<<" "<<
                   miny<<" "<< maxy<<" "<< mint<<" "<< maxt<<" "<<
                   Dx<<" "<< Dy<<" "<< Dt<<" "<< dist;
            return ostream.str();
        }
        void fromString(string s){
            istringstream istream(s);
            istream>>bt>> M>> lineCount>> trajCount>>
                   tl>> jt>> v>> minx>> maxx>>
                   miny>> maxy>> mint>> maxt>>
                   Dx>> Dy>> Dt>> dist;
            return;
        }
        friend SIDX_DLL std::ostream& operator<<(std::ostream& os,const trajStat &r);
    };
    SIDX_DLL std::ostream& operator<<(std::ostream& os, const trajStat& r);

    class SIDX_DLL XZ3Enocder{
    private:
        XZ3Enocder();
        ~XZ3Enocder(){delete singleton;}
        static XZ3Enocder* singleton;
    public:
        long encode(double x,double y,double z);
        double m_xmin,m_xmax,m_ymin,m_ymax,m_zmin,m_zmax;
        uint32_t m_length;
        static XZ3Enocder* instance();
    };
    namespace StorageManager {
        class tsExternalSorter
        {
        public:
            class Record
            {
            public:
                Record();
                Record(const Region& shape, id_type id,id_type pvId,id_type ntId, uint32_t len,
                        uint8_t* pData, uint32_t s,uint32_t level,MBC* mbc= nullptr);
                ~Record();

                bool operator<(const Record& r) const;

                void storeToFile(Tools::TemporaryFile& f);
                void loadFromFile(Tools::TemporaryFile& f);

                struct SortAscending : public std::binary_function<Record* const, Record* const, bool>
                {
                    bool operator()(Record* const r1, Record* const r2)
                    {
                        if (*r1 < *r2) return true;
                        else return false;
                    }
                };

            public:
                uint32_t m_level;
                Region m_r;
                MBC m_mbc;
                id_type m_id,m_pvId,m_ntId;
                uint32_t m_len;
                uint8_t* m_pData;
                uint32_t m_s;
            };

        public:
            tsExternalSorter(uint32_t u32PageSize, uint32_t u32BufferPages);
            virtual ~tsExternalSorter();

            void insert(Record* r);
            void sort();
            Record* getNextRecord();
            uint64_t getTotalEntries() const;

        private:
            class PQEntry
            {
            public:
                PQEntry(Record* r, uint32_t u32Index) : m_r(r), m_u32Index(u32Index) {}

                struct SortAscending : public std::binary_function<const PQEntry&, const PQEntry&, bool>
                {
                    bool operator()(const PQEntry& e1, const PQEntry& e2)
                    {
                        if (*(e1.m_r) < *(e2.m_r)) return true;
                        else return false;
                    }
                };

                Record* m_r;
                uint32_t m_u32Index;
            };

        private:
            bool m_bInsertionPhase;
            uint32_t m_u32PageSize;
            uint32_t m_u32BufferPages;
            Tools::SmartPointer<Tools::TemporaryFile> m_sortedFile;
            std::list<Tools::SmartPointer<Tools::TemporaryFile> > m_runs;
            std::vector<Record*> m_buffer;
            uint64_t m_u64TotalEntries;
            uint32_t m_stI;
        };//tsExternalSorter
        class ircData:public IData{
        public:
            id_type m_id;
            Region m_br;
            MBC m_bc;
            uint8_t* m_pData;
            uint32_t m_dataLength;
            ircData(id_type id,Region &br,MBC &bc)
                :m_id(id),m_br(br),m_bc(bc){}
            virtual ircData* clone(){
                return new ircData(m_id,m_br,m_bc);
            }
            virtual id_type getIdentifier() const{
                return m_id;
            }
            virtual void getShape(IShape** out) const{
                if(m_bc.m_dimension>0){
                    *out = new MBC(m_bc);
                }else{
                    *out = new Region(m_br);
                }
            }
            virtual void getData(uint32_t& len, uint8_t** data) const{
                len = m_dataLength;
                *data = 0;

                if (m_dataLength > 0)
                {
                    *data = new uint8_t[m_dataLength];
                    memcpy(*data, m_pData, m_dataLength);
                }
            }
        };
        class SBBStream:public IDataStream{
        public:
            std::fstream m_idFile;
            std::fstream m_brFile;
            std::fstream m_bcFile;
            bool m_binput=true;
            bool m_busembc=false;
            id_type m_size;
            id_type m_cur;
            string m_s;
            SBBStream(std::string &str){
                m_idFile.open(str+".id",std::ios::out);
                m_brFile.open(str+".mbr",std::ios::out);
                m_bcFile.open(str+".mbc",std::ios::out);
                m_size=0;
                m_cur=0;
                m_s=str;
            }
            bool hasNext(){
//                std::cerr<<m_cur<<m_size;
                return m_cur<m_size;
            }
            void inputSBB(id_type id, Region* br, MBC* bc=nullptr){
                m_size+=1;
//                return;
                if(!m_binput){
                    m_idFile.open(m_s+".id",std::ios::out);
                    m_brFile.open(m_s+".mbr",std::ios::out);
                    m_bcFile.open(m_s+".mbc",std::ios::out);
                    m_binput=true;
                    rewind();
                }
                m_idFile<<id<<"\n";
                m_brFile<<br->toString()<<"\n";
                if(m_busembc||bc!= nullptr){
                    m_busembc=true;
                    m_bcFile<<bc->toString()<<"\n";
                }
            }
            uint32_t size() override {
                return m_size;
            }
            void endInput(){
                m_binput=false;
                m_idFile<<"END"<<std::endl;
                m_idFile.flush();
                m_idFile.close();
                m_idFile.clear();
                m_idFile.open(m_s+".id",std::ios::in);
                m_brFile<<"END"<<std::endl;
                m_brFile.flush();
                m_brFile.close();
                m_brFile.clear();
                m_brFile.open(m_s+".mbr",std::ios::in);
                m_bcFile<<"END"<<std::endl;
                m_bcFile.flush();
                m_bcFile.close();
                m_bcFile.clear();
                m_bcFile.open(m_s+".mbc",std::ios::in);
                rewind();
                std::cerr<<m_s<<" input finished. We have "<<m_size<<" sub-bounding boxes.\n";
            }
            ircData* getNext() override {
                if(m_binput){
                    endInput();
                }
                id_type id;
                Region br;
                MBC bc;
                m_cur++;
                std::string s;
                std::getline(m_idFile,s);
                id=std::stoll(s);
                std::getline(m_brFile,s);
                br.loadFromString(s);
                if(m_busembc){
                    std::getline(m_bcFile,s);
                    bc.loadFromString(s);
                }
                if(m_size==m_cur){
                    std::cerr<<"SBB loading finished.\n";
                    m_idFile.flush();
                    m_idFile.close();
                    m_idFile.clear();
                    m_brFile.flush();
                    m_brFile.close();
                    m_brFile.clear();
                    m_bcFile.flush();
                    m_bcFile.close();
                    m_bcFile.clear();
                    m_size=0;
                }
//                if(m_size-m_cur<10){
//                    std::cerr<<m_cur<<m_size<<std::endl<<id<<std::endl<<br<<std::endl<<bc<<std::endl;
//                }
                return new ircData(id,br,bc);
            }
            void rewind() override{
                m_cur=0;
                m_idFile.clear();
                m_idFile.seekg(std::ios::beg);
                m_brFile.clear();
                m_brFile.seekg(std::ios::beg);
                m_bcFile.clear();
                m_bcFile.seekg(std::ios::beg);
            }
            ~SBBStream(){
                m_idFile.close();
                m_brFile.close();
                m_bcFile.close();
            }
        };

        class SIDX_DLL TrajStore:public IStorageManager{
        public:
            ~TrajStore();
            TrajStore(string &name,IStorageManager *store,uint32_t pageSize,int maxseg=100);
            void flush();
            void loadByteArray(const id_type page, uint32_t& len, uint8_t** data){
                auto start = std::chrono::system_clock::now();
                m_pStorageManager->loadByteArray(page,len,data);
                auto end = std::chrono::system_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                m_IOtime+=double(duration.count()) * std::chrono::microseconds::period::num
                        / std::chrono::microseconds::period::den;
            }//for inner nodes
            void storeByteArray(id_type& page, const uint32_t len, const uint8_t* const data){
                m_pStorageManager->storeByteArray(page,len,data);
            }//for inner nodes
            void deleteByteArray(const id_type page){
                m_pStorageManager->deleteByteArray(page);
            }//should not be used in bulkload mode i guess
            void releaseTmp(){
//                delete m_stream; //error on this line
                for(const auto e:m_entries) delete e.second;

                m_entries.clear();
                m_entryMbcs.clear();
                m_entryMbrs.clear();
                m_part2node.clear();
            };
            void cleanStatistic(){
                m_trajIO=0;m_indexIO=0;
                m_leaf1=0;m_leaf2=0;
                m_loadedTraj=0;
                m_boundingVisited=0;
                m_maxVelocity=0;
                m_IOtime=0;
            }
            class Entry{
            public:
                id_type m_page;
                uint32_t m_start;
                uint32_t m_len;
                id_type m_pvId;
                id_type m_ntId;
                Entry(id_type page,uint32_t start,uint32_t len,id_type pvId,id_type ntId);
                std::string toString();
                Entry(string &s);
            };
            id_type getTrajId(id_type id){return id/m_maxTrajSegs;}
            id_type getSegId(id_type id,uint32_t segnum){return id*m_maxTrajSegs+segnum;}
            void loadTrajs(vector<std::pair<id_type, Trajectory> > &trajs,double segpara,bool idFirst=false,bool output=true,int method=3);
            void loadSegments(vector<std::pair<id_type, vector<Trajectory>> > &trajs,bool idFirst=false,bool output=true);
            void loadSegments(std::string file=subTrajFile,bool idFirst=false,bool output=true);
            const Trajectory getTraj(id_type &id);
            const Trajectory getTrajByTime(id_type &id,double tstart,double tend, bool *hasPrev= nullptr,bool *hasNext=nullptr);
            const ShapeList getMBRsByTime(id_type &id,double tstart,double tend);
            const ShapeList getMBCsByTime(id_type &id,double tstart,double tend);
            const Region getMBR(id_type &id);
            const MBC getMBC(id_type &id);
            DiskMultiMap m_dentries;
            std::map<id_type, Entry*> m_entries;//map from seg id to entry
            std::map<id_type, MBC> m_entryMbcs;
            std::map<id_type, Region> m_entryMbrs;
            std::map<id_type,id_type> m_part2node;
            IStorageManager* m_pStorageManager;
            SBBStream *m_stream= nullptr;
            std::string m_name;
            uint32_t m_pageSize;
            uint32_t m_maxTrajSegs=100;
            uint32_t m_leaf1=0,m_leaf2=0;
            uint32_t m_trajIO=0,m_indexIO=0;
            uint32_t m_loadedTraj=0;
            uint32_t m_boundingVisited=0;
            uint32_t m_segCount=0;
            double m_timeCount=0;
            double m_maxVelocity=-1;
            double m_IOtime=0;
        };
    }//namespace StorageManager



    class baseSegmentStream:public IDataStream{
    public:
        std::map<id_type, MBC>::iterator iter;
        std::map<id_type, Region> *m_brs;
        std::map<id_type, MBC> *m_bcs;
        baseSegmentStream(TrajStore *ts)
            :m_brs(&ts->m_entryMbrs),m_bcs(&ts->m_entryMbcs),iter(ts->m_entryMbcs.begin()){}
        virtual bool hasNext() override
        {
            return iter!=m_bcs->end();
        }
        virtual IData* constructData(id_type id,Region mbr,MBC mbc)=0;
        virtual IData* getNext() override{
            uint8_t *data;
            uint32_t len;
            auto d=constructData(iter->first,(*m_brs)[iter->first],iter->second);
            iter++;
            return d;
        }
        virtual uint32_t size()
        {
            return m_bcs->size();
        }

        virtual void rewind(){iter=m_bcs->begin();}
    };


    inline double triangleArea(double d,double v,double t){
        if(d<v*t){
            return d*d/v;
        }
        else{
            return 0.5*d*(d-v*t)*t;
        }
    }
}