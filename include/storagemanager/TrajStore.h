//
// Created by Chuang on 2019/6/10.
//

#pragma once


#include <spatialindex/SpatialIndex.h>
#include <cstring>

using namespace SpatialIndex;
using namespace SpatialIndex::StorageManager;
using std::vector;

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
        double minx=1e300,maxx=-1e300,miny=1e300,maxy=-1e300,mint=1e300,maxt=-1e300;
        double Dx=0,Dy=0,Dt=0;
        double dist=0;
        static trajStat* instance();
        void clear(){
            M=lineCount=trajCount=jt=tl=Dx=Dy=Dt=0;
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
            std::fstream m_File;
            bool m_binput=true;
            bool m_busembc=false;
            id_type m_size;
            id_type m_cur;
            string m_s;
            SBBStream(std::string &str){
                m_File.open(str,std::ios::out);
                m_size=0;
                m_cur=0;
                m_s=str;
            }
            bool hasNext(){
                return m_cur<m_size;
            }
            void inputSBB(id_type id, Region* br, MBC* bc=nullptr){
                if(!m_binput){
                    m_File.open(m_s,std::ios::out);
                    m_binput=true;
                    rewind();
                }
                m_size+=1;
                m_File<<id<<std::endl<<br->toString()<<std::endl;
                if(bc!= nullptr){
                    m_busembc=true;
                    m_File<<bc->toString()<<std::endl;
                }
            }
            uint32_t size() override {
                return m_size;
            }
            void endInput(){
                m_File<<"END"<<std::endl;
                m_File.flush();
                m_File.close();
                m_File.clear();
                m_binput=false;
                m_File.open(m_s,std::ios::in);
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
                std::getline(m_File,s);
                id=std::stoll(s);
                std::getline(m_File,s);
                br.loadFromString(s);
                if(m_busembc){
                    std::getline(m_File,s);
                    bc.loadFromString(s);
                }
                if(m_size==m_cur){
                    std::cerr<<"SBB loading finished.\n";
                    m_File.close();
                    m_size=0;
                }
//                if(m_size-m_cur<10){
//                    std::cerr<<m_cur<<m_size<<std::endl<<id<<std::endl<<br<<std::endl<<bc<<std::endl;
//                }
                return new ircData(id,br,bc);
            }
            void rewind() override{
                m_cur=0;
                m_File.clear();
                m_File.seekg(std::ios::beg);
            }
            ~SBBStream(){
                m_File.close();
                std::remove(m_s.c_str());
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