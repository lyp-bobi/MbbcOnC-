//
// Created by chuang on 4/4/19.
//


#pragma once

namespace SpatialIndex
{
    namespace R2Tree
    {
        class ExternalSorter
        {
        public:
            class Record
            {
            public:
                Record();
                Record(const Mbbc& r, id_type id, uint32_t len, uint8_t* pData, uint32_t s);
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
                Mbbc m_Mbbc;
                id_type m_id;
                uint32_t m_len;
                uint8_t* m_pData;
                uint32_t m_s;
            };


            ExternalSorter(uint32_t u32PageSize, uint32_t u32BufferPages);
            virtual ~ExternalSorter();

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
            class cmpR{
            public:
                int dim;
                cmpR(int d):dim(d){}
                bool operator()(Record* const r1, Record* const r2){
                    if(dim==1) return r1->m_Mbbc.m_smbr.m_pLow[0]<r2->m_Mbbc.m_smbr.m_pLow[0];
                    else if(dim==2) return r1->m_Mbbc.m_embr.m_pLow[0]<r2->m_Mbbc.m_embr.m_pLow[0];
                    else if(dim==3) return r1->m_Mbbc.m_smbr.m_pLow[1]<r2->m_Mbbc.m_smbr.m_pLow[1];
                    else if(dim==4) return r1->m_Mbbc.m_embr.m_pLow[1]<r2->m_Mbbc.m_embr.m_pLow[1];
                    else throw Tools::IllegalArgumentException("dimension what?");
                }
            };
            class cmpE{
            public:
                int dim;
                cmpE(int d):dim(d){}
                bool operator()(const PQEntry &e1,const PQEntry &e2){
                    if(dim==1) return e1.m_r->m_Mbbc.m_smbr.m_pLow[0]<e2.m_r->m_Mbbc.m_smbr.m_pLow[0];
                    else if(dim==2) return e1.m_r->m_Mbbc.m_embr.m_pLow[0]<e2.m_r->m_Mbbc.m_embr.m_pLow[0];
                    else if(dim==3) return e1.m_r->m_Mbbc.m_smbr.m_pLow[1]<e2.m_r->m_Mbbc.m_smbr.m_pLow[1];
                    else if(dim==4) return e1.m_r->m_Mbbc.m_embr.m_pLow[1]<e2.m_r->m_Mbbc.m_embr.m_pLow[1];
                    else throw Tools::IllegalArgumentException("dimension what?");
                }
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
        };

        class BulkLoader
        {
        public:
            void bulkLoadUsingSTR(
                    R2Tree* pTree,
                    IDataStream& stream,
                    uint32_t bindex,
                    uint32_t bleaf,
                    uint32_t pageSize, // The number of node entries per page.
                    uint32_t numberOfPages // The total number of pages to use.
            );
            void bulkLoadUsingSTR2(
                    R2Tree* pTree,
                    IDataStream& stream,
                    uint32_t bindex,
                    uint32_t bleaf,
                    uint32_t pageSize, // The number of node entries per page.
                    uint32_t numberOfPages // The total number of pages to use.
            );
            void bulkLoadUsingSTR3(
                    R2Tree* pTree,
                    IDataStream& stream,
                    uint32_t bindex,
                    uint32_t bleaf,
                    uint32_t pageSize, // The number of node entries per page.
                    uint32_t numberOfPages // The total number of pages to use.
            );

        protected:
            void createLevel(
                    R2Tree* pTree,
                    Tools::SmartPointer<ExternalSorter> es,
                    uint32_t dimension,
                    uint32_t indexSize,
                    uint32_t leafSize,
                    uint32_t level,
                    Tools::SmartPointer<ExternalSorter> es2,
                    uint32_t pageSize,
                    uint32_t numberOfPages
            );
            void createLevel2(
                    R2Tree* pTree,
                    Tools::SmartPointer<ExternalSorter> es,
                    uint32_t dimension,
                    uint32_t indexSize,
                    uint32_t leafSize,
                    uint32_t level,
                    Tools::SmartPointer<ExternalSorter> es2,
                    uint32_t pageSize,
                    uint32_t numberOfPages
            );

            Node* createNode(
                    R2Tree* pTree,
                    std::vector<ExternalSorter::Record*>& e,
                    uint32_t level
            );
        };
    }
}
