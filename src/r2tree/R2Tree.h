//
// Created by chuang on 4/3/19.
//

#pragma once

#include "Statistics.h"
#include "Node.h"
#include "PointerPoolNode.h"

namespace SpatialIndex
{
    namespace R2Tree
    {
        class R2Tree : public ISpatialIndex
        {
        public:
            R2Tree(IStorageManager&, Tools::PropertySet&);
            // String                   Value     Description
            // ----------------------------------------------
            // IndexIndentifier         VT_LONG   If specified an existing index will be openened from the supplied
            //                          storage manager with the given index id. Behaviour is unspecified
            //                          if the index id or the storage manager are incorrect.
            // Dimension                VT_ULONG  Dimensionality of the data that will be inserted.
            // IndexCapacity            VT_ULONG  The index node capacity. Default is 100.
            // LeafCapactiy             VT_ULONG  The leaf node capacity. Default is 100.
            // EnsureTightMBRs          VT_BOOL   Default is true
            // IndexPoolCapacity        VT_LONG   Default is 100
            // LeafPoolCapacity         VT_LONG   Default is 100
            // RegionPoolCapacity       VT_LONG   Default is 1000
            // PointPoolCapacity        VT_LONG   Default is 500

            virtual ~R2Tree();

            //
            // ISpatialIndex interface
            //
            virtual void insertData(uint32_t len, const byte* pData, const IShape& shape, id_type shapeIdentifier);
            virtual bool deleteData(const IShape& shape, id_type id);
            virtual void containsWhatQuery(const IShape& query, IVisitor& v);
            virtual void intersectsWithQuery(const IShape& query, IVisitor& v);
            virtual void pointLocationQuery(const Point& query, IVisitor& v);
            virtual void nearestNeighborQuery(uint32_t k, const IShape& query, IVisitor& v, INearestNeighborComparator&);
            virtual void nearestNeighborQuery(uint32_t k, const IShape& query, IVisitor& v);
            virtual void selfJoinQuery(const IShape& s, IVisitor& v);
            virtual void queryStrategy(IQueryStrategy& qs);
            virtual void getIndexProperties(Tools::PropertySet& out) const;
            virtual void addCommand(ICommand* pCommand, CommandType ct);
            virtual bool isIndexValid();
            virtual void getStatistics(IStatistics** out) const;




        private:
            void initNew(Tools::PropertySet&);
            void initOld(Tools::PropertySet& ps);
            void storeHeader();
            void loadHeader();

            //void insertData_impl(uint32_t dataLength, byte* pData, Mbbc& mbbc, id_type id);
            //void insertData_impl(uint32_t dataLength, byte* pData, Mbbc& mbbc, id_type id, uint32_t level, byte* overflowTable);
            //bool deleteData_impl(const Region& mbr, id_type id);

            id_type writeNode(Node*);
            NodePtr readNode(id_type page);
            void deleteNode(Node*);

            void rangeQuery(RangeQueryType type, const IShape& query, IVisitor& v);

            IStorageManager* m_pStorageManager;

            id_type m_rootID, m_headerID;


            double m_fillFactor;

            uint32_t m_indexCapacity;

            uint32_t m_leafCapacity;

            uint32_t m_dimension;

            Mbbc m_infiniteMbbc;

            Statistics m_stats;

            bool m_bTightMBRs;


            Tools::PointerPool<Point> m_pointPool;
            Tools::PointerPool<Mbbc> m_MbbcPool;
            Tools::PointerPool<Node> m_indexPool;
            Tools::PointerPool<Node> m_leafPool;

#ifdef HAVE_PTHREAD_H
            pthread_mutex_t m_lock;
#endif

            friend class Node;
            friend class Leaf;
            friend class Index;
            friend class BulkLoader;

            friend std::ostream& operator<<(std::ostream& os, const R2Tree& t);
        };//R2Tree
    }
}
