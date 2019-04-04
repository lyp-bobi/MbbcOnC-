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
            IStorageManager* m_pStorageManager;

            id_type m_rootID, m_headerID;

        };//R2Tree
    }
}