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


            class SIDX_DLL Data : public IData, public Tools::ISerializable
            {
            public:
                Data(uint32_t len, byte* pData, Region& r, id_type id);
                virtual ~Data();

                virtual Data* clone();
                virtual id_type getIdentifier() const;
                virtual void getShape(IShape** out) const;
                virtual void getData(uint32_t& len, byte** data) const;
                virtual uint32_t getByteArraySize();
                virtual void loadFromByteArray(const byte* data);
                virtual void storeToByteArray(byte** data, uint32_t& len);

                id_type m_id;
                Region m_region;
                byte* m_pData;
                uint32_t m_dataLength;
            }; // Data
        };
    }
}