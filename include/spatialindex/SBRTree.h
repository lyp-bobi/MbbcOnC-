
#pragma once

namespace SpatialIndex
{
	namespace SBRTree
	{
		SIDX_DLL enum SBRTreeVariant
		{
			RV_LINEAR = 0x0,
			RV_QUADRATIC,
			RV_RSTAR
		};

		SIDX_DLL enum BulkLoadMethod
		{
			BLM_STR = 0x0
		};

		SIDX_DLL enum PersistenObjectIdentifier
		{
			PersistentIndex = 0x1,
			PersistentLeaf = 0x2
		};

		SIDX_DLL enum RangeQueryType
		{
			ContainmentQuery = 0x1,
			IntersectionQuery = 0x2
		};

		class SIDX_DLL Data : public IData, public Tools::ISerializable
		{
		public:
            Data(uint32_t len, uint8_t* pData, MBC& r, id_type id);
			Data(uint32_t len, uint8_t* pData, MBC& r,SBR &rg, id_type id);
			virtual ~Data();

			virtual Data* clone();
			virtual id_type getIdentifier() const;
			virtual void getShape(IShape** out) const;
			virtual void getData(uint32_t& len, uint8_t** data) const;
			virtual uint32_t getByteArraySize() const;
			virtual void loadFromByteArray(const uint8_t* data);
			virtual void storeToByteArray(uint8_t** data, uint32_t& len);

			id_type m_id;
			SBR m_sbr;
			MBC m_mbc;
			uint8_t* m_pData;
			uint32_t m_dataLength;
		}; // Data

		SIDX_DLL ISpatialIndex* returnSBRTree(IStorageManager& ind, Tools::PropertySet& in);
		SIDX_DLL ISpatialIndex* createNewSBRTree(
			IStorageManager& sm,
			double fillFactor,
			uint32_t indexCapacity,
			uint32_t leafCapacity,
			uint32_t dimension,
			SBRTreeVariant rv,
			id_type& indexIdentifier
		);
		SIDX_DLL ISpatialIndex* createAndBulkLoadNewSBRTree(
			BulkLoadMethod m,
			IDataStream& stream,
			IStorageManager& sm,
			double fillFactor,
			uint32_t indexCapacity,
			uint32_t leafCapacity,
			uint32_t dimension,
			SBRTreeVariant rv,
			id_type& indexIdentifier
		);
		SIDX_DLL ISpatialIndex* createAndBulkLoadNewSBRTree(
			BulkLoadMethod m,
			IDataStream& stream,
			IStorageManager& sm,
			Tools::PropertySet& ps,
			id_type& indexIdentifier
		);
		SIDX_DLL ISpatialIndex* loadSBRTree(IStorageManager& in, id_type indexIdentifier);
        SIDX_DLL ISpatialIndex* createAndBulkLoadNewSBRTreeWithTrajStore(
                IStorageManager *tsm,
                uint32_t indexCapacity,
                uint32_t dimension,
                id_type& indexIdentifier
        );
	}
}
