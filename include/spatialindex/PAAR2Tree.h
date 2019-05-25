
#pragma once


namespace SpatialIndex
{
	namespace PAAR2Tree
	{
		SIDX_DLL enum BulkLoadMethod
		{
            BLM_STR = 0x1,
            BLM_STR2 = 0x2,
            BLM_STR3 = 0x3
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
			Data(uint32_t len, uint8_t* pData, MBBCk& r, id_type id);
			virtual ~Data();

			virtual Data* clone();
			virtual id_type getIdentifier() const;
			virtual void getShape(IShape** out) const;
			virtual void getData(uint32_t& len, uint8_t** data) const;
			virtual uint32_t getByteArraySize() const;
			virtual void loadFromByteArray(const uint8_t* data);
			virtual void storeToByteArray(uint8_t** data, uint32_t& len);

			id_type m_id;
			MBBCk m_MBBCk;
			uint8_t* m_pData;
			uint32_t m_dataLength;
		}; // Data

		SIDX_DLL ISpatialIndex* returnPAAR2Tree(IStorageManager& ind, Tools::PropertySet& in);
		SIDX_DLL ISpatialIndex* createNewPAAR2Tree(
			IStorageManager& sm,
			double fillFactor,
			uint32_t indexCapacity,
			uint32_t leafCapacity,
			uint32_t dimension,
			int k,
			id_type& indexIdentifier
		);
		SIDX_DLL ISpatialIndex* createAndBulkLoadNewPAAR2Tree(
			BulkLoadMethod m,
			IDataStream& stream,
			IStorageManager& sm,
			double fillFactor,
			uint32_t indexCapacity,
			uint32_t leafCapacity,
			uint32_t dimension,
			int k,
			id_type& indexIdentifier
		);
		SIDX_DLL ISpatialIndex* createAndBulkLoadNewPAAR2Tree(
			BulkLoadMethod m,
			IDataStream& stream,
			IStorageManager& sm,
			Tools::PropertySet& ps,
			id_type& indexIdentifier
		);
		SIDX_DLL ISpatialIndex* loadPAAR2Tree(IStorageManager& in, id_type indexIdentifier);
	}
}
