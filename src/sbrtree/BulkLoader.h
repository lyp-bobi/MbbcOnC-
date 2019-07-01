
#pragma once

namespace SpatialIndex
{
	namespace SBRTree
	{
		class ExternalSorter
		{
		public:
			class Record
			{
			public:
				Record();
				Record(const SBR& shape, id_type id, uint32_t len, uint8_t* pData, uint32_t s,uint32_t level,MBC* mbc= nullptr);
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
				SBR m_sbr;
				MBC m_mbc;
				id_type m_id;
				uint32_t m_len;
				uint8_t* m_pData;
				uint32_t m_s;
			};

		public:
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
				PQEntry(Record* r, uint32_t u32Index) : m_sbr(r), m_u32Index(u32Index) {}

				struct SortAscending : public std::binary_function<const PQEntry&, const PQEntry&, bool>
				{
					bool operator()(const PQEntry& e1, const PQEntry& e2)
					{
						if (*(e1.m_sbr) < *(e2.m_sbr)) return true;
						else return false;
					}
				};

				Record* m_sbr;
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
		};

		class BulkLoader
		{
		public:
			void bulkLoadUsingSTR(
				SBRTree* pTree,
				IDataStream& stream,
				uint32_t bindex,
				uint32_t bleaf,
				uint32_t pageSize, // The number of node entries per page.
				uint32_t numberOfPages // The total number of pages to use.
			);

		protected:
			void createLevel(
				SBRTree* pTree,
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
				SBRTree* pTree,
				std::vector<ExternalSorter::Record*>& e,
				uint32_t level
			);
		};
	}
}
