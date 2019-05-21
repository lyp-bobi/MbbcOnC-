//
// Created by chuang on 4/8/19.
//


#pragma once


namespace SpatialIndex
{
    namespace PAAR2Tree
    {
        class Leaf : public Node
        {
        public:
            virtual ~Leaf();

        protected:
            Leaf(PAAR2Tree* pTree, id_type id);

            //virtual NodePtr chooseSubtree(const MBBCk& mbbc, uint32_t level, std::stack<id_type>& pathBuffer);
            virtual NodePtr findLeaf(const MBBCk& mbbc, id_type id, std::stack<id_type>& pathBuffer);

            //virtual void split(uint32_t dataLength, uint8_t* pData, MBBCk& mbbc, id_type id, NodePtr& left, NodePtr& right);

            //virtual void deleteData(id_type id, std::stack<id_type>& pathBuffer);

            friend class PAAR2Tree;
            friend class BulkLoader;
        }; // Leaf
    }
}
