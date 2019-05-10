//
// Created by chuang on 4/8/19.
//


#pragma once


namespace SpatialIndex
{
    namespace R2Tree
    {
        class Leaf : public Node
        {
        public:
            virtual ~Leaf();

        protected:
            Leaf(R2Tree* pTree, id_type id);

            //virtual NodePtr chooseSubtree(const Mbbc& mbbc, uint32_t level, std::stack<id_type>& pathBuffer);
            virtual NodePtr findLeaf(const Mbbc& mbbc, id_type id, std::stack<id_type>& pathBuffer);

            //virtual void split(uint32_t dataLength, uint8_t* pData, Mbbc& mbbc, id_type id, NodePtr& left, NodePtr& right);

            //virtual void deleteData(id_type id, std::stack<id_type>& pathBuffer);

            friend class R2Tree;
            friend class BulkLoader;
        }; // Leaf
    }
}
