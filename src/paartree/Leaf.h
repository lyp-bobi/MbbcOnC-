//
// Created by chuang on 4/8/19.
//


#pragma once


namespace SpatialIndex
{
    namespace PAARTree
    {
        class Leaf : public Node
        {
        public:
            virtual ~Leaf();

        protected:
            Leaf(PAARTree* pTree, id_type id);

            //virtual NodePtr chooseSubtree(const MBRk& mbbc, uint32_t level, std::stack<id_type>& pathBuffer);
            virtual NodePtr findLeaf(const MBRk& mbbc, id_type id, std::stack<id_type>& pathBuffer);

            //virtual void split(uint32_t dataLength, uint8_t* pData, MBRk& mbbc, id_type id, NodePtr& left, NodePtr& right);

            //virtual void deleteData(id_type id, std::stack<id_type>& pathBuffer);

            friend class PAARTree;
            friend class BulkLoader;
        }; // Leaf
    }
}
