//
// Created by chuang on 4/3/19.
//
#pragma once

namespace SpatialIndex {
    namespace R2Tree {
        class R2Tree;
        class Leaf;
        class Index;
        class Node;
        typedef Tools::PoolPointer<Node> NodePtr;

        class Node:public SpatialIndex::INode{
        public:
            virtual ~Node();

            //
            // Tools::IObject interface
            //
            virtual Tools::IObject* clone();

            //
            // Tools::ISerializable interface
            //
            virtual uint32_t getByteArraySize();
            virtual void loadFromByteArray(const byte* data);
            virtual void storeToByteArray(byte** data, uint32_t& len);

            //
            // SpatialIndex::IEntry interface
            //
            virtual id_type getIdentifier() const;
            virtual void getShape(IShape** out) const;

            //
            // SpatialIndex::INode interface
            //
            virtual uint32_t getChildrenCount() const;
            virtual id_type getChildIdentifier(uint32_t index)  const;
            virtual void getChildShape(uint32_t index, IShape** out)  const;
            virtual void getChildData(uint32_t index, uint32_t& length, byte** data) const;
            virtual uint32_t getLevel() const;
            virtual bool isIndex() const;
            virtual bool isLeaf() const;

        private:
            Node();
            Node(R2Tree* pTree, id_type id, uint32_t level, uint32_t capacity);

            virtual Node& operator=(const Node&);


            virtual void insertEntry(uint32_t dataLength, byte* pData, Mbbc& mbbc, id_type id);
            //virtual void deleteEntry(uint32_t index);

            //virtual bool insertData(uint32_t dataLength, byte* pData, Mbbc& mbbc, id_type id, std::stack<id_type>& pathBuffer, byte* overflowTable);
            //virtual void reinsertData(uint32_t dataLength, byte* pData, Mbbc& mbbc, id_type id, std::vector<uint32_t>& reinsert, std::vector<uint32_t>& keep);

            //virtual NodePtr chooseSubtree(const Mbbc& mbbc, uint32_t level, std::stack<id_type>& pathBuffer) = 0;
            virtual NodePtr findLeaf(const Mbbc& mbbc, id_type id, std::stack<id_type>& pathBuffer) = 0;

            R2Tree* m_pTree;
            // Parent of all nodes.

            uint32_t m_level;
            // The level of the node in the tree.
            // Leaves are always at level 0.

            id_type m_identifier;
            // The unique ID of this node.

            uint32_t m_children;
            // The number of children pointed by this node.

            uint32_t m_capacity;
            // Specifies the node capacity.

            Mbbc m_nodeMbbc;
            // The minimum bounding region enclosing all data contained in the node.

            byte** m_pData;
            // The data stored in the node.

            MbbcPtr* m_ptrMbbc;
            // The corresponding data MBBCs.

            id_type* m_pIdentifier;
            // The corresponding data identifiers.

            uint32_t* m_pDataLength;

            uint32_t m_totalDataLength;

            // Needed to access protected members without having to cast from Node.
            // It is more efficient than using member functions to access protected members.
            friend class R2Tree;
            friend class Leaf;
            friend class Index;
            friend class Tools::PointerPool<Node>;
            friend class BulkLoader;
        };
    }
}