#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <map>
#include <utility>

template<typename T>
class CGraph
{

public:

    typedef T TType;

    class CEdge;
    class CNode
    {

    public:
        typedef std::vector<const CEdge *> TEdges;

    private:
        T m_Data;
        TEdges m_Edges;

    private:
        CNode(const CNode &);
        CNode& operator=(const CNode &);

    public:
        CNode(const T &data);
        void AddEdge(const CEdge *edge);
        const T& GetData() const;
        const TEdges& GetEdges() const;
    };

    typedef std::pair<CNode *, CNode *> TEdgeBase;
    class CEdge : public TEdgeBase
    {

    private:
        CEdge(const CEdge &);
        CEdge& operator=(const CEdge &);

    public:
        CEdge(CNode *node1, CNode *node2);
        const CNode* GetFirstNode() const;
        const CNode* GetSecondNode() const;
    };

public:
    enum EGraphType
    {
        eGraphType_Undirected = 0,
        eGraphType_Directed
    };

public:
    typedef std::map<T, CNode *> TNodes;
    typedef std::vector<CEdge *> TEdges;
    typedef std::vector<const CNode *> TNodePtrs;

private:
    EGraphType m_GraphType;
    TNodes m_Nodes;
    TEdges m_Edges;

private:
    CGraph(const CGraph &);
    CGraph& operator=(const CGraph &);

public:
    CGraph(EGraphType graphType = eGraphType_Undirected);
    virtual ~CGraph();
    void AddData(const T &data1, const T &data2);
    void GetNodePtrs(TNodePtrs &nodePtrs) const;
    const CNode* GetNodePtr(const T &data) const;
};

template<typename T>
inline
CGraph<T>::CGraph(EGraphType graphType) : m_GraphType(graphType)
{
}

template<typename T>
inline
CGraph<T>::~CGraph()
{
    for (typename TNodes::const_iterator i = m_Nodes.begin();
         i != m_Nodes.end(); ++i) {
        delete i->second;
    }
    for (typename TEdges::const_iterator i = m_Edges.begin();
         i != m_Edges.end(); ++i) {
        delete *i;
    }
}

template<typename T>
inline
void
CGraph<T>::AddData(const T &data1, const T &data2)
{
    CNode **node1 = &m_Nodes[data1];
    CNode **node2 = &m_Nodes[data2];
    if (!*node1) {
        *node1 = new CNode(data1);
    }
    if (!*node2) {
        *node2 = new CNode(data2);
    }
    CEdge *edge = new CEdge(*node1, *node2);
    m_Edges.push_back(edge);
    (*node1)->AddEdge(edge);
    if (m_GraphType == eGraphType_Undirected) {
        (*node2)->AddEdge(edge);
    }
}

template<typename T>
inline
void
CGraph<T>::GetNodePtrs(TNodePtrs &nodePtrs) const
{
    nodePtrs.clear();
    for (typename TNodes::const_iterator i = m_Nodes.begin();
         i != m_Nodes.end(); ++i) {
        nodePtrs.push_back(i->second);
    }
}

template<typename T>
inline
const typename CGraph<T>::CNode*
CGraph<T>::GetNodePtr(const T &data) const
{
    typename TNodes::const_iterator i = m_Nodes.find(data);
    return i != m_Nodes.end() ? i->second : 0;
}

template<typename T>
inline
CGraph<T>::CNode::CNode(const T &data) : m_Data(data)
{
}

template<typename T>
inline
void
CGraph<T>::CNode::AddEdge(const CEdge *edge)
{
    m_Edges.push_back(edge);
}

template<typename T>
inline
const T&
CGraph<T>::CNode::GetData() const
{
    return m_Data;
}

template<typename T>
inline
const typename CGraph<T>::CNode::TEdges&
CGraph<T>::CNode::GetEdges() const
{
    return m_Edges;
}

template<typename T>
inline
CGraph<T>::CEdge::CEdge(CNode *node1, CNode *node2) :
    CGraph<T>::TEdgeBase(node1, node2)
{
}

template<typename T>
inline
const typename CGraph<T>::CNode*
CGraph<T>::CEdge::GetFirstNode() const
{
    return this->first;
}

template<typename T>
inline
const typename CGraph<T>::CNode*
CGraph<T>::CEdge::GetSecondNode() const
{
    return this->second;
}

#endif
