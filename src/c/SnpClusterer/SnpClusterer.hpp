#ifndef SNP_CLUSTERER_HPP
#define SNP_CLUSTERER_HPP

#include "Graph.hpp"
#include "SnpDataType.hpp"

#include <utility>
#include <string>
#include <vector>
#include <map>

class CSnpClusterer
{

public:
    typedef CSnpDataType* TSnp;
    typedef std::vector<TSnp> TSnps;
    typedef std::vector<TSnps> TClusters;
    typedef std::vector<std::string> TSnpData;

private:
    typedef CGraph<TSnp> TGraph;
    typedef std::pair<int, int> TInterval;
    typedef std::map<std::pair<std::string, TInterval>,
                     CSnpDataType *> TId2Snp;
    typedef std::map<TSnp, TSnp> TIndels;

private:
    class CSnpComparator
    {
    public:
        bool operator()(const TGraph::CNode *node1,
                        const TGraph::CNode *node2);
    };

private:
    TGraph m_Graph;
    TId2Snp m_Id2Snp;
    bool m_DeleteData;
    TIndels m_Indels;

private:
    CSnpClusterer(const CSnpClusterer &);
    CSnpClusterer& operator=(const CSnpClusterer &);
    CSnpDataType* x_GetSnp(const std::string &, bool &);
    void x_ComplementSnp(CSnpDataType *);
    void x_AddToDecomposedCluster(TClusters &, CSnpDataType *);

public:
    CSnpClusterer(bool deleteData = true);
    ~CSnpClusterer();
    void AddSnp(const TSnp &snp1, const TSnp &snp2);
    void AddIndel(const TSnp &snp1, const TSnp &snp2);
    void AddSnps(const TSnpData &snpData);
    void GetClusters(TClusters &clusters);
};

inline
CSnpClusterer::CSnpClusterer(bool deleteData) : m_DeleteData(deleteData)
{
}

inline
CSnpClusterer::~CSnpClusterer()
{
    if (m_DeleteData) {
        for (TId2Snp::iterator i = m_Id2Snp.begin(); i != m_Id2Snp.end();
             ++i) {
            delete i->second;
        }
    }
}

inline
void
CSnpClusterer::AddSnp(const TSnp &snp1, const TSnp &snp2)
{
    m_Graph.AddData(snp1, snp2);
}

inline
void
CSnpClusterer::AddIndel(const TSnp &snp1, const TSnp &snp2)
{
    m_Indels.insert(TIndels::value_type(snp1, snp2));
}

inline
void
CSnpClusterer::x_ComplementSnp(CSnpDataType *snp)
{
    snp->SetPlus(!snp->IsPlus());
}

inline
bool
CSnpClusterer::CSnpComparator::operator()(const TGraph::CNode *node1,
                                          const TGraph::CNode *node2)
{
    const TSnp &snp1 = node1->GetData();
    const TSnp &snp2 = node2->GetData();
    if (snp1->GetId() == snp2->GetId()) {
        return snp1->GetFrom() < snp2->GetFrom();
    }
    return snp1->GetId() < snp2->GetId();
}

#endif
