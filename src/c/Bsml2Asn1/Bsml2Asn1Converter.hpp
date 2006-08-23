#ifndef BSML_2_ASN1_CONVERTER_HPP
#define BSML_2_ASN1_CONVERTER_HPP

#include "Registry.hpp"

#include <string>
#include <vector>
#include <map>
#include <set>
#include <xalanc/Include/PlatformDefinitions.hpp>
#include <corelib/ncbistl.hpp>
#include <objects/seq/Seq_annot.hpp>
#include <objects/seq/MolInfo.hpp>
#include <objects/seqloc/Packed_seqint.hpp>
#include <objmgr/object_manager.hpp>

XALAN_CPP_NAMESPACE_BEGIN

class XalanSourceTreeInit;
class XalanSourceTreeDOMSupport;
class XalanSourceTreeParserLiaison;
class XPathEvaluator;
class XalanDocumentPrefixResolver;
class XalanNode;
class XalanElement;
class NodeRefList;
class XalanDOMString;

XALAN_CPP_NAMESPACE_END

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(objects)

class CSeq_loc;
class CSeq_feat;
class CSeq_id;
class CBioseq;
class CSeqdesc;
class CContact_info;
class CCit_sub;
class CAuthor;
class CAffil;
class CPerson_id;
class CTaxon1;
class CObjectManager;
class CSeq_entry;
class CScope;

END_SCOPE(objects)
END_NCBI_SCOPE

class CFastaIndexer;

class CBsml2Asn1Converter
{

public:
    CBsml2Asn1Converter();
    ~CBsml2Asn1Converter();
    void Convert(const std::string &bsmlFile, const std::string &asn1File,
                 std::istream *contact = 0, bool update = false,
                 bool genProdSet = false);

private:
    typedef XALAN_CPP_NAMESPACE_QUALIFIER XalanNode TNode;
    typedef XALAN_CPP_NAMESPACE_QUALIFIER XalanElement TElement;
    typedef XALAN_CPP_NAMESPACE_QUALIFIER NodeRefList TNodeList;
    typedef XALAN_CPP_NAMESPACE_QUALIFIER XalanDocumentPrefixResolver
            TResolver;
    typedef XALAN_CPP_NAMESPACE_QUALIFIER XalanDOMString TXalanString;
    typedef ncbi::objects::CSeq_loc TSeqLoc;
    typedef ncbi::objects::CSeq_feat TSeqFeat;
    typedef ncbi::objects::CSeq_id TSeqId;
    typedef ncbi::objects::CSeq_annot::TData::TFtable TFeats;
    typedef ncbi::objects::CBioseq TBioseq;
    typedef ncbi::objects::CSeq_entry TEntry;
    typedef ncbi::objects::CSeqdesc TSeqDesc;
    typedef ncbi::objects::CContact_info TContactInfo;
    typedef ncbi::objects::CCit_sub TCitation;
    typedef ncbi::objects::CTaxon1 TTaxon;
    typedef ncbi::objects::CAuthor TAuthor;
    typedef ncbi::objects::CAffil TAffil;
    typedef ncbi::objects::CMolInfo TMolInfo;
    typedef ncbi::objects::CPacked_seqint TPackedInt;
    typedef ncbi::objects::CPerson_id TName;
    typedef ncbi::objects::CObjectManager TObjectManager;
    typedef ncbi::objects::CScope TScope;
    typedef std::vector<TNode *> TNodes;
    typedef std::map<std::string, CFastaIndexer *> TIndexes;
    typedef std::map<std::string, std::string> TTitles;
    typedef std::map<std::string, std::string> TIds;
    typedef std::map<std::string, TNode *> TFeatures;
    typedef std::map<std::string, TMolInfo::ECompleteness> TCompleteness;
    typedef std::set<std::string> TWrittenIds;
    typedef std::set<std::string> TClassStrings;
    typedef CRegistry TRegistry;
    typedef TRegistry::TRegistryData TRegistryData;

private:
    struct SFeatGroup
    {
        TSeqFeat *m_Gene;
        TSeqFeat *m_Rna;
        TSeqFeat *m_Cds;
    };

private:
    typedef std::vector<SFeatGroup> TFeatGroups;

private:
    class CSeqIntComparator
    {

    private:
        typedef TPackedInt::Tdata::value_type TData;

    private:
        bool m_Forward;

    public:
        CSeqIntComparator(bool forward = true);
        bool operator()(const TData &int1, const TData &int2);
    };

private:
    ncbi::CRef<TObjectManager> m_ObjMgr;
    TClassStrings m_ValidClasses;

private:
    CBsml2Asn1Converter(const CBsml2Asn1Converter &);
    CBsml2Asn1Converter& operator=(const CBsml2Asn1Converter &);
    void x_Init();
    TBioseq* x_ProcessSequence(TNode *, TIndexes &, TTitles &, TIds &,
                               TCompleteness &,
                               int, TFeatures &, bool, bool, TWrittenIds &,
                               TWrittenIds &, TFeatGroups &, TResolver &);
    void x_ProcessFeatureGroup(TNode *, TNode *, TFeats &, TTitles &,
                               TIds &,
                               TCompleteness &, int,
                               TFeatures &, TWrittenIds &, TWrittenIds &,
                               TFeatGroups &, bool, TResolver &);
    void x_ProcessMiscFeature(TNode *, TNode *, TFeats &, TResolver &);
    TNodeList& x_GetNodeList(TNodeList &, TNode *,
                             const std::string &, TResolver &);
    TNode* x_GetNode(TNode *, const std::string &, TResolver &);
    TElement* x_GetElement(TNode *, const std::string &, TResolver &);
    TSeqFeat* x_CreateGeneFeat(TNode *, TNode *, TResolver &);
    TSeqFeat *x_CreateRnaFeat(TNode *, TNode *, TNodes &, TIds &, bool,
                              TResolver &);
    TSeqFeat* x_CreateCdsFeat(TNode *, TNode *, TNode *, TNode *, TNodes &,
                              int, TIds &, TResolver &);
    TSeqLoc* x_CreateSeqLoc(TNode *, TNode *);
    TSeqId* x_CreateSeqId(const std::string &);
    TSeqDesc* x_ProcessGenome(TNode *, TResolver &);
    TContactInfo* x_ProcessContact(TNode *, TResolver &);
    TCitation* x_ProcessCitation(TNode *, TResolver &);
    TSeqDesc* x_ProcessPublication(TNode *, TResolver &);
    TAuthor* x_CreateAuthor(TNode *, TResolver &, TNode * = 0);
    TAffil* x_CreateAffil(TNode *, TResolver &);
    std::string x_ConvertToStdStr(const TXalanString &);
    void x_CleanupIndexes(TIndexes &);
    bool x_IsCanonical(TNode *, TResolver &);
    void x_SortPackedInt(TPackedInt &);

    TContactInfo* x_ProcessContact(TRegistry &);
    TAffil* x_CreateAffil(const TRegistryData &);
    TName* x_CreateName(const TRegistryData &);
    TCitation* x_ProcessCitation(const TRegistry &);

    TEntry* x_CreateRnaProtSet(SFeatGroup &, TScope &);
    TEntry* x_CreateProteinSequence(TSeqFeat &, TSeqFeat &, TScope &);
    TEntry* x_CreateRnaSequence(TSeqFeat &, TSeqFeat &, TScope &);
    TEntry* x_CreateSequence(const std::string &, TSeqFeat &, bool);
    void x_SetTitles(TEntry &);
    void x_InitTaxon();

private:
    XALAN_CPP_NAMESPACE_QUALIFIER XalanSourceTreeInit *m_SrcTreeInit;
    XALAN_CPP_NAMESPACE_QUALIFIER XalanSourceTreeDOMSupport *m_DomSupport;
    XALAN_CPP_NAMESPACE_QUALIFIER XalanSourceTreeParserLiaison *m_Liaison;
    XALAN_CPP_NAMESPACE_QUALIFIER XPathEvaluator *m_XPathEval;
    TTaxon *m_Taxon;
};

inline
CBsml2Asn1Converter::CSeqIntComparator::CSeqIntComparator(bool forward) :
    m_Forward(forward)
{
}

#endif
