#ifndef FASTA_INDERXER_HPP
#define FASTA_INDERXER_HPP

#include <map>
#include <string>
#include <fstream>

class CFastaIndexer
{

public:
    typedef std::string TId;

private:
    typedef std::map<TId, std::istream::pos_type> TOffsets;

private:
    TOffsets m_Offsets;
    std::ifstream m_Data;

private:
    void x_Init();
    CFastaIndexer(const CFastaIndexer &);
    CFastaIndexer& operator=(const CFastaIndexer &);

public:
    CFastaIndexer(const std::string &fileName);
    std::string& GetSequence(const TId &id, std::string &seq);
};

#endif
