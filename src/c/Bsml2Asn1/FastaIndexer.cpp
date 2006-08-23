#include "FastaIndexer.hpp"

#include <iostream>
#include <stdexcept>

using namespace std;

CFastaIndexer::CFastaIndexer(const string &fileName) :
    m_Data(fileName.c_str())
{
    if (!m_Data.is_open()) {
        throw runtime_error("Error opening FastA data: " + fileName);
    }
    x_Init();
}

string&
CFastaIndexer::GetSequence(const TId &id, string &seq)
{
    seq.clear();
    TOffsets::const_iterator offset = m_Offsets.find(id);
    if (offset != m_Offsets.end()) {
        m_Data.clear();
        m_Data.seekg(offset->second);
        string buf;
        while (getline(m_Data, buf)) {
            if (buf[0] == '>') {
                break;
            }
            seq += buf;
        }
    }
    return seq;
}

void
CFastaIndexer::x_Init()
{
    string buf;
    while (getline(m_Data, buf)) {
        if (buf[0] == '>') {
            m_Offsets.insert(TOffsets::value_type(buf.substr(1),
                                                  m_Data.tellg()));
        }
    }
}
