#ifndef SNP_DATA_TYPE_HPP
#define SNP_DATA_TYPE_HPP

#include <string>
#include <ostream>

class CSnpDataType
{

private:
    std::string m_Id;
    int m_From;
    int m_To;
    bool m_Plus;

public:
    CSnpDataType(const std::string &data);
    CSnpDataType(const std::string &id, int from, int to, bool plus);
    const std::string& GetId() const;
    int GetFrom() const;
    int GetTo() const;
    bool IsPlus() const;
    void SetPlus(bool plus);
    int GetLength() const;
    friend std::ostream& operator<<(std::ostream &os, const CSnpDataType &snp);

private:
    void x_SetData(const std::string &id, int from, int to, bool plus);
};

inline
CSnpDataType::CSnpDataType(const std::string &id, int from, int to, bool plus)
{
    x_SetData(id, from, to, plus);
}

inline
const std::string&
CSnpDataType::GetId() const
{
    return m_Id;
}

inline
int
CSnpDataType::GetFrom() const
{
    return m_From;
}

inline
int
CSnpDataType::GetTo() const
{
    return m_To;
}

inline
bool
CSnpDataType::IsPlus() const
{
    return m_Plus;
}

inline
void
CSnpDataType::SetPlus(bool plus)
{
    m_Plus = plus;
}

inline
int
CSnpDataType::GetLength() const
{
    return GetTo() - GetFrom();
}

inline
void
CSnpDataType::x_SetData(const std::string &id, int from, int to, bool plus)
{
    m_Id = id;
    m_From = from;
    m_To = to;
    SetPlus(plus);
}

inline
std::ostream&
operator<<(std::ostream &os, const CSnpDataType &snp)
{
    os << snp.GetId() << "\t" << snp.GetFrom() << "\t" << snp.GetTo() << "\t"
       << (snp.IsPlus() ? "+" : "-");
    return os;
}

#endif
