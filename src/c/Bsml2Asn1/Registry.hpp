#ifndef REGISTRY_HPP
#define REGISTRY_HPP

#include <map>
#include <string>
#include <iostream>

class CRegistry
{

public:

    class CRegistryData
    {

    public:
        const static std::string sm_NotFound;

    public:
        typedef std::map<std::string, std::string> TData;
        typedef TData::const_iterator TConstIterator;
        typedef TConstIterator const_iterator;

    private:
        TData m_Data;

    public:
        const std::string& GetValue(const std::string &key) const;
        void AddValue(const std::string &key, const std::string &val);
        TConstIterator begin() const;
        TConstIterator end() const;
    };

public:
    typedef CRegistryData TRegistryData;
    typedef std::multimap<std::string, const TRegistryData *> TEntries;

private:
    typedef std::multimap<std::string, TRegistryData *> TData;

private:
    TData m_Data;

private:
    CRegistry(const CRegistry &);
    CRegistry& operator=(const CRegistry &);

public:
    CRegistry();
    CRegistry(std::istream &data);
    ~CRegistry();
    void ParseEntries(std::istream &data);
    void GetEntriesByName(const std::string &section, TEntries &entries) const;
    void GetAllEntries(TEntries &entries) const;
};

inline
CRegistry::CRegistry()
{
}

inline
CRegistry::CRegistry(std::istream &data)
{
    ParseEntries(data);
}

inline
CRegistry::~CRegistry()
{
    for (TData::iterator i = m_Data.begin(); i != m_Data.end(); ++i) {
        delete i->second;
    }
}

inline
void
CRegistry::GetEntriesByName(const std::string &section,
                            TEntries &entries) const
{
    entries.clear();
    for (TData::const_iterator i = m_Data.lower_bound(section);
         i != m_Data.upper_bound(section); ++i) {
        entries.insert(TEntries::value_type(i->first, i->second));
    }
}

inline
void
CRegistry::GetAllEntries(TEntries &entries) const
{
    entries.clear();
    for (TData::const_iterator i = m_Data.begin(); i != m_Data.end(); ++i) {
        entries.insert(TEntries::value_type(i->first, i->second));
    }
}

inline
const std::string&
CRegistry::CRegistryData::GetValue(const std::string &key) const
{
    TData::const_iterator i = m_Data.find(key);
    return i != m_Data.end() ? i->second : sm_NotFound;
}

inline
void
CRegistry::CRegistryData::AddValue(const std::string &key,
                                   const std::string &value)
{
    m_Data[key] = value;
}

inline
CRegistry::CRegistryData::TConstIterator
CRegistry::CRegistryData::begin() const
{
    return m_Data.begin();
}

inline
CRegistry::CRegistryData::TConstIterator
CRegistry::CRegistryData::end() const
{
    return m_Data.end();
}

#endif
