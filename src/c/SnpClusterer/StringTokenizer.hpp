#ifndef STRING_TOKENIZER_HPP
#define STRING_TOKENIZER_HPP

#include <string>
#include <cstring>
#include <stdexcept>
#include <vector>

class CStringTokenizer
{

public:
    typedef std::vector<std::string> TTokens;

public:
    class CIterator
    {

    private:
        friend class CStringTokenizer;

    private:
        CStringTokenizer *m_Tokenizer;

    private:
        CIterator(CStringTokenizer *);

    public:
        CIterator(const CIterator &);
        CIterator& operator=(const CIterator &);
        const char* operator*() const;
        const CIterator& operator++();
        operator bool() const;
    };
    typedef CIterator const_iterator;

private:
    char *m_OrigStr;
    char *m_TokenizedStr;
    char *m_Token;
    char *m_Pos;
    char *m_Delim;

private:
    CStringTokenizer(const CStringTokenizer &);
    CStringTokenizer& operator=(const CStringTokenizer &);

public:
    CStringTokenizer(const std::string &data, const std::string &delim);
    ~CStringTokenizer();
    bool HasMoreTokens() const;
    const char* GetToken() const;
    const char* NextToken();
    void Reset();
    const_iterator begin();
    TTokens& GetTokens(TTokens &);
};

inline
CStringTokenizer::CStringTokenizer(const std::string &data,
                                   const std::string &delim)
{
    m_OrigStr = strdup(data.c_str());
    m_TokenizedStr = 0;
    m_Delim = strdup(delim.c_str());
    if (!m_OrigStr || !m_Delim) {
        if (m_OrigStr) {
            free(m_OrigStr);
        }
        if (m_Delim) {
            free(m_Delim);
        }
        throw std::runtime_error("Out of memory");
    }
    Reset();
}

inline
CStringTokenizer::~CStringTokenizer()
{
    if (m_OrigStr) {
        free(m_OrigStr);
    }
    if (m_TokenizedStr) {
        free(m_TokenizedStr);
    }
    if (m_Delim) {
        free(m_Delim);
    }
}

inline
bool
CStringTokenizer::HasMoreTokens() const
{
    return m_Token;
}

inline
const char*
CStringTokenizer::GetToken() const
{
    return m_Token;
}

inline
const char*
CStringTokenizer::NextToken()
{
    m_Token = strtok_r(0, m_Delim, &m_Pos);
    return m_Token;
}

inline
void
CStringTokenizer::Reset()
{
    if (m_TokenizedStr) {
        free(m_TokenizedStr);
    }
    m_TokenizedStr = strdup(m_OrigStr);
    if (!m_TokenizedStr) {
        throw std::runtime_error("Out of memory");
    }
    m_Token = strtok_r(m_TokenizedStr, m_Delim, &m_Pos);
}

inline
CStringTokenizer::const_iterator
CStringTokenizer::begin()
{
    Reset();
    return const_iterator(this);
}

inline
CStringTokenizer::CIterator::CIterator(CStringTokenizer *tokenizer) :
    m_Tokenizer(tokenizer)
{
}

inline
CStringTokenizer::CIterator::CIterator(const CIterator &copy)
{
    m_Tokenizer = copy.m_Tokenizer;
}

inline
CStringTokenizer::CIterator&
CStringTokenizer::CIterator::operator=(const CIterator &copy)
{
    if (this != &copy) {
        m_Tokenizer = copy.m_Tokenizer;
    }
    return *this;
}

inline
const char*
CStringTokenizer::CIterator::operator*() const
{
    return m_Tokenizer->GetToken();
}

inline
const CStringTokenizer::CIterator&
CStringTokenizer::CIterator::operator++()
{
    m_Tokenizer->NextToken();
    return *this;
}

inline
CStringTokenizer::CIterator::operator bool() const
{
    return m_Tokenizer->HasMoreTokens();
}

inline
CStringTokenizer::TTokens&
CStringTokenizer::GetTokens(TTokens &tokens)
{
    while (HasMoreTokens()) {
        tokens.push_back(GetToken());
        NextToken();
    }
    return tokens;
}

#endif
