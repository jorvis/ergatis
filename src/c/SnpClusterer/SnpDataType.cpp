#include "SnpDataType.hpp"

#include <StringTokenizer.hpp>

using namespace std;

CSnpDataType::CSnpDataType(const string &data)
{
    CStringTokenizer stoken(data, "\t");
    CStringTokenizer::TTokens tokens;
    stoken.GetTokens(tokens);
    x_SetData(tokens[0], atoi(tokens[1].c_str()), atoi(tokens[2].c_str()),
            tokens[3] == "+");
}
