#include "Registry.hpp"

using namespace std;

const string CRegistry::CRegistryData::sm_NotFound = "";

void
CRegistry::ParseEntries(istream &data)
{
    string line;
    TRegistryData *reg = 0;
    while (getline(data, line)) {
        string::size_type idx;
        if (line[0] == ';') {
            continue;
        }
        if (line[0] == '[' && line[line.length() - 1] == ']') {
            string section = line.substr(1, line.length() - 2);
            reg = new TRegistryData;
            m_Data.insert(TData::value_type(section, reg));
        }
        else if ((idx = line.find_first_of('=')) != string::npos) {
            string key = line.substr(0, idx);
            string val = line.substr(idx + 1);
            reg->AddValue(key, val);
        }
    }
}
