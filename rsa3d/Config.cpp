//--------------------------------------------------------------------------------------------
// Simple configuration file parser
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "Config.h"
#include "Utils.h"


// Parses configuration parameters from given ifstream
//--------------------------------------------------------------------------------------------
Config::ptr Config::parse(std::ifstream & _file, char _delim)
{
    std::string line;
    std::string field, value;
    std::size_t line_num = 0;
    std::size_t pos;
    
    auto result = std::make_unique<Config>();
    
    // Parse line by line
    while (std::getline(_file, line)) {
        line_num++;
        if (line.empty())
            continue;
        
        // Split key from value
        pos = line.find(_delim);
        if (pos == std::string::npos)
            throw ConfigParseException("No '" + std::to_string(_delim) + "' sign in line "
                 + std::to_string(line_num));
        
        field = line.substr(0, pos);
        rtrim(field);
        if (result->hasParam(field))
            throw ConfigParseException("Redefinition of field \"" + field + "\" in line "
                + std::to_string(line_num));
        
        if (pos == line.length() - 1)
            value = "";
        else
            value = line.substr(pos + 1);
        ltrim(value);
        
        // Store key and value
        result->fieldMap[field] = value;
        result->keys.push_back(field);
    }

    return result;   
}

// Returns parameter value for given field name, or throws ConfigNoFieldException if not
// defined
//--------------------------------------------------------------------------------------------    
std::string Config::getParam(const std::string & _field) const
{
    auto iter = this->fieldMap.find(_field);
    if (iter == this->fieldMap.end())
        throw ConfigNoFieldException("No \"" + _field + "\" field in config");
    return (*iter).second;
}

// Returns true if a field of given name exists, false otherwise
//--------------------------------------------------------------------------------------------
bool Config::hasParam(const std::string & _field) const
{
    return this->fieldMap.find(_field) != this->fieldMap.end();
}

// Same as Config::getParam
//--------------------------------------------------------------------------------------------
std::string Config::operator[](const std::string & _field) const
{
    return this->getParam(_field);
}

// Returns a vector containing all keys in config
//--------------------------------------------------------------------------------------------
std::vector<std::string> Config::getKeys() const
{
    return this->keys;
}

