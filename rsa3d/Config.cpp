//--------------------------------------------------------------------------------------------
// Simple configuration file parser
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "Config.h"
#include "Utils.h"


// Parses configuration parameters from given istream
//--------------------------------------------------------------------------------------------
Config * Config::parse(std::istream & _file, char _delim)
{
    std::string line;
    std::string field, value;
    std::size_t line_num = 0;
    std::size_t pos;
    
    Config * result = new Config();
    
    // Parse line by line
    while (std::getline(_file, line)) {
        line_num++;
        if (line.empty() || line[0] == '#')
            continue;
        
        // Split key from value
        pos = line.find(_delim);
        if (pos == std::string::npos)
            throw ConfigParseException("No '" + std::to_string(_delim) + "' sign in line " + std::to_string(line_num));
        
        field = line.substr(0, pos);
        trim(field);
        if (result->hasParam(field))
            throw ConfigParseException("Redefinition of field \"" + field + "\" in line " + std::to_string(line_num));
        
        value = (pos == line.length() - 1) ? "" : line.substr(pos + 1);
        trim(value);
        
        // Store key and value
        result->fieldMap[field] = value;
        result->keys.push_back(field);
    }

    return result;   
}

// Returns parameter value for given field name as string, or throws ConfigNoFieldException
// if not defined
//--------------------------------------------------------------------------------------------    
std::string Config::getString(const std::string & _field) const
{
    auto iter = this->fieldMap.find(_field);
    if (iter == this->fieldMap.end())
        throw ConfigNoFieldException("No \"" + _field + "\" field in config");
    return (*iter).second;
}

// Other param getters parsing string to a specific type. They may throw some std:: exceptions
// while parsing
//-------------------------------------------------------------------------------------------- 
int Config::getInt(const std::string & _field) const
{
    return std::stoi(this->getString(_field));
}

unsigned int Config::getUnsignedInt(const std::string & _field) const
{
    return std::stoul(this->getString(_field));
}

double Config::getDouble(const std::string & _field) const
{
    return std::stod(this->getString(_field));
}

float Config::getFloat(const std::string & _field) const
{
    return std::stof(this->getString(_field));
}

// Returns true if a field of given name exists, false otherwise
//--------------------------------------------------------------------------------------------
bool Config::hasParam(const std::string & _field) const
{
    return this->fieldMap.find(_field) != this->fieldMap.end();
}

// Returns a vector containing all keys in config
//--------------------------------------------------------------------------------------------
std::vector<std::string> Config::getKeys() const
{
    return this->keys;
}

