//--------------------------------------------------------------------------------------------
// Simple configuration file parser
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CONFIG_H
    #define _CONFIG_H


#include <map>
#include <string>    
#include <istream>
#include <stdexcept>
#include <vector>


// Parse exception class
class ConfigParseException : public std::runtime_error
{
public:
    ConfigParseException(const std::string & _what) : std::runtime_error(_what)
    {}
};


// No field exception class
class ConfigNoFieldException : public std::runtime_error
{
public:
    ConfigNoFieldException(const std::string & _what) : std::runtime_error(_what)
    {}
};


class Config
{
private:
    std::map<std::string, std::string>  fieldMap;
    std::vector<std::string>            keys;
    
public:
    static Config * parse(std::istream & _file, char _delim = '=');
    
    std::string getString(const std::string & _field) const;
    int getInt(const std::string & _field) const;
    unsigned int getUnsignedInt(const std::string & _field) const;
    double getDouble(const std::string & _field) const;
    float getFloat(const std::string & _field) const;
    bool hasParam(const std::string & _field) const;
    std::vector<std::string> getKeys() const;
};

    
#endif // _CONFIG_H
