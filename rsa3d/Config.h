//--------------------------------------------------------------------------------------------
// Simple configuration file parser
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CONFIG_H
    #define _CONFIG_H


#include <map>
#include <string>    
#include <fstream>
#include <stdexcept>
#include <memory>
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
    typedef std::pair<std::string, std::string> map_entry;
    typedef std::unique_ptr<Config> ptr;
    
    std::map<std::string, std::string>  fieldMap;
    std::vector<std::string>            keys;
    
public:
    static ptr parse(std::ifstream & _file, char _delim = '=');
    
    std::string getParam(const std::string & _field) const;
    bool        hasParam(const std::string & _field) const;
    std::string operator[](const std::string & _field) const;
    std::vector<std::string> getKeys() const;
};

    
#endif // _CONFIG_H
