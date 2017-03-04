//============================================================================//
// A string wrapper class that convert string to other data types             //
//                                                                            //
// Chao Peng                                                                  //
// 06/07/2016                                                                 //
//============================================================================//

#include "ConfigValue.h"
#include "ConfigParser.h"
#include <climits>
#include <algorithm>

using namespace std;



//============================================================================//
// Constructors, Destructor and Assignment Operators                          //
//============================================================================//

ConfigValue::ConfigValue(const string &value)
: _value(value)
{}

ConfigValue::ConfigValue(string &&value)
: _value(move(value))
{}

ConfigValue::ConfigValue(const char *value)
: _value(value)
{}

ConfigValue::ConfigValue(const bool &value)
{
    if(value)
        _value = "1";
    else
        _value = "0";
}

ConfigValue::ConfigValue(const int &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const long &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const long long &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const unsigned &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const unsigned long &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const unsigned long long &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const float &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const double &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const long double &value)
: _value(to_string(value))
{}

ConfigValue &ConfigValue::operator =(const string &str)
{
    (*this)._value = str;
    return *this;
}

ConfigValue &ConfigValue::operator =(string &&str)
{
    (*this)._value = move(str);
    return *this;
}


//============================================================================//
// Public Member functions                                                    //
//============================================================================//

bool ConfigValue::Bool()
const
{
    if((_value == "1") ||
       (ConfigParser::strcmp_case_insensitive(_value, "T")) ||
       (ConfigParser::strcmp_case_insensitive(_value, "True")) ||
       (ConfigParser::strcmp_case_insensitive(_value, "Y")) ||
       (ConfigParser::strcmp_case_insensitive(_value, "Yes")))
        return true;

    if((_value == "0") ||
       (ConfigParser::strcmp_case_insensitive(_value, "F")) ||
       (ConfigParser::strcmp_case_insensitive(_value, "False")) ||
       (ConfigParser::strcmp_case_insensitive(_value, "N")) ||
       (ConfigParser::strcmp_case_insensitive(_value, "No")))
        return false;

    cout << "Config Value: Failed to convert "
         << _value << " to bool type. Return false."
         << endl;
    return false;
}

char ConfigValue::Char()
const
{
    try {
       int value = stoi(_value);
       if(value > CHAR_MAX)
           cout << "Config Value: Limit exceeded while converting "
                << _value << " to char." << endl;
       return (char) value;
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to char. 0 returned." << endl;
             return 0;
    }
}

unsigned char ConfigValue::UChar()
const
{
    try {
       unsigned long value = stoul(_value);
       if(value > UCHAR_MAX)
           cout << "Config Value: Limit exceeded while converting "
                << _value << " to unsigned char." << endl;
       return (unsigned char) value;
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned char. 0 returned." << endl;
             return 0;
    }
}

short ConfigValue::Short()
const
{
    try {
       int value = stoi(_value);
       if(value > SHRT_MAX)
           cout << "Config Value: Limit exceeded while converting "
                << _value << " to short." << endl;
       return (short) value;
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to short. 0 returned." << endl;
             return 0;
    }
}

unsigned short ConfigValue::UShort()
const
{
    try {
       unsigned long value = stoul(_value);
       if(value > USHRT_MAX)
           cout << "Config Value: Limit exceeded while converting "
                << _value << " to unsigned short." << endl;
       return (unsigned short) value;
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned short. 0 returned." << endl;
             return 0;
    }
}

int ConfigValue::Int()
const
{
    try {
       return stoi(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to int. 0 returned." << endl;
             return 0;
    }
}

unsigned int ConfigValue::UInt()
const
{
    try {
        return (unsigned int)stoul(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned int. 0 returned." << endl;
             return 0;
    }
}

long ConfigValue::Long()
const
{
    try {
        return stol(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to long. 0 returned." << endl;
             return 0;
    }
}

long long ConfigValue::LongLong()
const
{
    try {
        return stoll(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to long long. 0 returned." << endl;
             return 0;
    }
}

unsigned long ConfigValue::ULong()
const
{
    try {
        return stoul(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned long. 0 returned." << endl;
             return 0;
    }
}

unsigned long long ConfigValue::ULongLong()
const
{
    try {
        return stoull(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned long long. 0 returned." << endl;
             return 0;
    }
}

float ConfigValue::Float()
const
{
    try {
        return stof(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to float. 0 returned." << endl;
             return 0;
    }
}

double ConfigValue::Double()
const
{
    try {
        return stod(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to double. 0 returned." << endl;
             return 0;
    }
}

long double ConfigValue::LongDouble()
const
{
    try {
        return stold(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to long double. 0 returned." << endl;
             return 0;
    }
}

const char *ConfigValue::c_str()
const
{
    return _value.c_str();
}



//============================================================================//
// Other functions                                                            //
//============================================================================//

ostream &operator << (ostream &os, const ConfigValue &b)
{
    return  os << b.c_str();
};
