#ifndef PRAD_EXCEPTION_H
#define PRAD_EXCEPTION_H


#include <stdlib.h>
#include <string.h>
#include <exception>
#include <string>
#include <sstream>

class PRadException : public std::exception
{
public:
    enum PRadExceptionType
    {
        UNKNOWN_ERROR,
        ET_CONNECT_ERROR,
        ET_CONFIG_ERROR,
        ET_STATION_CONFIG_ERROR,
        ET_STATION_CREATE_ERROR,
        ET_STATION_ATTACH_ERROR,
        ET_READ_ERROR,
        ET_PUT_ERROR,
        HIGH_VOLTAGE_ERROR,
    };
    PRadException(const std::string &typ, const std::string &txt = "", const std::string &aux = "");
    PRadException(PRadExceptionType typ = UNKNOWN_ERROR, const std::string &txt = "", const std::string &aux = "");
    PRadException(PRadExceptionType typ, const std::string &txt, const std::string &file, const std::string &func, int line);
    virtual ~PRadException(void) throw() {};
    virtual std::string FailureDesc(void) const throw();
    virtual std::string FailureType(void) const throw();
    const char *what() const throw();

public:
    PRadExceptionType type;             // exception type
    std::string title;
    std::string text;     // primary text
    std::string auxText;  // auxiliary text
};

#endif


