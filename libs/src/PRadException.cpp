//============================================================================//
// A exception class for PRad Event Viewer                                    //
//                                                                            //
// Chao Peng                                                                  //
// 02/27/2016                                                                 //
//============================================================================//

#include "PRadException.h"

using namespace std;

PRadException::PRadException(const string &typ, const string &txt, const string &aux)
: title(typ), text(txt), auxText(aux)
{
}

PRadException::PRadException(PRadExceptionType typ, const string &txt, const string &aux)
: type(typ), text(txt), auxText(aux)
{
}

PRadException::PRadException(PRadExceptionType typ, const string &txt, const string &file, const string &func, int line)
: type(typ), text(txt)
{
    ostringstream oss;
    oss <<  "    evioException occured in file " << file << ", function " << func << ", line " << line;
    auxText=oss.str();
}


string PRadException::FailureDesc(void) const throw()
{
    ostringstream oss;
    oss << text << endl
        << auxText;
    return(oss.str());
}

string PRadException::FailureType(void) const throw()
{
    if(!title.empty())
        return title;

    string oss;
    switch(type)
    {
    case ET_CONNECT_ERROR:
        oss = "ET CONNECT ERROR";
        break;
    case ET_CONFIG_ERROR:
        oss = "ET CONFIG ERROR";
        break;
    case ET_STATION_CONFIG_ERROR:
        oss = "ET STATION CONFIG ERROR";
        break;
    case ET_STATION_CREATE_ERROR:
        oss = "ET STATION CREATE ERROR";
        break;
    case ET_STATION_ATTACH_ERROR:
        oss = "ET ATTACH ERROR";
        break;
    case ET_READ_ERROR:
        oss = "ET READ ERROR";
        break;
    case ET_PUT_ERROR:
        oss = "ET PUT ERROR";
        break;
    case HIGH_VOLTAGE_ERROR:
        oss = "HIGH VOLTAGE SYSTEM ERROR";
        break;
    default:
        oss = "UNKNOWN ERROR";
        break;
    }

    return(oss);
}

const char* PRadException::what() const throw()
{
    string failure = FailureType() + ": " + FailureDesc();
    return failure.c_str();
}
