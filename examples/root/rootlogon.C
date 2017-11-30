// loading the library while using root
{
    // load library
    string prad_lib = getenv("PRAD_LIB");
    string prad_inc = getenv("PRAD_INC");

    if(prad_lib.empty() || prad_inc.empty()) {
        cout << "WARNING: Environment variable PRAD_LIB and PRAD_INC is not defined." << endl;
        prad_lib = ".";
        prad_inc = ".";
    }

    gSystem->Load((prad_lib + "/libprana.so").c_str());
    gSystem->Load((prad_lib + "/libprconf.so").c_str());
    gSystem->Load((prad_lib + "/libcana.so").c_str());
    gSystem->Load((prad_lib + "/libcneural.so").c_str());
    // load include folder
    gSystem->AddIncludePath(prad_inc.c_str());
    // add include folder for the interpreter
    gInterpreter->AddIncludePath(prad_inc.c_str());
}

