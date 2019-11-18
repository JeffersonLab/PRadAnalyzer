// loading the library while using root
{
    // load library
    string prad_root = getenv("PRAD_PATH");
    string prad_lib = prad_root + "/lib";
    string prad_inc = prad_root + "/include";

    if (prad_root.empty()) {
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

