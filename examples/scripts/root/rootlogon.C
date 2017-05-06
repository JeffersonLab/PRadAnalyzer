// loading the library while using root
{
    // load library
    const char* prad_env = getenv("PRAD_LIB");
    string prad_lib_path;
    if(prad_env == nullptr) {
        cout << "WARNING: Environment variable PRAD_LIB is not defined." << endl;
        prad_lib_path = ".";
    } else {
        prad_lib_path = prad_env;
    }

    gSystem->Load((prad_lib_path + "/libprana.so").c_str());
    gSystem->Load((prad_lib_path + "/libcneural.so").c_str());
    // load include folder
    gSystem->AddIncludePath((prad_lib_path + "/include").c_str());
    // add include folder for the interpreter
    gInterpreter->AddIncludePath((prad_lib_path + "/include").c_str());
}

