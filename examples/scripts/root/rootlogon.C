// loading the library while using root
{
    // load library
    string prad_lib_path = getenv("PRAD_LIB");
    gSystem->Load((prad_lib_path + "/libprana.so").c_str());
    gSystem->Load((prad_lib_path + "/libcneural.so").c_str());
    // load include folder
    gSystem->AddIncludePath((prad_lib_path + "/include").c_str());
    // add include folder for the interpreter
    gInterpreter->AddIncludePath((prad_lib_path + "/include").c_str());
}

