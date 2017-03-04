g++ -std=c++11 -I../../include -c -o fileOrganizer.o fileOrganizer.cpp
g++ -Wl,-O1 -Wl,-z,relro -o fileOrganizer fileOrganizer.o ../../obj/ConfigParser.o ../../obj/ConfigValue.o

