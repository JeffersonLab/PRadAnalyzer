g++ -std=c++11 -I../../include -c -o mapCombiner.o mapCombiner.cpp
g++ -Wl,-O1 -Wl,-z,relro -o mapCombiner mapCombiner.o ../../obj/ConfigParser.o ../../obj/ConfigValue.o

