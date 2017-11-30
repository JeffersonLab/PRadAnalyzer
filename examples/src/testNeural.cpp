//============================================================================//
// An example showing how to use the neural network                           //
//                                                                            //
// Chao Peng                                                                  //
// 01/24/2017                                                                 //
//============================================================================//

#include "CNeuralNetwork.h"
#include <iostream>

#define NEW_NET

using namespace std;

int main(int /*argc*/, char * /*argv*/ [])
{
    CNeuralNetwork my_net;

#ifdef NEW_NET
    // create net with dimensions
    // 2 hidden layers have 5 and 3 neurons respectively
    std::vector<unsigned int> hidden = {15, 5, 3};
    // 5 inputs and 3 outputs with hidden layers
    my_net.CreateNet(5, 3, hidden);
    // initialize the weights with random values
    my_net.InitializeWeights();
#else
    // or create net from saved network data
    my_net.CreateNet("saved.net");
#endif

    // set input and expected output
    std::vector<double> input = {1, 2, 3, 4, 5};
    std::vector<double> expect = {0.3, 0.6, 0.9};
    std::vector<double> input2 = {5, 4, 3, 2, 1};
    std::vector<double> expect2 = {0.9, 0.6, 0.3};

    // output with random weights
    my_net.Update(input);
    cout << my_net.GetOutput().at(0) << ", "
         << my_net.GetOutput().at(1) << ", "
         << my_net.GetOutput().at(2)
         << endl;

    my_net.Update(input2);
    cout << my_net.GetOutput().at(0) << ", "
         << my_net.GetOutput().at(1) << ", "
         << my_net.GetOutput().at(2)
         << endl;

    // train it 20000 times
    int count = 20000;
    while(count-- > 0)
    {
        my_net.BP_Train(input, expect);
        my_net.BP_Train(input2, expect2);
    }

    // output after training
    my_net.Update(input);
    cout << my_net.GetOutput().at(0) << ", "
         << my_net.GetOutput().at(1) << ", "
         << my_net.GetOutput().at(2)
         << endl;

    my_net.Update(input2);
    cout << my_net.GetOutput().at(0) << ", "
         << my_net.GetOutput().at(1) << ", "
         << my_net.GetOutput().at(2)
         << endl;
    // save the result
    my_net.SaveNet("saved.net");
    return 0;
}
