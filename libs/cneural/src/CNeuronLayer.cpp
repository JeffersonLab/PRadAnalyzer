//============================================================================//
// A layer of neurons                                                         //
//                                                                            //
// Chao Peng                                                                  //
// 01/21/2017                                                                 //
//============================================================================//

#include "CNeuronLayer.h"
#include <iostream>



// constructor
CNeuronLayer::CNeuronLayer(unsigned int input, unsigned int size)
: input_size(input)
{
	for(unsigned int i = 0; i < size; ++i)
    {
		neurons.emplace_back(input);
    }
}

CNeuronLayer::CNeuronLayer(CNeuronLayer &prev_layer, unsigned int size)
: input_size(prev_layer.GetNeurons().size())
{
    for(unsigned int i = 0; i < size; ++i)
    {
        neurons.emplace_back(prev_layer.GetNeurons());
    }
}

void CNeuronLayer::Connect(CNeuronLayer &prev_layer)
{
    for(auto &neuron : neurons)
    {
        neuron.Connect(prev_layer.GetNeurons());
    }
}

// give the layer an input array and get its output array
void CNeuronLayer::Update(const std::vector<double> &input)
{
    if(input.size() != input_size)
    {
        std::cerr << "Input size " << input.size() << " unmatches expected size "
                  << input_size << ", abort layer signals update."
                  << std::endl;
        return;
    }

    for(auto &neuron : neurons)
    {
        neuron.Update(input);
    }
}

void CNeuronLayer::Update()
{
    for(auto &neuron : neurons)
    {
        neuron.Update();
    }
}

