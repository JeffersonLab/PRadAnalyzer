//============================================================================//
// A Neural Network Class                                                     //
//                                                                            //
// Chao Peng                                                                  //
// 01/21/2017                                                                 //
//============================================================================//

#include "CNeuralNetwork.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>



// constructor
CNeuralNetwork::CNeuralNetwork(double factor)
: learn_factor(factor), output_norm(1.0), output_shift(0.0)
{
    // place holder
}

// print out the structure
inline void __cnn_print_structure(const std::vector<CNeuronLayer> &layers)
{
    std::cout << "Neural network structure: " << std::endl;
    if(layers.empty()) {
        std::cout << "EMPTY!" << std::endl;
    } else {
        std::cout << std::setw(6) << " " << "Input dimension: "
                  << layers.front().GetInputSize() << std::endl
                  << std::setw(6) << " " << "Output dimension: "
                  << layers.back().GetNeurons().size() << std::endl
                  << std::setw(6) << " " << "Hidden layers: "
                  << layers.front().GetNeurons().size();
        for(unsigned int i = 1; i < layers.size() - 1; ++i)
        {
            std::cout << ", " << layers.at(i).GetNeurons().size();
        }
        std::cout << std::endl;
    }

}

// create a network
unsigned int CNeuralNetwork::CreateNet(unsigned int input_size,
                                       unsigned int output_size,
                                       const std::vector<unsigned int> &hidden_layers)
{
    // erase original layers
    layers.clear();

    // sanity check
    if(input_size == 0 || output_size == 0)
        return 0;

    // hidden layers
    for(auto &ilayer : hidden_layers)
    {
        if(ilayer <= 0)
            continue;
        layers.emplace_back(input_size, ilayer);
        input_size = ilayer;
    }

    // output layer
    layers.emplace_back(input_size, output_size);

    // build connections
    for(unsigned int i = 1; i < layers.size(); ++i)
    {
        layers.at(i).Connect(layers.at(i-1));
    }

    std::cout << "Create a new neural network" << std::endl;
    __cnn_print_structure(layers);

    return layers.size();
}

// helper functions for reading binary file
inline void __cnn_read_uint32(std::ifstream &ifs, uint32_t &word)
{
    if(!ifs.eof()) ifs.read((char*) &word, sizeof(uint32_t));
}

inline void __cnn_read_real64(std::ifstream &ifs, double &word)
{
    if(!ifs.eof()) ifs.read((char*) &word, sizeof(double));
}

// create a network from saved data file
unsigned int CNeuralNetwork::CreateNet(const char *path)
{
    std::ifstream inf(path, std::ios::in|std::ios::binary);
    if(!inf.is_open())
    {
        std::cerr << "Cannot open file " << path
                  << ", failed to create network from file."
                  << std::endl;
        return 0;
    }

    layers.clear();

    uint32_t uint_word;
    double real_word;

    __cnn_read_uint32(inf, uint_word);
    unsigned int input_size = uint_word;

    __cnn_read_uint32(inf, uint_word);
    unsigned int layers_size = uint_word;

    for(unsigned int i = 0; i < layers_size; ++i)
    {
        // create a layer
        __cnn_read_uint32(inf, uint_word);
        unsigned int neurons_size = uint_word;
        CNeuronLayer new_layer(input_size, neurons_size);
        input_size = neurons_size;

        // read weights from the layer
        for(unsigned int i = 0; i < neurons_size; ++i)
        {
            CNeuron &neuron = new_layer.GetNeurons().at(i);
            __cnn_read_uint32(inf, uint_word);
            unsigned int weights_size = uint_word;
            std::vector<double> weights(weights_size);
            for(unsigned int j = 0; j < weights_size; ++j)
            {
                __cnn_read_real64(inf, real_word);
                weights[j] = real_word;
            }
            neuron.SetWeights(weights);
        }

        layers.emplace_back(std::move(new_layer));
    }

    // read in normalization factor
    __cnn_read_real64(inf, output_norm);
    __cnn_read_real64(inf, output_shift);

    // build connections
    for(unsigned int i = 1; i < layers.size(); ++i)
    {
        layers.at(i).Connect(layers.at(i-1));
    }

    std::cout << "Create neural network from file "
              << "\"" << path << "\"\n"
              << "Output Normalization Factor: " << output_norm << "\n"
              << "Output Shift: " << output_shift
              << std::endl;

    __cnn_print_structure(layers);

    return layers.size();
}

// initialize all the weights, assign them with random numbers
void CNeuralNetwork::InitializeWeights()
{
    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::uniform_real_distribution<double> uni_dist(-1., 1.);

    for(auto &layer : layers)
    {
        for(auto &neuron : layer.GetNeurons())
        {
            std::vector<double> weights;
            weights.reserve(neuron.GetWeightSize());
            for(unsigned int i = 0; i < neuron.GetWeightSize(); ++i)
            {
                weights.push_back(uni_dist(rng));
            }
            neuron.SetWeights(weights);
        }
    }
}

// helper functions for writing binary file
inline void __cnn_write_uint32(std::ofstream &ofs, uint32_t word)
{
    ofs.write((char*) &word, sizeof(uint32_t));
}

inline void __cnn_write_real64(std::ofstream &ofs, double word)
{
    ofs.write((char*) &word, sizeof(double));
}

// save network to a data file
void CNeuralNetwork::SaveNet(const char *path)
const
{
    if(layers.empty())
    {
        std::cout << "This network has 0 layer, abort saving to " << path
                  << std::endl;
        return;
    }

    std::ofstream outf(path, std::ios::out | std::ios::binary);

    if(!outf.is_open())
    {
        std::cerr << "Cannot open file " << path
                  << ", save aborted!"
                  << std::endl;
        return;
    }

    // input size
    __cnn_write_uint32(outf, layers.front().GetInputSize());
    // number of total layers
    __cnn_write_uint32(outf, layers.size());
    for(auto &layer : layers)
    {
        // layer size
        __cnn_write_uint32(outf, layer.GetNeurons().size());
        // weights
        for(auto &neuron : layer.GetNeurons())
        {
            // number of weights
            auto weights = neuron.GetWeights();
            __cnn_write_uint32(outf, weights.size());
            for(auto &weight : weights)
            {
                __cnn_write_real64(outf, weight);
            }
        }
    }

    // output normalization factor
    __cnn_write_real64(outf, output_norm);
    __cnn_write_real64(outf, output_shift);

    std::cout << "Neural network data saved to "
              << "\"" << path << "\""
              << std::endl;
}

// give the network an input array and get the output array
void CNeuralNetwork::Update(const std::vector<double> &input)
{
    // first layer
    layers[0].Update(input);

    // the other layers will look for its previous layer's signals
    for(unsigned int i = 1; i < layers.size(); ++i)
    {
        layers[i].Update();
    }

    // save output
    auto &outn = layers.back().GetNeurons();
    output.clear();
    for(auto &neuron : outn)
        output.push_back((neuron.signal + output_shift)*output_norm);
}

// Training with erro back propagation
void CNeuralNetwork::BP_Train(const std::vector<double> &input,
                              const std::vector<double> &expect)
{
    // sanity check
    if(layers.empty())
        return;

    Update(input);
    std::vector<double> output = GetOutput();

    if(output.size() != expect.size())
    {
        std::cerr << "Unmatched dimension between expectation values ("
                  << expect.size() << ") and output values ("
                  << output.size() << "), abort back propagation training."
                  << std::endl;
        return;
    }

    // start from output layer
    auto &out_neurons = layers.back().GetNeurons();

    for(unsigned int i = 0; i < out_neurons.size(); ++i)
    {
        double dE = (output.at(i) - expect.at(i))/output_norm;
        // set error
        out_neurons.at(i).BP_Init(dE);
    }

    // initialize the response for other layers
    for(unsigned int i = 0; i < layers.size() - 1;  ++i)
    {
        auto &layer = layers.at(i);

        for(auto &neuron : layer.GetNeurons())
        {
            neuron.BP_Init(0.);
        }
    }

    // propagate responses backwardly
    for(auto it = layers.rbegin(); it != layers.rend(); ++it)
    {
        for(auto &neuron : it->GetNeurons())
        {
            neuron.BP_Propagate();
        }
    }

    // update weights for all the neurons
    // input layer is specially treated
    auto &in_neurons = layers.front().GetNeurons();
    for(auto &neuron : in_neurons)
    {
        neuron.BP_Learn(input, learn_factor);
    }

    for(unsigned int i = 1; i < layers.size(); ++i)
    {
        auto &layer = layers.at(i);

        for(auto &neuron : layer.GetNeurons())
        {
            neuron.BP_Learn(learn_factor);
        }
    }
}
