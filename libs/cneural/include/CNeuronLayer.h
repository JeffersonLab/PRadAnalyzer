#ifndef C_NEURON_LAYER_H
#define C_NEURON_LAYER_H

#include <vector>
#include "CNeuron.h"

class CNeuronLayer
{
public:
    CNeuronLayer(unsigned int con_size, unsigned int size);
    CNeuronLayer(CNeuronLayer &prev_layer, unsigned int size);

    void Connect(CNeuronLayer &prev_layer);
    void Update(const std::vector<double> &input);
    void Update();
    std::vector<CNeuron> &GetNeurons() {return neurons;}
    const std::vector<CNeuron> &GetNeurons() const {return neurons;}
    unsigned int GetInputSize() const {return input_size;}

private:
    unsigned int input_size;
    std::vector<CNeuron> neurons;
};

#endif
