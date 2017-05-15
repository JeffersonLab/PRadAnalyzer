#ifndef C_NEURON_H
#define C_NEURON_H

#include <vector>

class CNeuron
{
public:
    friend class CNeuronLayer;
    friend class CNeuralNetwork;

    struct Connection
    {
        double weight;
        CNeuron *neuron;

        Connection()
        : weight(0.), neuron(nullptr)
        {}

        Connection(double w, CNeuron *n)
        : weight(w), neuron(n)
        {}
    };

public:
    CNeuron(unsigned int size);
    CNeuron(std::vector<CNeuron> &neurons);

    void Connect(unsigned int idx, CNeuron *n);
    void Connect(std::vector<CNeuron> &neurons);
    void SetWeights(const std::vector<double> &w);
    void Update(const std::vector<double> &input);
    void Update();
    std::vector<double> GetWeights() const;
    unsigned int GetWeightSize() const {return connections.size() + 1;}

    // Back Propagation Trainning functions
    void BP_Init(const double &dE);
    void BP_Propagate();
    void BP_Learn(const double &factor);
    void BP_Learn(const std::vector<double> &input, const double &factor);

private:
    double sigmoid(const double &a, const double &p) const;
    double dsigmoid(const double &a, const double &p) const;

private:
    std::vector<Connection>	connections;
    double bias;
    double signal;
    double sigma;
};

#endif
