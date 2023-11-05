#include <iostream>
#include <vector>
#include <utility>
using namespace std;
class Evento{
    private: 
    pair<double, double> pi_;
    pair<double,double> k_;
    public:
    Evento(){};
    Evento(double p_phi, double p_theta, double k_phi, double k_theta);
    void set_pi(double phi, double theta);
    void set_k(double phi, double theta);
    pair<double, double> get_pi();
    pair<double, double> get_k();
    ~Evento(){};

};