#include <iostream>
#include <vector>
#include <utility>
#include "particle.h"

Evento::Evento(double p_phi, double p_theta, double k_phi, double k_theta){
    pi_.first=p_phi;
    pi_.second=p_theta;
    k_.first=k_phi;
    k_.second=k_theta;

}
void Evento::set_pi(double phi, double theta){
    pi_.first=phi;
    pi_.second=theta;
}
void Evento::set_k(double phi, double theta){
    k_.first=phi;
    k_.second=theta;

}
pair<double, double> Evento::get_pi(){
    return pi_;
}
pair<double, double> Evento::get_k(){
    return k_;
}