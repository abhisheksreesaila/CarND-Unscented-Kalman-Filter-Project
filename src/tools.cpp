#include "tools.h"

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    VectorXd rmse(4);
    rmse <<  0,0,0,0;
    if(estimations.size() !=0)
        if(estimations.size() == ground_truth.size())
        {
            for (int i=0; i < estimations.size(); ++i) {
                VectorXd diff =  estimations[i] - ground_truth[i];
                VectorXd power =  diff.array() * diff.array();
                rmse = rmse + power;
            }
        }
    rmse = rmse/estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;


}