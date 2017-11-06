#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {
}

Tools::~Tools() {
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		const vector<VectorXd> &ground_truth) {

	VectorXd rmse(4);
	rmse.setZero();

	if (!estimations.size()) {
		cout << "CalculateRMSE Error - estimation size is zero." << endl;
		return rmse;
	}
	if (estimations.size() != ground_truth.size()) {
		cout
				<< "CalculateRMSE Error - estimation size does not match ground_truth size."
				<< endl;
		return rmse;
	}

	//accumulate squared residuals
	for (int i = 0; i < estimations.size(); i++) {
		// Calculate residuals
		VectorXd residual = estimations[i] - ground_truth[i];
		// Square residuals coeffs
		VectorXd squared = residual.array().square();
		rmse += squared;
	}

	//calculate the mean
	rmse /= estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}
