#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/* Constant representing the zero threshold of the algorithm. */
#define ZERO_THRESH 0.0001

/* Converts time from milliseconds to seconds. */
#define CONVERT_USEC_TIME_TO_SECONDS(usec_time) usec_time *= 1e-6;

/* Normalize an angle to be in range from -PI to PI. */
#define NORMALIZE_ANGLE(angle) 		             \
		 {                                         \
			 while (angle > M_PI) angle -= 2.*M_PI;  \
			 while (angle < -M_PI) angle += 2.*M_PI; \
		 }

/* set x to the chi-squere 0.95 value for r levels of freedom. */
#define CHI_SQUARE_95(r, x)    \
		{					   \
			if (r == 2)        \
				x = 5.991;     \
			else if (r == 3)   \
				x = 7.815;     \
			else               \
				x = 0;         \
		}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// if this is true, laser NIS measurements will be calculated and printed.
	calc_laser_nis_ = false;

	// if this is true, radar NIS measurements will be calculated and printed.
	calc_radar_nis_ = false;

	// Number of elements in x vector
	n_x_ = 5;

	// Number of noise elements
	int n_noise = 2;

	// Number elements in an augmented vector.
	n_aug_ = n_x_ + n_noise;

	// initial state vector
	x_ = VectorXd(n_x_);

	// initial covariance matrix
	P_ = MatrixXd::Identity(n_x_, n_x_);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 0.1;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.2;

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

	// Makes sure that the first call to `ProcessMeasurement()` will initialize the filter.
	is_initialized_ = false;

	// The time in milliseconds of the latest measurement.
	time_us_ = 0;

	// The lambda value for the sigma points calculation.
	lambda_ = 3 - n_aug_;

	// Initialize weights
	weights_ = VectorXd(2 * n_aug_ + 1);

	// Number of laser samples that were processed by the filter.
	laser_n_samples_ = 0;

	// Number of radar samples that were processed by the filter.
	radar_n_samples_ = 0;

	// Number of laser samples that has scored below the threshold in NIS test.
	laser_samples_below_thold_ = 0;

	// Number of radar samples that has scored below the threshold in NIS test.
	radar_samples_below_thold_ = 0;

}

UKF::~UKF() {
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

	if (!is_initialized_) {
		is_initialized_ = true;

		// Set the first measurement.
		x_.setZero();
		if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) {
			x_(0) = meas_package.raw_measurements_(0);
			x_(1) = meas_package.raw_measurements_(1);

			P_(0, 0) = std_laspx_ * std_laspx_;
			P_(1, 1) = std_laspy_ * std_laspy_;
		} else if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
			double rho = meas_package.raw_measurements_(0);
			double phi = meas_package.raw_measurements_(1);

			x_(0) = rho * cos(phi);
			x_(1) = rho * sin(phi);


			// Set the px, py std of the P matrix to std_radr_^2
			// because var[r*cos(phi)] <= var[r] and var[r*sin(phi)] <= var[r]
			P_(0, 0) = std_radr_ * std_radr_;
			P_(1, 1) = std_radr_ * std_radr_;
		}

		time_us_ = meas_package.timestamp_;

		// calculate weights
		double weight_denominator = lambda_ + n_aug_;
		weights_(0) = lambda_ / weight_denominator;
		for (int weightIdx = 1; weightIdx < weights_.size(); weightIdx++) {
			weights_(weightIdx) = 1 / (2 * weight_denominator);
		}

		return;
	}

	double delta_t = meas_package.timestamp_ - time_us_;
	// convert deta_t into seconds.
	CONVERT_USEC_TIME_TO_SECONDS(delta_t);

	time_us_ = meas_package.timestamp_;

	if (delta_t > ZERO_THRESH) {
		Prediction(delta_t);
	} else {
		cout
				<< "Delta t is too small for anther prediction. Skipping prediction step."
				<< endl;
	}

	if (use_radar_
			&& meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
		UpdateRadar(meas_package);
	} else if (use_laser_
			&& meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) {
		UpdateLidar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

	/* ******************************************
	 * Create augmented covariance sigma points
	 * *****************************************
	 */

	int n_noise = n_aug_ - n_x_;
	double delta_t2 = delta_t * delta_t;
	double sigCoeff = sqrt(lambda_ + n_aug_);

	// Create the augmented x vector
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.setZero();
	x_aug.head(n_x_) = x_;
	/* no need to set mean values for the noises because it
	 * assumes that their mean is zero.
	 */

	// Create the augmented P matrix
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	P_aug.setZero();
	P_aug.topLeftCorner(n_x_, n_x_) << P_;
	P_aug.bottomRightCorner(n_noise, n_noise) << std_a_, 0, 0, std_yawdd_;
	//create square root matrix.
	MatrixXd A = P_aug.llt().matrixL();

	// Create sigma points matrix for x_
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	// Generate sigma points
	Xsig_aug.col(0) << x_aug;
	for (int colIdx = 1; colIdx <= n_aug_; colIdx++) {
		Xsig_aug.col(colIdx) << x_aug + (sigCoeff * A.col(colIdx - 1));
		Xsig_aug.col(n_aug_ + colIdx) << x_aug - (sigCoeff * A.col(colIdx - 1));
	}

	/* *************************************************************************
	 * Pass sigma points through the model equations and calculate the predicted
	 * covariance sigma points.
	 * *************************************************************************
	 */

	Xsig = MatrixXd(n_x_, Xsig_aug.cols());
	// Apply the non-linear model on the sigma points.
	for (int sigIdx = 0; sigIdx < Xsig.cols(); sigIdx++) {
		VectorXd x = Xsig_aug.col(sigIdx).head(n_x_);
		VectorXd nu = Xsig_aug.col(sigIdx).tail(n_noise);

		double v = x(2);
		double psi = x(3);
		double psi_d = x(4);

		VectorXd xpred_d = VectorXd(n_x_);
		if (psi_d > ZERO_THRESH) {
			double factor = v / psi_d;
			xpred_d << factor * (sin(psi + psi_d * delta_t) - sin(psi)), factor
					* (-cos(psi + psi_d * delta_t) + cos(psi)), 0, psi_d * delta_t, 0;
		} else {
			xpred_d << v * cos(psi) * delta_t, v * sin(psi) * delta_t, 0, psi_d
					* delta_t, 0;
		}

		double nu_a = nu(0);
		double nu_psi_dd = nu(1);

		VectorXd nu_vec = VectorXd(n_x_);
		nu_vec << 0.5 * delta_t2 * cos(psi) * nu_a, 0.5 * delta_t2 * sin(psi)
				* nu_a, delta_t * nu_a, 0.5 * delta_t2 * nu_psi_dd, delta_t * nu_psi_dd;

		Xsig.col(sigIdx) = x + xpred_d + nu_vec;
	}

	/* *************************************************************************
	 * Use the predicated sigma points to evaluate the predicated state
	 * mean and covariance.
	 * *************************************************************************
	 */

	//predict state mean
	x_.setZero();
	for (int sigIdx = 0; sigIdx < weights_.size(); sigIdx++) {
		x_ += weights_(sigIdx) * Xsig.col(sigIdx);
	}

	//predict state covariance matrix
	P_.setZero();
	for (int sigIdx = 0; sigIdx < weights_.size(); sigIdx++) {
		VectorXd diff = Xsig.col(sigIdx) - x_;
		NORMALIZE_ANGLE(diff(3));

		P_ += weights_(sigIdx) * diff * diff.transpose();
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

	VectorXd z = meas_package.raw_measurements_;
	int n_z = z.size();

	/* *************************************************************************
	 * Transform the predicted sigma points into the radar measurement space.
	 * *************************************************************************
	 */

	//transform sigma points into measurement space
	MatrixXd Zsig = MatrixXd(n_z, Xsig.cols());

	for (int sigIdx = 0; sigIdx < Xsig.cols(); sigIdx++) {
		double px = Xsig(0, sigIdx);
		double py = Xsig(1, sigIdx);

		Zsig.col(sigIdx) << px, py;
	}

	//calculate mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.setZero();
	for (int sigIdx = 0; sigIdx < Xsig.cols(); sigIdx++) {
		z_pred += weights_(sigIdx) * Zsig.col(sigIdx);
	}

	//calculate measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.setZero();
	for (int sigIdx = 0; sigIdx < Xsig.cols(); sigIdx++) {
		VectorXd diff = Zsig.col(sigIdx) - z_pred;

		S += weights_(sigIdx) * diff * diff.transpose();
	}

	// Add measurement noise R
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_;

	S += R;

	/* *****************************
	 * Update mean and covariance.
	 * *****************************
	 */

	//calculate cross correlation matrix
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.setZero();
	for (int sigIdx = 0; sigIdx < Xsig.cols(); sigIdx++) {
		VectorXd x_diff = Xsig.col(sigIdx) - x_;
		NORMALIZE_ANGLE(x_diff(3));

		VectorXd z_diff = Zsig.col(sigIdx) - z_pred;

		Tc += weights_(sigIdx) * x_diff * z_diff.transpose();
	}

	//calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//update state mean and covariance matrix
	VectorXd z_diff = z - z_pred;

	x_ += K * z_diff;
	P_ -= K * S * K.transpose();

	/* ********************
	 * Calc and print NIS.
	 * ********************
	 */

	if (calc_laser_nis_) {
		double nis = z_diff.transpose() * S.inverse() * z_diff;
		double nis_thold;

		CHI_SQUARE_95(n_z, nis_thold);

		if (nis_thold) {
			laser_n_samples_++;
			if (nis < nis_thold) {
				laser_samples_below_thold_++;
			}

			cout << "laser_consistency_precent = "
					<< laser_samples_below_thold_ * 100 / laser_n_samples_ << "%"
					<< endl;
		} else {
			cout << "Chi-square value for r=" << n_z
					<< " is unknown. Can't calculate NIS." << endl;
			calc_laser_nis_ = false;
		}
	}
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	VectorXd z = meas_package.raw_measurements_;
	int n_z = z.size();

	/* *************************************************************************
	 * Transform the predicted sigma points into the radar measurement space.
	 * *************************************************************************
	 */

	//transform sigma points into measurement space
	MatrixXd Zsig = MatrixXd(n_z, Xsig.cols());

	for (int sigIdx = 0; sigIdx < Xsig.cols(); sigIdx++) {
		double px = Xsig(0, sigIdx);
		double py = Xsig(1, sigIdx);
		double v = Xsig(2, sigIdx);
		double psi = Xsig(3, sigIdx);

		double rho = sqrt(px * px + py * py);
		double phi = atan2(py, px);
		double rho_dot = v * (px * cos(psi) + py * sin(psi)) / rho;

		Zsig.col(sigIdx) << rho, phi, rho_dot;
	}

	//calculate mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.setZero();
	for (int sigIdx = 0; sigIdx < Xsig.cols(); sigIdx++) {
		z_pred += weights_(sigIdx) * Zsig.col(sigIdx);
	}

	//calculate measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.setZero();
	for (int sigIdx = 0; sigIdx < Xsig.cols(); sigIdx++) {
		VectorXd diff = Zsig.col(sigIdx) - z_pred;
		NORMALIZE_ANGLE(diff(1));

		S += weights_(sigIdx) * diff * diff.transpose();
	}

	// Add measurement noise R
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_radr_ * std_radr_, 0, 0, 0, std_radphi_ * std_radphi_, 0, 0, 0, std_radrd_
			* std_radrd_;

	S += R;

	/* *****************************
	 * Update mean and covariance.
	 * *****************************
	 */

	//calculate cross correlation matrix
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.setZero();
	for (int sigIdx = 0; sigIdx < Xsig.cols(); sigIdx++) {
		VectorXd x_diff = Xsig.col(sigIdx) - x_;
		NORMALIZE_ANGLE(x_diff(3));

		VectorXd z_diff = Zsig.col(sigIdx) - z_pred;
		NORMALIZE_ANGLE(z_diff(1));

		Tc += weights_(sigIdx) * x_diff * z_diff.transpose();
	}

	//calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//update state mean and covariance matrix
	VectorXd z_diff = z - z_pred;
	NORMALIZE_ANGLE(z(1));

	x_ += K * z_diff;
	P_ -= K * S * K.transpose();

	/* ********************
	 * Calc and print NIS.
	 * ********************
	 */

	if (calc_radar_nis_) {
		double nis = z_diff.transpose() * S.inverse() * z_diff;
		double nis_thold;

		CHI_SQUARE_95(n_z, nis_thold);

		if (nis_thold) {
			radar_n_samples_++;
			if (nis < nis_thold) {
				radar_samples_below_thold_++;
			}

			cout << "radar_consistency_precent = "
					<< radar_samples_below_thold_ * 100 / radar_n_samples_ << "%"
					<< endl;
		} else {
			cout << "Chi-square value for r=" << n_z
					<< " is unknown. Can't calculate NIS." << endl;
			calc_radar_nis_ = false;
		}
	}
}
