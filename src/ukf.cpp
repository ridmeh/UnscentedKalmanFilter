
#include "Eigen/Dense"
#include "measurement_package.h"
#include "tools.h"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;
#include "ukf.h"

using namespace std;
// Process noise standard deviation longitudinal acceleration in m/s^2
double UKF::std_a_ = 0.2;

// Process noise standard deviation yaw acceleration in rad/s^2
double UKF::std_yawdd_ = 0.2;

// Laser measurement noise standard deviation position1 in m
double UKF::std_laspx_ = 0.15;

// Laser measurement noise standard deviation position2 in m
double UKF::std_laspy_ = 0.15;

// Radar measurement noise standard deviation radius in m
double UKF::std_radr_ = 0.3;

// Radar measurement noise standard deviation angle in rad
double UKF::std_radphi_ = 0.03;

// Radar measurement noise standard deviation radius change in m/s
double UKF::std_radrd_ = 0.3;

// State dimension
int UKF::n_x_ =5;

// Augmented state dimension
int UKF::n_aug_=7;

// Sigma point spreading parameter
double UKF::lambda_ = 3 - n_aug_;

// Initialize weights.
VectorXd UKF::weights_ = VectorXd(2 * UKF::n_aug_ + 1);


/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;
	is_initialized_ = false;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);

	weights_.fill(0.5 / (n_aug_ + lambda_));
	weights_(0) = lambda_ / (lambda_ + n_aug_);

	R_radar = MatrixXd(3, 3);
	R_radar << std_radr_ * std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;

	R_lidar = MatrixXd(2, 2);
	R_lidar << std_laspx_ * std_laspx_, 0,
		0, std_laspy_*std_laspy_;
	// Process noise standard deviation longitudinal acceleration in m/s^2
	//UKF::std_a_ = 30;

	// Process noise standard deviation yaw acceleration in rad/s^2
	//std_yawdd_ = 30;

	/**
	 * DO NOT MODIFY measurement noise values below.
	 * These are provided by the sensor manufacturer.
	 */

	 // Laser measurement noise standard deviation position1 in m
	//std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	//std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	//std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	//std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	//std_radrd_ = 0.3;

	/**
	 * End DO NOT MODIFY section for measurement noise values
	 */

	 /**
	  * TODO: Complete the initialization. See ukf.h for other member properties.
	  * Hint: one or more values initialized above might be wildly off...
	  */
	  // set state dimension
	//n_x_ = 5;
	// set augmented dimension
	//n_aug_ = 7;

	// Process noise standard deviation longitudinal acceleration in m/s^2
	//std_a_ = 0.2;

	// Process noise standard deviation yaw acceleration in rad/s^2
	//std_yawdd_ = 0.2;

	// define spreading parameter
	//lambda_ = 3 - n_aug_;

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
	/**
	 * TODO: Complete this function! Make sure you switch between lidar and radar
	 * measurements.
	 */


	// Initialize

	if (!is_initialized_) {
		/**
		 * TODO: Initialize the state ekf_.x_ with the first measurement.
		 * TODO: Create the covariance matrix.
		 * You'll need to convert radar from polar to cartesian coordinates.
		 */

		 // first measurement
		x_ = VectorXd::Zero(n_x_);

		if (measurement_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
			// TODO: Convert radar from polar to cartesian coordinates 
			//         and initialize state.
			x_(0) = measurement_pack.raw_measurements_(0);
			x_(1) = measurement_pack.raw_measurements_(1);
			x_(2) = 10;
			x_(3) = 0;
			x_(4) = .5;

		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
			// TODO: Initialize state.
			float ro = measurement_pack.raw_measurements_(0);
			float theta = measurement_pack.raw_measurements_(1);
			float ro_dot = measurement_pack.raw_measurements_(2);

			// Refer : https://www.mathsisfun.com/polar-cartesian-coordinates.html

			// Conver r,0 to x,y co-ordindate
			x_(0) = ro * cos(theta);
			x_(1) = ro * sin(theta);
			
			//ekf_.x_(2) = ro_dot * cos(theta);
			//ekf_.x_(3) = ro_dot * sin(theta);
			x_(2) = 10;
			x_(3) = 0;
			x_(4) = .5;

		}
		previous_timestamp_ = measurement_pack.timestamp_;
		P_ << 0.2, 0, 0, 0, 0,
			0, 0.2, 0, 0, 0,
			0, 0, 1, 0, 0,
			0, 0, 0, .1, 0,
			0, 0, 0, 0, .1;
		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}

	// 
	float delta_t = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;

	Prediction(delta_t);

	Tools tools = Tools();
	MatrixXd Xsig_outASP;
	//MatrixXd Xsig_outSPP;
	// Generate sigma points
	Xsig_outASP = UKF::AugmentedSigmaPoints(x_, P_);
	// Predict sigma points
	Xsig_pred_ = UKF::SigmaPointPrediction(Xsig_outASP, delta_t);
	// Predict Mean and co-variance
	UKF::PredictMeanAndCovariance(Xsig_pred_, &x_, &P_);

	
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
	{
		UpdateRadar(measurement_pack);
		 
	}
	else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_)
	{
		UpdateLidar(measurement_pack);
	}
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
}

// Predict 
// 1. Generate Sigma points
// 2. Predict Sigma points
// 3. Predict Mean and Co-variance



void UKF::UpdateRadar(MeasurementPackage meas_package) {
	// set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 3; // or 3


	// set vector for weights
	//VectorXd weights = VectorXd(2 * n_aug_ + 1);
	//double weight_0 = lambda_ / (lambda_ + n_aug_);
	//double weight = 0.5 / (lambda_ + n_aug_);
	//weights(0) = weight_0;

	//for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
	//	weights(i) = weight;
	//}

	// create example matrix with predicted sigma points in state space
	MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
	Xsig_pred = Xsig_pred_;

	// create example vector for predicted state mean
	VectorXd x = VectorXd(n_x_);
	x = x_;


	// create example matrix for predicted state covariance
	//MatrixXd P = MatrixXd(n_x_, n_x_);
	MatrixXd P = P_;


	// create example matrix with sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	Zsig.fill(0);
	// create example vector for mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		//transform sigma points into measurement space
		VectorXd state_vec = Xsig_pred.col(i);
		double px = state_vec(0);
		double py = state_vec(1);
		double v = state_vec(2);
		double yaw = state_vec(3);
		//double yaw_d = state_vec(4);

		double rho = sqrt(px*px + py * py);
		double phi = atan2(py, px);
		double  rho_d = (px*cos(yaw)*v + py * sin(yaw)*v) / rho;

		Zsig.col(i) << rho,
			phi,
			rho_d;

		//calculate mean predicted measurement
		z_pred += weights_(i) * Zsig.col(i);
	}



	// create example matrix for predicted measurement covariance
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd z_diff = Zsig.col(i) - z_pred;
		if (z_diff(1) > M_PI) {
			z_diff(1) -= 2. * M_PI;
		}
		else if (z_diff(1) < -M_PI) {
			z_diff(1) += 2. * M_PI;
		}
		S += weights_(i) * z_diff * z_diff.transpose();
	}
	// Add R to S

	S += R_radar;
	// create example vector for incoming radar measurement
	VectorXd z = VectorXd(n_z);

	z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), meas_package.raw_measurements_(2);

	 // create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	/**
	 * ---
	 */

	 // calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
	  // residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		// angle normalization
		while (z_diff(1) > M_PI) {z_diff(1) -= 2.*M_PI;}
		while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred.col(i) - x;
		// angle normalization
		while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	// Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	// residual
	VectorXd z_diff = z - z_pred;

	// angle normalization
	while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

	// update state mean and covariance matrix
	x = x + K * z_diff;
	P = P - K * S*K.transpose();

	/**
	 * Student part end
	 */

	 // print result
	cout << "Updated state x: " << std::endl << x << std::endl;
	cout << "Updated state covariance P: " << std::endl << P << std::endl;

	// write result
	x_ = x;
	P_ = P;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
	// set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 2; // or 3


	// set vector for weights
	//VectorXd weights = VectorXd(2 * n_aug_ + 1);
	//double weight_0 = lambda_ / (lambda_ + n_aug_);
	//double weight = 0.5 / (lambda_ + n_aug_);
	//weights(0) = weight_0;

	//for (int i = 1; i < 2 * n_aug_ + 1; ++i) {
	//	weights(i) = weight;
	//}

	// create example matrix with predicted sigma points in state space
	MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
	Xsig_pred = Xsig_pred_;
	// create vector for predicted state mean
	VectorXd x = VectorXd(n_x_);
	x = x_;
	// create matrix for predicted state covariance
	MatrixXd P = MatrixXd(n_x_, n_x_);
	P = P_;

	// create matrix with sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	// create vector for mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		//transform sigma points into measurement space
		VectorXd state_vec = Xsig_pred.col(i);
		double px = state_vec(0);
		double py = state_vec(1);

		Zsig.col(i) << px,
			py;

		//calculate mean predicted measurement
		z_pred += weights_(i) * Zsig.col(i);
	}
	// create example matrix for predicted measurement covariance
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S += weights_(i) * z_diff * z_diff.transpose();
	}

	// Add R to S
	S += R_lidar;

	// create example vector for incoming radar measurement
	VectorXd z = VectorXd(n_z);

	z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);
	// create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	/**
	 * ---
	 */

	 // calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
	  // residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		// state difference
		VectorXd x_diff = Xsig_pred.col(i) - x;


		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	// Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	// residual
	VectorXd z_diff = z - z_pred;

	// update state mean and covariance matrix
	x = x + K * z_diff;
	P = P - K * S*K.transpose();

	/**
	 * ---
	 */

	 // print result
	cout << "Updated state x: " << std::endl << x << std::endl;
	cout << "Updated state covariance P: " << std::endl << P << std::endl;

	// write result
	x_= x;
	P_ = P;
}

/**
 * Programming assignment functions:
 */
MatrixXd UKF::AugmentedSigmaPoints( VectorXd& x,  MatrixXd& P) {

	// set example state
	// create augmented mean vector
	VectorXd x_aug = VectorXd(7);

	// create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);

	// create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	/**
	 *---
	 */

	 // create augmented mean state
	x_aug.head(5) = x;
	x_aug(5) = 0;
	x_aug(6) = 0;

	// create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	// create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	// create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i < n_aug_; ++i) {
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
	}

	/**
	 * ---
	 */

	 // print result
	cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

	// write result
	//*Xsig_out = Xsig_aug;
	return Xsig_aug;
}

MatrixXd UKF::SigmaPointPrediction( MatrixXd& Xsig_aug, double delta_t) {

	// create example sigma point matrix
	//MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

	// create matrix with predicted sigma points as columns
	MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

	//double delta_t = 0.1; // time diff in sec

	/**
	 * --
	 */

	 // predict sigma points
	for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
		// extract values for better readability
		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		// predicted state values
		double px_p, py_p;

		// avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
		}
		else {
			px_p = p_x + v * delta_t*cos(yaw);
			py_p = p_y + v * delta_t*sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd * delta_t;
		double yawd_p = yawd;

		// add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a * delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd * delta_t;

		// write predicted sigma point into right column
		Xsig_pred(0, i) = px_p;
		Xsig_pred(1, i) = py_p;
		Xsig_pred(2, i) = v_p;
		Xsig_pred(3, i) = yaw_p;
		Xsig_pred(4, i) = yawd_p;
	}

	/**
	 *--
	 */

	 // print result
	cout << "Xsig_pred = " << endl << Xsig_pred << endl;

	// write result
	return Xsig_pred;
}


void UKF::PredictMeanAndCovariance( MatrixXd& Xsig_pred, VectorXd* x_out, MatrixXd* P_out) {


	// create example matrix with predicted sigma points
	//MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

	// create vector for weights
	//VectorXd weights = VectorXd(2 * n_aug_ + 1);

	// create vector for predicted state
	VectorXd x = VectorXd(n_x_);

	// create covariance matrix for prediction
	MatrixXd P = MatrixXd(n_x_, n_x_);


	/**
	 * ----
	 */

	 // set weights
	//double weight_0 = lambda_ / (lambda_ + n_aug_);
	//weights(0) = weight_0;
	//for (int i = 1; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 weights
	//	double weight = 0.5 / (n_aug_ + lambda_);
	//	weights(i) = weight;
	//}

	// predicted state mean
	x.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
		x = x + weights_(i) * Xsig_pred.col(i);
	}

	// predicted state covariance matrix
	P.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
	  // state difference
		VectorXd x_diff = Xsig_pred.col(i) - x;
		// angle normalization
		while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

		P = P + weights_(i) * x_diff * x_diff.transpose();
	}

	/**
	 * ----
	 */

	 // print result

	cout << "Predicted state" << endl;
	cout << x << std::endl;
	cout << "Predicted covariance matrix" << endl;
	cout << P << endl;

	// write result
	*x_out = x;
	*P_out = P;
}
