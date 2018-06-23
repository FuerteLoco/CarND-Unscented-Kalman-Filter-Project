#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // time when the state is true, in us
  //long long time_us_;

  // Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/

  if (!is_initialized_)
  {
    // first measurement
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      // radar measurement
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      x_ << rho * cos(phi), rho * sin(phi), 0, 0, 0;
    }
    else
    {
      // lidar measurement
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }
    
    P_ = MatrixXd::Identity(5, 5);
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }

  if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && (use_laser_ == false)) return;
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && (use_radar_ == false)) return;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // compute the time elapsed between the current and previous measurements
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;  //delta_t - expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else
  {
    UpdateLidar(meas_package);
  }
  
  return;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  // -------------------------------------------------------------------
  // Augmentation Assignment (17)
  // -------------------------------------------------------------------
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  //create augmented mean state
  x_aug.setZero();
  x_aug.head(n_x_) = x_;

  //create augmented covariance matrix
  P_aug.setZero();
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_,   n_x_  ) = std_a_     * std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;
  
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  
  //create augmented sigma points
  int col = 0;
  Xsig_aug.col(col++) = x_aug;
  A = sqrt(lambda_ + n_aug_) * A;
  for (int i=0; i<n_aug_; i++) Xsig_aug.col(col++) = x_aug + A.col(i);
  for (int i=0; i<n_aug_; i++) Xsig_aug.col(col++) = x_aug - A.col(i);
  
  // -------------------------------------------------------------------
  // Sigma Point Prediction (20)
  // -------------------------------------------------------------------
  
  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  double delta_t2 = delta_t * delta_t / 2.0;
  for (int col = 0; col < 2 * n_aug_ + 1; col++)
  {
    double v       = Xsig_aug(2, col);
    double psi     = Xsig_aug(3, col);
    double psi_dot = Xsig_aug(4, col);
    double nu_a    = Xsig_aug(5, col);
    double nu_psi  = Xsig_aug(6, col);

    VectorXd x_k = Xsig_aug.col(col).head(n_x_);
    VectorXd motion = VectorXd(n_x_);
    VectorXd noise = VectorXd(n_x_);

    noise << delta_t2 * cos(psi) * nu_a,
             delta_t2 * sin(psi) * nu_a,
             delta_t  * nu_a,
             delta_t2 * nu_psi,
             delta_t  * nu_psi;

    if (fabs(psi_dot) < 0.00001)
    {
      motion << v * cos(psi) * delta_t,
                v * sin(psi) * delta_t,
                0,
                psi_dot * delta_t,
                0;
    }
    else
    {
      motion << v / psi_dot * (sin(psi + psi_dot * delta_t) - sin(psi)),
                v / psi_dot * (-cos(psi + psi_dot * delta_t) + cos(psi)),
                0,
                psi_dot * delta_t,
                0;
    }
    Xsig_pred_.col(col) = x_k + motion + noise;
  }
  
  // -------------------------------------------------------------------
  // Predicted Mean and Covariance Assignment (23)
  // -------------------------------------------------------------------
  
  //set weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++) weights_(i) = 1 / (2 * (lambda_ + n_aug_));
  
  //predict state mean
  x_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; i++) x_ += weights_(i) * Xsig_pred_.col(i);
  
  //predict state covariance matrix
  P_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd diff = Xsig_pred_.col(i) - x_;
    while (diff(3) > M_PI) diff(3) -= 2.0 * M_PI;
    while (diff(3) < -M_PI) diff(3) += 2.0 * M_PI;
    P_ += weights_(i) * diff * diff.transpose();
  }
  
  return;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //set measurement dimension, lidar can measure x and y
  int n_z = 2;
    
  // -------------------------------------------------------------------
  // Predict Lidar Measurements Assignment (26)
  // -------------------------------------------------------------------

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  
  //transform sigma points into measurement space
  for (int col = 0; col < 2 * n_aug_ + 1; col++)
  {
      double p_x     = Xsig_pred_(0, col);
      double p_y     = Xsig_pred_(1, col);
      
      Zsig.col(col) << p_x, p_y;
  }

  //calculate mean predicted measurement
  z_pred.setZero();
  for (int col = 0; col < 2 * n_aug_ + 1; col++) z_pred += weights_(col) * Zsig.col(col);
  
  //calculate innovation covariance matrix S
  MatrixXd L = MatrixXd(n_z, n_z);
  L.setZero();
  L(0, 0) = std_laspx_ * std_laspx_;
  L(1, 1) = std_laspy_ * std_laspy_;

  S.setZero();
  for (int col = 0; col < 2 * n_aug_ + 1; col++)
  {
    VectorXd diff = Zsig.col(col) - z_pred;
    while (diff(1) > M_PI) diff(1) -= 2.0 * M_PI;
    while (diff(1) <- M_PI) diff(1) += 2.0 * M_PI;
    S += weights_(col) * diff * diff.transpose();
  }
  S += L;

  // -------------------------------------------------------------------
  // UKF Update Assignment (29)
  // -------------------------------------------------------------------

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  //calculate cross correlation matrix
  Tc.setZero();
  for (int col = 0; col < 2 * n_aug_ + 1; col++)
  {
    VectorXd z_diff = Zsig.col(col) - z_pred;
    while (z_diff(1) > M_PI) z_diff(1) -= 2.0 * M_PI;
    while (z_diff(1) <- M_PI) z_diff(1) += 2.0 * M_PI;

    VectorXd x_diff = Xsig_pred_.col(col) - x_;
    while (x_diff(3) > M_PI) x_diff(3) -= 2.0 * M_PI;
    while (x_diff(3) <- M_PI) x_diff(3) += 2.0 * M_PI;
    
    Tc += weights_(col) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  VectorXd z_diff = z - z_pred;
  while (z_diff(1) > M_PI) z_diff(1) -= 2.0 * M_PI;
  while (z_diff(1) <- M_PI) z_diff(1) += 2.0 * M_PI;
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  double nis_value = z_diff.transpose() * S.inverse() * z_diff;
  ofstream nis_file;
  nis_file.open("nis-lidar.txt", ios::app);
  nis_file << nis_value;
  nis_file << "\n";
  nis_file.close();
  
  return;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
    
  // -------------------------------------------------------------------
  // Predict Radar Measurements Assignment (26)
  // -------------------------------------------------------------------

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  
  //transform sigma points into measurement space
  for (int col = 0; col < 2 * n_aug_ + 1; col++)
  {
      double p_x     = Xsig_pred_(0, col);
      double p_y     = Xsig_pred_(1, col);
      double v       = Xsig_pred_(2, col);
      double psi     = Xsig_pred_(3, col);
      // double psi_dot = Xsig_pred_(4, col); // not needed
      
      double rho = sqrt(p_x * p_x + p_y * p_y);
      double phi = atan2(p_y, p_x);
      double rho_dot = (p_x * cos(psi) * v + p_y * sin(psi) * v) / rho;
      
      Zsig.col(col) << rho, phi, rho_dot;
  }

  //calculate mean predicted measurement
  z_pred.setZero();
  for (int col = 0; col < 2 * n_aug_ + 1; col++) z_pred += weights_(col) * Zsig.col(col);
  
  //calculate innovation covariance matrix S
  MatrixXd R = MatrixXd(n_z, n_z);
  R.setZero();
  R(0, 0) = std_radr_   * std_radr_;
  R(1, 1) = std_radphi_ * std_radphi_;
  R(2, 2) = std_radrd_  * std_radrd_;

  S.setZero();
  for (int col = 0; col < 2 * n_aug_ + 1; col++)
  {
    VectorXd diff = Zsig.col(col) - z_pred;
    while (diff(1) > M_PI) diff(1) -= 2.0 * M_PI;
    while (diff(1) <- M_PI) diff(1) += 2.0 * M_PI;
    S += weights_(col) * diff * diff.transpose();
  }
  S += R;

  // -------------------------------------------------------------------
  // UKF Update Assignment (29)
  // -------------------------------------------------------------------

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  //calculate cross correlation matrix
  Tc.setZero();
  for (int col = 0; col < 2 * n_aug_ + 1; col++)
  {
    VectorXd z_diff = Zsig.col(col) - z_pred;
    while (z_diff(1) > M_PI) z_diff(1) -= 2.0 * M_PI;
    while (z_diff(1) <- M_PI) z_diff(1) += 2.0 * M_PI;

    VectorXd x_diff = Xsig_pred_.col(col) - x_;
    while (x_diff(3) > M_PI) x_diff(3) -= 2.0 * M_PI;
    while (x_diff(3) <- M_PI) x_diff(3) += 2.0 * M_PI;
    
    Tc += weights_(col) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
  VectorXd z_diff = z - z_pred;
  while (z_diff(1) > M_PI) z_diff(1) -= 2.0 * M_PI;
  while (z_diff(1) <- M_PI) z_diff(1) += 2.0 * M_PI;
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  double nis_value = z_diff.transpose() * S.inverse() * z_diff;
  ofstream nis_file;
  nis_file.open("nis-radar.txt", ios::app);
  nis_file << nis_value;
  nis_file << "\n";
  nis_file.close();
  
  return;
}
