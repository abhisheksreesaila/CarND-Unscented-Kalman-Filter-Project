#include "ukf.h"
#include <iostream>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  previous_timestamp_ = 0;

  is_initialized_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  

      // State dimension
    n_x_ = 5;

    // Augmented state dimension
    n_aug_ = 7;

    // Sigma point spreading parameter
    lambda_ = 3-n_x_;

    //create vector for weights
    weights_ = VectorXd(2 * n_aug_ + 1);

    // state covariance matrix P_
    P_ = MatrixXd(n_x_, n_x_);
    P_ <<      1, 0, 0, 0, 0,
               0, 1, 0, 0, 0,
               0, 0, 1, 0, 0,
               0, 0, 0, 1, 0,
               0, 0, 0, 0, 1;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
 
    if (!is_initialized_) {
    
      x_ << 1, 1, 1, 1, 1;

      if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {

        x_(0) = meas_package.raw_measurements_(0);
        x_(1) = meas_package.raw_measurements_(1);

      }
      else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
  
        float ro = meas_package.raw_measurements_(0);
        float phi = meas_package.raw_measurements_(1);
        x_(0) = ro * cos(phi);
        x_(1) = ro * sin(phi);
      }

      previous_timestamp_ = meas_package.timestamp_;
      is_initialized_ = true;

      return;
    }

    /*****************************************************************************
    *  Prediction
    ****************************************************************************/
    
    float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	
    previous_timestamp_ = meas_package.timestamp_;

    Prediction(dt);

    /*****************************************************************************
    *  Update
    ****************************************************************************/

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      UpdateLidar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
    }
  

}

void UKF::Prediction(double delta_t) {


  /**************** creating augmented sigma matrices ******************/

  // create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  
   MatrixXd Q = MatrixXd(2,2);
   Q  << (std_a_*std_a_), 0, 0, (std_yawdd_*std_yawdd_);
   
   MatrixXd A = MatrixXd(5,2);
   A  << 0,0,0,0,0,0,0,0,0,0;
   
   MatrixXd B = MatrixXd(2,5);
   B  << 0,0,0,0,0,0,0,0,0,0;
  
  
  // create augmented covariance matrix
  P_aug.topLeftCorner(5,5) = P_;
  P_aug.topRightCorner(5,2) = A;
  P_aug.bottomLeftCorner(2,5) = B;
  P_aug.bottomRightCorner(2,2) = Q;
  

  // create square root matrix
  MatrixXd SQ = P_aug.llt().matrixL();  

  Xsig_aug.col(0) = x_aug;
  
  // create augmented sigma points
  for (int i = 0; i < n_aug_; ++i)
  {
  Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_)*SQ.col(i);
  Xsig_aug.col(i+1+n_aug_) =  x_aug  - sqrt(lambda_+n_aug_)*SQ.col(i);
  
  }
  
 
  /**********************Predict Sigma Points**********************************/

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // predict sigma points
       for(int i = 0; i< (2 * n_aug_ + 1); i++){

            VectorXd colValues =  Xsig_aug.col(i);
            double v = colValues(2); 
            double a = colValues(3); 
            double b = colValues(4);
            double c = b*delta_t;
            double nua = colValues(5); 
            double nub = colValues(6); 
            VectorXd colValues1 = VectorXd(5);
            VectorXd colValues2 = VectorXd(5);
            VectorXd colValues3 = VectorXd(5);
            colValues1 = colValues.head(5);
            if (b!=0)
            {
              colValues2(0) = (v/b) * (sin(a + c) - sin(a));
              colValues2(1)  = (v/b) * (-cos(a + c) + cos(a));
              colValues2(2) = 0;
              colValues2(3) = c;  
              colValues2(4) = 0;
    
              colValues3(0) = 0.5 * (delta_t*delta_t) * cos(a) * nua;
              colValues3(1) = 0.5  * (delta_t*delta_t) * sin(a) * nua;
              colValues3(2) = delta_t * nua;
              colValues3(3) = 0.5  * (delta_t*delta_t) * nub; 
              colValues3(4) = delta_t * nub;
         
              
            } 
              else {
                  
                  
              colValues2(0) = v * cos(a) * delta_t;
              colValues2(1) = v * sin(a) * delta_t;
              colValues2(2) = 0;
              colValues2(3) = 0;  
              colValues2(4) = 0;
            
              colValues3(0) = 0.5 * (delta_t*delta_t) * cos(a) *  nua ;
              colValues3(1) = 0.5 * (delta_t*delta_t) * sin(a) * nua;
              colValues3(2) = delta_t * nua;
              colValues3(3) = 0.5 * (delta_t*delta_t) * nub; 
              colValues3(4) = delta_t * nub;
              
                  
              }
          
          Xsig_pred_.col(i)=colValues1 + colValues2 + colValues3; 
          
       }
  
  /******************Convert Predicted Sigma Points to Mean/Covariance***********************/
  // set weights
  for( int i=0; i < weights_.size(); i++)
 {
    if(i == 0)
        weights_(i) = lambda_/(lambda_ + n_aug_);
    else
        weights_(i) = 1/(2 * (lambda_ + n_aug_));
 }  
 
  // predict state mean
   for(int i=0; i < n_x_; i++)
   {
    double totalSum = 0;
           for(int j=0; j< 2 * n_aug_ + 1; j++)
             totalSum = totalSum + (weights_(j) * Xsig_pred_(i, j));
     x_(i) = totalSum;
   }
   
  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

     

 
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {

    //extract measurement as VectorXd
  VectorXd z = meas_package.raw_measurements_;

  //set measurement dimension, lidar can measure p_x and p_y
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    // measurement model
    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();



}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  

  //extract measurement as VectorXd
  VectorXd z = meas_package.raw_measurements_;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1, i) = atan2(p_y, p_x);                                 //phi
    Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { 
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_,                       0,                     0,
                         0, std_radphi_*std_radphi_,                     0,
                         0,                       0, std_radrd_*std_radrd_;
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  /*****************************************************************************
  *  UKF Update for Radar
  ****************************************************************************/
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

}