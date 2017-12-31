#include <iostream>
#include "ukf.h"

UKF::UKF() {
  //TODO Auto-generated constructor stub
  Init();
}

UKF::~UKF() {
  //TODO Auto-generated destructor stub
}

void UKF::Init() {

}

/*******************************************************************************
* Programming assignment functions: 
*******************************************************************************/

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out, MatrixXd Xsig_pred) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  
  //create vector for predicted state
  VectorXd x = VectorXd(n_x);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x, n_x);


/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  // set weights
  double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
    x = x+ weights(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }


/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;
}