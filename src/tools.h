#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"
using Eigen::VectorXd;
using std::vector;
using Eigen::MatrixXd;


class Tools {
 public:
	 /*
	 // set state dimension
	 static int n_x;

	 // set augmented dimension
	 static int n_aug;

	 // Process noise standard deviation longitudinal acceleration in m/s^2
	 static double std_a;

	 // Process noise standard deviation yaw acceleration in rad/s^2
	 static double std_yawdd;

	 // define spreading parameter
	 static double lambda;
	 */

  /**
   * Constructor.
   */
  Tools();

  /**
   * Destructor.
   */
  virtual ~Tools();

  /**
   * A helper method to calculate RMSE.
   */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations,
                                const std::vector<Eigen::VectorXd> &ground_truth);



  static void UpdateState(VectorXd* x_out, MatrixXd* P_out);


};

#endif  // TOOLS_H_