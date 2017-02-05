Eigen::Affine3d Find3DAffineTransformScale(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out);
Eigen::Affine3d Find3DAffineTransform(Matrix3Xd P, Matrix3Xd Q);
void QuatfToEuler(Eigen::Vector3d& euler, const Eigen::Quaterniond& quat);
Eigen::Matrix3Xd transform_frame1_to_frame2(Eigen::Affine3d A , Eigen::Matrix3Xd coordIn,Eigen::Matrix3Xd coordOut);
