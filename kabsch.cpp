// This code is released in public domain
#define STANDALONE_TEST 1
#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
using namespace Eigen;
using namespace std;

Eigen::Affine3d Find3DAffineTransformScale(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out);
Eigen::Affine3d Find3DAffineTransform(Matrix3Xd P, Matrix3Xd Q); 
void QuatfToEuler(Eigen::Vector3d& euler, const Eigen::Quaterniond& quat);
#define MAXBUFSIZE  ((int) 1000)


std::string sep = "\n----------------------------------------\n";

Matrix3Xd readMatrix(const char *filename);

// Given two sets of 3D points, find the rotation + translation + scale
// which best maps the first set to the second.
// Source: http://en.wikipedia.org/wiki/Kabsch_algorithm


// The input 3D points are stored as columns.
Eigen::Affine3d Find3DAffineTransform(Matrix3Xd P, Matrix3Xd Q) {

  // Default output
  Affine3d A;
  A.linear() = Matrix3d::Identity(3, 3);
  A.translation() = Vector3d::Zero();

  if (P.cols() != Q.cols())
    throw "Find3DAffineTransform(): input data mis-match";

  // Center the data
  Vector3d p = P.rowwise().mean();
  Vector3d q = Q.rowwise().mean();

  Matrix3Xd X = P.colwise() - p;
  Matrix3Xd Y = Q.colwise() - q;

  // SVD
  MatrixXd Cov = X*Y.transpose();
  JacobiSVD<MatrixXd> svd(Cov, ComputeThinU | ComputeThinV);

  // Find the rotation, and prevent reflections
  Matrix3d I = Matrix3d::Identity(3, 3);
  double d = (svd.matrixV()*svd.matrixU().transpose()).determinant();
  (d > 0.0) ? d = 1.0 : d = -1.0;
  I(2, 2) = d;

  Matrix3d R = svd.matrixV()*I*svd.matrixU().transpose();

  // The final transform
  A.linear() = R;
  A.translation() = q - R*p;

  return A;
}


// The input 3D points are stored as columns.
Eigen::Affine3d Find3DAffineTransformScale(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out) {

  std::cout << "performing transform" << sep;
  // Default output
  Eigen::Affine3d A;
  A.linear() = Eigen::Matrix3d::Identity(3, 3);
  A.translation() = Eigen::Vector3d::Zero();

  if (in.cols() != out.cols())
    throw "Find3DAffineTransform(): input data mis-matchusing namespace std;";

  // First find the scale, by finding the ratio of sums of some distances,
  // then bring the datasets to the same scale.using namespace std;
  double dist_in = 0, dist_out = 0;
  for (int col = 0; col < in.cols()-1; col++) {
    dist_in  += (in.col(col+1) - in.col(col)).norm();
    dist_out += (out.col(col+1) - out.col(col)).norm();
  }
  if (dist_in <= 0 || dist_out <= 0)
  {
    std::cout << "same matrix : uity rot and trans" << sep;
    return A;
  }

  std::cout << "dist_out:"<< dist_out << ":dist_in:" << dist_in  << ":" << sep;
  double scale = dist_out/dist_in;
  out /= scale;

  // Find the centroids then shift to the origin
  Eigen::Vector3d in_ctr = Eigen::Vector3d::Zero();
  Eigen::Vector3d out_ctr = Eigen::Vector3d::Zero();
  for (int col = 0; col < in.cols(); col++) {
    in_ctr  += in.col(col);
    out_ctr += out.col(col);
  }
  in_ctr /= in.cols();
  out_ctr /= out.cols();
  for (int col = 0; col < in.cols(); col++) {
    in.col(col)  -= in_ctr;
    out.col(col) -= out_ctr;
  }

  // SVD
  Eigen::MatrixXd Cov = in * out.transpose();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // Find the rotation
  double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
  if (d > 0)
    d = 1.0;
  else
    d = -1.0;
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
  I(2, 2) = d;
  Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();

  // The final transform
  std::cout << ">scale<:" <<scale << std::endl;

  A.linear() = scale * R;
  A.translation() = scale*(out_ctr - R*in_ctr);

  return A;
}

// A function to test Find3DAffineTransform()
//IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
//IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
//IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
//IOFormat HeavyFmt(FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

void TestFind3DAffineTransform(){

  // Create datasets with known transform
  Eigen::Matrix3Xd in(3, 100), out(3, 100);
  Eigen::Quaternion<double> Q(1, 3, 5, 2);
  Q.normalize();
  Eigen::Matrix3d R = Q.toRotationMatrix();




  std::cout  << "IN.row:"<< in.rows() << sep;
  std::cout << "IN.col:" << in.cols() << sep;

  double scale = 1.0;// for rigid this value is 1.0
  for (int row = 0; row < in.rows(); row++) {
    for (int col = 0; col < in.cols(); col++) {
      in(row, col) = log(2*row + 10.0)/sqrt(1.0*col + 4.0) + sqrt(col*1.0)/(row + 1.0);
    }
  }

  std::cout  << "OUT.row:"<< out.rows() << sep;
  std::cout << "OUT.col:" << out.cols() << sep;

  Eigen::Vector3d S;
  S << -5, 6, -27;
  for (int col = 0; col < in.cols(); col++)
    out.col(col) = scale*R*in.col(col) + S;

  std::cout << in << std::endl;

  std::cout << out << std::endl;

  Eigen::Affine3d A = Find3DAffineTransform(in, out);

  // See if we got the transform we expected
  if ( (scale*R-A.linear()).cwiseAbs().maxCoeff() > 1e-13 ||
       (S-A.translation()).cwiseAbs().maxCoeff() > 1e-13)
    throw "Could not determine the affine transform accurately enough";

  std::cout << "Rotation matrix" << std::endl;
  std::cout << A.linear() << std::endl;


  std::cout << "Translation matrix" << std::endl;
  std::cout << A.translation() << std::endl;



}
static void show_usage(char *argv[])
{
    std::cerr << "Usage: " << argv[0] << " <option(s)> SOURCES"
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
              << "\t-w,--worldCoordinateFile\n"
              << "\t-l,--localCoordinateFile\n"
              << "\t-t,--test\n"
              << "\t-s,--scale\n"
              << std::endl;
}
#ifdef STANDALONE_TEST
int main(int argc, char* argv[])
{
    int allowScale = 0;
    int testNumber = -1;
#ifdef DEBUG_TEST
#else
    if (argc < 2) {
        show_usage(&argv[0]);
        return 1;
    }
#endif


    const char *worldCoordinateFile,*localCoordinateFile;
    worldCoordinateFile = "world.txt";
    localCoordinateFile = "local.txt";
    allowScale = 0;
#ifdef DEBUG_TEST
#else
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        std::cout <<"arg>"<<arg<<"<"<< sep;
        if ((arg == "-h") || (arg == "--help"))
        {
            show_usage(&argv[0]);
            return 0;
        }
        else if ((arg == "-w") || (arg == "--worldCoordinateFile"))
        {
            if (i + 1 < argc)
            { // Make sure we aren't at the end of argv!
                worldCoordinateFile = argv[i+1]; // Increment 'i' so we don't get the argument as the next argv[i].
                i++;
            }
            else
            { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--worldCoordinateFile option requires one argument." << std::endl;
                return 1;
            }
        }
        else if ((arg == "-l") || (arg == "--localCoordinateFile"))
        {
            if (i + 1 < argc)
            { // Make sure we aren't at the end of argv!
                localCoordinateFile = argv[i+1]; // Increment 'i' so we don't get the argument as the next argv[i].
                i++;
            }
            else
            { // Uh-oh, there was no argument to the localCoordinate1File option.
                std::cerr << "--localCoordinateFile option requires one argument." << std::endl;
                return 1;
            }
        }
        else if ((arg == "-s") || (arg == "--scale"))
        {
            if (i + 1 < argc)
            { // Make sure we aren't at the end of argv!
                allowScale = std::atoi(argv[i+1]); // Increment 'i' so we don't get the argument as the next argv[i].
                i++;
            }
            else
            { // Uh-oh, there was no argument to the localCoordinate1File option.
                std::cerr << "--allowscale 0|1 option requires one argument." << std::endl;
                return 1;
            }
        }
        else if ((arg == "-t") || (arg == "--test"))
        {
            if (i + 1 < argc)
            { // Make sure we aren't at the end of argv!
                testNumber = std::atoi(argv[i+1]); // Increment 'i' so we don't get the argument as the next argv[i].
                i++;
            }
            else
            { // Uh-oh, there was no argument to the localCoordinate1File option.
                std::cerr << "--allowscale 0|1 option requires one argument." << std::endl;
               return 1; 
            }
        }
        else
        {
             return 0;
        }
    }
#endif

    if(testNumber > 0)
    {
       TestFind3DAffineTransform();
       return 0;
    }

    // read in world coordinate
    Matrix3Xd worldCoord; 
    Matrix3Xd localCoord;
    std::cout << "worldCoordinateFile:"<<worldCoordinateFile << sep;
    std::cout << "localCoordinateFile:"<<localCoordinateFile << sep;

    std::cout << "worldCOORD" << sep;
    worldCoord = readMatrix(worldCoordinateFile);

    // read in local coordinate
     std::cout << "localCOORD" << sep;
    localCoord = readMatrix(localCoordinateFile);
    Eigen::Affine3d A;
    if(allowScale == 1)
    {
      A = Find3DAffineTransformScale(worldCoord,localCoord );
    }
    else
    {
      A = Find3DAffineTransform(worldCoord,localCoord );
    }

    std::cout << "Rotation matrix" << std::endl;
    std::cout << A.linear() << std::endl;


    std::cout << "Translation matrix" << std::endl;
    std::cout << A.translation() << std::endl;

  std::cout << "remappedpoints from world to local" << sep;
  
  Quaternion<double> q;
  Eigen::Vector3d euler;
  Matrix3d fixed;
 

  fixed = A.linear();// change from a dynamic to fix matrix
  q = fixed; 

  QuatfToEuler(euler, q);
  double frmRadian2Degree = (double)(180.0/3.141592653589793);  
  std::cout <<"deg:pitch:"<<  euler[0]* frmRadian2Degree << ":roll:" << euler[1]*frmRadian2Degree << ":pitch:"<< euler[2]*frmRadian2Degree << sep;


 
  Matrix3Xd  remappedpoints;
  remappedpoints = (A.linear()*worldCoord).colwise() + A.translation();
  std::cout << remappedpoints << sep;

  std::cout << "remappedpoints from local to world" << sep;
  Matrix3Xd  remappedpoints2;
  Matrix3Xd  atranspose,newLocal;
  atranspose = A.linear().transpose();// the rotation matrix is singular thus to invert the matrix we transpose it
  newLocal = localCoord.colwise() - A.translation(); // perform translation then the rotation
  
  fixed = atranspose;// change from a dynamic to fix matrix
  q = fixed; 
  QuatfToEuler(euler, q);

  std::cout <<"deg:pitch:"<<  euler[0]* frmRadian2Degree << ":roll:" << euler[1]*frmRadian2Degree << ":pitch:"<< euler[2]*frmRadian2Degree << sep;


  remappedpoints2 = atranspose * newLocal;
  std::cout << remappedpoints2 << sep;


    return 0;
}
#endif
#define F_PI (3.14159265359)
void QuatfToEuler(Eigen::Vector3d& euler, const Eigen::Quaterniond& quat)
{
  float sqw = quat.w()*quat.w();
  float sqx = quat.x()*quat.x();
  float sqy = quat.y()*quat.y();
  float sqz = quat.z()*quat.z();
  float unit = sqx + sqy + sqz + sqw; // if normalised is one, otherwise is correction factor
  float test = quat.x()*quat.y() + quat.z()*quat.w();
  if (test > 0.499*unit) { // singularity at north pole
    euler[1] = 2 * atan2(quat.x(),quat.w());
    euler[2] = F_PI/2;
    euler[0] = 0;
    return;
  }
  if (test < -0.499*unit) { // singularity at south pole
    euler[1] = -2 * atan2(quat.x(),quat.w());
    euler[2] = -F_PI/2;
    euler[0] = 0;
    return;
  }
  euler[1] = atan2(2*quat.y()*quat.w()-2*quat.x()*quat.z() , sqx - sqy - sqz + sqw);
  euler[2] = asin(2*test/unit);
  euler[0] = atan2(2*quat.x()*quat.w()-2*quat.y()*quat.z() , -sqx + sqy - sqz + sqw);
}


Matrix3Xd readMatrix(const char *filename)
{
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    ifstream infile;
    infile.open(filename);
    while (! infile.eof())
    {
        string line;
        getline(infile, line);

        int temp_cols = 0;
        stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = 3;

        rows++;


    }

    infile.close();
    int dataROWS;

    rows--;
    dataROWS = rows;
    cols = rows;
    rows = 3;



    // Populate matrix with numbers.
    std::cout << "row:" << rows <<":cols:"<<cols<<sep;
    Matrix3Xd result(rows,cols);

//    std::cout << result << sep;

    for (int row = 0; row< result.rows(); row++)
    {
        for (int col = 0; col< result.cols(); col++)
        {
            double d;
            d = buff[ 3*col+row ];

            result(row,col) = d;
        }
    }

    std::cout << result << sep;
    return result;
}
