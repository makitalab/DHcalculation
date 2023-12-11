const double PI = std::acos(-1.0);

class DHParams
{
public:
	DHParams();
	~DHParams();
	std::vector<double> theta;
	std::vector<double> alpha;
	std::vector<double> a_shift;
	std::vector<double> d_shift;
	int num_joint;
private:

};

DHParams::DHParams()
{
}

DHParams::~DHParams()
{
}


/*****************************************************************
* For matrix calculation with DH parameters
******************************************************************/
class DHCalc
{
	static const int matrix_size = 4;
public:
	DHCalc();
	~DHCalc();
	static const int matsize = matrix_size;
	int matxMat(std::vector<std::vector<double>> m1, std::vector<std::vector<double>> m2, std::vector<std::vector<double>>& mdst);
	int matIdentify(std::vector<std::vector<double>>& mat);
	int matCopy(std::vector<std::vector<double>> msrc, std::vector<std::vector<double>>& mdst);
	int transZ(double d_shift, std::vector<std::vector<double>>& mat);
	int transX(double a_shift, std::vector<std::vector<double>>& mat);
	int rotZ(double theta, std::vector<std::vector<double>>& mat);
	int rotX(double alpha, std::vector<std::vector<double>>& mat);
	int dhStandard(std::vector<std::vector<double>>& htm, DHParams params, std::vector<double>& joint_pos);
	int dhModified(std::vector<std::vector<double>>& htm, DHParams params, std::vector<double>& joint_pos);
private:
private:

};

DHCalc::DHCalc()
{
}

DHCalc::~DHCalc()
{
}

int DHCalc::matxMat(std::vector<std::vector<double>> m1, std::vector<std::vector<double>> m2, std::vector<std::vector<double>>& mdst) {
	int i, j, k;

	for (i = 0; i < DHCalc::matrix_size; i++) {
		for (j = 0; j < DHCalc::matrix_size; j++) {
			mdst[i][j] = 0.0;
			for (k = 0; k < DHCalc::matrix_size; k++) {
				mdst[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}

	return 0;
}

int DHCalc::matIdentify(std::vector<std::vector<double>>& mat) {
	int i, j;
	for (i = 0; i < DHCalc::matrix_size; i++) {
		for (j = 0; j < DHCalc::matrix_size; j++) {
			if (i == j)
				mat[i][j] = 1.0;
			else
				mat[i][j] = 0.0;
		}
	}
	return 0;
}

int DHCalc::matCopy(std::vector<std::vector<double>> msrc, std::vector<std::vector<double>>& mdst) {
	for (int i = 0; i < DHCalc::matrix_size; i++) {
		for (int j = 0; j < DHCalc::matrix_size; j++) {
			mdst[i][j] = msrc[i][j];
		}
	}
	return 0;
}

int DHCalc::transZ(double d_shift, std::vector<std::vector<double>>& mat) {

	DHCalc::matIdentify(mat);
	mat[2][3] = d_shift;

	return 0;
}

int DHCalc::transX(double a_shift, std::vector<std::vector<double>>& mat) {

	DHCalc::matIdentify(mat);
	mat[0][3] = a_shift;

	return 0;
}

int DHCalc::rotZ(double theta, std::vector<std::vector<double>>& mat) {

	DHCalc::matIdentify(mat);
	mat[0][0] = cos(theta);
	mat[0][1] = -1.0 * sin(theta);
	mat[1][0] = sin(theta);
	mat[1][1] = cos(theta);

	return 0;
}

int DHCalc::rotX(double alpha, std::vector<std::vector<double>>& mat) {

	DHCalc::matIdentify(mat);
	mat[1][1] = cos(alpha);
	mat[1][2] = -1.0 * sin(alpha);
	mat[2][1] = sin(alpha);
	mat[2][2] = cos(alpha);

	return 0;
}

int DHCalc::dhStandard(std::vector<std::vector<double>>& htm, DHParams params, std::vector<double>& joint_pos) {
	std::vector<std::vector<double>> tmpm1(DHCalc::matsize, std::vector<double>(DHCalc::matsize)), tmpm2(DHCalc::matsize, std::vector<double>(DHCalc::matsize));
	for (int i = 0; i < params.num_joint; i++) {
		DHCalc::matCopy(htm, tmpm1);
		DHCalc::rotZ(params.theta[i], tmpm2);
		DHCalc::matxMat(tmpm1, tmpm2, htm);
		DHCalc::matCopy(htm, tmpm1);
		DHCalc::transZ(params.d_shift[i], tmpm2);
		DHCalc::matxMat(tmpm1, tmpm2, htm);
		DHCalc::matCopy(htm, tmpm1);
		DHCalc::rotX(params.alpha[i], tmpm2);
		DHCalc::matxMat(tmpm1, tmpm2, htm);
		DHCalc::matCopy(htm, tmpm1);
		DHCalc::transX(params.a_shift[i], tmpm2);
		DHCalc::matxMat(tmpm1, tmpm2, htm);
		for (int k = 0; k < DHCalc::matsize - 1; k++) {
			joint_pos.push_back(htm[k][DHCalc::matsize - 1]);
		}
	}
	return 0;
}

int DHCalc::dhModified(std::vector<std::vector<double>>& htm, DHParams params, std::vector<double>& joint_pos) {
	std::vector<std::vector<double>> tmpm1(DHCalc::matsize, std::vector<double>(DHCalc::matsize)), tmpm2(DHCalc::matsize, std::vector<double>(DHCalc::matsize));
	joint_pos.clear();
	for (int i = 0; i < params.num_joint; i++) {
		DHCalc::matCopy(htm, tmpm1);
		DHCalc::rotX(params.alpha[i], tmpm2);
		DHCalc::matxMat(tmpm1, tmpm2, htm);
		DHCalc::matCopy(htm, tmpm1);
		DHCalc::transX(params.a_shift[i], tmpm2);
		DHCalc::matxMat(tmpm1, tmpm2, htm);
		DHCalc::matCopy(htm, tmpm1);
		DHCalc::rotZ(params.theta[i], tmpm2);
		DHCalc::matxMat(tmpm1, tmpm2, htm);
		DHCalc::matCopy(htm, tmpm1);
		DHCalc::transZ(params.d_shift[i], tmpm2);
		DHCalc::matxMat(tmpm1, tmpm2, htm);
		for (int k = 0; k < DHCalc::matsize - 1; k++) {
			joint_pos.push_back(htm[k][DHCalc::matsize - 1]);
		}
	}
	return 0;
}
