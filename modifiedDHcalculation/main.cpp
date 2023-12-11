// -*- coding:utf-8 -*-
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stddef.h>
#include <vector>
#include "./myconfig.hpp"

/***********************************************
Calculate Homogeneous Translate Matrix with standard/modified Denavit-Hartenberg (DH) parameters

************************************************/

const int xArm6ModifiedDH = 1;
const int xArm6StandardDH = 0;

int main(void) {
	DHCalc dh;
	int num_joint = 0;
	double t2_offset = 0.0;
	double t3_offset = 0.0;
	bool dh_standard = false;
	bool dh_modified = false;
/*
	if (xArm6ModifiedDH) {
		num_joint = 6;
		t2_offset = -1.0 * atan2(284.5, 53.5);
		t3_offset = -1.0 * t2_offset;
		dh_modified = 1;
	}
	else if (xArm6StandardDH) {
		num_joint = 6;
		t2_offset = -1.0 * atan2(284.5, 53.5);
		t3_offset = -1.0 * t2_offset;
		dh_standard = 1;
	}
	else {
		num_joint = 0;
	}
	*/
	const int dim = 3;
	DHParams params;
	std::vector<std::vector<double>> htm(dh.matsize, std::vector<double>(dh.matsize));
	std::vector<std::vector<double>> tmpm1(dh.matsize, std::vector<double>(dh.matsize)), tmpm2(dh.matsize, std::vector<double>(dh.matsize));
	int i, k, count;
	std::ifstream fs{};
	std::string str;

	if (xArm6ModifiedDH) {
		fs.open("xArm6ModifiedDHParameters.csv");
		if (fs.fail()) {
			std::cerr << "Failed to open DH parameter file." << std::endl;
			return -1;
		}
		else {
			std::getline(fs, str); // str == "num_joint"
			std::cout << "reading" + str << std::endl;
			std::getline(fs, str); // str == the value of num_joint
			num_joint = atoi(str.c_str());
			std::getline(fs, str); // str == "alpha"
			std::cout << "reading" + str << std::endl;
			for (i = 0; i < num_joint; i++) {
				std::getline(fs, str);
				params.alpha.push_back(atof(str.c_str()));
			}
			std::getline(fs, str); // str == "a_shift"
			std::cout << "reading" + str << std::endl;
			for (i = 0; i < num_joint; i++) {
				std::getline(fs, str);
				params.a_shift.push_back(atof(str.c_str()));
			}
			std::getline(fs, str); // str == "d_shift"
			std::cout << "reading" + str << std::endl;
			for (i = 0; i < num_joint; i++) {
				std::getline(fs, str);
				params.d_shift.push_back(atof(str.c_str()));
			}
		}
		dh_modified = true;
		t2_offset = -1.0 * std::atan2(284.5, 53.5);
		t3_offset = -1.0 * t2_offset;
	}
	else if (xArm6StandardDH) {
		dh_standard = true;
	}

	params.theta.push_back(0.0);
	params.theta.push_back(t2_offset);
	params.theta.push_back(t3_offset);
	params.theta.push_back(0.0);
	params.theta.push_back(0.0);
	params.theta.push_back(0.0);

	std::vector<double> j_pos(dim*params.num_joint);
	dh.matIdentify(htm);
	for (i = 0; i < num_joint; i++) {
		if (dh_modified) {
			dh.matCopy(htm, tmpm1);
			dh.rotX(params.alpha[i], tmpm2);
			dh.matxMat(tmpm1, tmpm2, htm);
			dh.matCopy(htm, tmpm1);
			dh.transX(params.a_shift[i], tmpm2);
			dh.matxMat(tmpm1, tmpm2, htm);
			dh.matCopy(htm, tmpm1);
			dh.rotZ(params.theta[i], tmpm2);
			dh.matxMat(tmpm1, tmpm2, htm);
			dh.matCopy(htm, tmpm1);
			dh.transZ(params.d_shift[i], tmpm2);
			dh.matxMat(tmpm1, tmpm2, htm);
		}
		else if (dh_standard) {
			dh.matCopy(htm, tmpm1);
			dh.rotZ(params.theta[i], tmpm2);
			dh.matxMat(tmpm1, tmpm2, htm);
			dh.matCopy(htm, tmpm1);
			dh.transZ(params.d_shift[i], tmpm2);
			dh.matxMat(tmpm1, tmpm2, htm);
			dh.matCopy(htm, tmpm1);
			dh.rotX(params.alpha[i], tmpm2);
			dh.matxMat(tmpm1, tmpm2, htm);
			dh.matCopy(htm, tmpm1);
			dh.transX(params.a_shift[i], tmpm2);
			dh.matxMat(tmpm1, tmpm2, htm);
		}
		for (k = 0; k < dim; k++) {
			j_pos.push_back(htm[k][3]);
		}
		printf_s("joint_position %d ( %7.3f, %7.3f, %7.3f )\n", i, j_pos[i*dim + 0], j_pos[i * dim + 1], j_pos[i * dim + 2]);
	}
	
	printf_s("Homogeneous Translation Matrix is\n");
	for (i = 0; i < dh.matsize; i++) {
		for (k = 0; k < dh.matsize; k++) {
			printf_s("%7.3f", htm[i][k]);
			if (k == dh.matsize - 1)
				printf_s("\n");
			else
				printf_s(", ");
		}
	}
	return 0;
}