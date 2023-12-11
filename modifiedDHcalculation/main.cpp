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
1) Select your robot and its DH parameters
2) Make a CSV file containing the DH parameters with the following format
	text: number of joints
	value
	text: alpha
	values
	text: a_shift
	values
	text: d_shift
	values

************************************************/

/*****************************************
Select your robot and its DH parameters
******************************************/
const bool xArm6ModifiedDH = false;
const bool xArm6StandardDH = true;

int main(void) {
	DHCalc dh;
	double t2_offset = 0.0;
	double t3_offset = 0.0;
	bool dh_standard = false;
	bool dh_modified = false;
	const int dim = 3;
	DHParams params;
	std::vector<std::vector<double>> htm(dh.matsize, std::vector<double>(dh.matsize));
	std::vector<std::vector<double>> tmpm1(dh.matsize, std::vector<double>(dh.matsize)), tmpm2(dh.matsize, std::vector<double>(dh.matsize));
	int i, k, count;
	std::string filename;
	
	/*************************************
	Read the file of DH parameters 
	**************************************/
	if (xArm6ModifiedDH) {
		filename = "xArm6ModifiedDHParameters.csv";
		readParams(params, filename);
		
		dh_modified = true;
		t2_offset = -1.0 * std::atan2(284.5, 53.5);
		t3_offset = -1.0 * t2_offset;
	}
	else if (xArm6StandardDH) {
		filename = "xArm6StandardDHParameters.csv";
		readParams(params, filename);
		
		dh_standard = true;
		t2_offset = -1.0 * std::atan2(284.5, 53.5);
		t3_offset = -1.0 * t2_offset;
	}

	/*************************************
	Set each joint variable
	**************************************/
	params.theta.clear();
	params.theta.push_back(0.0);
	params.theta.push_back(t2_offset);
	params.theta.push_back(t3_offset);
	params.theta.push_back(0.0);
	params.theta.push_back(0.0);
	params.theta.push_back(0.0);

	std::vector<double> j_pos(dim*params.num_joint);
	dh.matIdentify(htm);
	if (dh_modified) {
		dh.dhModified(htm, params, j_pos);
	}
	else if (dh_standard) {
		dh.dhStandard(htm, params, j_pos);
	}

	/*************************************
	Print every joint position
	**************************************/
	for (i = 0; i < params.num_joint; i++) {
		printf_s("joint_position %d (", i + 1);
		for (k = 0; k < dim; k++) {
			printf_s("%7.3f", j_pos[i * dim + k]);
			if (k == dim -1) {
				printf_s(")\n");
			}
			else {
				printf_s(", ");
			}
		}
	}

	/***************************************************************
	Print homogeneous translation matrix of end of the manipulator from the base
	***************************************************************/
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