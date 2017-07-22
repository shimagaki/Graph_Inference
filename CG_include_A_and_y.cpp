#include <iostream>
#include <fstream> 
#include <math.h>
#include <Eigen/Dense>
#include <vector>  //std::back_inserter
using namespace std;
using namespace Eigen;

int N ;
int n_iter = 100;
double epsilon_pAp = 1e-10;
double epsilon_rr = 1e-10;
/**** Ax = y ****/
MatrixXd A;
VectorXd y,x0,x,p0,p,r0,r;
double a,b;

vector<double> temp_y;
vector<double> temp_A;
vector<double> temp_x_fixed;
int size_y;
int size_A;
int size_x_fix;

void read_A_and_y(){

	double data_A,data_y,data_x_fix;

	string fname_A = "A_mat.dat";
	ifstream file_A(fname_A);
	string fname_y = "y_vec.dat";
	ifstream file_y(fname_y);
	string fname_x_fix = "config_S_step0.dat";
	ifstream file_x_fix(fname_x_fix);
	
	if(!file_A.is_open()){
		cout<< "Error: Couldn't open your file " + fname_A <<endl;	
		exit(1);
	} 
	if(!file_y.is_open()){
		cout<< "Error: Couldn't open your file " + fname_y <<endl;	
		exit(1);
	}
	if(!file_x_fix.is_open()){
		cout<< "Error: Couldn't open your file " + fname_x_fix <<endl;	
		exit(1);
	}

	while(file_A){
		file_A >> data_A;
		temp_A.push_back(data_A);
	} 
	while(file_y){
		file_y >> data_y;
		temp_y.push_back(data_y);
	} 
	while(file_x_fix){
		file_x_fix >> data_x_fix;
		temp_x_fixed.push_back(data_x_fix);
	} 
	size_y = temp_y.size()-1; // I don't know why last elemets would be a fake data.
	size_A = temp_A.size()-1;
	size_x_fix = temp_x_fixed.size()-1;
	//Update of N
	/*
	for(int i=0;i<size_y;++i)
	       cout<<i <<","<< temp_y[i]<<endl;	
	for(int j=0;j<size_x_fix;++j)
	       cout<<j <<","<< temp_x_fixed[j]<<endl;	
	*/
	N = size_y;
	}

void init_A_and_y(){
	y = VectorXd::Zero(N);
	A = MatrixXd::Zero(N,N);
}

void set_A_and_y(){

	for(int i=0;i<N;++i){
		y[i] = temp_y[i];
	for(int j=0;j<N;++j){
		A(i,j) = temp_A[i*N+j]; // should i use ato()? <= it seems not necessary. 
	}}

}

void init_vector(){
	x0= VectorXd::Zero(N);
	r0 = y;
	p0 = y;
}

void output_solution(){
	string fname = "config_S_solution.dat";
	ofstream file(fname);
	for(int j=0;j< size_x_fix-size_y ;++j){
		file<<temp_x_fixed[j]<<endl;
	}
	for(int i=0;i<N;++i){
		file<<x[i]<<endl;
	}
	file.close();
}

int main(int argc, char const* argv[]){
	read_A_and_y();
	init_A_and_y();
	set_A_and_y();
	//cout << "A = \n" << A << endl;
	//cout << "y = \n" << y << endl;
	init_vector();
	double r0r0, p0Ap0, rr;	
	for(int n=0;n<n_iter;++n){
		r0r0 = r0.dot(r0);
		p0Ap0 = p0.dot(A * p0);
		
		if(abs(p0Ap0)<epsilon_pAp){
			cout<< "A too small variables appears in p0.dot(A * p0) " <<endl;
			break;
		}	
		
		a = r0r0/ p0Ap0; 
		x = x0 + a * p0;
		r = r0 - a * (A*p0);
		rr = r.dot(r);
		
		cout<< n << " " << rr << endl;	
		if(rr < epsilon_rr){
			cout<< "The residual |Ax-b| is sufficiently small. " <<endl;
			break;
		}	
		
		b = rr / r0r0;
	       	p = r + b*p0;	
		// Update variables. 
		x0 = x; p0 = p; r0 = r;	
	}
	cout << "x = \n" << x << endl;
	output_solution();
	return 0;
}
