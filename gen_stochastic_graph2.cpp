#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <Eigen/Dense>
#include <vector>  
using namespace std;
using namespace Eigen;

//(N_l,N_u)=(3,7)
int N_l = 10;	//labeled 
int N_u= 90; 	//unlabeled 
int N = N_l+N_u;
int q = 3;
double p_conect = 0.5;
//double T = 0.01;
double T = 1.0;
double epsilon = 1e-3;
MatrixXd W = MatrixXd::Zero(N,N);
MatrixXd D = MatrixXd::Zero(N,N);
MatrixXd A(N_u,N_u);
VectorXd y(N_u);
vector< int > S(N);	//spin variables at each lattice site
vector< vector< int > > Neighbor(N);	//spin variables at each lattice site

inline double std_rand(){
	return rand() / (RAND_MAX + 1.0);
}

void set_randmat(double p){
	double r;
	for(int i=0;i<N;++i){
	for(int j=i+1; j<N; ++j){
		r = std_rand();
		if(r < p){
			W(i,j) = r; 
			W(j,i) = r;
		}
	}}
}

// i=0,1,.,l => Fixed
// i=l,..,l+u-1 => Variables, which should be optimize. 
void init_S(){
	for(int i=0;i<N_l;++i){
		S[i] = int(q*std_rand()); 
	}
	for(int j=N_l;j<N_l+N_u;++j){
		S[j] = int(q*std_rand()); 
	}
}

// A := D_uu - W_uu
void set_A_and_y(){
	VectorXd D(N);
	VectorXd one = VectorXd::Ones(N);
	VectorXd v(N_l);
	MatrixXd W_ul(N_u,N_l);
	
	D = W * one;
	
	for(int i=0;i<N_u;++i){
		A(i,i) += D[N_l+i]; 
		for(int j=0; j<N_u; ++j){
			A(i,j) -= W(N_l+i,N_l+j); 
	}}
	
	for(int i=0;i<N_l;++i){
		v[i] = S[i];
		for(int j=0;j<N_u;++j){
			W_ul(j,i) = W(j+N_l, i);	
		}	
	}
	y = W_ul * v;
}
void output_A_y(){
	string fname_A = "A_mat.dat";
	ofstream file_A(fname_A);
	string fname_y = "y_vec.dat";
	ofstream file_y(fname_y);
	for(int i=0;i<N_u;++i){
	for(int j=0;j<N_u;++j){
		file_A << A(i,j) <<" ";	
	}file_A<<endl;
		file_y << y[i]  << endl;	
	}
	file_A.close();
	file_y.close();
}

void set_neighbor(){
	double r;
	for(int i=0;i<N;++i){
	for(int j=0;j<N;++j){
		if(abs(W(i,j))>epsilon){
			//Want j (>0) appends to Neighbor[i]
			Neighbor[i].push_back(j);
		}
	}
	}
}

void output_W(){
	string fname = "adjacent.dat";
	ofstream file(fname);
	for(int i=0;i<N;++i){
	for(int j=0;j<N;++j){
		file<<W(i,j)<<" ";
	}file<<endl;
	}
	file.close();
}

void output_S(int step){
	string fname = "config_S_step"+to_string(step)+".dat";
	ofstream file(fname);
	for(int i=0;i<N;++i){
		file<<S[i]<<endl;
	}
	file.close();
}

double E(int i)
{

    // contribution to energy from 4 bonds at site (x,y)
    double E = 0;
    for(int j=0; j< Neighbor[i].size() ; ++j){
   	E += W(i,j)*pow((S[i]-S[j]),2); 
    }
    return E;
}

double E()
{
    // total energy
    double E_sum = 0;
    for (int i=0;i<N; ++i)
	    E_sum += E(i);
    return E_sum / 2;
}

int update_S(int state){
	return int(q*std_rand()); 
}

bool metropolis_step(int i){     
    // find the local energy at site (x,y)
    double E_i = E(i);

    // save the value of spin at (x,y)
    double S_save = S[i];

    // trial Metropolis step for spin at site (x,y)
    S[i]=update_S(S[i]);

    double E_i_trial = E(i);

    // Metropolis test
    bool accepted = false;
    double dE = E_i_trial - E_i;
    double w = exp(- dE / T);
    if (w > std_rand())
        accepted = true;
    else
        S[i] = S_save;
    return accepted;
}

int monte_carlo_sweep()
{
    int acceptances = 0;
    // randomized sweep
    for (int step = 0; step < N_u ; ++step) {
	// i is chosed from unlabeled set.
        int i = int(N_u * std_rand()) + N_l;
        if (metropolis_step(i))
            ++acceptances;
    }
    return acceptances;
}

int main(int argc, char const* argv[]){
	
	set_randmat(p_conect); //initialize of graph.
	init_S();
	output_W();
	
	set_A_and_y();
	output_A_y();
	
	set_neighbor();	
	//int equilibration_steps = 0;
	int equilibration_steps = 5000;
	//double dT = (T-0.000001)/equilibration_steps;
	for (int step = 0; step <= equilibration_steps; ++step) {
		monte_carlo_sweep();
		cout << step << "  " << E() << endl;
		if( step==0 || step==equilibration_steps ){
			output_S(step);
		}
	
		T = 1.0/pow(1+10*step,0.5); 
	}
	return 0;
}
