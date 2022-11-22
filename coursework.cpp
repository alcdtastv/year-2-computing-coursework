#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

//Runge-Kutta method for double pendulum
//This program takes as input the masses, lengths and initial angles of a double pendulum system and solves the equations of motion specified in equation 1
//of the handout through the use of the fourth order Runge-Kutta method.
//The initial conditions, timestep and final time have to be specified in a file called "parameters.txt"
//The output consists of a file called "output.txt" that contains the time and the x and y coordinates for both of the masses at each iteration

vector<long double> input(string name) {
	vector<long double> initialconditions(8);
	
	//this function takes as input the name of the input file and outputs a vector containing the initial conditions
	
	ifstream vIn(name, ios::in);  //creating an import variable
	
	if (vIn.good()) {      //checking that the file opened correctly
		
		while (true) {    //importing the initial conditions
			vIn >> initialconditions[0]
				>> initialconditions[1] 
				>> initialconditions[2] 
				>> initialconditions[3] 
				>> initialconditions[4] 
				>> initialconditions[5] 
				>> initialconditions[6] 
				>> initialconditions[7];
				
			if (vIn.eof()) {
				break;
			}
		}
		vIn.close();  //closing the file
	}
	else {
		cout << "File not opened correctly" << endl;  //outputting an error message in case it is impossible to open the file
	}
	return initialconditions;
} 

long double acceleration1(const long double m1, const long double m2, const long double l1, const long double l2,
						  const long double g, long double t1, long double t2, long double v1, long double v2) 
{
	
	//this function takes as input the masses, the lengths, the acceleration of gravity and the angles and angular velocities of 
	//the pendulums to compute the acceleration of the first mass 
	
	long double a1 = (-m2 * cos(t1-t2) * l1 * pow(v1,2) * sin(t1-t2)
				 +m2 * cos(t1-t2) * g * sin(t2)
				 -m2 * l2 * pow(v2,2) * sin(t1-t2)
				 -(m1+m2) * g * sin(t1))/
				 (l1 * (m1 + m2 - m2 * pow(cos(t1-t2),2)));
				 
	return a1;
}

long double acceleration2(const long double m1, const long double m2, const long double l1, const long double l2,
						  const long double g,long double t1, long double t2, long double v1, long double v2) 
{
	
	//this function takes as input the masses, the lengths, the acceleration of gravity and the angles and angular velocities of 
	//the pendulums to compute the acceleration of the second mass 
	
	long double a2 = ((m1 + m2) * l1 * pow(v1,2) * sin(t1-t2)
				 + pow(v2,2) * sin(t1-t2) * cos(t1-t2) *m2 *l2
				 + (m1 + m2) * cos(t1-t2) * g * sin(t1)
				 - (m1 + m2) * g * sin(t2))/
				 (l2 * (m1 + m2 * pow(sin(t1-t2), 2)));
				
	return a2;
}

long double velocity1(long double v1) {
	
	//this function takes as input the velocity of the first mass and returns it
	//it has been included only for consistency with the matlab code and ease of debugging
	
	long double function1 = v1;
	
	return function1;
}

long double velocity2(long double v2) {
	
	//this function takes as input the velocity of the second mass and returns it
	//it has been included only for consistency with the matlab code and ease of debugging
	
	long double function2 = v2;
	
	return function2;
}

void rungekuttapendulum(vector<long double> initialconditions) {
	
	//this function calculates the value of the four coefficients required to compute the fourth order runge-kutta approximation
	//of the double pendulum equations to subsequently calculate and output the angles of the two pendulums at each timestep
	
	//the input for this program is the vector of initial conditions produced by the function "input"
	
	//assigning the initial conditions to new variables to improve readability
	const long double m1 = initialconditions[0];
	const long double m2 = initialconditions[1];
	const long double l1 = initialconditions[2];
	const long double l2 = initialconditions[3];
	long double t1 = initialconditions[4];
	long double t2 = initialconditions[5];
	const long double timestep = initialconditions[6];
	const long double finaltime = initialconditions[7];
	long double v1 = 0;
	long double v2 = 0;
	const long double g = 9.81;
	
	//declaring variables for iteration
	long double k1t1, k2t1, k3t1, k4t1, k1t2, k2t2, k3t2, k4t2, k1v1, k2v1, k3v1, k4v1, k1v2, k2v2, k3v2, k4v2, x1, y1, x2, y2, t;
	
	//initialising output variable
	ofstream vOut("output.txt", ios::out | ios::trunc);
	
	//computing the initial positions of the masses from the angles of the pendulums
	x1 = l1 * sin(t1);
	y1 = -(l1 * cos(t1));
	x2 = l1 * sin(t1) + l2 * sin(t2);
	y2 = -(l1 * cos(t1) + l2 * cos(t2));
	
	t = 0;  //setting the initial time
	
	//writing to the output file the header, the initial time and the initial positions of the masses 
	vOut << setw(20) << setprecision(10) << "Time"
		 << setw(20) << setprecision(10) << "x1"
		 << setw(20) << setprecision(10) << "y1"
		 << setw(20) << setprecision(10) << "x2"
		 << setw(20) << setprecision(10) << "y2" << endl
		 << setw(20) << setprecision(10) << t
		 << setw(20) << setprecision(10) << x1
		 << setw(20) << setprecision(10) << y1
		 << setw(20) << setprecision(10) << x2
		 << setw(20) << setprecision(10) << y2 << endl;
	
	if (vOut.good()) {
		//performing runge kutta iteration and outputting the values of x1, y1, x2 and y2
		for (int i=0; i<finaltime/timestep; i++) {
			
			k1t1 = timestep * velocity1(v1);
			k1t2 = timestep * velocity2(v2);
			k1v1 = timestep * acceleration1(m1, m2, l1, l2, g, t1, t2, v1, v2);
			k1v2 = timestep * acceleration2(m1, m2, l1, l2, g, t1, t2, v1, v2);
			
			k2t1 = timestep * velocity1(v1 + k1v1/2);
			k2t2 = timestep * velocity2(v2 + k1v2/2);
			k2v1 = timestep * acceleration1(m1, m2, l1, l2, g, t1 + k1t1/2, t2 + k1t2/2, v1 + k1v1/2, v2 + k1v2/2);
			k2v2 = timestep * acceleration2(m1, m2, l1, l2, g, t1 + k1t1/2, t2 + k1t2/2, v1 + k1v1/2, v2 + k1v2/2);
			
			k3t1 = timestep * velocity1(v1 + k2v1/2);
			k3t2 = timestep * velocity2(v2 + k2v2/2);
			k3v1 = timestep * acceleration1(m1, m2, l1, l2, g, t1 + k2t1/2, t2 + k2t2/2, v1 + k2v1/2, v2 + k2v2/2);
			k3v2 = timestep * acceleration2(m1, m2, l1, l2, g, t1 + k2t1/2, t2 + k2t2/2, v1 + k2v1/2, v2 + k2v2/2);
			
			k4t1 = timestep * velocity1(v1 + k3v1);
			k4t2 = timestep * velocity2(v2 + k3v2);
			k4v1 = timestep * acceleration1(m1, m2, l1, l2, g, t1 + k3t1, t2 + k3t2, v1 + k3v1, v2 + k3v2);
			k4v2 = timestep * acceleration2(m1, m2, l1, l2, g, t1 + k3t1, t2 + k3t2, v1 + k3v1, v2 + k3v2);
			
			t1 += k1t1/6 + k2t1/3 + k3t1/3 + k4t1/6;
			t2 += k1t2/6 + k2t2/3 + k3t2/3 + k4t2/6;
			v1 += k1v1/6 + k2v1/3 + k3v1/3 + k4v1/6;
			v2 += k1v2/6 + k2v2/3 + k3v2/3 + k4v2/6;
			
			x1 = l1 * sin(t1);
			y1 = -(l1 * cos(t1));
			x2 = l1 * sin(t1) + l2 * sin(t2);
			y2 = -(l1 * cos(t1) + l2 * cos(t2));
			
			t += timestep;
			
			vOut << setw(20) << setprecision(10) << t
				 << setw(20) << setprecision(10) << x1
				 << setw(20) << setprecision(10) << y1
				 << setw(20) << setprecision(10) << x2
				 << setw(20) << setprecision(10) << y2 << endl;
		}
	}
}

int main () {
	rungekuttapendulum(input("parameters.txt"));
}