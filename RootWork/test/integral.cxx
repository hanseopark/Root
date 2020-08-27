#include<iostream>
#include<math.h>

double func(double x){
	
	double Func = 0;
	Func = exp(-x*x/10);

	return Func;
}

double integral(double start, double end){
	int num = 100;
	double ans = 0;
	double dx = ((end-start))/(num);

	for(int i =0; i< num; i++){
	double x_i = start + i*dx;
	ans += func(x_i)*dx;
	}

	return ans;

}

int main(){

	std::cout << integral(1,10) <<  std::endl;
}
