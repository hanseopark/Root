#include <iostream>
using namespace std;

int sum(int a, int b){ //受け取った引数をint型のa, bという名前で記憶する
	int sum = 0;
	sum = a + b;
	return sum;
}

int main(){
	cout << sum(1, 3) << endl; //sum(int a, int b)という関数に、a=1、b=3を代入する
	return 0;
}
