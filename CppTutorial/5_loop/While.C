#include <iostream>
using namespace std;

int main(){
	int start = 0;
	int end = 0;
	int sum = 0;
	
	cout << "Start = ";
	cin >> start;
	cout << "Start = ";
	cin >> end;

	while(start <= end){
		sum += start;
		start ++;
	}

	cout << "Result = " << sum << endl;
	return 0;
}
