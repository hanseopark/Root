#include <iostream>
#include <string>
using namespace std;

int main(){
	int start = 0;
	int end = 0;
	int sum = 0;
	cout << "Start = ";
	cin >> start;
	cout << "End = ";
	cin >> end;

	for(int i=start; i<=end; i++){
		sum += i;
	}

	cout << "Result = " << sum << endl;
	return 0;
}
