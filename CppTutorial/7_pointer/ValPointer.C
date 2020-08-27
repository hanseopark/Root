#include <iostream>
using namespace std;

int main(){
	int add = 10;
	int *p;
	p = &add;

	cout << "変数addの値は " << add << " です。" << endl;
	cout << "変数addのアドレスは " << p << " です。" << endl;
	return 0;
}
