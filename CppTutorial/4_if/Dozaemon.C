#include <iostream>
#include <string>
using namespace std;

string dozaemon(int say){
	string comment;

	if(say==0) comment = "ぼくタヌキ。";
	else comment = "3分間待ってやる。";

	return comment;
}

int main(){
	cout << dozaemon(0) << endl;
	cout << dozaemon(-1) << endl;
	return 0;
}
