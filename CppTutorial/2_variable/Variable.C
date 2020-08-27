#include <iostream>
using namespace std;

int main(){
	//変数を用意して任意の数値で初期化
	float height = 170;
	float weight = 65;
	int birth[2] = {0}; //配列のすべての要素を0で初期化。[2]は2個の要素を持つ配列という意味

	//配列の要素は[0]から始まる
	birth[0] = 1994;
	birth[1] = 6;

	cout << "私の身長は" << height << "fmです。" << endl;
	cout << "私の体重は" << weight << "tです。" << endl;
	cout << "私が生まれたのは、西暦" << birth[0] << "年、つまり令和" << birth[1] << "年です。" << endl;

	return 0;
}
