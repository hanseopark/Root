#include <iostream>
using namespace std;

class car{
	private:
		int weight = 1000; //kg
		float fuel = 50; //L
		float speed = 0; //km/h
	
	public:
		car(){};
		~car(){};
		void Accelerate(float level);
		void ShowInfo();
};

void car::Accelerate(float level){
	if(fuel<=0){
		cout << "Empty!" << endl;
	}
	if(fuel>0 && level>= 5*fuel){ //In this case accelerate with all fuel.
		speed += 5*fuel;
		fuel = 0;
	}
	if(fuel>0 && level<5*fuel){
		speed += level;
		fuel -= level/5;
	}
	return;
}

void car::ShowInfo(){
	cout << "Current speed = " << speed << " km/h" << endl;
	cout << "Curren fule   = " << fuel << " L" << endl;
	
	return;
}

int main(){
	car t8415;
	float aclevel = 0;
	for(int i=0; i<10; i++){
		cout << "Accelerate level? >>";
		cin >> aclevel;
		t8415.Accelerate(aclevel);
		t8415.ShowInfo();
	}
	return 0;
}
