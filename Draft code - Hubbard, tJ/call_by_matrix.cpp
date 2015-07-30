#include<iostream>
using namespace std;

void allocate(int mat[2][3]){
	mat[0][2]=1;
}

void allocate2(const int mat[2][3]){
//	mat[0][2]=3;
}

int main(){
	cout<<(-1+5)%5<<endl;
	int mat[2][3];
	cout<<"mat[0][2]="<<mat[0][2]<<endl;
	allocate(mat);
	cout<<"mat[0][2]="<<mat[0][2]<<endl;
	allocate2(mat);
	cout<<"mat[0][2]="<<mat[0][2]<<endl;
	
}
