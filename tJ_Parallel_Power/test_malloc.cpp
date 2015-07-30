#include<iostream>
#include<stdlib.h>

using namespace std;

int main(){
int ** array;
int i;
int *intptr;
int Size=4;
array = (int**) malloc(Size * sizeof(int*)+ Size*Size*sizeof(int) );
for(i=0, intptr=(int*)(array+Size); i<Size; i++, intptr+= Size){
	array[i]=intptr;
}

for(int l=0;l<4;l++){
	for(int m=0; m<4; m++){
		array[l][m]=l+m;
		cout<<array[l][m]<<endl;
	}
}




return 0;
}
