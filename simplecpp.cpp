/*Simple  C++ program for learning stuff*/

//#include <stdio>
#include <iostream>
//#include <stdlib>
//#include "stdlib.h"
#include <string>
#include <assert.h>
using namespace std;


void myfuncPrint(int num1, int num2);
void thisf();

/* This is where program runs */
int main( int argc,  char** argv )
{

 int number1 = atoi(argv[1]);
 int number2 = atoi(argv[2]);
 
 int numb1 = 0;
 int numb2 = 0; 

/*
for(int i=0; i < argc; i++)
 {  number1 = argv[i];
    number2 = argv[i+1];
     
// }
*/

myfuncPrint(number1, number2);
thisf();

return 0;
}



void myfuncPrint(int n1, int n2)
{
std::cout <<"This is C++ " << std::endl;
std::cout <<" I really Love it\n" << std::endl;

std::cout <<" The Sum of Numbers \t" << n1 << " and\t" << n2  << " is \n" << std::endl;

std::cout <<"Sum = " << n1 +  n2  << std::endl;
} 

void thisf()
{

std::cout <<" This funxtion is really nice can I change it!! " << std::endl;
}
