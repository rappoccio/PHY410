#include<iostream>
#include<string>
int main(){
   std::string order;
   std::cout << "Input items to be ordered" << std ::endl;
   getline(std::cin,order);
   std::cout << order[1,3] << std::endl;
   std::cout << sizeof(order) << std::endl;
   return 0;
}
