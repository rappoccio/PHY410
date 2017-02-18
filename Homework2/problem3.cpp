#include<iostream>
#include<string>
int main(){
   int i=0;
   std::string order, new_input;
   std::cout << "Input items to be ordered" << std ::endl;
   getline(std::cin,order);
   std::cout << order[0] << std::endl;
   std::cout << "\nItems Requested.....| Price:" << std::endl;
   std::cout <<   "------------------------" << std::endl;
   while (i < order.length()) {
      if (order[i] != ' ') {
         if (order[i] == 'e') {
            if (order[i+1] == '1') {
               std::cout << "Veggie Burger......| $7" << std::endl;
            } else if (order[i+1] == '2') {
               std::cout << "Falafel Wrap.......| $6" << std::endl;
            } else if (order[i+1] == '3'){
               std::cout << "Salami Sandwich....| $9" << std::endl;
            } else {
               std::cout << "Could not understand " << order.substr(i,i+1)
                         << "\nPlease input your selection again." << std::endl;
               std::cin  << new_input << std::endl;
         } else {
             std::cout << "crap" << i << " " << order[i] << std::endl;
         }
         i = i + 3;
      } else {
         i = i + 1;
      }
   }
   return 0;
}
