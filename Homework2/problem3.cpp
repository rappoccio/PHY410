#include<iostream>
#include<string>
int main(){
   bool falafel_special=false, veggie_special=false, salami_special=false 
   int i=0,j=0;
   std::string order;
   std::string new_input;
   std::cout << "Input items to be ordered" << std ::endl;
   getline(std::cin,order);
   std::cout << order[0] << std::endl;
   std::cout << "\nItems Requested........| Price:" << std::endl;
   std::cout <<   "-----------------------|-------" << std::endl; 
   while (i < order.length()) {
      if (order[i] != ' ') {
         if (order[i] == 'e') {
            if (order[i+1] == '1') {
               for (j=0;order.length();++j) {
                  if (order[j] == 'b'){
                     if ((order[j+1] == '1') || (order[j+1] == '2') 
                         || (order[j+1] == '7') || (order[j+1] == '8')) {
                        veggie_special = true;
                        break;
                     }
                  }
               }
               std::cout << "Veggie Burger..........| $7" << std::endl;
            } else if (order[i+1] == '2') {
               for (j=0;order.length();++j) {
                  if (order[j] == 's'){
                     if (order[j+1] == '2') {
                        falafel_special = true
                        break;
                     }
                  }
               }
               std::cout << "Falafel Wrap...........| $6" << std::endl;
            } else if (order[i+1] == '3'){
               for (j=0;order.length();++j) {
                  if (order[j] == 'b'){
                     if ((order[j+1] == '3') || (order[j+1] == '4')
                         || (order[j+1] == '5') || (order[j+1] == '6')) {
                        salami_special = true
                        break;
                     }
                  }
               }
               std::cout << "Salami Sandwich........| $9" << std::endl;
            } else {
               std::cout << "Could not parse selection " << order.substr(i,i+1)
                         << "\nPlease input your selection again." << std::endl;
               std::cin  >> new_input;
               
               order[i] = new_input[0];
               order[i+1] = new_input[1];
               i = 0;
               std::cout << "\nItems Requested........| Price:" << std::endl;
               std::cout <<   "-----------------------|-------" << std::endl;
               continue;
            }
         } else if (order[i] == 's') {
            if (order[i+1] == '1') {
               std::cout << "French Fries...........| $2" << std::endl;
            } else if (order[i+1] == '2') {
               std::cout << "Hummus w Pita Chips....| $3" << std::endl;
            } else if (order[i+1] == '3'){
               std::cout << "Celery and Carrots.....| $2" << std::endl;
            } else {
               std::cout << "Could not parse selection" << order.substr(i,i+1)
                         << "\nPlease input your selection again." << std::endl;
               std::cin  >> new_input;

               order[i] = new_input[0];
               order[i+1] = new_input[1];
               i = 0;
               std::cout << "\nItems Requested........| Price:" << std::endl;
               std::cout <<   "-----------------------|-------" << std::endl;
               continue;
            }
         } else if (order[i] == 'b') {
            if (order[i+1] == '1') {
               std::cout << "Tap Water..............| $0" << std::endl;
            } else if (order[i+1] == '2') {
               std::cout << "Sparkling Water........| $2" << std::endl;
            } else if (order[i+1] == '3'){
               std::cout << "Domestic Beer..........| $4" << std::endl;
            } else if (order[i+1] == '4') { 
               std::cout << "Imported Beer..........| $6" << std::endl;
            } else if (order[i+1] == '5') {
               std::cout << "Red Wine...............| $7" << std::endl;
            } else if (order[i+1] == '6') {
               std::cout << "White Wine.............| $7" << std::endl;
            } else if (order[i+1] == '7') {
               std::cout << "Coffee.................| $1" << std::endl;
            } else if (order[i+1] == '8') {
               std::cout << "Tea....................| $1" << std::endl;
            } else {
               std::cout << "Could not parse selection" << order.substr(i,i+1)
                         << "\nPlease input your selection again." << std::endl;
               std::cin  >> new_input;

               order[i] = new_input[0];
               order[i+1] = new_input[1];
               i = 0;
               std::cout << "\nItems Requested........| Price:" << std::endl;
               std::cout <<   "-----------------------|-------" << std::endl;
               continue;
            }

         } else {
             std::cout << "crap" << i << " " << order[i] << std::endl;
         }
         i = i + 2;
      } else {
         i = i + 1;
      }
   }
   if (veggie_special) {
      specials = specials + 1
   }
   return 0;
}
