#include<iostream>
#include<string>
int main(){
   int i=0,j=0,n=0,a;
   unsigned int total=0;
   std::string order;
   std::string new_input;
   std::cout << "Input items to be ordered" << std ::endl;
   getline(std::cin,order);
   //std::cout << order.length() << std::endl;
   //a = order.length();
   std::cout << "\nItems Requested........| Price:" << std::endl;
   std::cout <<   "-----------------------|-------" << std::endl;
   for (j=0;j<order.length();++j) {
      
                     //std::cout << "5" << std::endl;
      if (order[j] == 'e') {
         if (order[j+1] == '1') {
                     //std::cout << "6" <<n<< std::endl;
            for (n=0;n<order.length();++n) {
                     //std::cout << n<<" "<<j << std::endl;
               if (order[n] == 'b'){
                  if ((order[n+1] == '1') || (order[n+1] == '2') 
                      || (order[n+1] == '7') || (order[n+1] == '8')) {
                     order[j] = 'c';
                     order[n] = ' ';
                     order[n+1] = ' ';
                     //std::cout << "1" << std::endl;
                     break;
                  }
               }
            }
         continue;
         } else if (order[j+1] == '2') {
            for (n=0;n<order.length()-1;++n) {
               if (order[n] == 's'){
                  if (order[n+1] == '2') {
                     order[j] = 'c';
                     order[n] = ' ';
                     order[n+1] = ' ';
                     //std::cout << "2" << std::endl;
                     break;
                  }
               }
            }
         } else if (order[j+1] == '3'){
            for (n=0;n<order.length()-1;++n) {
               if (order[n] == 'b'){
                  if ((order[n+1] == '3') || (order[n+1] == '4')
                      || (order[n+1] == '5') || (order[n+1] == '6')) {
                     order[j] = 'c';
                     order[n] = ' ';
                     order[n+1] = ' ';
                     //std::cout << "3" << std::endl;
                     break;
                  }
               }
            }
         }
      }
   }
   i=0;
   while (i < order.length()) {
      if (order[i] != ' ') {
         //std::cout << "4" << std::endl;
         if (order[i] == 'e') {
            if (order[i+1] == '1') {
               std::cout << "Veggie Burger..........| $7" << std::endl;
               total = total + 7;
            } else if (order[i+1] == '2') {
               std::cout << "Falafel Wrap...........| $6" << std::endl;
               total = total + 6;
            } else if (order[i+1] == '3'){
               std::cout << "Salami Sandwich........| $9" << std::endl;
               total = total + 9;
            } else {
               std::cout << "Could not parse selection " 
                         << order.substr(i,i+1)
                         << "\nPlease input your selection again." 
                         << std::endl;
               std::cin  >> new_input;
               
               order[i] = new_input[0];
               order[i+1] = new_input[1];
               total = 0;
               i = 0;
               std::cout << "\nItems Requested........| Price:" << std::endl;
               std::cout <<   "-----------------------|-------" << std::endl;
               continue;
            }
         } else if (order[i] == 's') {
            if (order[i+1] == '1') {
               std::cout << "French Fries...........| $2" << std::endl;
               total = total + 2;
            } else if (order[i+1] == '2') {
               std::cout << "Hummus w Pita Chips....| $3" << std::endl;
               total = total + 3;
            } else if (order[i+1] == '3'){
               std::cout << "Celery and Carrots.....| $2" << std::endl;
               total = total + 2;
            } else {
               std::cout << "Could not parse selection" 
                         << order.substr(i,i+1)
                         << "\nPlease input your selection again." 
                         << std::endl;
               std::cin  >> new_input;

               order[i] = new_input[0];
               order[i+1] = new_input[1];
               total = 0;
               i = 0;
               std::cout << "\nItems Requested........| Price:" << std::endl;
               std::cout <<   "-----------------------|-------" << std::endl;
               continue;
            }
         } else if (order[i] == 'b') {
            if (order[i+1] == '1') {
               std::cout << "Tap Water..............| $0" << std::endl;
               total = total + 0;
            } else if (order[i+1] == '2') {
               std::cout << "Sparkling Water........| $2" << std::endl;
               total = total + 2;
            } else if (order[i+1] == '3'){
               std::cout << "Domestic Beer..........| $4" << std::endl;
               total = total + 4;
            } else if (order[i+1] == '4') { 
               std::cout << "Imported Beer..........| $6" << std::endl;
               total = total + 6;
            } else if (order[i+1] == '5') {
               std::cout << "Red Wine...............| $7" << std::endl;
               total = total + 7;
            } else if (order[i+1] == '6') {
               std::cout << "White Wine.............| $7" << std::endl;
               total = total + 7;
            } else if (order[i+1] == '7') {
               std::cout << "Coffee.................| $1" << std::endl;
               total = total + 7;
            } else if (order[i+1] == '8') {
               std::cout << "Tea....................| $1" << std::endl;
               total = total + 8;
            } else {
               std::cout << "Could not parse selection" 
                         << order.substr(i,i+1)
                         << "\nPlease input your selection again." 
                         << std::endl;
               std::cin  >> new_input;

               order[i] = new_input[0];
               order[i+1] = new_input[1];
               total = 0;
               i = 0;
               std::cout << "\nItems Requested........| Price:" << std::endl;
               std::cout <<   "-----------------------|-------" << std::endl;
               continue;
            }
         } else if (order[i] == 'c') {
            if (order[i+1] == '1') {
               std::cout << "Veggie Burger Special..| $8" << std::endl;
               total = total + 8;
            } else if (order[i+1] == '2') {
               std::cout << "Falafel Wrap Special...| $7" << std::endl;
               total = total + 7;
            } else if (order[i+1] == '3') {
               std::cout << "Salami Sandwich Special| $13" << std::endl;
               total = total + 13;
            }
         } else {
             i = i + 1;
         }
         i = i + 2;
      } else {
         i = i + 1;
      }
   }
   std::cout << "TOTAL = $" << total << std::endl;
   return 0;
}
