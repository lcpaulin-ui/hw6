
    #include <iostream>
    #include "hash.h"
    using namespace std; 
    int main()
   {
        MyStringHash hashk(true);
        string k("abc");
        size_t hk = hashk(k);
        size_t exp = 9953503400;
        cout << hk << endl;
        return 0; 

    }