#ifndef HASH_H
#define HASH_H

#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <ctime>
#include <cstdlib>.

typedef std::size_t HASH_INDEX_T;

struct MyStringHash {
    HASH_INDEX_T rValues[5] { 983132572, 1468777056, 552714139, 984953261, 261934300 };
    MyStringHash(bool debug = true)
    {
        srand(time(0)); 
        if(false == debug){
            generateRValues();
        }

        else {
            rValues[0] = 983132572;
            rValues[1] = 1468777056;
            rValues[2] = 552714139;
            rValues[3] = 984953261;
            rValues[4]=261934300; 
        }
    }
    // hash function entry point (i.e. this is h(k))
    HASH_INDEX_T operator()(const std::string& k) const
    {

        int str_size = k.size(); 
        int hash_idx = 4; 

        
        // loop over string, six chars at a time (or however big, taking 6steps down each loop )
        unsigned long long currhash[5] ={0,0,0,0,0}; 
        for (int i = k.size(); i > 0; i-=6) {


            unsigned long long temp = 0; 

            int start_idx = i - 6; // will start the string 6 steps down 
            if (start_idx < 0) { start_idx = 0; }

            for (int j = start_idx; j < i; j++){
                temp*= 36;
                temp += letterDigitToNumber(k[j]); 
        
            }

            currhash[hash_idx] = temp; 
            hash_idx--; 
    }

        // hash the string
        unsigned long long res = 0;  
        for (int i = 0; i < 5; i++){
            res+=(rValues[i]*currhash[i] ); 
        }

        return res; 

}
    // A likely helper function is to convert a-z,0-9 to an integral value 0-35
    HASH_INDEX_T letterDigitToNumber(char c) const
    {
        if (c >= '0' && c <= '9'){
            return (c - '0') + 26; 
        }

        else if (c >= 'A' && c <= 'Z'){
            c = c + 'a' - 'A'; 
        } 

        return c - 97; 


    }

    // Code to generate the random R values
    void generateRValues()
    {
        // obtain a seed from the system clock:
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937 generator (seed);  // mt19937 is a standard random number generator

        // Simply call generator() [it has an operator()] to get another random number
        for(int i{ 0 }; i < 5; ++i)
        {
            rValues[i] = generator();
        }
    }
};

#endif
