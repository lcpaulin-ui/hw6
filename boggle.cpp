#ifndef RECCHECK
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <random>
#include <iomanip>
#include <fstream>
#include <exception>
#endif

#include "boggle.h"

std::vector<std::vector<char> > genBoard(unsigned int n, int seed)
{
	//random number generator
	std::mt19937 r(seed);

	//scrabble letter frequencies
	//A-9, B-2, C-2, D-4, E-12, F-2, G-3, H-2, I-9, J-1, K-1, L-4, M-2, 
	//N-6, O-8, P-2, Q-1, R-6, S-4, T-6, U-4, V-2, W-2, X-1, Y-2, Z-1
	int freq[26] = {9,2,2,4,12,2,3,2,9,1,1,4,2,6,8,2,1,6,4,6,4,2,2,1,2,1};
	std::vector<char> letters;
	for(char c='A'; c<='Z';c++)
	{
		for(int i=0;i<freq[c-'A'];i++)
		{
			letters.push_back(c);
		}
	}
	std::vector<std::vector<char> > board(n);
	for(unsigned int i=0;i<n;i++)
	{
		board[i].resize(n);
		for(unsigned  int j=0;j<n;j++)
		{
			board[i][j] = letters[(r() % letters.size())];
		}
	}
	return board;
}

void printBoard(const std::vector<std::vector<char> >& board)
{
	unsigned int n = board.size();
	for(unsigned int i=0;i<n;i++)
	{
		for(unsigned int j=0;j<n;j++)
		{
			std::cout << std::setw(2) << board[i][j];
		}
		std::cout << std::endl;
	}
}

std::pair<std::set<std::string>, std::set<std::string> > parseDict(std::string fname)
{
	std::ifstream dictfs(fname.c_str());
	if(dictfs.fail())
	{
		throw std::invalid_argument("unable to open dictionary file");
	} 
	std::set<std::string> dict;
	std::set<std::string> prefix;
	std::string word;
	while(dictfs >> word)
	{
		dict.insert(word);
		for(unsigned int i=word.size()-1;i>=1;i--)
		{
			prefix.insert(word.substr(0,i));
		}
	}
	prefix.insert("");
	return make_pair(dict, prefix);
}

std::set<std::string> boggle(const std::set<std::string>& dict, const std::set<std::string>& prefix, const std::vector<std::vector<char> >& board)
{
	std::set<std::string> result;
	for(unsigned int i=0;i<board.size();i++)
	{
		for(unsigned int j=0;j<board.size();j++)
		{
			boggleHelper(dict, prefix, board, "", result, i, j, 0, 1);
			boggleHelper(dict, prefix, board, "", result, i, j, 1, 0);
			boggleHelper(dict, prefix, board, "", result, i, j, 1, 1);
		}
	}
	
	return result;
}


// r,c are curr loc in the board 
// dr , dc search direction
// / dont change dr and dc bnut use them to change r and c 
// build up a word along the direction of search and insert only the longest word found into the std::set<string> result set. 
bool boggleHelper(const std::set<std::string>& dict, const std::set<std::string>& prefix, const std::vector<std::vector<char> >& board, 
								   std::string word, std::set<std::string>& result, unsigned int r, unsigned int c, int dr, int dc)
{
	//add your solution here!

	// is word in dict? 
	
	// base case : finished going down this direction 
	if (r >=board.size() || c >= board[0].size()){
		// can i add it to the result? 
		if (word.size() != 0 && (dict.find(word) != dict.end())){ // i have agone as far as i could, so can i add it to the result?? this will be the longest 
			result.insert(word); 
			return true; 
		}
		return false; // becaude i cant keep going anymore , return false means i have to add to result 
	}

	word = word + board[r][c];  // add next letter 
	if (prefix.find(word) == prefix.end()) {
		if (dict.find(word) != dict.end() ){
			result.insert(word); 
		} 

		return false; 
	} // backtrack before doing anything else because this ins not a valid prefix ie no reason to keep going !! 

	// can i continue or should i stop 


	bool add = false; 

	// recurse 
	if (dr == 0 && dc == 1){

		add = boggleHelper(dict, prefix, board, word, result, r, c+1, 0,1); 
		
		// if (!add){ // if i cant continue, only add to res if its a word 
		// 	if (dict.find(word) != dict.end()){
		// 		result.insert(word);
		// 		return true;
		// 	}
		// 	else {
		// 		return false; 
		// 	}

	}
	// 	}
	// 	else if (add){
	// 		return true; 
	// 	}
	// }

	// else if (dr == 1 && dc == 1) {

	// bool add = boggleHelper(dict, prefix, board, word, result, r + 1, c + 1, 1,1);

	if (dr == 1 && dc == 0){

		 add = boggleHelper(dict, prefix, board, word, result, r+1, c, 1,0); 
		
		// if (!add){ // if i cant continue, only add to res if its a word 
		// 	if (dict.find(word) != dict.end()){
		// 		result.insert(word);
		// 		return true;
		// 	}
		// 	else {
		// 		return false; 
		// 	}

		// }
	}


	if (dr == 1 && dc == 1){

		 add = boggleHelper(dict, prefix, board, word, result, r+1, c+1, 1,1); 
		
	// 	if (!add){ // if i cant continue, only add to res if its a word 
	// 		if (dict.find(word) != dict.end()){
	// 			result.insert(word);
	// 			return true;
	// 		}
	// 		else {
	// 			return false; 
	// 		}

	// 	}
	// }
	}

	if (add == false){
		if (dict.find(word) != dict.end()){
			result.insert(word); 
			return true; 
		}
	}

	return add; 
}

// if not in prefix set backtrack 
// if in prefix set keep going 
// if in dictionary set return true so u can add to set , but keep going there might eb alonger one 
/// from guide  In our version of Boggle we also want to search for and return only the longest word that starts at a particular position (return EARS not EAR). 