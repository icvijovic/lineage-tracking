// ##############################################################################
// Performs error correction on parsed barcode files. Not for distribution
// Written and designed by Alex N Nguyen Ba and Jose Rojas Echenique
// nnguyenba@fas.harvard.edu and rojasechenique@fas.harvard.edu
// 2016
//
// Version 1.1. Modified from the original python script by Steve Hanov
//
// To run, use clean script on the Fastq file.
// Then run the error aware python parser to spit out a file that contains barcodes to be corrected.
// One column must indicate the population that the barcode is present in.
// Another column must be barcodes of length > 2
//
// LICENCE 
//
// ##############################################################################

// Compilation: 
// g++ error_correct.cpp -o error_correct -O3 -std=c++0x -lpthread -fopenmp -D_GLIBCXX_PARALLEL

// Typical use: error_correct -pop m -bc n -cpu o file
// Can also accept cin: cat file | error_correct -pop m -bc n -cpu o

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <functional>
#include <thread>
#include <parallel/algorithm>

#include <time.h>
#include <chrono> // Timing functions

#include <unistd.h> // read
#include <cstring>  // memchr
#include <mutex>

#include <sys/stat.h> // Check file exists

#include <signal.h> // Deal with interrupted events to remove temporary file.
 
using namespace std;

// The minimum cost of a given word to be changed to a word of the dictionary
int min_cost;
std::mutex mtx; 
int pop_id = 0;
int bc_id = 0;
int err_max = 2;
char tempFileName[] = "temporary_confirmed.XXXXXX";

static double diffclock(clock_t clock1,clock_t clock2)
{
    double diffticks=clock1-clock2;
    double diffms=(diffticks)/(CLOCKS_PER_SEC/1000);
    return diffms;
}
void signal_callback_handler(int signum){
	// Cleanup and close up stuff here
	unlink(tempFileName);
	// Terminate program
	exit(signum);
}

// Trie's node
struct trie{
    typedef std::unordered_map<char, trie*> next_t;
 
    // The set with all the letters which this node is prefix
    next_t next;
    // If word is equal to "" is because there is no word in the
    //  dictionary which ends here.
    string word;
 
    trie() : next(std::unordered_map<char, trie*>()) {}

    void insert(string w){
        w = "$" + w;
         
        int sz = w.size();
         
        trie* n = this;
        for (int i = 0; i < sz; ++i) {
            if (n->next.find(w[i]) == n->next.end()) {
                n->next[w[i]] = new trie();
            }
 
            n = n->next[w[i]];
        }
 
        n->word = w;
    }

};

void search_impl_ukk2(trie* tree, char ch, vector<int> &last_row, const string& word, int max_cost, std::unordered_map<string,int> &return_array, int half_depth, int half_cost, int left, int right, int depth, int reversed){
	int sz = last_row.size();
 
 	int min_element = 0x3f3f3f3f;
	// Calculate the min cost of insertion, deletion, match or substution
	int insert_or_del, replace;

	vector<int> current_row(sz,max_cost+1);
	current_row[0] = last_row[0] + 1;

	int min_right = -1;
	int max_left = sz;

	for (int i = left; i < right; ++i) {
		insert_or_del = min(current_row[i-1], last_row[i])+1;
		replace = (word[i-1] == ch) ? last_row[i-1] : (last_row[i-1] + 1);
 
		current_row[i] = min(insert_or_del, replace);
		min_element = min(current_row[i],min_element);

		// Keep track of where left and right becomes max cost
		if(current_row[i] <= max_cost){
			max_left = min(max_left,i);
			min_right = max(min_right,i);
		}
	}

	min_right = min(sz,min_right+2);
	if(max_left != 1){
		max_left += 1;
	}

	//int max_right = min(sz,right+1);
	//int min_left = max(1,right-max_cost*2-1);

	// When we find a cost that is less than the min_cost, is because
	// it is the minimum until the current row, so we update
	if ((current_row[sz-1] <= max_cost) && (tree->word != "")) {
        	//min_cost = current_row[sz-1]; use this to only return values that have the minimum cost.
		// If we found something lower than the maximum cost, then we populate the return array.
		string return_value = tree->word;
		if(reversed == 1){
			reverse(return_value.begin()+1,return_value.end());
		}

		return_array.insert({return_value,current_row[sz-1]});

  	}
 
	// If there is an element which is smaller than the current minimum cost,
	//  we can have another cost smaller than the current minimum cost, and we should continue
	if(depth <= half_depth && min_element <= half_cost){
		int next_depth = depth+1;
		for (trie::next_t::iterator it = tree->next.begin(); it != tree->next.end(); ++it) {
			search_impl_ukk2(it->second, it->first, current_row, word, max_cost, return_array, half_depth, half_cost, max_left, min_right, next_depth,reversed);
		}
	}
	else if (min_element <= max_cost && depth>half_depth) {
		int next_depth = depth+1;
		for (trie::next_t::iterator it = tree->next.begin(); it != tree->next.end(); ++it) {
			search_impl_ukk2(it->second, it->first, current_row, word, max_cost, return_array, half_depth, half_cost, max_left, min_right, next_depth,reversed);
		}
	}
}
std::unordered_map<string,int> search_ukk2(string word,trie &tree, trie &reverse_tree, int max_cost = 0){

	unsigned int half_depth_left = word.size();
	half_depth_left /= 2;
	int half_depth_right = word.size() - half_depth_left;

	string reversed_word = word;
	reverse(reversed_word.begin(),reversed_word.end());

	word = "$" + word;
	reversed_word = "$" + reversed_word;
     
	int sz = word.size();	

	vector<int> current_row(sz + 1);
	// Naive DP initialization
	for (int i = 0; i <= sz; ++i) current_row[i] = i;
	
	
	std::unordered_map<string,int> return_array;     

	// Search the forward trie
	// Search is done by halves.
	// For example, given a maximum cost of 2, a word has several ways of having mismatches:
	// 0 mismatches in the first half, 0 mismatches in the second half.
	// 1, 0
	// 2, 0
	// 1, 1
	// 0, 2
	// 0, 1
	// Therefore, we can completely find all the words with 2 mismatches by first looking for anyword that has 0 mismatches in the 1rst half, allowing 2 mismatches total.
	// We then search from the back, allowing 1 mismatch in the 2nd half, and allowing 2 total.
	// For 3, it can be 1 in first half, 2 in second half.


	for (trie::next_t::iterator it = tree.next.begin(); it != tree.next.end(); ++it) {
		// Here, we search by the node.
        	search_impl_ukk2(it->second, it->first, current_row, word, max_cost, return_array, half_depth_left, max(0,(max_cost+1)/2-1), 1,1+2+max_cost, 1,0);
    	}

	// Search the reverse trie
	for (trie::next_t::iterator it = reverse_tree.next.begin();it != reverse_tree.next.end();++it){
		search_impl_ukk2(it->second, it->first, current_row, reversed_word, max_cost, return_array, half_depth_right, max_cost/2, 1,1+2+max_cost, 1,1);

	}
	

 	// Return the matches, and the cost.
	return return_array;

	

}
void delete_from_trie(string word, trie &tree){
	
	// Strategy is to create a vector with pointers to the internal trie nodes until the final match
	word = "$" + word;
	
	int sz = word.size();

	trie * n;
	n = &tree;

	vector<trie*> path;
	path.reserve(sz);
	for(int i = 0;i<sz;++i){		
		if (n->next.find(word[i]) != n->next.end()) {		
			n = n->next[word[i]];
			path.push_back(n);
		}
		else{
			// Word not present in trie
			return;
		}
	}

	// Now see if the path was properly saved
	// GO from the last path up until and delete as you go until a termination condition arises.

	// There are a few simple termination conditions that we can check straight away
	if(path[path.size()-1]->next.size() != 0){
		path[path.size()-1]->word = ""; // Delete the word, done!
		return;
	}

	for(int i = path.size()-2;i > 0;--i){ // > 0 because we never delete the root
		// We've already checked that the next node is deletable. Now we just check when to stop.
		path[i]->next.erase(word[i+1]);

		if(path[i]->next.size() == 0 && path[i]->word == ""){
			//keep going
		}
		else{
			return;
		}	
	}
}

int string_to_int(const char *p) {
    int x = 0;
    bool neg = false;
    if (*p == '-') {
        neg = true;
        ++p;
    }
    while (*p >= '0' && *p <= '9') {
        x = (x*10) + (*p - '0');
        ++p;
    }
    if (neg) {
        x = -x;
    }
    return x;
}

class StringRef
{
private:
    char const*     begin_;
    int             size_;

public:
    int size() const { return size_; }
    char const* begin() const { return begin_; }
    char const* end() const { return begin_ + size_; }

    StringRef( char const* const begin, int const size )
        : begin_( begin )
        , size_( size )
    {}
};

vector<StringRef> split4(const char *str, int size, char delimiter = '	'){
	vector<StringRef> results;
	const char *p = str;
	const char *q = str;
	
	while(p = (char*) memchr(p,delimiter,(str+size) - p)){
		results.emplace_back(StringRef(q,p-q));

		++p;
		q = p;
	}
	results.emplace_back(StringRef(q,(str+size)-q));

	return results;
}
int correct(int pop,int start_c, int stop_c, std::unordered_map<int,vector<string>> &bcs, std::unordered_map<int,trie> &final_confirmed, std::unordered_map<int,trie> &reversed_final_confirmed, std::unordered_map<int,std::unordered_map<string,vector<string>>> &reads,int &corrected, int &ndiscarded){
	std::unordered_map<string,int> matches;
	string bc;

	for(int it2 = start_c; it2 < stop_c;++it2){
		bc = bcs[pop][it2];
		
		matches = search_ukk2(bc,final_confirmed[pop],reversed_final_confirmed[pop],err_max);

		if(matches.size() > 0){
			// Found at least one match. Correct.
			// Correct using the best match	

			auto it3_begin = matches.begin();
			auto it3_end = matches.end();
			auto best_match = it3_begin;
			int prev = best_match->second;		

			while(++it3_begin != it3_end){
				if(prev > it3_begin->second){
					best_match = it3_begin;
					prev =  best_match->second;
				}
			}

			string m = best_match->first.substr(1);

			int size = reads[pop][bc].size();
			mtx.lock();
			corrected += size;
			mtx.unlock();
			
			for(int i = 0;i < size;++i){
				vector<StringRef> const fields = split4(reads[pop][bc][i].c_str(),reads[pop][bc][i].length(),'	');

				string output;
				output.reserve(reads[pop][bc][i].length()+16);

				output.append(fields[0].begin(),fields[0].end() - fields[0].begin());
				for(unsigned int j = 1;j < fields.size();++j){
					if(j == bc_id){
						output.append("	").append(m);
					}
					else{
						output.append("	").append(fields[j].begin(),fields[j].end() - fields[j].begin());

					}
				}
				output.append("\n");
				
				cout << output;
			}
		}
		else{
			mtx.lock();
			ndiscarded += reads[pop][bc].size();
			mtx.unlock();
		}

	}
}
int correct_OMP(int pop,int start_c, int stop_c, std::unordered_map<int,vector<string>> &bcs, std::unordered_map<int,trie> &final_confirmed, std::unordered_map<int,trie> &reversed_final_confirmed, std::unordered_map<int,std::unordered_map<string,vector<string>>> &reads,int &corrected, int &ndiscarded){
	#pragma omp parallel for reduction(+:corrected, ndiscarded)
	for(int it2 = start_c; it2 < stop_c;++it2){
		string bc = bcs[pop][it2];
		
		std::unordered_map<string,int> matches = search_ukk2(bc,final_confirmed[pop],reversed_final_confirmed[pop],err_max);

		if(matches.size() > 0){
			// Found at least one match. Correct.
			// Correct using the best match

			auto it3_begin = matches.begin();
			auto it3_end = matches.end();
			auto best_match = it3_begin;
			int prev = best_match->second;		

			while(++it3_begin != it3_end){
				if(prev > it3_begin->second){
					best_match = it3_begin;
					prev =  best_match->second;
				}
			}

			string m = best_match->first.substr(1);

			int size = reads[pop][bc].size();
			corrected += size;
		
			for(int i = 0;i < size;++i){
				vector<StringRef> const fields = split4(reads[pop][bc][i].c_str(),reads[pop][bc][i].length(),'	');

				string output;
				output.reserve(reads[pop][bc][i].length()+16);

				output.append(fields[0].begin(),fields[0].end() - fields[0].begin());
				for(unsigned int j = 1;j < fields.size();++j){
					if(j == bc_id){
						output.append("	").append(m);
					}
					else{
						output.append("	").append(fields[j].begin(),fields[j].end() - fields[j].begin());

					}
				}
				output.append("\n");
				
				cout << output;
			}
		}
		else{
			ndiscarded += reads[pop][bc].size();
		}
	}
}

int correct_first_pass(int pop,std::unordered_map<int,std::unordered_map<string,int>> &confirmed_reads,std::unordered_map<int,trie> &confirmed, std::unordered_map<int,trie> &reversed_confirmed, std::unordered_map<int,unordered_map<string,string>> &correction_map,int &corrected_in_first_pass){
	// Make a vector of pairs for this population, we will sort this
		std::vector<pair<int, string>> pairs;
		pairs.reserve(confirmed_reads[pop].size());
		for(auto it2 = confirmed_reads[pop].begin();it2 != confirmed_reads[pop].end();++it2){
			// Make a vector of pairs based on these values
			pairs.emplace_back(pair<int,string>(confirmed_reads[pop][it2->first],it2->first));
		}

		// Now sort the vector of pairs
		sort(pairs.begin(), pairs.end(), greater<pair<int,string> >());

		// Make a new trie.
		// Now loop through the vector and create a new trie based on error correction
		for(auto it2 = pairs.begin();it2 != pairs.end();++it2){
			// The barcode is it2->second
			// We check if that barcode is in the confirmed_reads array, else, it has been error corrected.
			
			if(confirmed_reads[pop][it2->second] == 0){
				continue;
			}
			else{
				// This is a real barcode, so we add it to the new trie. 
				//final_confirmed[pop].insert(it2->second);

				// Add it to the reversed trie
				//string reversed_bc = it2->second;
				//reverse(reversed_bc.begin(),reversed_bc.end());
				//reversed_final_confirmed[pop].insert(reversed_bc);

				// Add to the correction map
				correction_map[pop].insert({it2->second,it2->second});
				
				// We search the old trie for all barcodes within 2 errors, and make sure they do not appear in the new trie.
				std::unordered_map<string,int> matches;
				matches = search_ukk2(it2->second,confirmed[pop],reversed_confirmed[pop], err_max);

				// If there's a match, then we should make a map so that the next search can be done fast.
				// We delete the match from the trie, so that it can never be used as a correct barcode in the following steps.
				for(auto it3 = matches.begin();it3 != matches.end();++it3){
					
					string matched_bc = it3->first.substr(1);
					
					// Check if the matched barcode is less than 5% of the current count. Errors should not be a significant portion of the cluster.
					if(confirmed_reads[pop][matched_bc] > confirmed_reads[pop][it2->second]/32 || matched_bc == it2->second){
						continue;
					} 

					
					string reversed_matched_bc = matched_bc;
					reverse(reversed_matched_bc.begin(),reversed_matched_bc.end());
					// Remove from trie
					delete_from_trie(matched_bc,confirmed[pop]);
					delete_from_trie(reversed_matched_bc,reversed_confirmed[pop]);

					correction_map[pop].insert({matched_bc,it2->second});
					mtx.lock();
					corrected_in_first_pass += confirmed_reads[pop][matched_bc];
					mtx.unlock();
					confirmed_reads[pop].erase(matched_bc);
				}
							
			}
			
		}
		// Done populating the new trie for this population
}
int file_exists(char *name){
  struct stat   buffer;
  return (stat (name, &buffer) == 0);
}
int main(int argc, char * argv[]) {
	cin.tie(NULL);
	//omp_set_num_threads(6);
	istream *in;
	// Must be declared here for scope reasons
	ifstream ifn;

	// Parse arguments
	// Typical use: error_correct -pop n -bc m file
	// Can also accept cin: cat file | error_correct -pop n -bc m


	int is_file = 0;
	int cpu_count = 1;
	int no_pop = 0;

	if(argc < 5){
		cerr << "error_correct [-cpu i] [-err i] -pop n -bc m file" << endl;
		cerr << "cat file | error_correct [-cpu i] [-err i] -pop n -bc m" << endl;
		return 1;
	}
	else{
		for(int i = 1;i < argc;++i){
			if(string(argv[i]) == "-pop"){
				pop_id = string_to_int(argv[i+1]);
				++i;
				continue;
			}
			if(string(argv[i]) == "-bc"){
				bc_id = string_to_int(argv[i+1]);
				++i;
				continue;
			}
			if(string(argv[i]) == "-err"){
				err_max = string_to_int(argv[i+1]);
				++i;
				continue;
			}
			if(string(argv[i]) == "-cpu"){
				cpu_count = string_to_int(argv[i+1]);
				omp_set_num_threads(cpu_count);
				++i;
				continue;
			}
			if(i == argc-1){
				is_file = 1;
			}
		}
	}
	if(pop_id == bc_id){
		cerr << "Column for the population and the barcodes should be different." << endl;
		return 1;
	}
	if(pop_id == 0){
		no_pop = 1;
	}
	if(is_file == 1){
		if(file_exists(argv[argc-1])){
			ifn.open(argv[argc-1]);
			in = &ifn;
		}
		else{
			cerr << "Cannot access file " << argv[argc-1] << endl;
			return 1;
		}
	}
	else{
		in = &cin;
	}
	--pop_id;
	--bc_id;


	// Read the file
	std::unordered_map<int,trie> confirmed;
	std::unordered_map<int,trie> reversed_confirmed;

	std::unordered_map<int,std::unordered_map<string,vector<string>>> reads;
	std::unordered_map<int,std::unordered_map<string,int>> confirmed_reads;

	// Make a temporary file with all the confirmed barcodes.
	int fd = mkstemp(tempFileName);
	if (fd == -1){
		cerr << "Could not make a temporary file." << endl;
		return 1;
	}
	
	ofstream temporary_confirmed;
	temporary_confirmed.open(tempFileName);
	signal(SIGINT, signal_callback_handler);
	int number_of_reads_total = 0;
	int number_of_reads_temporary_confirmed = 0;

	int time_0, time_1;
	time_0 = clock();
	for(string read; getline( *in, read); ){
		vector<StringRef> const fields = split4(read.c_str(),read.length());
	
		++number_of_reads_total;
		int pop;
		if(no_pop == 1){
			pop = 0;
		}
		else{
			pop = string_to_int(string(fields[pop_id].begin(),fields[pop_id].end()).c_str());
		}
		string bc(fields[bc_id].begin(),fields[bc_id].end());
		
		// Search the trie for the bc
		// Has the tree been created for this population?
		if(confirmed_reads[pop].find(bc) != confirmed_reads[pop].end()){
			++confirmed_reads[pop][bc];
			++number_of_reads_temporary_confirmed;

			// Print the read to the temporary_confirmed read file
			temporary_confirmed << read << "\n";

			// cout << read << endl;
			continue;
		}		
		// Once you reach 10, then move the reads from reads to confirmed_reads. Add the confirmed to the trie.
		if(reads[pop][bc].size() == 9){
			confirmed[pop].insert(bc);

			// Insert into the reverse trie the reversed string
			string reversed_bc = bc;
			reverse(reversed_bc.begin(),reversed_bc.end());
			reversed_confirmed[pop].insert(reversed_bc);

			confirmed_reads[pop][bc] = 10;
			number_of_reads_temporary_confirmed += 10;

			temporary_confirmed << read << "\n";

			for(int i = 0;i<9;++i){
				temporary_confirmed << reads[pop][bc][i] << "\n";
			}
			reads[pop].erase(bc);
			continue;
		}

		// None of the above passed the test, add the read to the vector.
		reads[pop][bc].push_back(read);

		
		
	}
	time_1 = clock();
	cerr << "Done reading reads. There were " << number_of_reads_total << " reads that had more than 3 columns." << endl;
	cerr << "Moved " << number_of_reads_temporary_confirmed << " to the temporary confirmed file." << endl;
	cerr << "This took : " << diffclock(time_1,time_0) << " milliseconds." << endl << endl;
	
	// close the temporary confirmed file
	temporary_confirmed.close();

	// At this point, all the reads are within the reads vector and into the temporary file.
	// We now want a list of reads within the temporary file sorted by their highest count. Can we get this?
	// Loop through the confirmed reads.

	std::unordered_map<int,unordered_map<string,string>> correction_map;
	int corrected_in_first_pass = 0;

	vector<thread> v;
	time_0 = clock();

	for(auto it=confirmed_reads.begin(); it != confirmed_reads.end();++it){
		int pop = it->first;

		v.push_back(thread(correct_first_pass,pop,ref(confirmed_reads),ref(confirmed),ref(reversed_confirmed),ref(correction_map),ref(corrected_in_first_pass)));
	}
	for(int i = 0;i < v.size();++i){
		v[i].join();
	}
	time_1 = clock();
	cerr << "Done creating clusters. This took " << diffclock(time_1,time_0) << " milliseconds." << endl << endl;

	time_0 = clock();
	// New trie can now be used to error correct. First, we need to correct all the initially thought to be confirmed barcodes. 
	ifstream ifn2;
	ifn2.open(tempFileName);
	in = &ifn2;
	//chars = 0;
	
	
	string prefix;
	for(string read; getline( *in, read); ){
		vector<StringRef> const fields = split4(read.c_str(),read.length());
		int pop;
		if(no_pop == 1){
			pop = 0;
		}
		else{
			pop = string_to_int(string(fields[pop_id].begin(),fields[pop_id].end()).c_str());
		}

		string bc = string(fields[bc_id].begin(),fields[bc_id].end());

		// If the correction is equal, just print the read, because it's faster maybe.

		if(bc == correction_map[pop][bc]){
			cout << read.append("\n");
		}
		else{
			// print the corrected read
			int size = fields.size();
			string output;
			output.reserve(read.length()+16);

			output.append(fields[0].begin(),fields[0].end() - fields[0].begin());
			for(unsigned int j = 1;j < size;++j){
				if(j == bc_id){
					output.append("	").append(correction_map[pop][bc]);
				}
				else{
					output.append("	").append(fields[j].begin(),fields[j].end() - fields[j].begin());
				}
			}
			output.append("\n");
			
			cout << output;
		}
	}
	time_1 = clock();
	cerr << "Corrected " << corrected_in_first_pass << " in the first pass." << endl;
	cerr << "This took : " << diffclock(time_1,time_0) << " milliseconds." << endl << endl;
	unlink(tempFileName);
	// Now correct the remainers.
	int ndiscarded = 0;
	int corrected = 0;

	// We're going to try to thread the correction within a single population.
	// Problem, the reads[pop] is an unordered map.
	// Unordered maps only have forward iterators, and so it's not so simple to just separate the whole map into chunks of equal size to be threaded.
	// Although in practice it is possible to forward the iterator to the right location, this is O(n).
	// An alternative, memory heavy, solution is to put all the barcode IDs into a vector, which has a random access iterator.
	// We then thread through the vector instead.
	std::unordered_map<int, vector<string>> bcs;
	int k = 0;
	time_0 = clock();
	auto t_start = std::chrono::high_resolution_clock::now();
	for(auto it=confirmed.begin(); it != confirmed.end();++it){
		int pop = it->first;
		int size_of_pop = reads[pop].size();
	
		bcs[pop].reserve(size_of_pop);
		for(auto it2 = reads[pop].begin();it2 != reads[pop].end();++it2){
			bcs[pop].push_back(it2->first);
		}

		correct_OMP(pop,0,size_of_pop,bcs, confirmed, reversed_confirmed, reads, corrected, ndiscarded);
	}
	time_1 = clock();
	auto t_end = std::chrono::high_resolution_clock::now();

	cerr << "Corrected " << corrected << endl;
	cerr << "Discarded " << ndiscarded << endl;
	cerr << "This took : " << diffclock(time_1,time_0) << " milliseconds." << endl;
	cerr << "Wallclock : " << std::chrono::duration<double, std::milli>(t_end-t_start).count() << " milliseconds." << endl;



}