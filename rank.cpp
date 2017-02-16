#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

vector <double> get_rank(vector <double> & data, string method = "average") {
	int len = data.size();
	vector <int> index(len, 0);
	vector <double> rank(len, 0);
	for (int i = 0; i < len; ++i) {
		index[i] = i;
	}
	sort(index.begin(), index.end(),
	    [&](const int & a, const int & b) {
	    	return (data[a] > data[b]);
		}
	);
	if (method == "average") {
		int sumranks = 0;
        int dupcount = 0;
        for (int i = 0; i < len; ++i) {
    	    sumranks = sumranks + i;
    	    dupcount = dupcount + 1;
    	    if ((i == len -1) || (data[index[i]] != data[index[i + 1]]) ) {
    		    double averank = double(sumranks) / double(dupcount) + 1;
    		    for (int j = i - dupcount + 1; j < i + 1; ++j) {
    			    rank[index[j]] = averank;
			    }
			    sumranks = 0;
			    dupcount = 0;
		    }
	    }
	} else if (method == "min") {
		int current = 1, skip = 1;
        rank[index[0]] = current;
        for (int i = 1; i < len; i++) {
            if (data[index[i]] == data[index[i - 1]]) {
                skip++;
            } else {
                current += skip;
                skip = 1;
            }
            rank[index[i]] = current;
        }		
	} else if (method == "dense") {
		int current = 1;
        rank[index[0]] = current;
        for (int i = 1; i < len; i++) {
            if (data[index[i]] != data[index[i - 1]]) {
                current++;
            }
            rank[index[i]] = current;
        }
	}
	return rank;
}

int main() {
    vector<double> data = {3.4,  6.2, 8.8, 3.4, 2.3};
    vector <double> rank_average = get_rank(data, "average");
    vector <double> rank_min = get_rank(data, "min");
    vector <double> rank_dense = get_rank(data, "dense");
    int len = data.size();
    cout << "original" << '\t' << "rank_average" << '\t' << "rank_min" << '\t' << "rank_dense" << endl;
    for (int i = 0; i < len; i++) {
    	cout << data[i] << '\t' << rank_average[i] << '\t' << rank_min[i] << '\t' << rank_dense[i] << endl;
    }
    return 0;
}
