#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
int main() {
	ofstream eta_file ("eta_values.txt");
	double eta_gap_start, eta_gap_end, eta_gap_step, eta_g;
	double eta_start_diff, eta_end_diff, eta_diff_step, eta_d; 
	
	cout << "Enter start eta_gap value:" << endl;
	cin >> eta_gap_start;
	cout << "Enter end eta_gap value:" << endl;
	cin >> eta_gap_end;
	cout << "Enter eta_gap stepsize:" << endl;
	cin >> eta_gap_step;
	int eta_gap = (int)((eta_gap_end - eta_gap_start)/eta_gap_step);
        cout << eta_gap << endl;

	cout << "Enter start eta_diff value:" << endl;
	cin >> eta_start_diff;
	cout << "Enter end eta_diff value:" << endl;
	cin >> eta_end_diff;
	cout << "Enter end eta_diff value:" << endl;
	cin >> eta_diff_step;
	int eta_diff = (int)((eta_end_diff - eta_start_diff)/eta_diff_step);
	cout << eta_diff << endl;	
        
	if (eta_file.is_open()) {
		eta_g = eta_gap_start;
		for (int i = 0; i <= eta_gap + 1; i++) {
			eta_d = eta_start_diff;
			for (int j = 0; j <= eta_diff; j++) {
				eta_file << setprecision(1) << fixed << eta_g <<" "<< eta_d<<"\n";
				cout << setprecision(1) << fixed << eta_g <<" "<< eta_d<<"\n";
				eta_d += eta_diff_step;
			}
			eta_g += eta_gap_step;		
		}
	} else cout << "Unable to Open file";
	eta_file.close();
	return 0;
}

