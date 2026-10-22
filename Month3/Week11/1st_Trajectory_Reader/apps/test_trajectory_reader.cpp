#include <iostream>
#include <iomanip>
#include "trajectory_reader.h"
#include "utils.h"

using namespace std;

int main()
{
	cout << "=== TRAJECTORY READER TEST ===" << endl;

	string csv_file = "../../data/trajectory_crctd.csv";
	auto epochs = readTrajectoryCsv(csv_file);

	cout << "Loaded " << epochs.size() << " epochs" << endl;

	if (epochs.empty())
	{
		cerr << "No Epochs loaded" << endl;
		return 1;
	}

	// Show first + last epoch

	// Lambda function
	// This creates a lambda function and stores it in a variable called show.
	// Then you can call it like any function: show(epochs.front(), "First epoch:");

	// auto: Let the compiler figure out the exact type.
	// show: The name of the variable holding the lambda.
	// []: The CAPTURE clause — which outside variables can I use? [] means no variable from main()
	// (const TrajectoryEpoch& e, const char* label): The parameters — like any function.
	// {}: function body
	auto show = [](const TrajectoryEpoch& e, const char* label)
		{
			double x, y, z;
			llaToEcef(e.lat_deg, e.lon_deg, e.height_m, x, y, z);

			cout << label << endl;
			cout << fixed << setprecision(3);
			cout << " t_gps = " << e.t_gps << " s" << endl;
			cout << " lat/lon/h = " << e.lat_deg << ", " << e.lon_deg << ", " << e.height_m << " m" << endl;
			cout << " ECEF x,y,z = " << x << ", " << y << ", " << z << endl;
			cout << " sigma0 = " << e.sigma0_m << " m" << endl;
			cout << " PDOP = " << e.pdop << endl;
			cout << " n_sats = " << e.n_sats << endl;

		};

	show(epochs.front(), "First Epoch:");
	show(epochs.back(), "\nLast Epoch:");

	return 0;
}

