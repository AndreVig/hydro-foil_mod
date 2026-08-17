#include "integrals.h"
#include "particle_data_group.h"
#include "utils.h"
#include "surface.h"
#include <filesystem>
#include <functional>
#include <map>
#include <string>

/*
Run './foil <surface_file> <output_folder> <options>'

TO TURN PARALLEL INTEGRATION OFF(ON) COMMENT(UN-COMMENT) THE OPEN_MP_FLAG LINE OF THE MAKEFILE

Default behaviour:
computation of Lambda polarization at midrapidity using the improved formula and the pysr-modified beta.dat vHLLE file.

Options that can be added in the execution command and their explanantion:
	-D --> to include feed-down corrections (deacys)
	-R --> to include rapidity dependence in the polarization calculation
	--iso --> to use the isothermal approximation (old) formula for the polarization calculation
	--impr-std --> to use the improved formula for the polarization calculation, and the original beta.dat vHLLE file
	--impr-hybrid --> to use the improved formula for the polarization calculation, and the pysr-hybrid beta.dat vHLLE file
*/

using namespace std;

//enum class surface_mode {standard, modified, isothermal};	// "modified" refers to the beta.dat file (output of vHLLE) with the inclusion of the unit normal vector
//enum class formula_mode {improved, improved_mod, old};	// "improved_mod" refers to the improved formula adapted to the modified beta.dat file
enum class polarization_mode {isothermal, improved_std, improved_pysr, improved_hybrid};

// Use "hypersurface_reader" as an alias for the function type used to read the hypersurface file
using hypersurface_reader = std::function<void(string, vector<element>&)>;
// Use "polarization_func_rapidity" as an alias for the function type used to compute the polarization with rapidity dependence
using polarization_func_rapidity = std::function<void(double, double, double, pdg_particle, vector<element>&, ofstream&)>;
// Use "polarization_func_midrapidity" as an alias for the function type used to compute the polarization at midrapidity
using polarization_func_midrapidity = std::function<void(double, double, pdg_particle, vector<element>&, ofstream&)>;



int main(int argc, char** argv){

bool decay = false;
bool rapidity = false;	// to compute polarization as a function of rapidity (instead of at midrapidity)
//surface_mode selected_surface_mode = surface_mode::standard;	// initialize the surface_mode variable
//formula_mode selected_formula_mode = formula_mode::improved;	// initialize the formula_mode variable
polarization_mode selected_polarization_mode = polarization_mode::improved_pysr;

/*
// Map for the hypersurface reading functions
map<surface_mode, hypersurface_reader> reader_map = {
	{surface_mode::standard, read_hypersrface},
	{surface_mode::modified, read_hypersrface_with_normal},
	{surface_mode::isothermal, read_hypersrface_iso}
};

// Map for the polarization functions with rapidity dependence
map<formula_mode, polarization_func_rapidity> formula_map_rapidity = {
	{formula_mode::improved, modified_polarization_rapidity_linear},
	{formula_mode::improved_mod, modified_polarization_rapidity_linear_mod},
	{formula_mode::old, polarization_exact_rapidity}
};

// Map for the polarization functions at midrapidity
map<formula_mode, polarization_func_midrapidity> formula_map_midrapidity = {
	{formula_mode::improved, modified_polarization_midrapidity_linear},
	{formula_mode::improved_mod, modified_polarization_midrapidity_linear_mod},
	{formula_mode::old, polarization_midrapidity_linear}
};
*/

// Map for the hypersurface reading functions
map<polarization_mode, hypersurface_reader> reader_map = {
	{polarization_mode::improved_std, read_hypersrface},
	{polarization_mode::improved_pysr, read_hypersrface},
	{polarization_mode::improved_hybrid, read_hypersrface_with_normal},
	{polarization_mode::isothermal, read_hypersrface_iso}
};

// Map for the polarization functions with rapidity dependence
map<polarization_mode, polarization_func_rapidity> formula_map_rapidity = {
	{polarization_mode::improved_std, improved_polarization_rapidity_linear},
	{polarization_mode::improved_pysr, improved_polarization_rapidity_linear},
	{polarization_mode::improved_hybrid, improved_polarization_rapidity_linear_mod},
	{polarization_mode::isothermal, polarization_exact_rapidity}
};

// Map for the polarization functions at midrapidity
map<polarization_mode, polarization_func_midrapidity> formula_map_midrapidity = {
	{polarization_mode::improved_std, improved_polarization_midrapidity_linear},
	{polarization_mode::improved_pysr, improved_polarization_midrapidity_linear},
	{polarization_mode::improved_hybrid, improved_polarization_midrapidity_linear_mod},
	{polarization_mode::isothermal, polarization_midrapidity_linear}
};

// Map for the output filename
map<polarization_mode, string> output_fname_map = {
	{polarization_mode::improved_std, "/primary_improved_std"},
	{polarization_mode::improved_pysr, "/primary_improved_pysr"},
	{polarization_mode::improved_hybrid, "/primary_improved_hybrid"},
	{polarization_mode::isothermal, "/primary_isothermal"}
};

if(argc<3){
    cout << "INVALID SINTAX!" << endl;
	cout << "use './foil <surface_file> <output_folder> <flags>' to compute Lambda polarization at decoupling." << endl;
	exit(1);
}

if(argc>3){
	for (int i=3; i<argc; i++) {
		if (argv[i]=="-D"s) {
			decay = true;
			cout << "Including calculations for the feed-down corrections!" << endl;
		} else if (argv[i]=="-R"s) {
			rapidity = true;
			cout << "Calculating polarization in the rapidity window [-1,1]!" << endl;
		} else if (argv[i]=="--iso"s) {
			cout << "Using the isothermal approximation (old) formula for the polarization calculation!" << endl;
			selected_polarization_mode = polarization_mode::isothermal;
		} else if (argv[i]=="--impr-std"s) {
			cout << "Using the improved formula and the original beta.dat vHLLE file for the polarization calculation!" << endl;
			selected_polarization_mode = polarization_mode::improved_std;
		} else if (argv[i]=="--impr-hybrid"s) {
			cout << "Using the improved formula and the pysr-hybrid beta.dat vHLLE file for the polarization calculation!" << endl;
			selected_polarization_mode = polarization_mode::improved_hybrid;
		} else {
			cout << "Unknown flag ignored: " << argv[i] << endl;
		}
	}
}
if(selected_polarization_mode == polarization_mode::improved_pysr) {
	cout << "Using the improved formula and the pysr-modified beta.dat vHLLE file for the polarization calculation!" << endl;
}

string surface_file = argv[1];
string output_folder = argv[2];
filesystem::create_directories(output_folder);

vector<element> hypersup = {};
reader_map[selected_polarization_mode](surface_file, hypersup);	// read the hypersurface input file

string output_filename = output_fname_map[selected_polarization_mode];
if(rapidity){
	output_filename = output_filename + "_rapidity";
} else {
	output_filename = output_filename + "_midrapidity";
}

int size_pt = 30;
int size_phi = 30;
int size_y = 20;
vector<double> pT = linspace(0,6,size_pt);
vector<double> phi =  linspace(0,2*PI,size_phi);
vector<double> y_rap =  linspace(-1,1,size_y);

string name_file_primary = output_folder + output_filename;
std::filesystem::path f{name_file_primary};
bool primary_exists = std::filesystem::exists(f);
if(!primary_exists){
	ofstream fout(name_file_primary);
	if (!fout) {
		cout << "I/O error with " << name_file_primary << endl;
		exit(1);
	}

	pdg_particle Lambda(3122);
	Lambda.print();

	if (rapidity) {
		polarization_func_rapidity polarization_calc = formula_map_rapidity[selected_polarization_mode];
		for(double ipt : pT){
			for(double iphi : phi){
				for(double iy : y_rap){
					polarization_calc(ipt, iphi, iy, Lambda, hypersup, fout);
				}
			}
		}
	} else {
		polarization_func_midrapidity polarization_calc = formula_map_midrapidity[selected_polarization_mode];
		for(double ipt : pT){
			for(double iphi : phi){
				polarization_calc(ipt, iphi, Lambda, hypersup, fout);
			}
		}
	}
} else {
	cout<< "Primary file already exists! Skipping calculation..." <<endl;
}





////////////////////////////////////DECAYS///////////////////////////////	 
if(decay){
	string table_file_sigma0 = output_folder + "/TableSigma0";
	string table_file_sigmastar = output_folder + "/TableSigmastar";
	
	pdg_particle Sigma0(3212);
	pdg_particle SigmaStar(3224); //NB: there are three sigma* decaying to lambdas
	Sigma0.print();
	SigmaStar.print();

	std::filesystem::path fS0{table_file_sigma0};
	bool S0table_exists = std::filesystem::exists(fS0);
	std::filesystem::path fSs{table_file_sigmastar};
	bool Sstable_exists = std::filesystem::exists(fSs);
	if(!S0table_exists || !Sstable_exists){
		ofstream fout_sigma0(table_file_sigma0);
			if (!fout_sigma0) {
				cout << "I/O error with " << table_file_sigma0 << endl;
				exit(1);
			}
		ofstream fout_sigmastar(table_file_sigmastar);
			if (!fout_sigmastar) {
				cout << "I/O error with " << table_file_sigmastar << endl;
				exit(1);
			}

		vector<double> pT_table = linspace(0,6.2,2*size_pt);
		vector<double> phi_table =  linspace(0,2*PI,2*size_phi);
		vector<double> y_rap_table =  linspace(-1,1,2*size_y);
		for(double ipt : pT_table){
			for(double iphi : phi_table){
				for(double iy : y_rap_table){
					polarization_exact_rapidity(ipt, iphi, iy, Sigma0, hypersup, fout_sigma0);
					polarization_exact_rapidity(ipt, iphi, iy, SigmaStar, hypersup, fout_sigmastar);
				}
			}
		}
	}
	else{
		cout<<"Tables for interpolation already exist! Skipping calculation..."<<endl;
	}

	string FD_file_sigma0 = output_folder + "/FeedDown_Sigma0";
	string FD_file_sigmastar = output_folder + "/FeedDown_Sigmastar";
	ofstream FDoutsigma0(FD_file_sigma0);
		if (!FDoutsigma0) {
			cout << "I/O error with " << FD_file_sigma0 << endl;
			exit(1);
		}
	ofstream FDoutsigmastar(FD_file_sigmastar);
		if (!FDoutsigmastar) {
			cout << "I/O error with " << FD_file_sigmastar << endl;
			exit(1);
		}

	interpolator spectrum_interpolatorS0(table_file_sigma0,3);
    array<interpolator,4> S_vorticity_interpolatorS0{ {{table_file_sigma0,4},
                                {table_file_sigma0,5},
                                {table_file_sigma0,6},
                                {table_file_sigma0,7}} };
    array<interpolator,4> S_shear_interpolatorS0{ {{table_file_sigma0, 8},
                            {table_file_sigma0, 9},
                            {table_file_sigma0, 10},
                            {table_file_sigma0, 11}} };

	interpolator spectrum_interpolatorSs(table_file_sigmastar,3);
    array<interpolator,4> S_vorticity_interpolatorSs{ {{table_file_sigmastar,4},
                                {table_file_sigmastar,5},
                                {table_file_sigmastar,6},
                                {table_file_sigmastar,7}} };
    array<interpolator,4> S_shear_interpolatorSs{ {{table_file_sigmastar, 8},
                            {table_file_sigmastar, 9},
                            {table_file_sigmastar, 10},
                            {table_file_sigmastar, 11}} };

	for(double ipt : pT){
		for(double iphi : phi){
			for(double iy : y_rap){
				Lambda_polarization_FeedDown(ipt, iphi, iy, Sigma0, 
						spectrum_interpolatorS0, S_vorticity_interpolatorS0, S_shear_interpolatorS0, FDoutsigma0);
				Lambda_polarization_FeedDown(ipt, iphi, iy, SigmaStar, 
						spectrum_interpolatorSs, S_vorticity_interpolatorSs, S_shear_interpolatorSs, FDoutsigmastar);
			}
		}
	}
}

cout<<"The calculation is done!"<<endl;
return 0;
}
