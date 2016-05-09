#include <iostream>
#include <fstream>
//#define debs
#include "plot.h"

using std::cout;
using std::endl;
using std::string;
using std::ofstream;


int main (int argc, char *argv[]){

	cout << "******************************** "<<PACKAGE_STRING<<" ***********************************" << endl;
	LHAPDF::initLHAPDF();
	LHAPDF::setPDFPath(LHAPDF_PDF_DIR);
	cout << "PDF path = " << LHAPDF::pdfsetsPath() << endl;




	Plot plot;
	gsl_set_error_handler_off();

	ofstream of;
	of.open("debug.out");
	of.close();


	if(argc>=2){
		cout << "argv[" << 0 << "] = " << argv[0] << endl;
		for(int i=1; i<argc; ++i){
		cout << "argv[" << i << "] = " << argv[i] << endl;
			plot.init(argv[i]);
			plot.createPlot();
		}
	}else{
		cout << " ! NO FILE ARGUMENT ! " << endl;
		plot.print_standard_in();
		plot.init("example.in");
		plot.createPlot();
	}


	if(debug.size()){
		of.open("debug.out");
		for(uint i=0; i<debug[0].size(); ++i){
			for(uint j=0; j<debug.size(); ++j){
				of << debug[j][i] << "\t";
			}
			of << endl;
		}
		of.close();
	}
	cout << "******************************** "<<PACKAGE_STRING<<" ***********************************" << endl;
	return EXIT_SUCCESS;
}
