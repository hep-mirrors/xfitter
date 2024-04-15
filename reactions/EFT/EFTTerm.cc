   /*
     @file EFTTerm.cc
     @date 2023-03
     @author X.M. Shen <xmshen137@gmail.com>
   */

#include <string>
#include <cstring>
#include "yaml-cpp/yaml.h"
#include "EFTTerm.h"

using namespace std;

//------------------------------------------------------------------------------------
void EFTTerm::initParamName(vector<string> name_EFT_param_in){
  size_t i = 0;
  num_param = name_EFT_param_in.size();

  if (num_param > MAX_NUM_PARAM)
    hf_errlog(23040302, "F: too many EFT parameters");

  for (string name : name_EFT_param_in) {
    name_EFT_param.push_back(name);
    // find_EFT_param_id0.insert(std::make_pair(name, i));
    find_EFT_param_id1.insert(std::make_pair(name, ++i));
  }

  assert(i == num_param);
}

//------------------------------------------------------------------------------------
void EFTTerm::readInput(){
  num_bin = 0;  
  if (input_type == "fixed") 
    readFixedInput();
  else
    readMixedInput();
}

//------------------------------------------------------------------------------------
void EFTTerm::readFixedInput() {

  // read coefficients for all files
  for (string fname : filename_list) {

    YAML::Node coeff_node = YAML::LoadFile(fname);

    size_t num_bin_one_file = 0;

    // read linear coefficients
    for (size_t i=1; i <= num_param; i++) {
      string param_name = name_EFT_param[i-1];

      if (coeff_node[param_name]) {
	// check the number of bins
	if (num_bin_one_file == 0) {
	  num_bin_one_file = coeff_node[param_name].size();
	  num_bin += num_bin_one_file;
	} 
	else if (coeff_node[param_name].size() != num_bin_one_file) {
	  hf_errlog(23032903, "F: # coefficients != # bins");
	}

	// read the coeff.
	if (coeff.count(i) > 0) {
	  for (double val: coeff_node[param_name].as<std::vector<double> >() ) {
	    (*coeff[i]).push_back(val);
	  }
	} 
	else {
	  coeff.insert(std::make_pair(i, new vector<double>(coeff_node[param_name].as<std::vector<double> >() )));

	  // if (debug > 0) {
	  //   std::cout << "=======================================================" << std::endl;
	  //   std::cout << "EFTTerm.init: size of map coeff: " << coeff.size() << std::endl;
	  // }
	}
      } 
      else {
	hf_errlog(23032901, "F: EFT coefficients missing for: " + param_name);
      }
    } // end of reading linear coeff.

    std::cout << "EFT reaction: num_bin (from fixed input) = " << num_bin << std::endl;

    // read quadratic coefficients
    for (size_t i=1; i <= num_param; i++) {
      for (size_t j=i; j <= num_param; j++) {
	string param_name1 = name_EFT_param[i-1] + "*" + name_EFT_param[j-1];
	string param_name2 = name_EFT_param[j-1] + "*" + name_EFT_param[i-1];
	vector<double>* pvd;
	bool found = false;
	if (coeff_node[param_name1]) {
	  pvd = new vector<double>(coeff_node[param_name1].as<std::vector<double> >());
	  found = true;
	} else if (coeff_node[param_name2]) {
	  pvd = new vector<double>(coeff_node[param_name2].as<std::vector<double> >());
	  found = true;
	} else {
	  hf_errlog(23032904, "I: EFT coefficients missing for: " + param_name1);
	}

	if (found) {
	  if (coeff.count(i*100 + j) > 0) {
	    for (double val : (*pvd)) (*coeff[i*100 + j]).push_back(val);
	  } else {
	    coeff.insert(std::make_pair(i*100+j, pvd));
	  }
	}
      }
    } // end of reading quadratic coeff.

  } // loop over files
}

//------------------------------------------------------------------------------------
void EFTTerm::readMixedInput(){

  assert (filename_list.size() == 1);
  string fname = filename_list[0];
  YAML::Node node = YAML::LoadFile(fname);

  assert (node.Type() == YAML::NodeType::Map);

  string grid_dir = "/";
  bool save_grid_Q = true;

  // the info entry
  if (node["info"]) {

    // num_bin is now optional, can be inferred from the xsec entry (arrays or grids)
    if (node["info"]["num_bin"]) {
      num_bin = node["info"]["num_bin"].as<size_t>();
      std::cout << "EFT reaction: num_bin (from YAML file) = " << num_bin << std::endl;
    }
    // else
    //   hf_errlog(23060102, "F: EFT: `num_bin` not found in the `info` entry");

    // more arguments:
    if (node["info"]["grid_dir"])
      grid_dir = node["info"]["grid_dir"].as<string>();

    if (node["info"]["save_grid_in_memory"]) {
      save_grid_Q = node["info"]["save_grid_in_memory"].as<bool>();
      if (save_grid_Q)
	hf_errlog(23061507, "I: EFT: save grids in memory");
      else
	hf_errlog(23061508, "I: EFT: read grids for each iteration");
    }
    
    if (node["info"]["rows_before_transpose"]) {
      rows_before_transpose = node["info"]["rows_before_transpose"].as<int>();
      // if ( num_bin % rows_before_transpose != 0 )
      // 	hf_errlog(23062001, "S: rows_before_transpose does not divide num_bin");	
    }

    ///////////////////////////////////////////////////////
    // more arguments:
    if (node["info"]["scaleQ1"]) {
      scaleQ1 = node["info"]["scaleQ1"].as<bool>();      
    }

    if (node["info"]["scaleQ2"]) {
      scaleQ2 = node["info"]["scaleQ2"].as<bool>();      
    }

    if (scaleQ1) {
      if (node["info"]["scaling1"]) {
	vector<double> vec = node["info"]["scaling1"].as<vector<double> >();
	scaling1.resize(vec.size());
	std::copy(vec.begin(), vec.end(), std::begin(scaling1));
      }
      else 
	hf_errlog(23061901, "F: scaling xsec asked, but array `scaling1` not provided");
    }

    if (scaleQ2) {
      if (node["info"]["scaling2"]) {
	vector<double> vec = node["info"]["scaling2"].as<vector<double> >();
	scaling2.resize(vec.size());
	std::copy(vec.begin(), vec.end(), std::begin(scaling2));
      }
      else 
	hf_errlog(23070301, "F: scaling xsec asked, but array `scaling2` not provided");
    }

    ///////////////////////////////////////////////////////
    // more optional arguements:
    if (normQ) {
      if (node["info"]["binning_for_norm"]) {
	vector<double> vecN = node["info"]["binning_for_norm"].as<vector<double> >();
	binning_for_norm.resize(vecN.size());
	std::copy(vecN.begin(), vecN.end(), std::begin(binning_for_norm));
      }
      else {
	hf_errlog(24041503, "F: `binning_for_norm` is mandatory if normQ = True");
	// hf_errlog(23062601, "I: `binning_for_norm` not provided, assumed to be 1.");
	// if (num_bin > 0) {
	//   binning_for_norm.resize(num_bin);
	//   binning_for_norm = 1.0;
	// }
	// else
	//   hf_errlog(23062602, "F: num_bin should > 0.");
      }
    }

    ///////////////////////////////////////////////////////
    // only available through TermInfo now
    // if (node["info"]["debug"])
    //   debug = node["info"]["debug"].as<int>();

  }
  else {
    hf_errlog(23060101, "F: `info` entry missing in the EFT YAML file");
  } // end of info node

  /////////////////////////////////////////////////////////////////////////////
  // check all other entries
  for (YAML::const_iterator it=node.begin(); it!=node.end(); ++it ) {

    string entry_name = it->first.as<string>() ;
    if (entry_name == "info")
      continue;

    YAML::Node entry = it->second;
    assert (entry.Type() == YAML::NodeType::Map);

    if (entry["mask"])
      if (entry["mask"].as<bool>() == true)
	continue;

    string type;
    string param_name;
    if (entry["type"])
      type = entry["type"].as<string>(); // use string instead of char for future extension
    else
      hf_errlog(23052303, "F: type not given for entry " + entry_name );

    // 
    // if (num_bin <= 0) {
    //   num_bin = readNumBin(entry, entry_name);
    // }

    // C term
    if (type == "C") {
      assert (prvec_C == nullptr);
      prvec_C = new RawVec(entry, entry_name, num_bin, grid_dir, xi_ren, xi_fac, save_grid_Q);

      assert (prvec_C->format != "FR");
    }
    else if (type=="l" || type=="L" || type=="q" || type=="Q") {
      if (entry["param"]) {
	param_name = entry["param"].as<string>();

	if (find_EFT_param_id1.count(param_name) == 0) {
	  if (type == "l" || type == "L") {
	    hf_errlog(24041108, "I: EFT reaction: " + param_name + " is not fitted" );
	  }
	  continue;
	}

	size_t id = find_EFT_param_id1[param_name];
	if (basis.count(id) == 0) {
	  basis.insert(make_pair(id, new Vec(1)));
	}

	RawVec* prvec = new RawVec(entry, entry_name, num_bin, grid_dir, xi_ren, xi_fac, save_grid_Q);

	raw_basis.push_back(prvec);
	basis[id]->addIng(prvec, 1.0);
      }
      else {
	hf_errlog(23052401, "F: param name not given for " + entry_name );
      }
    } // end of l, q
    else if (type=="m" || type=="M") {
      if (entry["param"]) {
	vector<string> param_names = entry["param"].as<vector<string> >();
	
	assert (param_names.size() == 2);
	if (find_EFT_param_id1.count(param_names[0]) * find_EFT_param_id1.count(param_names[1]) == 0 ) 
	  continue;

	size_t i1 = find_EFT_param_id1[param_names[0]];
	size_t i2 = find_EFT_param_id1[param_names[1]];

	if (i1 > i2) {
	  i1 += i2;
	  i2 = i1 - i2;
	  i1 -= i1;
	}

	assert(basis.count(i1*100+i2) == 0);
	basis.insert(make_pair(i1*100+i2, new Vec(2)));
	
	RawVec* prvec = new RawVec(entry, entry_name, num_bin, grid_dir, xi_ren, xi_fac, save_grid_Q);

	raw_basis.push_back(prvec);
	basis[i1*100+i2]->addIng(prvec, 1.0);
      }
    } // end of m, M
    else {
      hf_errlog(23052402, "F: EFT: unknown type for entry: " + entry_name );
    }
  } // end of loop over entries

  //-------------------------------------------------------
  // initialize each vec and raw vec
  for (size_t i=1; i<=num_param; ++i) {
    if (basis.count(i) == 0) {
      hf_errlog(23052601, "F: linear term not found for " + name_EFT_param[i-1]);
    }
    initlq(i);
  }

  for (size_t i=1; i<num_param; ++i) {
    for (size_t j=i+1; j<=num_param; ++j) {
      initm(i, j);
    }
  }

  initrvec();

}
//------------------------------------------------------------------------------------


void EFTTerm::initlq(size_t i){
  Vec* pvecl = basis[i];
  assert (pvecl->ingredients.size() == 1 || pvecl->ingredients.size() == 2);

  if (pvecl->ingredients.size() == 1) {
    hf_errlog(23052602, "I: quadrtic term for " + name_EFT_param[i-1] + " not found");
    ingredient* ing = pvecl->ingredients[0];

    assert (ing->prvec->type == 1 || ing->prvec->type == -1);
    
    if (ing->prvec->type == -1) {
      double val = ing->prvec->param_val1;
      ing->coeff = 1.0 / val;
      pvecl->addIng(prvec_C, -1.0/val);
    }
  }
  else { // len==2
    Vec* pvecq = new Vec(2);
    basis.insert(make_pair(100*i+i, pvecq));

    RawVec* prvec1 = pvecl->ingredients[0]->prvec;
    RawVec* prvec2 = pvecl->ingredients[1]->prvec;    
    pvecl->ingredients.clear();

    int type1 = prvec1->type;
    int type2 = prvec2->type;

    // assert (solvablelq(type1, type2));
    if (type1 == 1) {
      if (type2 == 2)
	solvelq(pvecl, pvecq, prvec1, prvec2);
      else if (type2 == -2) 
	solvelQ(pvecl, pvecq, prvec1, prvec2);
    }
    else if (type2 == 1) {
      if (type1 == 2)
	solvelq(pvecl, pvecq, prvec2, prvec1);
      else if (type1 == -2) 
	solvelQ(pvecl, pvecq, prvec2, prvec1);
    }
    else
      solveNol(pvecl, pvecq, prvec1, prvec2);
  }
  
}

void EFTTerm::solvelq(Vec* pvecl, Vec* pvecq, RawVec* prvec1, RawVec* prvec2) {
  pvecl->addIng(prvec1, 1.0);
  pvecq->addIng(prvec2, 1.0);
}

void EFTTerm::solvelQ(Vec* pvecl, Vec* pvecq, RawVec* prvec1, RawVec* prvec2) {
  double val = prvec2->param_val1;

  pvecl->addIng(prvec1, 1.0);

  pvecq->addIng(prvec2, 1.0/val/val);
  pvecq->addIng(prvec1, -1.0/val);
  pvecq->addIng(prvec_C, -1.0/val/val);
}

void EFTTerm::solveNol(Vec* pvecl, Vec* pvecq, RawVec* prvec1, RawVec* prvec2) {
  int type1 = prvec1->type;
  int type2 = prvec2->type;

  double val1 = prvec1->param_val1;
  double val2 = prvec2->param_val1;

  vector<double> c1;
  vector<double> c2;
  lqQCoeff(c1, type1, val1);
  lqQCoeff(c2, type2, val2);

  double det = c1[1]*c2[2] - c2[1]*c1[2];
  if (det == 0.0)
    hf_errlog(23053101, "F: EFT: can not solve l/q for " + prvec1->param_name1);
  else {
    if (det < 0.0001 && det > -0.0001 ) {
      hf_errlog(23053102, "W: small det in solving l/q for " + prvec1->param_name1);
    }
    pvecl->addIng(prvec_C, (-c2[2]*c1[0] + c1[2]*c2[0]) / det );
    pvecl->addIng(prvec1, c2[2]/det);
    pvecl->addIng(prvec2, -c1[2]/det);

    pvecq->addIng(prvec_C, (c2[1]*c1[0] - c1[1]*c2[0]) / det );
    pvecq->addIng(prvec1, -c2[1]/det);
    pvecq->addIng(prvec2, c1[1]/det);
  }
}

void EFTTerm::lqQCoeff(vector<double>& c, int type, double val){
  if (type == -1) {
    c.push_back(1.0);
    c.push_back(val);
    c.push_back(0.0);
  }
  else if (type == 2){
    c.push_back(0.0);
    c.push_back(0.0);
    c.push_back(val*val);
  }
  else if (type == -2){
    c.push_back(1.0);
    c.push_back(val);
    c.push_back(val*val);
  }
}

//------------------------------------------------------------------------------------

int EFTTerm::initm(size_t i1, size_t i2){
  size_t im = i1 * 100 + i2;
  if (basis.count(im) == 0) {
    hf_errlog(23052603, "I: mixed term not found: " 
	      + name_EFT_param[i1-1] + "*" + name_EFT_param[i2-1]);
    return 0;
  }

  Vec* pvecm = basis[im];
  if (pvecm->ingredients[0]->prvec->type == 3)
    return 1;

  if (basis.count(i1*101) * basis.count(i2*101) == 0) {
    return 0;
  }
    
  RawVec* prvecM = pvecm->ingredients[0]->prvec;
  double c1 = prvecM->param_val1;
  double c2 = prvecM->param_val2;

  // M
  pvecm->ingredients[0]->coeff = 1.0/c1/c2;
  
  // C
  pvecm->addIng(prvec_C, -1.0/c1/c2);

  // l
  for (auto ing: basis[i1]->ingredients)
    pvecm->addIng(ing->prvec, ing->coeff * (-1.0 / c2));

  for (auto ing: basis[i2]->ingredients)
    pvecm->addIng(ing->prvec, ing->coeff * (-1.0 / c1));

  // q
  for (auto ing: basis[i1*101]->ingredients)
    pvecm->addIng(ing->prvec, ing->coeff * (-1.0 * c1 / c2));

  for (auto ing: basis[i2*101]->ingredients)
    pvecm->addIng(ing->prvec, ing->coeff * (-1.0 * c2 / c1));

  return 1;
}

//------------------------------------------------------------------------------------
void EFTTerm::initrvec(){

  if (prvec_C == nullptr)
    hf_errlog(23053104, "F: central value not found in the mixed input file ");

  if (prvec_C->format == "FA") {
    for (auto prvec: raw_basis) {

      if (prvec->format == "FR") {
	prvec->FR2FA(prvec_C->value_list);
      }
    }
  }

}

//------------------------------------------------------------------------------------
void EFTTerm::initIter(valarray<double>& list_val){
  // initialization for each term
  setValEFT(list_val);
  
  if (input_type == "mixed") {
    updatervec();
    book();
  }
}

//------------------------------------------------------------------------------------
void EFTTerm::updatervec() {
  // FR -> FA
  if (prvec_C->format != "FR" && prvec_C->format != "FA") {
    prvec_C->convolute();

    for (auto prvec: raw_basis) {
      if (prvec->format == "FR") {
	prvec->FR2FA(prvec_C->value_list);
      }
    }
  }

  // convolution
  for (auto prvec: raw_basis) {
    if (prvec->format != "FR" && prvec->format != "FA") {
      prvec->convolute();
    }
  }
} 

//------------------------------------------------------------------------------------
void EFTTerm::book() {

  if (no_central)
    prvec_C->setCoeff(0.0);
  else
    prvec_C->setCoeff(1.0);

  for (auto prvec: raw_basis)
    prvec->setCoeff(0.0);

  for (size_t i=0; i < num_param; ++i) {
    // l
    double vi = val_EFT_param[i];
    basis[i+1]->book(vi);

    // q
    if (basis.count((i+1)*101) > 0)  {
      basis[(i+1)*101]->book(vi*vi);
    }
    
    // m
    for (size_t j=i+1; j < num_param; ++j) {
      if (basis.count((i+1)*100+j+1) > 0) {
	basis[100*(i+1)+j+1]->book(vi * val_EFT_param[j]);
      }
    }
  }
}

//------------------------------------------------------------------------------------
void EFTTerm::setValEFT(valarray<double>& list_val) {
  // executed for each computation

  for (size_t i=0; i<num_param; i++)
    val_EFT_param[i] = list_val[i];

  if (debug > 2) {
    std::cout << "=======================================================" << std::endl;
    std::cout << "EFTTerm.setValEFT: current values of EFT param:" << std::endl;
    for (size_t i=0; i<num_param; i++) {
      std::cout << name_EFT_param[i] << " = " <<  val_EFT_param[i] << std::endl;
    }
  }
};

//------------------------------------------------------------------------------------
void EFTTerm::calcXSec(valarray<double>& xsec) {

  assert (xsec.size() == num_bin);

  if (input_type == "fixed")
    calcXSecFixed(xsec);
  else
    calcXSecMixed(xsec);

  if (scaleQ1)
    scaleXSec1(xsec);

  if (scaleQ2)
    scaleXSec2(xsec);

  if (normQ)
    normXSec(xsec);

  if (rows_before_transpose > 1) {
    transpose(xsec);
  }

};

//------------------------------------------------------------------------------------
void EFTTerm::transpose(valarray<double>& xsec) {
  // tranpose xsec as if it is a matrix with "rows_before_transpose" rows
  int m = rows_before_transpose;
  int n = num_bin / rows_before_transpose;

  valarray<double> tmp(xsec.size());

  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      tmp[j*m + i] = xsec[i*n + j];

  xsec = std::move(tmp);
}

//------------------------------------------------------------------------------------
void EFTTerm::scaleXSec1(valarray<double>& xsec) {
  assert (scaling1.size() == num_bin);

  xsec = xsec * scaling1;
}

void EFTTerm::scaleXSec2(valarray<double>& xsec) {
  assert (scaling2.size() == num_bin);

  xsec = xsec * scaling2;
}

//------------------------------------------------------------------------------------
void EFTTerm::normXSec(valarray<double>& xsec) {
  assert (binning_for_norm.size() == num_bin);

  xsec /= (xsec * binning_for_norm).sum();
}

//------------------------------------------------------------------------------------
void EFTTerm::calcXSecMixed(valarray<double>& xsec) {

  xsec = 0.0;

  prvec_C->increaseXSecInPlace(xsec);

  for (auto prvec : raw_basis) {
    prvec->increaseXSecInPlace(xsec);
  }

  if (abs_output == false)
    for (size_t i=0; i<num_bin; ++i)
      xsec[i] /= prvec_C->value_list[i]; // value_list is a vector, not a valarray; so can only be divided bin by bin

};

//------------------------------------------------------------------------------------
void EFTTerm::calcXSecFixed(valarray<double>& xsec) {

  // C
  if (no_central) {
    xsec = 0.0;
  }
  else {
    xsec = 1.0;
  }

  // l
  for (size_t i=0; i < num_param; i++) {
    if (coeff.count(i+1) > 0) {

      for (size_t k=0; k < num_bin; k++) {
	xsec[k] += (*(coeff[i+1]))[k] * val_EFT_param[i];
      }
      /////////////////////////////////
      // if (debug > 10) {
      // 	std::cout << "=======================================================" << std::endl;
      // 	std::cout << "EFTTerm.calcXSec: l:" + name_EFT_param[i] + " = " << val_EFT_param[i] << std::endl;
      // 	std::cout << "EFTTerm.calcXSec: new xsec=" << std::endl;
      // 	for (double v : xsec) {
      // 	  std::cout << v << ", ";
      // 	}
      // 	std::cout << std::endl;
      // }
      /////////////////////////////////
    } 
    else {
      hf_errlog(23051601, "F: linear coefficients missing for: " + name_EFT_param[i]);
    }
  
  }

  // q&m
  for (size_t i=0; i < num_param; i++) {
    for (size_t j=i; j < num_param; j++) {
      if (coeff.count((i+1)*100+j+1) > 0) {

	vector<double>* pvec = coeff[(i+1)*100+j+1];
	for (size_t k=0; k < num_bin; k++)
	  xsec[k] += (*pvec)[k] * val_EFT_param[i] * val_EFT_param[j];
	/////////////////////////////////
	// if (debug > 10) {
	//   std::cout << "=======================================================" << std::endl;
	//   std::cout << "EFTTerm.calcXSec: q/m: " << name_EFT_param[i] +  "*" + name_EFT_param[j] 
	// 	    << " = " << val_EFT_param[i] * val_EFT_param[j] << std::endl;
	//   std::cout << "EFTTerm.calcXSec:" << std::endl;
	//   for (double v : xsec) {
	//     std::cout << v << ", ";
	//   }
	//   std::cout << std::endl;
	// }
	/////////////////////////////////
      }
    }
  } // end of q&m

  // return xsec;
}

