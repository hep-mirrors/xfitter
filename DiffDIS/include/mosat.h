/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
  \version 2.11.01
_____________________________________________________________*/

#ifndef MOSAT_H_
#define MOSAT_H_

// oooooooooooooooooooooooooo
class MoSat_t {
	double glu1[6], qua[6];
	double x0G,x0q, QQ0, lambda, b, sigma0;
	double Fac_x0G,Fac_x0q;
	double NNq, NormG, NormFull;
	const int Nc; ///< # colours
	const int Nf; ///< # flavours
  double AlphaQCD;
	
	// =============================================
	double SumChSqr(int nf) {
    double sq = 0;
    for(int i=1; i <= nf; i++) sq += (i & 1)? 1 : 4;
    return sq/9;
	}
  
public:
	
  // GBW fit:
  // x0=3.04e-4; 
  static const double x0_GBW = 3.04e-4;
  
	// =============================================
	MoSat_t(double alphas=0.4, double sxG=2.) : Nc(3), Nf(3) {
    Fac_x0q = 1;
		x0q = Fac_x0q*x0_GBW; 
    
    QQ0=1; lambda=0.288; 
		// b=6;
		b=7;
		sigma0 = 23.03/0.389;
    
		double NN0 = Nc*Nc*sigma0*sigma0/(4.0*pow(2*M_PI,6)*b);
    // --- this is N_G:
    NormG = NN0 *SumChSqr(Nf);
		NNq = NormG/(2*M_PI*Nc);
    NormFull = NN0*2;
    
    SetParams(alphas, sxG);
  }
	
	// =============================================
	void SetParams(double alphas, double sxG) {
    Fac_x0G = sxG;
		x0G = Fac_x0G*x0_GBW; 
		AlphaQCD = alphas;
	}
	
	// =============================================
	void SetAlphaQCD(double alphas) {
		AlphaQCD = alphas;
    // cout << "AlphaQCD = " << AlphaQCD << endl;
	}
	
	// =============================================
	void Setx0(double Vx0G=x0_GBW) {
		x0G = Vx0G; 
		// x0q = Vx0q; 
	}
	
	// =============================================
	void Setx0(double Vx0G, double Vx0q) {
    // if(Vx0q <= 0.) Vx0q = Vx0G;
		x0G = Vx0G; 
		x0q = Vx0q; 
	}
	
private:

	// =============================================
	double Phi_G(double beta) {
    // return 1/(1-beta);
    double blog = beta*log(beta);
    return 1 - 0.802*beta + (1.902 + 1.121*blog)*blog;
  }

	// ======================================
	void FillGlu1(double QQ, double xi) {
		double lnRQ2,RQ2, a;
		RQ2 = pow(xi/x0G,lambda)*QQ/QQ0;
		a = 1/RQ2;
		lnRQ2 = log(RQ2);
		glu1[0] = (-0.76663057e1 + 0.58068961e1 * lnRQ2) * a;
		glu1[1] = 0.58068962e1 * a;
		a /= RQ2;
		glu1[2] = (0.10071675e2 + 0.50265482e1 * lnRQ2) * a;
		glu1[3] = -(0.12158956e1 + (0.56098253e1 + 0.33510320e1 * lnRQ2) * lnRQ2) * a;
		a /= RQ2;
		glu1[4] = -(-0.11990366e2 + 0.46026271e1 * lnRQ2) * a;
		glu1[5] = -(-0.30217814e2 + 0.13807881e2 * lnRQ2) * a;
	}
	
	// ======================================
	void FillQua(double QQ, double xi, double beta) {
		double lnRQ2,RQ2, beta2, Lb;
		RQ2 = pow(xi/x0q,lambda)*QQ/QQ0;
		lnRQ2 = log(RQ2);
		beta2 = beta*beta;
		Lb = log(beta);
		double Cq = 1/RQ2;
		// Cq[1] = 0.1e1;
		// Cq[2] = beta * beta * pow(RQ2, -2);
		// Cq[3] = pow(beta, 0.3e1) * pow(RQ2, -2);
		// Cq[4] = pow(beta, 0.4e1) * pow(RQ2, -3);
		// Cq[5] = pow(beta, 0.4e1) * pow(RQ2, -3);
		
		qua[0] = Cq * (0.8751972557e2 * beta * Lb * Lb
									 + 0.1984401710e4 * (beta - 0.1004779751e1) * (beta * beta + 0.4779751104e-2 * beta + 0.4802597123e-2) * Lb 
									 - 0.9171958778e2 * (beta + 0.3738847857e1) * (beta - 0.8122589727e0) * (beta - 0.9999997186e0) 
									   * (beta2 + 0.1667152178e0 * beta + 0.2180196914e-1) * (beta2 - 0.4884241252e0 * beta + 0.2201927927e1) 
										 * (beta2 - 0.3824270028e1 * beta + 0.4391375051e1)
									);
		qua[1] = 0;
		Cq *= beta2/RQ2;
		qua[2] = 0.1984401708e4 * Cq * beta2 * (-1 + beta);
		Cq *= beta;
		qua[3] = Cq * (0.9922008534e3 * (beta - 0.4999929012e0) * (beta - 0.5000070988e0) * lnRQ2
               	   + 0.1991006016e5 * (beta + 0.6361585135e0) * (beta - 0.3233057258e-1) * Lb 
									 + 0.1523572274e4 * (beta + 0.1905660707e1) * (beta + 0.1569068433e-1) * (beta - 0.9887699436e0) * (beta - 0.8482298686e1)
									);
		Cq *= beta/RQ2;
		qua[4] = Cq * (0.2381282049e5 * beta * (beta - 0.4999999994e0) * (beta - 1) * lnRQ2 - 0.1151186161e5 * Lb * Lb 
									 - 0.5458995858e6 * Lb * beta 
									 + 0.3444545763e5 * (beta + 0.1075046457e1) * (beta + 0.1491985948e0) * (beta - 0.9999998675e0) * (beta - 0.2523348113e1) * (beta - 0.5438506575e1)
									);
		qua[5] = Cq * (-0.2381282048e5 * (beta - 0.2113248652e0) * (beta - 0.4999999997e0) * (beta - 0.7886751363e0) * lnRQ2
									 + (0.6521020294e6 - 0.8758869400e7 * beta) * Lb * Lb 
									 - 0.2381282048e6 * (beta - 0.3592637121e1) * (beta2 + 0.2092637120e1 * beta + 0.8118085799e1) * Lb 
									 - 0.2427821461e6 * (beta - 0.1000678465e1) * (beta - 0.1612517231e1) * (beta - 0.4637784774e1) * (beta2 + 0.3015142716e0 * beta + 0.1138825405e2)
									);
	}
	
	// =============================================
	double YL(double y) {
		double yb = 1 - y;
		return 2*yb/(1+yb*yb);
	}

public:

  /// Contributions to x sigma_red
	// =============================================
	double Ghi(double y, double QQ, double xi, double beta, bool all=false) {
		double cL = YL(y);
		FillGlu1(QQ,xi);
		double tw = glu1[2]+glu1[4] + cL*(glu1[3]+glu1[5]);
		if(all) tw += glu1[0] +cL*glu1[1];
		// return tw*AlphaQCD.val(QQ)*NormG*QQ *Phi_G(beta);
		return tw*AlphaQCD*NormG*QQ *Phi_G(beta);
		// return (tw > 0 && beta < 0.5) ? tw*AlphaQCD.val(QQ)*NormG*QQ/(1-beta) : 0;
	}

	// =============================================
	double Qhi(double y, double QQ, double xi, double beta, bool all=false) {
		double cL = YL(y);
		FillQua(QQ,xi,beta);
		double tw = qua[2]+qua[4] + 4*cL*(qua[3]+qua[5]); //--- 4 is from twist_L normalization
		if(all) tw += qua[0];
		return tw*NNq*QQ;
		// return (tw > 0 && beta < 0.5) ? tw*NNq*QQ : 0;
	}

	// =============================================
	double HiTwist(double y, double QQ, double xi, double beta) {
		return Ghi(y, QQ, xi, beta) + Qhi(y, QQ, xi, beta);
	}

  #if 0
	// =============================================
	void GetGlu(double QQ, double xi, double beta, double val[]) {
		double yb = 1 - QQ/(xi*sHERA)/beta;
		double cL = 2*yb/(1+yb*yb);
		// double A = AlphaQCD.val(QQ)*NormG*QQ *Phi_G(beta);
		double A = AlphaQCD*NormG*QQ *Phi_G(beta);
		// double A = 0.2*NormG*QQ *Phi_G(beta);
		FillGlu1(QQ,xi);
		val[0] = A*(glu1[0] + cL*glu1[1]);
		val[1] = A*(glu1[2] + cL*glu1[3]);
		val[2] = A*(glu1[4] + cL*glu1[5]);
	}

	// =============================================
	void GetQua(double QQ, double xi, double beta, double val[]) {
		double yb = 1 - QQ/(xi*sHERA)/beta;
		double cL = 2*yb/(1+yb*yb);
		double A = NNq*QQ;
		FillQua(QQ,xi,beta);
		val[0] = A*qua[0];
		val[1] = A*(qua[2] + 4*cL*qua[3]); //--- 4 is from twist_L normalization
		val[2] = A*(qua[4] + 4*cL*qua[5]); //--- 4 is from twist_L normalization
	}

	// =============================================
	void GetAll(double QQ, double xi, double beta, double val[]) {
		double gv[3];
    GetGlu(QQ,xi,beta,gv);
    GetQua(QQ,xi,beta,val);
		for(int j=0; j < 3; j++) val[j] += gv[j];
	}


	// =============================================
	double Full(double QQ, double xi, double beta) {
		double yb = 1 - QQ/(xi*sHERA)/beta;
		double cL = 2*yb/(1+yb*yb);
		// double A = NormFull*AlphaQCD.val(QQ)*QQ *Phi_G(beta);
		double A = NormFull*AlphaQCD*QQ *Phi_G(beta);
		double lnRQ,RQ;
		RQ = sqrt(pow(xi/x0G,lambda)*QQ/QQ0);
		lnRQ = log(RQ);
    double cF2sat = 0.67706 *exp(-.40307*RQ)+( 0.249485 + 0.194740*lnRQ)/RQ;
    double cFLsat = 0.223532*exp(-.40929*RQ)+(-0.012976 + 0.0604318*lnRQ)/RQ;
		double Qbox[3];
    GetQua(QQ,xi,beta,Qbox);
    return Qbox[0]+Qbox[1]+Qbox[2] + A*(cF2sat - (1-cL)*cFLsat);
	}
  #endif
	
  friend class dDIS_t;
  
};

#endif
