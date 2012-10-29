/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2011--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include "matrix.h"

//========================
double Vector::ScalarProduct(Vector& v) {
  if(v.len != len) throw("Vector::ScalarProduct: diff. lengths");
  double s=0;
  for(int r=first; r <= last; r++) s += vp[r]*v[r-first+v.first];
  return s;
}

//========================
void Vector::Lmul(SqMatrix_t& MM) {
  SqMatrix_t M(first, len);
  Vector V(first,len);
  M.Set(MM);
  for(int r=first; r <= last; r++) {
    for(int j=first; j <= last; j++) V[r] += M[r][j]*vp[j];
  }
  Set(V);
}

//==================================
void SqMatrix_t::GetCol(int k, Vector& c) {
  for(int row=first; row <= last; row++) c.vp[row] = mp[row][k];
}

//==================================
Vector& SqMatrix_t::GetCol(int k) {
  Vector* c = new Vector(first, len);
  for(int row=first; row <= last; row++) c->vp[row] = mp[row][k];
  return *c;
}

//==================================
double SqMatrix_t::MxElem(Vector& v) {
  Vector c(v);
  c.Lmul(*this);
  return v.ScalarProduct(c);
}

//==================================================
void MM_mul(const SqMatrix_t& A, const SqMatrix_t& B, SqMatrix_t& C) {
  //Vector* b=new Vector(a.first, a.last);
  for(int r=A.first; r <= A.last; r++)
  for(int k=A.first; k <= A.last; k++) {
    double s=0;
    for(int j=A.first; j <= A.last; j++) s += A[r][j]*B[j][k];
    C[r][k]=s;
  }
}

//==================================================
SqMatrix_t operator+(const SqMatrix_t& A, const SqMatrix_t& B) {
  SqMatrix_t C(A);
  for(int r=A.first; r <= A.last; r++)
  for(int k=A.first; k <= A.last; k++) {
    C[r][k] += B[r][k];
  }
  return C;
}

//============================================
SqMatrix_t& SqMatrix_t::Rmul(const SqMatrix_t& B) {
  SqMatrix_t C(first, len);
  MM_mul(*this, B, C);
  Set(C);
  return *this;
}

//=============================================
SqMatrix_t& SqMatrix_t::Inverse() {
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
	int astart=first;
  StartIndex(1);
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;
	int n = len;

	indxc=new int[n+1];
	indxr=new int[n+1];
	ipiv=new int[n+1];
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(mp[j][k]) >= big) {
							big=fabs(mp[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) throw("gaussj: Singular Matrix-1");
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(mp[irow][l],mp[icol][l])
			//for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (mp[icol][icol] == 0.0) throw("gaussj: Singular Matrix-2");
		pivinv=1.0/mp[icol][icol];
		mp[icol][icol]=1.0;
		for (l=1;l<=n;l++) mp[icol][l] *= pivinv;
		//for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=mp[ll][icol];
				mp[ll][icol]=0.0;
				for (l=1;l<=n;l++) mp[ll][l] -= mp[icol][l]*dum;
				//for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(mp[k][indxr[l]],mp[k][indxc[l]]);
	}
	delete[] ipiv;
	delete[] indxr;
	delete[] indxc;
  StartIndex(astart);
	return *this;
#undef SWAP
}

//=============================================
int SqMatrix_t::Diagonalize(Vector& d) {
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

	int j,iq,ip,i, nrot;
	//int n = len, astart=first, dstart=d.first;
	int n = len, dstart=d.first;
	double tresh,theta,tau,t,sm,s,h,g,c;
	Vector b(1,n);
	Vector z(1,n);
  SqMatrix_t v(1,n);
  SqMatrix_t a(1,n); a.Set(*this);
  //StartIndex(1);
  d.StartIndex(1);

	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
      //--- a is now lower triangular
      //StartIndex(astart);
      //v.StartIndex(vstart);
      Set(v);
      d.StartIndex(dstart);
			return nrot;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
					&& (fabs(d[iq])+g) == fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((fabs(h)+g) == fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++nrot;
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	cout << "Too many iterations in routine jacobi"<<endl;
#undef ROTATE
}

  //====================================
  ostream& operator<<(ostream &ostr, const SqMatrix_t& a) {
    a.Write(ostr);
    return ostr << endl;
  }

