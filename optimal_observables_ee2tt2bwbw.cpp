const Double_t pi = 3.141592653589793;
Double_t pi2 = pi*pi;

Double_t mod(TComplex z){
	Double_t r,i;
	r = z.Re();
	i = z.Im();
	return r*r+i*i;
};

Double_t Lam = 1000.;
Double_t m = 172.5;
Double_t gz = 2.4952  ;
Double_t gw = 2.085  ;
Double_t mw = 79.82436 ;
Double_t mz = 91.1876 ;
Double_t gf = 1.16637e-5 ;
Double_t sW2 = 1-pow(mw/mz,2) ;
Double_t sW= sqrt(sW2); Double_t cW2 = 1-sW2; Double_t cW = sqrt(cW2);
Double_t apM1 = pi/(sqrt(2)*gf*mw*mw*sW2) ;
Double_t ap =1./apM1  ;
Double_t aps = 0.1076705;
Double_t pb =pow(197.3269631e-18,2)/1e-40 ;

void constants(){
	cout << "mt  = " << m   << endl;
	cout << "mz  = " << mz  << endl;
	cout << "mw  = " << mw  << endl;
	cout << "sw2 = " << sW2 << endl;
	cout << "gf  = " << gf  << endl;
	cout << "apM1= " << apM1 << endl;
}


Double_t pre = ap*pi/sW2/cW2;
Double_t pre2 = pre*pre;

Double_t LamM2 = 1./Lam/Lam;
Double_t LamM4 = LamM2*LamM2;

Double_t m2 = m*m;
Double_t w = mw/m; Double_t w2 = w*w;
Double_t ap2 = ap*ap;

TComplex i = TComplex(0.,1.);
const Double_t Sleft(Int_t k, Double_t sqs, const Double_t t0, const Double_t t1, const Double_t f1, const Double_t t2, const Double_t f2) {

Double_t s, s2, s4, r, r2, g2, g, b2, b, ipropz, rpropz, mpropz, prefac,
	mymodleft, mymodright, myreleft, myreleftbis, myreright,
	expr;
TComplex propz;

s = sqs*sqs; s2 = s*s; s4 = s2*s2;
r = m2/s; r2 = r*r; g2 =1./(4*r); g=sqrt(g2);
b2 = 1-4*r; b = sqrt(b2);
prefac = pb/(2*s) /(8*pi) *b /4. *1e3;

propz = 1./(s-mz*mz +i*gz*mz);
ipropz = propz.Im();
rpropz = propz.Re();
mpropz = rpropz*rpropz+ipropz*ipropz;

mymodleft = mod(16*sW2*cW2 + (3-8*sW2)*(1-2*sW2)*s*propz);
myreleft  = 16*rpropz*sW2*cW2 + (3-8*sW2)*(1-2*sW2)*s*mpropz;
myreleftbis = 16*sW2*cW2 + (3-8*sW2)*(1-2*sW2)*s*rpropz;

mymodright = mod(8*cW2 - (3-8*sW2)*propz*s);
myreright = 8*cW2*rpropz - (3-8*sW2)*s*mpropz;

//sm
if (k==0){
	expr =
	- 1./6.*(mymodleft - 9*b2*pow((1-2*sW2),2)*mpropz*s2)*pow(sin(t0),2)*sin(t1)*sin(t2)*cos(f1 - f2)*pow((1-2*w2),2)*pow((1+2*w2),-2)*pre2
	+ 1./6.*pow(sin(t0),2)*mymodleft*pow(g2,-1)*pre2
	- 16*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 - f2)*pre*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*s*b*ipropz*pi
	+ 2*cos(t0)*(1-2*sW2)*s*b*pre2*myreleft
	+ 2*cos(t0)*cos(t1)*cos(t2)*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*s*b*pre2*myreleft
	+ 1./6.*(1 + pow(cos(t0),2))*(mymodleft + 9*b2*pow((1-2*sW2),2)*mpropz*s2)*pre2
	+ 1./6.*(1 + pow(cos(t0),2))*(mymodleft + 9*b2*pow((1-2*sW2),2)*mpropz*s2)*cos(t1)*cos(t2)*pow((1-2*w2),2)*pow((1+2*w2),-2)*pre2
	- (1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*s*b*pre2*myreleft
	+ 16*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(t0)*pre*ap*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*s*pow(g,-1)*b*ipropz*pi
	- 8*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(2*t0)*pre*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*s*pow(g,-1)*b*ipropz*pi
	- 1./6.*(sin(t1)*sin(t2)*cos(f1 + f2) + cos(t1)*cos(t2))*pow(sin(t0),2)*pow((1-2*w2),2)*pow((1+2*w2),-2)*mymodleft*pow(g2,-1)*pre2
	+ 1./2.*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(2*t0)*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*s*pow(g,-1)*b*pre2*myreleft
	+ 1./3.*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(t0)*(1-2*w2)*pow((1+2*w2),-1)*pow(g,-1)*mymodleft*pre2
	- 1./6.*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(2*t0)*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow(g,-1)*mymodleft*pre2
	;
	expr +=  - (sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(t0)*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*s*pow(g,-1)*b*pre2*myreleft
	- 1./3.*(cos(t1) + cos(t2))*(mymodleft + 9*b2*pow((1-2*sW2),2)*mpropz*s2)*cos(t0)*(1-2*w2)*pow((1+2*w2),-1)*pre2
	;
}//Re(lqA)
else if (k==1){
	expr =
	+ pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 - f2)*pre*pow((1-2*w2),2)*pow((1+2*w2),-2)*(3-8*sW2)*(1-2*sW2)*b*ipropz*LamM2*s2
	+ 3*pow(sin(t0),2)*sin(t1)*sin(t2)*cos(f1 - f2)*pre*b2*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*rpropz*LamM2*s2
	+ 2*cos(t0)*cos(t1)*cos(t2)*pre*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*b*LamM2*myreleftbis
	+ 2*cos(t0)*pre*s*b*LamM2*myreleftbis
	+ 3*(1 + pow(cos(t0),2))*cos(t1)*cos(t2)*pre*b2*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*rpropz*LamM2*s2
	+ 3*(1 + pow(cos(t0),2))*pre*b2*(1-2*sW2)*rpropz*LamM2*s2
	- (1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*pre*(1-2*w2)*pow((1+2*w2),-1)*s*b*LamM2*myreleftbis
	- (sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(t0)*pre*(1-2*w2)*pow((1+2*w2),-1)*(3-8*sW2)*(1-2*sW2)*pow(g,-1)*b*ipropz*LamM2*s2
	+ 1./2.*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(2*t0)*pre*pow((1-2*w2),2)*pow((1+2*w2),-2)*(3-8*sW2)*(1-2*sW2)*pow(g,-1)*b*ipropz*LamM2*s2
	+ 1./2.*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(2*t0)*pre*(1-2*w2)*pow((1+2*w2),-1)*s*pow(g,-1)*b*LamM2*myreleftbis
	- (sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(t0)*pre*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*pow(g,-1)*b*LamM2*myreleftbis
	- 6*(cos(t1) + cos(t2))*cos(t0)*pre*b2*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*rpropz*LamM2*s2
	;
}//Re(eqA)
else if (k==2){
	expr = 0;

}//Re(pqA)
else if (k==3){
	expr =
	- 32*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 - f2)*pre*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*b*ipropz*pi*r*LamM2*s2
	+ 6*pow(sin(t0),2)*sin(t1)*sin(t2)*cos(f1 - f2)*b2*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow((1-2*sW2),2)*s*mpropz*r*LamM2*s2*pre2
	+ 4*cos(t0)*(1-2*sW2)*b*r*LamM2*s2*pre2*myreleft
	+ 4*cos(t0)*cos(t1)*cos(t2)*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*b*r*LamM2*s2*pre2*myreleft
	+ 6*(1 + pow(cos(t0),2))*cos(t1)*cos(t2)*b2*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow((1-2*sW2),2)*s*mpropz*r*LamM2*s2*pre2
	+ 6*(1 + pow(cos(t0),2))*b2*pow((1-2*sW2),2)*s*mpropz*r*LamM2*s2*pre2
	- 2*(1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*b*r*LamM2*s2*pre2*myreleft
	+ 32*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(t0)*pre*ap*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*pow(g,-1)*b*ipropz*pi*r*LamM2*s2
	- 16*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(2*t0)*pre*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*pow(g,-1)*b*ipropz*pi*r*LamM2*s2
	+ (sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(2*t0)*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*pow(g,-1)*b*r*LamM2*s2*pre2*myreleft
	- 2*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(t0)*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*pow(g,-1)*b*r*LamM2*s2*pre2*myreleft
	- 12*(cos(t1) + cos(t2))*cos(t0)*b2*(1-2*w2)*pow((1+2*w2),-1)*pow((1-2*sW2),2)*s*mpropz*r*LamM2*s2*pre2
	;
}//Re(lqV)
else if (k==4){
	expr =
	+ 3*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 - f2)*pre*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*b*ipropz*LamM2*s2
	+ pow(sin(t0),2)*sin(t1)*sin(t2)*cos(f1 - f2)*pre*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*LamM2*myreleftbis
	- pow(sin(t0),2)*pre*s*LamM2*pow(g2,-1)*myreleftbis
	- 6*cos(t0)*cos(t1)*cos(t2)*pre*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*b*rpropz*LamM2*s2
	- 6*cos(t0)*pre*(1-2*sW2)*b*rpropz*LamM2*s2
	- (1 + pow(cos(t0),2))*cos(t1)*cos(t2)*pre*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*LamM2*myreleftbis
	- (1 + pow(cos(t0),2))*pre*s*LamM2*myreleftbis
	+ 3*(1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*pre*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*b*rpropz*LamM2*s2
	- 3*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(t0)*pre*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*pow(g,-1)*b*ipropz*LamM2*s2
	+ 3./2.*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(2*t0)*pre*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*pow(g,-1)*b*ipropz*LamM2*s2
	+ (sin(t1)*sin(t2)*cos(f1 + f2) + cos(t1)*cos(t2))*pow(sin(t0),2)*pre*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*LamM2*pow(g2,-1)*myreleftbis
	- 3./2.*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(2*t0)*pre*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*pow(g,-1)*b*rpropz*LamM2*s2
	- 2*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(t0)*pre*(1-2*w2)*pow((1+2*w2),-1)*s*pow(g,-1)*LamM2*myreleftbis
	+ (sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(2*t0)*pre*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*pow(g,-1)*LamM2*myreleftbis
	;
	expr +=  + 3*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(t0)*pre*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*pow(g,-1)*b*rpropz*LamM2*s2
	+ 2*(cos(t1) + cos(t2))*cos(t0)*pre*(1-2*w2)*pow((1+2*w2),-1)*s*LamM2*myreleftbis
	;
}//Re(eqV)
else if (k==5){
	expr = 0;

}//Re(pqV)
else if (k==6){
	expr =
	- 2*pow(sin(t0),2)*(1-2*sW2)*r*LamM2*s2*pow(g2,-1)*pre2*myreleft
	+ 2*pow(sin(t0),2)*sin(t1)*sin(t2)*cos(f1 - f2)*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*r*LamM2*s2*pre2*myreleft
	- 12*cos(t0)*pow((1-2*sW2),2)*s*b*mpropz*r*LamM2*s2*pre2
	- 12*cos(t0)*cos(t1)*cos(t2)*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow((1-2*sW2),2)*s*b*mpropz*r*LamM2*s2*pre2
	- 2*(1 + pow(cos(t0),2))*(1-2*sW2)*r*LamM2*s2*pre2*myreleft
	- 2*(1 + pow(cos(t0),2))*cos(t1)*cos(t2)*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*r*LamM2*s2*pre2*myreleft
	+ 6*(1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*(1-2*w2)*pow((1+2*w2),-1)*pow((1-2*sW2),2)*s*b*mpropz*r*LamM2*s2*pre2
	+ 2*(sin(t1)*sin(t2)*cos(f1 + f2) + cos(t1)*cos(t2))*pow(sin(t0),2)*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*r*LamM2*s2*pow(g2,-1)*pre2*myreleft
	- 3*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(2*t0)*(1-2*w2)*pow((1+2*w2),-1)*pow((1-2*sW2),2)*s*pow(g,-1)*b*mpropz*r*LamM2*s2*pre2
	- 4*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(t0)*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*pow(g,-1)*r*LamM2*s2*pre2*myreleft
	+ 2*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(2*t0)*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*pow(g,-1)*r*LamM2*s2*pre2*myreleft
	+ 6*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(t0)*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow((1-2*sW2),2)*s*pow(g,-1)*b*mpropz*r*LamM2*s2*pre2
	+ 4*(cos(t1) + cos(t2))*cos(t0)*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*r*LamM2*s2*pre2*myreleft
	;
}//Re(uZ)
else if (k==7){
	expr =
	+ 16*(1 + 2*w2)*pow(sin(t0),2)*pre*ap*pow((1+2*w2),-1)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*pi*r*LamM2*s2*myreleft
	+ 96*(1 + 2*w2)*cos(t0)*pre*ap*pow((1+2*w2),-1)*pow((1-2*sW2),2)*pow(sW,-1)*pow(cW,-1)*s*b*mpropz*pi*r*LamM2*s2
	+ 32*(64*rpropz*sW2*cW2*w2 - 64*rpropz*sW2*cW2*pow(w2,2) + 3*(1-2*sW2)*s*mpropz - 12*(1-2*sW2)*s*mpropz*pow(w2,2) + 4*(3-8*sW2)*(1-2*sW2)*s*mpropz*w2 - 4*(3-8*sW2)*(1-2*sW2)*s*mpropz*pow(w2,2))*cos(t0)*cos(t1)*cos(t2)*pre*ap*(1-2*w2)*pow(
	(1+2*w2),-3)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*b*pi*r*LamM2*s2
	- 16./3.*(2*mymodleft*w2 - 2*mymodleft*pow(w2,2) + 48*(1-2*sW2)*s*rpropz*sW2*cW2 - 192*(1-2*sW2)*s*rpropz*sW2*cW2*pow(w2,2) + 3*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2 - 12*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2*pow(w2,2) - 18*b2*pow((1-2*sW2),2)*
	mpropz*s2*w2 + 18*b2*pow((1-2*sW2),2)*mpropz*s2*pow(w2,2))*pow(sin(t0),2)*sin(t1)*sin(t2)*cos(f1 - f2)*pre*ap*(1-2*w2)*pow((1+2*w2),-3)*pow(sW,-1)*pow(cW,-1)*s*pi*r*LamM2
	- 1024*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 - f2)*ap2*(1-w2)*(1-2*w2)*pow((1+2*w2),-3)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*b*ipropz*r*LamM2*pi2*s2*w2
	+ 16*(1 + pow(cos(t0),2))*(1 + 2*w2)*pre*ap*pow((1+2*w2),-1)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*pi*r*LamM2*s2*myreleft
	+ 16./3.*(1 + pow(cos(t0),2))*(2*mymodleft*w2 - 2*mymodleft*pow(w2,2) + 48*(1-2*sW2)*s*rpropz*sW2*cW2 - 192*(1-2*sW2)*s*rpropz*sW2*cW2*pow(w2,2) + 3*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2 - 12*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2*pow(w2,2) + 18*b2*
	pow((1-2*sW2),2)*mpropz*s2*w2 - 18*b2*pow((1-2*sW2),2)*mpropz*s2*pow(w2,2))*cos(t1)*cos(t2)*pre*ap*(1-2*w2)*pow((1+2*w2),-3)*pow(sW,-1)*pow(cW,-1)*s*pi*r*LamM2
	- 16*(1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*(32*rpropz*sW2*cW2*w2 - 32*rpropz*sW2*cW2*pow(w2,2) + 3*(1-2*sW2)*s*mpropz - 12*(1-2*sW2)*s*mpropz*pow(w2,2) + 2*(3-8*sW2)*(1-2*sW2)*s*mpropz*w2 - 2*(3-8*sW2)*(1-2*sW2)*s*mpropz*pow(w2,2))*pre*ap*pow(
	(1+2*w2),-2)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*b*pi*r*LamM2*s2
	+ 128*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(2*t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*(1-g2)*pow(sW,-1)*pow(cW,-1)*pow(g,-1)*ipropz*r*LamM2*pi2*s2
	+ 512*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(t0)*ap2*(1-w2)*pow((1+2*w2),-2)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*pow(g,-1)*b*ipropz*r*LamM2*pi2*s2*w2
	;
	expr +=  - 512*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(2*t0)*ap2*(1-w2)*(1-2*w2)*pow((1+2*w2),-3)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*pow(g,-1)*b*ipropz*r*LamM2*pi2*s2*w2
	- 256*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*(1-g2)*pow(sW,-1)*pow(cW,-1)*pow(g,-1)*ipropz*r*LamM2*pi2*s2
	- 16./3.*(sin(t1)*sin(t2)*cos(f1 + f2) + cos(t1)*cos(t2))*(2*mymodleft*w2 - 2*mymodleft*pow(w2,2) + 48*(1-2*sW2)*s*rpropz*sW2*cW2*g2 - 192*(1-2*sW2)*s*rpropz*sW2*cW2*g2*pow(w2,2) + 3*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2*g2 - 12*(3-8*sW2)*pow(
	(1-2*sW2),2)*mpropz*s2*g2*pow(w2,2))*pow(sin(t0),2)*pre*ap*(1-2*w2)*pow((1+2*w2),-3)*pow(sW,-1)*pow(cW,-1)*s*pi*r*LamM2*pow(g2,-1)
	+ 8*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*(32*rpropz*sW2*cW2*w2 - 32*rpropz*sW2*cW2*pow(w2,2) + 3*(1-2*sW2)*s*mpropz*g2 - 12*(1-2*sW2)*s*mpropz*g2*pow(w2,2) + 2*(3-8*sW2)*(1-2*sW2)*s*mpropz*w2 - 2*(3-8*sW2)*(1-2*sW2)*s*mpropz*pow(w2,2))*sin(2*t0)*
	pre*ap*pow((1+2*w2),-2)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*pow(g,-1)*b*pi*r*LamM2*s2
	+ 16./3.*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*(2*mymodleft*w2 - 2*mymodleft*pow(w2,2) + 48*(1-2*sW2)*s*rpropz*sW2*cW2 - 192*(1-2*sW2)*s*rpropz*sW2*cW2*pow(w2,2) + 48*(1-2*sW2)*s*rpropz*sW2*cW2*g2 - 192*(1-2*sW2)*s*rpropz*sW2*cW2*g2*pow(w2,2) + 3*
	(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2 - 12*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2*pow(w2,2) + 3*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2*g2 - 12*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2*g2*pow(w2,2))*sin(t0)*pre*ap*pow((1+2*w2),-2)*pow(sW,-1)*pow(cW,-1)*s*
	pow(g,-1)*pi*r*LamM2
	- 16*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*(64*rpropz*sW2*cW2*w2 - 64*rpropz*sW2*cW2*pow(w2,2) + 3*(1-2*sW2)*s*mpropz*g2 - 12*(1-2*sW2)*s*mpropz*g2*pow(w2,2) + 4*(3-8*sW2)*(1-2*sW2)*s*mpropz*w2 - 4*(3-8*sW2)*(1-2*sW2)*s*mpropz*pow(
	w2,2))*sin(t0)*pre*ap*(1-2*w2)*pow((1+2*w2),-3)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*pow(g,-1)*b*pi*r*LamM2*s2
	- 8./3.*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*(4*mymodleft*w2 - 4*mymodleft*pow(w2,2) + 48*(1-2*sW2)*s*rpropz*sW2*cW2 - 192*(1-2*sW2)*s*rpropz*sW2*cW2*pow(w2,2) + 48*(1-2*sW2)*s*rpropz*sW2*cW2*g2 - 192*(1-2*sW2)*s*rpropz*sW2*cW2*g2*
	pow(w2,2) + 3*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2 - 12*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2*pow(w2,2) + 3*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2*g2 - 12*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2*g2*pow(w2,2))*sin(2*t0)*pre*ap*(1-2*w2)*pow((1+2*w2),-3)*
	pow(sW,-1)*pow(cW,-1)*s*pow(g,-1)*pi*r*LamM2
	- 32./3.*(cos(t1) + cos(t2))*(mymodleft*w2 - mymodleft*pow(w2,2) + 48*(1-2*sW2)*s*rpropz*sW2*cW2 - 192*(1-2*sW2)*s*rpropz*sW2*cW2*pow(w2,2) + 3*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2 - 12*(3-8*sW2)*pow((1-2*sW2),2)*mpropz*s2*pow(w2,2) + 9*b2*pow(
	(1-2*sW2),2)*mpropz*s2*w2 - 9*b2*pow((1-2*sW2),2)*mpropz*s2*pow(w2,2))*cos(t0)*pre*ap*pow((1+2*w2),-2)*pow(sW,-1)*pow(cW,-1)*s*pi*r*LamM2
	;
}//Re(uA)
else if (k==8){
	expr =
	+ 32*(1 + 2*w2)*pow(sin(t0),2)*pre*ap*pow((1+2*w2),-1)*s*pi*r*LamM2*myreleftbis
	+ 192*(1 + 2*w2)*cos(t0)*pre*ap*pow((1+2*w2),-1)*(1-2*sW2)*b*rpropz*pi*r*LamM2*s2
	- 32*(3 - 12*pow(w2,2) + 32*sW2*w2 - 32*sW2*pow(w2,2))*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 - f2)*pre*ap*(1-2*w2)*pow((1+2*w2),-3)*(1-2*sW2)*b*ipropz*pi*r*LamM2*s2
	- 32./3.*(48*sW2*pow(cW2,2) - 192*sW2*pow(cW2,2)*pow(w2,2) + mymodleft*w2 - mymodleft*pow(w2,2) + 3*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2 - 12*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2*pow(w2,2) - 9*b2*pow((1-2*sW2),2)*mpropz*s2*w2 + 9*b2*pow((1-2*sW2),2)*
	mpropz*s2*pow(w2,2))*pow(sin(t0),2)*sin(t1)*sin(t2)*cos(f1 - f2)*pre*ap*(1-2*w2)*pow((1+2*w2),-3)*s*pi*r*pow(cW2,-1)*LamM2
	+ 64*(3*rpropz*cW2 - 12*rpropz*cW2*pow(w2,2) + 32*rpropz*sW2*cW2*w2 - 32*rpropz*sW2*cW2*pow(w2,2) + 2*(3-8*sW2)*(1-2*sW2)*s*mpropz*w2 - 2*(3-8*sW2)*(1-2*sW2)*s*mpropz*pow(w2,2))*cos(t0)*cos(t1)*cos(t2)*pre*ap*(1-2*w2)*pow((1+2*w2),-3)*(1-2*sW2)*b*
	pi*r*pow(cW2,-1)*LamM2*s2
	+ 32*(1 + pow(cos(t0),2))*(1 + 2*w2)*pre*ap*pow((1+2*w2),-1)*s*pi*r*LamM2*myreleftbis
	+ 32./3.*(1 + pow(cos(t0),2))*(48*sW2*pow(cW2,2) - 192*sW2*pow(cW2,2)*pow(w2,2) + mymodleft*w2 - mymodleft*pow(w2,2) + 3*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2 - 12*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2*pow(w2,2) + 9*b2*pow((1-2*sW2),2)*mpropz*s2*w2 - 9*b2*
	pow((1-2*sW2),2)*mpropz*s2*pow(w2,2))*cos(t1)*cos(t2)*pre*ap*(1-2*w2)*pow((1+2*w2),-3)*s*pi*r*pow(cW2,-1)*LamM2
	- 32*(1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*(3*rpropz*cW2 - 12*rpropz*cW2*pow(w2,2) + 16*rpropz*sW2*cW2*w2 - 16*rpropz*sW2*cW2*pow(w2,2) + (3-8*sW2)*(1-2*sW2)*s*mpropz*w2 - (3-8*sW2)*(1-2*sW2)*s*mpropz*pow(w2,2))*pre*ap*pow((1+2*w2),-2)*
	(1-2*sW2)*b*pi*r*pow(cW2,-1)*LamM2*s2
	+ 32*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*(3*g2 - 12*g2*pow(w2,2) + 16*sW2*w2 - 16*sW2*pow(w2,2))*sin(t0)*pre*ap*pow((1+2*w2),-2)*(1-2*sW2)*pow(g,-1)*b*ipropz*pi*r*LamM2*s2
	- 16*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(2*t0)*pre*ap*(1-2*w2)*pow((1+2*w2),-1)*(3-8*sW2)*(1-2*sW2)*(1-g2)*pow(g,-1)*ipropz*pi*r*LamM2*s2
	;
	expr +=  - 16*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*(3*g2 - 12*g2*pow(w2,2) + 32*sW2*w2 - 32*sW2*pow(w2,2))*sin(2*t0)*pre*ap*(1-2*w2)*pow((1+2*w2),-3)*(1-2*sW2)*pow(g,-1)*b*ipropz*pi*r*LamM2*s2
	+ 32*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(t0)*pre*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*(3-8*sW2)*(1-2*sW2)*(1-g2)*pow(g,-1)*ipropz*pi*r*LamM2*s2
	- 32./3.*(sin(t1)*sin(t2)*cos(f1 + f2) + cos(t1)*cos(t2))*(48*sW2*pow(cW2,2)*g2 - 192*sW2*pow(cW2,2)*g2*pow(w2,2) + mymodleft*w2 - mymodleft*pow(w2,2) + 3*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2*g2 - 12*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2*g2*pow(w2,2))*pow(
	sin(t0),2)*pre*ap*(1-2*w2)*pow((1+2*w2),-3)*s*pi*r*pow(cW2,-1)*LamM2*pow(g2,-1)
	+ 32./3.*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*(48*sW2*pow(cW2,2) - 192*sW2*pow(cW2,2)*pow(w2,2) + 48*sW2*pow(cW2,2)*g2 - 192*sW2*pow(cW2,2)*g2*pow(w2,2) + mymodleft*w2 - mymodleft*pow(w2,2) + 3*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2 - 12*(3-8*sW2)*
	(1-2*sW2)*s*rpropz*cW2*pow(w2,2) + 3*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2*g2 - 12*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2*g2*pow(w2,2))*sin(t0)*pre*ap*pow((1+2*w2),-2)*s*pow(g,-1)*pi*r*pow(cW2,-1)*LamM2
	+ 16*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*(3*rpropz*cW2*g2 - 12*rpropz*cW2*g2*pow(w2,2) + 16*rpropz*sW2*cW2*w2 - 16*rpropz*sW2*cW2*pow(w2,2) + (3-8*sW2)*(1-2*sW2)*s*mpropz*w2 - (3-8*sW2)*(1-2*sW2)*s*mpropz*pow(w2,2))*sin(2*t0)*pre*ap*pow(
	(1+2*w2),-2)*(1-2*sW2)*pow(g,-1)*b*pi*r*pow(cW2,-1)*LamM2*s2
	- 16./3.*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*(48*sW2*pow(cW2,2) - 192*sW2*pow(cW2,2)*pow(w2,2) + 48*sW2*pow(cW2,2)*g2 - 192*sW2*pow(cW2,2)*g2*pow(w2,2) + 2*mymodleft*w2 - 2*mymodleft*pow(w2,2) + 3*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2
	- 12*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2*pow(w2,2) + 3*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2*g2 - 12*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2*g2*pow(w2,2))*sin(2*t0)*pre*ap*(1-2*w2)*pow((1+2*w2),-3)*s*pow(g,-1)*pi*r*pow(cW2,-1)*LamM2
	- 32*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*(3*rpropz*cW2*g2 - 12*rpropz*cW2*g2*pow(w2,2) + 32*rpropz*sW2*cW2*w2 - 32*rpropz*sW2*cW2*pow(w2,2) + 2*(3-8*sW2)*(1-2*sW2)*s*mpropz*w2 - 2*(3-8*sW2)*(1-2*sW2)*s*mpropz*pow(w2,2))*sin(t0)*pre
	*ap*(1-2*w2)*pow((1+2*w2),-3)*(1-2*sW2)*pow(g,-1)*b*pi*r*pow(cW2,-1)*LamM2*s2
	- 32./3.*(cos(t1) + cos(t2))*(96*sW2*pow(cW2,2) - 384*sW2*pow(cW2,2)*pow(w2,2) + mymodleft*w2 - mymodleft*pow(w2,2) + 6*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2 - 24*(3-8*sW2)*(1-2*sW2)*s*rpropz*cW2*pow(w2,2) + 9*b2*pow((1-2*sW2),2)*mpropz*s2*w2 - 9*b2*
	pow((1-2*sW2),2)*mpropz*s2*pow(w2,2))*cos(t0)*pre*ap*pow((1+2*w2),-2)*s*pi*r*pow(cW2,-1)*LamM2
	;
}//Im(uZ)
else if (k==9){
	expr =
	- 16*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 + f2)*pre*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*b*pi*r*LamM2*s2*myreleft
	- 24*(sin(f2)*sin(t2) - sin(f1)*sin(t1))*sin(2*t0)*pre*b2*ap*(1-2*w2)*pow((1+2*w2),-1)*pow((1-2*sW2),2)*pow(sW,-1)*pow(cW,-1)*s*g*mpropz*pi*r*LamM2*s2
	- 16*(sin(f2)*sin(t2) - sin(f1)*sin(t1))*sin(t0)*pre*ap*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*g*b*pi*r*LamM2*s2*myreleft
	+ 8*(sin(f2)*sin(t2)*cos(t1) - sin(f1)*sin(t1)*cos(t2))*sin(2*t0)*pre*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*g*b*pi*r*LamM2*s2*myreleft
	+ 48*(sin(f2)*sin(t2)*cos(t1) - sin(f1)*sin(t1)*cos(t2))*sin(t0)*pre*b2*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow((1-2*sW2),2)*pow(sW,-1)*pow(cW,-1)*s*g*mpropz*pi*r*LamM2*s2
	+ 128*(sin(t1)*cos(f1) + sin(t2)*cos(f2))*sin(2*t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*g*b*ipropz*r*LamM2*pi2*s2
	- 256*(sin(t1)*cos(f1)*cos(t2) + sin(t2)*cos(f2)*cos(t1))*sin(t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*g*b*ipropz*r*LamM2*pi2*s2
	+ 256*(cos(t1) - cos(t2))*pow(sin(t0),2)*ap2*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*pow(sW,-1)*pow(cW,-1)*b*ipropz*r*LamM2*pi2*s2
	;
}//Im(uA)
else if (k==10){
	expr =
	- 32*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 + f2)*pre*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*b*pi*r*LamM2*myreleftbis
	- 48*(sin(f2)*sin(t2) - sin(f1)*sin(t1))*sin(2*t0)*pre*b2*ap*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*g*rpropz*pi*r*LamM2*s2
	- 32*(sin(f2)*sin(t2) - sin(f1)*sin(t1))*sin(t0)*pre*ap*(1-2*w2)*pow((1+2*w2),-1)*s*g*b*pi*r*LamM2*myreleftbis
	+ 16*(sin(f2)*sin(t2)*cos(t1) - sin(f1)*sin(t1)*cos(t2))*sin(2*t0)*pre*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*g*b*pi*r*LamM2*myreleftbis
	+ 96*(sin(f2)*sin(t2)*cos(t1) - sin(f1)*sin(t1)*cos(t2))*sin(t0)*pre*b2*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*g*rpropz*pi*r*LamM2*s2
	- 16*(sin(t1)*cos(f1) + sin(t2)*cos(f2))*sin(2*t0)*pre*ap*(1-2*w2)*pow((1+2*w2),-1)*(3-8*sW2)*(1-2*sW2)*g*b*ipropz*pi*r*LamM2*s2
	- 96*(sin(t1)*cos(f1) + sin(t2)*cos(f2))*sin(t0)*pre*b2*ap*(1-2*w2)*pow((1+2*w2),-1)*(1-2*sW2)*g*ipropz*pi*r*LamM2*s2
	+ 48*(sin(t1)*cos(f1)*cos(t2) + sin(t2)*cos(f2)*cos(t1))*sin(2*t0)*pre*b2*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-2*sW2)*g*ipropz*pi*r*LamM2*s2
	+ 32*(sin(t1)*cos(f1)*cos(t2) + sin(t2)*cos(f2)*cos(t1))*sin(t0)*pre*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*(3-8*sW2)*(1-2*sW2)*g*b*ipropz*pi*r*LamM2*s2
	- 32*(cos(t1) - cos(t2))*pow(sin(t0),2)*pre*ap*(1-2*w2)*pow((1+2*w2),-1)*(3-8*sW2)*(1-2*sW2)*b*ipropz*pi*r*LamM2*s2
	;
}//Re(dw)
else if (k==11){
	expr = 0;

}//Re(pqP)
else if (k==12){
	expr = 0;

}
return expr*prefac;
}
const Double_t Sright(Int_t k, Double_t sqs, const Double_t t0, const Double_t t1, const Double_t f1, const Double_t t2, const Double_t f2) {

Double_t s, s2, s4, r, r2, g2, g, b2, b, ipropz, rpropz, mpropz, prefac,
	mymodleft, mymodright, myreleft, myreleftbis, myreright,
	expr;
TComplex propz;

s = sqs*sqs; s2 = s*s; s4 = s2*s2;
r = m2/s; r2 = r*r; g2 =1./(4*r); g=sqrt(g2);
b2 = 1-4*r; b = sqrt(b2);
prefac = pb/(2*s) /(8*pi) *b /4. *1e3;

propz = 1./(s-mz*mz +i*gz*mz);
ipropz = propz.Im();
rpropz = propz.Re();
mpropz = rpropz*rpropz+ipropz*ipropz;

mymodleft = mod(16*sW2*cW2 + (3-8*sW2)*(1-2*sW2)*s*propz);
myreleft  = 16*rpropz*sW2*cW2 + (3-8*sW2)*(1-2*sW2)*s*mpropz;
myreleftbis = 16*sW2*cW2 + (3-8*sW2)*(1-2*sW2)*s*rpropz;

mymodright = mod(8*cW2 - (3-8*sW2)*propz*s);
myreright = 8*cW2*rpropz - (3-8*sW2)*s*mpropz;

//sm
if (k==0){
	expr =
	- 2./3.*(mymodright - 9*b2*mpropz*s2)*pow(sin(t0),2)*sin(t1)*sin(t2)*cos(f1 - f2)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow(cW2,-2)*pi2
	+ 32*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 - f2)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*b*ipropz*pow(cW2,-1)*pi2
	+ 2./3.*pow(sin(t0),2)*ap2*mymodright*pow(cW2,-2)*pi2*pow(g2,-1)
	+ 8*cos(t0)*cos(t1)*cos(t2)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*b*pow(cW2,-2)*pi2*myreright
	+ 8*cos(t0)*ap2*s*b*pow(cW2,-2)*pi2*myreright
	+ 2./3.*(1 + pow(cos(t0),2))*(mymodright + 9*b2*mpropz*s2)*cos(t1)*cos(t2)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow(cW2,-2)*pi2
	+ 2./3.*(1 + pow(cos(t0),2))*(mymodright + 9*b2*mpropz*s2)*ap2*pow(cW2,-2)*pi2
	+ 4*(1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*ap2*(1-2*w2)*pow((1+2*w2),-1)*s*b*pow(cW2,-2)*pi2*myreright
	- 32*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*s*pow(g,-1)*b*ipropz*pow(cW2,-1)*pi2
	- 16*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(2*t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*pow(g,-1)*b*ipropz*pow(cW2,-1)*pi2
	- 2./3.*(sin(t1)*sin(t2)*cos(f1 + f2) + cos(t1)*cos(t2))*pow(sin(t0),2)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*mymodright*pow(cW2,-2)*pi2*pow(g2,-1)
	+ 2*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(2*t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*s*pow(g,-1)*b*pow(cW2,-2)*pi2*myreright
	+ 4./3.*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*pow(g,-1)*mymodright*pow(cW2,-2)*pi2
	+ 2./3.*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(2*t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow(g,-1)*mymodright*pow(cW2,-2)*pi2
	;
	expr +=  + 4*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*pow(g,-1)*b*pow(cW2,-2)*pi2*myreright
	+ 4./3.*(cos(t1) + cos(t2))*(mymodright + 9*b2*mpropz*s2)*cos(t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*pow(cW2,-2)*pi2
	;
}//Re(lqA)
else if (k==1){
	expr = 0;

}//Re(eqA)
else if (k==2){
	expr =
	- 4*(8*cW2 - (3-8*sW2)*s*rpropz)*cos(t0)*cos(t1)*cos(t2)*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*b*pi*pow(cW2,-1)*LamM2
	- 4*(8*cW2 - (3-8*sW2)*s*rpropz)*cos(t0)*ap*s*b*pi*pow(cW2,-1)*LamM2
	- 2*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 - f2)*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*(3-8*sW2)*b*ipropz*pi*pow(cW2,-1)*LamM2*s2
	- 6*pow(sin(t0),2)*sin(t1)*sin(t2)*cos(f1 - f2)*b2*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*rpropz*pi*pow(cW2,-1)*LamM2*s2
	- 6*(1 + pow(cos(t0),2))*cos(t1)*cos(t2)*b2*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*rpropz*pi*pow(cW2,-1)*LamM2*s2
	- 6*(1 + pow(cos(t0),2))*b2*ap*rpropz*pi*pow(cW2,-1)*LamM2*s2
	- 2*(1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*(8*cW2 - (3-8*sW2)*s*rpropz)*ap*(1-2*w2)*pow((1+2*w2),-1)*s*b*pi*pow(cW2,-1)*LamM2
	+ 2*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(t0)*ap*(1-2*w2)*pow((1+2*w2),-1)*(3-8*sW2)*pow(g,-1)*b*ipropz*pi*pow(cW2,-1)*LamM2*s2
	+ (sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(2*t0)*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*(3-8*sW2)*pow(g,-1)*b*ipropz*pi*pow(cW2,-1)*LamM2*s2
	- (sin(t1)*cos(f1) - sin(t2)*cos(f2))*(8*cW2 - (3-8*sW2)*s*rpropz)*sin(2*t0)*ap*(1-2*w2)*pow((1+2*w2),-1)*s*pow(g,-1)*b*pi*pow(cW2,-1)*LamM2
	- 2*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*(8*cW2 - (3-8*sW2)*s*rpropz)*sin(t0)*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*pow(g,-1)*b*pi*pow(cW2,-1)*LamM2
	- 12*(cos(t1) + cos(t2))*cos(t0)*b2*ap*(1-2*w2)*pow((1+2*w2),-1)*rpropz*pi*pow(cW2,-1)*LamM2*s2
	;
}//Re(pqA)
else if (k==3){
	expr =
	+ 64*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 - f2)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*b*ipropz*r*pow(cW2,-1)*LamM2*pi2*s2
	+ 24*pow(sin(t0),2)*sin(t1)*sin(t2)*cos(f1 - f2)*b2*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*mpropz*r*pow(cW2,-2)*LamM2*pi2*s2
	+ 16*cos(t0)*cos(t1)*cos(t2)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*b*r*pow(cW2,-2)*LamM2*pi2*s2*myreright
	+ 16*cos(t0)*ap2*b*r*pow(cW2,-2)*LamM2*pi2*s2*myreright
	+ 24*(1 + pow(cos(t0),2))*cos(t1)*cos(t2)*b2*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*mpropz*r*pow(cW2,-2)*LamM2*pi2*s2
	+ 24*(1 + pow(cos(t0),2))*b2*ap2*s*mpropz*r*pow(cW2,-2)*LamM2*pi2*s2
	+ 8*(1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*ap2*(1-2*w2)*pow((1+2*w2),-1)*b*r*pow(cW2,-2)*LamM2*pi2*s2*myreright
	- 64*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*pow(g,-1)*b*ipropz*r*pow(cW2,-1)*LamM2*pi2*s2
	- 32*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(2*t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow(g,-1)*b*ipropz*r*pow(cW2,-1)*LamM2*pi2*s2
	+ 4*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(2*t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*pow(g,-1)*b*r*pow(cW2,-2)*LamM2*pi2*s2*myreright
	+ 8*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow(g,-1)*b*r*pow(cW2,-2)*LamM2*pi2*s2*myreright
	+ 48*(cos(t1) + cos(t2))*cos(t0)*b2*ap2*(1-2*w2)*pow((1+2*w2),-1)*s*mpropz*r*pow(cW2,-2)*LamM2*pi2*s2
	;
}//Re(lqV)
else if (k==4){
	expr = 0;

}//Re(eqV)
else if (k==5){
	expr =
	+ 2*(8*cW2 - (3-8*sW2)*s*rpropz)*pow(sin(t0),2)*sin(t1)*sin(t2)*cos(f1 - f2)*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*pi*pow(cW2,-1)*LamM2
	- 2*(8*cW2 - (3-8*sW2)*s*rpropz)*pow(sin(t0),2)*ap*s*pi*pow(cW2,-1)*LamM2*pow(g2,-1)
	- 6*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 - f2)*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*b*ipropz*pi*pow(cW2,-1)*LamM2*s2
	- 12*cos(t0)*cos(t1)*cos(t2)*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*b*rpropz*pi*pow(cW2,-1)*LamM2*s2
	- 12*cos(t0)*ap*b*rpropz*pi*pow(cW2,-1)*LamM2*s2
	- 2*(1 + pow(cos(t0),2))*(8*cW2 - (3-8*sW2)*s*rpropz)*cos(t1)*cos(t2)*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*pi*pow(cW2,-1)*LamM2
	- 2*(1 + pow(cos(t0),2))*(8*cW2 - (3-8*sW2)*s*rpropz)*ap*s*pi*pow(cW2,-1)*LamM2
	- 6*(1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*ap*(1-2*w2)*pow((1+2*w2),-1)*b*rpropz*pi*pow(cW2,-1)*LamM2*s2
	+ 6*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(t0)*ap*(1-2*w2)*pow((1+2*w2),-1)*pow(g,-1)*b*ipropz*pi*pow(cW2,-1)*LamM2*s2
	+ 3*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(2*t0)*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow(g,-1)*b*ipropz*pi*pow(cW2,-1)*LamM2*s2
	+ 2*(sin(t1)*sin(t2)*cos(f1 + f2) + cos(t1)*cos(t2))*(8*cW2 - (3-8*sW2)*s*rpropz)*pow(sin(t0),2)*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*pi*pow(cW2,-1)*LamM2*pow(g2,-1)
	- 4*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*(8*cW2 - (3-8*sW2)*s*rpropz)*sin(t0)*ap*(1-2*w2)*pow((1+2*w2),-1)*s*pow(g,-1)*pi*pow(cW2,-1)*LamM2
	- 3*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(2*t0)*ap*(1-2*w2)*pow((1+2*w2),-1)*pow(g,-1)*b*rpropz*pi*pow(cW2,-1)*LamM2*s2
	- 2*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*(8*cW2 - (3-8*sW2)*s*rpropz)*sin(2*t0)*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*pow(g,-1)*pi*pow(cW2,-1)*LamM2
	;
	expr +=  - 6*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(t0)*ap*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow(g,-1)*b*rpropz*pi*pow(cW2,-1)*LamM2*s2
	- 4*(cos(t1) + cos(t2))*(8*cW2 - (3-8*sW2)*s*rpropz)*cos(t0)*ap*(1-2*w2)*pow((1+2*w2),-1)*s*pi*pow(cW2,-1)*LamM2
	;
}//Re(pqV)
else if (k==6){
	expr =
	- 8*pow(sin(t0),2)*sin(t1)*sin(t2)*cos(f1 - f2)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*r*pow(cW2,-2)*LamM2*pi2*s2*myreright
	+ 8*pow(sin(t0),2)*ap2*r*pow(cW2,-2)*LamM2*pi2*s2*pow(g2,-1)*myreright
	+ 48*cos(t0)*cos(t1)*cos(t2)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*b*mpropz*r*pow(cW2,-2)*LamM2*pi2*s2
	+ 48*cos(t0)*ap2*s*b*mpropz*r*pow(cW2,-2)*LamM2*pi2*s2
	+ 8*(1 + pow(cos(t0),2))*cos(t1)*cos(t2)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*r*pow(cW2,-2)*LamM2*pi2*s2*myreright
	+ 8*(1 + pow(cos(t0),2))*ap2*r*pow(cW2,-2)*LamM2*pi2*s2*myreright
	+ 24*(1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*ap2*(1-2*w2)*pow((1+2*w2),-1)*s*b*mpropz*r*pow(cW2,-2)*LamM2*pi2*s2
	- 8*(sin(t1)*sin(t2)*cos(f1 + f2) + cos(t1)*cos(t2))*pow(sin(t0),2)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*r*pow(cW2,-2)*LamM2*pi2*s2*pow(g2,-1)*myreright
	+ 12*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(2*t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*s*pow(g,-1)*b*mpropz*r*pow(cW2,-2)*LamM2*pi2*s2
	+ 16*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*sin(t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*pow(g,-1)*r*pow(cW2,-2)*LamM2*pi2*s2*myreright
	+ 8*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(2*t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*pow(g,-1)*r*pow(cW2,-2)*LamM2*pi2*s2*myreright
	+ 24*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*sin(t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*pow(g,-1)*b*mpropz*r*pow(cW2,-2)*LamM2*pi2*s2
	+ 16*(cos(t1) + cos(t2))*cos(t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*r*pow(cW2,-2)*LamM2*pi2*s2*myreright
	;
}//Re(uZ)
else if (k==7){
	expr =
	- 64*(1 + 2*w2)*pow(sin(t0),2)*ap2*pow((1+2*w2),-1)*sW*pow(cW,-1)*r*pow(cW2,-1)*LamM2*pi2*s2*myreright
	- 384*(1 + 2*w2)*cos(t0)*ap2*pow((1+2*w2),-1)*sW*pow(cW,-1)*s*b*mpropz*r*pow(cW2,-1)*LamM2*pi2*s2
	- 128*( - 32*rpropz*cW2*w2 + 32*rpropz*cW2*pow(w2,2) + 3*s*mpropz - 12*s*mpropz*pow(w2,2) + 4*(3-8*sW2)*s*mpropz*w2 - 4*(3-8*sW2)*s*mpropz*pow(w2,2))*cos(t0)*cos(t1)*cos(t2)*ap2*(1-2*w2)*pow((1+2*w2),-3)*sW*pow(cW,-1)*b*r*pow(cW2,-1)*LamM2*pi2*s2
	- 64./3.*(2*mymodright*w2 - 2*mymodright*pow(w2,2) - 24*s*rpropz*cW2 + 96*s*rpropz*cW2*pow(w2,2) + 3*(3-8*sW2)*mpropz*s2 - 12*(3-8*sW2)*mpropz*s2*pow(w2,2) - 18*b2*mpropz*s2*w2 + 18*b2*mpropz*s2*pow(w2,2))*pow(sin(t0),2)*sin(t1)*sin(t2)*cos(f1 - 
	f2)*ap2*(1-2*w2)*pow((1+2*w2),-3)*sW*pow(cW,-1)*s*r*pow(cW2,-1)*LamM2*pi2
	+ 2048*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 - f2)*ap2*(1-w2)*(1-2*w2)*pow((1+2*w2),-3)*sW*pow(cW,-1)*b*ipropz*r*LamM2*pi2*s2*w2
	- 64*(1 + pow(cos(t0),2))*(1 + 2*w2)*ap2*pow((1+2*w2),-1)*sW*pow(cW,-1)*r*pow(cW2,-1)*LamM2*pi2*s2*myreright
	+ 64./3.*(1 + pow(cos(t0),2))*(2*mymodright*w2 - 2*mymodright*pow(w2,2) - 24*s*rpropz*cW2 + 96*s*rpropz*cW2*pow(w2,2) + 3*(3-8*sW2)*mpropz*s2 - 12*(3-8*sW2)*mpropz*s2*pow(w2,2) + 18*b2*mpropz*s2*w2 - 18*b2*mpropz*s2*pow(w2,2))*cos(t1)*cos(t2)*ap2*
	(1-2*w2)*pow((1+2*w2),-3)*sW*pow(cW,-1)*s*r*pow(cW2,-1)*LamM2*pi2
	- 64*(1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*( - 16*rpropz*cW2*w2 + 16*rpropz*cW2*pow(w2,2) + 3*s*mpropz - 12*s*mpropz*pow(w2,2) + 2*(3-8*sW2)*s*mpropz*w2 - 2*(3-8*sW2)*s*mpropz*pow(w2,2))*ap2*pow((1+2*w2),-2)*sW*pow(cW,-1)*b*r*pow(cW2,-1)*LamM2*
	pi2*s2
	+ 256*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(2*t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*(1-g2)*sW*pow(cW,-1)*pow(g,-1)*ipropz*r*LamM2*pi2*s2
	- 1024*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(t0)*ap2*(1-w2)*pow((1+2*w2),-2)*sW*pow(cW,-1)*pow(g,-1)*b*ipropz*r*LamM2*pi2*s2*w2
	- 1024*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(2*t0)*ap2*(1-w2)*(1-2*w2)*pow((1+2*w2),-3)*sW*pow(cW,-1)*pow(g,-1)*b*ipropz*r*LamM2*pi2*s2*w2
	;
	expr +=  + 512*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*(1-g2)*sW*pow(cW,-1)*pow(g,-1)*ipropz*r*LamM2*pi2*s2
	- 64./3.*(sin(t1)*sin(t2)*cos(f1 + f2) + cos(t1)*cos(t2))*(2*mymodright*w2 - 2*mymodright*pow(w2,2) - 24*s*rpropz*cW2*g2 + 96*s*rpropz*cW2*g2*pow(w2,2) + 3*(3-8*sW2)*mpropz*s2*g2 - 12*(3-8*sW2)*mpropz*s2*g2*pow(w2,2))*pow(sin(t0),2)*ap2*(1-2*w2)*
	pow((1+2*w2),-3)*sW*pow(cW,-1)*s*r*pow(cW2,-1)*LamM2*pi2*pow(g2,-1)
	- 32*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*( - 16*rpropz*cW2*w2 + 16*rpropz*cW2*pow(w2,2) + 3*s*mpropz*g2 - 12*s*mpropz*g2*pow(w2,2) + 2*(3-8*sW2)*s*mpropz*w2 - 2*(3-8*sW2)*s*mpropz*pow(w2,2))*sin(2*t0)*ap2*pow((1+2*w2),-2)*sW*pow(cW,-1)*pow(g,-1)*b
	*r*pow(cW2,-1)*LamM2*pi2*s2
	+ 64./3.*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*(2*mymodright*w2 - 2*mymodright*pow(w2,2) - 24*s*rpropz*cW2 + 96*s*rpropz*cW2*pow(w2,2) - 24*s*rpropz*cW2*g2 + 96*s*rpropz*cW2*g2*pow(w2,2) + 3*(3-8*sW2)*mpropz*s2 - 12*(3-8*sW2)*mpropz*s2*pow(w2,2) + 3
	*(3-8*sW2)*mpropz*s2*g2 - 12*(3-8*sW2)*mpropz*s2*g2*pow(w2,2))*sin(t0)*ap2*pow((1+2*w2),-2)*sW*pow(cW,-1)*s*pow(g,-1)*r*pow(cW2,-1)*LamM2*pi2
	- 64*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*( - 32*rpropz*cW2*w2 + 32*rpropz*cW2*pow(w2,2) + 3*s*mpropz*g2 - 12*s*mpropz*g2*pow(w2,2) + 4*(3-8*sW2)*s*mpropz*w2 - 4*(3-8*sW2)*s*mpropz*pow(w2,2))*sin(t0)*ap2*(1-2*w2)*pow((1+2*w2),-3)*sW
	*pow(cW,-1)*pow(g,-1)*b*r*pow(cW2,-1)*LamM2*pi2*s2
	+ 32./3.*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*(4*mymodright*w2 - 4*mymodright*pow(w2,2) - 24*s*rpropz*cW2 + 96*s*rpropz*cW2*pow(w2,2) - 24*s*rpropz*cW2*g2 + 96*s*rpropz*cW2*g2*pow(w2,2) + 3*(3-8*sW2)*mpropz*s2 - 12*(3-8*sW2)*mpropz*
	s2*pow(w2,2) + 3*(3-8*sW2)*mpropz*s2*g2 - 12*(3-8*sW2)*mpropz*s2*g2*pow(w2,2))*sin(2*t0)*ap2*(1-2*w2)*pow((1+2*w2),-3)*sW*pow(cW,-1)*s*pow(g,-1)*r*pow(cW2,-1)*LamM2*pi2
	+ 128./3.*(cos(t1) + cos(t2))*(mymodright*w2 - mymodright*pow(w2,2) - 24*s*rpropz*cW2 + 96*s*rpropz*cW2*pow(w2,2) + 3*(3-8*sW2)*mpropz*s2 - 12*(3-8*sW2)*mpropz*s2*pow(w2,2) + 9*b2*mpropz*s2*w2 - 9*b2*mpropz*s2*pow(w2,2))*cos(t0)*ap2*pow(
	(1+2*w2),-2)*sW*pow(cW,-1)*s*r*pow(cW2,-1)*LamM2*pi2
	;
}//Re(uA)
else if (k==8){
	expr =
	+ 64*(1 + 2*w2)*(8*cW2 - (3-8*sW2)*s*rpropz)*pow(sin(t0),2)*ap2*pow((1+2*w2),-1)*s*r*pow(cW2,-1)*LamM2*pi2
	+ 384*(1 + 2*w2)*cos(t0)*ap2*pow((1+2*w2),-1)*b*rpropz*r*pow(cW2,-1)*LamM2*pi2*s2
	+ 64*(3 - 12*pow(w2,2) + 32*sW2*w2 - 32*sW2*pow(w2,2))*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 - f2)*ap2*(1-2*w2)*pow((1+2*w2),-3)*b*ipropz*r*pow(cW2,-1)*LamM2*pi2*s2
	- 64./3.*(24*pow(cW2,2) - 96*pow(cW2,2)*pow(w2,2) + 2*mymodright*sW2*w2 - 2*mymodright*sW2*pow(w2,2) - 3*(3-8*sW2)*s*rpropz*cW2 + 12*(3-8*sW2)*s*rpropz*cW2*pow(w2,2) - 18*b2*mpropz*sW2*s2*w2 + 18*b2*mpropz*sW2*s2*pow(w2,2))*pow(sin(t0),2)*sin(t1)*
	sin(t2)*cos(f1 - f2)*ap2*(1-2*w2)*pow((1+2*w2),-3)*s*r*pow(cW2,-2)*LamM2*pi2
	+ 128*(3*rpropz*cW2 - 12*rpropz*cW2*pow(w2,2) + 32*rpropz*sW2*cW2*w2 - 32*rpropz*sW2*cW2*pow(w2,2) - 4*(3-8*sW2)*s*mpropz*sW2*w2 + 4*(3-8*sW2)*s*mpropz*sW2*pow(w2,2))*cos(t0)*cos(t1)*cos(t2)*ap2*(1-2*w2)*pow((1+2*w2),-3)*b*r*pow(cW2,-2)*LamM2*pi2*
	s2
	+ 64*(1 + pow(cos(t0),2))*(1 + 2*w2)*(8*cW2 - (3-8*sW2)*s*rpropz)*ap2*pow((1+2*w2),-1)*s*r*pow(cW2,-1)*LamM2*pi2
	+ 64./3.*(1 + pow(cos(t0),2))*(24*pow(cW2,2) - 96*pow(cW2,2)*pow(w2,2) + 2*mymodright*sW2*w2 - 2*mymodright*sW2*pow(w2,2) - 3*(3-8*sW2)*s*rpropz*cW2 + 12*(3-8*sW2)*s*rpropz*cW2*pow(w2,2) + 18*b2*mpropz*sW2*s2*w2 - 18*b2*mpropz*sW2*s2*pow(w2,2))*
	cos(t1)*cos(t2)*ap2*(1-2*w2)*pow((1+2*w2),-3)*s*r*pow(cW2,-2)*LamM2*pi2
	+ 64*(1 + pow(cos(t0),2))*(cos(t1) + cos(t2))*(3*rpropz*cW2 - 12*rpropz*cW2*pow(w2,2) + 16*rpropz*sW2*cW2*w2 - 16*rpropz*sW2*cW2*pow(w2,2) - 2*(3-8*sW2)*s*mpropz*sW2*w2 + 2*(3-8*sW2)*s*mpropz*sW2*pow(w2,2))*ap2*pow((1+2*w2),-2)*b*r*pow(cW2,-2)*
	LamM2*pi2*s2
	- 64*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*(3*g2 - 12*g2*pow(w2,2) + 16*sW2*w2 - 16*sW2*pow(w2,2))*sin(t0)*ap2*pow((1+2*w2),-2)*pow(g,-1)*b*ipropz*r*pow(cW2,-1)*LamM2*pi2*s2
	- 32*(sin(f2)*sin(t2) + sin(f1)*sin(t1))*sin(2*t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*(3-8*sW2)*(1-g2)*pow(g,-1)*ipropz*r*pow(cW2,-1)*LamM2*pi2*s2
	;
	expr +=  - 32*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*(3*g2 - 12*g2*pow(w2,2) + 32*sW2*w2 - 32*sW2*pow(w2,2))*sin(2*t0)*ap2*(1-2*w2)*pow((1+2*w2),-3)*pow(g,-1)*b*ipropz*r*pow(cW2,-1)*LamM2*pi2*s2
	- 64*(sin(f2)*sin(t2)*cos(t1) + sin(f1)*sin(t1)*cos(t2))*sin(t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*(3-8*sW2)*(1-g2)*pow(g,-1)*ipropz*r*pow(cW2,-1)*LamM2*pi2*s2
	- 64./3.*(sin(t1)*sin(t2)*cos(f1 + f2) + cos(t1)*cos(t2))*(24*pow(cW2,2)*g2 - 96*pow(cW2,2)*g2*pow(w2,2) + 2*mymodright*sW2*w2 - 2*mymodright*sW2*pow(w2,2) - 3*(3-8*sW2)*s*rpropz*cW2*g2 + 12*(3-8*sW2)*s*rpropz*cW2*g2*pow(w2,2))*pow(sin(t0),2)*ap2*
	(1-2*w2)*pow((1+2*w2),-3)*s*r*pow(cW2,-2)*LamM2*pi2*pow(g2,-1)
	+ 64./3.*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*(24*pow(cW2,2) - 96*pow(cW2,2)*pow(w2,2) + 24*pow(cW2,2)*g2 - 96*pow(cW2,2)*g2*pow(w2,2) + 2*mymodright*sW2*w2 - 2*mymodright*sW2*pow(w2,2) - 3*(3-8*sW2)*s*rpropz*cW2 + 12*(3-8*sW2)*s*rpropz*cW2*pow(
	w2,2) - 3*(3-8*sW2)*s*rpropz*cW2*g2 + 12*(3-8*sW2)*s*rpropz*cW2*g2*pow(w2,2))*sin(t0)*ap2*pow((1+2*w2),-2)*s*pow(g,-1)*r*pow(cW2,-2)*LamM2*pi2
	+ 32*(sin(t1)*cos(f1) - sin(t2)*cos(f2))*(3*rpropz*cW2*g2 - 12*rpropz*cW2*g2*pow(w2,2) + 16*rpropz*sW2*cW2*w2 - 16*rpropz*sW2*cW2*pow(w2,2) - 2*(3-8*sW2)*s*mpropz*sW2*w2 + 2*(3-8*sW2)*s*mpropz*sW2*pow(w2,2))*sin(2*t0)*ap2*pow((1+2*w2),-2)*pow(
	g,-1)*b*r*pow(cW2,-2)*LamM2*pi2*s2
	+ 32./3.*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*(24*pow(cW2,2) - 96*pow(cW2,2)*pow(w2,2) + 24*pow(cW2,2)*g2 - 96*pow(cW2,2)*g2*pow(w2,2) + 4*mymodright*sW2*w2 - 4*mymodright*sW2*pow(w2,2) - 3*(3-8*sW2)*s*rpropz*cW2 + 12*(3-8*sW2)*s*
	rpropz*cW2*pow(w2,2) - 3*(3-8*sW2)*s*rpropz*cW2*g2 + 12*(3-8*sW2)*s*rpropz*cW2*g2*pow(w2,2))*sin(2*t0)*ap2*(1-2*w2)*pow((1+2*w2),-3)*s*pow(g,-1)*r*pow(cW2,-2)*LamM2*pi2
	+ 64*(sin(t1)*cos(f1)*cos(t2) - sin(t2)*cos(f2)*cos(t1))*(3*rpropz*cW2*g2 - 12*rpropz*cW2*g2*pow(w2,2) + 32*rpropz*sW2*cW2*w2 - 32*rpropz*sW2*cW2*pow(w2,2) - 4*(3-8*sW2)*s*mpropz*sW2*w2 + 4*(3-8*sW2)*s*mpropz*sW2*pow(w2,2))*sin(t0)*ap2*(1-2*w2)*
	pow((1+2*w2),-3)*pow(g,-1)*b*r*pow(cW2,-2)*LamM2*pi2*s2
	+ 128./3.*(cos(t1) + cos(t2))*(24*pow(cW2,2) - 96*pow(cW2,2)*pow(w2,2) + mymodright*sW2*w2 - mymodright*sW2*pow(w2,2) - 3*(3-8*sW2)*s*rpropz*cW2 + 12*(3-8*sW2)*s*rpropz*cW2*pow(w2,2) + 9*b2*mpropz*sW2*s2*w2 - 9*b2*mpropz*sW2*s2*pow(w2,2))*cos(t0)*
	ap2*pow((1+2*w2),-2)*s*r*pow(cW2,-2)*LamM2*pi2
	;
}//Im(uZ)
else if (k==9){
	expr =
	+ 64*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 + f2)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*sW*pow(cW,-1)*b*r*pow(cW2,-1)*LamM2*pi2*s2*myreright
	+ 96*(sin(f2)*sin(t2) - sin(f1)*sin(t1))*sin(2*t0)*b2*ap2*(1-2*w2)*pow((1+2*w2),-1)*sW*pow(cW,-1)*s*g*mpropz*r*pow(cW2,-1)*LamM2*pi2*s2
	+ 64*(sin(f2)*sin(t2) - sin(f1)*sin(t1))*sin(t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*sW*pow(cW,-1)*g*b*r*pow(cW2,-1)*LamM2*pi2*s2*myreright
	+ 32*(sin(f2)*sin(t2)*cos(t1) - sin(f1)*sin(t1)*cos(t2))*sin(2*t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*sW*pow(cW,-1)*g*b*r*pow(cW2,-1)*LamM2*pi2*s2*myreright
	+ 192*(sin(f2)*sin(t2)*cos(t1) - sin(f1)*sin(t1)*cos(t2))*sin(t0)*b2*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*sW*pow(cW,-1)*s*g*mpropz*r*pow(cW2,-1)*LamM2*pi2*s2
	+ 256*(sin(t1)*cos(f1) + sin(t2)*cos(f2))*sin(2*t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*sW*pow(cW,-1)*g*b*ipropz*r*LamM2*pi2*s2
	+ 512*(sin(t1)*cos(f1)*cos(t2) + sin(t2)*cos(f2)*cos(t1))*sin(t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*sW*pow(cW,-1)*g*b*ipropz*r*LamM2*pi2*s2
	- 512*(cos(t1) - cos(t2))*pow(sin(t0),2)*ap2*(1-2*w2)*pow((1+2*w2),-1)*sW*pow(cW,-1)*b*ipropz*r*LamM2*pi2*s2
	;
}//Im(uA)
else if (k==10){
	expr =
	- 64*(8*cW2 - (3-8*sW2)*s*rpropz)*pow(sin(t0),2)*sin(t1)*sin(t2)*sin(f1 + f2)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*b*r*pow(cW2,-1)*LamM2*pi2
	- 64*(sin(f2)*sin(t2) - sin(f1)*sin(t1))*(8*cW2 - (3-8*sW2)*s*rpropz)*sin(t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*s*g*b*r*pow(cW2,-1)*LamM2*pi2
	- 96*(sin(f2)*sin(t2) - sin(f1)*sin(t1))*sin(2*t0)*b2*ap2*(1-2*w2)*pow((1+2*w2),-1)*g*rpropz*r*pow(cW2,-1)*LamM2*pi2*s2
	- 32*(sin(f2)*sin(t2)*cos(t1) - sin(f1)*sin(t1)*cos(t2))*(8*cW2 - (3-8*sW2)*s*rpropz)*sin(2*t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*s*g*b*r*pow(cW2,-1)*LamM2*pi2
	- 192*(sin(f2)*sin(t2)*cos(t1) - sin(f1)*sin(t1)*cos(t2))*sin(t0)*b2*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*g*rpropz*r*pow(cW2,-1)*LamM2*pi2*s2
	- 32*(sin(t1)*cos(f1) + sin(t2)*cos(f2))*sin(2*t0)*ap2*(1-2*w2)*pow((1+2*w2),-1)*(3-8*sW2)*g*b*ipropz*r*pow(cW2,-1)*LamM2*pi2*s2
	+ 192*(sin(t1)*cos(f1) + sin(t2)*cos(f2))*sin(t0)*b2*ap2*(1-2*w2)*pow((1+2*w2),-1)*g*ipropz*r*pow(cW2,-1)*LamM2*pi2*s2
	+ 96*(sin(t1)*cos(f1)*cos(t2) + sin(t2)*cos(f2)*cos(t1))*sin(2*t0)*b2*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*g*ipropz*r*pow(cW2,-1)*LamM2*pi2*s2
	- 64*(sin(t1)*cos(f1)*cos(t2) + sin(t2)*cos(f2)*cos(t1))*sin(t0)*ap2*pow((1-2*w2),2)*pow((1+2*w2),-2)*(3-8*sW2)*g*b*ipropz*r*pow(cW2,-1)*LamM2*pi2*s2
	+ 64*(cos(t1) - cos(t2))*pow(sin(t0),2)*ap2*(1-2*w2)*pow((1+2*w2),-1)*(3-8*sW2)*b*ipropz*r*pow(cW2,-1)*LamM2*pi2*s2
	;
}//Re(dw)
else if (k==11){
	expr = 0;

}//Re(pqP)
else if (k==12){
	expr = 0;

}
return expr*prefac;
}


class LinearIntegrand {
	public:
	Double_t sqs, pp, pe, k, l;
	
	LinearIntegrand( Double_t input_sqs,
		Double_t input_pp, Double_t input_pe,
		Int_t input_k,
		Int_t input_l){
		
		sqs = input_sqs;
		
		pp = input_pp;
		pe = input_pe;

		k = input_k;
		l = input_l;

	}

	Double_t operator()(const Double_t *x) const{
		Double_t t0,t1,f1,t2,f2,expr;

		t0 = x[0];
		t1 = x[1];
		f1 = x[2];
		t2 = x[3];
		f2 = x[4];

		if( k==0 ){
			expr    = (1+pp)*(1-pe)*Sleft (l,sqs,t0,t1,f1,t2,f2)
				+ (1-pp)*(1+pe)*Sright(l,sqs,t0,t1,f1,t2,f2);
		} else {
			expr    =((1+pp)*(1-pe)*Sleft (k,sqs,t0,t1,f1,t2,f2)
				 +(1-pp)*(1+pe)*Sright(k,sqs,t0,t1,f1,t2,f2))
				*((1+pp)*(1-pe)*Sleft (l,sqs,t0,t1,f1,t2,f2)
				 +(1-pp)*(1+pe)*Sright(l,sqs,t0,t1,f1,t2,f2))
				/((1+pp)*(1-pe)*Sleft (0,sqs,t0,t1,f1,t2,f2)
				 +(1-pp)*(1+pe)*Sright(0,sqs,t0,t1,f1,t2,f2))
				;
		}

	return expr*sin(t0)*sin(t1)*sin(t2)/(32.*pi2);
	}
	Double_t operator()(const Double_t t0, const Double_t t1, const Double_t f1, const Double_t t2, const Double_t f2){
	const Double_t x[5] = {t0,t1,f1,t2,f2};
	return operator()(x);
	}
};


array<Double_t,2> integral(LinearIntegrand integrand){
    const Double_t ini[5] = {0.,0.,-pi,0.,-pi};
    const Double_t fin[5] = {pi,pi, pi,pi, pi};
    array<Double_t,2> res;
    
    // kDEFAULT = -1, kADAPTIVE, kVEGAS, kMISER, kPLAIN
    ROOT::Math::Functor func( integrand, 5 );
    ROOT::Math::IntegratorMultiDim inte( func ); //, ROOT::Math::IntegrationMultiDim::kADAPTIVE );
    res[0] = inte.Integral(ini,fin);
    res[1] = inte.Error();
    return res;
}


