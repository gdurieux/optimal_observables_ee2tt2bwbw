// compile with g++ optimal_observables_ee2tt2bwbw_from_lhe.cpp `root-config --libs --cflags`
#include "string.h"
#include <iostream>
#include <fstream>

#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include <stdio.h>
#include <complex.h>

#define TAM 1000
#define PI 3.1415
#define Nobs 6

using namespace std;

const float sqs = 500.;
const float polp = 0.;
const float pole = 0.;
const string lhe_name = "unweighted_events.lhe";
const string root_name = "output.root";

float mod(float _Complex z){
    float r,i;
    r = __real__ z;
    i = __imag__ z;
    return r*r+i*i;
};

const float pi = 3.141592653589793;
const float Lam = 1000.;
const float m = 172.5;
const float gz = 2.4952  ;
const float gw = 2.085  ;
const float mw = 79.82436 ;
const float mz = 91.1876 ;
const float gf = 1.16637e-5 ;
const float sW2 = 1-pow(mw/mz,2) ;
const float sW= sqrt(sW2);
const float cW2 = 1-sW2;
const float cW = sqrt(cW2);
const float apM1 = pi/(sqrt(2)*gf*mw*mw*sW2) ;
const float ap =1./apM1  ;
const float aps = 0.1076705;
const float pb =pow(197.3269631e-18,2)/1e-40 ;
const float LamM2 = 1./Lam/Lam;
const float m2 = m*m;
const float w = mw/m;
const float w2 = w*w;
const float ap2 = ap*ap;
const float pi2 = pi*pi;
const float _Complex i = {0,1};
const float s = sqs*sqs;
const float s2 = s*s;
const float s4 = s2*s2;
const float r = m2/s;
const float r2 = r*r;
const float g2 =1./(4*r);
const float g=sqrt(g2);
const float b2 = 1-4*r;
const float b = sqrt(b2);
const float pre = ap*pi/sW2/cW2;
const float pre2 = pre*pre;
const float _Complex propz = 1./(s-mz*mz +i*gz*mz);
const float ipropz = __imag__ propz;
const float rpropz = __real__ propz;
const float mpropz = rpropz*rpropz+ipropz*ipropz;

/*
	Compute the cosine and sine of the azimuthal helicity angles φW⁺,
	as well as the cosine of the polar helicity angle θW⁺
	
	cos(θt ) = (e⁺∙t) / sqrt{ (t∙t) (e⁺∙e⁺) }
	cos(θW⁺) = (W⁺∙t) / sqrt{ (t∙t) (W⁺∙W⁺) }
	cos(θW⁻) = (W⁻∙t̅) / sqrt{ (t̅∙t̅) (W⁻∙W⁻) }
	cos(φW⁺) = [ (e⁻∙t) (W⁺∙t) - (e⁻∙W⁺) (t∙t) ] /sqrt{ (t∙t) (e⁻∙e⁻) - (t∙e⁻)² } /sqrt{ (t∙t) (W⁺∙W⁺) - (t∙W⁺)² }
	cos(φW⁻) = [ (e⁺∙t̅) (W⁻∙t̅) - (e⁺∙W⁻) (t̅∙t̅) ] /sqrt{ (t̅∙t̅) (e⁺∙e⁺) - (t̅∙e⁺)² } /sqrt{ (t̅∙t̅) (W⁻∙W⁻) - (t̅∙W⁻)² }
	sin(φW⁺) =        e⁻∙(t⨯W⁺) sqrt{t∙t}       /sqrt{ (t∙t) (e⁻∙e⁻) - (t∙e⁻)² } /sqrt{ (t∙t) (W⁺∙W⁺) - (t∙W⁺)² }
	sin(φW⁻) =        e⁺∙(t̅⨯W⁻) sqrt{t̅∙t̅}       /sqrt{ (t̅∙t̅) (e⁺∙e⁺) - (t̅∙e⁺)² } /sqrt{ (t̅∙t̅) (W⁻∙W⁻) - (t̅∙W⁻)² }
	.    e⁻ is the three-momentum of the beam electron in the centre-of-mass frame,
	.    e⁺ is the three-momentum of the beam positron in the centre-of-mass frame,
	.    t  is the three-momentum of the     top       in the centre-of-mass frame,
	.    t̅  is the three-momentum of the antitop       in the centre-of-mass frame,
	.    W⁺ is the three-momentum of the W⁺ in the restframe of the     top
	.    W⁻ is the three-momentum of the W⁻ in the restframe of the antitop
*/
std::vector<float> angle_calc(TLorentzVector top, TLorentzVector bosonW, TLorentzVector direction){

	float tmp_den, tmp_num;
	std::vector<float> hel(3);

	TLorentzVector boostW(bosonW);
	boostW.Boost(-(top).BoostVector());    // brings W boson to the top restframe.
	TVector3 dirVect    = direction.Vect();// e⁻
	TVector3 topVect    = top.Vect();      // t
	TVector3 boostWVect = boostW.Vect();   // W⁺

	// denominator common to cos(φW⁺) and sin(φW⁺)
	tmp_den = sqrt(topVect.Mag2()*   dirVect.Mag2() - pow(topVect.Dot(   dirVect),2))*
	          sqrt(topVect.Mag2()*boostWVect.Mag2() - pow(topVect.Dot(boostWVect),2));
	// numerator of cos(φW⁺)
	tmp_num = (dirVect.Dot(   topVect))*(boostWVect.Dot(topVect)) 
	        - (dirVect.Dot(boostWVect))*topVect.Mag2();
	// the cosine of the azimuthal helicity angle cos(φW⁺)
	hel[0] = tmp_num/tmp_den;
	
	// the numerator of sin(φW⁺)
	tmp_num = dirVect.Dot(topVect.Cross(boostWVect))*topVect.Mag();
	// the sine of the azimuthal helicity angle sin(φW⁺)
	hel[1] = tmp_num/tmp_den;
	
	// the cosine of the polar helicity angle cos(θW⁺)
	tmp_num = topVect.Dot(boostWVect);
	tmp_den = (topVect.Mag())*(boostWVect.Mag());
	hel[2] = tmp_num/tmp_den;

	return hel;
}


/*
	Linear EFT terms for a right-handed beam	
*/
float Sright(int j, float t0, float t1, float f1, float t2, float f2){
	
	f1 = f1 + pi;
	f2 = f2 + pi;
	
	float expr;
	float prefac = pb/(2*s) /(8*pi) *b /4. *1e3;


	float mymodleft, mymodright, myreleft, myreleftbis, myreright;

	mymodleft = mod(16*sW2*cW2 + (3-8*sW2)*(1-2*sW2)*s*propz);
	myreleft  = 16*rpropz*sW2*cW2 + (3-8*sW2)*(1-2*sW2)*s*mpropz;
	myreleftbis = 16*sW2*cW2 + (3-8*sW2)*(1-2*sW2)*s*rpropz;
	mymodright = mod(8*cW2 - (3-8*sW2)*propz*s);
	myreright = 8*cW2*rpropz - (3-8*sW2)*s*mpropz;
    
        switch( j ){
                //S['sm']
            case 0:
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
                break;
                
                //S['Re(lqA)']
            case 1:
                expr = 0;
                break;
                
                //S['Re(eqA)']
            case 2:
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
                break;
                
                //S['Re(pqA)']
            case 3:
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
                break;
                
                //S['Re(lqV)']
            case 4:
                expr = 0;
                break;
                
                //S['Re(eqV)']
            case 5:
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
                break;
                
                //S['Re(pqV)']
            case 6:
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
                break;
                
                //S['Re(uZ)']
            case 7:
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
                break;
                
                //S['Re(uA)']
            case 8:
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
                break;
                
                //S['Im(uZ)']
            case 9:
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
                break;
                
                //S['Im(uA)']
            case 10:
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
                break;
                
                //S['Re(dw)']
            case 11:
                expr = 0;
                break;
                
            case 12:
                expr = 0;
                break;
                
            default:
                std::cout << "No valid integer" << std::endl;
                break;
        }
        return expr*prefac;
}
    
/*
	Linear EFT terms for a left-handed beam	
*/
float Sleft(Int_t j, float t0, float t1, float f1, float t2, float f2){
	float expr;
	float prefac = pb/(2*s) /(8*pi) *b /4. *1e3;

	float mymodleft, mymodright, myreleft, myreleftbis, myreright;
	mymodleft = mod(16*sW2*cW2 + (3-8*sW2)*(1-2*sW2)*s*propz);
	myreleft  = 16*rpropz*sW2*cW2 + (3-8*sW2)*(1-2*sW2)*s*mpropz;
	myreleftbis = 16*sW2*cW2 + (3-8*sW2)*(1-2*sW2)*s*rpropz;
	mymodright = mod(8*cW2 - (3-8*sW2)*propz*s);
	myreright = 8*cW2*rpropz - (3-8*sW2)*s*mpropz;
    
        switch( j ){
            case 0:
                //	S['sm']
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
                break;
                
            case 1:
                //	S['Re(lqA)']
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
                break;
                
            case 2:
                //	S['Re(eqA)']
                expr = 0;
                break;
                
            case 3:
                //	S['Re(pqA)']
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
                break;
                
            case 4:
                //	S['Re(lqV)']
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
                break;
                
            case 5:
                //	S['Re(eqV)']
                expr = 0;
                break;
                
            case 6:
                //	S['Re(pqV)']
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
                break;
                
            case 7:
                //	S['Re(uZ)']
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
                break;
                
            case 8:
                //	S['Re(uA)']
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
                break;
                
            case 9:
                //	S['Im(uZ)']
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
                break;
                
            case 10:
                //	S['Im(uA)']
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
                break;
                
                //S['Re(dw)']
            case 11:
                expr = 0;
                break;
                
            case 12:
                expr = 0;
                break;
                
            default:
                std::cout << "No valid integer" << std::endl;
                break;
                
        }
        return expr*prefac;
}

/*
	Compute the optimal observable 
	Compute covariance matrix elements for k=1,...10, l=1,...10	
*/
float Operator(float pp, float pe, int k, float t0, float t1, float f1, float t2, float f2){
	float expr;
	expr =((1+pp)*(1-pe)*Sleft (k,t0,t1,f1,t2,f2)
	      +(1-pp)*(1+pe)*Sright(k,t0,t1,f1,t2,f2))
	     /((1+pp)*(1-pe)*Sleft (0,t0,t1,f1,t2,f2)
	      +(1-pp)*(1+pe)*Sright(0,t0,t1,f1,t2,f2));
    
	return expr;
}

/*
	Read LHE file and compute everything
	Store values in a root file as ntuple
*/
void opt_obs(){

	ifstream file;
	int n;

	file.open(lhe_name.c_str());

	TNtuple *ntuple = new TNtuple("ntuple","ntuple","t1t:t1eta:t1y:t1phi:t1theta:t1m:t1e:t2t:t2eta:t2y:t2phi:t2theta:t2m:t2e:mtt:cosphiWp:sinphiWp:costhetaWp:cosphiWm:sinphiWm:costhetaWm:phiWp:phiWp2:thetaWp:phiWm:phiWm2:thetaWm:SM:lqA:eqA:pqA:lqV:eqV:pqV:ReuZ:ReuA:ImuZ:ImuA:dW:pqP");
	
	// definitions and initializations
	float crosssection = 0.;				// from <init> tag
	float uncertainty = 0.;					// its uncertainty
	float weightssum = 0.;					// sum of event weights
	vector<float> means(10);				// values of the optimal observables
	vector<vector<float>> covariance(10, vector<float>(10));// stripped inverse covariance matrix: cov(Ci,Cj)⁻¹/(Lumi*eff) in fb
	for( int index=0; index<10; ++index){
		means[index] = 0.;
		for( int jndex=0; jndex<10; ++jndex){
			covariance[index][jndex] = 0.;
		}
	}
	int counter=0; // count file lines
	int nevt=0;    // count events
	float afb_top, err_afb_top;
	float afb_atop, err_afb_atop;


	if (file) {

		int no = 0;
		std::cout <<std::endl<<"Running on "+lhe_name+" assuming: sqrt(s) = "<< sqs << " GeV,  P(e+,e-) = ("<< polp << "," << pole << ")" << std::endl<< std::endl;

		while (!file.eof() && counter < 400000000) {

			string store;
			file >> store;
			
			// get total cross section from <init>
			if (store.find("<init>") == 0 && no == 0){
				float y[10];
				float z[4];
				file
					>> y[0] >> y[1] // beam ids
					>> y[2] >> y[3]	// beam energies
					>> y[4] >> y[5] >> y[6] >> y[7] // pdf type and ids
					>> y[8] // weight strategy
					>> y[9];// number of subprocesses
				for(int index=0; index<y[9]; index++){
					file
						>> z[0] // cross section for this subprocess
						>> z[1] // uncertaincy on the cross section
						>> z[2] // maximum weight of events
						>> z[3];// process index
					crosssection += z[0];
					uncertainty += z[1]*z[1];
				}
				//std::cout <<"Cross-section [pb]: " << z[0] << " +- " << z[1] << "." << std::endl;
				no = 1;
			}
		  
			// get event information
			if (store.find("<event>") == 0) {

				nevt++;

				int np=0;	// number of particles in the event
				int nn=0;	// id of the process
				float x[10];	// weight, c.o.m. energy, ...
				file >> np >> nn >> x[0] >> x[1] >> x[2] >> x[3];
				weightssum += x[0];

				int pdg;	// particle PDG id
				int id[4];	// status (-1=initial, 1=final, 2=intermediate), mother particles, and colour flows
				float p[4];	// particle four momentum
				float mm;	// particle mass
				float xx[2];	// 
				
				// four vectors
				TLorentzVector ELECBEAM(0,0,0,0);
				TLorentzVector  POSBEAM(0,0,0,0);

				TLorentzVector  TOP(0,0,0,0);
				TLorentzVector ATOP(0,0,0,0);

				TLorentzVector WP(0,0,0,0);
				TLorentzVector WM(0,0,0,0);

				TLorentzVector  BQUARK(0,0,0,0);
				TLorentzVector ABQUARK(0,0,0,0);

				// W⁻ decay products
				TLorentzVector LEPM(0,0,0,0);
				TLorentzVector ANEU(0,0,0,0);
				TLorentzVector QUARKU(0,0,0,0);
				TLorentzVector AQUARKD(0,0,0,0);

				// W⁺ decay products
				TLorentzVector AQUARKU(0,0,0,0);
				TLorentzVector QUARKD(0,0,0,0);
				TLorentzVector LEPP(0,0,0,0);
				TLorentzVector NEU(0,0,0,0);
				
				// loop over particles
				for (int index=0;index<np;index++) {
					file >> pdg >> id[0] >> id[1] >> id[2] >> id[3] >> id[4] >> p[0] >> p[1] >> p[2] >> p[3] >> mm >> xx[0] >> xx[1];
					//1-->down	2-->up		3-->strange
					//4-->charm	5-->bottom	6-->top
					//7-->	8-->	9-->	10-->
					//11-->electron	12-->electron neutrino
					//13-->muon	14-->muon neutrino
					//15-->tau	16-->tau neutrino
					//17-->	18-->	19-->	20-->
					//21-->gluón	22-->fotón	23-->	24-->W boson
					if ( pdg== 6)                                         TOP.SetPxPyPzE(p[0],p[1],p[2],p[3]);
					if ( pdg==-6)                                        ATOP.SetPxPyPzE(p[0],p[1],p[2],p[3]);
					if ( pdg== 5 && id[0]==1)                          BQUARK.SetPxPyPzE(p[0],p[1],p[2],p[3]);
					if ( pdg==-5 && id[0]==1)                         ABQUARK.SetPxPyPzE(p[0],p[1],p[2],p[3]);
					if ( pdg== 24)                                         WP.SetPxPyPzE(p[0],p[1],p[2],p[3]);
					if ( pdg==-24)                                         WM.SetPxPyPzE(p[0],p[1],p[2],p[3]);
					if ((pdg== 11 || pdg== 13 || pdg== 15) && id[0]==1)  LEPM.SetPxPyPzE(p[0],p[1],p[2],p[3]);
					if ((pdg==-11 || pdg==-13 || pdg==-15) && id[0]==1)  LEPP.SetPxPyPzE(p[0],p[1],p[2],p[3]);
					if ( pdg== 11 && id[0]==-1)                      ELECBEAM.SetPxPyPzE(p[0],p[1],p[2],p[3]);
					if ( pdg==-11 && id[0]==-1)                       POSBEAM.SetPxPyPzE(p[0],p[1],p[2],p[3]);
					if ( pdg== 1 || pdg== 3)                           QUARKD.SetPxPyPzE(p[0],p[1],p[2],p[3]);
					if ( pdg==-1 || pdg==-3)                          AQUARKD.SetPxPyPzE(p[0],p[1],p[2],p[3]);
					if ( pdg== 2 || pdg== 4)                           QUARKU.SetPxPyPzE(p[0],p[1],p[2],p[3]);
					if ( pdg==-2 || pdg==-4)                          AQUARKU.SetPxPyPzE(p[0],p[1],p[2],p[3]);


				}
				 
				// store kinematic variables
				// assume information about the top and W were present in the event
				float t[100];
				t[0]  = TOP.Pt();
				t[1]  = TOP.Eta();
				t[2]  = TOP.Rapidity();
				t[3]  = TOP.Phi();
				t[4]  = TOP.Theta();
				t[5]  = TOP.M();
				t[6]  = TOP.E();

				t[7]  = ATOP.Pt();
				t[8]  = ATOP.Eta();
				t[9]  = ATOP.Rapidity();
				t[10] = ATOP.Phi();
				t[11] = ATOP.Theta();
				t[12] = ATOP.M();
				t[13] = ATOP.E();

				t[14] = (TOP+ATOP).M();

				 
				std::vector<float> helP, helM;
				 
				helP = angle_calc( TOP, WP, ELECBEAM);
				helM = angle_calc(ATOP, WM,  POSBEAM);
				helM[0] *= -1.;
				helM[1] *= -1.;


				t[15] = helP[0];
				t[16] = helP[1];
				t[17] = helP[2];

				t[18] = helM[0];
				t[19] = helM[1];
				t[20] = helM[2];

				float phi0M, phi1M, phi0P, phi1P;
				 
				if(helP[0]>0 && helP[1]>0){
					phi0P=TMath::ACos(helP[0]);
					phi1P=TMath::ASin(helP[1]);
				}
				if(helP[0]<0 && helP[1]>0){
					phi0P=TMath::ACos(helP[0]);
					phi1P=pi-TMath::ASin(helP[1]);
				}
				if(helP[0]>0 && helP[1]<0){
					phi0P=-1*TMath::ACos(helP[0]);
					phi1P=TMath::ASin(helP[1]);
				}
				if(helP[0]<0 && helP[1]<0){
					phi0P=-1*TMath::ACos(helP[0]);
					phi1P=-1*(TMath::ASin(helP[1])+pi);
				}

				if(helM[0]>0 && helM[1]>0){
					phi0M=TMath::ACos(helM[0]);
					phi1M=TMath::ASin(helM[1]);
				}
				if(helM[0]<0 && helM[1]>0){
					phi0M=TMath::ACos(helM[0]);
					phi1M=pi-TMath::ASin(helM[1]);
				}
				if(helM[0]>0 && helM[1]<0){
					phi0M=-1*TMath::ACos(helM[0]);
					phi1M=TMath::ASin(helM[1]);
				}
				if(helM[0]<0 && helM[1]<0){
					phi0M=-1*TMath::ACos(helM[0]);
					phi1M=-1*(TMath::ASin(helM[1])+pi);
				}

				t[21] = phi0P;
				t[22] = phi1P;
				t[23] = TMath::ACos(helP[2]);

				t[24] = phi0M;
				t[25] = phi1M;
				t[26] = TMath::ACos(helM[2]);
				
				// the top production polar angle
				float costheta0 = TOP.Vect().Dot( ELECBEAM.Vect() )/ TOP.Vect().Mag() / ELECBEAM.Vect().Mag();
				
				// value of the optimal observable for that event	
				for(int index=0;index<13;index++){
					t[27+index] = Operator(polp, pole, index, TMath::ACos(costheta0), TMath::ACos(helP[2]), phi0P, TMath::ACos(helM[2]), phi0M);
				}
				// optimal observable sample mean and covariance matrix, weighted by the event weight x[0]
				for( int index=0; index<10; ++index){
					means[index] += t[28+index]*x[0];
					for( int jndex=0; jndex<10; ++jndex){
						covariance[index][jndex] += t[28+index]*t[28+jndex]*x[0];
					}
				}
				 
			ntuple->Fill(&t[0]);
			}

		counter++;
		}

	}
	
	// print results
	uncertainty = sqrt(uncertainty);
	cout << nevt << " events,\t" << crosssection << " pb +- " << uncertainty << " (weights' sum / nevt = " << weightssum/nevt << " pb)" << endl;
	cout << "--- optimal observable covariance matrix divided by luminosity cov(Oi,Oj)/Lumi [pb] (= cov(Ci,Cj)^-1/Lumi ) ---" << endl;
	for( int index=0; index<10; ++index){
		for( int jndex=0; jndex<10; ++jndex){
			covariance[index][jndex] *= 1./nevt;
			cout << covariance[index][jndex] << "	";
		}
		cout << endl;
	}
	cout << "--- optimal observable means divided by luminosity [pb] (= linear EFT dependence) ---" << endl;
	for( int index=0; index<10; ++index){
		means[index] *= 1./nevt;
		cout << means[index] << endl;
	}
	cout << "---" << endl;

    
	TFile *f = new TFile(root_name.c_str(),"RECREATE");
	ntuple->Write();
	f->Close();
}

int main(){
	opt_obs();
}
