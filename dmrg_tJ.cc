#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <math.h>       /* exp */
using namespace std;
using namespace itensor;

int main()
    {
int Nx;

double t;

double Jx;

double Jy;

double Jz;

int jp;


int MSW;

cout << "Hopping amplitude t";
cin >> t;
cout << "Number of X unit cells Nx";
cin >> Nx;
cout << "The magnetic coupling Jx";
cin >> Jx;
cout << "The magnetic coupling Jy";
cin >> Jy;
cout << "The magnetic coupling Jz";
cin >> Jz;

cout << "The initial position";
cin >> jp;

cout << "The number of sweep MSW";
cin >> MSW;

Jx=Jx*4;
Jy=Jy*4;
Jz=Jz*4;

// read in the parameters

int N=2;

int N_=4*N*Nx;
        
int NDD;


auto sites_ = tJ(N_);



auto    H = MPO(sites_);

//ofstream myfile;
//myfile.open ("density_t"+ std::to_string(t) +"_Jx_"+std::to_string(Jx)+"_Jy_"+std::to_string(Jy)+"_Jz_"+std::to_string(Jz)+"_Nx_"+std::to_string(Nx)+".dat");


std::vector<Index> links(N_+1);



int NDX=4*Nx-1;
int NDY=4*Nx+1;
int NDZ=4*Nx-1;

    for(int l = 0; l <= N_; ++l) 
        {
        links.at(l) = Index(nameint("hl",l),5*NDY+NDZ+9);
        }

    Index const& last = (links.at(0));

    for(int n = 1; n <= N_; ++n)
        {

        auto& W = H.Anc(n);
        auto row = dag(links.at(n-1));
        auto col = (n==N_ ? last : links.at(n));

        W = ITensor(dag(sites_(n)),prime(sites_(n)),row,col);

        W += sites_.op("Id",n) * row(1) * col(1); //ending state
        W += sites_.op("Id",n) * row(2) * col(2); //starting state

        W += sites_.op("Sp",n) * row(3) * col(1)*0.5;	        	//X
        W += sites_.op("Sm",n) * row(3) * col(1)*0.5;	        	//X

        W += sites_.op("Sm",n) * row(4) * col(1)*0.5;			//Y
        W +=-sites_.op("Sp",n) * row(4) * col(1)*0.5;			//Y

        W += sites_.op("Sm",n) * row(4+NDY) * col(1)*0.5;		//Y
        W +=-sites_.op("Sp",n) * row(4+NDY) * col(1)*0.5;		//Y

        W += sites_.op("Sz",n) * row(5+NDY) * col(1);			//Z

        W += sites_.op("Cup",n) * row(NDY+NDZ+6) * col(1);			//Cup
        W += sites_.op("Cup",n) * row(2*NDY+NDZ+6) * col(1);			//Cup

        W += sites_.op("Cdagup",n) * row(2*NDY+NDZ+7) * col(1);			//Cdagup
        W += sites_.op("Cdagup",n) * row(3*NDY+NDZ+7) * col(1);			//Cdagup

        W += sites_.op("Cdn",n) * row(3*NDY+NDZ+8) * col(1);			//Cdn
        W += sites_.op("Cdn",n) * row(4*NDY+NDZ+8) * col(1);			//Cdn

        W += sites_.op("Cdagdn",n) * row(4*NDY+NDZ+9) * col(1);			//Cdagdn
        W += sites_.op("Cdagdn",n) * row(5*NDY+NDZ+9) * col(1);			//Cdagdn

       for (int jj1 = 5; jj1 <=NDY+3; ++jj1)				{W += sites_.op("Id",n) * row(jj1) * col(jj1-1); } //ok
       for (int jj2 = 6+NDY; jj2 <=NDY+NDZ+5; ++jj2)			{W += sites_.op("Id",n) * row(jj2) * col(jj2-1); }

       for (int jj1 = NDY+NDZ+7; jj1 <=2*NDY+NDZ+5; ++jj1)		{W += sites_.op("Id",n) * row(jj1) * col(jj1-1); } //ok
       for (int jj2 = 2*NDY+NDZ+8; jj2 <=3*NDY+NDZ+6; ++jj2)		{W += sites_.op("Id",n) * row(jj2) * col(jj2-1); }


       for (int jj1 = 3*NDY+NDZ+9; jj1 <=4*NDY+NDZ+7; ++jj1)		{W += sites_.op("Id",n) * row(jj1) * col(jj1-1); } //ok
       for (int jj2 = 4*NDY+NDZ+10; jj2 <=5*NDY+NDZ+8; ++jj2)		{W += sites_.op("Id",n) * row(jj2) * col(jj2-1); }


        if (n%2==1)
{
        W += sites_.op("Sm",n) * row(2) * col(3) * Jx *0.5;		//NN	
        W += sites_.op("Sp",n) * row(2) * col(3) * Jx *0.5;		//NN

        W += -t*sites_.op("Cup",n) * row(2) * col(3*NDY+NDZ+7) ;		//NN	
        W += -t*sites_.op("Cdn",n) * row(2) * col(5*NDY+NDZ+9) ;		//NN
        W += t*sites_.op("Cdagup",n) * row(2) * col(2*NDY+NDZ+6) ;		//NN	
        W += t*sites_.op("Cdagdn",n) * row(2) * col(4*NDY+NDZ+8) ;		//NN


}

        if (n%2==0 && n%(4*Nx)!=0)
{
        W += -sites_.op("Sm",n) * row(2) * col(4+NDY) * Jy *0.5;		//NN
        W += sites_.op("Sp",n) * row(2) * col(4+NDY) * Jy *0.5;			//NN

        W += -t*sites_.op("Cup",n) * row(2) * col(3*NDY+NDZ+7) ;		//NN	
        W += -t*sites_.op("Cdn",n) * row(2) * col(5*NDY+NDZ+9) ;		//NN
        W += t*sites_.op("Cdagup",n) * row(2) * col(2*NDY+NDZ+6) ;		//NN	
        W += t*sites_.op("Cdagdn",n) * row(2) * col(4*NDY+NDZ+8) ;		//NN

    }



        if (n%(4*Nx)==1)
{
        W += -sites_.op("Sm",n) * row(2) * col(3+NDY) * Jy *0.5;
        W += sites_.op("Sp",n) * row(2) * col(3+NDY) * Jy*0.5;


        W += -t*sites_.op("Cup",n) * row(2) * col(NDY+2*NDZ+5) ;		//NN	
        W += -t*sites_.op("Cdn",n) * row(2) * col(2*NDY+2*NDZ+6) ;		//NN
        W += t*sites_.op("Cdagup",n) * row(2) * col(3*NDY+2*NDZ+7) ;		//NN	
        W += t*sites_.op("Cdagdn",n) * row(2) * col(4*NDY+2*NDZ+8) ;		//NN
   
}


	if (n%2==1 && (n<=4*Nx) && n>1)
{

NDD=N_-2*n+2;   
W += Jz*sites_.op("Sz",n) * row(2)* col(4+NDY+NDD);

W += -t*sites_.op("Cup",n) * row(2) * col(NDY+NDZ+5+NDD) ;		//NN	
W += -t*sites_.op("Cdn",n) * row(2) * col(2*NDY+NDZ+6+NDD) ;		//NN
W += t*sites_.op("Cdagup",n) * row(2) * col(3*NDY+NDZ+7+NDD) ;		//NN	
W += t*sites_.op("Cdagdn",n) * row(2) * col(4*NDY+NDZ+8+NDD) ;		//NN
}  //ok

        if (n==1)				
{
  W += Jz*sites_.op("Sz",n) * row(2)* col(4+NDY+4*Nx);
 
  W += -t*sites_.op("Cup",n) * row(2) * col(NDY+NDZ+5+4*Nx) ;		//NN	
  W += -t*sites_.op("Cdn",n) * row(2) * col(2*NDY+NDZ+6+4*Nx) ;		//NN
  W += t*sites_.op("Cdagup",n) * row(2) * col(3*NDY+NDZ+7+4*Nx) ;		//NN	
  W += t*sites_.op("Cdagdn",n) * row(2) * col(4*NDY+NDZ+8+4*Nx) ;		//NN
}   //ok

	if (n%2==0 && (n>4*Nx))
{ 
NDD=2*(N_-n)+2; 
W += Jz*sites_.op("Sz",n) *row(2)* col(4+NDY+NDD);

W += -t*sites_.op("Cup",n) * row(2) * col(NDY+NDZ+5+NDD) ;		//NN	
W += -t*sites_.op("Cdn",n) * row(2) * col(2*NDY+NDZ+6+NDD) ;		//NN
W += t*sites_.op("Cdagup",n) * row(2) * col(3*NDY+NDZ+7+NDD) ;		//NN	
W += t*sites_.op("Cdagdn",n) * row(2) * col(4*NDY+NDZ+8+NDD) ;		//NN
}	//ok


}


    auto sweeps = Sweeps(MSW);
    sweeps.maxm() = 20,40,80,120,200,300;
    sweeps.cutoff() = 1E-15,Args("Repeat",20),1E-15;
    sweeps.niter() = 3,2;

    sweeps.noise() = 1E-7,1E-8,0.0;

    auto state = InitState(sites_);
    for(int i = 1; i <= N_; ++i) 
        {
          state.set(i,"Up");
        }
 state.set(jp,"Emp");

    auto psi = MPS(state);


    auto LH = setElt(links.at(0)(2));
    auto RH = setElt(dag(last)(1));

    H.Anc(0) = LH;
    H.Anc(N_+1) = RH;


auto res = idmrg(psi,H,sweeps,{"OutputLevel",1});



normalize(psi);


    return 0;
    }

