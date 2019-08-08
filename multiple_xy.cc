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

int MSW;


int ip;
cout << "The position ip";
cin >> ip;

cout << "Number of X unit cells Nx";
cin >> Nx;
cout << "The magnetic coupling Jx";
cin >> Jx;
cout << "The magnetic coupling Jy";
cin >> Jy;
cout << "The magnetic coupling Jz";
cin >> Jz;
cout << "The number of sweep MSW";
cin >> MSW;
// read in the parameters

Jx=Jx*4;
Jy=Jy*4; // The spin operator here gives eigenvalues of +/- 1/2 rather than +/- 1
Jz=Jz*4;

int N=2;

int N_=4*N*Nx;
        
int NDD;

int NDY=4*Nx-1;
int NDZ=8*Nx-2;

auto sites_ = tJ(N_);

auto    H = MPO(sites_);

//ofstream myfile1;
//myfile1.open ("density_t"+ std::to_string(t) +"_Jx_"+std::to_string(Jx)+"_Jy_"+std::to_string(Jy)+"_Jz_"+std::to_string(Jz)+"_Nx_"+std::to_string(Nx)+".dat");

std::vector<Index> links(N_+1); // This is the vector of bond dimension


    for(int l = 0; l <= N_; ++l) 
        {
        links.at(l) = Index(nameint("hl",l),5+NDY+NDZ);
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


       for (int jj1 = 5; jj1 <=NDY+3; ++jj1)		{W += sites_.op("Id",n) * row(jj1) * col(jj1-1); } //ok
       for (int jj2 = 6+NDY; jj2 <=NDY+NDZ+5; ++jj2)	{W += sites_.op("Id",n) * row(jj2) * col(jj2-1); }

        if (n%2==1)
{
        W += sites_.op("Sm",n) * row(2) * col(3) * Jx *0.5;		//NN	
        W += sites_.op("Sp",n) * row(2) * col(3) * Jx *0.5;		//NN
}

        if (n%2==0 && n%(4*Nx)!=0)
{
        W += -sites_.op("Sm",n) * row(2) * col(4+NDY) * Jy *0.5;		//NN
        W += sites_.op("Sp",n) * row(2) * col(4+NDY) * Jy *0.5;			//NN
    }



        if (n%(4*Nx)==1)
{
        W += -sites_.op("Sm",n) * row(2) * col(3+NDY) * Jy *0.5;
        W += sites_.op("Sp",n) * row(2) * col(3+NDY) * Jy*0.5;    
}


	if (n%2==1 && (n<=4*Nx) && n>1)		{ NDD=N_-2*n+2;   W += Jz*sites_.op("Sz",n) * row(2)* col(4+NDY+NDD);}  //ok

        if (n==1)				{ W += Jz*sites_.op("Sz",n) * row(2)* col(4+NDY+4*Nx); }   //ok

	if (n%2==0 && (n>4*Nx))			{ NDD=2*(N_-n)+2; W += Jz*sites_.op("Sz",n) *row(2)* col(4+NDY+NDD);}	//ok


}


    auto sweeps = Sweeps(MSW);
    sweeps.maxm() = 20,40,80,120,200,300;
    sweeps.cutoff() = 1E-15,Args("Repeat",20),1E-15;
    sweeps.niter() = 3,2;

    auto state = InitState(sites_);
    for(int i = 1; i <= N_; ++i) 
        {
          state.set(i,"Up");
        }
    auto psi = MPS(state);


    auto LH = setElt(links.at(0)(2));
    auto RH = setElt(dag(last)(1));

    H.Anc(0) = LH;
    H.Anc(N_+1) = RH;


auto res = idmrg(psi,H,sweeps,{"OutputLevel",1});



normalize(psi);


for (int j=ip+1;j<=N_;++j)
{
//Given an MPS or IQMPS called "psi",
//constructed from a SiteSet "sites"


//Replace "Op1" and "Op2" with the actual names
//of the operators you want to measure
auto op_i = sites_.op("Sz",ip);
auto op_j = sites_.op("Sz",j);

//below we will assume j > i

//'gauge' the MPS to site i
//any 'position' between i and j, inclusive, would work here
psi.position(ip); 

psi.Anc(1) *= psi.A(0); //Uncomment if doing iDMRG calculation

//index linking i to i+1:
auto ir = commonIndex(psi.A(ip),psi.A(ip+1),Link);

auto C = psi.A(ip)*op_i*dag(prime(psi.A(ip),Site,ir));
for(int k = ip+1; k < j; ++k)
    {
    C *= psi.A(k);
    C *= dag(prime(psi.A(k),Link));
    }
C *= psi.A(j);
C *= op_j;
//index linking j to j-1:
auto jl = commonIndex(psi.A(j),psi.A(j-1),Link);
C *= dag(prime(psi.A(j),jl,Site));

auto result = C.real(); //or C.cplx() if expecting complex

println(j,result);



}


    return 0;
    }
