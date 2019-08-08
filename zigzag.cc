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
;

int N=2;

int N_=4*N*Nx;
        
int NDX=4*Nx-1;
int NDY=4*Nx+1;
int NDZ=4*Nx-1;

auto sites_ = tJ(N_);

int NC=0;


auto    H = MPO(sites_);

ofstream myfile;
myfile.open ("Sweep_NO"+ std::to_string(MSW) +"_Jx_"+std::to_string(Jx)+"_Jy_"+std::to_string(Jy)+"_Jz_"+std::to_string(Jz)+"_Nx_"+std::to_string(Nx)+".dat");


Jx=Jx*4;
Jy=Jy*4;
Jz=Jz*4;

std::vector<Index> links(N_+1);

//The names of these indices refer to their Sz quantum numbers
std::vector<Index> q0(N_+1),
                   qP(N_+1),
                   qM(N_+1);

    for(int l = 0; l <= N_; ++l) 
        {
        links.at(l) = Index(nameint("hl",l),NDX+NDY+NDZ+5);
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
        W += sites_.op("Sz",n) * row(3) * col(1);				//Z
        W += sites_.op("Sz",n) * row(3+NDZ) * col(1);				//Z

        W += sites_.op("Sp",n) * row(4+NDZ) * col(1)*0.5;	       		//X
        W += sites_.op("Sm",n) * row(4+NDZ) * col(1)*0.5;	    		//X

        W += sites_.op("Sp",n) * row(NDX+NDZ+4) * col(1)*0.5;			//X
        W += sites_.op("Sm",n) * row(NDX+NDZ+4) * col(1)*0.5;			//X

        W += sites_.op("Sm",n) * row(NDX+NDZ+5) * col(1)*0.5;			//Y
        W +=-sites_.op("Sp",n) * row(NDX+NDZ+5) * col(1)*0.5;			//Y
        W += sites_.op("Sm",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5;		//Y
        W +=-sites_.op("Sp",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5;		//Y

        for (int jj = 4; jj <=NDZ+2; ++jj)          {W += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 5+NDZ; jj <=NDX+NDZ+3; ++jj)          {W += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 6+NDX+NDZ; jj <=NDX+NDY+NDZ+4; ++jj)  {W += sites_.op("Id",n) * row(jj) * col(jj-1);}


        if (n%4==1 )
{
        W += sites_.op("Sm",n) * row(2) * col(NDX+NDZ+4) * Jx *0.5;
        W += sites_.op("Sp",n) * row(2) * col(NDX+NDZ+4) * Jx *0.5;

        W += -sites_.op("Sm",n) * row(2) * col(NDX+NDY+NDZ+4) * Jy *0.5;		// pb
        W +=  sites_.op("Sp",n) * row(2) * col(NDX+NDY+NDZ+4) * Jy *0.5;		// pb
}


        if (n%4==2)
{
       W += sites_.op("Sz",n) * row(2) * col(3+NDZ) * Jz;
}



        if (n%4==3)
{
        W += -sites_.op("Sm",n) * row(2) * col(NDX+NDY+NDZ+5) * Jy *0.5;
        W +=  sites_.op("Sp",n) * row(2) * col(NDX+NDY+NDZ+5) * Jy *0.5;
}


        if (n%4==0 && ( n%(4*Nx)!=0 ) )
{
        W += sites_.op("Sz",n) * row(2) * col(3+NDZ) * Jz;
        W += sites_.op("Sm",n) * row(2) * col(NDX+NDZ+3) * Jx*0.5;		// pb
        W += sites_.op("Sp",n) * row(2) * col(NDX+NDZ+3) * Jx*0.5;		// pb
}


        if ((n%(4*Nx)==0))
{
        W += sites_.op("Sz",n) * row(2) * col(2+NDZ) * Jz;
	}

}
    auto sweeps = Sweeps(MSW);
    sweeps.maxm() = 40,120,400,800,2000;
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

