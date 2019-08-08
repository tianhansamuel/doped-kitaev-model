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

// read in the parameters

int N=1;

int N_=8*N*Nx;
        
int NDX=4*Nx+1;
int NDY=4*Nx-1;


auto sites_ = tJ(N_);



auto    HtJ = MPO(sites_);

ofstream myfile;
myfile.open ("density_t"+ std::to_string(t) +"_Jx_"+std::to_string(Jx)+"_Jy_"+std::to_string(Jy)+"_Jz_"+std::to_string(Jz)+"_Nx_"+std::to_string(Nx)+".dat");


std::vector<Index> links1(N_+1);




    auto sweeps = Sweeps(10);
    sweeps.maxm() = 20,80,140,200;
    sweeps.cutoff() = 1E-10,Args("Repeat",10),1E-14;
    sweeps.niter() = 3,2;

    auto state = InitState(sites_);
    for(int i = 1; i <= N_; ++i) 
        {
        if(i%2 == 1)
            state.set(i,"Up");
        else
            state.set(i,"Dn");
        }

    state.set(jp,"Emp");
    auto psi = MPS(state);



/////////////////////// The tJ Hamiltonian

for(int l = 0; l <= N_; ++l) 
        {
        links1.at(l) = Index(nameint("hl",l),5*NDX+NDY+9);
        }



    Index const& last1 = (links1.at(0));

    for(int n = 1; n <= N_; ++n)
        {

        auto& W1 = HtJ.Anc(n);
        auto row = dag(links1.at(n-1));
        auto col = (n==N_ ? last1 : links1.at(n));

        W1 = ITensor(dag(sites_(n)),prime(sites_(n)),row,col);
    for(int l = 0; l <= N_; ++l) 
        {
        links1.at(l) = Index(nameint("hl",l),NDX+5*NDY+NDZ+9);
        }

    Index const& last = (links1.at(0));

    for(int n = 1; n <= N_; ++n)
        {

        auto& W1 = H.Anc(n);
        auto row = dag(links1.at(n-1));
        auto col = (n==N_ ? last : links1.at(n));

        W1 = ITensor(dag(sites_(n)),prime(sites_(n)),row,col);

        if (n%4==1 && (n%(4*Nx)!=1))
{
        W1 += sites_.op("Id",n) * row(1) * col(1); //ending state
        W1 += sites_.op("Id",n) * row(2) * col(2); //starting state
        W1 += sites_.op("Sz",n) * row(3) * col(1);					//Z
        W1 += sites_.op("Sz",n) * row(3+NDZ) * col(1);					//Z
        W1 += sites_.op("Sp",n) * row(4+NDZ) * col(1)*0.5;	        		//X
        W1 += sites_.op("Sm",n) * row(4+NDZ) * col(1)*0.5;	    			//X
        W1 += sites_.op("Sp",n) * row(NDX+NDZ+4) * col(1)*0.5;				//X
        W1 += sites_.op("Sm",n) * row(NDX+NDZ+4) * col(1)*0.5;				//X
        W1 += sites_.op("Sm",n) * row(NDX+NDZ+5) * col(1)*0.5*Cplx_i;			//Y
        W1 +=-sites_.op("Sp",n) * row(NDX+NDZ+5) * col(1)*0.5*Cplx_i;			//Y
        W1 += sites_.op("Sm",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5*Cplx_i;		//Y
        W1 +=-sites_.op("Sp",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5*Cplx_i;		//Y

        W1 += sites_.op("Sz",n) * row(2) * col(3+NDZ)*Jz;					//Z

        for (int jj = 4; jj <=NDZ+2; ++jj)          {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 5+NDZ; jj <=NDX+NDZ+3; ++jj)          {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 6+NDX+NDZ; jj <=NDX+NDY+NDZ+4; ++jj)  {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}

        W1 += sites_.op("Cup",n) * row(NDX+NDY+NDZ+6) * col(1);			//Cup
        W1 += sites_.op("Cup",n) * row(NDX+2*NDY+NDZ+6) * col(1);		//Cup		NN
        W1 += sites_.op("Cdn",n) * row(NDX+2*NDY+NDZ+7) * col(1);		//Cdn
        W1 += sites_.op("Cdn",n) * row(NDX+3*NDY+NDZ+7) * col(1);		//Cdn		NN

        W1 += sites_.op("Cdagup",n) * row(NDX+3*NDY+NDZ+8) * col(1);		//Cup
        W1 += sites_.op("Cdagup",n) * row(NDX+4*NDY+NDZ+8) * col(1);		//Cup		NN
        W1 += sites_.op("Cdagdn",n) * row(NDX+4*NDY+NDZ+9) * col(1);		//Cdn
        W1 += sites_.op("Cdagdn",n) * row(NDX+5*NDY+NDZ+9) * col(1);		//Cdn		NN

        for (int jj = NDX+NDY+NDZ+7; jj <=NDX+2*NDY+NDZ+5; ++jj)                {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+2*NDY+NDZ+8; jj <=NDX+3*NDY+NDZ+6; ++jj)              {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+3*NDY+NDZ+9; jj <=NDX+4*NDY+NDZ+7; ++jj)	        {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+4*NDY+NDZ+10; jj <=NDX+5*NDY+NDZ+8; ++jj)  		{W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}

        W1 += -t*sites_.op("Cdagup",n) * row(2) * col(NDX+2*NDY+NDZ+6);					//Z
        W1 += -t*sites_.op("Cdagdn",n) * row(2) * col(NDX+3*NDY+NDZ+7);					//Z
        W1 += t*sites_.op("Cdagup",n) * row(2) * col(NDX+4*NDY+NDZ+8);					//Z
        W1 += t*sites_.op("Cdagdn",n) * row(2) * col(NDX+5*NDY+NDZ+9);					//Z


}


        if (n%(4*Nx)==1)
{
        W1 += sites_.op("Id",n) * row(1) * col(1); //ending state
        W1 += sites_.op("Id",n) * row(2) * col(2); //starting state
        W1 += sites_.op("Sz",n) * row(3) * col(1);					//Z
        W1 += sites_.op("Sz",n) * row(3+NDZ) * col(1);					//Z
        W1 += sites_.op("Sp",n) * row(4+NDZ) * col(1)*0.5;	        		//X
        W1 += sites_.op("Sm",n) * row(4+NDZ) * col(1)*0.5;	    			//X
        W1 += sites_.op("Sp",n) * row(NDX+NDZ+4) * col(1)*0.5;				//X
        W1 += sites_.op("Sm",n) * row(NDX+NDZ+4) * col(1)*0.5;				//X
        W1 += sites_.op("Sm",n) * row(NDX+NDZ+5) * col(1)*0.5*Cplx_i;			//Y
        W1 +=-sites_.op("Sp",n) * row(NDX+NDZ+5) * col(1)*0.5*Cplx_i;			//Y
        W1 += sites_.op("Sm",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5*Cplx_i;		//Y
        W1 +=-sites_.op("Sp",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5*Cplx_i;		//Y


        for (int jj = 4; jj <=NDZ+2; ++jj)          {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 5+NDZ; jj <=NDX+NDZ+3; ++jj)          {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 6+NDX+NDZ; jj <=NDX+NDY+NDZ+4; ++jj)  {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}


        W1 += sites_.op("Cup",n) * row(NDX+NDY+NDZ+6) * col(1);			//Cup
        W1 += sites_.op("Cup",n) * row(NDX+2*NDY+NDZ+6) * col(1);		//Cup		NN
        W1 += sites_.op("Cdn",n) * row(NDX+2*NDY+NDZ+7) * col(1);		//Cdn
        W1 += sites_.op("Cdn",n) * row(NDX+3*NDY+NDZ+7) * col(1);		//Cdn		NN

        W1 += sites_.op("Cdagup",n) * row(NDX+3*NDY+NDZ+8) * col(1);		//Cup
        W1 += sites_.op("Cdagup",n) * row(NDX+4*NDY+NDZ+8) * col(1);		//Cup		NN
        W1 += sites_.op("Cdagdn",n) * row(NDX+4*NDY+NDZ+9) * col(1);		//Cdn
        W1 += sites_.op("Cdagdn",n) * row(NDX+5*NDY+NDZ+9) * col(1);		//Cdn		NN

        for (int jj = NDX+NDY+NDZ+7; jj <=NDX+2*NDY+NDZ+5; ++jj)                {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+2*NDY+NDZ+8; jj <=NDX+3*NDY+NDZ+6; ++jj)              {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+3*NDY+NDZ+9; jj <=NDX+4*NDY+NDZ+7; ++jj)	        {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+4*NDY+NDZ+10; jj <=NDX+5*NDY+NDZ+8; ++jj)  		{W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}




}



        if (n%4==2)
{
        W1 += sites_.op("Id",n) * row(1) * col(1); //ending state
        W1 += sites_.op("Id",n) * row(2) * col(2); //starting state
        W1 += sites_.op("Sz",n) * row(3) * col(1);					//Z
        W1 += sites_.op("Sz",n) * row(3+NDZ) * col(1);					//Z

        W1 += sites_.op("Sp",n) * row(4+NDZ) * col(1)*0.5;	        		//X
        W1 += sites_.op("Sm",n) * row(4+NDZ) * col(1)*0.5;	    			//X

        W1 += sites_.op("Sp",n) * row(NDX+NDZ+4) * col(1)*0.5;				//X
        W1 += sites_.op("Sm",n) * row(NDX+NDZ+4) * col(1)*0.5;				//X

        W1 += sites_.op("Sm",n) * row(NDX+NDZ+5) * col(1)*0.5*Cplx_i;			//Y
        W1 +=-sites_.op("Sp",n) * row(NDX+NDZ+5) * col(1)*0.5*Cplx_i;			//Y
        W1 += sites_.op("Sm",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5*Cplx_i;		//Y
        W1 +=-sites_.op("Sp",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5*Cplx_i;		//Y

        for (int jj = 4; jj <=NDZ+2; ++jj)          {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 5+NDZ; jj <=NDX+NDZ+3; ++jj)          {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 6+NDX+NDZ; jj <=NDX+NDY+NDZ+4; ++jj)  {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}

        W1 += sites_.op("Sm",n) * row(2) * col(NDX+NDZ+4) * Jx*0.5;
        W1 += sites_.op("Sp",n) * row(2) * col(NDX+NDZ+4) * Jx*0.5;

        W1 += sites_.op("Sm",n) * row(2) * col(NDX+NDY+NDZ+4) * Jy*Cplx_i *0.5;
        W1 +=-sites_.op("Sp",n) * row(2) * col(NDX+NDY+NDZ+4) * Jy*Cplx_i *0.5;


        W1 += sites_.op("Cup",n) * row(NDX+NDY+NDZ+6) * col(1);			//Cup
        W1 += sites_.op("Cup",n) * row(NDX+2*NDY+NDZ+6) * col(1);		//Cup		NN
        W1 += sites_.op("Cdn",n) * row(NDX+2*NDY+NDZ+7) * col(1);		//Cdn
        W1 += sites_.op("Cdn",n) * row(NDX+3*NDY+NDZ+7) * col(1);		//Cdn		NN

        W1 += sites_.op("Cdagup",n) * row(NDX+3*NDY+NDZ+8) * col(1);		//Cup
        W1 += sites_.op("Cdagup",n) * row(NDX+4*NDY+NDZ+8) * col(1);		//Cup		NN
        W1 += sites_.op("Cdagdn",n) * row(NDX+4*NDY+NDZ+9) * col(1);		//Cdn
        W1 += sites_.op("Cdagdn",n) * row(NDX+5*NDY+NDZ+9) * col(1);		//Cdn		NN

        for (int jj = NDX+NDY+NDZ+7; jj <=NDX+2*NDY+NDZ+5; ++jj)                {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+2*NDY+NDZ+8; jj <=NDX+3*NDY+NDZ+6; ++jj)              {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+3*NDY+NDZ+9; jj <=NDX+4*NDY+NDZ+7; ++jj)	        {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+4*NDY+NDZ+10; jj <=NDX+5*NDY+NDZ+8; ++jj)  		{W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}

        W1 += -t*sites_.op("Cdagup",n) * row(2) * col(NDX+2*NDY+NDZ+6);					//Z
        W1 += -t*sites_.op("Cdagdn",n) * row(2) * col(NDX+3*NDY+NDZ+7);					//Z
        W1 += t*sites_.op("Cdagup",n) * row(2) * col(NDX+4*NDY+NDZ+8);					//Z
        W1 += t*sites_.op("Cdagdn",n) * row(2) * col(NDX+5*NDY+NDZ+9);					//Z

        W1 += -t*sites_.op("Cdagup",n) * row(2) * col(NDX+2*NDY+NDZ+5);					//Z
        W1 += -t*sites_.op("Cdagdn",n) * row(2) * col(NDX+3*NDY+NDZ+6);					//Z
        W1 += t*sites_.op("Cdagup",n) * row(2) * col(NDX+4*NDY+NDZ+7);					//Z
        W1 += t*sites_.op("Cdagdn",n) * row(2) * col(NDX+5*NDY+NDZ+8);					//Z


}


        if (n%4==3)
{
        W1 += sites_.op("Id",n) * row(1) * col(1); //ending state
        W1 += sites_.op("Id",n) * row(2) * col(2); //starting state
        W1 += sites_.op("Sz",n) * row(3) * col(1);					//Z
        W1 += sites_.op("Sz",n) * row(3+NDZ) * col(1);					//Z

        W1 += sites_.op("Sp",n) * row(4+NDZ) * col(1)*0.5;	        		//X
        W1 += sites_.op("Sm",n) * row(4+NDZ) * col(1)*0.5;	    			//X

        W1 += sites_.op("Sp",n) * row(NDX+NDZ+4) * col(1)*0.5;				//X
        W1 += sites_.op("Sm",n) * row(NDX+NDZ+4) * col(1)*0.5;				//X

        W1 += sites_.op("Sm",n) * row(NDX+NDZ+5) * col(1)*0.5*Cplx_i;			//Y
        W1 +=-sites_.op("Sp",n) * row(NDX+NDZ+5) * col(1)*0.5*Cplx_i;			//Y
        W1 += sites_.op("Sm",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5*Cplx_i;		//Y
        W1 +=-sites_.op("Sp",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5*Cplx_i;		//Y

        for (int jj = 4; jj <=NDZ+2; ++jj)			{W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 5+NDZ; jj <=NDX+NDZ+3; ++jj)		{W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 6+NDX+NDZ; jj <=NDX+NDY+NDZ+4; ++jj)	{W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}

        W1 += sites_.op("Sz",n) * row(2) * col(3+NDZ) * Jz;
        W1 += sites_.op("Sm",n) * row(2) * col(NDX+NDZ+3) * Jx *0.5;
        W1 += sites_.op("Sp",n) * row(2) * col(NDX+NDZ+3) * Jx *0.5


        W1 += sites_.op("Cup",n) * row(NDX+NDY+NDZ+6) * col(1);			//Cup
        W1 += sites_.op("Cup",n) * row(NDX+2*NDY+NDZ+6) * col(1);		//Cup		NN
        W1 += sites_.op("Cdn",n) * row(NDX+2*NDY+NDZ+7) * col(1);		//Cdn
        W1 += sites_.op("Cdn",n) * row(NDX+3*NDY+NDZ+7) * col(1);		//Cdn		NN

        W1 += sites_.op("Cdagup",n) * row(NDX+3*NDY+NDZ+8) * col(1);		//Cup
        W1 += sites_.op("Cdagup",n) * row(NDX+4*NDY+NDZ+8) * col(1);		//Cup		NN
        W1 += sites_.op("Cdagdn",n) * row(NDX+4*NDY+NDZ+9) * col(1);		//Cdn
        W1 += sites_.op("Cdagdn",n) * row(NDX+5*NDY+NDZ+9) * col(1);		//Cdn		NN

        for (int jj = NDX+NDY+NDZ+7; jj <=NDX+2*NDY+NDZ+5; ++jj)                {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+2*NDY+NDZ+8; jj <=NDX+3*NDY+NDZ+6; ++jj)              {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+3*NDY+NDZ+9; jj <=NDX+4*NDY+NDZ+7; ++jj)	        {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+4*NDY+NDZ+10; jj <=NDX+5*NDY+NDZ+8; ++jj)  		{W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}

        W1 += -t*sites_.op("Cdagup",n) * row(2) * col(NDX+2*NDY+NDZ+6);					//Z
        W1 += -t*sites_.op("Cdagdn",n) * row(2) * col(NDX+3*NDY+NDZ+7);					//Z
        W1 += t*sites_.op("Cdagup",n) * row(2) * col(NDX+4*NDY+NDZ+8);					//Z
        W1 += t*sites_.op("Cdagdn",n) * row(2) * col(NDX+5*NDY+NDZ+9);					//Z

        W1 += -t*sites_.op("Cdagup",n) * row(2) * col(NDX+2*NDY+NDZ+3);					//Z
        W1 += -t*sites_.op("Cdagdn",n) * row(2) * col(NDX+3*NDY+NDZ+4);					//Z
        W1 += t*sites_.op("Cdagup",n) * row(2) * col(NDX+4*NDY+NDZ+5);					//Z
        W1 += t*sites_.op("Cdagdn",n) * row(2) * col(NDX+5*NDY+NDZ+6);					//Z


}


        if (n%4==0 && (n%(4*Nx)!=0))
{
        W1 += sites_.op("Id",n) * row(1) * col(1); //ending state
        W1 += sites_.op("Id",n) * row(2) * col(2); //starting state
        W1 += sites_.op("Sz",n) * row(3) * col(1);					//Z
        W1 += sites_.op("Sz",n) * row(3+NDZ) * col(1);					//Z
        W1 += sites_.op("Sp",n) * row(4+NDZ) * col(1)*0.5;	        		//X
        W1 += sites_.op("Sm",n) * row(4+NDZ) * col(1)*0.5;	    			//X
        W1 += sites_.op("Sp",n) * row(NDX+NDZ+4) * col(1)*0.5;				//X
        W1 += sites_.op("Sm",n) * row(NDX+NDZ+4) * col(1)*0.5;				//X
        W1 += sites_.op("Sm",n) * row(NDX+NDZ+5) * col(1)*0.5*Cplx_i;			//Y
        W1 +=-sites_.op("Sp",n) * row(NDX+NDZ+5) * col(1)*0.5*Cplx_i;			//Y
        W1 += sites_.op("Sm",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5*Cplx_i;		//Y
        W1 +=-sites_.op("Sp",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5*Cplx_i;		//Y

        for (int jj = 4; jj <=NDZ+2; ++jj)          {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 5+NDZ; jj <=NDX+NDZ+3; ++jj)          {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 6+NDX+NDZ; jj <=NDX+NDY+NDZ+4; ++jj)  {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}

        W1 += sites_.op("Sm",n) * row(2) * col(NDX+NDY+NDZ+5) * Jy*Cplx_i*0.5;
        W1 +=-sites_.op("Sp",n) * row(2) * col(NDX+NDY+NDZ+5) * Jy*Cplx_i*0.5;

        W1 += sites_.op("Cup",n) * row(NDX+NDY+NDZ+6) * col(1);			//Cup
        W1 += sites_.op("Cup",n) * row(NDX+2*NDY+NDZ+6) * col(1);		//Cup		NN
        W1 += sites_.op("Cdn",n) * row(NDX+2*NDY+NDZ+7) * col(1);		//Cdn
        W1 += sites_.op("Cdn",n) * row(NDX+3*NDY+NDZ+7) * col(1);		//Cdn		NN

        W1 += sites_.op("Cdagup",n) * row(NDX+3*NDY+NDZ+8) * col(1);		//Cup
        W1 += sites_.op("Cdagup",n) * row(NDX+4*NDY+NDZ+8) * col(1);		//Cup		NN
        W1 += sites_.op("Cdagdn",n) * row(NDX+4*NDY+NDZ+9) * col(1);		//Cdn
        W1 += sites_.op("Cdagdn",n) * row(NDX+5*NDY+NDZ+9) * col(1);		//Cdn		NN

        for (int jj = NDX+NDY+NDZ+7; jj <=NDX+2*NDY+NDZ+5; ++jj)                {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+2*NDY+NDZ+8; jj <=NDX+3*NDY+NDZ+6; ++jj)              {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+3*NDY+NDZ+9; jj <=NDX+4*NDY+NDZ+7; ++jj)	        {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+4*NDY+NDZ+10; jj <=NDX+5*NDY+NDZ+8; ++jj)  		{W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}

        W1 += -t*sites_.op("Cdagup",n) * row(2) * col(NDX+2*NDY+NDZ+6);					//Z
        W1 += -t*sites_.op("Cdagdn",n) * row(2) * col(NDX+3*NDY+NDZ+7);					//Z
        W1 += t*sites_.op("Cdagup",n) * row(2) * col(NDX+4*NDY+NDZ+8);					//Z
        W1 += t*sites_.op("Cdagdn",n) * row(2) * col(NDX+5*NDY+NDZ+9);					//Z




	}

        if ((n%(4*Nx)==0))
{
        W1 += sites_.op("Id",n) * row(1) * col(1); //ending state
        W1 += sites_.op("Id",n) * row(2) * col(2); //starting state
        W1 += sites_.op("Sz",n) * row(3) * col(1);					//Z
        W1 += sites_.op("Sz",n) * row(3+NDZ) * col(1);					//Z
        W1 += sites_.op("Sp",n) * row(4+NDZ) * col(1)*0.5;	        		//X
        W1 += sites_.op("Sm",n) * row(4+NDZ) * col(1)*0.5;	    			//X
        W1 += sites_.op("Sp",n) * row(NDX+NDZ+4) * col(1)*0.5;				//X
        W1 += sites_.op("Sm",n) * row(NDX+NDZ+4) * col(1)*0.5;				//X
        W1 += sites_.op("Sm",n) * row(NDX+NDZ+5) * col(1)*0.5*Cplx_i;			//Y
        W1 +=-sites_.op("Sp",n) * row(NDX+NDZ+5) * col(1)*0.5*Cplx_i;			//Y
        W1 += sites_.op("Sm",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5*Cplx_i;		//Y
        W1 +=-sites_.op("Sp",n) * row(NDX+NDY+NDZ+5) * col(1)*0.5*Cplx_i;		//Y

        for (int jj = 4; jj <=NDZ+2; ++jj)          {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 5+NDZ; jj <=NDX+NDZ+3; ++jj)          {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = 6+NDX+NDZ; jj <=NDX+NDY+NDZ+4; ++jj)  {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}

        W1 += sites_.op("Sz",n) * row(2) * col(2+NDZ) * Jz;
        W1 += sites_.op("Sm",n) * row(2) * col(NDX+NDY+NDZ+5) * Jy*Cplx_i *0.5;
        W1 +=-sites_.op("Sp",n) * row(2) * col(NDX+NDY+NDZ+5) * Jy*Cplx_i *0.5;


        W1 += sites_.op("Cup",n) * row(NDX+NDY+NDZ+6) * col(1);			//Cup
        W1 += sites_.op("Cup",n) * row(NDX+2*NDY+NDZ+6) * col(1);		//Cup		NN
        W1 += sites_.op("Cdn",n) * row(NDX+2*NDY+NDZ+7) * col(1);		//Cdn
        W1 += sites_.op("Cdn",n) * row(NDX+3*NDY+NDZ+7) * col(1);		//Cdn		NN

        W1 += sites_.op("Cdagup",n) * row(NDX+3*NDY+NDZ+8) * col(1);		//Cup
        W1 += sites_.op("Cdagup",n) * row(NDX+4*NDY+NDZ+8) * col(1);		//Cup		NN
        W1 += sites_.op("Cdagdn",n) * row(NDX+4*NDY+NDZ+9) * col(1);		//Cdn
        W1 += sites_.op("Cdagdn",n) * row(NDX+5*NDY+NDZ+9) * col(1);		//Cdn		NN

        for (int jj = NDX+NDY+NDZ+7; jj <=NDX+2*NDY+NDZ+5; ++jj)                {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+2*NDY+NDZ+8; jj <=NDX+3*NDY+NDZ+6; ++jj)              {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+3*NDY+NDZ+9; jj <=NDX+4*NDY+NDZ+7; ++jj)	        {W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}
        for (int jj = NDX+4*NDY+NDZ+10; jj <=NDX+5*NDY+NDZ+8; ++jj)  		{W1 += sites_.op("Id",n) * row(jj) * col(jj-1);}

        W1 += -t*sites_.op("Cdagup",n) * row(2) * col(NDX+2*NDY+NDZ+6);					//Z
        W1 += -t*sites_.op("Cdagdn",n) * row(2) * col(NDX+3*NDY+NDZ+7);					//Z

        W1 += t*sites_.op("Cdagup",n) * row(2) * col(NDX+4*NDY+NDZ+8);					//Z
        W1 += t*sites_.op("Cdagdn",n) * row(2) * col(NDX+5*NDY+NDZ+9);					//Z



        W1 += -t*sites_.op("Cdagup",n) * row(2) * col(NDX+2*NDY+NDZ+3);					//Z
        W1 += -t*sites_.op("Cdagdn",n) * row(2) * col(NDX+3*NDY+NDZ+4);					//Z
        W1 += t*sites_.op("Cdagup",n) * row(2) * col(NDX+4*NDY+NDZ+5);					//Z
        W1 += t*sites_.op("Cdagdn",n) * row(2) * col(NDX+5*NDY+NDZ+6);					//Z


	}


}
    auto LH1 = setElt(links1.at(0)(2));
    auto RH1 = setElt(dag(last1)(1));

    HtJ.Anc(0) = LH1;
    HtJ.Anc(N_+1) = RH1;

    auto res = idmrg(psi,HtJ,sweeps,{"OutputLevel",1});



    return 0;
    }

