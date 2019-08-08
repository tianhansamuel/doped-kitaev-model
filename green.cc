#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <math.h>       /* exp */
using namespace std;
using namespace itensor;

int main()
    {
int Nx;

int Ny;

double t;

double U;

double ttotal;

int jp;

double energy;




cout << "Hopping amplitude t";
cin >> t;
cout << "Number of X unit cells Nx";
cin >> Nx;
cout << "Number of Y unit cells Ny";
cin >> Ny;
cout << "The Hubbard interaction U";
cin >> U;
cout << "The initial position";
cin >> jp;
cout << "Total time evolution";
cin >> ttotal;

Hubbard sites;

readFromFile("sites_file_1",sites);

ifstream EE;

EE.open("energy_Jx_4_Jy_4_Jz_4_Nx_1.dat");

EE>>energy;

int N_=4*Nx;

ofstream myfile;
myfile.open ("density_t_"+ std::to_string(t) +"_U_"+std::to_string(U)+"_Nx_"+std::to_string(Nx)+"_Ny_"+std::to_string(Ny)+".dat");


MPS psi(sites);
readFromFile("psi_Jx_4_Jy_4_Jz_4_Nx_1",psi);

auto psi_i = vector<MPS>(N_+1);

for (int ii=1; ii<=N_;++ii)
{
MPS psi1(sites);
readFromFile("psi_Jx_4_Jy_4_Jz_4_Nx_1_Hole_"+std::to_string(ii),psi1);
psi_i.at(ii)=psi1;
}

auto sites1=Hubbard(Ny*N_);

auto tau=0.02;

auto args1 = Args("Cutoff=",1e-10,"Maxm=",2000);

auto psi_hh=psi_i.at(jp);		//hole MPS

auto psi_start=MPS(sites1);
auto psi_final=MPS(sites1);

// psi		half-filled MPS

///////////////////////////////////////////////////////////////////////////////////////


auto psi_fh=MPS(sites);
auto ampo0 = AutoMPO(sites);
ampo0 += "Ntot", 1;
auto Adnop0=MPO(ampo0);
fitApplyMPO(psi,Adnop0,psi_fh,args1);

normalize(psi_fh);

/////////////////////////////////////////////////////////////////////////////////////////

//	regauging the tensor


for (int iy=0; iy<Ny;++iy)
{
	for (int ix=1;ix<=N_;++ix)
	{
	int is=ix+iy*N_;
	println(is);
	auto II = ITensor(dag(sites(ix)),sites1(is));
	II.set(dag(sites(ix,"Up")),sites1(is,"Up"), 1);
	II.set(dag(sites(ix,"Dn")),sites1(is,"Dn"), 1);
	II.set(dag(sites(ix,"Emp")),sites1(is,"Emp"), 1);

	if (iy==(Ny-1)/2)
	{
	psi_start.setA(is,psi_hh.A(ix));
	psi_start.Anc(is) *=II;
	}
	else	
	{
	psi_start.setA(is,psi_fh.A(ix));
	psi_start.Anc(is) *=II;
	}
	psi_final.setA(is,psi_fh.A(ix));
	psi_final.Anc(is) *=II;
	}
}


/////////////////////////////////////////////////////////////////////////////////////////


auto psi_f = vector<MPS>(N_*Ny+1);

for (int is=1; is<=N_*Ny;++is)
{
	auto psi_fh=MPS(sites1);
	auto ampoj = AutoMPO(sites1);
	ampoj += "Cdn", is;
	auto Adnopj=MPO(ampoj);
	fitApplyMPO(psi_final,Adnopj,psi_fh,args1);
	psi_f.at(is)=psi_fh;
}



auto nt = int(ttotal/tau+(1e-9*(ttotal/tau)));

for (int it=0; it<nt; ++it)
{
	double tt=it*tau;
	println(tt);
for(int jjy = 0; jjy <Ny; ++jjy)
{
	auto ampo1 = AutoMPO(sites1);
	ampo1 +=   0.5*U,  "Nupdn",Nx+jjy*N_;
	ampo1 +=   -t,  "Cdagup",Nx+jjy*N_,"Cup",jjy*N_+1;
	ampo1 +=   -t,  "Cdagdn",Nx+jjy*N_,"Cdn",jjy*N_+1;
	ampo1 +=   -t,  "Cdagup",jjy*N_+1,"Cup",Nx+jjy*N_;
	ampo1 +=   -t,  "Cdagdn",jjy*N_+1,"Cdn",Nx+jjy*N_;
	auto expH1 = toExpH<ITensor>(ampo1,tau*Cplx_i);
	fitApplyMPO(psi_start,expH1,psi_start,args1);
	normalize(psi_start);

	auto ampo1b = AutoMPO(sites1);
	ampo1b +=   0.5*U,  "Nupdn",2*Nx+jjy*N_;
	ampo1b +=   -t,  "Cdagup",2*Nx+jjy*N_,"Cup",jjy*N_+1+Nx;
	ampo1b +=   -t,  "Cdagdn",2*Nx+jjy*N_,"Cdn",jjy*N_+1+Nx;
	ampo1b +=   -t,  "Cdagup",jjy*N_+1+Nx,"Cup",2*Nx+jjy*N_;
	ampo1b +=   -t,  "Cdagdn",jjy*N_+1+Nx,"Cdn",2*Nx+jjy*N_;
	auto expH1b = toExpH<ITensor>(ampo1b,tau*Cplx_i);
	fitApplyMPO(psi_start,expH1b,psi_start,args1);
	normalize(psi_start);

	auto ampo2 = AutoMPO(sites1);
	ampo2 +=   0.5*U,  "Nupdn",3*Nx+jjy*N_;
	ampo2 +=   -t,  "Cdagup",3*Nx+jjy*N_,"Cup",jjy*N_+1+2*Nx;
	ampo2 +=   -t,  "Cdagdn",3*Nx+jjy*N_,"Cdn",jjy*N_+1+2*Nx;
	ampo2 +=   -t,  "Cdagup",jjy*N_+1+2*Nx,"Cup",3*Nx+jjy*N_;
	ampo2 +=   -t,  "Cdagdn",jjy*N_+1+2*Nx,"Cdn",3*Nx+jjy*N_;
	auto expH2 = toExpH<ITensor>(ampo2,tau*Cplx_i);
	fitApplyMPO(psi_start,expH2,psi_start,args1);
	normalize(psi_start);


	auto ampo2b = AutoMPO(sites1);
	ampo2b +=   0.5*U,  "Nupdn",4*Nx+jjy*N_;
	ampo2b +=   -t,  "Cdagup",4*Nx+jjy*N_,"Cup",jjy*N_+1+3*Nx;
	ampo2b +=   -t,  "Cdagdn",4*Nx+jjy*N_,"Cdn",jjy*N_+1+3*Nx;
	ampo2b +=   -t,  "Cdagup",jjy*N_+1+3*Nx,"Cup",4*Nx+jjy*N_;
	ampo2b +=   -t,  "Cdagdn",jjy*N_+1+3*Nx,"Cdn",4*Nx+jjy*N_;
	auto expH2b = toExpH<ITensor>(ampo2b,tau*Cplx_i);
	fitApplyMPO(psi_start,expH2b,psi_start,args1);
	normalize(psi_start);

  for (int jjx=1; jjx<N_; ++jjx)
  {
	auto ampo5 = AutoMPO(sites1);

	ampo5 +=   0.5*U,  "Nupdn",jjx+jjy*N_;

	ampo5 +=   -t,  "Cdagup",jjx+jjy*N_,"Cup",jjx+jjy*N_+1;
	ampo5 +=   -t,  "Cdagdn",jjx+jjy*N_,"Cdn",jjx+jjy*N_+1;
	ampo5 +=   -t,  "Cdagup",jjx+jjy*N_+1,"Cup",jjx+jjy*N_;
	ampo5 +=   -t,  "Cdagdn",jjx+jjy*N_+1,"Cdn",jjx+jjy*N_;
	auto expH5 = toExpH<ITensor>(ampo5,tau*Cplx_i);
	fitApplyMPO(psi_start,expH5,psi_start,args1);
	normalize(psi_start);

	if(jjx>=1 && jjx<Nx)	//ok
	{
	auto ampo6 = AutoMPO(sites1);

	ampo6 +=   0.5*U,  "Nupdn",jjx+jjy*N_;

	ampo6 +=   -t,  "Cdagup",jjx+jjy*N_,"Cup",2*Nx+1-jjx+jjy*N_;
	ampo6 +=   -t,  "Cdagdn",jjx+jjy*N_,"Cdn",2*Nx+1-jjx+jjy*N_;
	ampo6 +=   -t,  "Cdagup",2*Nx-jjx+jjy*N_+1,"Cup",jjx+jjy*N_;
	ampo6 +=   -t,  "Cdagdn",2*Nx-jjx+jjy*N_+1,"Cdn",jjx+jjy*N_;
	auto expH6 = toExpH<ITensor>(ampo6,tau*Cplx_i);
	fitApplyMPO(psi_start,expH6,psi_start,args1);
	normalize(psi_start);
	}

	if (jjx>=1+Nx && jjx<2*Nx)	//ok
	{
	auto ampo7 = AutoMPO(sites1);

	ampo7 +=   0.5*U,  "Nupdn",jjx+jjy*N_;

	ampo7 +=   -t,  "Cdagup",jjx+jjy*N_,"Cup",4*Nx+1-jjx+jjy*N_;
	ampo7 +=   -t,  "Cdagdn",jjx+jjy*N_,"Cdn",4*Nx+1-jjx+jjy*N_;
	ampo7 +=   -t,  "Cdagup",4*Nx+1-jjx+jjy*N_,"Cup",jjx+jjy*N_;
	ampo7 +=   -t,  "Cdagdn",4*Nx+1-jjx+jjy*N_,"Cdn",jjx+jjy*N_;
	auto expH7 = toExpH<ITensor>(ampo7,tau*Cplx_i);
	fitApplyMPO(psi_start,expH7,psi_start,args1);
	normalize(psi_start);
   }


	if (jjx>=1+2*Nx && jjx<3*Nx)	//ok
	{
	auto ampo8 = AutoMPO(sites1);

	ampo8 +=   0.5*U,  "Nupdn",jjx+jjy*N_;

	ampo8 +=   -t,  "Cdagup",jjx+jjy*N_,"Cup",6*Nx+1-jjx+jjy*N_;
	ampo8 +=   -t,  "Cdagdn",jjx+jjy*N_,"Cdn",6*Nx+1-jjx+jjy*N_;
	ampo8 +=   -t,  "Cdagup",6*Nx+1-jjx+jjy*N_,"Cup",jjx+jjy*N_;
	ampo8 +=   -t,  "Cdagdn",6*Nx+1-jjx+jjy*N_,"Cdn",jjx+jjy*N_;
	auto expH8 = toExpH<ITensor>(ampo8,tau*Cplx_i);
	fitApplyMPO(psi_start,expH8,psi_start,args1);
	normalize(psi_start);
   }

   if (jjx>=3*Nx+1 && jjx<4*Nx)	//ok
  {
	auto ampo4 = AutoMPO(sites1);

	ampo4 +=   0.5*U,  "Nupdn",jjx+jjy*N_;

	ampo4 +=   -t,  "Cdagup",jjx+jjy*N_,"Cup",4*Nx-jjx+1+(jjy+1)%Ny*N_;
	ampo4 +=   -t,  "Cdagdn",jjx+jjy*N_,"Cdn",4*Nx-jjx+1+(jjy+1)%Ny*N_;
	ampo4 +=   -t,  "Cdagup",4*Nx+1-jjx+(jjy+1)%Ny*N_,"Cup",jjx+jjy*N_;
	ampo4 +=   -t,  "Cdagdn",4*Nx+1-jjx+(jjy+1)%Ny*N_,"Cdn",jjx+jjy*N_;
	auto expH4 = toExpH<ITensor>(ampo4,tau*Cplx_i);
	fitApplyMPO(psi_start,expH4,psi_start,args1);
	normalize(psi_start);
   }


}
}


for (int is=1; is<=N_;++is)
{

	auto psi_ff=psi_f.at(is);
        auto GR1 = overlapC(psi_ff, psi_start)*exp (Cplx_i*tau*energy);
        myfile << tt<<"\t";
        myfile << is<<"\t";
        myfile << GR1.real()<<"\t";
        myfile << GR1.imag()<<"\t";
        myfile << abs(GR1)<<"\n";
        println(abs(GR1));
}


}

    return 0;
    }
