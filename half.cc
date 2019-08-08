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

double U;

int MSW;

int n_s;

cout << "The hopping amplitude";
cin >> t;
cout << "The number of plaquettes in the X direction";
cin >> Nx;
cout << "The Hubbard interaction U";
cin >> U;
cout << "The number of extra holes n_s";
cin >> n_s;
cout << "The number of sweep MSW";
cin >> MSW;
// read in the parameters


ofstream myfile;
myfile.open ("hopping_t_"+std::to_string(t)+"_Hubbard_interaction_U_"+std::to_string(U)+"_Nx_"+std::to_string(Nx)+"_Ns_"+std::to_string(n_s)+".dat");



int N_=4*Nx;
        
int NDD=2*Nx-1;

auto sites_ = Hubbard(N_);

auto    H = MPO(sites_);


std::vector<Index> links(N_+1); // This is the vector of bond dimension

std::vector<Index> links1(N_+1); // This is the vector of bond dimension

    for(int l = 0; l <= N_; ++l) 
        {
        links.at(l) = Index(nameint("hl",l),6+4*NDD);
        }

    Index const& last = (links.at(0));

    for(int n = 1; n <= N_; ++n)
        {

        auto& W = H.Anc(n);
        auto row = dag(links.at(n-1));
        auto col = (n==N_ ? last : links.at(n));

        W = ITensor(dag(sites_(n)),prime(sites_(n)),row,col);

        W += sites_.op("Id",n) * row(1) * col(1); //starting state

        W += sites_.op("Adagup",n) * row(2) * col(1);	 
        W += sites_.op("Adagup",n) * row(2+NDD) * col(1);

        W += sites_.op("Aup",n) * row(3+NDD) * col(1);	
        W +=-sites_.op("Aup",n) * row(3+2*NDD) * col(1);

        W += sites_.op("Adagdn",n) * row(4+2*NDD) * col(1);
        W +=-sites_.op("Adagdn",n) * row(4+3*NDD) * col(1);

        W += sites_.op("Adn",n) * row(5+3*NDD) * col(1);
        W += sites_.op("Adn",n) * row(5+4*NDD) * col(1);


        W += sites_.op("Id",n) * row(6+4*NDD) * col(6+4*NDD);	//starting state

        W += U*sites_.op("Nupdn",n) * row(6+4*NDD) * col(1);	//Hubbard interaction

        for (int jj1 = 3; jj1 <=NDD+2; ++jj1)		{W += sites_.op("Id",n) * row(jj1) * col(jj1-1); }		//ok
        for (int jj2 = 4+NDD; jj2 <=2*NDD+3; ++jj2)	{W += sites_.op("Id",n) * row(jj2) * col(jj2-1); }


        for (int jj3 = 5+2*NDD; jj3 <=3*NDD+4; ++jj3)	{W += sites_.op("Id",n) * row(jj3) * col(jj3-1); }		//ok
        for (int jj4 = 6+3*NDD; jj4 <=4*NDD+5; ++jj4)	{W += sites_.op("Id",n) * row(jj4) * col(jj4-1); }

        W += t*sites_.op("Aup",n) * col(2+NDD) * row(6+4*NDD) ;        //NN    
        W += -t*sites_.op("Adagup",n) * col(3+2*NDD) * row(6+4*NDD) ;        //NN
        W += t*sites_.op("Adn",n) * col(4+3*NDD) * row(6+4*NDD) ;      //NN    
        W += -t*sites_.op("Adagdn",n) * col(5+4*NDD) * row(6+4*NDD) ;      //NN

	if(n>=1 && n<Nx)
	{
	int NDH=2*Nx+1-2*n;
        W += t*sites_.op("Aup",n) * col(2+NDH) * row(6+4*NDD) ;        //NN    
        W += -t*sites_.op("Adagup",n) * col(3+NDD+NDH) * row(6+4*NDD) ;        //NN
        W += t*sites_.op("Adn",n) * col(4+2*NDD+NDH) * row(6+4*NDD) ;      //NN    
        W += -t*sites_.op("Adagdn",n) * col(5+3*NDD+NDH) * row(6+4*NDD) ;      //NN
	}

	if(n>=2*Nx+1 && n<3*Nx)
	{ 
	int NDH1=2*Nx+1-2*(n-2*Nx);
        W += t*sites_.op("Aup",n) * col(2+NDH1) * row(6+4*NDD) ;        //NN    
        W += -t*sites_.op("Adagup",n) * col(3+NDD+NDH1) * row(6+4*NDD) ;        //NN
        W += t*sites_.op("Adn",n) * col(4+2*NDD+NDH1) * row(6+4*NDD) ;      //NN    
        W += -t*sites_.op("Adagdn",n) * col(5+3*NDD+NDH1) * row(6+4*NDD) ;      //NN
	}

	if(n>=Nx+1 && n<2*Nx)
	{
	int NDH2=2*(2*Nx-n)+1;
        W += t*sites_.op("Aup",n) * col(2+NDH2) * row(6+4*NDD) ;        //NN    
        W += -t*sites_.op("Adagup",n) * col(3+NDD+NDH2) * row(6+4*NDD) ;        //NN
        W += t*sites_.op("Adn",n) * col(4+2*NDD+NDH2) * row(6+4*NDD) ;      //NN    
        W += -t*sites_.op("Adagdn",n) * col(5+3*NDD+NDH2) * row(6+4*NDD) ;      //NN
	}

	if(n>=3*Nx+1 && n<4*Nx)
	{
	int NDH3=2*(4*Nx-n)+1;
        W += t*sites_.op("Aup",n) * col(2+NDH3) * row(6+4*NDD) ;        //NN    
        W += -t*sites_.op("Adagup",n) * col(3+NDD+NDH3) * row(6+4*NDD) ;        //NN
        W += t*sites_.op("Adn",n) * col(4+2*NDD+NDH3) * row(6+4*NDD) ;      //NN    
        W += -t*sites_.op("Adagdn",n) * col(5+3*NDD+NDH3) * row(6+4*NDD) ;      //NN
	}
	if(n%(2*Nx)==1)
	{
	int NDH4=Nx-1;
        W += t*sites_.op("Aup",n) * col(2+NDH4) * row(6+4*NDD) ;        //NN    
        W += -t*sites_.op("Adagup",n) * col(3+NDD+NDH4) * row(6+4*NDD) ;        //NN
        W += t*sites_.op("Adn",n) * col(4+2*NDD+NDH4) * row(6+4*NDD) ;      //NN    
        W += -t*sites_.op("Adagdn",n) * col(5+3*NDD+NDH4) * row(6+4*NDD) ;      //NN
	}
	if(n%(2*Nx)==Nx+1)
	{
	int NDH5=Nx-1;
        W += t*sites_.op("Aup",n) * col(2+NDH5) * row(6+4*NDD) ;        //NN    
        W += -t*sites_.op("Adagup",n) * col(3+NDD+NDH5) * row(6+4*NDD) ;        //NN
        W += t*sites_.op("Adn",n) * col(4+2*NDD+NDH5) * row(6+4*NDD) ;      //NN    
        W += -t*sites_.op("Adagdn",n) * col(5+3*NDD+NDH5) * row(6+4*NDD) ;      //NN
	}

}


    auto sweeps = Sweeps(MSW);
    sweeps.maxm() = 20,40,80,120,200,300;
    sweeps.cutoff() = 1E-15,Args("Repeat",20),1E-15;
    sweeps.niter() = 3,2;

    sweeps.noise() = 1E-7,1E-8,0.0;

    auto state = InitState(sites_);
    for(int i = n_s+1; i <= N_; ++i) 
        {
          state.set(i,"Up");
        }

    for(int i = 1; i <= n_s; ++i) 
        {
          state.set(i,"Emp");
        }

    auto psi = MPS(state);


    auto LH = setElt(links.at(0)(6+4*NDD));
    auto RH = setElt(dag(last)(1));

    H.Anc(0) = LH;
    H.Anc(N_+1) = RH;


auto res = idmrg(psi,H,sweeps,{"OutputLevel",1});

normalize(psi);

double energy=res.energy;


myfile << energy;

// create a hole


auto args = Args("Cutoff=",1e-12,"Maxm=",5000);

auto psi_i = vector<MPS>(N_+1);
for (int js=1; js<=N_; ++js)
{
auto ampo1 = AutoMPO(sites_);
ampo1 += "Cdn", js;
auto Adnop=MPO(ampo1);
auto psi_1=MPS(state);
fitApplyMPO(psi, Adnop,psi_1,args);
psi_i.at(js)=psi_1;
}




writeToFile(format("sites_file_%d",Nx),sites_);

writeToFile(format("psi_t_%d_U_%d_Nx_%d_Ns_%d",t,U,Nx,n_s),psi);     //file name will be psi_100

for (int i1=1; i1<=N_; ++i1)
{
println(i1);
writeToFile(format("psi_t_%d_U_%d_Nx_%d_Ns%d_Hole_%d",t,U,Nx, n_s, i1),psi_i.at(i1));     //file name will be psi_100
}



    return 0;
    }
