# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <algorithm>    

using namespace std;

int main()
{
	//Domain
	
	double a,b,rho1,rho2,mu,un,us,ve,vw,time,rad,xd,yd,gx,gy,lrg;
	a=1;
	b=1;
	rho1=1;
	rho2=2;
	mu=0.01;
	un=0;
	us=0;
	ve=0;
	vw=0;
	time=1.5;
	rad=0.15;
	xd=0.5;
	yd=0.7;
	gx=0;
	gy=100;
	lrg=1000;

	//Numerical
	
	int nx,ny,nt;
	double dt,maxiter,maxeps,beta,max;
	nx=32;
	ny=32;
	nt=100;
	dt=0.00125;
	maxiter=200;
	maxeps=0.001;
	beta=1.2;
	max=0;
	
	//Initialization
	
	double **u,**v,**p,**ut,**uu,**vt,**vv,**tmp1,**tmp2,**rt,**p_check,**err,**rho;
	u=new double *[nx+1];
	v=new double *[nx+2];
	p=new double *[nx+2];
	ut=new double *[nx+1];
	uu=new double *[nx+1];
	vt=new double *[nx+2];
	vv=new double *[nx+1];
	tmp1=new double *[nx+2];	
	tmp2=new double *[nx+2];	
	rt=new double *[nx+2];
	p_check=new double *[nx+2];
	err=new double *[nx+2];
	rho=new double *[nx+2];

	for (int i=0;i<=nx+1;i++)
	{
		v[i]=new double[ny+1];
		p[i]=new double[ny+2];
		p_check[i]=new double [ny+2];
		vt[i]=new double[ny+1];
		tmp1[i]=new double[ny+2];
		tmp2[i]=new double[ny+2];	
		rt[i]=new double [ny+2];
		err[i]=new double [ny+2];
		rho[i]=new double [ny+2];
	}

	for (int i=0;i<=nx;i++)
	{
		u[i]=new double[ny+2];
		ut[i]=new double[ny+2];
		uu[i]=new double[ny+1];
		vv[i]=new double[ny+1];
	}

	//Grid
	
	double *x,*y,dx,dy;
	x=new double [nx+1];
	y=new double [ny+1];
	dx=a/nx;
	dy=b/ny;
	for (int i=0;i<=nx;i++)
	{
		x[i]=(i)*dx;

	}
	for (int j=0;j<=ny;j++)
	{
		y[j]=(j)*dy;

	}
	cout<<dt<<endl;
	cout<<dx<<endl;
	cout<<dy<<endl;
	//Density
	
	double **r;
	r=new double *[nx+2];
	
	for (int i=0;i<=nx+1;i++)
	{
		r[i]=new double[ny+2];
	}
	
	//Initialization 
	
	for (int i=0;i<=nx+1;i++)
	{	for (int j=0;j<=ny+1;j++)
		{
   			if ((pow(x[i]-xd,2)+pow(y[j]-yd,2))<pow(rad,2))
			{
				r[i][j]=rho2;
			}
			else
			{
				r[i][j]=rho1;
			}
				rt[i][j]=r[i][j];
		}
	}
	dt=0.00125;

	//velocity advection
	ofstream output;
	string file="results.txt";
	output.open(file.c_str());
	for(int i=0;i<=nx+1;i++)
		{
			for(int j=0;j<=ny+1;j++)
			{
				output<<r[i][j]<<" ";	
			}	
			output<<endl;
		}	
	for (int  n=0;n<1000;n++)
	{

		cout<<"Iteration "<<n<<endl;
		for(int i=0;i<=nx;i++)
			{
				u[i][0]=2*us-u[i][1];
				u[i][ny+1]=2*un-u[i][ny];
			}
		for(int j=0;j<=ny;j++)
			{
				v[0][j]=2*vw-v[1][j];
				v[nx+1][j]=2*ve-v[nx][j];
			}
		
		cout<<"computing x velocity contribution"<<endl;
		
		for(int i=1;i<nx;i++)
		{
			for(int j=1;j<ny+1;j++)
			{
			
				ut[i][j]=u[i][j]+dt*(-0.25*((pow(u[i+1][j]+u[i][j],2)-pow(u[i][j]+u[i-1][j],2))/dx+((u[i][j+1]+u[i][j])*(v[i+1][j]+v[i][j])-(u[i][j]+u[i][j-1])*(v[i+1][j-1]+v[i][j-1]))/dy)+mu/(0.5*(r[i+1][j]+r[i][j]))*((u[i+1][j]-2*u[i][j]+u[i-1][j])/pow(dx,2)+(u[i][j+1]-2*u[i][j]+u[i][j-1])/pow(dy,2))+gx);

			}

		}

		cout<<"computing y velocity contribution"<<endl;
		
		for(int i=1;i<nx+1;i++)
		{
			for(int j=1;j<ny;j++)
			{
				vt[i][j]=v[i][j]+dt*(-0.25*(((u[i][j+1]+u[i][j])*(v[i+1][j]+v[i][j])-(u[i-1][j+1]+u[i-1][j])*(v[i][j]+v[i-1][j]))/dx+(pow(v[i][j+1]+v[i][j],2)-pow(v[i][j]+v[i][j-1],2))/dy)+mu/(0.5*(r[i][j+1]+r[i][j]))*((v[i+1][j]-2*v[i][j]+v[i-1][j])/pow(dx,2)+(v[i][j+1]-2*v[i][j]+v[i][j-1])/pow(dy,2))+gy);						
			}
		}
	for (int i=0;i<=nx+1;i++)
	{
		for (int j=0;j<=ny+1;j++)
		{
		rt[i][j]=r[i][j];
		}
	}
		for(int i=0;i<=nx+1;i++)
			{
				rt[i][0]=lrg;
				rt[i][ny+1]=lrg;
			}
		for(int j=0;j<=ny+1;j++)
			{
				rt[0][j]=lrg;
				rt[nx+1][j]=lrg;
			}
		
		cout<<"computing tmp1"<<endl;

	
		for(int i=1;i<nx+1;i++)
		{
			for(int j=1;j<ny+1;j++)
			{

				tmp1[i][j]=0.5*((ut[i][j]-ut[i-1][j])/dx+(vt[i][j]-vt[i][j-1])/dy)/dt;
				//tmp2[i][j]=pow((1/dx)*(1/(dx*(rt[i+1][j]+rt[i][j]))+1/(dx*(rt[i-1][j]+rt[i][j])) )+(1/dy)*(1/(dy*(rt[i][j+1]+rt[i][j]))+1/(dy*(rt[i][j-1]+rt[i][j]))),-1);
				tmp2[i][j]=1.0/( (1./dx)*( 1./(dx*(rt[i+1][j]+rt[i][j]))+1./(dx*(rt[i-1][j]+rt[i][j])) )+(1./dy)*(1./(dy*(rt[i][j+1]+rt[i][j]))+1./(dy*(rt[i][j-1]+rt[i][j])) ) );
				//cout<<tmp1[i][j]<<" ";

			}	//cout<<endl;
		}
		//pressure solution SOR method	

		for (int n=0;n<=maxiter;n++)
		{
	
		max=0;

		for(int i=0;i<=nx+1;i++)
		{
			for(int j=0;j<=ny+1;j++)
			{
				p_check[i][j]=p[i][j];
			}
		}
		for(int i=1;i<nx+1;i++)
		{
			for(int j=1;j<ny+1;j++)
			{
				p[i][j]=(1-beta)*p[i][j]+beta*tmp2[i][j]*((1/dx)*(p[i+1][j]/(dx*(rt[i+1][j]+rt[i][j]))+p[i-1][j]/(dx*(rt[i-1][j]+rt[i][j])))+(1/dy)*( p[i][j+1]/(dy*(rt[i][j+1]+rt[i][j]))+p[i][j-1]/(dy*(rt[i][j-1]+rt[i][j])))- tmp1[i][j]);
			}
		}
		
		for(int i=0;i<=nx+1;i++)
		{
			for(int j=0;j<=ny+1;j++)
			{
				err[i][j]=abs(p_check[i][j]-p[i][j]);
			}
			
		}
		for(int i=0;i<=nx+1;i++)
		{
			for(int j=0;j<=ny+1;j++)
			{
				if(err[i][j]>max)
				{
					max=err[i][j];
				}
			}
			
		}
		if(max<maxeps) {break;}
		}
		for(int i=1;i<nx;i++)
		{
			for(int j=1;j<=ny;j++) 
			{
			u[i][j]=ut[i][j]-dt*(2/dx)*(p[i+1][j]-p[i][j])/(r[i+1][j]+r[i][j]);
			}
			
		}
		for(int i=1;i<=nx;i++)
		{	
			for(int j=1;j<ny;j++) 
			{
			v[i][j]=vt[i][j]-dt*(2/dy)*(p[i][j+1]-p[i][j])/(r[i][j+1]+r[i][j]);		
			}
			
		}
	for (int i=0;i<=nx+1;i++)
	{
		for (int j=0;j<=ny+1;j++)
		{
		rho[i][j]=r[i][j];
		}
	}
	for(int i=1;i<nx+1;i++)
		{
			for(int j=1;j<ny+1;j++)
			{
				
				r[i][j]=rho[i][j]-(0.5*dt/dx)*(u[i][j]*(rho[i+1][j]
+rho[i][j])-u[i-1][j]*(rho[i-1][j]+rho[i][j]))-(0.5* dt/dy)*(v[i][j]*(rho[i][j+1]+rho[i][j])-v[i][j-1]*(rho[i][j-1]+rho[i][j]) )+(mu*dt/dx/dx)*(rho[i+1][j]-2.0*rho[i][j]+rho[i-1][j])+(mu*dt/dy/dy)*(rho[i][j+1]-2.0*rho[i][j]+rho[i][j-1]);
				
			}
		}
		//output<<"ITER "<<n<<endl;
		for(int i=0;i<=nx+1;i++)
		{
			for(int j=0;j<=ny+1;j++)
			{
				output<<r[i][j]<<" ";	
			}	
			output<<endl;
		}	
	}
			output.close();
return 0;
}
