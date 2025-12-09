/********************************************************************************/
/***********Elastic Velocity-Stress Finite-difference Modeling with PML**********/
/**************High Order (Highest: sixteen order) Finite-difference*************/
/***************************Written By Zhang Jianming,2022.10.10*****************/
/***************************************CopyRight********************************/


#include<stdio.h>   
#include<math.h>   
#include<iostream>   
#include<fstream>   
#include<iomanip>   
using namespace std;

#define PI 3.1415926   
#define dx 10
#define dz 10   
#define pml 50
#define NX (1000+2*pml)   
#define NZ (1000+2*pml)   
#define NT 10001   
#define dt 0.005   
#define N 6 
float fi=0*(PI/4.0);


int main()
{
    //cout<<"sin(pi/6)="<<sin(PI/6)<<endl;
    float **Txx=new float*[NX],**Txx_x=new float*[NX],**Txx_z=new float*[NX];
    float **Tzz=new float*[NX],**Tzz_x=new float*[NX],**Tzz_z=new float*[NX];
    float **Txz=new float*[NX],**Txz_x=new float*[NX],**Txz_z=new float*[NX];
    float **Vx=new float*[NX],**Vx_x=new float*[NX],**Vx_z=new float*[NX];
    //A wierd question: icpc compeler broke down using the below Vz, why?
    float **Vz=new float*[NX],**Vz_x=new float*[NX],**Vz_z=new float*[NX];
    float **L=new float*[NX],**M=new float*[NX],**e=new float*[NX];
    float **C11=new float*[NX],**C33=new float*[NX],**C44=new float*[NX],**C66=new float*[NX],**C13=new float*[NX],**Rou=new float*[NX];
    float **C_11=new float*[NX],**C_13=new float*[NX],**C_15=new float*[NX],**C_33=new float*[NX],**C_35=new float*[NX],**C_55=new float*[NX],**O=new float*[NX];
    float **Vp=new float*[NX],**Eps=new float*[NX],**Del=new float*[NX],**Vs=new float*[NX],**Gam=new float*[NX],**F=new float*[NX];
    float **data_vx=new float*[NX],**data_vz=new float*[NX];//Observed Vx and Vz components 
    for(int i=0;i<NX;i++)
    {
        Txx[i]=new float[NZ];Txx_x[i]=new float[NZ];Txx_z[i]=new float[NZ];
        Tzz[i]=new float[NZ];Tzz_x[i]=new float[NZ];Tzz_z[i]=new float[NZ];
        Txz[i]=new float[NZ];Txz_x[i]=new float[NZ];Txz_z[i]=new float[NZ];
        Vx[i]=new float[NZ];Vx_x[i]=new float[NZ];Vx_z[i]=new float[NZ];
        Vz[i]=new float[NZ];Vz_x[i]=new float[NZ];Vz_z[i]=new float[NZ];
        L[i]=new float[NZ];M[i]=new float[NZ];e[i]=new float[NZ];
        C11[i]=new float[NZ];C33[i]=new float[NZ];C44[i]=new float[NZ];C66[i]=new float[NZ];C13[i]=new float[NZ];Rou[i]=new float[NZ];
        Vp[i]=new float[NZ];Eps[i]=new float[NZ];Del[i]=new float[NZ];Vs[i]=new float[NZ];Gam[i]=new float[NZ];F[i]=new float[NZ];
        C_11[i]=new float[NZ];C_13[i]=new float[NZ];C_15[i]=new float[NZ];C_33[i]=new float[NZ];C_35[i]=new float[NZ];C_55[i]=new float[NZ];O[i]=new float[NZ];
        for(int j=0;j<NZ;j++)
        {
            Txx[i][j]=0;Txx_x[i][j]=0;Txx_z[i][j]=0;
            Tzz[i][j]=0;Tzz_x[i][j]=0;Tzz_z[i][j]=0;
            Txz[i][j]=0;Txz_x[i][j]=0;Txz_z[i][j]=0;
            Vx[i][j]=0;Vx_x[i][j]=0;Vx_z[i][j]=0;
            Vz[i][j]=0;Vz_x[i][j]=0;Vz_z[i][j]=0;
        }
        data_vx[i]=new float[NT];data_vz[i]=new float[NT];
        for(int j=0;j<NT;j++)
        {
            data_vx[i][j]=0;
            data_vz[i][j]=0;
        }
    }
    
    //FILE *filename;
    //filename=fopen("vp.dat","rb");
    ////for(int i=0;i<NX;i++)
    ////    for(int j=0;j<NZ;j++)
    //for (int i = pml; i < NX-pml; i++)
    //    for (int j = pml; j < NZ-pml; j++)
    //        fread(&Vp[i][j], sizeof(float), 1, filename);
    //fclose(filename);

    //FILE *filenameE;
    //filenameE=fopen("eps.dat","rb");
    ////for(int i=0;i<NX;i++)
    ////    for(int j=0;j<NZ;j++)
    //for (int i = pml; i < NX-pml; i++)
    //    for (int j = pml; j < NZ-pml; j++)
    //        fread(&Eps[i][j], sizeof(float), 1, filenameE);
    //fclose(filenameE);

    //FILE *filenameD;
    //filenameD=fopen("delta.dat","rb");
    ////for(int i=0;i<NX;i++)
    ////    for(int j=0;j<NZ;j++)
    //for (int i = pml; i < NX-pml; i++)
    //    for (int j = pml; j < NZ-pml; j++)
    //        fread(&Del[i][j], sizeof(float), 1, filenameD);
    //fclose(filenameD);


    ////Model extroplate
    //for(int i=0;i<pml;i++)
    //    for(int j=0;j<NZ;j++)
    //    {
    //        Vp[i][j]=Vp[pml][j];
    //        Eps[i][j]=Eps[pml][j];
    //        Del[i][j]=Del[pml][j];
    //    }
    //for(int i=NX-pml;i<NX;i++)
    //    for(int j=0;j<NZ;j++)
    //    {
    //        Vp[i][j]=Vp[NX-pml-1][j];
    //        Eps[i][j]=Eps[NX-pml-1][j];
    //        Del[i][j]=Del[NX-pml-1][j];
    //    }
    //for(int i=0;i<NX;i++)
    //    for(int j=0;j<pml;j++)
    //    {
    //        Vp[i][j]=Vp[i][pml];
    //        Eps[i][j]=Eps[i][pml];
    //        Del[i][j]=Del[i][pml];
    //    }
    //for(int i=0;i<NX;i++)
    //    for(int j=NZ-pml;j<NZ;j++)
    //    {
    //        Vp[i][j]=Vp[i][NZ-pml-1];
    //        Eps[i][j]=Eps[i][NZ-pml-1];
    //        Del[i][j]=Del[i][NZ-pml-1];
    //    }


    int i,j,k,f0=15;float t0=1.2/f0;
    //int SX=NX/2;int SZ=pml+10;
    int SX=NX/2;int SZ=NZ/2;
    //int SX=2*pml;int SZ=2*pml;
    FILE *fp1,*fp11,*fp12,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7;
    float ddx[NX][NZ];float ddz[NX][NZ];
    float R=0.001,Vmax=7500;
    int x,z,plx,plz;
    for(i=0;i<NX;i++)
        for(j=0;j<NZ;j++)
        {
            //Iso
            //L[i][j]=1.19*pow(10, 10);M[i][j]=5.4*pow(10, 9);e[i][j]=2.0*pow(10, 3);
            //C11[i][j]=L[i][j]+2*M[i][j];
            //C33[i][j]=C11[i][j];
            //C44[i][j]=M[i][j];
            //C66[i][j]=M[i][j];
            //C13[i][j]=L[i][j];
            //Rou[i][j]=e[i][j];
            //Ice 
            //C11[i][j]=13.33*pow(10, 9);
            //C33[i][j]=14.28*pow(10, 9);
            //C44[i][j]=3.26*pow(10, 9);
            //C66[i][j]=3.65*pow(10, 9);
            //C13[i][j]=5.08*pow(10, 9);
            //Rou[i][j]=0.91*pow(10, 3);
            //Quartz 
            //C11[i][j]=116.6*pow(10, 9);
            //C33[i][j]=110.40*pow(10, 9);
            //C44[i][j]=36.06*pow(10, 9);
            //C66[i][j]=49.95*pow(10, 9);
            //C13[i][j]=32.8*pow(10, 9);
            //Rou[i][j]=2.65*pow(10, 3);
            //Mantle 
            //C11[i][j]=212.8*pow(10, 9);
            //C33[i][j]=249.18*pow(10, 9);
            //C44[i][j]=70.47*pow(10, 9);
            //C66[i][j]=66.72*pow(10, 9);
            //C13[i][j]=76.24*pow(10, 9);
            //Rou[i][j]=3.22*pow(10, 3);
            //Ref 
            //C11[i][j]=25.5*pow(10, 9);
            //C33[i][j]=18.4*pow(10, 9);
            //C44[i][j]=5.6*pow(10, 9);
            //C66[i][j]=5.6*pow(10, 9);
            //C13[i][j]=10.4*pow(10, 9);
            //Rou[i][j]=2.44*pow(10, 3);

            //Thomsen 
            Vp[i][j]=3000;
            Eps[i][j]*=0;
            Del[i][j]=0.0;
            Vs[i][j]=2000;//Vp[i][j]/sqrt(3);
            Gam[i][j]=0.1;
            if(j>7.0/10.0*NZ)
            //if(j>-200*i/1000.0+800)
            //if(j>-100*cos(2*PI*i/500.0)+700)
            {
                //Vp[i][j]=4000;
                Vs[i][j]=2800;
                //Vs[i][j]=2800;
            }

            if(j>8.0/10.0*NZ)
            {
                Vp[i][j]=4000;
                //Vs[i][j]=2800;
            }
            //if(j>80) {Vp[i][j]=1000;Vs[i][j]=800;}

            Rou[i][j]=3.44*pow(10, 3);
            C11[i][j]=Rou[i][j]*(1.0+2*Eps[i][j])*pow(Vp[i][j],2);
            C33[i][j]=Rou[i][j]*pow(Vp[i][j],2);
            C44[i][j]=Rou[i][j]*pow(Vs[i][j],2);
            C66[i][j]=Rou[i][j]*(1.0+2*Gam[i][j])*pow(Vs[i][j],2);
            F[i][j]=1.0-pow(Vs[i][j],2)/pow(Vp[i][j],2);
            C13[i][j]=Rou[i][j]*pow(Vp[i][j],2)*sqrt(F[i][j]*(F[i][j]+2*Del[i][j]))-Rou[i][j]*pow(Vs[i][j],2);


            O[i][j]=fi;
            C_11[i][j]=C11[i][j]*pow(cos(O[i][j]),4)+C33[i][j]*pow(sin(O[i][j]),4)+(2*C13[i][j]+4*C44[i][j])*pow(sin(O[i][j]),2)*pow(cos(O[i][j]),2);
            C_13[i][j]=(C11[i][j]+C33[i][j]-4*C44[i][j])*pow(sin(O[i][j]),2)*pow(cos(O[i][j]),2)+C13[i][j]*(pow(sin(O[i][j]),4)+pow(cos(O[i][j]),4));
            C_15[i][j]=(C13[i][j]-C11[i][j]+2*C44[i][j])*pow(sin(O[i][j]),1)*pow(cos(O[i][j]),3)-(C13[i][j]-C33[i][j]+2*C44[i][j])*pow(sin(O[i][j]),3)*pow(cos(O[i][j]),1);
            C_33[i][j]=C11[i][j]*pow(sin(O[i][j]),4)+C33[i][j]*pow(cos(O[i][j]),4)+(2*C13[i][j]+4*C44[i][j])*pow(sin(O[i][j]),2)*pow(cos(O[i][j]),2);
            C_35[i][j]=(C13[i][j]-C11[i][j]+2*C44[i][j])*pow(sin(O[i][j]),3)*pow(cos(O[i][j]),1)-(C13[i][j]-C33[i][j]+2*C44[i][j])*pow(sin(O[i][j]),1)*pow(cos(O[i][j]),3);
            C_55[i][j]=(C11[i][j]+C33[i][j]-2*C13[i][j])*pow(sin(O[i][j]),2)*pow(cos(O[i][j]),2)+C44[i][j]*pow(pow(cos(O[i][j]),2)-pow(sin(O[i][j]),2),2);
        }

    //PML
    plx=pml*dx;plz=pml*dz;
    for(i=0;i<NX;i++)
        for(j=0;j<NZ;j++)
        {
            if(i>=0&&i<pml&&j>=0&&j<pml)
            {
                x=pml-i;z=pml-j;
                ddx[i][j]=-log(R)*3*Vmax*x*x/(2*plx*plx);
                ddz[i][j]=-log(R)*3*Vmax*z*z/(2*plz*plz);
            }

            else if(i>=0&&i<pml&&j>NZ-pml&&j<NZ)
            {
                x=pml-i;z=j-(NZ-pml);
                ddx[i][j]=-log(R)*3*Vmax*x*x/(2*plx*plx);
                ddz[i][j]=-log(R)*3*Vmax*z*z/(2*plz*plz);
            }
            else if(i>NX-pml&&i<NX&&j>=0&&j<pml)
            {
                x=i-(NX-pml);z=pml-j;
                ddx[i][j]=-log(R)*3*Vmax*x*x/(2*plx*plx);
                ddz[i][j]=-log(R)*3*Vmax*z*z/(2*plz*plz);
            }
            else if(i>NX-pml&&i<NX&&j>NZ-pml&&j<NZ)
            {
                x=i-(NX-pml);z=j-(NZ-pml);
                ddx[i][j]=-log(R)*3*Vmax*x*x/(2*plx*plx);
                ddz[i][j]=-log(R)*3*Vmax*z*z/(2*plz*plz);
            }
            else if(i>=pml&&i<=NX-pml&&j>=0&&j<pml)
            {
                x=i-pml;z=pml-j;
                ddx[i][j]=0;
                ddz[i][j]=-log(R)*3*Vmax*z*z/(2*plz*plz);
            }
            else if(i>=pml&&i<=NX-pml&&j>NZ-pml&&j<NZ)
            {
                x=i-pml;z=j-(NZ-pml);
                ddx[i][j]=0;
                ddz[i][j]=-log(R)*3*Vmax*z*z/(2*plz*plz);
            }
            else if(i>=0&&i<pml&&j>=pml&&j<=NZ-pml)
            {
                x=pml-i;z=j-pml;
                ddx[i][j]=-log(R)*3*Vmax*x*x/(2*plx*plx);
                ddz[i][j]=0;
            }
            else if(i>NX-pml&&i<NX&&j>=pml&&j<=NZ-pml)
            {
                x=i-(NX-pml);z=j-pml;
                ddx[i][j]=-log(R)*3*Vmax*x*x/(2*plx*plx);
                ddz[i][j]=0;
            }
            else if(i>=pml&&i<=NX-pml&&j>=pml&&j<=NZ-pml)
            {
                x=i-pml;z=j-pml;
                ddx[i][j]=0;
                ddz[i][j]=0;
            }
        }

    float cof[N];
    if(N==1)
    {    cof[0]=1.0;}
    if(N==2)
    {    cof[0]=1.125;cof[1]=-0.041666667;}
    if(N==3)
    {    cof[0]=1.171875;cof[1]=-0.065104167;cof[2]=0.0046875;}
    if(N==4)
    {    cof[0]=1.196289;cof[1]=-0.0797526;cof[2]=0.009570313;cof[3]=-0.0006975447;}
    if(N==5)
    {    cof[0]=1.2112427;cof[1]=-0.08972168;cof[2]=0.013842773;cof[3]=-0.0017656599;cof[4]=0.00011867947;}
    if(N==6)
    {    cof[0]=1.2213364;cof[1]=-0.096931458;cof[2]=0.017447662;cof[3]=-0.0029672895;cof[4]=0.00035900540;cof[5]=-0.000021847812;}
    if(N==7)
    {    cof[0]=1.2286062;cof[1]=-0.10238385;cof[2]=0.02047677;cof[3]=-0.0041789327;cof[4]=0.00068945355;cof[5]=-0.000076922503;cof[6]=0.0000042365148;}
    if(N==8)
    {    cof[0]=1.2340911;cof[1]=-0.10664985;cof[2]=0.023036367;cof[3]=-0.0053423856;cof[4]=0.0010772712;cof[5]=-0.00016641888;cof[6]=0.000017021711;cof[7]=-0.00000085234642;}

    for(k=0;k<NT;k++)
    {
        for(i=N;i<NX-N;i++)
        {
            for(j=N;j<NZ-N;j++)
            {
                //high order TTI
                float pxVx=0;float pzVx=0;float pxVz=0;float pzVz=0;
                for(int in=0;in<N;in++)
                {
                    pxVx+=cof[in]*(Vx[i+in+1][j]-Vx[i-in][j]);
                    pxVz+=cof[in]*(Vz[i+in][j]-Vz[i-in-1][j]);
                    pzVx+=cof[in]*(Vx[i][j+in+1]-Vx[i][j-in]);
                    pzVz+=cof[in]*(Vz[i][j+in]-Vz[i][j-in-1]);
                }
                Txx_x[i][j]=1.0/(1+0.5*dt*ddx[i][j])*((1-0.5*ddx[i][j]*dt)*Txx_x[i][j]+C_11[i][j]*dt/dx*pxVx+C_15[i][j]*dt/dx*pxVz);
                Txx_z[i][j]=1.0/(1+0.5*dt*ddz[i][j])*((1-0.5*ddz[i][j]*dt)*Txx_z[i][j]+C_13[i][j]*dt/dz*pzVz+C_15[i][j]*dt/dz*pzVx);
                Tzz_x[i][j]=1.0/(1+0.5*dt*ddx[i][j])*((1-0.5*ddx[i][j]*dt)*Tzz_x[i][j]+C_13[i][j]*dt/dx*pxVx+C_35[i][j]*dt/dx*pxVz);
                Tzz_z[i][j]=1.0/(1+0.5*dt*ddz[i][j])*((1-0.5*ddz[i][j]*dt)*Tzz_z[i][j]+C_33[i][j]*dt/dz*pzVz+C_35[i][j]*dt/dz*pzVx);
                Txx[i][j]=Txx_x[i][j]+Txx_z[i][j];
                Tzz[i][j]=Tzz_x[i][j]+Tzz_z[i][j];
                if(i==SX&&j==SZ)
                {
                    float tmp=pow(PI*f0*(k*dt-t0),2);
                    //Ricker
                    Txx_x[i][j]+=1*exp(-1.0*tmp)*(1-2.0*tmp);
                    Txx_z[i][j]+=1*exp(-1.0*tmp)*(1-2.0*tmp);
                    Tzz_x[i][j]+=1*exp(-1.0*tmp)*(1-2.0*tmp);
                    Tzz_z[i][j]+=1*exp(-1.0*tmp)*(1-2.0*tmp);
                    //Gauss
                    //Tzz[i][j]+=-1.0*PI*PI*f0*f0*(k*dt-t0)*exp(-1.0*tmp)*(3.0-2.0*tmp);
                    //Txx[i][j]+=-1.0*PI*PI*f0*f0*(k*dt-t0)*exp(-1.0*tmp)*(3.0-2.0*tmp);
                }
                pxVx=0;pzVx=0;pxVz=0;pzVz=0;
                for(int in=0;in<N;in++)
                {
                    pxVx+=cof[in]*(Vx[i+in+1][j]-Vx[i-in][j]);
                    pxVz+=cof[in]*(Vz[i+in][j]-Vz[i-in-1][j]);
                    pzVx+=cof[in]*(Vx[i][j+in+1]-Vx[i][j-in]);
                    pzVz+=cof[in]*(Vz[i][j+in]-Vz[i][j-in-1]);
                }
                Txz_x[i][j]=1.0/(1+0.5*dt*ddx[i][j])*((1-0.5*ddx[i][j]*dt)*Txz_x[i][j]+C_15[i][j]*dt/dx*pxVx+C_55[i][j]*dt/dx*pxVz);
                Txz_z[i][j]=1.0/(1+0.5*dt*ddz[i][j])*((1-0.5*ddz[i][j]*dt)*Txz_z[i][j]+C_35[i][j]*dt/dz*pzVz+C_55[i][j]*dt/dz*pzVx);
                Txz[i][j]=Txz_x[i][j]+Txz_z[i][j];
            }
        }

        for(i=N;i<NX-N;i++)
        {
            for(j=N;j<NZ-N;j++)
            {
                float pxTxx=0;float pzTxz=0;
                for(int in=0;in<N;in++)
                {
                    pxTxx+=cof[in]*(Txx[i+in][j]-Txx[i-in-1][j]);
                    pzTxz+=cof[in]*(Txz[i][j+in]-Txz[i][j-in-1]);
                }
                Vx_x[i][j]=1.0/(1+0.5*dt*ddx[i][j])*((1-0.5*ddx[i][j]*dt)*Vx_x[i][j]+dt/Rou[i][j]/dx*pxTxx);
                Vx_z[i][j]=1.0/(1+0.5*dt*ddz[i][j])*((1-0.5*ddz[i][j]*dt)*Vx_z[i][j]+dt/Rou[i][j]/dz*pzTxz);
                Vx[i][j]=Vx_x[i][j]+Vx_z[i][j];
                //if(i==SX&&j==SZ)
                //{
                //    float tmp=pow(PI*f0*(k*dt-t0),2);
                //    //Ricker
                //    Vz_x[i][j]+=1*exp(-1.0*tmp)*(1-2.0*tmp);
                //    Vz_z[i][j]+=1*exp(-1.0*tmp)*(1-2.0*tmp);
                //}
                float pxTxz=0;float pzTzz=0;
                for(int in=0;in<N;in++)
                {
                    pxTxz+=cof[in]*(Txz[i+in+1][j]-Txz[i-in][j]);
                    pzTzz+=cof[in]*(Tzz[i][j+in+1]-Tzz[i][j-in]);
                }
                Vz_x[i][j]=1.0/(1+0.5*dt*ddx[i][j])*((1-0.5*ddx[i][j]*dt)*Vz_x[i][j]+dt/Rou[i][j]/dx*pxTxz);
                Vz_z[i][j]=1.0/(1+0.5*dt*ddz[i][j])*((1-0.5*ddz[i][j]*dt)*Vz_z[i][j]+dt/Rou[i][j]/dz*pzTzz);
                Vz[i][j]=Vz_x[i][j]+Vz_z[i][j];
            }
        }


        //for (i = 0; i < NX; i++)
        for (i = pml; i < NX-pml; i++)
        {
            int rz=SZ;
            data_vx[i][k]=Txx[i][rz]+Tzz[i][rz];
            //data_vx[i][k]=Vx[i][rz];
            data_vz[i][k]=Vz[i][rz];
        }


        if(k%200==0)
        {
            //FILE *fp121;
            //char name[256];
            //sprintf(name,"T_TTI2.dat%d",k);
            //fp121=fopen(name, "wb+");
            ////for (i = 0; i < NX; i++)
            ////    for (j = 0; j < NZ; j++)
            //for (i = pml; i < NX-pml; i++)
            //    for (j = pml; j < NZ-pml; j++)
            //    {
            //        float T=1.0/2*(Txx[i][j]+Tzz[i][j]);
            //        fwrite(&T, sizeof(float), 1, fp121);
            //    }
            //fclose(fp121);
            //FILE *fp122;
            char name[256];
            //sprintf(name,"Vx.dat%d",k);
            //fp122=fopen(name, "wb+");
            ////for (i = 0; i < NX; i++)
            ////    for (j = 0; j < NZ; j++)
            //for (i = pml; i < NX-pml; i++)
            //    for (j = pml; j < NZ-pml; j++)
            //    {
            //        fwrite(&Vx[i][j], sizeof(float), 1, fp122);
            //    }
            //fclose(fp122);
            FILE *fp123;
            sprintf(name,"u.dat%d",k);
            fp123=fopen(name, "wb+");
            //for (i = 0; i < NX; i++)
            //    for (j = 0; j < NZ; j++)
            for (i = pml; i < NX-pml; i++)
                for (j = pml; j < NZ-pml; j++)
                {
                    float u=Txx[i][j]+Tzz[i][j];
                    //fwrite(&Vz[i][j], sizeof(float), 1, fp123);
                    fwrite(&u, sizeof(float), 1, fp123);
                }
            fclose(fp123);
        }

    }


    fp1=fopen("Txx.dat", "wb+");
    for (i = pml; i < NX-pml; i++)
        for (j = pml; j < NZ-pml; j++)
        {
            fwrite(&Txx[i][j], sizeof(float), 1, fp1);
        }
    fp2=fopen("Tzz.dat", "wb+");
    for (i = pml; i < NX-pml; i++)
        for (j = pml; j < NZ-pml; j++)
        {
            fwrite(&Tzz[i][j], sizeof(float), 1, fp2);
        }
    fp12=fopen("T.dat", "wb+");
    for (i = pml; i < NX-pml; i++)
        for (j = pml; j < NZ-pml; j++)
        {
            float T=1.0/2*(Txx[i][j]+Tzz[i][j]);
            fwrite(&T, sizeof(float), 1, fp12);
        }
    fp3=fopen("Txz.dat", "wb+");
    for (i = pml; i < NX-pml; i++)
        for (j = pml; j < NZ-pml; j++)
        {
            fwrite(&Txz[i][j], sizeof(float), 1, fp3);
        }
    fp4=fopen("Vx.dat", "wb+");
    for (i = pml; i < NX-pml; i++)
        for (j = pml; j < NZ-pml; j++)
        {
            fwrite(&Vx[i][j], sizeof(float), 1, fp4);
        }
    fp5=fopen("Vz.dat", "wb+");
    for (i = pml; i < NX-pml; i++)
        for (j = pml; j < NZ-pml; j++)
        {
            fwrite(&Vz[i][j], sizeof(float), 1, fp5);
        }
    //fp6=fopen("data_vx.dat", "wb+");
    fp6=fopen("uobs.dat", "wb+");
    //for (i = 0; i < NX; i++)
    for (i = pml; i < NX-pml; i++)
        for (j = 0; j < NT; j++)
        {
            fwrite(&data_vx[i][j], sizeof(float), 1, fp6);
        }
    //    fp7=fopen("data_vz.dat", "wb+");
    //    for (i = 0; i < NX; i++)
    //        for (j = 0; j < NT; j++)
    //        {
    //            fwrite(&data_vz[i][j], sizeof(float), 1, fp7);
    //        }
    fp1=fopen("Vp.dat", "wb+");
    for (i = pml; i < NX-pml; i++)
        for (j = pml; j < NZ-pml; j++)
        {
            fwrite(&Vp[i][j], sizeof(float), 1, fp1);
        }
    //We should delete memory here!!!
}
