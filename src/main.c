# include <stdio.h>
# include <math.h>

#define MAX_GPSSTA 1000
#define R 6371.393
#define PI 3.14159265358979323846
#define RAD2DEG (180.0/PI)
#define DEG2RAD	(PI/180.0)

typedef struct
{
	double b0;
	double l0;
	double d_b0;
	double d_l0;
	double step_b0;
	double step_l0;
	double depth_st;
	double depth_ed;
	double depth_step;
	double velo_st;
	double velo_ed;
	double velo_step;  
	double b0_ref;
	double l0_ref;
	double depth_ref;                        
}CFG;

typedef struct
{
	char staName[256];
	double staB;
	double staL;
	double staH;
	double t; // seconds after ot
	double t_residual;
}GPSSTA;

typedef struct
{
	double B;
	double L;
	double depth;
	double velo;
	double ot_mean;
	double ot_max;
	double ot_min;
}HYPCENTER;


double calDistance(double b1,double l1,double b2,double l2)
{
	b1=b1*DEG2RAD;
	b2=b2*DEG2RAD;
	l1=l1*DEG2RAD;
	l2=l2*DEG2RAD;
	return acos(sin(b1)*sin(b2)+cos(b1)*cos(b2)*cos(l1-l2))*R;
}

void makeList(double* list,int* num,double st,double ed,double step)
{
	int MaxNum=100000;
	double element=st;
	for(int i=0;i<MaxNum;i++)
	{
		if(element>ed) 
		{
			break;
		}
		list[*num]=element;
		*num=*num+1;
		element=element+step;
	}
}

void cal_hypCenter_GRID(HYPCENTER* hypcenter,char* logFileName,GPSSTA* gpsSta,int gpsStaNum,CFG* cfg)
{
	double misfit_min=999999;
	int count=0;

	int coupleNum=0;
	for(int i=0;i<gpsStaNum;i++)
	{
		coupleNum=coupleNum+gpsStaNum-1-i;
	}

	FILE* fp_log=fopen(logFileName,"w");

	double b_list[10000],l_list[10000],depth_list[10000],velo_list[10000];
	int b_num=0,l_num=0,depth_num=0,velo_num=0;

	if(cfg->b0==0)
	{
		double b0_mean=0;
		double l0_mean=0;
		for(int i=0;i<gpsStaNum;i++)
		{
			b0_mean=b0_mean+gpsSta[i].staB;
			l0_mean=l0_mean+gpsSta[i].staL;
		}
		b0_mean=b0_mean/gpsStaNum;
		l0_mean=l0_mean/gpsStaNum;
		cfg->b0=b0_mean;
		cfg->l0=l0_mean;
	}

	makeList(b_list,&b_num,cfg->b0-cfg->d_b0,cfg->b0+cfg->d_b0,cfg->step_b0);
	makeList(l_list,&l_num,cfg->l0-cfg->d_l0,cfg->l0+cfg->d_l0,cfg->step_l0);
	makeList(depth_list,&depth_num,cfg->depth_st,cfg->depth_ed,cfg->depth_step);
	makeList(velo_list,&velo_num,cfg->velo_st,cfg->velo_ed,cfg->velo_step);

	double misfit=0;
	double epidist_i=0,hypdist_i=0;
	double epidist_j=0,hypdist_j=0;

	for(int i=0;i<b_num;i++)
	{
		printf("\n%6d  %6.3f\n",i,b_list[i]);
    	fprintf(fp_log,"\n%6d  %6.3f\n",i,b_list[i]);
		for(int j=0;j<l_num;j++)
		{
			for(int k=0;k<depth_num;k++)
			{
				for(int v=0;v<velo_num;v++)
				{
					misfit=0;
					for(int ista=0;ista<gpsStaNum-1;ista++)
					{
						epidist_i=calDistance(b_list[i],l_list[j],gpsSta[ista].staB,gpsSta[ista].staL);
						hypdist_i=sqrt(depth_list[k]*depth_list[k]+epidist_i*epidist_i);
						for(int jsta=ista+1;jsta<gpsStaNum;jsta++)
						{
							epidist_j=calDistance(b_list[i],l_list[j],gpsSta[jsta].staB,gpsSta[jsta].staL);
						    hypdist_j=sqrt(depth_list[k]*depth_list[k]+epidist_j*epidist_j);
							double tmp=(hypdist_i-hypdist_j)/velo_list[v]-(gpsSta[ista].t-gpsSta[jsta].t);
							if(tmp<0)
							{
								tmp=tmp*(-1.0);
							}
							misfit=misfit+tmp;
						}
					}
					if(misfit<misfit_min)
					{
						misfit_min=misfit;
						hypcenter->B=b_list[i];
						hypcenter->L=l_list[j];
						hypcenter->depth=depth_list[k];
						hypcenter->velo=velo_list[v];
						count=count+1;

						printf("%6d epicenter=[%6.3f,%6.3f],depth=%6.3fkm,velo=%6.3fkm/s,misfit(sum)=%6.3fs,misfit(average)=%6.3fs\n",count,hypcenter->B,hypcenter->L,hypcenter->depth,hypcenter->velo,misfit_min,misfit_min/coupleNum);
                        fprintf(fp_log,"%6d epicenter=[%6.3f,%6.3f],depth=%6.3fkm,velo=%6.3fkm/s,misfit(sum)=%6.3fs,misfit(average)=%6.3fs\n",count,hypcenter->B,hypcenter->L,hypcenter->depth,hypcenter->velo,misfit_min,misfit_min/coupleNum);
					}
				}
			}
		}
	}
	fclose(fp_log);
}

void cal_OT(HYPCENTER* hypcenter,GPSSTA* gpsSta,int gpsStaNum)
{
	double average_ot=0,max_ot=-999999,min_ot=999999;

	for(int i=0;i<gpsStaNum;i++)
	{
		double epiDist=calDistance(hypcenter->B,hypcenter->L,gpsSta[i].staB,gpsSta[i].staL);
		double hypDist=sqrt(epiDist*epiDist+hypcenter->depth*hypcenter->depth);
		double travel_t=hypDist/hypcenter->velo;
		gpsSta[i].t_residual=gpsSta[i].t-travel_t;
		average_ot=average_ot+gpsSta[i].t_residual;
		if(gpsSta[i].t_residual>max_ot){
			max_ot=gpsSta[i].t_residual;
		}
		if(gpsSta[i].t_residual<min_ot){
			min_ot=gpsSta[i].t_residual;
		}
	}
	average_ot=average_ot/gpsStaNum;

	hypcenter->ot_mean=average_ot;
	hypcenter->ot_max=max_ot;
	hypcenter->ot_min=min_ot;
}

void readGPSstaFile(char* gpsSta_filename,GPSSTA* gpsSta,int* gpsStaNum)
{
    FILE* fp_gps=fopen(gpsSta_filename,"r");
    char line[256]="\0";

	while(!feof(fp_gps))
	{
		fgets(line,256,fp_gps);
		sscanf(line,"%s%lf%lf%lf%lf",gpsSta[*gpsStaNum].staName,&(gpsSta[*gpsStaNum].staB),&(gpsSta[*gpsStaNum].staL),&(gpsSta[*gpsStaNum].staH),&(gpsSta[*gpsStaNum].t));				
	    *gpsStaNum=*gpsStaNum+1;
	}
	fclose(fp_gps);
}

void readCFGfile(char* cfg_filename,CFG* cfg)
{
	FILE* fp_cfg=fopen(cfg_filename,"r");
    char line[256]="\0";

	fgets(line,256,fp_cfg); 
	fgets(line,256,fp_cfg); 

	fgets(line,256,fp_cfg); 
	fgets(line,256,fp_cfg);
	sscanf(line,"%lf%lf%lf",&(cfg->b0_ref),&(cfg->l0_ref),&(cfg->depth_ref));

	fgets(line,256,fp_cfg); 
	fgets(line,256,fp_cfg);
	sscanf(line,"%lf%lf",&(cfg->b0),&(cfg->l0));

	fgets(line,256,fp_cfg); 
	fgets(line,256,fp_cfg); 
	sscanf(line,"%lf%lf",&(cfg->d_b0),&(cfg->step_b0));

	fgets(line,256,fp_cfg); 
    fgets(line,256,fp_cfg); 
	sscanf(line,"%lf%lf",&(cfg->d_l0),&(cfg->step_l0));

	fgets(line,256,fp_cfg); 
    fgets(line,256,fp_cfg); 
	sscanf(line,"%lf%lf%lf",&(cfg->depth_st),&(cfg->depth_ed),&(cfg->depth_step));

	fgets(line,256,fp_cfg); 
    fgets(line,256,fp_cfg); 
	sscanf(line,"%lf%lf%lf",&(cfg->velo_st),&(cfg->velo_ed),&(cfg->velo_step));

	fclose(fp_cfg);
}

void writeHypRes(HYPCENTER* hypcenter,CFG* cfg,char* hypResFileName,GPSSTA* gpsSta,int gpsStaNum)
{
	FILE* fp_res=fopen(hypResFileName,"w");
	fprintf(fp_res,"\n\n******************************Results********************************\n");
	fprintf(fp_res,">>epicenter_b=%12.4f[%12.4f,%12.4f]\n",hypcenter->B,cfg->b0-cfg->d_b0,cfg->b0+cfg->d_b0);
	fprintf(fp_res,">>epicenter_l=%12.4f[%12.4f,%12.4f]\n",hypcenter->L,cfg->l0-cfg->d_l0,cfg->l0+cfg->d_l0);
	fprintf(fp_res,">>      depth=%12.4f[%12.4f,%12.4f]\n",hypcenter->depth,cfg->depth_st,cfg->depth_ed);
	fprintf(fp_res,">>       velo=%12.4f[%12.4f,%12.4f]\n",hypcenter->velo,cfg->velo_st,cfg->velo_ed);
	fprintf(fp_res,">>         ot=%12.4f[%12.4f,%12.4f]\n",hypcenter->ot_mean,hypcenter->ot_min,hypcenter->ot_max);
	fprintf(fp_res,">>>hypCenterRefe=(%12.4f,%12.4f,%12.4f)\n",cfg->b0_ref,cfg->l0_ref,cfg->depth_ref);
	fprintf(fp_res,">>>epiCenterBias=%12.4fkm\n",calDistance(cfg->b0_ref,cfg->l0_ref,hypcenter->B,hypcenter->L));
	fprintf(fp_res,"*****Number,staCode,Latitude,Longitude,Height,t_residual(or OT)*****\n");
	for(int i=0;i<gpsStaNum;i++)
	{
		fprintf(fp_res,"No.%5d %6s %12.4f %12.4f %12.4f %12.4f |\n",i+1,gpsSta[i].staName,gpsSta[i].staB,gpsSta[i].staL,gpsSta[i].staH,gpsSta[i].t_residual);
	}
	fclose(fp_res);

	printf("\n\n******************************Results********************************\n");
	printf(">>epicenter_b=%12.4f[%12.4f,%12.4f]\n",hypcenter->B,cfg->b0-cfg->d_b0,cfg->b0+cfg->d_b0);
	printf(">>epicenter_l=%12.4f[%12.4f,%12.4f]\n",hypcenter->L,cfg->l0-cfg->d_l0,cfg->l0+cfg->d_l0);
	printf(">>      depth=%12.4f[%12.4f,%12.4f]\n",hypcenter->depth,cfg->depth_st,cfg->depth_ed);
	printf(">>       velo=%12.4f[%12.4f,%12.4f]\n",hypcenter->velo,cfg->velo_st,cfg->velo_ed);
	printf(">>         ot=%12.4f[%12.4f,%12.4f]\n",hypcenter->ot_mean,hypcenter->ot_min,hypcenter->ot_max);
	printf(">>>hypCenterRefe=(%12.4f,%12.4f,%12.4f)\n",cfg->b0_ref,cfg->l0_ref,cfg->depth_ref);
	printf(">>>epiCenterBias=%12.4fkm\n",calDistance(cfg->b0_ref,cfg->l0_ref,hypcenter->B,hypcenter->L));
	printf("*****Number,staCode,Latitude,Longitude,Height,t_residual(or OT)*****\n");
	for(int i=0;i<gpsStaNum;i++)
	{
		printf("No.%5d %6s %12.4f %12.4f %12.4f %12.4f |\n",i+1,gpsSta[i].staName,gpsSta[i].staB,gpsSta[i].staL,gpsSta[i].staH,gpsSta[i].t_residual);
	}

}

int main(int argc, char* argv[])
{
	char cfg_filename[256]="../workspace/inputs/cfg_hyp.txt";
	char gpsSta_filename[256]="../workspace/inputs/gpsSta.dat";
	char logFileName[256]="../workspace/results/logHyp.txt";
	char hypResFileName[256]="../workspace/results/hypRes.txt";

	CFG cfg;
	GPSSTA gpsSta[MAX_GPSSTA];
	HYPCENTER hypcenter;
	int gpsStaNum=0;

    readCFGfile(cfg_filename,&cfg);
	readGPSstaFile(gpsSta_filename,gpsSta,&gpsStaNum);
	cal_hypCenter_GRID(&hypcenter,logFileName,gpsSta,gpsStaNum,&cfg);
	cal_OT(&hypcenter,gpsSta,gpsStaNum);
    writeHypRes(&hypcenter,&cfg,hypResFileName,gpsSta,gpsStaNum);

	return 0;
}