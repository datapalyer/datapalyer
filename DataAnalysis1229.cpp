/*****************************
2PI数据处理程程序2.3版
(1)DataAnalysisx.x.cpp
版本1.2
1.定义了参数结构体，事例结构体，粒子结构体；
2.读入处理程序需要的一切参数，处理数据。

版本1.3更新
1.增加了找角度圆心的代码(analysis_mode=4)。

版本2.0更新
1.添加注释；
2.增加了从命令行读入参数文件名，从参数文件中读入参数的代码，这样可以在项目属性->调试中加入命令参数调试程序。

版本3.0更新
1.符合增加了单路角度范围参数Phi1_min,Phi1_max,Phi2_min,Phi2_max，输入文件需要修改。

版本3.1更新：
1.符合增加了能谱角度二维分布图。

版本4.0更新：
1.由于uv，uw，vw图像不能完全重合，三种组合分开刻度求能量，参数结构体parameter和输入文件发生变化
2.弹性分角度刻度
3.能量的计算分uv，uw,vw来写，getenergy()函数重写

版本4.1更新
1.稍修改符合判断函数的限定范围
2.增加了模式5，输出能量和角度分布

版本4.3更新
1.增加模式6，直接(x,y)输出
2.整合所有函数在common文件中

版本5.0更新：
1.增加了每个文件入射能量可以修正的内容，E0_modify[num_files];
2.输入文件有变动，在name后面要加上修改的值。

版本5.1更新：
1.角度刻度,角度刻度采用20160524的结果
2.修改bug，ems模式中ow忘记写入（+param.ow）

版本5.2更新：
1.偶然符合时间窗对称取min<fabs(tdiff)<max,即对accidentalcoincidence()函数修改了下
版本5.2s更新：
个别地方（+param.ow）忘加

版本20161031更新：
角度刻度为20161031结果

版本20161229_ver2
角度刻度为20161229结果
角度圆心可以改变

版本2PIDataAnalysis20190108_ver2_2
20190106刻度 角度16~160 192~340（原194~340）
******************************/
#include <stdio.h>
#include<stdlib.h>
#include<string.h>
#include<iostream>
#include<fstream>
#include<iomanip>

#define MAX_LINE 200
typedef unsigned int uint32;  //32位无符号整型
typedef unsigned char uint8;  //8位无符号整型
#define max_hit 3  //最大多击数目
#define max_ch 16  //最大通道数
#define u1ch 0  //u1,u2,v1,v2,w1,w2,MCP对应的通道
#define u2ch 2
#define v1ch 4
#define v2ch 6
#define w1ch 8
#define w2ch 10
#define MCPch 12       //位置信号和时间信号对应的TDC编号
#define MAX_origin_num 5000   //输出原始数据最大个数，用于查看原始数据是否正确
using namespace std;

//参数结构体，所有参数都记录在这个结构体中
typedef struct
{
	char workspace[MAX_LINE];   //处理结果文件存放的位置，注意要以“/”结尾
	int uTsum_Lower;            //u路Tsum的最小值
	int uTsum_Up;               //u路Tsum的最大值
	int vTsum_Lower;	//v路
	int vTsum_Up;		
	int wTsum_Lower;	//w路
	int wTsum_Up;
	double du;    //Half of the signal transmission speed on u-layer
	double dv;    //Half of the signal transmission speed on v-layer
	double dw;    //Half of the signal transmission speed on w-layer
	double ow;    //Offset of w-layer
	double rotate_angle; //rotate angle
	double x1_energy_offset; //能量圆心(x<0)
	double y1_energy_offset; 
	double x2_energy_offset; //能量圆心(x>0)
	double y2_energy_offset;
	double x1_angle_offset;  //求角分布时圆心的偏移，注意这个是在能量圆心的基础上继续偏移
	double y1_angle_offset;
	double x2_angle_offset;
	double y2_angle_offset;
	int analysis_mode;      //数据处理模式，1->Tsum谱，2->弹性刻度,3->符合，4->角分布找圆心
	int num_files;          //处理文件的数目
	char **file;            //处理文件的路径及文件名
	char **name;            //结果文件的前缀
	double *E0_modify;
	int time_min,time_max,time_step;  //Tsum时间谱的最小值，最大值和步长
	double position_min,position_max,position_step;  //位置分布的最小值，最大值和步长
	double radius_min,radius_max,radius_step;        //半径分布的最小值，最大值和步长
	double phi_min,phi_max,phi_step;  //角度分布的最小值，最大值和步长
	double tdiff_min,tdiff_max,tdiff_step;  ////符合时间谱的最小值，最大值和步长
	double *cuv0,*cuv1,*cuv2;  //E=c0+c1*r+c2*r^2
	double *cuw0,*cuw1,*cuw2;
	double *cvw0,*cvw1,*cvw2;
	int cali_angle_num;  //the number of angles for esnergy calibration
	double E0,E0_shift,energy_min,energy_max,energy_step; //入射能量，能量修正，能谱的最小值，最大值和步长
	double angle_min,angle_max,angle_step; //方位角分布的最小值，最大值和步长
	double E12_min,E12max; //单路E1或者E2的的最小值，最大值和步长
	double Phi1_min,Phi1_max,Phi2_min,Phi2_max;
	double truecoincidence1,truecoincidence2; //真符合窗的范围
	double accidentalcoincidence1,accidentalcoincidence2;  //偶然符合窗的范围
	int state_num;   //作角分布时能量区间的数目
	double *state_energy_min,*state_energy_max,*state_energy_mid;//作角分布时能量区间的的最小值，最大值和时间能峰的位置
	double *elastic_energy; //单路刻度每个数据文件对应的能量值
//	double auv1[4][4],auw1[4][4],avw1[4][4];//角度刻度系数
//	double auv2[4][4],auw2[4][4],avw2[4][4];//角度刻度系数
}parameter;  

//事例结构体，TDC记录事件的信息
typedef struct
{
	int count;  //事例编号
	int u1[max_hit],u2[max_hit],v1[max_hit],v2[max_hit],w1[max_hit],w2[max_hit],MCP[max_hit];
	int nhit[max_ch]; //每个通道多击的数目
	int uTsum,vTsum,wTsum;
	//pattern:第n击通道是否全的判断，u1->000001(1),u2->000010(2),u3->000100(4),u4->001000(8),w1->010000(16),w2->100000(32)
	int pattern[max_hit];  
	int status;  //status shown in global trailer
}event;

//粒子结构体，单个粒子的信息
typedef struct
{
	double u,v,w;  //u1-u2,v1-v2...
	double uTsum,vTsum,wTsum; //u1+u2,...
	int uflag,vflag,wflag;   //Tsum是否在真实Tsum范围内
	int whichtwolayer;//哪两层计算的(x,y)和(r,phi),1->uv,2->uw,3->vw
	double xuv,yuv,xuw,yuw,xvw,yvw,x,y; 
	double ruv,ruw,rvw,phiuv,phiuw,phivw,r,phi;
	double euv,euw,evw,auv,auw,avw; //能量和角度刻度的结果
	double energy,angle; 
}particle;

int read_param(parameter &param);
int read_param_fromfile(parameter &param,char *name);
int print_param(parameter param);

#include"common.h"
#include"getTsum.h"
#include"elastic.h"
#include"ems.h"

int main(int argc,char* argv[])
{
	parameter param;
	if(argc>1)
		read_param_fromfile(param,argv[1]);
	else
		read_param(param);
//	print_param(param);
	if(param.analysis_mode==1)
		getTsum(param);
	else if(param.analysis_mode==2)
		elastic(param);
	else if(param.analysis_mode==3)
		ems(param);
	else if(param.analysis_mode==4)
	{
		double scale,step;
		cout<<"Input the varation scale and step:"<<endl;
		cin>>scale>>step;
		int i,j,n;
		n=(int)(scale/step+1.5);
		cout<<"Scale: "<<scale<<"\t"<<"Step: "<<step<<"\t"<<n<<endl;
		
		double x_start=param.x2_angle_offset-scale/2;
		param.x1_angle_offset=param.x2_angle_offset=x_start;
		double y_start=param.y2_angle_offset-scale/2;
		param.y1_angle_offset=param.y2_angle_offset=y_start;
		for(i=0;i<n;i++)
		{
			param.y1_angle_offset=param.y2_angle_offset=y_start;	
			for(j=0;j<n;j++)
			{
				cout<<"Center: "<<param.x1_angle_offset<<"\t"<<param.y1_angle_offset<<endl;
				sprintf(param.name[0],"_%.1f_%.1f",param.x1_angle_offset,param.y1_angle_offset);
				ems(param);
				param.y2_angle_offset+=step;
				param.y1_angle_offset=param.y2_angle_offset;	
			}
			param.x2_angle_offset+=step;
			param.x1_angle_offset=param.x2_angle_offset;
		}
	}
	else if(param.analysis_mode==5)
		noncoin(param);
	else if(param.analysis_mode==6)
		directxy(param);
	else
		return 0;
	return 1;

}

//从输入流中读参数信息
int read_param(parameter &param)
{
	char buff[MAX_LINE];
	scanf("%s%s",buff,param.workspace);   
	printf("workspace: %s\n",param.workspace); 
	scanf("%s%d%d",buff,&param.uTsum_Lower,&param.uTsum_Up); 
	printf("uTsum %d %d\n",param.uTsum_Lower,param.uTsum_Up); 
	scanf("%s%d%d",buff,&param.vTsum_Lower,&param.vTsum_Up); 
	printf("vTsum %d %d\n",param.vTsum_Lower,param.vTsum_Up);  
	scanf("%s%d%d",buff,&param.wTsum_Lower,&param.wTsum_Up);
	printf("wTsum %d %d\n",param.wTsum_Lower,param.wTsum_Up); 
	scanf("%s%lf%lf%lf",buff,&param.du,&param.dv,&param.dw); 
	printf("velocity: %f %f %f \n",param.du,param.dv,param.dw); 
	scanf("%s%lf",buff,&param.ow); 
	printf("w offset %f \n",param.ow);
	scanf("%s%lf",buff,&param.rotate_angle); 
	printf("rotate angle: %f \n",param.rotate_angle);
	scanf("%s%lf%lf",buff,&param.x1_energy_offset,&param.y1_energy_offset); 
	printf("xy energy offset for x<0: %f %f \n",param.x1_energy_offset,param.y1_energy_offset);
	scanf("%s%lf%lf",buff,&param.x2_energy_offset,&param.y2_energy_offset);
	printf("xy energy offset for x>0: %f %f \n",param.x2_energy_offset,param.y2_energy_offset);
	scanf("%s%lf%lf",buff,&param.x1_angle_offset,&param.y1_angle_offset);
	printf("xy angle offset for x<0: %f %f \n",param.x1_angle_offset,param.y1_angle_offset);
	scanf("%s%lf%lf",buff,&param.x2_angle_offset,&param.y2_angle_offset); 
	printf("xy angle offset for x>0: %f %f \n",param.x2_angle_offset,param.y2_angle_offset);
	scanf("%s%d%d%d",buff,&param.time_min,&param.time_max,&param.time_step); 
	printf("Time statistics: %d %d %d\n",param.time_min,param.time_max,param.time_step); 
	scanf("%s%lf%lf%lf",buff,&param.position_min,&param.position_max,&param.position_step);
	printf("Position statistics: %f %f %f \n",param.position_min,param.position_max,param.position_step); 
	scanf("%s%lf%lf%lf",buff,&param.radius_min,&param.radius_max,&param.radius_step); 
	printf("Radius statistics: %f %f %f \n",param.radius_min,param.radius_max,param.radius_step);
	scanf("%s%lf%lf%lf",buff,&param.phi_min,&param.phi_max,&param.phi_step);
	printf("Phi_angle statistics: %f %f %f \n",param.phi_min,param.phi_max,param.phi_step); 
	scanf("%s%lf%lf%lf",buff,&param.tdiff_min,&param.tdiff_max,&param.tdiff_step);  
	printf("Time difference statistics: %f %f %f \n",param.tdiff_min,param.tdiff_max,param.tdiff_step);
	scanf("%s%lf",buff,&param.E0);  
	printf("Incident energy: %f\n",param.E0);
	scanf("%s%lf",buff,&param.E0_shift);  
	printf("Energy shift: %f\n",param.E0_shift); 
	scanf("%s%lf%lf%lf",buff,&param.energy_min,&param.energy_max,&param.energy_step);  
	printf("Binding energy step: %f %f %f\n",param.energy_min,param.energy_max,param.energy_step); 
	scanf("%s%lf%lf",buff,&param.E12_min,&param.E12max);
	printf("E12 min max: %f %f\n",param.E12_min,param.E12max);
	scanf("%s%lf%lf%lf%lf",buff,&param.Phi1_min,&param.Phi1_max,&param.Phi2_min,&param.Phi2_max);
	printf("Phi12 min max: %f %f %f %f\n",param.Phi1_min,param.Phi1_max,param.Phi2_min,param.Phi2_max);
	scanf("%s%lf%lf%lf",buff,&param.angle_min,&param.angle_max,&param.angle_step); 
	printf("Azimuth angle step: %f %f %f\n",param.angle_min,param.angle_max,param.angle_step);
	scanf("%s%lf%lf",buff,&param.truecoincidence1,&param.truecoincidence2);
	printf("True coincidence gate: %f %f\n",param.truecoincidence1,param.truecoincidence2); 
	scanf("%s%lf%lf",buff,&param.accidentalcoincidence1,&param.accidentalcoincidence2);
	printf("Accidental coincidence gate: %f %f\n",param.accidentalcoincidence1,param.accidentalcoincidence2); 	
	scanf("%s%d",buff,&param.state_num);
	param.state_energy_mid=(double *)malloc(sizeof(double)*param.state_num);
	param.state_energy_min=(double *)malloc(sizeof(double)*param.state_num);
	param.state_energy_max=(double *)malloc(sizeof(double)*param.state_num);
	for(int i=0;i<param.state_num;i++) 	
		scanf("%s%lf%lf%lf",buff,&param.state_energy_mid[i],&param.state_energy_min[i],&param.state_energy_max[i]);
	printf("NO. of state: %d\n",param.state_num);   
	for(int i=0;i<param.state_num;i++) 
		printf("state %d: %f %f %f\n",i+1,param.state_energy_mid[i],param.state_energy_min[i],param.state_energy_max[i]);  
	
	scanf("%s%d",buff,&param.cali_angle_num); 
	param.cuv0=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cuv1=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cuv2=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cuw0=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cuw1=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cuw2=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cvw0=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cvw1=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cvw2=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	for(int i=0;i<param.cali_angle_num;i++)
	{
		scanf("%s%lf%lf%lf",buff,&param.cuv0[i],&param.cuv1[i],&param.cuv2[i]);
		scanf("%lf%lf%lf",&param.cuw0[i],&param.cuw1[i],&param.cuw2[i]);
		scanf("%lf%lf%lf",&param.cvw0[i],&param.cvw1[i],&param.cvw2[i]);
	}
	printf("Energy calibration for each angle for uv:\n");
	for(int i=0;i<param.cali_angle_num;i++) 
		printf("Angle %.1f~%.1f: %f %f %f\n",(360.0/param.cali_angle_num)*i,(360.0/param.cali_angle_num)*(i+1),param.cuv0[i],param.cuv1[i],param.cuv2[i]); 
	printf("Energy calibration for each angle for uw:\n");
	for(int i=0;i<param.cali_angle_num;i++) 
		printf("Angle %.1f~%.1f: %f %f %f\n",(360.0/param.cali_angle_num)*i,(360.0/param.cali_angle_num)*(i+1),param.cuw0[i],param.cuw1[i],param.cuw2[i]);
	printf("Energy calibration for each angle for vw:\n");
	for(int i=0;i<param.cali_angle_num;i++) 
		printf("Angle %.1f~%.1f: %f %f %f\n",(360.0/param.cali_angle_num)*i,(360.0/param.cali_angle_num)*(i+1),param.cvw0[i],param.cvw1[i],param.cvw2[i]); 

	scanf("%s%d",buff,&param.analysis_mode); printf("Analysis mode: %d->%s\n",param.analysis_mode,mode_name(param.analysis_mode)); 
	scanf("%s%d",buff,&param.num_files);     printf("Number of files: %d\n",param.num_files);  
	param.file=(char **)malloc(sizeof(char *)*param.num_files);  
	param.name=(char **)malloc(sizeof(char *)*param.num_files);
	param.E0_modify=(double *)malloc(sizeof(double)*param.num_files);
	for(int i=0;i<param.num_files;i++)
	{
		param.file[i]=(char *)malloc(sizeof(char *)*MAX_LINE);   
		param.name[i]=(char *)malloc(sizeof(char *)*MAX_LINE);   
		scanf("%s%s",buff,param.file[i]);
		scanf("%s%s%lf",buff,param.name[i],&param.E0_modify[i]);
	}
	for(int i=0;i<param.num_files;i++) 
	{ 
		printf("File: %s\n",param.file[i]);  
		printf("Name: %s %f\n",param.name[i],param.E0_modify[i]); 
	}
	if(param.analysis_mode==2)
	{
		param.elastic_energy=(double *)malloc(sizeof(double)*param.num_files);
		scanf("%s",buff);
		for(int i=0;i<param.num_files;i++)
		{
			scanf("%lf",&param.elastic_energy[i]);
		}
	}
	return 1;
}

//从参数文件中读入参数
int read_param_fromfile(parameter &param,char *name)
{
	FILE *fp;
	if( (fp=fopen(name,"r"))==NULL )
	{
		cout<<"Parameter file error!";
		return 0;
	}
	char buff[128];
	fscanf(fp,"%s%s",buff,param.workspace);   
	printf("workspace: %s\n",param.workspace); 
	fscanf(fp,"%s%d%d",buff,&param.uTsum_Lower,&param.uTsum_Up); 
	printf("uTsum %d %d\n",param.uTsum_Lower,param.uTsum_Up); 
	fscanf(fp,"%s%d%d",buff,&param.vTsum_Lower,&param.vTsum_Up); 
	printf("vTsum %d %d\n",param.vTsum_Lower,param.vTsum_Up);  
	fscanf(fp,"%s%d%d",buff,&param.wTsum_Lower,&param.wTsum_Up);
	printf("wTsum %d %d\n",param.wTsum_Lower,param.wTsum_Up); 
	fscanf(fp,"%s%lf%lf%lf",buff,&param.du,&param.dv,&param.dw); 
	printf("velocity: %f %f %f \n",param.du,param.dv,param.dw); 
	fscanf(fp,"%s%lf",buff,&param.ow); 
	printf("w offset %f \n",param.ow);
	fscanf(fp,"%s%lf",buff,&param.rotate_angle); 
	printf("rotate angle: %f \n",param.rotate_angle);
	fscanf(fp,"%s%lf%lf",buff,&param.x1_energy_offset,&param.y1_energy_offset); 
	printf("xy energy offset for x<0: %f %f \n",param.x1_energy_offset,param.y1_energy_offset);
	fscanf(fp,"%s%lf%lf",buff,&param.x2_energy_offset,&param.y2_energy_offset);
	printf("xy energy offset for x>0: %f %f \n",param.x2_energy_offset,param.y2_energy_offset);
	fscanf(fp,"%s%lf%lf",buff,&param.x1_angle_offset,&param.y1_angle_offset);
	printf("xy angle offset for x<0: %f %f \n",param.x1_angle_offset,param.y1_angle_offset);
	fscanf(fp,"%s%lf%lf",buff,&param.x2_angle_offset,&param.y2_angle_offset); 
	printf("xy angle offset for x>0: %f %f \n",param.x2_angle_offset,param.y2_angle_offset);
	fscanf(fp,"%s%d%d%d",buff,&param.time_min,&param.time_max,&param.time_step); 
	printf("Time statistics: %d %d %d\n",param.time_min,param.time_max,param.time_step); 
	fscanf(fp,"%s%lf%lf%lf",buff,&param.position_min,&param.position_max,&param.position_step);
	printf("Position statistics: %f %f %f \n",param.position_min,param.position_max,param.position_step); 
	fscanf(fp,"%s%lf%lf%lf",buff,&param.radius_min,&param.radius_max,&param.radius_step); 
	printf("Radius statistics: %f %f %f \n",param.radius_min,param.radius_max,param.radius_step);
	fscanf(fp,"%s%lf%lf%lf",buff,&param.phi_min,&param.phi_max,&param.phi_step);
	printf("Phi_angle statistics: %f %f %f \n",param.phi_min,param.phi_max,param.phi_step); 
	fscanf(fp,"%s%lf%lf%lf",buff,&param.tdiff_min,&param.tdiff_max,&param.tdiff_step);  
	printf("Time difference statistics: %f %f %f \n",param.tdiff_min,param.tdiff_max,param.tdiff_step);
	fscanf(fp,"%s%lf",buff,&param.E0);  
	printf("Incident energy: %f\n",param.E0);
	fscanf(fp,"%s%lf",buff,&param.E0_shift);  
	printf("Energy shift: %f\n",param.E0_shift); 
	fscanf(fp,"%s%lf%lf%lf",buff,&param.energy_min,&param.energy_max,&param.energy_step);  
	printf("Binding energy step: %f %f %f\n",param.energy_min,param.energy_max,param.energy_step); 
	fscanf(fp,"%s%lf%lf",buff,&param.E12_min,&param.E12max);
	printf("E12 min max: %f %f\n",param.E12_min,param.E12max);
	fscanf(fp,"%s%lf%lf%lf%lf",buff,&param.Phi1_min,&param.Phi1_max,&param.Phi2_min,&param.Phi2_max);
	printf("Phi12 min max: %f %f %f %f\n",param.Phi1_min,param.Phi1_max,param.Phi2_min,param.Phi2_max);
	fscanf(fp,"%s%lf%lf%lf",buff,&param.angle_min,&param.angle_max,&param.angle_step); 
	printf("Azimuth angle step: %f %f %f\n",param.angle_min,param.angle_max,param.angle_step);
	fscanf(fp,"%s%lf%lf",buff,&param.truecoincidence1,&param.truecoincidence2);
	printf("True coincidence gate: %f %f\n",param.truecoincidence1,param.truecoincidence2); 
	fscanf(fp,"%s%lf%lf",buff,&param.accidentalcoincidence1,&param.accidentalcoincidence2);
	printf("Accidental coincidence gate: %f %f\n",param.accidentalcoincidence1,param.accidentalcoincidence2); 	
	fscanf(fp,"%s%d",buff,&param.state_num);
	param.state_energy_mid=(double *)malloc(sizeof(double)*param.state_num);
	param.state_energy_min=(double *)malloc(sizeof(double)*param.state_num);
	param.state_energy_max=(double *)malloc(sizeof(double)*param.state_num);
	for(int i=0;i<param.state_num;i++) 	
		fscanf(fp,"%s%lf%lf%lf",buff,&param.state_energy_mid[i],&param.state_energy_min[i],&param.state_energy_max[i]);
	printf("NO. of state: %d\n",param.state_num);   
	for(int i=0;i<param.state_num;i++) 
		printf("state %d: %f %f %f\n",i+1,param.state_energy_mid[i],param.state_energy_min[i],param.state_energy_max[i]);  
	
	fscanf(fp,"%s%d",buff,&param.cali_angle_num); 
	param.cuv0=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cuv1=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cuv2=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cuw0=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cuw1=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cuw2=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cvw0=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cvw1=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	param.cvw2=(double *)malloc(sizeof(double)*param.cali_angle_num); 
	for(int i=0;i<param.cali_angle_num;i++)
	{
		fscanf(fp,"%s%lf%lf%lf",buff,&param.cuv0[i],&param.cuv1[i],&param.cuv2[i]);
		fscanf(fp,"%lf%lf%lf",&param.cuw0[i],&param.cuw1[i],&param.cuw2[i]);
		fscanf(fp,"%lf%lf%lf",&param.cvw0[i],&param.cvw1[i],&param.cvw2[i]);
	}
	printf("Energy calibration for each angle for uv:\n");
	for(int i=0;i<param.cali_angle_num;i++) 
		printf("Angle %.1f~%.1f: %f %f %f\n",(360.0/param.cali_angle_num)*i,(360.0/param.cali_angle_num)*(i+1),param.cuv0[i],param.cuv1[i],param.cuv2[i]); 
	printf("Energy calibration for each angle for uw:\n");
	for(int i=0;i<param.cali_angle_num;i++) 
		printf("Angle %.1f~%.1f: %f %f %f\n",(360.0/param.cali_angle_num)*i,(360.0/param.cali_angle_num)*(i+1),param.cuw0[i],param.cuw1[i],param.cuw2[i]);
	printf("Energy calibration for each angle for vw:\n");
	for(int i=0;i<param.cali_angle_num;i++) 
		printf("Angle %.1f~%.1f: %f %f %f\n",(360.0/param.cali_angle_num)*i,(360.0/param.cali_angle_num)*(i+1),param.cvw0[i],param.cvw1[i],param.cvw2[i]); 

	fscanf(fp,"%s%d",buff,&param.analysis_mode); printf("Analysis mode: %d->%s\n",param.analysis_mode,mode_name(param.analysis_mode)); 
	fscanf(fp,"%s%d",buff,&param.num_files);     printf("Number of files: %d\n",param.num_files);  
	param.file=(char **)malloc(sizeof(char *)*param.num_files);  
	param.name=(char **)malloc(sizeof(char *)*param.num_files);
	param.E0_modify=(double *)malloc(sizeof(double)*param.num_files);
	for(int i=0;i<param.num_files;i++)
	{
		param.file[i]=(char *)malloc(sizeof(char *)*MAX_LINE);   
		param.name[i]=(char *)malloc(sizeof(char *)*MAX_LINE);   
		fscanf(fp,"%s%s",buff,param.file[i]);
		fscanf(fp,"%s%s%lf",buff,param.name[i],&param.E0_modify[i]);
	}
	for(int i=0;i<param.num_files;i++) 
	{ 
		printf("File: %s\n",param.file[i]);  
		printf("Name: %s %f\n",param.name[i],param.E0_modify[i]); 
	}
	if(param.analysis_mode==2)
	{
		param.elastic_energy=(double *)malloc(sizeof(double)*param.num_files);
		fscanf(fp,"%s",buff);
		for(int i=0;i<param.num_files;i++)
		{
			fscanf(fp,"%lf",&param.elastic_energy[i]);
		}
	}
	fclose(fp);  
	return 1;
}
int print_param(parameter param)
{
	FILE *fp;
	char name[MAX_LINE];
	sprintf(name,"%s%s_sum_input.txt",param.workspace,param.name[0]);
	if( (fp=fopen(name,"w"))==NULL )
	{
		cout<<"Parameter file error!";
		return 0;
	}
	fprintf(fp,"workspace: %s\n",param.workspace); 
	fprintf(fp,"uTsum: %d %d\n",param.uTsum_Lower,param.uTsum_Up); 
	fprintf(fp,"vTsum: %d %d\n",param.vTsum_Lower,param.vTsum_Up);  
	fprintf(fp,"wTsum: %d %d\n",param.wTsum_Lower,param.wTsum_Up); 
	fprintf(fp,"HalfVelocity(mm/ns): %f %f %f \n",param.du,param.dv,param.dw); 
	fprintf(fp,"w_offset: %.2f \n",param.ow);
	fprintf(fp,"RotateAngle: %.1f \n",param.rotate_angle);
	fprintf(fp,"xy_energy_offset1: %.2f %.2f \n",param.x1_energy_offset,param.y1_energy_offset);
	fprintf(fp,"xy_energy_offset2: %.2f %.2f \n",param.x2_energy_offset,param.y2_energy_offset);
	fprintf(fp,"xy_angle_offset1: %.2f %.2f \n",param.x1_angle_offset,param.y1_angle_offset);
	fprintf(fp,"xy_angle_offset2: %.2f %.2f \n",param.x2_angle_offset,param.y2_angle_offset);
	fprintf(fp,"Tsum_step(ch): %d %d %d\n",param.time_min,param.time_max,param.time_step); 
	fprintf(fp,"Postion_step(mm): %.1f %.1f %.1f \n",param.position_min,param.position_max,param.position_step); 
	fprintf(fp,"Radius_step(mm): %.1f %.1f %.1f \n",param.radius_min,param.radius_max,param.radius_step);
	fprintf(fp,"Phi_step(degree): %.1f %.1f %.1f \n",param.phi_min,param.phi_max,param.phi_step); 
	fprintf(fp,"Tdiff_step(ns): %.1f %.1f %.1f \n",param.tdiff_min,param.tdiff_max,param.tdiff_step);
	fprintf(fp,"E0(eV): %.3f\n",param.E0);
	fprintf(fp,"Energy_shift(eV): %.3f\n",param.E0_shift); 
	fprintf(fp,"E_step(eV): %.3f %.3f %.3f\n",param.energy_min,param.energy_max,param.energy_step); 
	fprintf(fp,"E12_min_max: %.1f %.1f\n",param.E12_min,param.E12max);
	fprintf(fp,"Phi12_min_max: %.1f %.1f %.1f %.1f\n",param.Phi1_min,param.Phi1_max,param.Phi2_min,param.Phi2_max);
	fprintf(fp,"Azimuth_angle_step: %.1f %.1f %.1f\n",param.angle_min,param.angle_max,param.angle_step);
	fprintf(fp,"True_coincidence_gate(ns): %.1f %.1f\n",param.truecoincidence1,param.truecoincidence2); 
	fprintf(fp,"Accidental_coincidence_gate(ns): %.1f %.1f\n",param.accidentalcoincidence1,param.accidentalcoincidence2); 	
	fprintf(fp,"NO_of_States: %d\n",param.state_num);   
	for(int i=0;i<param.state_num;i++) 
		fprintf(fp,"state%d: %.2f %.2f %.2f\n",i+1,param.state_energy_mid[i],param.state_energy_min[i],param.state_energy_max[i]);  
	fprintf(fp,"\n");
	fprintf(fp,"Energy_calibration_for_uv_uw_vw %d\n",param.cali_angle_num);
	for(int i=0;i<param.cali_angle_num;i++)
	{
		fprintf(fp,"Angle_%.1f~%.1f: %f %f %f ",(360.0/param.cali_angle_num)*i,(360.0/param.cali_angle_num)*(i+1),param.cuv0[i],param.cuv1[i],param.cuv2[i]);  
		fprintf(fp,"%f %f %f ",param.cuw0[i],param.cuw1[i],param.cuw2[i]);
		fprintf(fp,"%f %f %f\n",param.cvw0[i],param.cvw1[i],param.cvw2[i]);
	}
	fprintf(fp,"\n");
	fprintf(fp,"Analysis_Mode: %d\n",param.analysis_mode); 
	fprintf(fp,"NO.ofFiles: %d\n",param.num_files);  
	for(int i=0;i<param.num_files;i++) 
	{ 
		fprintf(fp,"File: %s\n",param.file[i]);  
		fprintf(fp,"Name: %s  %.3f\n",param.name[i],param.E0_modify[i]); 
	}
	if(param.analysis_mode==2)
	{
		fprintf(fp,"Elastic_energy:");
		for(int i=0;i<param.num_files;i++)
		{
			fprintf(fp," %.2f",param.elastic_energy[i]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);  
	return 1;
}

/*****************************
 输入文件内容实例：
 ----------------------------------------------------------
Workspace D:\2PIData\ems\
uTsum	7000	7200
vTsum	6650	7100
wTsum	6700	7050
HalfVelocity(mm/ns)	0.326	0.326	0.333
w_offset(mm)	0.06
RotateAngle	0
xy_energy_offset1(mm)	-1.76	-0.06
xy_energy_offset2(mm)	-1.76	-0.06
xy_angle_offset1(mm)	1.20	-0.50
xy_angle_offset2(mm)	1.20	-0.50
Tsum_step(ch)	0	9000	10
Postion_step(mm)	-60.0	60.0	0.2
Radius_step(mm)	0	60.0	0.2
Phi_step(degree)	0.0	360.0	8.0
Tdiff_step(ns)	-300	300	1.0
E0(eV)	1230.0
Energy_shift(eV)	0.00
E_step(eV)	-50.125	100.125	0.25
E12_min_max	580	625
Phi12min_max 0 180 180 360
Azimuth_angle_step	-179.5	179.5	1.0
True_coincidence(ns)	-4.0	4.0
Accidental_coincidence(ns)	10.0	50.0
NO_of_States	2
3p	15.76	15.00	22.00
3s	29.24	28.00	35.00

Energy_calibration	1
AllAngles	655.060498	-2.081001	0.007549  655.060498	-2.081001	0.007549  655.060498	-2.081001	0.007549

Analysis_Mode 3
NO.ofFiles 2
File	D:\2PIData\201512\20151219_1230_1.dat
name	20151219_1230_1
File	D:\2PIData\201512\20151219_1230_2.dat
------------------------------------------------------------------
********************
Analysis mode: 1->getTsum, 2->elastic scattering, 3->ems measurement
4->圆心扫描找最佳角分布的圆心

*****************************/