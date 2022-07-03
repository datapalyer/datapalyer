/*******************
版本2.2更新
1.增加free函数，及时释放掉内存空间，防止溢出。

版本3.0更新
1.增加了p1/2.phi min，max的判断，去掉了p1.x*p2.x<0的条件;
2.修改bug：conincidence()函数，vw判断时，找u。
//2.增加了分角度能谱的输出，分解度解谱的结果。

版本3.1更新：
1.符合增加了能谱角度二维分布图。

版本4.0更新
1.更改了计算电子能量的方式，每两层的能量刻度公式独立
2.完善了distribution_ems结构体的应用。

版本4.1更新
1.修改bug：dis.EPhi1/2对uv，uw，vw都要有
2.加入了uvwTdiff，TsumRes改写了相关写法
***********************/

#include<math.h>
#define uvwdiff 2.0
#define uvwTdiff 10.0
#define TsumRes 10.0
using namespace std;

typedef struct
{
	int *uTdiff,*vTdiff,*wTdiff;  //u,v,w层测量的两个电子时间差
//	int *uTdiff1,*vTdiff1,*wTdiff1;  
//	int *uTdiff2,*vTdiff2,*wTdiff2;
	int **uvTdiff1,**uwTdiff1,**vwTdiff1; //u,v,w时间差的关联，根据剩余一层修正前和修正后
	int **uvTdiff2,**uwTdiff2,**vwTdiff2; 
	int *u12_u11,*u22_u11,*u22_u21,*u12_u21,*v12_v11,*v22_v11; //时间信号检查
	int *v22_v21,*v12_v21,*w12_w11,*w22_w11,*w22_w21,*w12_w21;
//	int *u12_u11_eff,*u22_u11_eff,*u22_u21_eff,*u12_u21_eff,*v12_v11_eff,*v22_v11_eff;
//	int *v22_v21_eff,*v12_v21_eff,*w12_w11_eff,*w22_w11_eff,*w22_w21_eff,*w12_w21_eff;
	int *bindingenergy1,**azimuthangle1; //真符合窗下的能谱和角分布
	int *bindingenergy2,**azimuthangle2; //偶然符合窗下的能谱和角分布
	int **EPhi1,**EPhi2;  //能量角度二维分布
}distribution_ems;

void initialize_dis_ems(distribution_ems &dis,int tdiff_num,int energy_num,int state_num,int angle_num)
{
	initialize_int1(dis.uTdiff,tdiff_num);	initialize_int1(dis.vTdiff,tdiff_num);	initialize_int1(dis.wTdiff,tdiff_num);	
//	initialize_int1(dis.uTdiff1,tdiff_num);	initialize_int1(dis.vTdiff1,tdiff_num);	initialize_int1(dis.wTdiff1,tdiff_num);	
//	initialize_int1(dis.uTdiff2,tdiff_num);	initialize_int1(dis.vTdiff2,tdiff_num);	initialize_int1(dis.wTdiff2,tdiff_num);
	initialize_int2(dis.uvTdiff1,tdiff_num,tdiff_num);
	initialize_int2(dis.uwTdiff1,tdiff_num,tdiff_num);
	initialize_int2(dis.vwTdiff1,tdiff_num,tdiff_num);
	initialize_int2(dis.uvTdiff2,tdiff_num,tdiff_num);
	initialize_int2(dis.uwTdiff2,tdiff_num,tdiff_num);
	initialize_int2(dis.vwTdiff2,tdiff_num,tdiff_num);
	initialize_int1(dis.u12_u11,tdiff_num);	initialize_int1(dis.u22_u11,tdiff_num);
	initialize_int1(dis.u22_u21,tdiff_num);	initialize_int1(dis.u12_u21,tdiff_num);
	initialize_int1(dis.v12_v11,tdiff_num);	initialize_int1(dis.v22_v11,tdiff_num);
	initialize_int1(dis.v22_v21,tdiff_num);	initialize_int1(dis.v12_v21,tdiff_num);
	initialize_int1(dis.w12_w11,tdiff_num);	initialize_int1(dis.w22_w11,tdiff_num);
	initialize_int1(dis.w22_w21,tdiff_num);	initialize_int1(dis.w12_w21,tdiff_num);
	initialize_int1(dis.bindingenergy1,energy_num);
	initialize_int2(dis.azimuthangle1,state_num,angle_num);
	initialize_int1(dis.bindingenergy2,energy_num);
	initialize_int2(dis.azimuthangle2,state_num,angle_num);
	initialize_int2(dis.EPhi1,energy_num,angle_num);
	initialize_int2(dis.EPhi2,energy_num,angle_num);
}

void free_ems(distribution_ems &dis,int n1,int n2,int n3) //n1=Tdiff_num,n2=state_num,n3=energy_num
{
	free(dis.uTdiff);	free(dis.vTdiff);	free(dis.wTdiff);
//	free(dis.uTdiff1);	free(dis.vTdiff1);	free(dis.wTdiff1);
//	free(dis.uTdiff2);	free(dis.vTdiff2);	free(dis.wTdiff2);
	for(int i=0;i<n1;i++)
	{
		free(dis.uvTdiff1[i]);free(dis.uvTdiff2[i]);
		free(dis.uwTdiff1[i]);free(dis.uwTdiff2[i]);
		free(dis.vwTdiff1[i]);free(dis.vwTdiff2[i]);
	}
	free(dis.uvTdiff1);free(dis.uvTdiff2);
	free(dis.uwTdiff1);free(dis.uwTdiff2);
	free(dis.vwTdiff1);free(dis.vwTdiff2);
	free(dis.u12_u11);free(dis.u22_u11);
	free(dis.u22_u21);free(dis.u12_u21);
	free(dis.v12_v11);free(dis.v22_v11);
	free(dis.v22_v21);free(dis.v12_v21);
	free(dis.w12_w11);free(dis.w22_w11);
	free(dis.w22_w21);free(dis.w12_w21);
	for(int i=0;i<n2;i++)
	{
		free(dis.azimuthangle1[i]);
		free(dis.azimuthangle2[i]);
	}
	free(dis.bindingenergy1);free(dis.azimuthangle1);
	free(dis.bindingenergy2);free(dis.azimuthangle2);
	for(int i=0;i<n3;i++)
	{
		free(dis.EPhi1[i]);
		free(dis.EPhi2[i]);
	}
	free(dis.EPhi1);
	free(dis.EPhi2);
}

int conincidence(event &e,particle &p1,particle &p2,parameter param,distribution_ems &dis)
{
	double uTdiff,vTdiff,wTdiff;
	int u_same=0,v_same=0,w_same=0;//判断是否两个粒子同时到达，即两个粒子的Tsum都在正确Tsum范围内

	//u路信号分配，u1[0]+u2[0]->Tsum
	if(e.nhit[u1ch]>=2&&e.nhit[u2ch]>=2)
	{
		fill_double((e.u1[1]-e.u1[0])*0.025,dis.u12_u11,param.tdiff_min,param.tdiff_max,param.tdiff_step);
		fill_double((e.u2[1]-e.u1[0])*0.025,dis.u22_u11,param.tdiff_min,param.tdiff_max,param.tdiff_step);
		fill_double((e.u2[1]-e.u2[0])*0.025,dis.u22_u21,param.tdiff_min,param.tdiff_max,param.tdiff_step);
		fill_double((e.u1[1]-e.u2[0])*0.025,dis.u12_u21,param.tdiff_min,param.tdiff_max,param.tdiff_step);

		if( testTsum(e.u1[0]+e.u2[1],'u',param) && !testTsum(e.u1[1]+e.u2[0],'u',param))	
			swap(e.u2[0],e.u2[1]);
		if( !testTsum(e.u1[0]+e.u2[1],'u',param) && testTsum(e.u1[1]+e.u2[0],'u',param))	
			swap(e.u1[0],e.u1[1]);
		if( testTsum(e.u1[0]+e.u2[1],'u',param) && testTsum(e.u1[1]+e.u2[0],'u',param))
		{
			if( (e.u1[0]+e.u2[1])<(e.u1[1]+e.u2[0]) )
				swap(e.u2[0],e.u2[1]); 
			else
				swap(e.u1[0],e.u1[1]);
			u_same=1;
		}
 //     if( testTsum(e.u1[0]+e.u2[0],'u',param))
//			fill_double( ((e.u1[1]+e.u2[1])-(e.u1[0]+e.u2[0]))/2,dis.uTdiff1,param.tdiff_min,param.tdiff_max,param.tdiff_step);
	}
	if(e.nhit[v1ch]>=2&&e.nhit[v2ch]>=2)
	{
		fill_double((e.v1[1]-e.v1[0])*0.025,dis.v12_v11,param.tdiff_min,param.tdiff_max,param.tdiff_step);
		fill_double((e.v2[1]-e.v1[0])*0.025,dis.v22_v11,param.tdiff_min,param.tdiff_max,param.tdiff_step);
		fill_double((e.v2[1]-e.v2[0])*0.025,dis.v22_v21,param.tdiff_min,param.tdiff_max,param.tdiff_step);
		fill_double((e.v1[1]-e.v2[0])*0.025,dis.v12_v21,param.tdiff_min,param.tdiff_max,param.tdiff_step);
		if( testTsum(e.v1[0]+e.v2[1],'v',param) && !testTsum(e.v1[1]+e.v2[0],'v',param))	
			swap(e.v2[0],e.v2[1]);
		if( !testTsum(e.v1[0]+e.v2[1],'v',param) && testTsum(e.v1[1]+e.v2[0],'v',param))	
			swap(e.v1[0],e.v1[1]);
		if( testTsum(e.v1[0]+e.v2[1],'v',param) && testTsum(e.v1[1]+e.v2[0],'v',param))
		{
			if( (e.v1[0]+e.v2[1])<(e.v1[1]+e.v2[0]) )
				swap(e.v2[0],e.v2[1]); 
			else
				swap(e.v1[0],e.v1[1]);
			v_same=1;
		}
//        if( testTsum(e.v1[0]+e.v2[0],'v',param))
//			fill_double( ((e.v1[1]+e.v2[1])-(e.v1[0]+e.v2[0]))/2,dis.vTdiff1,param.tdiff_min,param.tdiff_max,param.tdiff_step);
	}
	if(e.nhit[w1ch]>=2&&e.nhit[w2ch]>=2)
	{
		fill_double((e.w1[1]-e.w1[0])*0.025,dis.w12_w11,param.tdiff_min,param.tdiff_max,param.tdiff_step);
		fill_double((e.w2[1]-e.w1[0])*0.025,dis.w22_w11,param.tdiff_min,param.tdiff_max,param.tdiff_step);
		fill_double((e.w2[1]-e.w2[0])*0.025,dis.w22_w21,param.tdiff_min,param.tdiff_max,param.tdiff_step);
		fill_double((e.w1[1]-e.w2[0])*0.025,dis.w12_w21,param.tdiff_min,param.tdiff_max,param.tdiff_step);

		if( testTsum(e.w1[0]+e.w2[1],'w',param) && !testTsum(e.w1[1]+e.w2[0],'w',param))
			swap(e.w2[0],e.w2[1]);
		if( !testTsum(e.w1[0]+e.w2[1],'w',param) && testTsum(e.w1[1]+e.w2[0],'w',param))
			swap(e.w1[0],e.w1[1]);
		if( testTsum(e.w1[0]+e.w2[1],'w',param) && testTsum(e.w1[1]+e.w2[0],'w',param))	
		{
			if( (e.w1[0]+e.w2[1])<(e.w1[1]+e.w2[0]) )
				swap(e.w2[0],e.w2[1]); 
			else
				swap(e.w1[0],e.w1[1]);
			w_same=1;
		}
 //       if( testTsum(e.w1[0]+e.w2[0],'w',param))
//			fill_double( ((e.w1[1]+e.w2[1])-(e.w1[0]+e.w2[0]))/2,dis.wTdiff1,param.tdiff_min,param.tdiff_max,param.tdiff_step);
	}

	//两层组合判断双击
	//u,v两层有双击
	if(e.nhit[u1ch]>=2&&e.nhit[u2ch]>=2&&e.nhit[v1ch]>=2&&e.nhit[v2ch]>=2)
	{
		initialize_particle(p1);
		initialize_particle(p2);
		e.uTsum=(e.u1[0]+e.u2[0]);
		e.vTsum=(e.v1[0]+e.v2[0]);
		if( testTsum(e.uTsum,'u',param) && testTsum(e.vTsum,'v',param) )
		{
			p1.uTsum=e.uTsum*0.025;
			p2.uTsum=(e.u1[1]+e.u2[1])*0.025;
			p1.vTsum=e.vTsum*0.025;
			p2.vTsum=(e.v1[1]+e.v2[1])*0.025;
			uTdiff=(p2.uTsum-p1.uTsum)/2+1;
			vTdiff=(p2.vTsum-p1.vTsum)/2;
			p1.u=(e.u1[0]-e.u2[0])*0.025*param.du;
			p2.u=(e.u1[1]-e.u2[1])*0.025*param.du;
			p1.uflag=1;
			p2.uflag=1;
			p1.v=(e.v1[0]-e.v2[0])*0.025*param.dv;
			p2.v=(e.v1[1]-e.v2[1])*0.025*param.dv;
			p1.vflag=1;
			p2.vflag=1;

			fill2D_double(uTdiff,vTdiff,dis.uvTdiff1,param.tdiff_min,param.tdiff_max,param.tdiff_step,param.tdiff_min,param.tdiff_max,param.tdiff_step);

			if(uTdiff<TsumRes||vTdiff<TsumRes)
			{
				double w_test;
				int w_flag=0;
				for(int i=0;i<e.nhit[w1ch]&&w_flag==0;i++)
				{
					for(int j=0;j<e.nhit[w2ch]&&w_flag==0;j++)
					{
						if(testTsum(e.w1[i]+e.w2[j],'w',param))
						{
							w_test=(e.w1[i]-e.w2[j])*0.025*param.dw+param.ow;
							w_flag=1;
							break;
						}
					}
				}
				if(w_flag)
				{
					if(fabs(p1.u-p1.v-w_test)<uvwdiff || fabs(p2.u-p2.v-w_test)<uvwdiff)
						return 1;
					else if(fabs(p1.u-p2.v-w_test)<uvwdiff || fabs(p2.u-p1.v-w_test)<uvwdiff)
					{
						swap_double(p1.v,p2.v);
						swap_double(p1.vTsum,p2.vTsum);
						vTdiff=(p2.vTsum-p1.vTsum)/2;
						return 1;
					}
				}
			}//END of if(u_same==1||v_same==1)

			if(fabs(uTdiff-vTdiff)<uvwTdiff)
				return 1;
		}//END of if( testTsum(e.uTsum,'u',param) && testTsum(e.vTsum,'v',param) )
	}//END of if(e.nhit[u1ch]>=2&&e.nhit[u2ch]>=2&&e.nhit[v1ch]>=2&&e.nhit[v2ch]>=2)

	if(e.nhit[u1ch]>=2&&e.nhit[u2ch]>=2&&e.nhit[w1ch]>=2&&e.nhit[w2ch]>=2)
	{
        initialize_particle(p1);
    	initialize_particle(p2);
		e.uTsum=(e.u1[0]+e.u2[0]);
		e.wTsum=(e.w1[0]+e.w2[0]);
		if( testTsum(e.uTsum,'u',param) && testTsum(e.wTsum,'w',param) )
		{
			p1.uTsum=e.uTsum*0.025;
			p2.uTsum=(e.u1[1]+e.u2[1])*0.025;
			p1.wTsum=e.wTsum*0.025;
			p2.wTsum=(e.w1[1]+e.w2[1])*0.025;
			uTdiff=(p2.uTsum-p1.uTsum)/2+1;
			wTdiff=(p2.wTsum-p1.wTsum)/2+1;
			p1.u=(e.u1[0]-e.u2[0])*0.025*param.du;
			p2.u=(e.u1[1]-e.u2[1])*0.025*param.du;
			p1.uflag=1;
			p2.uflag=1;
			p1.w=(e.w1[0]-e.w2[0])*0.025*param.dw+param.ow;
			p2.w=(e.w1[1]-e.w2[1])*0.025*param.dw+param.ow;
			p1.wflag=1;
			p2.wflag=1;

			fill2D_double(uTdiff,wTdiff,dis.uwTdiff1,param.tdiff_min,param.tdiff_max,param.tdiff_step,param.tdiff_min,param.tdiff_max,param.tdiff_step);

			if(uTdiff<TsumRes||wTdiff<TsumRes)
			{
				double v_test;
				int v_flag=0;
				for(int i=0;i<e.nhit[v1ch]&&v_flag==0;i++)
				{
					for(int j=0;j<e.nhit[v2ch]&&v_flag==0;j++)
					{
						if(testTsum(e.v1[i]+e.v2[j],'v',param))
						{
							v_test=(e.v1[i]-e.v2[j])*0.025*param.dv;
							v_flag=1;
							break;
						}
					}
				}
				if(v_flag)
				{
					if(fabs(p1.u-p1.w-v_test)<uvwdiff || fabs(p2.u-p2.w-v_test)<uvwdiff)
						return 1;
					else if(fabs(p1.u-p2.w-v_test)<uvwdiff || fabs(p2.u-p1.w-v_test)<uvwdiff)
					{
						swap_double(p1.w,p2.w);
						swap_double(p1.wTsum,p2.wTsum);
						wTdiff=(p2.wTsum-p1.wTsum)/2+1;
						return 1;
					}
				}
			}//END of if(u_same==1||w_same==1)

			if(fabs(uTdiff-wTdiff)<uvwTdiff)
				return 1;
		}//END of if( testTsum(e.uTsum,'u',param) && testTsum(e.wTsum,'w',param) )
	}//END of if(e.nhit[u1ch]>=2&&e.nhit[u2ch]>=2&&e.nhit[w1ch]>=2&&e.nhit[w2ch]>=2)
	
	
    if(e.nhit[v1ch]>=2&&e.nhit[v2ch]>=2&&e.nhit[w1ch]>=2&&e.nhit[w2ch]>=2)
	{
        initialize_particle(p1);
    	initialize_particle(p2);
		e.vTsum=(e.v1[0]+e.v2[0]);
		e.wTsum=(e.w1[0]+e.w2[0]);
		if( testTsum(e.vTsum,'v',param) && testTsum(e.wTsum,'w',param) )
		{
            p1.vTsum=e.vTsum*0.025;
			p2.vTsum=(e.v1[1]+e.v2[1])*0.025;
			p1.wTsum=e.wTsum*0.025;
			p2.wTsum=(e.w1[1]+e.w2[1])*0.025;
			vTdiff=(p2.vTsum-p1.vTsum)/2;
			wTdiff=(p2.wTsum-p1.wTsum)/2+1;
			p1.v=(e.v1[0]-e.v2[0])*0.025*param.dv;
			p2.v=(e.v1[1]-e.v2[1])*0.025*param.dv;
			p1.vflag=1;
			p2.vflag=1;
			p1.w=(e.w1[0]-e.w2[0])*0.025*param.dw+param.ow;
			p2.w=(e.w1[1]-e.w2[1])*0.025*param.dw+param.ow;
			p1.wflag=1;
			p2.wflag=1;

			fill2D_double(vTdiff,wTdiff,dis.vwTdiff1,param.tdiff_min,param.tdiff_max,param.tdiff_step,param.tdiff_min,param.tdiff_max,param.tdiff_step);

			if(vTdiff<TsumRes||wTdiff<TsumRes)
			{
				double u_test;
				int u_flag=0;
				for(int i=0;i<e.nhit[u1ch]&&u_flag==0;i++)
				{
					for(int j=0;j<e.nhit[u2ch]&&u_flag==0;j++)
					{
						if(testTsum(e.u1[i]+e.u2[j],'u',param))
						{
							u_test=(e.u1[i]-e.u2[j])*0.025*param.du;
							u_flag=1;
							break;
						}
					}
				}
				if(u_flag)
				{
					if(fabs(u_test-p1.v-p1.w)<uvwdiff || fabs(u_test-p2.v-p2.w)<uvwdiff)
						return 1;
					else if(fabs(u_test-p1.v-p2.w)<uvwdiff || fabs(u_test-p2.v-p1.w)<uvwdiff)
					{
						swap_double(p1.w,p2.w);
						swap_double(p1.wTsum,p2.wTsum);
						wTdiff=(p2.wTsum-p1.wTsum)/2+1;
						return 1;
					}
				}
			}//END of if(v_same==1||w_same==1)

			if(fabs(vTdiff-wTdiff)<uvwTdiff)
				return 1;

		}//END of if( testTsum(e.uTsum,'v',param) && testTsum(e.wTsum,'w',param) )
	}//END of if(e.nhit[v1ch]>=2&&e.nhit[v2ch]>=2&&e.nhit[w1ch]>=2&&e.nhit[w2ch]>=2)

	return 0;
}

int ems(parameter param)
{
	int time_num=(int)((param.time_max-param.time_min)/param.time_step+1);
	int position_num=(int)((param.position_max-param.position_min)/param.position_step+1);
	int radius_num=(int)((param.radius_max-param.radius_min)/param.radius_step+1);
	int phi_num=(int)((param.phi_max-param.phi_min)/param.phi_step+0.99);
	int tdiff_num=(int)((param.tdiff_max-param.tdiff_min)/param.tdiff_step+1);
	int energy_num=(int)((param.energy_max-param.energy_min)/param.energy_step+1);
	int angle_num=(int)((param.angle_max-param.angle_min)/param.angle_step+1);
	event e; 
	particle p1,p2;

	distribution_ems dis_sum;
	initialize_dis_ems(dis_sum,tdiff_num,energy_num,param.state_num,angle_num);

	int num_file;
	for(num_file=0;num_file<param.num_files;num_file++)
	{
		FILE *fp_in;
		fp_in=fopen(param.file[num_file],"rb");
		if(fp_in==NULL)
		{
			cout<<"Input file error!"<<endl;
			return 0;
		}

		FILE *fp_timespectra,*fp_bes,*fp_angledis,*fp_timecor,*fp_EPhi;
		char timespectra[MAX_LINE],BES[MAX_LINE],angledis[MAX_LINE],time_correlation[MAX_LINE],EPhi[MAX_LINE];
		sprintf(timespectra,"%s%s_timespectra.txt",param.workspace,param.name[num_file]);  //时间谱
		fp_timespectra=fopen(timespectra,"w+");
		sprintf(time_correlation,"%s%s_timecorrelation.txt",param.workspace,param.name[num_file]);  //三层时间差关联图
		fp_timecor=fopen(time_correlation,"w+");
		sprintf(BES,"%s%s_BES.txt",param.workspace,param.name[num_file]);  //能谱
		fp_bes=fopen(BES,"w+");
		sprintf(angledis,"%s%s_AngleDistribution.txt",param.workspace,param.name[num_file]);  //方位角分布
		fp_angledis=fopen(angledis,"w+");
		sprintf(EPhi,"%s%s_EPhi2D.txt",param.workspace,param.name[num_file]);  //能量角度二维谱
		fp_EPhi=fopen(EPhi,"w+");
		
		//初始化统计分布
		distribution_ems dis;
		initialize_dis_ems(dis,tdiff_num,energy_num,param.state_num,angle_num);

		double bindingenergy;
		double azimuthangle;
		int count_total=0,count_eff1=0,count_eff2=0;
		int count3=0,count4=0;
		while(!feof(fp_in))
		{
			initialize_event(e);
			initialize_particle(p1);
			initialize_particle(p2);
			readevent(e,fp_in);
			count_total++;
			if(count_total%500000==0)
				printf("total eff1 eff2: %d %d %d %d %d\n",count_total,count_eff1,count_eff2,count3,count4);
			if(!completesignal(e.pattern[1])) continue;	
			count_eff1++;
			if(!conincidence(e,p1,p2,param,dis))	continue;
			count_eff2++;
			
			if(p1.uflag*p2.uflag*p1.vflag*p2.vflag)
			{
				getposition(p1,param);
				getposition(p2,param);
				if(p1.phi>p2.phi)
					swap_particle(p1,p2);

				if(p1.angle<param.Phi1_min||p1.angle>param.Phi1_max || p2.angle<param.Phi2_min||p2.angle>param.Phi2_max)
					continue;

				double uTdiff=(p2.uTsum-p1.uTsum)/2+1;
				double vTdiff=(p2.vTsum-p1.vTsum)/2;

				getenergy(p1,param);
				getenergy(p2,param);
				if(p1.energy<param.E12_min||p1.energy>param.E12max||p2.energy<param.E12_min||p2.energy>param.E12max)
					continue;
				bindingenergy=getIP(p1,p2,param)+param.E0_modify[num_file];
				azimuthangle=getAzimuth(p1,p2,param);

				fill2D_double(uTdiff,vTdiff,dis.uvTdiff2,param.tdiff_min,param.tdiff_max,param.tdiff_step,param.tdiff_min,param.tdiff_max,param.tdiff_step);
				fill_double(uTdiff,dis.uTdiff,param.tdiff_min,param.tdiff_max,param.tdiff_step);
				fill_double(vTdiff,dis.vTdiff,param.tdiff_min,param.tdiff_max,param.tdiff_step);
				count3++;
				if(truecoincidence(uTdiff,param)&&truecoincidence(vTdiff,param))
				{
					count4++;
					fill_double(bindingenergy,dis.bindingenergy1,param.energy_min,param.energy_max,param.energy_step);
					for(int istate=0;istate<param.state_num;istate++)
					{
						if(bindingenergy>=param.state_energy_min[istate]&&bindingenergy<=param.state_energy_max[istate])
							fill_double(azimuthangle,dis.azimuthangle1[istate],param.angle_min,param.angle_max,param.angle_step);

					}
					fill2D_double(bindingenergy,azimuthangle,dis.EPhi1,param.energy_min,param.energy_max,param.energy_step,param.angle_min,param.angle_max,param.angle_step);
				}
				if(accidentalcoincidence(uTdiff,param)&&accidentalcoincidence(vTdiff,param))
				{
					fill_double(bindingenergy,dis.bindingenergy2,param.energy_min,param.energy_max,param.energy_step);
					for(int istate=0;istate<param.state_num;istate++)
					{
						if(bindingenergy>=param.state_energy_min[istate]&&bindingenergy<=param.state_energy_max[istate])
							fill_double(azimuthangle,dis.azimuthangle2[istate],param.angle_min,param.angle_max,param.angle_step);

					}
					fill2D_double(bindingenergy,azimuthangle,dis.EPhi2,param.energy_min,param.energy_max,param.energy_step,param.angle_min,param.angle_max,param.angle_step);
				}
			}
			else if(p1.uflag*p2.uflag*p1.wflag*p2.wflag)
			{
                getposition(p1,param);
				getposition(p2,param);
				if(p1.phi>p2.phi)
					swap_particle(p1,p2);
				if(p1.angle<param.Phi1_min||p1.angle>param.Phi1_max || p2.angle<param.Phi2_min||p2.angle>param.Phi2_max)
					continue;

				double uTdiff=(p2.uTsum-p1.uTsum)/2+1;
				double wTdiff=(p2.wTsum-p1.wTsum)/2+1;
				getenergy(p1,param);
				getenergy(p2,param);
				if(p1.energy<param.E12_min||p1.energy>param.E12max||p2.energy<param.E12_min||p2.energy>param.E12max)
					continue;
				bindingenergy=getIP(p1,p2,param)+param.E0_modify[num_file];
				azimuthangle=getAzimuth(p1,p2,param);

				fill2D_double(uTdiff,wTdiff,dis.uwTdiff2,param.tdiff_min,param.tdiff_max,param.tdiff_step,param.tdiff_min,param.tdiff_max,param.tdiff_step);
				fill_double(uTdiff,dis.uTdiff,param.tdiff_min,param.tdiff_max,param.tdiff_step);
				fill_double(wTdiff,dis.wTdiff,param.tdiff_min,param.tdiff_max,param.tdiff_step);
				count3++;
				if(truecoincidence(uTdiff,param)&&truecoincidence(wTdiff,param))
				{
					count4++;
					fill_double(bindingenergy,dis.bindingenergy1,param.energy_min,param.energy_max,param.energy_step);
					for(int istate=0;istate<param.state_num;istate++)
					{
						if(bindingenergy>=param.state_energy_min[istate]&&bindingenergy<=param.state_energy_max[istate])
							fill_double(azimuthangle,dis.azimuthangle1[istate],param.angle_min,param.angle_max,param.angle_step);

					}
					fill2D_double(bindingenergy,azimuthangle,dis.EPhi1,param.energy_min,param.energy_max,param.energy_step,param.angle_min,param.angle_max,param.angle_step);

				}
				if(accidentalcoincidence(uTdiff,param)&&accidentalcoincidence(wTdiff,param))
				{
					fill_double(bindingenergy,dis.bindingenergy2,param.energy_min,param.energy_max,param.energy_step);
					for(int istate=0;istate<param.state_num;istate++)
					{
						if(bindingenergy>=param.state_energy_min[istate]&&bindingenergy<=param.state_energy_max[istate])
							fill_double(azimuthangle,dis.azimuthangle2[istate],param.angle_min,param.angle_max,param.angle_step);

					}
					fill2D_double(bindingenergy,azimuthangle,dis.EPhi2,param.energy_min,param.energy_max,param.energy_step,param.angle_min,param.angle_max,param.angle_step);
				}
			}
			else if(p1.vflag*p2.vflag*p1.wflag*p2.wflag)
			{
				getposition(p1,param);
				getposition(p2,param);
				if(p1.phi>p2.phi)
					swap_particle(p1,p2);
				if(p1.angle<param.Phi1_min||p1.angle>param.Phi1_max || p2.angle<param.Phi2_min||p2.angle>param.Phi2_max)
					continue;
				double vTdiff=(p2.vTsum-p1.vTsum)/2;
				double wTdiff=(p2.wTsum-p1.wTsum)/2+1;
				getenergy(p1,param);
				getenergy(p2,param);
				if(p1.energy<param.E12_min||p1.energy>param.E12max||p2.energy<param.E12_min||p2.energy>param.E12max)
					continue;
				bindingenergy=getIP(p1,p2,param)+param.E0_modify[num_file];
				azimuthangle=getAzimuth(p1,p2,param);

				fill2D_double(vTdiff,wTdiff,dis.vwTdiff2,param.tdiff_min,param.tdiff_max,param.tdiff_step,param.tdiff_min,param.tdiff_max,param.tdiff_step);
				fill_double(vTdiff,dis.vTdiff,param.tdiff_min,param.tdiff_max,param.tdiff_step);
				fill_double(wTdiff,dis.wTdiff,param.tdiff_min,param.tdiff_max,param.tdiff_step);
				count3++;
				if(truecoincidence(vTdiff,param)&&truecoincidence(wTdiff,param))
				{
					count4++;
					fill_double(bindingenergy,dis.bindingenergy1,param.energy_min,param.energy_max,param.energy_step);
					for(int istate=0;istate<param.state_num;istate++)
					{
						if(bindingenergy>=param.state_energy_min[istate]&&bindingenergy<=param.state_energy_max[istate])
							fill_double(azimuthangle,dis.azimuthangle1[istate],param.angle_min,param.angle_max,param.angle_step);

					}
					fill2D_double(bindingenergy,azimuthangle,dis.EPhi1,param.energy_min,param.energy_max,param.energy_step,param.angle_min,param.angle_max,param.angle_step);

				}
				if(accidentalcoincidence(vTdiff,param)&&accidentalcoincidence(wTdiff,param))
				{
					fill_double(bindingenergy,dis.bindingenergy2,param.energy_min,param.energy_max,param.energy_step);
					for(int istate=0;istate<param.state_num;istate++)
					{
						if(bindingenergy>=param.state_energy_min[istate]&&bindingenergy<=param.state_energy_max[istate])
							fill_double(azimuthangle,dis.azimuthangle2[istate],param.angle_min,param.angle_max,param.angle_step);

					}
					fill2D_double(bindingenergy,azimuthangle,dis.EPhi2,param.energy_min,param.energy_max,param.energy_step,param.angle_min,param.angle_max,param.angle_step);
				}

			}//End of else if(p1.vflag*p2.vflag*p1.wflag*p2.wflag)

		}//End of while(!feof(fp_in))
		printf("total eff1 eff2: %d %d %d %d %d\n",count_total,count_eff1,count_eff2,count3,count4);

		fprintf(fp_timespectra,"tdiff uTdiff vTdiff wTdiff u12_u11_dis u22_u11_dis u22_u21_dis u12_u21_dis v12_v11_dis v22_v11_dis v22_v21_dis v12_v21_dis w12_w11_dis w22_w11_dis w22_w21_dis w12_w21_dis\n");
		for(int i=0;i<tdiff_num;i++)
		{
			double diff=param.tdiff_min+(i+0.5)*param.tdiff_step;
			fprintf(fp_timespectra,"%.2f %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",diff,dis.uTdiff[i],dis.vTdiff[i],dis.wTdiff[i],
				dis.u12_u11[i],dis.u22_u11[i],dis.u22_u21[i],dis.u12_u21[i],dis.v12_v11[i],dis.v22_v11[i],dis.v22_v21[i],dis.v12_v21[i],
				dis.w12_w11[i],dis.w22_w11[i],dis.w22_w21[i],dis.w12_w21[i]);
		}

		fprintf(fp_timecor,"tdiff1 tdiff2 uvdiff1 uwdiff1 vwdiff1 uvdiff2 uwdiff2 vwdiff2\n");
		for(int i=0;i<tdiff_num;i++)
		{
			for(int j=0;j<tdiff_num;j++)
			{
				double diff1=param.tdiff_min+(i+0.5)*param.tdiff_step;
				double diff2=param.tdiff_min+(j+0.5)*param.tdiff_step;
				fprintf(fp_timecor,"%.2f %.2f %d %d %d %d %d %d\n",diff1,diff2,dis.uvTdiff1[i][j],
					dis.uwTdiff1[i][j],dis.vwTdiff1[i][j],dis.uvTdiff2[i][j],dis.uwTdiff2[i][j],dis.vwTdiff2[i][j]);
			}
		}

		fprintf(fp_bes,"energy truecoin accidentcoin\n");
		for(int i=0;i<energy_num;i++)
		{
			fprintf(fp_bes,"%.2f %d %d\n",param.energy_min+(i+0.5)*param.energy_step,dis.bindingenergy1[i],dis.bindingenergy2[i]);
		}

		fprintf(fp_angledis,"azimuth_angle truecoin accidentcoin\n");
		for(int i=0;i<angle_num;i++)
		{
			fprintf(fp_angledis,"%.2f ",param.angle_min+(i+0.5)*param.angle_step);
			for(int istate=0;istate<param.state_num;istate++)
			{
				fprintf(fp_angledis,"%d %d ",dis.azimuthangle1[istate][i],dis.azimuthangle2[istate][i]);
			}
			fprintf(fp_angledis,"\n");
		}

		fprintf(fp_EPhi,"energy angle truecoin accidentcoin");
		for(int i=0;i<energy_num;i++)
		{
			for(int j=0;j<angle_num;j++)
			{
				fprintf(fp_EPhi,"%.2f %.2f %d %d\n",param.energy_min+(i+0.5)*param.energy_step,
					param.angle_min+(j+0.5)*param.angle_step,dis.EPhi1[i][j],dis.EPhi2[i][j]);
			}
		}

		fclose(fp_in);
		fclose(fp_timespectra);
		fclose(fp_bes);
		fclose(fp_angledis);
		fclose(fp_timecor);
		fclose(fp_EPhi);
		cout<<param.file[num_file]<<" is over"<<endl;

		for(int i=0;i<tdiff_num;i++)
		{
			dis_sum.uTdiff[i]+=dis.uTdiff[i];
			dis_sum.vTdiff[i]+=dis.vTdiff[i];
			dis_sum.wTdiff[i]+=dis.wTdiff[i];
			dis_sum.u12_u11[i]+=dis.u12_u11[i];
			dis_sum.u22_u11[i]+=dis.u22_u11[i];
			dis_sum.u22_u21[i]+=dis.u22_u21[i];
			dis_sum.u12_u21[i]+=dis.u12_u21[i];
			dis_sum.v12_v11[i]+=dis.v12_v11[i];
			dis_sum.v22_v11[i]+=dis.v22_v11[i];
			dis_sum.v22_v21[i]+=dis.v22_v21[i];
			dis_sum.v12_v21[i]+=dis.v12_v21[i];
			dis_sum.w12_w11[i]+=dis.w12_w11[i];
			dis_sum.w22_w11[i]+=dis.w22_w11[i];
			dis_sum.w22_w21[i]+=dis.w22_w21[i];
			dis_sum.w12_w21[i]+=dis.w12_w21[i];
			for(int j=0;j<tdiff_num;j++)
			{
				dis_sum.uvTdiff1[i][j]+=dis.uvTdiff1[i][j];
				dis_sum.uwTdiff1[i][j]+=dis.uwTdiff1[i][j];
				dis_sum.vwTdiff1[i][j]+=dis.vwTdiff1[i][j];
				dis_sum.uvTdiff2[i][j]+=dis.uvTdiff2[i][j];
				dis_sum.uwTdiff2[i][j]+=dis.uwTdiff2[i][j];
				dis_sum.vwTdiff2[i][j]+=dis.vwTdiff2[i][j];
			}
		}
//cout<<"1"<<endl;
		for(int i=0;i<energy_num;i++)
		{
			dis_sum.bindingenergy1[i]+=dis.bindingenergy1[i];
			dis_sum.bindingenergy2[i]+=dis.bindingenergy2[i];
		}
//		cout<<"2"<<endl;
		for(int i=0;i<angle_num;i++)
		{
			for(int istate=0;istate<param.state_num;istate++)
			{
				dis_sum.azimuthangle1[istate][i]+=dis.azimuthangle1[istate][i];
				dis_sum.azimuthangle2[istate][i]+=dis.azimuthangle2[istate][i];
			}
		}
//		cout<<"3"<<endl;
		for(int i=0;i<energy_num;i++)
		{
			for(int j=0;j<angle_num;j++)
			{
				dis_sum.EPhi1[i][j]+=dis.EPhi1[i][j];
				dis_sum.EPhi2[i][j]+=dis.EPhi2[i][j];
			}
		}
//		cout<<"4"<<endl;
		free_ems(dis,tdiff_num,param.state_num,energy_num);

	}//End of for(num_file=0;num_file<param.num_files;num_file++)

	if(param.num_files>1)
	{
		FILE *fp_timespectra,*fp_bes,*fp_angledis,*fp_timecor,*fp_EPhi;
		char timespectra[MAX_LINE],BES[MAX_LINE],angledis[MAX_LINE],time_correlation[MAX_LINE],EPhi[MAX_LINE];
		sprintf(timespectra,"%s%s_sum_timespectra.txt",param.workspace,param.name[0]);
		sprintf(time_correlation,"%s%s_sum_timecorrelation.txt",param.workspace,param.name[0]);
		sprintf(BES,"%s%s_sum_BES.txt",param.workspace,param.name[0]);
		sprintf(angledis,"%s%s_sum_AngleDistribution.txt",param.workspace,param.name[0]);
		fp_timespectra=fopen(timespectra,"w+");
		fp_bes=fopen(BES,"w+");
		fp_angledis=fopen(angledis,"w+");
		fp_timecor=fopen(time_correlation,"w+");
		sprintf(EPhi,"%s%s_sum_EPhi2D.txt",param.workspace,param.name[0]);  //能量角度二维谱
		fp_EPhi=fopen(EPhi,"w+");

		fprintf(fp_timespectra,"tdiff uTdiff vTdiff wTdiff u12_u11_dis u22_u11_dis u22_u21_dis u12_u21_dis v12_v11_dis v22_v11_dis v22_v21_dis v12_v21_dis w12_w11_dis w22_w11_dis w22_w21_dis w12_w21_dis\n");
		for(int i=0;i<tdiff_num;i++)
		{
			double diff=param.tdiff_min+(i+0.5)*param.tdiff_step;
			fprintf(fp_timespectra,"%.2f %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",diff,dis_sum.uTdiff[i],dis_sum.vTdiff[i],dis_sum.wTdiff[i],
				dis_sum.u12_u11[i],dis_sum.u22_u11[i],dis_sum.u22_u21[i],dis_sum.u12_u21[i],dis_sum.v12_v11[i],dis_sum.v22_v11[i],dis_sum.v22_v21[i],
				dis_sum.v12_v21[i],dis_sum.w12_w11[i],dis_sum.w22_w11[i],dis_sum.w22_w21[i],dis_sum.w12_w21[i]);
		}

		fprintf(fp_timecor,"tdiff1 tdiff2 uvdiff1 uwdiff1 vwdiff1 uvdiff2 uwdiff2 vwdiff2\n");
		for(int i=0;i<tdiff_num;i++)
		{
			for(int j=0;j<tdiff_num;j++)
			{
				double diff1=param.tdiff_min+(i+0.5)*param.tdiff_step;
				double diff2=param.tdiff_min+(j+0.5)*param.tdiff_step;
				fprintf(fp_timecor,"%.2f %.2f %d %d %d %d %d %d\n",diff1,diff2,dis_sum.uvTdiff1[i][j],
					dis_sum.uwTdiff1[i][j],dis_sum.vwTdiff1[i][j],dis_sum.uvTdiff2[i][j],dis_sum.uwTdiff2[i][j],dis_sum.vwTdiff2[i][j]);
			}
		}

		fprintf(fp_bes,"energy truecoin accidentcoin\n");
		for(int i=0;i<energy_num;i++)
		{
			fprintf(fp_bes,"%.2f %d %d\n",param.energy_min+(i+0.5)*param.energy_step,dis_sum.bindingenergy1[i],dis_sum.bindingenergy2[i]);
		}

		fprintf(fp_angledis,"azimuth_angle truecoin accidentcoin\n");
		for(int i=0;i<angle_num;i++)
		{
			fprintf(fp_angledis,"%.2f ",param.angle_min+(i+0.5)*param.angle_step);
			for(int istate=0;istate<param.state_num;istate++)
			{
				fprintf(fp_angledis,"%d %d ",dis_sum.azimuthangle1[istate][i],dis_sum.azimuthangle2[istate][i]);
			}
			fprintf(fp_angledis,"\n");
		}

		fprintf(fp_EPhi,"energy angle truecoin accidentcoin");
		for(int i=0;i<energy_num;i++)
		{
			for(int j=0;j<angle_num;j++)
			{
				fprintf(fp_EPhi,"%.2f %.2f %d %d\n",param.energy_min+(i+0.5)*param.energy_step,
					param.angle_min+(j+0.5)*param.angle_step,dis_sum.EPhi1[i][j],dis_sum.EPhi2[i][j]);
			}
		}

		fclose(fp_timespectra);
		fclose(fp_bes);
		fclose(fp_angledis);
		fclose(fp_timecor);
		fclose(fp_EPhi);
	}

	free_ems(dis_sum,tdiff_num,param.state_num,energy_num);
	return 1;
}

typedef struct
{
	int *energy1,*energy2;
	int *angle;
	int *auv,*auw,*avw;
	int *auv_uw,*auv_vw,*auw_vw;
}distribution_noncoin;

int noncoin(parameter param)
{
	double energy_min=param.E12_min-5;
	double energy_max=param.E12max+5;
	double energy_step=0.1;
	int energy_num=(int)((energy_max-energy_min)/energy_step+1);
	double angle_step=param.phi_step;
	int phi_num=(int)(360.0/angle_step+1);
	double adiff_min,adiff_max,adiff_step;
	adiff_min=-9.9;
	adiff_max=9.9;
	adiff_step=0.2;
	int adiff_num=(int)((adiff_max-adiff_min)/adiff_step+1);

	distribution_noncoin dis_sum;
	dis_sum.energy1=(int *)malloc(sizeof(int)*energy_num);
	dis_sum.energy2=(int *)malloc(sizeof(int)*energy_num);
	dis_sum.angle=(int *)malloc(sizeof(int)*phi_num);
	dis_sum.auv=(int *)malloc(sizeof(int)*phi_num);
	dis_sum.auw=(int *)malloc(sizeof(int)*phi_num);
	dis_sum.avw=(int *)malloc(sizeof(int)*phi_num);
	dis_sum.auv_uw=(int *)malloc(sizeof(int)*adiff_num);
	dis_sum.auv_vw=(int *)malloc(sizeof(int)*adiff_num);
	dis_sum.auw_vw=(int *)malloc(sizeof(int)*adiff_num);
	for(int i=0;i<energy_num;i++)
	{	
		dis_sum.energy1[i]=0;
		dis_sum.energy2[i]=0;
	}
	for(int i=0;i<phi_num;i++) 
	{
		dis_sum.angle[i]=0;
		dis_sum.auv[i]=0;dis_sum.auw[i]=0;dis_sum.avw[i]=0;
	}
	for(int i=0;i<adiff_num;i++) 
	{
		dis_sum.auv_uw[i]=0;
		dis_sum.auv_vw[i]=0;
		dis_sum.auw_vw[i]=0;
	}

	event e; 
	particle p1;

	int num_file;
	for(num_file=0;num_file<param.num_files;num_file++)
	{
		FILE *fp_in;
		fp_in=fopen(param.file[num_file],"rb");
		if(fp_in==NULL)
		{
			cout<<"Input file error!"<<endl;
			return 0;
		}

		FILE *fp_energy,*fp_angle,*fp_adiff;
		char out_energy[MAX_LINE],out_angle[MAX_LINE],out_adiff[MAX_LINE];	
		sprintf(out_energy,"%s%s_energy.txt",param.workspace,param.name[num_file]);
		sprintf(out_angle,"%s%s_angle.txt",param.workspace,param.name[num_file]);	
		sprintf(out_adiff,"%s%s_angle_uvwdiff.txt",param.workspace,param.name[num_file]);	
		fp_energy=fopen(out_energy,"w+");
		fp_angle=fopen(out_angle,"w+");
		fp_adiff=fopen(out_adiff,"w+");

		distribution_noncoin dis;
		dis.energy1=(int *)malloc(sizeof(int)*energy_num);
		dis.energy2=(int *)malloc(sizeof(int)*energy_num);
		dis.angle=(int *)malloc(sizeof(int)*phi_num);
		dis.auv=(int *)malloc(sizeof(int)*phi_num);
		dis.auw=(int *)malloc(sizeof(int)*phi_num);
		dis.avw=(int *)malloc(sizeof(int)*phi_num);
		dis.auv_uw=(int *)malloc(sizeof(int)*adiff_num);
		dis.auv_vw=(int *)malloc(sizeof(int)*adiff_num);
		dis.auw_vw=(int *)malloc(sizeof(int)*adiff_num);
		for(int i=0;i<energy_num;i++)
		{	
			dis.energy1[i]=0;
			dis.energy2[i]=0;
		}
		for(int i=0;i<phi_num;i++) 
		{
			dis.angle[i]=0;
			dis.auv[i]=0;dis.auw[i]=0;dis.avw[i]=0;
		}
		for(int i=0;i<adiff_num;i++) 
		{
			dis.auv_uw[i]=0;
			dis.auv_vw[i]=0;
			dis.auw_vw[i]=0;
		}
		while(!feof(fp_in))
		{
			initialize_event(e);
			initialize_particle(p1);
			readevent(e,fp_in);
			if(e.nhit[u1ch]>0&&e.nhit[u2ch]>0)
			{
				e.uTsum=(e.u1[0]+e.u2[0]);
				if(e.uTsum>=param.uTsum_Lower&&e.uTsum<=param.uTsum_Up)
				{
					p1.u=(e.u1[0]-e.u2[0])*0.025*param.du;
					p1.uTsum=e.uTsum*0.025;
					p1.uflag=1;
				}
			}
			if(e.nhit[v1ch]>0&&e.nhit[v2ch]>0)
			{
				e.vTsum=(e.v1[0]+e.v2[0]);
				if(e.vTsum>=param.vTsum_Lower&&e.vTsum<=param.vTsum_Up)
				{
					p1.v=(e.v1[0]-e.v2[0])*0.025*param.dv;
					p1.vTsum=e.vTsum*0.025;
					p1.vflag=1;
				}
			}
			if(e.nhit[w1ch]>0&&e.nhit[w2ch]>0)
			{
				e.wTsum=(e.w1[0]+e.w2[0]);			
				if(e.wTsum>=param.wTsum_Lower&&e.wTsum<=param.wTsum_Up)
				{
					p1.w=(e.w1[0]-e.w2[0])*0.025*param.dw+param.ow;
					p1.wTsum=e.wTsum*0.025;
					p1.wflag=1;					
				}
			}
			getposition(p1,param);
			getenergy(p1,param);
			
			if(p1.angle>=param.Phi1_min&&p1.angle<=param.Phi1_max)
				fill_double(p1.energy,dis.energy1,energy_min,energy_max,energy_step);
			if(p1.angle>=param.Phi2_min&&p1.angle<=param.Phi2_max)
				fill_double(p1.energy,dis.energy2,energy_min,energy_max,energy_step);
			if(p1.energy>=param.E12_min&&p1.energy<=param.E12max)
			{
				fill_double(p1.angle,dis.angle,0,360,angle_step);
				fill_double(p1.auv,dis.auv,0,360,angle_step);
				fill_double(p1.auw,dis.auw,0,360,angle_step);
				fill_double(p1.avw,dis.avw,0,360,angle_step);
				if(p1.auv>0&&p1.auw>0&&p1.avw>0)
				{
					fill_double(p1.auv-p1.auw,dis.auv_uw,adiff_min,adiff_max,adiff_step);
					fill_double(p1.auv-p1.avw,dis.auv_vw,adiff_min,adiff_max,adiff_step);
					fill_double(p1.auw-p1.avw,dis.auw_vw,adiff_min,adiff_max,adiff_step);
				}
			}

		}
		add_dis(dis_sum.energy1,dis.energy1,energy_num);
		add_dis(dis_sum.energy2,dis.energy2,energy_num);
		add_dis(dis_sum.angle,dis.angle,phi_num);
		add_dis(dis_sum.auv,dis.auv,phi_num);
		add_dis(dis_sum.auw,dis.auw,phi_num);
		add_dis(dis_sum.avw,dis.avw,phi_num);
		add_dis(dis_sum.auv_uw,dis.auv_uw,adiff_num);
		add_dis(dis_sum.auv_vw,dis.auv_vw,adiff_num);
		add_dis(dis_sum.auw_vw,dis.auw_vw,adiff_num);

		fprintf(fp_energy,"energy count1 count2\n");
		for(int j=0;j<energy_num;j++)
		{
			double energy=energy_min+(j+0.5)*energy_step;
			fprintf(fp_energy,"%.3f %d %d\n",energy,dis.energy1[j],dis.energy2[j]);
		}
		fprintf(fp_angle,"a_cord angle auv auw avw\n");
		for(int j=0;j<phi_num;j++)
		{
			double phi=(j+0.5)*angle_step;
			fprintf(fp_angle,"%.1f %d %d %d %d\n",phi,dis.angle[j],dis.auv[j],dis.auw[j],dis.avw[j]);
		}
		fprintf(fp_adiff,"a_cord auv-auw auv-avw auw-avw\n");
		for(int j=0;j<adiff_num;j++)
		{
			double adiff=adiff_min+(j+0.5)*adiff_step;
			fprintf(fp_adiff,"%.2f %d %d %d\n",adiff,dis.auv_uw[j],dis.auv_vw[j],dis.auw_vw[j]);
		}
		fclose(fp_in);
		fclose(fp_energy);
		fclose(fp_angle);
		fclose(fp_adiff);
		cout<<param.file[num_file]<<" is over"<<endl;
	}
	//将分布的总和输出
	if(param.num_files>1)
	{
		FILE *fp_energy,*fp_angle,*fp_adiff;
		char out_energy[MAX_LINE],out_angle[MAX_LINE],out_adiff[MAX_LINE];	
		sprintf(out_energy,"%s%s_sum_energy.txt",param.workspace,param.name[0]);
		sprintf(out_angle,"%s%s_sum_angle.txt",param.workspace,param.name[0]);
		sprintf(out_adiff,"%s%s_sum_angle_uvwdiff.txt",param.workspace,param.name[0]);
		fp_energy=fopen(out_energy,"w+");
		fp_angle=fopen(out_angle,"w+");
		fp_adiff=fopen(out_adiff,"w+");

		fprintf(fp_energy,"energy count1 count2\n");
		for(int j=0;j<energy_num;j++)
		{
			double energy=energy_min+(j+0.5)*energy_step;
			fprintf(fp_energy,"%.3f %d %d\n",energy,dis_sum.energy1[j],dis_sum.energy2[j]);
		}
		fprintf(fp_angle,"a_cord angle auv auw avw\n");
		for(int j=0;j<phi_num;j++)
		{
			double phi=(j+0.5)*angle_step;
			fprintf(fp_angle,"%.2f %d %d %d %d\n",phi,dis_sum.angle[j],dis_sum.auv[j],dis_sum.auw[j],dis_sum.avw[j]);
		}

		fprintf(fp_adiff,"a_cord auv-auw auv-avw auw-avw\n");
		for(int j=0;j<adiff_num;j++)
		{
			double adiff=adiff_min+(j+0.5)*adiff_step;
			fprintf(fp_adiff,"%.2f %d %d %d\n",adiff,dis_sum.auv_uw[j],dis_sum.auv_vw[j],dis_sum.auw_vw[j]);
		}
		fclose(fp_energy);
		fclose(fp_angle);
		fclose(fp_adiff);

	}//End of "if(param.num_files>1)"
	return 1;
}
