/**************************************
版本2.0更新
1.重写分布，引入结构体distribution_elastic;
2.分角度拟合半径中心，二次多项式拟合
**************************************/

#include<math.h>
using namespace std;

typedef struct
{
	int **xyuv,**rphiuv,**xyuw,**rphiuw,**xyvw,**rphivw,**xy,**rphi;
	int *ruv,*ruw,*rvw,*r,*phiuv,*phiuw,*phivw,*phi;
	int *xuv_uw,*xuv_vw,*xuw_vw,*yuv_uw,*yuv_vw,*yuw_vw; //xuv-xuw,...
}distribution_elastic;

void initialize_dis_elastic(distribution_elastic &dis,int position_num,int radius_num,int phi_num)
{
	initialize_int2(dis.xyuv,position_num,position_num);
	initialize_int2(dis.xyuw,position_num,position_num);
	initialize_int2(dis.xyvw,position_num,position_num);
	initialize_int2(dis.xy,position_num,position_num);
	initialize_int2(dis.rphiuv,radius_num,phi_num);
	initialize_int2(dis.rphiuw,radius_num,phi_num);
	initialize_int2(dis.rphivw,radius_num,phi_num);
	initialize_int2(dis.rphi,radius_num,phi_num);
	initialize_int1(dis.ruv,radius_num);	initialize_int1(dis.phiuv,phi_num);
	initialize_int1(dis.ruw,radius_num);	initialize_int1(dis.phiuw,phi_num);	
	initialize_int1(dis.rvw,radius_num);	initialize_int1(dis.phivw,phi_num);
	initialize_int1(dis.r,radius_num);	initialize_int1(dis.phi,phi_num);
	initialize_int1(dis.xuv_uw,position_num);
	initialize_int1(dis.xuv_vw,position_num);
	initialize_int1(dis.xuw_vw,position_num);
	initialize_int1(dis.yuv_uw,position_num);
	initialize_int1(dis.yuv_vw,position_num);
	initialize_int1(dis.yuw_vw,position_num);
}
void free_dis_elastic(distribution_elastic &dis,int position_num,int radius_num)
{
	free(dis.r);free(dis.phi);
	free(dis.ruv);free(dis.phiuv);
	free(dis.ruw);free(dis.phiuw);
	free(dis.rvw);free(dis.phivw);
	free(dis.xuv_uw);free(dis.xuv_vw);free(dis.xuw_vw);
	free(dis.yuv_uw);free(dis.yuv_vw);free(dis.yuw_vw);
	int i;
	for(i=0;i<position_num;i++)
	{	free(dis.xyuv[i]);free(dis.xyuw[i]);free(dis.xyvw[i]);free(dis.xy[i]);	}
	for(i=0;i<radius_num;i++)
	{	free(dis.rphiuv[i]);free(dis.rphiuw[i]);free(dis.rphivw[i]);free(dis.rphi[i]);	}
}

int fit_radius(int **rphi,double *r0,double *w0,double r_min,double r_max,double r_step,double a_min,double a_max,double a_step)
{
	int i,j;
	int num_angle=(int)((a_max-a_min)/a_step+0.99);
	int num_radius=(int)((r_max-r_min)/r_step+1);
	int **dis;
	dis=(int **)malloc(sizeof(int *)*num_angle);
	for(i=0;i<num_angle;i++)
		dis[i]=(int *)malloc(sizeof(int)*num_radius);
	for(i=0;i<num_angle;i++)
		for(j=0;j<num_radius;j++)
			dis[i][j]=rphi[j][i];
	for(i=0;i<num_angle;i++)
	{
		double area=0,x0=0,width=0;
		gauss_fit(area,x0,width,dis[i],r_min,r_max,r_step);
		r0[i]=x0;
		w0[i]=width;
	}
	for(i=0;i<num_angle;i++)
		free(dis[i]);
	free(dis);
	return 1;
}

int elastic(parameter param)
{
	int time_num=(int)((param.time_max-param.time_min)/param.time_step+1);
	int position_num=(int)((param.position_max-param.position_min)/param.position_step+1);
	int radius_num=(int)((param.radius_max-param.radius_min)/param.radius_step+1);
	int phi_num=(int)((param.phi_max-param.phi_min)/param.phi_step+0.99);

	distribution_elastic dis_sum;
	initialize_dis_elastic(dis_sum,position_num,radius_num,phi_num);

	event e; 
	particle p1;

	double *r0,*x0,*y0,x0_avg=0,y0_avg=0;//各个圆弧的半径和圆心
	r0=(double *)malloc(sizeof(double)*param.num_files);
	x0=(double *)malloc(sizeof(double)*param.num_files);
	y0=(double *)malloc(sizeof(double)*param.num_files);

	FILE *fp_summary;
	char summary[MAX_LINE];
	sprintf(summary,"%s%s_sum_elastic.txt",param.workspace,param.name[0]);
	fp_summary=fopen(summary,"w+");
	fprintf(fp_summary,"FILE total_count elastic_count x_center y_center radius area rpeak rwidth\n");

	double *area,*rpeak,*rwidth; //每段圆弧所有角度半径分布高斯拟合值
	area=(double *)malloc(sizeof(double)*param.num_files);
	rpeak=(double *)malloc(sizeof(double)*param.num_files);
	rwidth=(double *)malloc(sizeof(double)*param.num_files);

	double **rauv,**wauv,**rauw,**wauw,**ravw,**wavw;//分角度刻度的中心值和半宽
	rauv=(double **)malloc(sizeof(double *)*phi_num);
	wauv=(double **)malloc(sizeof(double *)*phi_num);
	rauw=(double **)malloc(sizeof(double *)*phi_num);
	wauw=(double **)malloc(sizeof(double *)*phi_num);
	ravw=(double **)malloc(sizeof(double *)*phi_num);
	wavw=(double **)malloc(sizeof(double *)*phi_num);
	for(int i=0;i<phi_num;i++)
	{
		rauv[i]=(double *)malloc(sizeof(double)*param.num_files);
		wauv[i]=(double *)malloc(sizeof(double)*param.num_files);
		rauw[i]=(double *)malloc(sizeof(double)*param.num_files);
		wauw[i]=(double *)malloc(sizeof(double)*param.num_files);
		ravw[i]=(double *)malloc(sizeof(double)*param.num_files);
		wavw[i]=(double *)malloc(sizeof(double)*param.num_files);
	}

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

		FILE *fp_out,*fp_outxyDis,*fp_rphiDis,*fp_radius,*fp_phi,*fp_xydiff;
		char output_name[MAX_LINE],outxyDis_name[MAX_LINE],rphiDis[MAX_LINE],radius[MAX_LINE],phi[MAX_LINE],xtdiff[MAX_LINE];	
		sprintf(output_name,"%s%s_uvwxy.txt",param.workspace,param.name[num_file]);
		sprintf(outxyDis_name,"%s%s_xy.txt",param.workspace,param.name[num_file]);
		sprintf(rphiDis,"%s%s_rphi.txt",param.workspace,param.name[num_file]);
		sprintf(radius,"%s%s_radius.txt",param.workspace,param.name[num_file]);
		sprintf(phi,"%s%s_phi.txt",param.workspace,param.name[num_file]);
		sprintf(xtdiff,"%s%s_xydifference.txt",param.workspace,param.name[num_file]);		

		fp_out=fopen(output_name,"w+");
		fp_outxyDis=fopen(outxyDis_name,"w+");
		fp_rphiDis=fopen(rphiDis,"w+");
		fp_radius=fopen(radius,"w+");
		fp_phi=fopen(phi,"w+");
		fp_xydiff=fopen(xtdiff,"w+");
		
		distribution_elastic dis;
		initialize_dis_elastic(dis,position_num,radius_num,phi_num);
		int elastic_count=0,total_count=0;

		fprintf(fp_out,"count u v w xuv yuv xuw yuw xvw yvw x y uTsum vTsum wTsum\n");
		while(!feof(fp_in))
		{
			initialize_event(e);
			initialize_particle(p1);
			readevent(e,fp_in);
			total_count++;
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
			if(getposition(p1,param))
			{
				elastic_count++;
			}

			fill_double(p1.ruv,dis.ruv,param.radius_min,param.radius_max,param.radius_step);
			fill_double(p1.ruw,dis.ruw,param.radius_min,param.radius_max,param.radius_step);
			fill_double(p1.rvw,dis.rvw,param.radius_min,param.radius_max,param.radius_step);
			fill_double(p1.r,dis.r,param.radius_min,param.radius_max,param.radius_step);
			fill_double(p1.phiuv,dis.phiuv,param.phi_min,param.phi_max,param.phi_step);
			fill_double(p1.phiuw,dis.phiuw,param.phi_min,param.phi_max,param.phi_step);
			fill_double(p1.phivw,dis.phivw,param.phi_min,param.phi_max,param.phi_step);
			fill_double(p1.phi,dis.phi,param.phi_min,param.phi_max,param.phi_step);
			fill2D_double(p1.xuv,p1.yuv,dis.xyuv,param.position_min,param.position_max,param.position_step,param.position_min,param.position_max,param.position_step);
			fill2D_double(p1.xuw,p1.yuw,dis.xyuw,param.position_min,param.position_max,param.position_step,param.position_min,param.position_max,param.position_step);
			fill2D_double(p1.xvw,p1.yvw,dis.xyvw,param.position_min,param.position_max,param.position_step,param.position_min,param.position_max,param.position_step);
			fill2D_double(p1.x,p1.y,dis.xy,param.position_min,param.position_max,param.position_step,param.position_min,param.position_max,param.position_step);
			fill2D_double(p1.ruv,p1.phiuv,dis.rphiuv,param.radius_min,param.radius_max,param.radius_step,param.phi_min,param.phi_max,param.phi_step);
			fill2D_double(p1.ruw,p1.phiuw,dis.rphiuw,param.radius_min,param.radius_max,param.radius_step,param.phi_min,param.phi_max,param.phi_step);
			fill2D_double(p1.rvw,p1.phivw,dis.rphivw,param.radius_min,param.radius_max,param.radius_step,param.phi_min,param.phi_max,param.phi_step);
			fill2D_double(p1.r,p1.phi,dis.rphi,param.radius_min,param.radius_max,param.radius_step,param.phi_min,param.phi_max,param.phi_step);	

			if(p1.uflag&&p1.vflag&&p1.wflag)
			{
				fill_double(p1.xuv-p1.xuw,dis.xuv_uw,param.position_min,param.position_max,param.position_step);
				fill_double(p1.xuv-p1.xvw,dis.xuv_vw,param.position_min,param.position_max,param.position_step);
				fill_double(p1.xuw-p1.xvw,dis.xuw_vw,param.position_min,param.position_max,param.position_step);
				fill_double(p1.yuv-p1.yuw,dis.yuv_uw,param.position_min,param.position_max,param.position_step);
				fill_double(p1.yuv-p1.yvw,dis.yuv_vw,param.position_min,param.position_max,param.position_step);
				fill_double(p1.yuw-p1.yvw,dis.yuw_vw,param.position_min,param.position_max,param.position_step);
			}
			if(e.count<MAX_origin_num)
			fprintf(fp_out,"%d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",e.count,
				p1.u,p1.v,p1.w,p1.xuv,p1.yuv,p1.xuw,p1.yuw,p1.xvw,p1.yvw,p1.x,p1.y,p1.uTsum,p1.vTsum,p1.wTsum);
		}

		add_dis(dis_sum.ruv,dis.ruv,radius_num);
		add_dis(dis_sum.ruw,dis.ruw,radius_num);
		add_dis(dis_sum.rvw,dis.rvw,radius_num);
		add_dis(dis_sum.r,dis.r,radius_num);
		add_dis(dis_sum.phiuv,dis.phiuv,phi_num);
		add_dis(dis_sum.phiuw,dis.phiuw,phi_num);
		add_dis(dis_sum.phivw,dis.phivw,phi_num);
		add_dis(dis_sum.phi,dis.phi,phi_num);
		add_dis(dis_sum.xuv_uw,dis.xuv_uw,position_num);
		add_dis(dis_sum.xuv_vw,dis.xuv_vw,position_num);
		add_dis(dis_sum.xuw_vw,dis.xuw_vw,position_num);
		add_dis(dis_sum.yuv_uw,dis.yuv_uw,position_num);
		add_dis(dis_sum.yuv_vw,dis.yuv_vw,position_num);
		add_dis(dis_sum.yuw_vw,dis.yuw_vw,position_num);
		add_dis2D(dis_sum.xyuv,dis.xyuv,position_num,position_num);
		add_dis2D(dis_sum.xyuw,dis.xyuw,position_num,position_num);
		add_dis2D(dis_sum.xyvw,dis.xyvw,position_num,position_num);
		add_dis2D(dis_sum.xy,dis.xy,position_num,position_num);
		add_dis2D(dis_sum.rphiuv,dis.rphiuv,radius_num,phi_num);
		add_dis2D(dis_sum.rphiuw,dis.rphiuw,radius_num,phi_num);
		add_dis2D(dis_sum.rphivw,dis.rphivw,radius_num,phi_num);		
		add_dis2D(dis_sum.rphi,dis.rphi,radius_num,phi_num);

		fprintf(fp_outxyDis,"x y xyuv xyuw xyvw xy\n");
		for(int i=0;i<position_num;i++)
		{
			for(int j=0;j<position_num;j++)
			{
				double x=param.position_min+(i+0.5)*param.position_step;
				double y=param.position_min+(j+0.5)*param.position_step;
				fprintf(fp_outxyDis,"%.2f %.2f %d %d %d %d\n",x,y,dis.xyuv[i][j],dis.xyuw[i][j],dis.xyvw[i][j],dis.xy[i][j]);
			}
		}
		fprintf(fp_rphiDis,"r phi rphiuv rphiuw rphivw rphi\n");
		for(int i=0;i<radius_num;i++)
		{
			for(int j=0;j<phi_num;j++)
			{
				double r=param.radius_min+(i+0.5)*param.radius_step;
				double phi=param.phi_min+(j+0.5)*param.phi_step;
				fprintf(fp_rphiDis,"%.2f %.1f %d %d %d %d\n",r,phi,dis.rphiuv[i][j],dis.rphiuw[i][j],dis.rphivw[i][j],dis.rphi[i][j]);
			}
		}

		fprintf(fp_radius,"cord ruv ruw rvw r\n");
		for(int i=0;i<radius_num;i++)
		{
			double r=param.radius_min+(i+0.5)*param.radius_step;
			fprintf(fp_radius,"%.2f %d %d %d %d\n",r,dis.ruv[i],dis.ruw[i],dis.rvw[i],dis.r[i]);
		}

		fprintf(fp_phi,"cord phiuv phiuw phivw phi\n");
		for(int j=0;j<phi_num;j++)
		{
			double phi=param.phi_min+(j+0.5)*param.phi_step;
			fprintf(fp_phi,"%.1f %d %d %d %d\n",phi,dis.phiuv[j],dis.phiuw[j],dis.phivw[j],dis.phi[j]);
		}

		fprintf(fp_xydiff,"position xuv-xuw xuv-xvw xuw-xvw yuv-yuw yuv-yvw yuw-yvw\n");
		for(int i=0;i<position_num;i++)
		{
			double x=param.position_min+(i+0.5)*param.position_step;
			fprintf(fp_xydiff,"%.2f %d %d %d %d %d %d\n",x,dis.xuv_uw[i],dis.xuv_vw[i],dis.xuw_vw[i],dis.yuv_uw[i],dis.yuv_vw[i],dis.yuw_vw[i]);
		}

		fclose(fp_in);
		fclose(fp_out);
		fclose(fp_outxyDis);
		fclose(fp_rphiDis);
		fclose(fp_radius);
		fclose(fp_phi);
		fclose(fp_xydiff);
		
		getcenter(x0[num_file],y0[num_file],r0[num_file],dis.xy,param);
		x0_avg+=x0[num_file];
		y0_avg+=y0[num_file];
		fprintf(fp_summary,"%s %d %d %.3f %.3f %.3f ",param.name[num_file],total_count,elastic_count,x0[num_file],y0[num_file],r0[num_file]);
		if(gauss_fit(area[num_file],rpeak[num_file],rwidth[num_file],dis.r,param.radius_min,param.radius_max,param.radius_step))
			fprintf(fp_summary,"%.3f %.3f %.3f\n",area[num_file],rpeak[num_file],rwidth[num_file]);
		else
			fprintf(fp_summary,"\n");

		double *rt,*wt;
		rt=(double *)malloc(sizeof(double)*phi_num);
		wt=(double *)malloc(sizeof(double)*phi_num);
		fit_radius(dis.rphiuv,rt,wt,param.radius_min,param.radius_max,param.radius_step,param.phi_min,param.phi_max,param.phi_step);
		for(int i=0;i<phi_num;i++)
		{
			rauv[i][num_file]=rt[i];
			wauv[i][num_file]=wt[i];
		}
		fit_radius(dis.rphiuw,rt,wt,param.radius_min,param.radius_max,param.radius_step,param.phi_min,param.phi_max,param.phi_step);
		for(int i=0;i<phi_num;i++)
		{
			rauw[i][num_file]=rt[i];
			wauw[i][num_file]=wt[i];
		}
		fit_radius(dis.rphivw,rt,wt,param.radius_min,param.radius_max,param.radius_step,param.phi_min,param.phi_max,param.phi_step);
		for(int i=0;i<phi_num;i++)
		{
			ravw[i][num_file]=rt[i];
			wavw[i][num_file]=wt[i];
		}
		free(rt);free(wt);
		cout<<param.file[num_file]<<" is over"<<endl;
	}

	//将分布的总和输出
	if(param.num_files>1)
	{
		x0_avg/=param.num_files;
		y0_avg/=param.num_files;
		fprintf(fp_summary,"\nAveraged center: %.3f %.3f\n",x0_avg,y0_avg);

		FILE *fp_outxyDis_sum,*fp_rphiDis_sum,*fp_radius_sum,*fp_phi_sum,*fp_xydiff_sum;
		char outxyDis_sum[MAX_LINE],rphiDis_sum[MAX_LINE],radius_sum[MAX_LINE],phi_sum[MAX_LINE],xydiff_sum[MAX_LINE];
	
		sprintf(outxyDis_sum,"%s%s_sum_xy.txt",param.workspace,param.name[0]);
		sprintf(rphiDis_sum,"%s%s_sum_rphi.txt",param.workspace,param.name[0]);
		sprintf(radius_sum,"%s%s_sum_radius.txt",param.workspace,param.name[0]);
		sprintf(phi_sum,"%s%s_sum_phi.txt",param.workspace,param.name[0]);
		sprintf(xydiff_sum,"%s%s_sum_xydifference.txt",param.workspace,param.name[0]);

		fp_outxyDis_sum=fopen(outxyDis_sum,"w+");
		fp_rphiDis_sum=fopen(rphiDis_sum,"w+");
		fp_radius_sum=fopen(radius_sum,"w+");
		fp_phi_sum=fopen(phi_sum,"w+");
		fp_xydiff_sum=fopen(xydiff_sum,"w+");

		fprintf(fp_outxyDis_sum,"x y xyuv xyuw xyvw xy\n");
		for(int i=0;i<position_num;i++)
		{
			for(int j=0;j<position_num;j++)
			{
				double x=param.position_min+(i+0.5)*param.position_step;
				double y=param.position_min+(j+0.5)*param.position_step;
				fprintf(fp_outxyDis_sum,"%.2f %.2f %d %d %d %d\n",x,y,dis_sum.xyuv[i][j],dis_sum.xyuw[i][j],dis_sum.xyvw[i][j],dis_sum.xy[i][j]);
			}
		}
		fprintf(fp_rphiDis_sum,"r phi rphiuv rphiuw rphivw rphi\n");
		for(int i=0;i<radius_num;i++)
		{
			for(int j=0;j<phi_num;j++)
			{
				double r=param.radius_min+(i+0.5)*param.radius_step;
				double phi=param.phi_min+(j+0.5)*param.phi_step;
				fprintf(fp_rphiDis_sum,"%.2f %.1f %d %d %d %d\n",r,phi,dis_sum.rphiuv[i][j],dis_sum.rphiuw[i][j],dis_sum.rphivw[i][j],dis_sum.rphi[i][j]);
			}
		}

		fprintf(fp_radius_sum,"cord ruv ruw rvw r\n");
		for(int i=0;i<radius_num;i++)
		{
			double r=param.radius_min+(i+0.5)*param.radius_step;
			fprintf(fp_radius_sum,"%.2f %d %d %d %d\n",r,dis_sum.ruv[i],dis_sum.ruw[i],dis_sum.rvw[i],dis_sum.r[i]);
		}

		fprintf(fp_phi_sum,"cord phiuv phiuw phivw phi\n");
		for(int j=0;j<phi_num;j++)
		{
			double phi=param.phi_min+(j+0.5)*param.phi_step;
			fprintf(fp_phi_sum,"%.1f %d %d %d %d\n",phi,dis_sum.phiuv[j],dis_sum.phiuw[j],dis_sum.phivw[j],dis_sum.phi[j]);
		}

		fprintf(fp_xydiff_sum,"position xuv-xuw xuv-xvw xuw-xvw yuv-yuw yuv-yvw yuw-yvw\n");
		for(int i=0;i<position_num;i++)
		{
			double x=param.position_min+(i+0.5)*param.position_step;
			fprintf(fp_xydiff_sum,"%.2f %d %d %d %d %d %d\n",x,dis_sum.xuv_uw[i],dis_sum.xuv_vw[i],dis_sum.xuw_vw[i],dis_sum.yuv_uw[i],dis_sum.yuv_vw[i],dis_sum.yuw_vw[i]);
		}

		fclose(fp_outxyDis_sum);
		fclose(fp_rphiDis_sum);
		fclose(fp_radius_sum);
		fclose(fp_phi_sum);
		fclose(fp_xydiff_sum);
	}//End of "if(param.num_files>1)"
	fclose(fp_summary);

	//输出拟合的分角度半径值和半宽
	FILE *fp_calib;
	char calib[MAX_LINE];
	sprintf(calib,"%s%s_sum_calib.txt",param.workspace,param.name[0]);
	fp_calib=fopen(calib,"w+");
	fprintf(fp_calib,"Angle Energy ruv wuv ruw wuv rvw uvw\n");
	for(int i=0;i<phi_num;i++)
		for(int j=0;j<param.num_files;j++)
		{
			fprintf(fp_calib,"%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",
				param.phi_min+param.phi_step*(i+0.5),param.elastic_energy[j],
				rauv[i][j],wauv[i][j],rauw[i][j],wauw[i][j],ravw[i][j],wavw[i][j]);
		}
	fclose(fp_calib);

	if(param.num_files>2)
	{
		double *calib_r,*calib_e;//多项式拟合的r,E值
		double calib_c[3];
		calib_r=(double *)malloc(sizeof(double)*param.num_files);
		calib_e=(double *)malloc(sizeof(double)*param.num_files);
		
		if(param.cali_angle_num<phi_num)
		{
			param.cuv0=(double *)realloc(param.cuv0,sizeof(double)*phi_num); 
			param.cuv1=(double *)realloc(param.cuv1,sizeof(double)*phi_num); 
			param.cuv2=(double *)realloc(param.cuv2,sizeof(double)*phi_num); 
			param.cuw0=(double *)realloc(param.cuw0,sizeof(double)*phi_num); 
			param.cuw1=(double *)realloc(param.cuw1,sizeof(double)*phi_num); 
			param.cuw2=(double *)realloc(param.cuw2,sizeof(double)*phi_num); 
			param.cvw0=(double *)realloc(param.cvw0,sizeof(double)*phi_num); 
			param.cvw1=(double *)realloc(param.cvw1,sizeof(double)*phi_num); 
			param.cvw2=(double *)realloc(param.cvw2,sizeof(double)*phi_num);
		}
		param.cali_angle_num=phi_num;

		for(int i=0;i<phi_num;i++)
		{
			int k=0;
			for(int j=0;j<param.num_files;j++)
			{
				if(rauv[i][j]>1.0)
				{
					calib_r[k]=rauv[i][j];
					calib_e[k++]=param.elastic_energy[j];
				}
			}
			polynomial_fit(calib_r,calib_e,2,calib_c,k);
			param.cuv0[i]=calib_c[0];
			param.cuv1[i]=calib_c[1];
			param.cuv2[i]=calib_c[2];
		}

		for(int i=0;i<phi_num;i++)
		{
			int k=0;
			for(int j=0;j<param.num_files;j++)
			{
				if(rauw[i][j]>1.0)
				{
					calib_r[k]=rauw[i][j];
					calib_e[k++]=param.elastic_energy[j];
				}
			}
			polynomial_fit(calib_r,calib_e,2,calib_c,k);
			param.cuw0[i]=calib_c[0];
			param.cuw1[i]=calib_c[1];
			param.cuw2[i]=calib_c[2];
		}

		for(int i=0;i<phi_num;i++)
		{
			int k=0;
			for(int j=0;j<param.num_files;j++)
			{
				if(ravw[i][j]>1.0)
				{
					calib_r[k]=ravw[i][j];
					calib_e[k++]=param.elastic_energy[j];
				}
			}
			polynomial_fit(calib_r,calib_e,2,calib_c,k);
			param.cvw0[i]=calib_c[0];
			param.cvw1[i]=calib_c[1];
			param.cvw2[i]=calib_c[2];
		}
		print_param(param);
	}
	return 1;
}


int directxy(parameter param)
{
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

		FILE *fp_uv,*fp_uw,*fp_vw;
		char out_uv[MAX_LINE],out_uw[MAX_LINE],out_vw[MAX_LINE];	
		sprintf(out_uv,"%s%s_xyuv.txt",param.workspace,param.name[num_file]);
		sprintf(out_uw,"%s%s_xyuw.txt",param.workspace,param.name[num_file]);
		sprintf(out_vw,"%s%s_xyvw.txt",param.workspace,param.name[num_file]);

		fp_uv=fopen(out_uv,"w+");
		fp_uw=fopen(out_uw,"w+");
		fp_vw=fopen(out_vw,"w+");

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
			if(p1.ruv<param.radius_max)
				fprintf(fp_uv,"%.3f %.3f\n",p1.xuv,p1.yuv);
			if(p1.ruw<param.radius_max)
				fprintf(fp_uw,"%.3f %.3f\n",p1.xuw,p1.yuw);
			if(p1.ruv<param.radius_max)
				fprintf(fp_vw,"%.3f %.3f\n",p1.xvw,p1.yvw);
		}

		fclose(fp_in);
		fclose(fp_uv);
		fclose(fp_uw);
		fclose(fp_vw);
		cout<<param.file[num_file]<<" is over"<<endl;
	}
	//将分布的总和输出
	return 1;
}