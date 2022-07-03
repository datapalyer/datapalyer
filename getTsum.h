/********************************
getTsumx.x.cpp：getTsum函数int getTsum(parameter param)
版本1.1
输出文件：
1.每个数据文件对应的输出
(1) "*_pattern.txt"  pattern分布，用来判断信号丢失；
(2)"*_Tsum.txt"  Tsum分布图；
(3)"*_Tsum2D.txt" Tsum二维分布图；
(4)"*_origin.txt" 输出MAX_origin_num个转化的原始数据。
2.多个个数据文件的总结
(1) "*_sum_pattern.txt"  pattern分布，用来判断信号丢失；
(2)"*_sum_Tsum.txt"  Tsum分布图；
(3)"*_sum_Tsum2D.txt" Tsum二维分布图；
(4)"*_sum_TsumSummary.txt" getTsum总结，包括u1,u2，u1&u2,u1&u2&uTsum的计数，Tsum分布Gassian拟合的结果。

版本2.0更新：
1.完全更改了分布的写法，定义了Tsum分布结构体distribution_Tsum，提高程序可读性；
2.去掉了在Tsum范围内的eff分布，这个在弹性和符合数据处理中再体现；
3.修改相关bug,如开辟的内存空间及时释放等。
*************************************/

using namespace std;
int readevent(event &e,FILE *fp_in);

//Tsum谱结构体
typedef struct
{
		int *uTsum,*vTsum,*wTsum; //Tsum分布
		int *u11,*u21,*v11,*v21,*w11,*w21;  //第一击的分布
		int *u12,*u22,*v12,*v22,*w12,*w22,*MCP2;  //第二击的分布
		int *u1diff,*u2diff,*v1diff,*v2diff,*w1diff,*w2diff;  //第一击和第二击时间差的分布
		int **uTsum2D,**vTsum2D,**wTsum2D; //Tsum二维分布
		int **pattern;	//pattern的分布
}distribution_Tsum;

void initialize_dis_Tsum(distribution_Tsum &dis,int time_num,int position_num)
{
	initialize_int1(dis.uTsum,time_num);initialize_int1(dis.vTsum,time_num);initialize_int1(dis.wTsum,time_num);
	initialize_int1(dis.u11,time_num);initialize_int1(dis.u21,time_num);
	initialize_int1(dis.v11,time_num);initialize_int1(dis.v21,time_num);
	initialize_int1(dis.w11,time_num);initialize_int1(dis.w21,time_num);
	initialize_int1(dis.u12,time_num);initialize_int1(dis.u22,time_num);
	initialize_int1(dis.v12,time_num);initialize_int1(dis.v22,time_num);
	initialize_int1(dis.w12,time_num);initialize_int1(dis.w22,time_num);
	initialize_int1(dis.MCP2,time_num);
	initialize_int1(dis.u1diff,time_num);initialize_int1(dis.u2diff,time_num);
	initialize_int1(dis.v1diff,time_num);initialize_int1(dis.v2diff,time_num);
	initialize_int1(dis.w1diff,time_num);initialize_int1(dis.w2diff,time_num);
	initialize_int2(dis.uTsum2D,time_num,position_num);
	initialize_int2(dis.vTsum2D,time_num,position_num);
	initialize_int2(dis.wTsum2D,time_num,position_num);
	initialize_int2(dis.pattern,64,max_hit);
}

void free_dis_Tsum(distribution_Tsum &dis,int time_num)
{
	free(dis.uTsum);free(dis.vTsum);free(dis.wTsum);
	free(dis.u11);free(dis.u21);free(dis.v11);free(dis.v21);free(dis.w11);free(dis.w21);
	free(dis.u12);free(dis.u22);free(dis.v12);free(dis.v22);free(dis.w12);free(dis.w22);free(dis.MCP2);
	free(dis.u1diff);free(dis.u2diff);free(dis.v1diff);free(dis.v2diff);free(dis.w1diff);free(dis.w2diff);
	for(int i=0;i<time_num;i++)
	{
		free(dis.uTsum2D[i]);free(dis.vTsum2D[i]);free(dis.wTsum2D[i]);
	}
	free(dis.uTsum2D);free(dis.vTsum2D);free(dis.wTsum2D);
	for(int i=0;i<64;i++)
		free(dis.pattern[i]);
	free(dis.pattern);
}

int getTsum(parameter param)
{
	int time_num=(int)((param.time_max-param.time_min)/param.time_step+1);
	int position_num=(int)((param.position_max-param.position_min)/param.position_step+1);

	distribution_Tsum dis_sum;
	initialize_dis_Tsum(dis_sum,time_num,position_num);

	int *totalcount,*u1count,*u2count,*v1count,*v2count,*w1count,*w2count;
	int *u12count,*u12effcount,*v12count,*v12effcount,*w12count,*w12effcount;
	totalcount=(int *)malloc(sizeof(int)*(param.num_files+1));
	u1count=(int *)malloc(sizeof(int)*(param.num_files+1));
	u2count=(int *)malloc(sizeof(int)*(param.num_files+1));
	v1count=(int *)malloc(sizeof(int)*(param.num_files+1));
	v2count=(int *)malloc(sizeof(int)*(param.num_files+1));
	w1count=(int *)malloc(sizeof(int)*(param.num_files+1));
	w2count=(int *)malloc(sizeof(int)*(param.num_files+1));
	u12count=(int *)malloc(sizeof(int)*(param.num_files+1));
	u12effcount=(int *)malloc(sizeof(int)*(param.num_files+1));
	v12count=(int *)malloc(sizeof(int)*(param.num_files+1));
	v12effcount=(int *)malloc(sizeof(int)*(param.num_files+1));
	w12count=(int *)malloc(sizeof(int)*(param.num_files+1));
	w12effcount=(int *)malloc(sizeof(int)*(param.num_files+1));
	for(int icount=0;icount<(param.num_files+1);icount++)
	{
		totalcount[icount]=0;u1count[icount]=0;u2count[icount]=0;v1count[icount]=0;v2count[icount]=0;w1count[icount]=0;w2count[icount]=0;
		u12count[icount]=0;u12effcount[icount]=0;v12count[icount]=0;v12effcount[icount]=0;w12count[icount]=0;w12effcount[icount]=0;
	}
	event e; 
	int num_file;
	for(num_file=0;num_file<param.num_files;num_file++)
	{
		FILE *fp_in,*fp_out,*fp_outpattern,*fp_outTsum,*fp_outTsumDis;
		fp_in=fopen(param.file[num_file],"rb");
		if(fp_in==NULL)
		{
			cout<<"Input file error!"<<endl;
			return 0;
		}
		char output_name[MAX_LINE],outpattern_name[MAX_LINE],outTsum_name[MAX_LINE],outTsumDis_name[MAX_LINE];
		sprintf(outpattern_name,"%s%s_pattern.txt",param.workspace,param.name[num_file]);
		sprintf(outTsum_name,"%s%s_Tsum.txt",param.workspace,param.name[num_file]);
		sprintf(outTsumDis_name,"%s%s_Tsum2D.txt",param.workspace,param.name[num_file]);
		sprintf(output_name,"%s%s_origin.txt",param.workspace,param.name[num_file]);

		fp_out=fopen(output_name,"w+");
		fp_outpattern=fopen(outpattern_name,"w+");
		fp_outTsum=fopen(outTsum_name,"w+");
		fp_outTsumDis=fopen(outTsumDis_name,"w+");

		distribution_Tsum dis;
		initialize_dis_Tsum(dis,time_num,position_num);

		fprintf(fp_out,"count u1_1 u1_2 u2_1 u2_2 v1_1 v1_2 v2_1 v2_2 w1_1 w1_2 w2_1 w2_2 MCP0 MCP1 status uTsum vTsum wTsum\n");
		while(!feof(fp_in))
		{
			initialize_event(e);
			readevent(e,fp_in);//读一个事例

			totalcount[num_file]++; totalcount[param.num_files]++;
			u1count[num_file]+=fill_int(e.u1[0],dis.u11,param.time_min,param.time_max,param.time_step);
			u2count[num_file]+=fill_int(e.u2[0],dis.u21,param.time_min,param.time_max,param.time_step);
			v1count[num_file]+=fill_int(e.v1[0],dis.v11,param.time_min,param.time_max,param.time_step);
			v2count[num_file]+=fill_int(e.v2[0],dis.v21,param.time_min,param.time_max,param.time_step);
			w1count[num_file]+=fill_int(e.w1[0],dis.w11,param.time_min,param.time_max,param.time_step);
			w2count[num_file]+=fill_int(e.w2[0],dis.w21,param.time_min,param.time_max,param.time_step);

			if(e.nhit[u1ch]>1){ 
				fill_int(e.u1[1]-e.u1[0],dis.u1diff,param.time_min,param.time_max,param.time_step);
				fill_int(e.u1[1],dis.u12,param.time_min,param.time_max,param.time_step);
			}
			if(e.nhit[u2ch]>1){
				fill_int(e.u2[1]-e.u2[0],dis.u2diff,param.time_min,param.time_max,param.time_step);
				fill_int(e.u2[1],dis.u22,param.time_min,param.time_max,param.time_step);
			}
			if(e.nhit[v1ch]>1){
				fill_int(e.v1[1]-e.v1[0],dis.v1diff,param.time_min,param.time_max,param.time_step);
				fill_int(e.v1[1],dis.v12,param.time_min,param.time_max,param.time_step);
			}
			if(e.nhit[v2ch]>1){
				fill_int(e.v2[1]-e.v2[0],dis.v2diff,param.time_min,param.time_max,param.time_step);
				fill_int(e.v2[1],dis.v22,param.time_min,param.time_max,param.time_step);
			}
			if(e.nhit[w1ch]>1){
				fill_int(e.w1[1]-e.w1[0],dis.w1diff,param.time_min,param.time_max,param.time_step);
				fill_int(e.w1[1],dis.w12,param.time_min,param.time_max,param.time_step);
			}
			if(e.nhit[w2ch]>1){
				fill_int(e.w2[1]-e.w2[0],dis.w2diff,param.time_min,param.time_max,param.time_step);
				fill_int(e.w2[1],dis.w22,param.time_min,param.time_max,param.time_step);
			}

			if(e.nhit[u1ch]>0&&e.nhit[u2ch]>0)
			{
				e.uTsum=(e.u1[0]+e.u2[0]);
				double u=(e.u1[0]-e.u2[0])*0.025*param.du;
				u12count[num_file]+=fill_int(e.uTsum,dis.uTsum,param.time_min,param.time_max,param.time_step);
				fill2D_double(e.uTsum,u,dis.uTsum2D,param.time_min,param.time_max,param.time_step,param.position_min,param.position_max,param.position_step);
				if(e.uTsum>=param.uTsum_Lower&&e.uTsum<=param.uTsum_Up)
				{
					u12effcount[num_file]++;
				}
			}
			else
				e.uTsum=-1;

			if(e.nhit[v1ch]>0&&e.nhit[v2ch]>0)
			{
				e.vTsum=(e.v1[0]+e.v2[0]);
				double v=(e.v1[0]-e.v2[0])*0.025*param.dv;
				v12count[num_file]+=fill_int(e.vTsum,dis.vTsum,param.time_min,param.time_max,param.time_step);
				fill2D_double(e.vTsum,v,dis.vTsum2D,param.time_min,param.time_max,param.time_step,param.position_min,param.position_max,param.position_step);

				if(e.vTsum>=param.vTsum_Lower&&e.vTsum<=param.vTsum_Up)
				{
					v12effcount[num_file]++;
				}
			}
			else
				e.vTsum=-1;

			if(e.nhit[w1ch]>0&&e.nhit[w2ch]>0)
			{
				e.wTsum=(e.w1[0]+e.w2[0]);
				double w=(e.w1[0]-e.w2[0])*0.025*param.dw+param.ow;
				w12count[num_file]+=fill_int(e.wTsum,dis.wTsum,param.time_min,param.time_max,param.time_step);
				fill2D_double(e.wTsum,w,dis.wTsum2D,
					param.time_min,param.time_max,param.time_step,param.position_min,param.position_max,param.position_step);

				if(e.wTsum>=param.wTsum_Lower&&e.wTsum<=param.wTsum_Up)
				{
					w12effcount[num_file]++;
				}
			}
			else
				e.vTsum=-1;

			if(e.nhit[MCPch]>1)
			{
				fill_int(e.MCP[1],dis.MCP2,param.time_min,param.time_max,param.time_step);
			}
			dis.pattern[e.pattern[0]][0]++;
			dis.pattern[e.pattern[1]][1]++;
			dis.pattern[e.pattern[2]][2]++;

			if(e.count<MAX_origin_num)
				fprintf(fp_out,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",e.count,e.u1[0],e.u1[1],e.u2[0],e.u2[1],
					e.v1[0],e.v1[1],e.w1[0],e.w1[1],e.w2[0],e.w2[1],e.MCP[0],e.MCP[1],e.status,e.uTsum,e.vTsum,e.wTsum);
		}//End of while(!feof(fp_in))

		for(int i=0;i<64;i++)
		{
			fprintf(fp_outpattern,"%d ",i);
			for(int j=0;j<max_hit;j++)
			{
				fprintf(fp_outpattern,"%d ",dis.pattern[i][j]);
			}
			fprintf(fp_outpattern,"\n");
		}
		fprintf(fp_outTsum,"time_mid t uTsum vTsum wTsum u1 u2 v1 v2 w1 w2 u12 u22 v12 v22 w12 w22 MCP2 u1diff u2diff v1diff v1diff w1diff w2diff\n");
		for(int j=0;j<time_num;j++)
		{
			double time_mid=param.time_min+(j+0.5)*param.time_step;
			double t=time_mid*0.025;
			fprintf(fp_outTsum,"%.1f %.3f %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
				time_mid,t,dis.uTsum[j],dis.vTsum[j],dis.wTsum[j],
				dis.u11[j],dis.u21[j],dis.v11[j],dis.v21[j],dis.w11[j],dis.w21[j],
				dis.u12[j],dis.u22[j],dis.v12[j],dis.v22[j],dis.w12[j],dis.w22[j],dis.MCP2[j],
				dis.u1diff[j],dis.u2diff[j],dis.v1diff[j],dis.v2diff[j],dis.w1diff[j],dis.w2diff[j]);
		}
		for(int iT=0;iT<time_num;iT++)
		{
			for(int k=0;k<position_num;k++)
			{
				double time_mid=param.time_min+(iT+0.5)*param.time_step;
				double position_mid=param.position_min+(k+0.5)*param.position_step;
				fprintf(fp_outTsumDis,"%.1f %.1f %d %d %d\n",time_mid,position_mid,dis.uTsum2D[iT][k],dis.vTsum2D[iT][k],dis.wTsum2D[iT][k]);
			}
		}

		fclose(fp_in);
		fclose(fp_out);
		fclose(fp_outpattern);
		fclose(fp_outTsum);
		fclose(fp_outTsumDis);
		cout<<param.file[num_file]<<" is over"<<endl;

		u1count[param.num_files]+=u1count[num_file];
		u2count[param.num_files]+=u2count[num_file];
		v1count[param.num_files]+=v1count[num_file];
		v2count[param.num_files]+=v2count[num_file];
		w1count[param.num_files]+=w1count[num_file];
		w2count[param.num_files]+=w2count[num_file];
		u12count[param.num_files]+=u12count[num_file];
		v12count[param.num_files]+=v12count[num_file];
		w12count[param.num_files]+=w12count[num_file];
		u12effcount[param.num_files]+=u12effcount[num_file];
		v12effcount[param.num_files]+=v12effcount[num_file];
		w12effcount[param.num_files]+=w12effcount[num_file];

		add_dis(dis_sum.uTsum,dis.uTsum,time_num);add_dis(dis_sum.vTsum,dis.vTsum,time_num);add_dis(dis_sum.wTsum,dis.wTsum,time_num);
		add_dis(dis_sum.u11,dis.u11,time_num);add_dis(dis_sum.u21,dis.u21,time_num);
		add_dis(dis_sum.v11,dis.v11,time_num);add_dis(dis_sum.v21,dis.v21,time_num);
		add_dis(dis_sum.w11,dis.w11,time_num);add_dis(dis_sum.w21,dis.w21,time_num);
		add_dis(dis_sum.u12,dis.u12,time_num);add_dis(dis_sum.u22,dis.u22,time_num);
		add_dis(dis_sum.v12,dis.v12,time_num);add_dis(dis_sum.v22,dis.v22,time_num);
		add_dis(dis_sum.w12,dis.w12,time_num);add_dis(dis_sum.w22,dis.w22,time_num);add_dis(dis_sum.MCP2,dis.MCP2,time_num);
		add_dis(dis_sum.u1diff,dis.u1diff,time_num);add_dis(dis_sum.u2diff,dis.u2diff,time_num);
		add_dis(dis_sum.v1diff,dis.v1diff,time_num);add_dis(dis_sum.v2diff,dis.v2diff,time_num);
		add_dis(dis_sum.w1diff,dis.w1diff,time_num);add_dis(dis_sum.w2diff,dis.w2diff,time_num);
		add_dis2D(dis_sum.uTsum2D,dis.uTsum2D,time_num,position_num);
		add_dis2D(dis_sum.vTsum2D,dis.vTsum2D,time_num,position_num);
		add_dis2D(dis_sum.wTsum2D,dis.wTsum2D,time_num,position_num);
		add_dis2D(dis_sum.pattern,dis.pattern,64,max_hit);

		free_dis_Tsum(dis,time_num);

	}//End of for(num_file=0;num_file<param.num_files;num_file++)

	FILE *fp_summary;
	char summary[MAX_LINE];
	sprintf(summary,"%s%s_sum_TsumSummary.txt",param.workspace,param.name[0]);
	fp_summary=fopen(summary,"w+");
	fprintf(fp_summary,"name totalcount u1count u2count v1count v2count w1count w2count u12count u12effcount v12count v12effcount w12count w12effcount\n");
	for(int icount=0;icount<param.num_files;icount++)
	{
		fprintf(fp_summary,"%s %d %d %d %d %d %d %d %d %d %d %d %d %d\n",param.name[icount],
			totalcount[icount],u1count[icount],u2count[icount],v1count[icount],v2count[icount],w1count[icount],w2count[icount],
			u12count[icount],u12effcount[icount],v12count[icount],v12effcount[icount],w12count[icount],w12effcount[icount]);
	}
	fprintf(fp_summary,"Total %d %d %d %d %d %d %d %d %d %d %d %d %d\n",totalcount[param.num_files],u1count[param.num_files],u2count[param.num_files],v1count[param.num_files],
		v2count[param.num_files],w1count[param.num_files],w2count[param.num_files],u12count[param.num_files],u12effcount[param.num_files],v12count[param.num_files],
		v12effcount[param.num_files],w12count[param.num_files],w12effcount[param.num_files]);

	double area[3],Tsum[3],delat_T[3];
	fprintf(fp_summary,"\n");
	if(gauss_fit(area[0],Tsum[0],delat_T[0],dis_sum.uTsum,param.time_min,param.time_max,param.time_step))
		fprintf(fp_summary,"uTsum %.2f %.2f %.2f %.4f\n",area[0],Tsum[0],delat_T[0],104/Tsum[0]*20);
	if(gauss_fit(area[1],Tsum[1],delat_T[1],dis_sum.vTsum,param.time_min,param.time_max,param.time_step))
		fprintf(fp_summary,"vTsum %.2f %.2f %.2f %.4f\n",area[1],Tsum[1],delat_T[1],104/Tsum[1]*20);
	if(gauss_fit(area[2],Tsum[2],delat_T[2],dis_sum.wTsum,param.time_min,param.time_max,param.time_step))
		fprintf(fp_summary,"wTsum %.2f %.2f %.2f %.4f\n",area[2],Tsum[2],delat_T[2],104/Tsum[2]*20);
	
	fclose(fp_summary);

	if(param.num_files>1)
	{
		FILE *fp_sumpattern,*fp_sumTsum,*fp_sumTsumDis;
		char sumpattern_name[MAX_LINE],sumTsum_name[MAX_LINE],sumTsumDis_name[MAX_LINE];
		sprintf(sumpattern_name,"%s%s_sum_pattern.txt",param.workspace,param.name[0]);
		sprintf(sumTsum_name,"%s%s_sum_Tsum.txt",param.workspace,param.name[0]);
		sprintf(sumTsumDis_name,"%s%s_sum_Tsum2D.txt",param.workspace,param.name[0]);		

		fp_sumpattern=fopen(sumpattern_name,"w+");
		fp_sumTsum=fopen(sumTsum_name,"w+");
		fp_sumTsumDis=fopen(sumTsumDis_name,"w+");		

		for(int i=0;i<64;i++)
		{
			fprintf(fp_sumpattern,"%d %d %d %d\n",i,dis_sum.pattern[i][0],dis_sum.pattern[i][1],dis_sum.pattern[i][2]);
		}

		fprintf(fp_sumTsum,"time_mid t uTsum vTsum wTsum u1 u2 v1 v2 w1 w2 u12 u22 v12 v22 w12 w22 MCP2 u1diff u2diff v1diff v1diff w1diff w2diff\n");
		for(int j=0;j<time_num;j++)
		{
			double time_mid=param.time_min+(j+0.5)*param.time_step;
			double t=time_mid*0.025;
			fprintf(fp_sumTsum,"%.1f %.3f %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
				time_mid,t,dis_sum.uTsum[j],dis_sum.vTsum[j],dis_sum.wTsum[j],
				dis_sum.u11[j],dis_sum.u21[j],dis_sum.v11[j],dis_sum.v21[j],dis_sum.w11[j],dis_sum.w21[j],
				dis_sum.u12[j],dis_sum.u22[j],dis_sum.v12[j],dis_sum.v22[j],dis_sum.w12[j],dis_sum.w22[j],dis_sum.MCP2[j],
				dis_sum.u1diff[j],dis_sum.u2diff[j],dis_sum.v1diff[j],dis_sum.v2diff[j],dis_sum.w1diff[j],dis_sum.w2diff[j]);
		}
		for(int iT=0;iT<time_num;iT++)
		{
			for(int k=0;k<position_num;k++)
			{
				double time_mid=param.time_min+(iT+0.5)*param.time_step;
				double position_mid=param.position_min+(k+0.5)*param.position_step;
				fprintf(fp_sumTsumDis,"%.1f %.1f %d %d %d\n",time_mid,position_mid,dis_sum.uTsum2D[iT][k],dis_sum.vTsum2D[iT][k],dis_sum.wTsum2D[iT][k]);
			}
		}
		fclose(fp_sumpattern);
		fclose(fp_sumTsum);
		fclose(fp_sumTsumDis);
	}//End of if(param.num_files>1)

	free_dis_Tsum(dis_sum,position_num);

	return 1;
}
