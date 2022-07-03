/**********************
版本1.4更新
1.readevent()函数：max_hit限制优化，防止ch[channel][e.nhit[channel]]溢出.
2.readevent()函数：bug修正，judge[max_hit]->judge[6].
3.gauss_fit()函数修改bug
****************************/

#define ch_none -1

int readevent(event &e,FILE *fp_in); //从文件中读出一个事例，赋值到e中
int fill_int(int x,int *dis,int min,int max,int step); //把数x添加到dis分布中
int fill_double(double x,int *dis,double min,double max,double step); //把数x添加到dis分布中
int fill2D_int(int x,int y,int **dis,int xmin,int xmax,int xstep,int ymin,int ymax,int ystep);
int fill2D_double(double x,double y,int **dis,double xmin,double xmax,double xstep,double ymin,double ymax,double ystep);  //把点(x,y)添加到dis分布中
int initialize_event(event &e);  //初始化参数结构体e
int initialize_particle(particle &p);   //初始化粒子结构体p
int initialize_int2(int **&dis,int n1,int n2); //初始化二维分布dis (n1*n2)
int initialize_int1(int *&dis,int n1); //初始化一维分布dis
int rotatexy(double &x,double &y,double angle); //将(x,y)绕圆心转angle角度
double getangle(double x,double y);  //(x,y)算角度
double getradius(double x,double y);  //(x,y)算半径
int getposition(particle &p,parameter param); //由u,v,w算粒子的位置
char * mode_name(int mode); //由Analysis_Mode的值给出分析模式的名称
int getenergy(particle &p,parameter param);  //由位置算能量
double getIP(particle p1,particle p2,parameter param);  //计算电离能
double getAzimuth(particle p1,particle p2,parameter param); //计算相对方位角
int testTsum(int Tsum,char layer,parameter param);  //判断Tsum值是否在正确范围内
int completesignal(int pattern); //判断信号是否完整
int truecoincidence(double t,parameter param); //判断时间差是否在真符合窗内
int accidentalcoincidence(double t,parameter param); //判断时间差是否在偶然符合窗内
int getcenter(double &x0,double &y0,double &r0,int **xydis,parameter param); //由(x,y)的分布拟合圆心
int GaussXiaoYuan(double **A,double *X,double *B,int N); //高斯消元法求解n元一次方程组AX=B
double polynomial_fit(double *x,double *y,int m,double *c,int n); //n阶多项式拟合得到系数c
int gauss_newton_fit(double &A,double &B,double &C,double *y,double *x,int n);
//newton迭代法高斯拟合y=exp(Ax^2+Bx+C)，A,B,C的初始值要在正确值附近
int gauss_fit(double &area,double &x0,double &width,int *dis,double min,double max,double step);
//高斯拟合dis分布
void swap(int &A,int &B); //交换A与B
void swap_double(double &A,double &B); //交换A与B
void swap_particle(particle &p1,particle &p2); //交换p1与p2
int add_dis(int *des,int *org,int n); //将分布org添加到des中去
int add_dis2D(int **des,int **org,int n1,int n2); //将分布org添加到des中去

int readevent(event &e,FILE *fp_in)
{
		uint32 raw_dec=0;
		int type=0;       //word type: global header, global trailer...
		int flag=0;  //to show that a complete global header and trailer
		int channel=0;   //TDC channel
		int judge[max_ch]={0}; //judge whether there is hitting
		int overflow;//judge whether the count is over one period (52us)
		int ch[max_ch][max_hit]={0}; //每个通道的TDC计数

		while(fread(&raw_dec,4,1,fp_in))
		{
			//type=raw_dec/0x8000000;  //2^27=134217728=0x8000000,bits[31:27],word type;
			type=raw_dec>>27;
			if(type==8)              //Global header
			{
				//e.count=(raw_dec%0x8000000)/32;  //%2^27/2^5
			    e.count=(raw_dec&0x07FFFFE0)>>5;
				if(flag==1)
				{
					for(int j=0;j<max_ch;j++)
					{
						e.nhit[j]=0;
						for(int k=0;k<max_hit;k++)	ch[j][k]=0;
					}
					cout<<"data error: event count "<<e.count-1<<endl; //two helder together
				}    
				flag=1;         //there is a global header
			}
			else if(type==0)
			{
				if(flag==0)  continue;  //no helder
//				channel=(raw_dec%67108864)/2097152;
				//channel=(raw_dec%0x4000000)/0x200000;    //TDC channel
				channel=(raw_dec>>21)&0x0000001F;
				if(channel>=max_ch)
				{
					printf("event %d channel %d error\n",e.count,channel);
					continue;
				}
				if(e.nhit[channel]<max_hit)	ch[channel][e.nhit[channel]]=(raw_dec&0x001FFFFF);
				e.nhit[channel]++;
				if((raw_dec&0x001FFFFF)>(0x200000-60*1024))   //60*1024channel->1.5us
					overflow=1;
			}
			else if(type==16)
			{
				if(flag==0)	continue;
				e.status=(raw_dec>>24)&0x00000007;
				flag=0;            //global trailer, a EVENT is over.
				break;
			}
			if(feof(fp_in)) break;
		}

		if(overflow==1)
			for(int j=0;j<max_ch;j++)
			{
				for(int k=0;k<e.nhit[j]&&k<max_hit;k++)
					if(ch[j][k]<1024*120)
						ch[j][k]+=0x200000;
			}

		int t_zero=ch[MCPch][0];
		for(int j=0;j<max_ch;j++)
		{
			int ihit;
			for(ihit=0;ihit<e.nhit[j]&&ihit<max_hit;ihit++)
				ch[j][ihit]-=t_zero;
			for(;ihit<max_hit;ihit++)
				ch[j][ihit]=ch_none;
		}
		e.u1[0]=ch[u1ch][0];e.u1[1]=ch[u1ch][1];e.u2[0]=ch[u2ch][0];e.u2[1]=ch[u2ch][1];
		e.v1[0]=ch[v1ch][0];e.v1[1]=ch[v1ch][1];e.v2[0]=ch[v2ch][0];e.v2[1]=ch[v2ch][1];
		e.w1[0]=ch[w1ch][0];e.w1[1]=ch[w1ch][1];e.w2[0]=ch[w2ch][0];e.w2[1]=ch[w2ch][1];
		e.MCP[0]=ch[MCPch][0];e.MCP[1]=ch[MCPch][1];
		int j;
		for(int i=0;i<max_hit;i++)
		{
			j=u1ch;
			judge[j]=e.nhit[j]>i?(1<<0):0;	e.pattern[i]+=judge[j];	
			j=u2ch;
			judge[j]=e.nhit[j]>i?(1<<1):0;	e.pattern[i]+=judge[j];	
			j=v1ch;
			judge[j]=e.nhit[j]>i?(1<<2):0;	e.pattern[i]+=judge[j];	
			j=v2ch;
			judge[j]=e.nhit[j]>i?(1<<3):0;	e.pattern[i]+=judge[j];	
			j=w1ch;
			judge[j]=e.nhit[j]>i?(1<<4):0;	e.pattern[i]+=judge[j];	
			j=w2ch;
			judge[j]=e.nhit[j]>i?(1<<5):0;	e.pattern[i]+=judge[j];	
		}
		return 1;
}

int fill_int(int x,int *dis,int min,int max,int step)
{
	int i;
	if(x<min||x>max)	return 0;
	i=(x-min)/step;	dis[i]++;
	return 1;
}

int fill_double(double x,int *dis,double min,double max,double step)
{
	int i;
	if(x<min||x>max)	return 0;
	i=(int)((x-min)/step);
	dis[i]++;
	return 1;
}

int fill2D_int(int x,int y,int **dis,int xmin,int xmax,int xstep,int ymin,int ymax,int ystep)
{
	int i,j;
	if(x<xmin||x>xmax||y<ymin||y>ymax)	return 0;
	i=(x-xmin)/xstep;	j=(y-ymin)/ystep;
	dis[i][j]++;
	return 1;
}

int fill2D_double(double x,double y,int **dis,double xmin,double xmax,double xstep,double ymin,double ymax,double ystep)
{
	int i,j;
	if(x<xmin||x>xmax||y<ymin||y>ymax)	return 0;
	i=(int)((x-xmin)/xstep);	j=(int)((y-ymin)/ystep);
	dis[i][j]++;
	return 1;
}

int initialize_event(event &e)
{
	int i;
	e.count=0;
	for(i=0;i<3;i++)
	{
		e.u1[i]=0;e.u2[i]=0;e.v1[i]=0;e.v2[i]=0;e.w1[i]=0;e.w2[i]=0;e.MCP[i]=0;
		e.pattern[i]=0;
	}
	for(i=0;i<max_ch;i++)
	{
		e.nhit[i]=0;
	}
	e.uTsum=0;e.vTsum=0;e.wTsum=0;
	e.status; 
	return 1;
}

int initialize_particle(particle &p)
{
	p.u=p.v=p.w=100;
	p.uTsum=p.vTsum=p.wTsum=-1;
	p.uflag=p.vflag=p.wflag=0;
	p.x=p.xuv=p.xuw=p.xvw=100;
	p.y=p.yuv=p.yuw=p.yvw=100;
	p.r=p.ruv=p.ruw=p.rvw=100;
	p.phi=p.phiuv=p.phiuw=p.phivw=-1;
	p.energy=p.euv=p.euw=p.evw=-1;
	p.angle=p.auv=p.auw=p.avw=-1;
	p.whichtwolayer=0;
	return 1;
}

int initialize_int2(int **&dis,int n1,int n2)
{
	dis=(int **)malloc(sizeof(int *)*n1);
	for(int i=0;i<n1;i++)
	{
		dis[i]=(int *)malloc(sizeof(int)*n2);
		for(int j=0;j<n2;j++)
			dis[i][j]=0;
	}
	return 1;
}

int initialize_int1(int *&dis,int n1)
{
	dis=(int *)malloc(sizeof(int)*n1);
	for(int j=0;j<n1;j++)
		dis[j]=0;
	return 1;
}

int rotatexy(double &x,double &y,double angle)
{
	double x0,y0;
	x0=x*cos(angle*3.141592654/180)+y*sin(angle*3.141592654/180);
	y0=-x*sin(angle*3.141592654/180)+y*cos(angle*3.141592654/180);
	x=x0;
	y=y0;
	return 1;
}

double getangle(double x,double y)
{
	double angle;
	angle=atan(y/x)*180/3.141592654;
	if(x<0) angle+=180;
	if(x>=0&&y<0) angle+=360;
	angle-=90;
	if(angle<0)
		angle+=360;
	if(angle>=360)
		angle-=360;
	return angle;
}

double getradius(double x,double y)
{
	return sqrt(x*x+y*y);
}

int getposition(particle &p,parameter param)
{
	double angle=param.rotate_angle;
	int flag=0;
	if(p.vflag&&p.wflag)
	{
		p.xvw=p.v+p.w;
		p.yvw=(p.w-p.v)/sqrt(3);
		rotatexy(p.xvw,p.yvw,angle);
		if(p.xvw<0){ p.xvw+=param.x1_energy_offset;p.yvw+=param.y1_energy_offset;}
		else{ p.xvw+=param.x2_energy_offset;p.yvw+=param.y2_energy_offset;}
		p.rvw=getradius(p.xvw,p.yvw);
		p.phivw=getangle(p.xvw,p.yvw);
		double p_rvw=getradius(p.xvw+param.x1_angle_offset,p.yvw+param.y1_angle_offset);
		double p_phivw=getangle(p.xvw+param.x1_angle_offset,p.yvw+param.y1_angle_offset);
		if(p.rvw<=param.radius_max)
		{
			p.x=p.xvw;	p.y=p.yvw; 
			p.r=p.rvw; p.phi=p.phivw;
			if(p.xvw<0)
				p.avw=-9.41501 + 0.0487752*p_rvw + 0.0145034*pow(p_rvw,2) - 0.000519285*pow(p_rvw,3) + 4.87776e-6*pow(p_rvw,4) + 0.910576*p_phivw - 0.0100772*p_rvw*p_phivw + 
   0.00024438*pow(p_rvw,2)*p_phivw - 1.5545199999999998e-6*pow(p_rvw,3)*p_phivw + 0.00430255*pow(p_phivw,2) + 2.17247e-6*p_rvw*pow(p_phivw,2) - 
   4.0536e-7*pow(p_rvw,2)*pow(p_phivw,2) - 0.0000271525*pow(p_phivw,3) + 1.0256799999999998e-7*p_rvw*pow(p_phivw,3) + 4.36416e-8*pow(p_phivw,4);
			else
				p.avw=-260.824 + 0.158682*p_rvw + 0.0243335*pow(p_rvw,2) - 0.000688468*pow(p_rvw,3) + 5.39968e-6*pow(p_rvw,4) + 4.61658*p_phivw - 0.0102905*p_rvw*p_phivw + 
   0.0000725447*pow(p_rvw,2)*p_phivw - 8.11871e-8*pow(p_rvw,3)*p_phivw - 0.0180425*pow(p_phivw,2) + 0.00003709*p_rvw*pow(p_phivw,2) - 
   1.34555e-7*pow(p_rvw,2)*pow(p_phivw,2) + 0.0000380643*pow(p_phivw,3) - 3.90656e-8*p_rvw*pow(p_phivw,3) - 2.8867600000000002e-8*pow(p_phivw,4);
				/*-256.824 + 0.158682*p_rvw + 0.0243335*pow(p_rvw,2) - 0.000688468*pow(p_rvw,3) + 5.39968e-6*pow(p_rvw,4) + 4.61658*p_phivw - 0.0102905*p_rvw*p_phivw + 
   0.0000725447*pow(p_rvw,2)*p_phivw - 8.11871e-8*pow(p_rvw,3)*p_phivw - 0.0180425*pow(p_phivw,2) + 0.00003709*p_rvw*pow(p_phivw,2) - 
   1.34555e-7*pow(p_rvw,2)*pow(p_phivw,2) + 0.0000380643*pow(p_phivw,3) - 3.90656e-8*p_rvw*pow(p_phivw,3) - 2.8867600000000002e-8*pow(p_phivw,4);*/
			p.angle=p.avw;
			flag=3;
		}
	}
	if(p.uflag&&p.wflag)
	{
		p.xuw=p.u;
		p.yuw=(2*p.w-p.u)/sqrt(3);
		rotatexy(p.xuw,p.yuw,angle);
		if(p.xuw<0){ p.xuw+=param.x1_energy_offset;p.yuw+=param.y1_energy_offset;}
		else{ p.xuw+=param.x2_energy_offset;p.yuw+=param.y2_energy_offset;}
		p.ruw=getradius(p.xuw,p.yuw);
		p.phiuw=getangle(p.xuw,p.yuw);

		double p_ruw=getradius(p.xuw+param.x1_angle_offset,p.yuw+param.y1_angle_offset);
		double p_phiuw=getangle(p.xuw+param.x1_angle_offset,p.yuw+param.y1_angle_offset);

		if(p.ruw<=param.radius_max)
		{
			p.x=p.xuw; p.y=p.yuw; 
			p.r=p.ruw; p.phi=p.phiuw;
			if(p.xuw<0)
				p.auw=-1.40317 - 1.50247*p_ruw + 0.075774*pow(p_ruw,2) - 0.00147571*pow(p_ruw,3) + 0.0000101647*pow(p_ruw,4) + 1.11586*p_phiuw + 0.00255882*p_ruw*p_phiuw - 
   0.000215814*pow(p_ruw,2)*p_phiuw + 2.3264900000000002e-6*pow(p_ruw,3)*p_phiuw - 0.00160446*pow(p_phiuw,2) + 0.0000198*p_ruw*pow(p_phiuw,2) + 
   2.3640599999999998e-7*pow(p_ruw,2)*pow(p_phiuw,2) + 0.0000100525*pow(p_phiuw,3) - 1.16608e-7*p_ruw*pow(p_phiuw,3) - 1.79618e-8*pow(p_phiuw,4);
			else
				p.auw=102.377 - 0.451657*p_ruw + 0.0463457*pow(p_ruw,2) - 0.00106063*pow(p_ruw,3) + 6.46286e-6*pow(p_ruw,4) - 0.528989*p_phiuw - 0.00647507*p_ruw*p_phiuw - 
   0.0000111002*pow(p_ruw,2)*p_phiuw + 1.07914e-6*pow(p_ruw,3)*p_phiuw + 0.00892433*pow(p_phiuw,2) + 0.0000301017*p_ruw*pow(p_phiuw,2) - 
   1.30926e-7*pow(p_ruw,2)*pow(p_phiuw,2) - 0.0000237687*pow(p_phiuw,3) - 3.16261e-8*p_ruw*pow(p_phiuw,3) + 2.3231e-8*pow(p_phiuw,4);
				
				/*106.377 - 0.451657*p_ruw + 0.0463457*pow(p_ruw,2) - 0.00106063*pow(p_ruw,3) + 6.46286e-6*pow(p_ruw,4) - 0.528989*p_phiuw - 0.00647507*p_ruw*p_phiuw - 
   0.0000111002*pow(p_ruw,2)*p_phiuw + 1.07914e-6*pow(p_ruw,3)*p_phiuw + 0.00892433*pow(p_phiuw,2) + 0.0000301017*p_ruw*pow(p_phiuw,2) - 
   1.30926e-7*pow(p_ruw,2)*pow(p_phiuw,2) - 0.0000237687*pow(p_phiuw,3) - 3.16261e-8*p_ruw*pow(p_phiuw,3) + 2.3231e-8*pow(p_phiuw,4);*/
			p.angle=p.auw;
			flag=2;
		}
	}
	if(p.uflag&&p.vflag)
	{
		p.xuv=p.u;
		p.yuv=(p.u-2*p.v)/sqrt(3);
		rotatexy(p.xuv,p.yuv,angle);
		if(p.xuv<0){ p.xuv+=param.x1_energy_offset;p.yuv+=param.y1_energy_offset;}
		else{ p.xuv+=param.x2_energy_offset;p.yuv+=param.y2_energy_offset;}
		p.ruv=getradius(p.xuv,p.yuv);
		p.phiuv=getangle(p.xuv,p.yuv);

		double p_ruv=getradius(p.xuv+param.x1_angle_offset,p.yuv+param.y1_angle_offset);
		double p_phiuv=getangle(p.xuv+param.x1_angle_offset,p.yuv+param.y1_angle_offset);

		if(p.ruv<=param.radius_max)
		{
			p.x=p.xuv; p.y=p.yuv; 
			p.r=p.ruv; p.phi=p.phiuv;
			if(p.xuv<0)
				p.auv=-12.9852 - 0.114795*p_ruv + 0.0160625*pow(p_ruv,2) - 0.000341482*pow(p_ruv,3) + 2.22427e-6*pow(p_ruv,4) + 1.13288*p_phiuv - 0.00540993*p_ruv*p_phiuv - 
   7.692129999999999e-6*pow(p_ruv,2)*p_phiuv + 6.995269999999999e-7*pow(p_ruv,3)*p_phiuv + 0.000924862*pow(p_phiuv,2) + 0.0000345653*p_ruv*pow(p_phiuv,2) - 
   1.77194e-7*pow(p_ruv,2)*pow(p_phiuv,2) - 0.0000156238*pow(p_phiuv,3) - 5.4844e-8*p_ruv*pow(p_phiuv,3) + 4.7242300000000003e-8*pow(p_phiuv,4);
			else
				p.auv=-237.323 - 0.71412*p_ruv + 0.0281956*pow(p_ruv,2) - 0.000520176*pow(p_ruv,3) + 4.796579999999999e-6*pow(p_ruv,4) + 4.83136*p_phiuv - 0.000802809*p_ruv*p_phiuv - 
   0.0000484436*pow(p_ruv,2)*p_phiuv - 2.39067e-7*pow(p_ruv,3)*p_phiuv - 0.0221689*pow(p_phiuv,2) + 0.0000167148*p_ruv*pow(p_phiuv,2) + 
   1.38611e-7*pow(p_ruv,2)*pow(p_phiuv,2) + 0.0000543294*pow(p_phiuv,3) - 3.65796e-8*p_ruv*pow(p_phiuv,3) - 4.83401e-8*pow(p_phiuv,4);
				
/*				-233.323 - 0.71412*p_ruv + 0.0281956*pow(p_ruv,2) - 0.000520176*pow(p_ruv,3) + 4.796579999999999e-6*pow(p_ruv,4) + 4.83136*p_phiuv - 0.000802809*p_ruv*p_phiuv - 
   0.0000484436*pow(p_ruv,2)*p_phiuv - 2.39067e-7*pow(p_ruv,3)*p_phiuv - 0.0221689*pow(p_phiuv,2) + 0.0000167148*p_ruv*pow(p_phiuv,2) + 
   1.38611e-7*pow(p_ruv,2)*pow(p_phiuv,2) + 0.0000543294*pow(p_phiuv,3) - 3.65796e-8*p_ruv*pow(p_phiuv,3) - 4.83401e-8*pow(p_phiuv,4)   */
			p.angle=p.auv;
			flag=1;
		}
	}
	p.whichtwolayer=flag;
	return flag;
}
char * mode_name(int mode)
{
	static char str[50];
	switch(mode)
	{
	case 1: strcpy(str,"getTsum");break;
	case 2: strcpy(str,"elastic scattering");break;
	case 3: strcpy(str,"ems measurement");break;
	case 4: strcpy(str,"find the angle center");break;
	case 5: strcpy(str,"non coincidence mode");break;
	case 6: strcpy(str,"direct (x,y) output");break;
	default: strcpy(str,"error");
	}
	return str;
}

int getenergy(particle &p,parameter param)
{
	if(p.r>param.radius_max) 
		return 0;
	int i=(int)(p.phi/360.0*param.cali_angle_num);
	if(p.whichtwolayer==1)
		p.energy=param.cuv0[i]+param.cuv1[i]*p.ruv+param.cuv2[i]*p.ruv*p.ruv;
	else if(p.whichtwolayer==2)
		p.energy=param.cuw0[i]+param.cuw1[i]*p.ruw+param.cuw2[i]*p.ruw*p.ruw;
	else if(p.whichtwolayer==3)
		p.energy=param.cvw0[i]+param.cvw1[i]*p.rvw+param.cvw2[i]*p.rvw*p.rvw;
	else
		return 0;
	return 1;
}

double getIP(particle p1,particle p2,parameter param)
{
	double IP;
	IP=param.E0+param.E0_shift-p1.energy-p2.energy;
	return IP;
}

double getAzimuth(particle p1,particle p2,parameter param)
{
	double phi;	
//	int i,j;
/*	p1.phi=0;p2.phi=0;
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			if(p1.whichtwolayer==1) p1.phi+=param.auv1[i][j]*pow(p1.ruv,i)*pow(p1.phiuv,j);
			if(p1.whichtwolayer==2) p1.phi+=param.auw1[i][j]*pow(p1.ruw,i)*pow(p1.phiuw,j);
			if(p1.whichtwolayer==3) p1.phi+=param.avw1[i][j]*pow(p1.rvw,i)*pow(p1.phivw,j);
			if(p2.whichtwolayer==1) p2.phi+=param.auv2[i][j]*pow(p2.ruv,i)*pow(p2.phiuv,j);
			if(p2.whichtwolayer==2) p2.phi+=param.auw2[i][j]*pow(p2.ruw,i)*pow(p2.phiuw,j);
			if(p2.whichtwolayer==3) p2.phi+=param.avw2[i][j]*pow(p2.rvw,i)*pow(p2.phivw,j);
		}
	}
*/
//	phi=p2.angle*1.05-p1.angle*1.04-180-12;
	phi=p2.angle-p1.angle-180;
	return phi;
}

int testTsum(int Tsum,char layer,parameter param)
{
    int Tsum_Lower,Tsum_Upper;
    switch(layer)
	{
	case 'u': Tsum_Lower=param.uTsum_Lower;Tsum_Upper=param.uTsum_Up;break;
	case 'v': Tsum_Lower=param.vTsum_Lower;Tsum_Upper=param.vTsum_Up;break;
	case 'w': Tsum_Lower=param.wTsum_Lower;Tsum_Upper=param.wTsum_Up;break;
	default: return -1;
	}
	if(Tsum<=Tsum_Upper&&Tsum>=Tsum_Lower)
        return 1;
	else
        return 0;
}

int completesignal(int pattern)
{
    if(pattern==15 || pattern ==31 || pattern==47 || pattern==51 || pattern==55 || pattern>=59)
        return 1;
	else
        return 0;
}

int truecoincidence(double t,parameter param)
{
	if(t<=param.truecoincidence2&&t>=param.truecoincidence1)
        return 1;
	else
        return 0;
}

int accidentalcoincidence(double t,parameter param)
{
	t=fabs(t);
	if(t<=param.accidentalcoincidence2&&t>=param.accidentalcoincidence1)
        return 1;
	else
        return 0;
}

int getcenter(double &x0,double &y0,double &r0,int **xydis,parameter param)  //拟合圆心
{
		double x11,x12,x13,x21,x22,x23,x31,x32,x33,y1,y2,y3,a,b,c;//系数
		x11=0;x12=0;x13=0;x21=0;x22=0;x23=0;x31=0;x32=0;x33=0;y1=0;y2=0;y3=0;
		int count=0,max=0;

		int num=(int)((param.position_max-param.position_min)/param.position_step);
		double x,y;
		for(int i=0;i<num;i++)
		{
			for(int j=0;j<num;j++)
			{
				if(max<xydis[i][j])
					max=xydis[i][j];
			}
		}
		for(int i=0;i<num;i++)
		{
			for(int j=0;j<num;j++)
			{
				if(xydis[i][j]>max/10)
				{
					x=param.position_min+param.position_step*(i+0.5);
					y=param.position_min+param.position_step*(j+0.5);
					x11+=x*x*xydis[i][j];
					x12+=x*y*xydis[i][j];
					x13+=x*xydis[i][j];
					y1+=(x*x*x+x*y*y)*xydis[i][j];
					x22+=y*y*xydis[i][j];
					x23+=y*xydis[i][j];
					y2+=(x*x*y+y*y*y)*xydis[i][j];
					count+=xydis[i][j];
				}				
			}
		}
		x21=x12;
		x31=x13;
		x32=x23;
		x33=count;
		y3=x11+x22;
		a=-((-(x23*x32*y1)+x22*x33*y1+x13*x32*y2-x12*x33*y2-x13*x22*y3+x12*x23*y3)/(-(x13*x22*x31)+x12*x23*x31+x13*x21*x32-x11*x23*x32-x12*x21*x33+x11*x22*x33));
		b=-((-(x23*x31*y1)+x21*x33*y1+x13*x31*y2-x11*x33*y2-x13*x21*y3+x11*x23*y3)/(x13*x22*x31-x12*x23*x31-x13*x21*x32+x11*x23*x32+x12*x21*x33-x11*x22*x33));
		c=-((-(x22*x31*y1)+x21*x32*y1+x12*x31*y2-x11*x32*y2-x12*x21*y3+x11*x22*y3)/(-(x13*x22*x31)+x12*x23*x31+x13*x21*x32-x11*x23*x32-x12*x21*x33+x11*x22*x33));
		x0=-a/2;
		y0=-b/2;
		r0=sqrt(a*a/4+b*b/4-c);
		return 1;
}

int GaussXiaoYuan(double **A,double *X,double *B,int N)
{	//高斯列主元消元法求解方程组AX=B
	double **a,*b,temp,s;
	int i,j,k;
	a=(double **)malloc(sizeof(double *)*N);
	b=(double *)malloc(sizeof(double)*N);
	for(i=0;i<N;i++)
	{
		a[i]=(double *)malloc(sizeof(double)*N);
	}
	
	for(i=0;i<N;i++)
	{
		b[i]=B[i];
		for(j=0;j<N;j++)	a[i][j]=A[i][j];
	}
	for(i=0;i<N-1;i++)
	{
		s=fabs(a[i][i]);k=i;
		for(j=i+1;j<N;j++)
		{
			if(s<fabs(a[j][i]))
			{
				s=fabs(a[j][i]);k=j;
			}
		}
		if(k!=i)
		{
			for(j=i;j<N;j++)
			{
				temp=a[i][j];a[i][j]=a[k][j];a[k][j]=temp;
			}
			temp=b[i];b[i]=b[k];b[k]=temp;
		}
		for(j=i+1;j<N;j++)
		{
			temp=a[j][i]/a[i][i];
			for(k=i;k<N;k++)	a[j][k]-=a[i][k]*temp;
			b[j]-=temp*b[i];
		}
	}

	for(i=N-1;i>=0;i--)
	{
		temp=0;
		for(j=i+1;j<N;j++)	temp+=a[i][j]*X[j];
		X[i]=(b[i]-temp)/a[i][i];
	}
	return 1;
}

double polynomial_fit(double *x,double *y,int m,double *c,int n)
{ //m是多项式的阶数，n是数据个数。
	if(m>=n)
	{
		for(int i=0;i<m+1;i++)
			c[i]=0;
		return 0;
	}
	double **A,*B,R2;
	A=(double **)malloc(sizeof(double *)*(m+1));
	for (int i=0;i<m+1;i++)
		A[i]=(double *)malloc(sizeof(double)*(m+1));
	B=(double *)malloc(sizeof(double)*(m+1));

	for(int i=0;i<m+1;i++)
	{
		for(int j=0;j<m+1;j++)
		{
			A[i][j]=0;
			for(int k=0;k<n;k++)
				A[i][j]+=pow(x[k],i)*pow(x[k],j);
		}
		B[i]=0;
		for(int k=0;k<n;k++)
			B[i]+=pow(x[k],i)*y[k];
	}
	GaussXiaoYuan(A,c,B,m+1);
	
	double y2=0,y2r=0;
	for(int k=0;k<n;k++)
	{
		y2+=y[k]*y[k];
		double yr=y[k];
		for(int i=0;i<m+1;i++)
			yr-=c[i]*pow(x[k],i);
		y2r+=yr*yr;
	}
	R2=1-y2r/y2;
	return R2;
} 

int gauss_newton_fit(double &A,double &B,double &C,double *y,double *x,int n)
{//y=exp(A*x^2+B*x+C)
	double **J,delta[3]={0},f[3]={0};
	J=(double **)malloc(sizeof(double *)*3);
	J[0]=(double *)malloc(sizeof(double)*3);
	J[1]=(double *)malloc(sizeof(double)*3);
	J[2]=(double *)malloc(sizeof(double)*3);
	for(int count=0;count<128;count++)
	{
		
		J[0][0]=J[0][1]=J[0][2]=J[1][0]=J[1][1]=J[1][2]=J[2][0]=J[2][1]=J[2][2]=0;
		f[0]=f[1]=f[2]=0;
		double yi;
		for(int i=0;i<n;i++)
		{
			yi=exp(A*x[i]*x[i]+B*x[i]+C);
			f[0]-=(yi-y[i])*yi*x[i]*x[i];
			f[1]-=(yi-y[i])*yi*x[i];
			f[2]-=(yi-y[i])*yi;
			J[0][0]+=(2*yi-y[i])*yi*pow(x[i],4);
			J[0][1]+=(2*yi-y[i])*yi*pow(x[i],3);
			J[0][2]+=(2*yi-y[i])*yi*pow(x[i],2);
			J[1][2]+=(2*yi-y[i])*yi*pow(x[i],1);
			J[2][2]+=(2*yi-y[i])*yi;
		}
		J[1][0]=J[0][1];
		J[1][1]=J[2][0]=J[0][2];
		J[2][1]=J[1][2];
		GaussXiaoYuan(J,delta,f,3);

		A+=delta[0];
		B+=delta[1];
		C+=delta[2];
	
		if(fabs(delta[0]/A)<1e-4 && fabs(delta[1]/B)<1e-4 && fabs(delta[2]/C)<1e-4 )
		{
			return 1;
		}
	}
	free(J[0]);free(J[1]);free(J[2]);free(J);
	return 0;
}

int gauss_fit(double &area,double &x0,double &width,int *dis,double min,double max,double step)
{
	int num=(int)((max-min)/step);
	int i_max=0,dis_max=0,i_start=0,i_stop=0;
	double *dis_log,*xdis,*x_cord;
	double A[3];
	for(int i=0;i<num;i++)
	{
		if(dis_max<dis[i])
		{
			dis_max=dis[i];
			i_max=i;
		}
	}
	for(int i=i_max;i<num;i++)
	{		
		if(dis[i]<=(dis_max/20+1))	
			break;
		i_stop=i;
	}
	for(int i=i_max;i>=0;i--)
	{		
		if(dis[i]<=(dis_max/20+1))	
			break;
		i_start=i;
	}
	if((i_stop-i_start)<=3)
		return 0;
	dis_log=(double *)malloc(sizeof(double)*(i_stop-i_start+1));
	x_cord=(double *)malloc(sizeof(double)*(i_stop-i_start+1));
	xdis=(double *)malloc(sizeof(double)*(i_stop-i_start+1));
	for(int i=i_start;i<=i_stop;i++)	{dis_log[i-i_start]=log(dis[i]);xdis[i-i_start]=dis[i]; x_cord[i-i_start]=min+step*(i+0.5);}
	polynomial_fit(x_cord,dis_log,2,A,(i_stop-i_start+1));   //lny=A[2]x^2+A[1]*x+A[0];
	if(A[2]>=0) return 0;
	if(gauss_newton_fit(A[2],A[1],A[0],xdis,x_cord,(i_stop-i_start+1))==0)
		return 0;
	x0=-A[1]/2.0/A[2];
	width=sqrt(-0.5/A[2])*2.35482;
	area=exp(A[0]-A[1]*A[1]/A[2]/4.0)*sqrt(2*3.14159)*sqrt(-0.5/A[2]);
	return 1;
}

void swap(int &A,int &B)
{
	int temp;
	temp=A;
	A=B;
	B=temp;
}

void swap_double(double &A,double &B)
{
	double temp;
	temp=A;
	A=B;
	B=temp;
}

void swap_particle(particle &p1,particle &p2)
{
	swap_double(p1.u,p2.u);swap_double(p1.v,p2.v);swap_double(p1.w,p2.w);
	swap_double(p1.uTsum,p2.uTsum);	swap_double(p1.vTsum,p2.vTsum);	swap_double(p1.wTsum,p2.wTsum);
	swap(p1.uflag,p2.uflag);swap(p1.vflag,p2.vflag);swap(p1.wflag,p2.wflag);
	swap_double(p1.xuv,p2.xuv);swap_double(p1.yuv,p2.yuv);
	swap_double(p1.xuw,p2.xuw);swap_double(p1.yuw,p2.yuw);
	swap_double(p1.xvw,p2.xvw);swap_double(p1.yvw,p2.yvw);
	swap_double(p1.ruv,p2.ruv);swap_double(p1.phiuv,p2.phiuv);
	swap_double(p1.ruw,p2.ruw);swap_double(p1.phiuw,p2.phiuw);
	swap_double(p1.rvw,p2.rvw);swap_double(p1.phivw,p2.phivw);
	swap_double(p1.x,p2.x);	swap_double(p1.y,p2.y);
	swap_double(p1.r,p2.r);	swap_double(p1.phi,p2.phi);
	swap_double(p1.euv,p2.euv);swap_double(p1.euw,p2.euw);swap_double(p1.evw,p2.evw);
	swap_double(p1.auv,p2.auv);swap_double(p1.auw,p2.auw);swap_double(p1.avw,p2.avw);
	swap_double(p1.energy,p2.energy);
	swap_double(p1.angle,p2.angle);
	swap(p1.whichtwolayer,p2.whichtwolayer);
}

int add_dis(int *des,int *org,int n)
{
	int i;
	for(i=0;i<n;i++)
		des[i]+=org[i];
	return i;
}
int add_dis2D(int **des,int **org,int n1,int n2)
{
	int i,j;
	for(i=0;i<n1;i++)
		for(j=0;j<n2;j++)
			des[i][j]+=org[i][j];
	return i*j;
}
