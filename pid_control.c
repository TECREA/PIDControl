#include "pid_control.h"

/*===========================================================================================================*/
double __pid_controller(CPID *CT, const double w,const double y){
    static double e=0,u=0,v=0,de=0;
    static double td=0,ti=0,temp=0;
    if ((CT->iflag)!=-1){
        CT->e1=0;
        CT->ie=0;
        CT->u1=0;
        CT->iflag=-1;
        (CT->dt)=((CT->dt)<=0)? 1.0  : (CT->dt);
        (CT->kw)=((CT->kw)<=0)? 1.0  : (CT->dt);
    }

    if (CT->FORM==1){ //conversion serie-paralelo
        ti=(CT->kc)/(CT->ki);
        td=(CT->kd)/(CT->kc);
        temp=(1+(td/ti));
        (CT->kc)=(CT->kc)*temp;
        (CT->ki)=(CT->kc)/(ti*temp);
        (CT->kd)=(CT->kc)*(td/temp);
    }
    e=w-y;
    if(fabs(e)<=(CT->epsilon)) e=0;
    (CT->ie)+=(e+(CT->u1))*(CT->dt); //parte integral con antiwindup de seguimiento
    de = (e-(CT->e1))/(CT->dt);
    v  = (CT->kc)*e + (CT->ki)*(CT->ie) + (CT->kd)*de;
    u  = (v>(CT->umax))? (CT->umax) : ( (v<(CT->umin))? (CT->umin): v );
    (CT->u1)=(CT->kw)*(u-v);
    (CT->e1)=e;
    return u;
}

/*===========================================================================================================*/
double __dtf_filt__fcn(DTF *TF,const double uk){
	int k;
	double yk=0,normt=TF->den[0];
	
	if ((TF->iflag)!=-1){
	if (( (TF->ymin)==(TF->ymax) ) && ( (TF->ymin)==0 ) ){
		(TF->ymin)=-1E100;
		(TF->ymax)= 1E100;
	}
		(TF->iflag)=-1;
	
	}
	
	update_reg(TF->_ureg, uk);
	for(k=0;k<50;k++){ yk+=(((TF->_ureg[k])*(TF->num[k]  )) + ((TF->_yreg[k])*(TF->den[k+1])) ); }
	yk=(normt==1)? yk : yk/normt;
        update_reg(TF->_yreg,-yk);
	return (yk>(TF->ymax))? (TF->ymax) : ( (yk<(TF->ymin))? (TF->ymin) : yk );
}

/*===========================================================================================================*/
void __reg_update_fcn(double regarray[],double newp,int length){
	int i;
        for (i=length-1;i>=1;i=i-1)
            regarray[i]=regarray[i-1];
        regarray[0]=newp;
}
/*=================================================================================*/
double __ctf_tras_fcn_eval(CTF *TF,const double ut){
	double temp=0;
	int k,n;
	if ((TF->iflag)!=-1){
		if (( (TF->ymin)==(TF->ymax) ) && ( (TF->ymin)==0 ) ){
			(TF->ymin)=-1E100;
			(TF->ymax)= 1E100;
		}

		(TF->_a0)=1.0*(TF->den[0]);
		for(n=50;n>0;n--){
			if ((TF->den[n])!=0) break;
		}
		
		for(k=0;k<=n;k++){
			(TF->num[k])/=(double)(TF->_a0);
			(TF->den[k])/=(double)(TF->_a0);
		}
		(TF->order)=n;
		(TF->iflag)=-1;
		(TF->dt)=((TF->dt)<=0)? 0.01  : (TF->dt);
	}
	n=(TF->order);
	/*-------------------------*/
	for(k=0;k<n;k++){
		temp+=(TF->_x[k])*-(TF->den[k+1]);
		(TF->_dx[k])=(TF->_x[k]);
	}
	update_reg((TF->_dx),(temp+ut));
	temp=0.0;
	for(k=0;k<n;k++){
		(TF->_x[k])+=(TF->_dx[k])*(TF->dt);
		temp+=(TF->_x[k])*((TF->num[k+1])-(TF->den[k+1])*(TF->num[0]) );
	}
	temp+=(TF->num[0])*ut;
	return (temp>(TF->ymax))? (TF->ymax) : ( (temp<(TF->ymin))? (TF->ymin) : temp );
}


