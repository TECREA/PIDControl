/*******************************************************************************
 *	pid_control - A library for basic PID control and simulatión of continous 
 *                      and discrete LTI transfer functions
 * 
 *	Copyright (C) 2013 Eng. Juan Camilo Gómez C. MSc. (kmilo17pet@gmail.com)
 *              Jaime Isaza Cadavid Colombian Polytechnic 
 *
 *	pid_control is free software: you can redistribute it and/or modify it
 *	under the terms of the GNU Lesser General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	pid_control is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public License
 *	along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
*******************************************************************************/

#ifndef PID_CONTROL_H
#define	PID_CONTROL_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <math.h>
    
typedef struct {
    double kc,ki,kd,dt,umin,umax,epsilon,kw,e1,ie,u1;
    int iflag,FORM;
}CPID;

typedef struct {
	double num[50],den[51],ymax,ymin;
	double _yreg[50], _ureg[50];
	int iflag;
}DTF;

typedef struct{
	double num[51], den[51],ymax,ymin,dt;
	double _x[50], _dx[50], _a0;
	int iflag,order;
}CTF;


extern double __ctf_tras_fcn_eval(CTF *TF,const double ut);
#define ctransferfcn(TF,UT)		__ctf_tras_fcn_eval(&TF,UT)

extern void __reg_update_fcn(double regarray[],double newp,int length);
#define update_reg(R,NEWVAL)  __reg_update_fcn(R,NEWVAL,(sizeof(R)/sizeof(R[0])))   
extern double __pid_controller(CPID *CT, const double w,const double y);
#define     pid_control(C,W,Y)    __pid_controller(&C,W,Y);
extern double __dtf_filt__fcn(DTF *TF,const double uk);
#define dtransferfcn(TF,UK)		__dtf_filt__fcn(&TF,UK)



#ifdef	__cplusplus
}
#endif

#endif	/* PID_CONTROL_H */

