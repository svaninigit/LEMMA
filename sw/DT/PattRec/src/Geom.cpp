#include "Geom.h"


bool debug_geom  = false;


Geom::~Geom()
{   
  return;
}

float Geom::get_x_wire(int CH, int SL, int L, int W)
{
  float _x = 0;
  if( (CH==11 && (SL==1||SL==3)) 
      || CH==9 ||CH==8){
    if(L==1 || L==3 || L==5 || L==7)
      _x= (W-1)*4.2-148.6;
    if(L==2 || L==4 || L==6 || L==8)
      _x= (W-1)*4.2-150.7;
  }
  if( CH==10 && (SL==1||SL==3) ){
    if(L==1 || L==3 || L==5 || L==7)
      _x=-(W-1)*4.2+148.6;
    if(L==2 || L==4 || L==6 || L==8)
      _x=-(W-1)*4.2+150.7;
  }
  if( (CH==11||CH==10) && (SL==2) ){
    if(L==1 || L==3)
      _x= 117.35 - (W-1)*4.2;
    if(L==2 || L==4)
      _x= 119.45 - (W-1)*4.2;
  }
  return _x;
}
  
float Geom::get_y_wire(int CH, int SL, int L, int W)
{
  float _y = 0;
  float dy_ch1 = 0.1581;	//chamber 1 honeycomb correction

  float ay_phi_SL1[4] = {10.75, 9.45, 8.15, 6.85};
  float ay_phi_SL3[4] = {12.85, 14.15, 15.45, 16.75};
  float ay_theta[4] = {7.45, 8.75, 10.05, 11.35};
  float ay_phi_34[4] = {-1.95,-0.65,0.65,1.95};
  
  if(CH==11 && SL==1)
    for(int i=0;i<4;i++) _y=-ay_phi_SL1[L-1]-dy_ch1;
  if(CH==11 && SL==3)
    for(int i=0;i<4;i++) _y=ay_phi_SL3[L-1];
  if(CH==11 && SL==2)
    for(int i=0;i<4;i++) _y=ay_theta[L-1];

  if(CH==10 && SL==1)
    for(int i=0;i<4;i++) _y=ay_phi_SL1[L-1];
  if(CH==10 && SL==3)
    for(int i=0;i<4;i++) _y=-ay_phi_SL3[L-1];
  if(CH==10 && SL==2)
    for(int i=0;i<4;i++) _y=-ay_theta[L-1];

  if( CH==8||CH==9 )
    for(int i=0;i<4;i++) _y=ay_phi_34[L-1];
  
  return _y;
}

void Geom::printGeom()
{}

