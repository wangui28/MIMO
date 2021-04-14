function PopObj_up =func_pso_up(a1,b1,GAIN_C2BS,GAIN_C_UP,GAIN_C_DOWN,GAIN_D2BS,GAIN_BS2D,GAIN_D_UP,GAIN_D_DOWN,GAIN_C2D,GAIN_D2C)
global CUE DUE
N=size(a1,1);
PopObj_up=zeros(N,1);
for i=1:N
p_c=a1(i);
p_d=b1(i);
N0=3.981e-13;
PopObj_up(i)=log2((1+GAIN_C2BS(CUE)*p_c/(GAIN_D2BS(DUE)*p_d+N0)).*(1+GAIN_D_UP(DUE)*p_d/(GAIN_C2D(CUE,DUE)*p_c+N0)));
if (GAIN_C2BS(CUE)*p_c/(GAIN_D2BS(DUE)*p_d+N0)<5) || (GAIN_D_UP(DUE)*p_d/(GAIN_C2D(CUE,DUE)*p_c+N0)<5)  PopObj_up(i)=0;end
end

end
   
   
   
