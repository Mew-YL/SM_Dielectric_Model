function [dcsr,aa]=Wang(f,T,vwc,vsand,vclay,bd)

    if nargin<6
        error('Not enough input arguments.');
    elseif nargin>6
        error('Too many input arguments.');
    else
       wp = 0.06774 - 0.00064 * vsand + 0.004788 * vclay;
       mvt=0.49*wp+0.165;% unit% 最大结合水
       P=1-bd/2.65;
       r=0.481-0.57*wp;
       e0b=88.045-0.4147*T+6.295 * 10^(-4)*T^2+1.075 * 10^(-5)*T^3;
       t1b=(1.1109 * 10^(-10)-3.824*10^(-12)*T+6.938 * 10^(-14)*T^2-5.096 * 10^(-16)*T^3)/(2*pi);
       eb1=4.9+(e0b-4.9)/(1+(2*pi*f*t1b).^2);
       
       if vwc<=mvt
          ebw=3.2+(eb1-3.2)*vwc/mvt*r;
          dcsr=vwc*ebw+(P-vwc)+(1-P)*5.5;
       else
          ebw=3.2+(eb1-3.2)*r;
          dcsr=vwc*ebw+(P-vwc)+(1-P)*5.5+(vwc-mvt)*ebw;   
       end
       
     end
end