% Describing the dynamics
function nt= local_dynamics(t,ni,rho_a, rho_p, na, ntotal,theta, h,nt,num,den)
%{
Local dynamics: Lotka-Volterra equation with holling type II functional
response.The mutualistic and competitive interactions are proportional to the niche
overlaps.

ni: abundances of species
rho_a,rho_p: intrisic growth rates for especies guilds
na: number of animals (users in our context)
ntotal: total number of especies na+np (or nu+nh)
h: handling time
theta_ik: biadjacency matrix
gammaik: mutualistic matrix
betaik: competition matrix
%}

%     nt=zeros(ntotal,1);
%     num=zeros(ntotal,1);
%     den=zeros(ntotal,1);
    global gamma_ap;
    global beta_aa;
    global beta_pp;
    num(1:na) = gamma_ap*ni(na+1:ntotal);
    num(na+1:ntotal) = ((ni(1:na)')*gamma_ap)';
   
    
    den(1:na) = 1+h*theta*ni(na+1:ntotal);
    den(na+1:ntotal) = 1+h*(theta')*ni(1:na);
    mutualism = num./den;
   
    nt(1:na)=ni(1:na).*(rho_a - beta_aa*ni(1:na) + mutualism(1:na));
    nt(na+1:ntotal)=ni(na+1:ntotal).*( rho_p - beta_pp*ni(na+1:ntotal) + mutualism(na+1:end));
end
