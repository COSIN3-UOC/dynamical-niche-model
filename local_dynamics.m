% Describing the dynamics
function nt= local_dynamics(t,ni,rho_a, rho_p, na, ntotal, h, theta_ik, gammaik, betaik)
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

    nt=zeros(1,ntotal);
    for i=1:ntotal
        if (i <= na) % for users (animals) especies
            competition_a=0.;
            for j=1:na
                competition_a=competition_a+betaik(i,j)*ni(j);
            end
            num=0.;
            theta_nk=0.;
            for k=na+1:ntotal
                num=num+gammaik(i,k)*ni(k);
                theta_nk=theta_nk+theta_ik(i,k)*ni(k);
            end
            mutualism_a=num /(1. + h*theta_nk);
            n_animal=ni(i)*(rho_a - competition_a + mutualism_a);
            nt(i)=n_animal; 
        else %for hashtags (plants)
            competition_p=0.;
            for j=na+1:ntotal
                competition_p=competition_p+betaik(i,j)*ni(j);
            end
            nump=0.;
            theta_nkp=0.;
            for k=1:na
                nump=nump+gammaik(i,k)*ni(k);
                theta_nkp=theta_nkp+theta_ik(i,k)*ni(k);
            end
            %lotka-volterra equations
            mutualism_p=nump /(1. + h*theta_nkp);
            n_plant=ni(i)*(rho_p - competition_p + mutualism_p);
            nt(i)=n_plant;
        end
    end
    nt=nt';
end