mydir='./results/'; %path to save files
n_re=10;
parfor rep=1:n_re % to execute different realizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    omegam=.07;   % intensity of mutualism
    omegac=.07;  % intensity of competition
    llambda=.6; % inter-intra competition
    na=100; % number of animals
    np=100; % number of plants
    ntotal=na+np;
    Con=4/ntotal^0.8; %connectance

    %intrisic growth rate animals (users) and plants (hashtags). In case you want to assing different rho's to different guilds
    rho_a =1; 
    rho_p =1; 

    n0= 0.2; % initial abundances
    n_topics=4; %topics along the niche axis
    ht=0.1; %handling factor
    T=100; % integration time within time step
    t_max=30000;

    % Niche widths. General definition in case you want to assign a different width to each specie
    ssigmap= ones(1,np).*0.1; 
    ssigmaa=ones(1,na).*0.1;
    ssigma=[ssigmaa, ssigmap];

    % Creating the adjacency matrix theta
    theta=int16(rand(na,np)<Con);
    theta_ik=zeros(na+np,na+np);
    theta_ik(1:na,na+1:end)=theta;
    theta_ik(na+1:end,1:na)=theta';
    out_matrix=theta_ik(1:na,na+1:end);
    parsave_functions.parsave_matrices(fullfile(mydir,sprintf('matrix_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f_time_%08d.mat',rep,llambda,omegam,omegac,0)),out_matrix);

    % choosing the center positions of the topics
    H_n_u=[0.1,0.35,0.6,0.83]; % users
    H_n_h=[0.1,0.35,0.6,0.83]; % hashtags
    %assigning the niches positions
    Ht=niche_topics.pre_event_niche_topics(n_topics,na,np,H_n_u,H_n_h);
    % computing the pre event niches overlaps
    Hapt=niche_overlaps.pre_event_overlap(Ht,ssigma);

    % Creating the mutualistic matrix gamma
    gammaik=omegam.*theta_ik.*Hapt;
    % Creating the competitive matrix beta
    beta=zeros(na+np,na+np);
    beta(1:na,1:na)=1;
    beta(na+1:end,na+1:end)=1;
    betaik=omegac.*(llambda-llambda.*Hapt+(1-llambda).*Hapt).*beta;
    for i=1:na+np
        betaik(i,i)=1;
    end

    % integrating the initial dynamics
    fprintf('initial integration for realization %i\n',rep);
    ni0=ones(1,ntotal).*n0; % initial abundances
    [t0,y0] = ode45(@(t,ni)local_dynamics(t,ni,rho_a, rho_a, na, ntotal, ht, theta_ik, gammaik, betaik),[0 T],ni0);
    tTot=t0;
    yTot=y0;
    parsave_functions.parsave_abundances(fullfile(mydir,sprintf('results_abundances_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f_time_%08d.mat',rep,llambda,omegam,omegac,tTot(end))),tTot,yTot);

    % rewiring
    rewire_counts=0;
    iteration=1;
    fprintf('starting rewirings for realization %i\n',rep);
    tic;
    while (iteration<t_max)
        %fprintf('iteration %i for realization %i\n',iteration, rep);
        [ii,ii_abundance,jj_old,jj_new,degrees,theta_ik, gammaik]=link_rewire(yTot,na, ntotal, theta_ik, gammaik, omegam, Hapt); % rewire link function
        [ti,yi] = ode45(@(t,ni)local_dynamics(t,ni,rho_a, rho_p, na, ntotal, ht, theta_ik, gammaik, betaik),[0 T],yTot(end,:)); %integrate dynamics
        ti=ti+(T*iteration);
        tTot=ti(2:end);
        yTot=yi(2:end,:);
        new_abundance_ii=yi(end,ii);
        if (new_abundance_ii<=ii_abundance) %ineffective rewire attemp, recover old link      
            %fprintf('ineffective rewiring\n');
            theta_ik(ii,jj_old)=1;
            theta_ik(jj_old,ii)=1;
            theta_ik(ii,jj_new)=0;
            theta_ik(jj_new,ii)=0;
            gammaik=omegam.*theta_ik.*Hapt; % updating gamma_ik
            k_old=(degrees(jj_old) + 1 );
            degrees(jj_old)=k_old; % updating the degrees of j and j'
            k_new=(degrees(jj_new) - 1 );
            degrees(jj_new)=k_new;
        end

         if (rewire_counts==1000) % save data every 1000 successful rewires
            out_matrix=theta_ik(1:na,na+1:end);
            parsave_functions.parsave_matrices(fullfile(mydir,sprintf('matrix_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f_time_%08d.mat',rep,llambda,omegam,omegac,tTot(end))), out_matrix);
            parsave_functions.parsave_abundances(fullfile(mydir,sprintf('results_abundances_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f_time_%08d.mat',rep,llambda,omegam,omegac,tTot(end))),tTot,yTot);
            rewire_counts=0;
        end
        iteration=iteration+1;
        rewire_counts=rewire_counts+1;
    end
    %saving additional data needed to run post event dynamics
    parsave_functions.parsave_matrices(fullfile(mydir,sprintf('matrix_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f_time_%08d.mat',rep,llambda,omegam,omegac,tTot(end))), out_matrix);
    parsave_functions.parsave_abundances(fullfile(mydir,sprintf('results_abundances_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f_time_%08d.mat',rep,llambda,omegam,omegac,tTot(end))),tTot,yTot);
    parsave_functions.parsave_niche_overlaps(fullfile(mydir,sprintf('niche_matrix_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f.mat',rep,llambda,omegam,omegac)),Hapt);
    parsave_functions.parsave_niche_centers(fullfile(mydir,sprintf('positions_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f.mat',rep,llambda,omegam,omegac)),Ht(1:na),Ht(na+1:end));
    fprintf('end rewirings for realization %i\n',rep);
    elapsed_time=toc;
    fprintf('elapsed time %i for realization %i\n',elapsed_time, rep);
end