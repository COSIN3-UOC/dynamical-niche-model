mydir='./results/'; %path to save files
n_re=10;
parfor rep=1:n_re % to execute different realizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    omegam=.07;   % intensity of mutualism
    omegac=.07;  % intensity of competition
    llambda=.6; % inter-intra competition
    %intrisic growth rate animals (users) and plants (hashtags). In case you want to assing different rho's to different guilds
    rho_a =1; 
    rho_p =1; 
    ht=0.1; %handling factor
    T=100; % integration time within time step
    event_t_max=70000;
    t_max=30000;
    
    % loading the data t run simulation with external event
    % you need to load the pre event final niche overlap and adjacency
    % matrices, final pre event abundances and niche positions.
    niche_matrix=load(fullfile(mydir,sprintf('niche_matrix_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f.mat',rep,llambda,omegam,omegac)));
    niche_positions=load(fullfile(mydir,sprintf('positions_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f.mat',rep,llambda,omegam,omegac)));
    abundances=load(fullfile(mydir,sprintf('results_abundances_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f_time_%08d.mat',rep,llambda,omegam,omegac,t_max*T)));
    theta=load(fullfile(mydir,sprintf('matrix_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f_time_%08d.mat',rep,llambda,omegam,omegac,t_max*T)));
    Hapt=niche_matrix.niche_matrix; %niche overlaps at pre-event
    
    %niche positions at pre-event
    Ha=niche_positions.Ha;
    Hp=niche_positions.Hp;
    theta=theta.out_matrix;
    
    %final abundances at pre-event
    yTot=abundances.yTot(end,:);
    tTot=abundances.tTot(end);
    na=length(Ha);
    np=length(Hp);
    ntotal=na+np;
    
    % Niche widths. General definition in case you want to assign a different width to each specie
    ssigmap= ones(1,np).*0.1; % niche width
    ssigmaa=ones(1,na).*0.1;
    ssigma=[ssigmaa, ssigmap];
    
    %final adjacency matrix at pre-event
    theta_ik=zeros(na+np,na+np);
    theta_ik(1:na,na+1:end)=theta(1:na,1:np);
    theta_ik(na+1:end,1:na)=theta(1:na,1:np)';
    degrees=sum(theta_ik(:,:),1);
    
    %% creating the niche overlap with external event
    % in case of a sudden event
    %parameters
%     tau=0;
%     alpha=0.0001;
%     ev_position=0.6;
%     Hee=niche_topics.event_niche_topic(na,ev_position);
%     %overlap
%     Hij_event=niche_overlaps.sudden_event_overlap(Hapt,Hee',Ha,Hp,na,ssigmaa,tau,alpha);
%     Hapt=Hij_event;
    
    % in case of an expected event
    %parameters
    tau=0;
    alpha=10000;
    tau_o=15000;
    aa=2;
    ev_position=0.6;
    Hee=niche_topics.event_niche_topic(na,ev_position);
    %overlap
    Hij_event=niche_overlaps.expected_event_overlap(Hapt,Hee,Ha,Hp,na,ssigmaa,tau_o,tau,alpha,aa);
    Hapt=Hij_event;
    
    % creating the mutualism matrix gamma
    gammaik=omegam.*theta_ik.*Hapt;
    % Creating the competitive matrix beta
    beta=zeros(na+np,na+np);
    beta(1:na,1:na)=1;
    beta(na+1:end,na+1:end)=1;
    betaik=omegac.*(llambda-llambda.*Hapt+(1-llambda).*Hapt).*beta;
    for i=1:na+np
        betaik(i,i)=1;
    end
% start rewirings after event
    iteration=t_max+1;
    rewire_counts=0;
    tic;
    fprintf('starting rewirings for realization %i\n',rep);
    while (iteration<event_t_max)
        [ii,ii_abundance,jj_old,jj_new,degrees,theta_ik, gammaik]=link_rewire(yTot,na, ntotal, theta_ik, gammaik, omegam, Hapt); % rewire link function
        [ti,yi] = ode45(@(t,ni)local_dynamics(t,ni,rho_a, rho_p, na, ntotal, ht, theta_ik, gammaik, betaik),[0 T],yTot(end,:));
        ti=ti+(T*iteration);
        tTot=ti(2:end);
        yTot=yi(2:end,:);
        new_abundance_ii=yi(end,ii);
        if (new_abundance_ii<=ii_abundance) %ineffective rewire attemp          
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
        %updating tau and recomputing overlap of event
        tau=tau+1;
        %Hij_event=niche_overlaps.sudden_event_overlap(Hapt,Hee,Ha,Hp,na,ssigmaa,tau,alpha); %sudden
        Hij_event=niche_overlaps.expected_event_overlap(Hapt,Hee,Ha,Hp,na,ssigmaa,tau_o,tau,alpha,aa); %expected
        Hapt=Hij_event;
        
        % updating the mutualism matrix gamma
        gammaik=omegam.*theta_ik.*Hapt;
        % updating the competitive matrix beta
        beta=zeros(na+np,na+np);
        beta(1:na,1:na)=1;
        beta(na+1:end,na+1:end)=1;
        betaik=omegac.*(llambda-llambda.*Hapt+(1-llambda).*Hapt).*beta;
        for i=1:na+np
            betaik(i,i)=1;
        end                   
        if (rewire_counts==1000)
            out_matrix=theta_ik(1:na,na+1:end);
            parsave_functions.parsave_matrices(fullfile(mydir,sprintf('matrix_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f_time_%08d.mat',rep,llambda,omegam,omegac,tTot(end))), out_matrix);
            parsave_functions.parsave_abundances(fullfile(mydir,sprintf('results_abundances_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f_time_%08d.mat',rep,llambda,omegam,omegac,tTot(end))),tTot,yTot);
            rewire_counts=0;
        end
        iteration=iteration+1;
        rewire_counts=rewire_counts+1;
    end
    parsave_functions.parsave_matrices(fullfile(mydir,sprintf('matrix_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f_time_%08d.mat',rep,llambda,omegam,omegac,tTot(end))), out_matrix);
    parsave_functions.parsave_abundances(fullfile(mydir,sprintf('results_abundances_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f_time_%08d.mat',rep,llambda,omegam,omegac,tTot(end))),tTot,yTot);
    parsave_functions.parsave_niche_overlaps(fullfile(mydir,sprintf('post_event_niche_matrix_rep_%.2d_lambda_%.2f_mutualism_%.2f_competition_%.2f.mat',rep,llambda,omegam,omegac)),Hapt);
    fprintf('end rewirings for realization %i\n',rep);
    elapsed_time=toc;
    fprintf('elapsed time %i for realization %i\n',elapsed_time, rep);
end