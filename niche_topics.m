classdef niche_topics
%{
Assign the pre and post event niche topics along the niche axis
%}  
    methods(Static)

        function Ht= pre_event_niche_topics(n_topics,na,np,H_n_u,H_n_h)
        %{
        inputs:
        nt: number of topics
        na,np: number of users (animals) and hashtags (plants), respectively
        H_n_u, H_n_h: center positions of each topic

        outputs:
        Ht: a vector containing the center positions of the niches for each especie
        %}   
            %distribute the list of topics among the users randomly
            rng(1); topics=randi(n_topics,1,na);
            [especies_per_topic,~]=hist(topics,unique(topics));

            %distribute the list of topics among the hashtags randomly
            rng(5); topics_h=randi(n_topics,1,np);
            [hashtags_per_topic,~]=hist(topics_h,unique(topics_h));

            Ha_=[];
            Hp_=[];
            for i=1:n_topics
                Ha_=[Ha_, H_n_u(i).*ones(1,especies_per_topic(i))];
                Hp_=[Hp_, H_n_h(i).*ones(1,hashtags_per_topic(i))];
            end

            %Users: adding noise around the niche topics positions
            noise_niche = (.1 .* Ha_).* rand(1, length(Ha_)); 
            Ha = Ha_ + noise_niche;

            %Hashtags: % adding noise around the niche topics positions
            noise_niche_h = (.12 .* Hp_).* rand(1, length(Hp_)); 
            Hp = Hp_ + noise_niche_h;

            Ht=[Ha,Hp];
        end
        
        function Hee= event_niche_topic(na,ev_position)
            Hee=ones(1,na).*ev_position;
            noise_niche_event = (.1 .* Hee).* rand(1, length(Hee));
            Hee=Hee+noise_niche_event;
        end
        
    end
end