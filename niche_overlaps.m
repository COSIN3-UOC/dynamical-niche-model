%computing the niche overlaps
classdef niche_overlaps
%{
Compute the niche overlaps pre and post events
%}  
    methods(Static)
        
        function Hapt = pre_event_overlap(Ht,ssigma)
        %{
        Ht: niches center for all the especies (a+p) distributed in different
        ssigma: width of the gaussian niche profile of species. 

        Hapt: Niche overlap of all espcies. The overlap is computed in general form, assuming different width to each
        especie
        %}
            Hapt=exp(-((Ht-Ht').^2./(2.*(ssigma.^2+ssigma'.^2)))).*((2.*(ssigma.*ssigma'))./(ssigma.^2 + ssigma'.^2)).^(1/2);
        end
        
        function Hij_event = sudden_event_overlap(Hapt,Hee,Ha,Hp,na,ssigmaa,tau,alpha)
        %{
        Hapt: niche overlaps pre event
        
        Hij_event: Niche overlap of all especies with sudden event. The overlap is computed 
        in general form, assuming different width to each especie
        %}
            Hij_event=Hapt;
            f_event_uu=exp(-(2*alpha*tau))+exp(-(alpha*tau)).*( exp(-(Ha-Hee).^2./(4.*ssigmaa.^2)) + exp(-(Ha'-Hee').^2./(4.*ssigmaa.^2)) );
            f_event_uh=(exp(-(Ha-Hp').^2./(4.*ssigmaa.^2)))*(1-exp(-(alpha*tau))) + (exp(-(alpha*tau)).*( exp(-(Ha-Hee).^2./(4.*ssigmaa.^2))));
            Hij_event(1:na,1:na)=((Hapt(1:na,1:na))*(1-exp(-(alpha*tau)))+f_event_uu)/(exp(-(2*alpha*tau))+2*exp(-(alpha*tau))+1);
            Hij_event(1:na,na+1:end)=(f_event_uh)'/(exp(-(alpha*tau))+1);
            Hij_event(na+1:end,1:na)=(f_event_uh)/(exp(-(alpha*tau))+1);
        end
        
        function Hij_event = expected_event_overlap(Hapt,Hee,Ha,Hp,na,ssigmaa,tau_o,tau,alpha,aa)
        %{
        Hapt: niche overlaps pre event

        Hij_event: Niche overlap of all especies with expected event. The overlap is computed 
        in general form, assuming different width to each especie
        %}
            Hij_event=Hapt;
            f_event_uu=((1/(1 + ((tau-tau_o)/alpha)^(2*aa))^2) + (1/(1 + ((tau-tau_o)/alpha)^(2*aa))).*( exp(-(Ha-Hee).^2./(4.*ssigmaa.^2)) + exp(-(Ha'-Hee').^2./(4.*ssigmaa.^2))));
            f_event_uh=(exp(-(Ha-Hp').^2./(4.*ssigmaa.^2)))*(1-(1/(1 + ((tau-tau_o)/alpha)^(2*aa)))) + (1/(1 + ((tau-tau_o)/alpha)^(2*aa))).*( exp(-(Ha-Hee).^2./(4.*ssigmaa.^2)));
            Hij_event(1:na,1:na)=((Hapt(1:na,1:na))*(1-(1/(1 + ((tau-tau_o)/alpha)^(2*aa))))+f_event_uu)/((1/(1 + ((tau-tau_o)/alpha)^(2*aa)))^2+2*(1/(1 + ((tau-tau_o)/alpha)^(2*aa)))+1);
            Hij_event(1:na,na+1:end)=(f_event_uh)'/((1/(1 + ((tau-tau_o)/alpha)^(2*aa)))+1);
            Hij_event(na+1:end,1:na)=(f_event_uh)/((1/(1 + ((tau-tau_o)/alpha)^(2*aa)))+1);
        end
    
    end
end