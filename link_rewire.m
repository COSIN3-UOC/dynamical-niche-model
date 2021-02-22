function [ii,ii_abundance,jj_old,jj_new,degrees,theta_ik, gammaik]=link_rewire(abundances_vector,na, ntotal, theta_ik, gammaik, omegam, Hapt)
    %Randomly select a specie across guilds and perform a link rewire
    cont=false;
    degrees=sum(theta_ik(:,:),1);
    while (cont==false)

        ii=randi([1,ntotal],1); % select specie i to rewire
        ii_abundance=abundances_vector(end,ii); %check abundance in case of link recovery
        if ((degrees(ii)>1) && (length(find(theta_ik(ii,na+1:end) == 0))>1) && (length(find(theta_ik(ii,1:na) == 0))>1) )
            if ii<=na % spcecie i is animal
                possible_links_indices=find(theta_ik(ii,na+1:end) == 0) + na;
                jj_new=datasample(possible_links_indices,1); % specie j' to add new link
                existing_links_indices=find(theta_ik(ii,na+1:end) == 1)+ na;
                jj_old=datasample(existing_links_indices,1); % specie j to remove link with p_ij
            else % specie is plant
                possible_links_indices=find(theta_ik(ii,1:na) == 0);
                jj_new=datasample(possible_links_indices,1); % specie j' to add new link
                existing_links_indices=find(theta_ik(ii,1:na) == 1);
                jj_old=datasample(existing_links_indices,1); % specie j to remove link with p_ij
            end
            p_ij=(1 - (degrees(jj_old))^-1); % rewiring probability
            cont=p_ij>rand;
            if cont==true
                theta_ik(ii,jj_old)=0;
                theta_ik(jj_old,ii)=0;
                theta_ik(ii,jj_new)=1;
                theta_ik(jj_new,ii)=1;
                gammaik=omegam.*theta_ik.*Hapt; % updating gamma_ik
                k_old=degrees(jj_old) - 1;
                degrees(jj_old)=k_old; % updating the degrees of j and j'
                k_new=degrees(jj_new) + 1;
                degrees(jj_new)=k_new;
            end
        end
    end
end