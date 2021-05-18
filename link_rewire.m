function [ii,ii_abundance,jj_old,jj_new,degrees,theta, gamma_ap]=link_rewire(abundances_vector,na, ntotal, degrees, omegam,theta, Hap)
    %Randomly select a specie across guilds and perform a link rewire
    cont=false;
    while (cont==false)

        ii=randi([1,ntotal],1); % select specie i to rewire
        ii_abundance=abundances_vector(end,ii); %check abundance in case of link recovery
        if ((degrees(ii)>1) && (degrees(ii)<na) && (degrees(ii)<(ntotal-na)) )
            if ii<=na % spcecie i is animal
                possible_links_indices=find(theta(ii,:) == 0) + na;
                jj_new=datasample(possible_links_indices,1); % specie j' to add new link
                existing_links_indices=find(theta(ii,:) == 1)+ na;
                jj_old=datasample(existing_links_indices,1); % specie j to remove link with p_ij
            else % specie is plant
                possible_links_indices=find(theta(:,(ii-na)) == 0);
                jj_new=datasample(possible_links_indices,1); % specie j' to add new link
                existing_links_indices=find(theta(:,(ii-na)) == 1);
                jj_old=datasample(existing_links_indices,1); % specie j to remove link with p_ij
            end
            p_ij=(1 - (degrees(jj_old))^-1); % rewiring probability
            cont=p_ij>rand;
            if cont==true
                if ii<=na
                    theta(ii,(jj_old-na))=0;
                    theta(ii,(jj_new-na))=1;
                else
                    theta(jj_old,(ii-na))=0;
                    theta(jj_new,(ii-na))=1;
                end
                gamma_ap=omegam.*theta.*Hap; % updating gamma_ik
                k_old=degrees(jj_old) - 1;
                degrees(jj_old)=k_old; % updating the degrees of j and j'
                k_new=degrees(jj_new) + 1;
                degrees(jj_new)=k_new;
            end
        end
    end
end
