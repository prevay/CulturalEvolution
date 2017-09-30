function [potential_partners,preferred_ind_variant,divider] = get_potential_partners(i,agents,config,pref_locus,ind_locus,extra,n,range,adjmat)
%Function for selecting a set of candidates for interaction
%  Selects a set of agents who are potential interaction partners for the
%  given agent 'i' based on its preferences for tags.
    divider=[];
    preferred_ind_variant=[];
   if strcmp(config(7),'1')
      potential_partners=1:1:length(agents.indicators);

   elseif strcmp(config(2),'1') %ACT-R memory
      if n==1 || strcmp(config(7),'1') %just act randomly in the first generation or if unbiased config
        preferred_ind_variant=randi(range,1);
      else
        

        %calculate base-line activatin for each indicator assuming decay rate of 0.5
        pref_bla=log(2*agents.pref_n(i,:,pref_locus).*agents.pref_t(i,:,pref_locus).^(-0.5));
        %the 'recall anything above threshold' method:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        retrieved_indices=find(pref_bla>extra);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(retrieved_indices)
            preferred_ind_variant=randi(range,1); %pick a random tag if none are retrieved    
        else
            preferred_ind_variant=retrieved_indices(randsample(length(retrieved_indices),1)); %otherwise randomly select one from the retrieved ones
        end
      end

      %find all the agents with the selected tag in the right place:
      if strcmp(config(5),'1')|| strcmp(config(4),'0') 
          potential_partners=find(agents.indicators(:,ind_locus)==preferred_ind_variant);
      else
          potential_partners=find(agents.indicators(:,1)==preferred_ind_variant);
          divider=length(potential_partners);
          potential_partners=[potential_partners; find(agents.indicators(:,2)==preferred_ind_variant)];
      end
   else %naive memoery
      if strcmp(config(7),'1') %pick tags randomly in unbiased config
          preferred_ind_variant=randi(range,1);
          potential_partners=find(agents.indicators(:,ind_locus)==preferred_ind_variant);
      else
          memidx=randi(extra,1);%otherwise pick a random tag from those currently memorize
          preferred_ind_variant=agents.prefs(i,memidx,pref_locus);
          %find all those with the tag in the appropriate place:
          if strcmp(config(5),'1') || strcmp(config(4),'0')
              potential_partners=find(agents.indicators(:,ind_locus)==preferred_ind_variant);
          else
              potential_partners=find(agents.indicators(:,1)==preferred_ind_variant);
              divider=length(potential_partners);
              potential_partners=[potential_partners; find(agents.indicators(:,2)==preferred_ind_variant)];
          end
      end
   end

   %in the grid scenario only select those who are in the moore nbhd:
   if strcmp(config(3),'1')
      potential_partners=intersect(potential_partners,agents.nbhd(i)); 
   end
   
   %in the network scenario only select the subset of those who are linked to the agent:
   if strcmp(config(8),'1')
     potential_partners=intersect(potential_partners,find(adjmat(i,:)));   
   end



end

