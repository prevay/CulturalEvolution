function [partner,ind_locus,partner_ind_locus,agents] = get_actual_partner(i,agents,potential_partners,ind_locus,pref_locus,config,n,divider)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

       partner=[];
       partner_ind_locus=[];
       if ~isempty(potential_partners)
          potential_partner=potential_partners(randsample(length(potential_partners),1));  %randomly sample from the set of potential partners
          if strcmp(config(5),'0') && strcmp(config(4),'1')
              %special case for 1 preference and 2 tags--find if prefered
              %need to identify which tag we are looking at:
               idx=find(potential_partners==potential_partner);
               if idx<=divider
                   ind_locus=1;
               else
                   ind_locus=2;
               end 
          end
          
          %also need to identify which tag to look at in the 2 preference 2
          %tag scenario, using the agent's own preference to tag map:
          if strcmp(config(5),'1') && strcmp(config(4),'1')
              partner_ind_locus=agents.p2i(potential_partner,pref_locus);
          else    
              partner_ind_locus=1;
          end
                  
   	       if strcmp(config(7),'1')
		   partner=potential_partner; %in the unbiased config we are done
           %otherwise the chosen partner needs to accept:
	       elseif strcmp(config(2),'1') %ACT-R case:
               % so we check if the ego's tag is preferred by the chosen
               % partner:
              if strcmp(config(5),'1') || strcmp(config(4),'0') 
                  potential_partner_pref_bla=log(2*agents.pref_n(potential_partner,agents.indicators(i,partner_ind_locus),partner_ind_locus).*...
                      agents.pref_t(potential_partner,agents.indicators(i,partner_ind_locus)).^(-0.5));
                  potential_partner_bad_bla=log(2*agents.bad_n(potential_partner,agents.indicators(i,partner_ind_locus),partner_ind_locus).*...
                      agents.bad_t(potential_partner,agents.indicators(i,partner_ind_locus)).^(-0.5));
              else 
                  potential_partner_pref_bla1=log(2*agents.pref_n(potential_partner,agents.indicators(i,1)).*...
                      agents.pref_t(potential_partner,agents.indicators(i,1)).^(-0.5));
                  potential_partner_bad_bla=log(2*agents.bad_n(potential_partner,agents.indicators(i,1)).*...
                      agents.bad_t(potential_partner,agents.indicators(i,2)).^(-0.5));
                  potential_partner_pref_bla2=log(2*agents.pref_n(potential_partner,agents.indicators(i,2)).*...
                      agents.pref_t(potential_partner,agents.indicators(i,2)).^(-0.5));
                  potential_partner_bad_bla=[potential_partner_bad_bla log(2*agents.bad_n(potential_partner,agents.indicators(i,2)).*...
                      agents.bad_t(potential_partner,agents.indicators(i,2)).^(-0.5))];
                  [potential_partner_pref_bla,idx]=max([potential_partner_pref_bla1,potential_partner_pref_bla2]);
                  potential_partner_bad_bla=potential_partner_bad_bla(idx);
              end
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if potential_partner_pref_bla>potential_partner_bad_bla || n==1
                  partner=potential_partner; % if preference is mutual then interact
              end
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           else %naive memory case:
              if strcmp(config(5),'1') || strcmp(config(4),'0')
                  if ~isempty(find(agents.prefs(potential_partner,:,partner_ind_locus)==agents.indicators(i,partner_ind_locus),1))
                      partner=potential_partner; % if preference is mutual then interact
                  end
              else
                  if ~isempty(find(agents.prefs(potential_partner,:,1)==agents.indicators(i,1) || ...
                          agents.prefs(potential_partner,:,2)==agents.indicators(i,2),1))
                      partner=potential_partner; % if preference is mutual then interact
                  end
              end
           end
          
           %keep track of this interaction:
           agents.interactions(potential_partner)=agents.interactions(potential_partner)+1;
           if ~isempty(partner)
               agents.good_interactions(potential_partner)=agents.good_interactions(potential_partner)+1;
           end
      end

end

